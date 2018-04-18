library(cluster)
optimizer <- "mosek"

if (optimizer == "cplex") {
  source("lmkkmeans_efficient_train_cplex.R") 
}
if (optimizer == "mosek") {
  source("lmkkmeans_efficient_train_mosek.R") 
}

cohorts <- c("ACC",  "BLCA", "BRCA", "CESC", "CHOL",
             "COAD", "DLBC", "ESCA", "GBM",  "HNSC",
             "KICH", "KIRC", "KIRP", "LAML", "LGG",
             "LIHC", "LUAD", "LUSC", "MESO", "OV",
             "PAAD", "PCPG", "PRAD", "READ", "SARC",
             "SKCM", "STAD", "TGCT", "THCA", "THYM",
             "UCEC", "UCS",  "UVM")

data_path <- "./data"
result_path <- "./results"

pdist <- function(X1, X2) {
  if (identical(X1, X2) == TRUE) {
    D <- as.matrix(dist(X1))
  }
  else {
    D <- as.matrix(dist(rbind(X1, X2)))
    D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
  }
  return(D)
}

for (cohort in cohorts) {
  for (pathway in c("hallmark", "kegg", "pid")) {
    for (cluster_count in 2:10)
      if (file.exists(sprintf("%s/state_cohort_%s_pathway_%s_cluster_%s.RData", result_path, cohort, pathway, cluster_count)) == FALSE) {
        print(sprintf("running %s cohort with %d clusters ...", cohort, cluster_count))
        load(sprintf("%s/TCGA-%s.RData", data_path, cohort))
        load(sprintf("./msigdb/%s.RData", pathway))
        N_pathway <- length(pathways)
        
        X_train <- log2(TCGA$mrna[unique(rownames(TCGA$mrna)),] + 1)
        X_train <- scale(X_train)
        X_train <- X_train[,colMeans(is.na(X_train)) == 0]
        
        pathway_names <- sapply(1:length(pathways), FUN = function(m) pathways[[m]]$name)
        N_train <- nrow(X_train)
        K_train <- array(0, dim = c(N_train, N_train, N_pathway), dimnames = list(rownames(X_train), rownames(X_train), pathway_names))
        for (m in 1:N_pathway) {
          hit_features <- which(colnames(X_train) %in% pathways[[m]]$symbols)
          D_train <- pdist(X_train[, hit_features, drop = FALSE], X_train[, hit_features, drop = FALSE])
          sigma <- mean(D_train)
          K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
        }
        
        parameters <- list()
        parameters$cluster_count <- cluster_count
        parameters$iteration_count <- 1000
        
        state <- lmkkmeans_efficient_train(K_train, parameters)
        
        all_genes <- unique(unlist(sapply(1:length(pathways), FUN = function(m) pathways[[m]]$symbols)))
        result <- silhouette(state$clustering, as.matrix(dist(X_train[,colnames(X_train) %in% all_genes], "euclidean"))^2)
        state$silhouette_score <- as.numeric(colMeans(result)[3])
        
        save("state", file = sprintf("%s/state_cohort_%s_pathway_%s_cluster_%s.RData", result_path, cohort, pathway, cluster_count))
      }
  }
}