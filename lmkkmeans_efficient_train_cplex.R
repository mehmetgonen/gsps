library(Rcplex)

lmkkmeans_efficient_train <- function(Km, parameters) {
  state <- list()
  state$time <- system.time({
    N <- dim(Km)[2]
    P <- dim(Km)[3]
    Theta <- matrix(1 / P, N, P)
    rownames(Theta) <- dimnames(Km)[[1]]
    colnames(Theta) <- dimnames(Km)[[3]]
    K_Theta <- matrix(0, nrow(Km), ncol(Km))
    for (m in 1:P) {
      K_Theta <- K_Theta + (Theta[,m,drop = FALSE] %*% t(Theta[,m,drop = FALSE])) * Km[,,m]  
    }

    objective <- rep(0, parameters$iteration_count)
    dual_variables <- rep(0, N)
    for (iter in 1:parameters$iteration_count) {
      print(sprintf("running iteration %d...", iter))
      H <- eigen(K_Theta, symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
      HHT <- H %*% t(H)

      values <- sapply(1:P, function(m) {((diag(1, N, N) - HHT) * Km[,,m]) %*% Theta[,m] - dual_variables})
      values <- order(colMeans(abs(values)))
      m1 <- values[1]
      m2 <- values[length(values)]

      Q <- matrix(0, 2 * N, 2 * N)
      Q[1:N, 1:N] <- diag(1, N, N) * Km[,,m1] - HHT * Km[,,m1]
      Q[(N + 1):(2 * N), (N + 1):(2 * N)] <- diag(1, N, N) * Km[,,m2] - HHT * Km[,,m2]

      opts <- list()
      opts$trace <- 0
      result <- Rcplex(cvec = rep(0, 2 * N), 
                       Amat = matrix(rep(diag(1, N, N), 2), nrow = N, ncol = 2 * N), 
                       bvec = rep(1, N) - rowSums(Theta[, -c(m1, m2), drop = FALSE]), 
                       Qmat = Q, 
                       lb = rep(0, 2 * N), 
                       ub = rep(1, 2 * N), 
                       control = opts, 
                       objsense = "min", 
                       sense = rep("E", N))
      dual_variables <- result$extra$slack
      Theta[, c(m1, m2)] <- matrix(result$xopt, N, 2, byrow = FALSE)
      K_Theta <- matrix(0, nrow(Km), ncol(Km))
      for (m in 1:P) {
        K_Theta <- K_Theta + (Theta[,m,drop = FALSE] %*% t(Theta[,m,drop = FALSE])) * Km[,,m]  
      }

      objective[iter] <- sum(diag(t(H) %*% K_Theta %*% H)) - sum(diag(K_Theta))
    }
    H_normalized <- H / matrix(sqrt(rowSums(H^2, 2)), nrow(H), parameters$cluster_count, byrow = FALSE)

    set.seed(NULL)
    state$clustering <- kmeans(H_normalized, centers = parameters$cluster_count, iter.max = 1000, nstart = 10)$cluster
    names(state$clustering) <- dimnames(Km)[[1]]
    state$K_Theta <- K_Theta
    state$objective <- objective
    state$parameters <- parameters
    state$Theta <- Theta
  })
  return(state)
}