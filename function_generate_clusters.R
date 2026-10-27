source("function_generate_beta.R")
source("function_generate_intercept.R")
suppressPackageStartupMessages(library(MASS, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))



generate_observed_clusters <- function(N, # number of observations in each cluster
                             K, # number of clusters
                             tau, # tau
                             tau_beta, # tau beta
                             p, # expected event rate,
                             type, # predictor type
                             iter, # iteration
                             beta,
                             b_x1,
                             mean_x1
                             ) 
{
  # This function is for generating (number of K) clusters of the heterogeneous populations (observed causes)
  # for iteration i(ter)
    # Generate x2, x3, x4, x5
    x2 <- rnorm(n = N * K, mean = 0, sd = 1)
    x3 <- rnorm(n = N * K, mean = 0, sd = 1)
    x4 <- rnorm(n = N * K, mean = 0, sd = 1)
    x5 <- rnorm(n = N * K, mean = 0, sd = 1)
    
    # X1 according to X_mean from tau (for covariate shift)
    x1 <- matrix(nrow = N, ncol = K)
    for (k in 1:K) {
        # sample x1 from the corresponding mean x1
        x1[,k] <- rnorm(n = N, mean = mean_x1[k], sd = 1)
    }
    # combine Xs
    X <- data.frame(x1 = as.vector(x1), x2 = x2, x3 = x3, x4 = x4, x5 = x5, cluster = rep(1:K, each = N))

    b_x1_broadcasted <- rep(b_x1, each = N)
    beta_broadcasted <- matrix(rep(beta, each=N*K), ncol=length(beta))
    beta_broadcasted[,2] <- b_x1_broadcasted
    int_X <- cbind(1, X[, 1:5])
    lp <- as.numeric(rowSums(int_X * beta_broadcasted))

    #lp <- as.vector(as.matrix(cbind(1, X[,1:5])) %*% beta)
    p <- exp(lp)/(1+exp(lp))
    y <- rbinom(N*K, 1, p)
    clusters <- cbind(y, X)
    return(clusters)
}

generate_unobserved_clusters <- function(N, # number of observations in each cluster
                             K, # number of clusters
                             tau, # tau
                             tau_beta, # tau beta
                             p, # expected event rate,
                             type, # predictor type
                             iter, # iteration
                             beta,
                             b_x1,
                             mean_x1
                             ) 
{
  # This function is for generating (number of K) clusters of the heterogeneous populations (unobserved causes)
  # for iteration i(ter)
    # Generate x2, x3, x4, x5
    x2 <- rnorm(n = N * K, mean = 0, sd = 1)
    x3 <- rnorm(n = N * K, mean = 0, sd = 1)
    x4 <- rnorm(n = N * K, mean = 0, sd = 1)
    x5 <- rnorm(n = N * K, mean = 0, sd = 1)
    
    # X1 according to X_mean from tau (for covariate shift)
    x1 <- matrix(nrow = N, ncol = K)
    for (k in 1:K) {
        # sample x1 from the corresponding mean x1
        x1[,k] <- rnorm(n = N, mean = mean_x1[k], sd = 1)
    }
    # combine Xs
    # Now change x1 to u1
    X <- data.frame(u1 = as.vector(x1), x2 = x2, x3 = x3, x4 = x4, x5 = x5, cluster = rep(1:K, each = N))

    b_x1_broadcasted <- rep(b_x1, each = N)
    beta_broadcasted <- matrix(rep(beta, each=N*K), ncol=length(beta))
    beta_broadcasted[,2] <- b_x1_broadcasted
    int_X <- cbind(1, X[, 1:5])
    lp <- as.numeric(rowSums(int_X * beta_broadcasted))

    #lp <- as.vector(as.matrix(cbind(1, X[,1:5])) %*% beta)
    p <- exp(lp)/(1+exp(lp))
    y <- rbinom(N*K, 1, p)
    clusters <- cbind(y, X)
    return(clusters)
}


generate_proxy_clusters <- function(N, # number of observations in each cluster
                             K, # number of clusters
                             tau, # tau
                             tau_beta, # tau beta
                             tau_sd, # tau sd of L1
                             p, # expected event rate,
                             type, # predictor type
                             iter, # iteration
                             beta,
                             b_x1,
                             mean_x1
                             ) 
{
  # This function is for generating (number of K) clusters of the heterogeneous populations (proxy)
  # for iteration i(ter)
    # Generate x2, x3, x4, x5
    x2 <- rnorm(n = N * K, mean = 0, sd = 1)
    x3 <- rnorm(n = N * K, mean = 0, sd = 1)
    x4 <- rnorm(n = N * K, mean = 0, sd = 1)
    x5 <- rnorm(n = N * K, mean = 0, sd = 1)
    
    # X1 according to X_mean from tau (for covariate shift)
    x1 <- matrix(nrow = N, ncol = K)
    for (k in 1:K) {
        # sample x1 from the corresponding mean x1
        x1[,k] <- rnorm(n = N, mean = mean_x1[k], sd = 1)
    }
    # combine Xs
    # Now change x1 to u1
    X <- data.frame(u1 = as.vector(x1), x2 = x2, x3 = x3, x4 = x4, x5 = x5, cluster = rep(1:K, each = N))

    b_x1_broadcasted <- rep(b_x1, each = N)
    beta_broadcasted <- matrix(rep(beta, each=N*K), ncol=length(beta))
    beta_broadcasted[,2] <- b_x1_broadcasted
    int_X <- cbind(1, X[, 1:5])
    lp <- as.numeric(rowSums(int_X * beta_broadcasted))

    #lp <- as.vector(as.matrix(cbind(1, X[,1:5])) %*% beta)
    p <- exp(lp)/(1+exp(lp))
    y <- rbinom(N*K, 1, p)

    # add proxy l1
    l1 <- rnorm(n = N*K, mean = X$u1, sd = tau_sd)
    clusters <- cbind(y, X, l1 = as.vector(l1))
    return(clusters)
}