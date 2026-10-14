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
 
                event_rate <- p
                X_overarch <- readRDS(paste("~/scratch/simulated_data/overarching/iter_", iter,".rds", sep=""))
                # version 2026 - generate intercept once per iteration
                b0 <- generate_intercept(N_pop = N, prob = p, lp = X_overarch  %*% beta[2:6])
                x1 <- lapply(1:K, function(k) { rnorm(n = N, mean = mean_x1[k], sd = 1) })
                out <- lapply(1:K, function(k) {
                        beta[2] <- b_x1[k]
                        
                        # generate X
                        X <- X_overarch
                        X[,1] <- x1[[k]]
                        colnames(X) <- c("x1", "x2", "x3", "x4", "x5")
                        
                        
                        # generate intercept
                        #lp <- X %*% beta[2:6]
                        # version 2025 
                        #b0 <- generate_intercept(N_pop=N, prob=p, lp=lp)
                        
                        X <- cbind(rep(1, N), X)
                        
                        beta[1] <- b0
                        
                        # generate y
                        lp <- as.vector(as.matrix(X) %*% beta)
                        p <- exp(lp)/(1+exp(lp))
                        y <- rbinom(N, 1, p)
                        
                        # Dichotomize X if type = "b"
                        
                        X_dicho <- match(c("x1", "x2", "x3", "x4", "x5"), colnames(X))
                        if (type == "b") {
                          X[, X_dicho] <- sapply(X_dicho, function(x) { X[, x] >  qnorm(0.9, 0, 1) } )
                        }
                        
                        # combine all
                        cluster <- cbind(X, lp=lp, p=p,y=y, cluster=k)
                        out = list(cluster=as.data.frame(cluster), beta=beta, mu_x1=mean_x1[[k]])

                        return(out)
  })
  paramfile <- paste("N",N,"iter", iter, "tau", tau, "tau_beta", tau_beta, "p", event_rate)
  paramfile <- paste("~/scratch/parameters/observed/",paramfile, ".rds", sep="")
saveRDS(lapply(out, function(x) x[c("beta", "mu_x1")]), file = paramfile)
return(do.call(rbind, lapply(out, function(x) x$cluster)))
}


generate_unobserved_clusters <- function(N, # number of observations in each cluster
                             K, # number of clusters
                             tau, # tau
                             tau_beta, # tau beta
                             p ,# expected event rate,
                             type, # predictor type
                             iter, # iteration
                             beta,
                             b_x1,
                             mean_u1
)
{
  # This function is for generating (number of K) clusters of the heterogeneous populations (fixed proxy)
  event_rate <- p
  X_overarch <- readRDS(paste("~/scratch/simulated_data/overarching/iter_", iter,".rds", sep=""))
  b0 <- generate_intercept(N_pop = N, prob = p, lp = X_overarch  %*% beta[2:6])
  u1 <- lapply(1:K, function(k) { rnorm(n = N, mean = mean_u1[k], sd = 1) })
  #P <- 5  # total predictors = 5
  out <- lapply(1:K, function(k) {
    beta[2] <- b_x1[k]

    # generate u and l
    X <- X_overarch
    X[,1] <- u1[[k]]
    colnames(X) <- c("u1", "x2", "x3", "x4", "x5")

    #lp <- X %*% beta[2:6]
    
    # generate intercept
    #b0 <- generate_intercept(N_pop=N, prob=p, lp=lp)
    X <- cbind(rep(1, N), X)
    beta[1] <- b0
    
    # generate y
    lp <- as.vector(as.matrix(X) %*% beta)
    p <- exp(lp)/(1+exp(lp))
    y <- rbinom(N, 1, p)

    # Dichotomize X if type = "b"
    X_dicho <- match(c("x2", "x3", "x4", "x5"), colnames(X))
    if (type == "b") {
      X[, X_dicho] <- sapply(X_dicho, function(x) { X[, x] >  qnorm(0.9, 0, 1) } )
    }

    # combine all
    cluster <- cbind(X,lp=lp, p=p,y=y, cluster=k)
    out = list(cluster=as.data.frame(cluster), beta=beta, mu_u1=mean_u1[[k]])
    
    return( out)
  })
  paramfile <- paste("N",N,"iter", iter, "tau", tau, "tau_beta", tau_beta, "p", event_rate)
  paramfile <- paste("~/scratch/parameters/observed/",paramfile, ".rds", sep="")
saveRDS(lapply(out, function(x) x[c("beta", "mu_u1")]), file = paramfile)
return(do.call(rbind, lapply(out, function(x) x$cluster)))
}


generate_proxy_clusters <- function(N, # number of observations in each cluster
                                       K, # number of clusters
                                       tau, # tau mean of U1
                                       tau_beta, # tau beta of U1-X1
                                       tau_sd, # tau sd of L1
                                       p, # expected event rate)
                                       type, # predictor type
                                       iter, # iteration
                                       beta,
                                       b_x1,
                                       mean_u1
)
{
  # This function is for generating (number of K) clusters of the heterogeneous populations (fixed proxy)
  event_rate <- p
  X_overarch <- readRDS(paste("~/scratch/simulated_data/overarching/iter_", iter,".rds", sep=""))
  b0 <- generate_intercept(N_pop = N, prob = p, lp = X_overarch  %*% beta[2:6])
  u1 <- lapply(1:K, function(k) { rnorm(n = N, mean = mean_u1[k], sd = 1) })
  #P <- 5  # total predictors = 5
  out <- lapply(1:K, function(k) {
    beta[2] <- b_x1[k]

    # generate u and l
    X <- X_overarch
    X[,1] <- u1[[k]]
    l1 <- rnorm(n=N, mean=u1[[k]], sd=tau_sd)
    colnames(X) <- c("u1", "x2", "x3", "x4", "x5")

    #lp <- X %*% beta[2:6]
    
    # generate intercept
    #b0 <- generate_intercept(N_pop=N, prob=p, lp=lp)
    X <- cbind(rep(1, N), X)
    beta[1] <- b0
    
    # generate y
    lp <- as.vector(as.matrix(X) %*% beta)
    p <- exp(lp)/(1+exp(lp))
    y <- rbinom(N, 1, p)

    # add l1
    X <- cbind(X, l1)

    # Dichotomize X if type = "b"
    X_dicho <- match(c("l1", "x2", "x3", "x4", "x5"), colnames(X))
    if (type == "b") {
      X[, X_dicho] <- sapply(X_dicho, function(x) { X[, x] >  qnorm(0.9, 0, 1) } )
    }

    # combine all
    cluster <- cbind(X,lp=lp, p=p,y=y, cluster=k)
    out = list(cluster=as.data.frame(cluster), beta=beta, mu_u1=mean_u1[[k]])
    
    return( out)
  })
  paramfile <- paste("N",N,"iter", iter, "tau", tau, "tau_beta", tau_beta, "tau_sd", tau_sd, "p", event_rate)
  paramfile <- paste("~/scratch/parameters/proxy/",paramfile, ".rds", sep="")
saveRDS(lapply(out, function(x) x[c("beta", "mu_u1")]), file = paramfile)
return(do.call(rbind, lapply(out, function(x) x$cluster)))
}