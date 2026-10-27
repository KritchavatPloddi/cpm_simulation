setwd("~/simulation")
source("function_utils.R")
source("function_generate_intercept.R")
# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments with defaults
n_iter <- get_arg("n_iter", 1000)
N <- get_arg("N", 1000)
K <- get_arg("K", 10)
structure <- ifelse(any(grepl("--structure=", args)), sub("--structure=", "", args[grep("--structure=", args)]), 'observed')
tau <- get_arg("tau", 0)
tau_beta <- get_arg("tau_beta", 0)
tau_sd <- get_arg("tau_sd",0)
p <- get_arg("p", 0.2)
# Predictor type: c = continuous, b = binary
type <- ifelse(any(grepl("--type=", args)), sub("--type=", "", args[grep("--type=", args)]), "c")
seed <- get_arg("seed", 68250340)


source("function_generate_clusters.R")

start <- Sys.time()

title_txt <- list(N = N, tau = tau, tau_beta = tau_beta, tau_sd = tau_sd, p = p, K = K, seed = seed, type = type)

print(paste(Sys.time(), ":", title_txt))
# generate beta
#set.seed(seed)
#mean_beta <- 2
#mean_beta <- 1.2
#sd_beta <- 1
# simulate 5 betas plus 1 intercept with the uniform distribution (the intercept will be modified for the expected prevalence)
#beta <- rnorm(n = 6, mean = mean_beta, sd = sd_beta)
#beta <- runif(n = 6, min = 0, max = 2)

#### ### Generate overarching intercept
#### X_list <- lapply(1:n_iter, function(iter) {
####     readRDS(paste0("~/scratch/simulated_data/overarching/iter_", iter, ".rds"))
#### })
#### X <- do.call(rbind, X_list)
#### beta[1] <- generate_intercept(N_pop = dim(X)[1], prob = p, X %*% beta[2:6])
#### rm(X_list, X)
set.seed(seed)
beta <- readRDS("~/simulation/simulated_data/betas.rds")

### Old 2025 anaylsis
#b_x1 <- matrix(rnorm(n = K*n_iter, mean = beta[2], sd = tau_beta), ncol = K)

### New analysis Ranked b_x1  
### Target population would k = 1, 3, 5
### First sampling beta x1 from a normal distribution with mean = 1.35 (beta[2]), and sd = tau_beta) and sorted
#b_x1 <- sort(rnorm(K, mean = beta[2], sd = tau_beta))
### Broadcast sorted beta x1 to be N x K matrix
#b_x1 <- matrix(rep(b_x1,N), nrow = N, byrow = TRUE)

# Generate beta x1 vector with length of K
b_x1 <- sort(rnorm(K, mean = beta[2], sd = tau_beta))

# Generate X1 mean vector with length of K
mean_x1 <- sort(rnorm(K, mean=0, sd=tau))

saveRDS(mean_x1, paste("~/simulation/simulated_data/x1_tau",tau,".rds", sep="_"))
saveRDS(b_x1, paste("~/simulation/simulated_data/b1_tau_beta",tau_beta,".rds", sep="_"))


# Generate the n_iter x K matrix
#mean_x1 <- matrix(rnorm(K * n_iter, mean = 0, sd = tau), nrow = n_iter, ncol = K)

# 1. generate mean of cluster (subpopulation)-wise X1/U1/L1 matrix for KxN matrix, mean 0 and sd tau for covariate shift
#mean_x1_unsorted <- matrix(rnorm(K*n_iter, mean=0, sd=tau), ncol = K)
# 2. Sort each row individually from lowest to highest (left to right)
#mean_x1 <- t(apply(mean_x1_unsorted, 1, sort))
#set.seed(1)
for (iter in 1:n_iter) {
  #print(iter)
  fname_string <- paste("iter", iter, 
                        "N", title_txt$N, 
                        "tau", title_txt$tau,
                        "tau_beta", title_txt$tau_beta,
                        "tau_sd", title_txt$tau_sd,
                        "p", title_txt$p,
                        "K", title_txt$K,
                        "seed", title_txt$seed
                        )

  if (structure == "observed") {
    clusters <- generate_observed_clusters(iter = iter, N = N, K = K, tau = tau, tau_beta = tau_beta, p = p, type = type, 
                                          beta = beta, b_x1 = b_x1, mean_x1 = mean_x1)
  } else if(structure == "unobserved") {
    clusters <- generate_unobserved_clusters(iter = iter, N = N, K = K, tau = tau, tau_beta = tau_beta, p = p, type = type,
                                             beta = beta, b_x1 = b_x1, mean_x1 = mean_x1)
  } else if (structure == "proxy") {
    
    clusters <- generate_proxy_clusters(iter = iter, N = N, K = K, tau = tau, tau_beta = tau_beta, tau_sd = tau_sd, p = p, type = type, 
                                        beta = beta, b_x1 = b_x1, mean_x1 = mean_x1)
  } 
  
  #clusters$cluster <- ceiling(seq_len(nrow(clusters)) / (N * 10 / K))
  save_path <- paste("~/scratch/simulated_data/",type,"/", structure, "/", fname_string, ".rds", sep="")
  saveRDS(clusters, save_path)
  print(paste("Simulated data at:", save_path))
}

### proportion of development population for the benchmark model
dev_size <- 0.8
tpop_dev_index <- sample(1:N, size = dev_size * N)
saveRDS(tpop_dev_index, paste0("simulated_data/tpop_dev_index_N_", N, ".rds"))

end <- Sys.time()
total_time <- difftime(end, start, units = "mins")
print(paste("Total time:", round(total_time, 2), "minutes"))

