setwd("~/simulation")
source("function_utils.R")
source("function_generate_intercept.R")
source("function_generate_clusters.R")

### Initiate parameters

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


##############################
### Start generating datasets
##############################
start <- Sys.time()
title_txt <- list(N = N, tau = tau, tau_beta = tau_beta, tau_sd = tau_sd, p = p, K = K, seed = seed, type = type)
print(paste(Sys.time(), ":", title_txt))


set.seed(seed)
beta <- readRDS("~/simulation/simulated_data/betas.rds")

# Generate beta x1 vector with length of K
b_x1 <- sort(rnorm(K, mean = beta[2], sd = tau_beta))

# Generate X1 mean vector with length of K
mean_x1 <- sort(rnorm(K, mean=0, sd=tau))

saveRDS(mean_x1, paste("~/simulation/simulated_data/x1_tau",tau,".rds", sep="_"))
saveRDS(b_x1, paste("~/simulation/simulated_data/b1_tau_beta",tau_beta,".rds", sep="_"))


# Generate dataset iteration-by-iteration
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

### NO LONGER USED AS CURRENTLY IMPLEMENTING BOOTSTRAPPING
### proportion of development population for the benchmark model
dev_size <- 0.8
tpop_dev_index <- sample(1:N, size = dev_size * N)
saveRDS(tpop_dev_index, paste0("simulated_data/tpop_dev_index_N_", N, ".rds"))

end <- Sys.time()
total_time <- difftime(end, start, units = "mins")
print(paste("Total time:", round(total_time, 2), "minutes"))

