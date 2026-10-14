suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
setwd('~/simulation')
source("function_utils.R")
source("function_simulation_nrun.R")
source("function_simulation_nrun_2pop.R")

start_timestamp <- Sys.time()
total_cores = detectCores()
start_message <- paste(start_timestamp,": Simulation started, toal cores available -,", total_cores, "\n")
cat(start_message)
# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to get argument or set default
# Assign arguments with defaults
n_iter <- get_arg("n_iter", 500)
N <- get_arg("N", 1000)
K <- get_arg("K", 10)
structure <- ifelse(any(grepl("--structure=", args)), sub("--structure=", "", args[grep("--structure=", args)]), 'observed')
tau <- get_arg("tau", 0)
tau_beta <- get_arg("tau_beta", 0)
tau_sd <- get_arg("tau_sd", 0)
p <- get_arg("p", 0.2)
type <- ifelse(any(grepl("--type=", args)), sub("--type=", "", args[grep("--type=", args)]), "c")
seed <- get_arg("seed", 68250340)
cores <- get_arg("cores", 1)
dpop_num <- get_arg("dpop_num", 9)
run_all <- get_arg("run_all", 1)

if (structure == "observed") {
  f1 <- "y ~ x1 + x2 + x3 + x4 + x5"
} else if (structure == "unobserved") { 
  f1 <- "y ~ x2 + x3 + x4 + x5"
} else {
  f1 <- "y ~ l1 + x2 + x3 + x4 + x5 "
}
f0 <- "y ~ 1"

title_txt <- list(n_iter = n_iter, N = N, tau = tau, tau_beta = tau_beta, tau_sd = tau_sd, p = p, K = K, type = type, seed = seed)
title_txt_string <- paste("n_iter =", title_txt$n_iter, ", N =", title_txt$N, "\ntau =", title_txt$tau, ", tau_beta =", title_txt$tau_beta,  ", tau_sd =", title_txt$tau_sd, ", p =", title_txt$p, ", K =", title_txt$K, ", type =", title_txt$type, ", seed =", title_txt$seed)
fname_string <- paste("n_iter", title_txt$n_iter, "N", title_txt$N, "tau", title_txt$tau, "tau_beta", title_txt$tau_beta, "tau_sd", title_txt$tau_sd, "p", title_txt$p, "K", title_txt$K, "seed", title_txt$seed)


# Print all arguments
cat(paste("n_iter =", title_txt$n_iter, ", N =", title_txt$N, ", tau =", title_txt$tau, ", tau_beta =", title_txt$tau_beta,  ", tau_sd =", title_txt$tau_sd, ", p =", title_txt$p, ", K =", title_txt$K, ", type =", title_txt$type, ", seed =", title_txt$seed, '\n'))
cat(paste("num dpop =", dpop_num))

if (run_all == 1) {
  iters_run <- 1:n_iter
  cat(paste("\nrunning iters 1-1000\n"))
} else if (run_all == 0) {
  folder <- paste("iters-run-multiple/", type, "/", structure,"/", fname_string, sep = "")
  if (!dir.exists(folder)) {
    dir.create(folder)
  } 
  iters_run <- readRDS(paste(folder, ".rds", sep = ""))
  cat(paste("\nrunning iters", paste(iters_run, collapse=",") ,"\n"))
}



if (dpop_num == 9) {
  results <- simulation_nrun_fnc(n_iter = n_iter, N = N, K = K,
                                structure = structure,
                                tau = tau,
                                tau_beta = tau_beta,
                                tau_sd = tau_sd,
                                p = p,
                                f1 = f1,
                                f0 = f0,
                                perfFUN = 'mse',
                                save_data = FALSE,
                                seed = seed,
                                iters_run = iters_run,
                                cores = cores)
if (run_all == 1) {
  write_message <- paste(Sys.time(), ": Writing to", paste("output-multiple/", type, "/", structure, "/", fname_string, ".csv\n", sep = ""))
  cat(write_message)
  write.csv(results, paste("output-multiple/", type, "/", structure,"/", fname_string, ".csv", sep = ""), row.names = FALSE)
}
} else if (dpop_num == 1) {
  results <- simulation_nrun_fnc_2pop(n_iter = n_iter, N = N, K = K,
                                structure = structure,
                                tau = tau,
                                tau_beta = tau_beta,
                                tau_sd = tau_sd,
                                p = p,
                                f1 = f1,
                                f0 = f0,
                                perfFUN = 'mse',
                                save_data = FALSE,
                                seed = seed,
                                iters_run = iters_run,
                                cores = cores)
if (run_all == 1) {
  write_message <- paste(Sys.time(), ": Writing to", paste("output-single/", type, "/", structure, "/", fname_string, ".csv\n", sep = ""))
  cat(write_message)
  write.csv(results, paste("output-single/", type, "/", structure,"/", fname_string, ".csv", sep = ""), row.names = FALSE)
}
}

end_timestamp <- Sys.time()
end_message <- paste(end_timestamp,": Simulation ended \n")
time_message <- paste("Total time:", as.character(end_timestamp-start_timestamp))
cat(end_message)
cat(time_message)
