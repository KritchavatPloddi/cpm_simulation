suppressPackageStartupMessages(library(parallel, quietly = TRUE))

parallel_sapply <- function(X, FUN, mc.cores, ...) {
  res <- mclapply(X, FUN, mc.cores = mc.cores, ...)
  simplify2array(res)
}

parallel_lapply <- function(X, FUN, mc.cores, ...) {
  res <- mclapply(X, FUN, mc.cores = mc.cores, ...)
  return(res)
}


get_arg <- function(arg_name, default) {
  value <- args[grep(paste0("--", arg_name, "="), args)]
  if (length(value) == 0) {
    return(default)
  } else {
    return(as.numeric(sub(paste0("--", arg_name, "="), "", value)))
  }
}

read_data <- function(structure='observed', iter=1, tau=0, tau_beta=0, tau_sd=0, verbose=TRUE) {
  path <- paste('~/scratch/simulated_data/c/',structure,'/iter ',iter,' N 1000 tau ',tau,' tau_beta ',tau_beta,' tau_sd ',tau_sd,' p 0.2 K 10 seed 68250340.rds', sep="")
  if (verbose) {
    print(path)  
  }
  return(readRDS(path))
}