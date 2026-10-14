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