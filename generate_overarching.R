
set.seed(68250340)
# 5 Predictors
P <- 5

# 1000 population

N <- 1000

# 1000 iters
iters <- 1000
for (iter in 1:iters) {
    print(iter)
    X <- matrix(rnorm(n = P * N, mean = 0, sd = 1), ncol=5)
    saveRDS(X, paste("~/scratch/simulated_data/overarching/iter_", iter, ".rds", sep=""))
}