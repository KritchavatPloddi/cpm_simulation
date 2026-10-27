### Generate betas and intercept
source("~/simulation/function_generate_intercept.R")
K <- 10
N <- 1000
p <- 0.2

P <- 5
set.seed(1234)

OR_min <- 1.3
OR_max <- 1.5

beta <- runif(n = P, min = log(OR_min), max = log(OR_max))

# simulate Xs

# Generate x1, x2, x3, x4, x5
x1 <- rnorm(n = N * K, mean = 0, sd = 1)
x2 <- rnorm(n = N * K, mean = 0, sd = 1)
x3 <- rnorm(n = N * K, mean = 0, sd = 1)
x4 <- rnorm(n = N * K, mean = 0, sd = 1)
x5 <- rnorm(n = N * K, mean = 0, sd = 1)

X <- cbind(x1, x2, x3, x4, x5)

intercept <- generate_intercept(N_pop = N*K, prob = p, X %*% beta)

betas <- c(intercept, beta)

print("Intercept and beta")
print(betas)
print(exp(betas))


saveRDS(betas, "~/simulation/simulated_data/betas.rds")
print("Betas saved to ~/simulation/simulated_data/betas.rds")