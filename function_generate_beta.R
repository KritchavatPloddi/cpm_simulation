# This function generates beta (b0, b1, .., b5) coefficients (excluding intercept) which will be used for subsequent cluster generation.
# returns a matrix with each column represents each cluster, and each row represents each beta coefficients
# no intercept will be generated at this point


# This function generates beta (b0, b1, .., b5) coefficients which will be used for subsequent cluster generation.
# returns a matrix with each column represents each cluster, and each row represents each beta coefficients

generate_beta <- function(P = 5, mean_beta = 1, sd_beta = 1) { 
  # P = number of predictors

  P1 <- P + 1 # number of predictors + 1 intercept
  B <- rnorm(n=P, mean = mean_beta, sd = sd_beta)
  
  return(B)
}
