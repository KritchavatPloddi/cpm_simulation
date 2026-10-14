# this function generate intercept to control overall event rate

generate_intercept <- function(N_pop, # number of population
                               prob, #expected event rate) 
                               lp # linear predictor excluding intercept 
                              )
                               {
  y <- rbinom(N_pop, 1, prob=prob)
  X <- as.vector(lp)
  model <- glm(y ~ offset(X), family=binomial(link = 'logit'))
  intercept <- as.numeric(coef(model)[1])
  return(intercept)
}

