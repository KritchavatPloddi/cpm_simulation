

# y_hat = predicted outcome
# y = actual outcome

eva_brier <- function(y_hat, y) {
  return(   mean((y - y_hat)^2  ))
}

eva_auc <- function(y_hat, y) {
  return( auc(roc(y, y_hat, quiet=TRUE)) )
}

eva_citl <- function(LP, y) {
  # calibration-in-the-large
  model <- glm(y ~ offset(LP), family=binomial(link = 'logit'))
  citl <- as.numeric(coef(model)[1])
  return(citl)
  
}

eva_cslope <- function(LP, y) {
  # Calibration slope
  model <- glm(y ~ LP,  family=binomial(link = 'logit'))
  cslope <- as.numeric(coef(model)[2])
  return(cslope)
}



ext_validation <- function(model, target) {
  pred_y <- predict(model, newdata=target, type='response')
  LP <- log(pred_y / (1 - pred_y))
  
  
  result_auc <- eva_auc(y_hat=pred_y, y=target$y)
  result_brier <- eva_brier(y_hat=pred_y, y=target$y)
  result_citl <- eva_citl(LP=LP, y=target$y)
  result_cslope <- eva_cslope(LP=LP, y=target$y)
  
  result <- list(auc=result_auc, brier=result_brier, citl=result_citl, cslope=result_cslope)
  
  return(result)
}

ext_validation_pred <- function(y, y_pred) {
  if (mean(is.na(y_pred)) == 1 ) {
    return(list(auc=NA, brier=NA, citl=NA, cslope=NA))
  } else {
    LP <- log(y_pred / (1 - y_pred))
    result_auc <- eva_auc(y_hat=y_pred, y=y)
    result_brier <- eva_brier(y_hat=y_pred, y=y)
    result_citl <- eva_citl(LP=LP, y=y)
    result_cslope <- eva_cslope(LP=LP, y=y)
    result <- list(auc=result_auc, brier=result_brier, citl=result_citl, cslope=result_cslope)
    return(result)
  }


  
  
  
  
  
}

boot_validation <- function(boot_n, data, f1) {
  ## Bootstrap validation
  
  f1 <- as.formula(f1)
  
  ### 1. Apparent performance
  fit_appr <- glm(f1, data=data, family=binomial)
  pref_appr <- unlist(ext_validation(model=fit_appr, target=data))
  
  ### 2. Bootstrap sample
  optimism_list <- lapply(1:10, function(i) {
    boot_index <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    data_boot <- data[boot_index,]
    fit_boot <- glm(f1, data=data_boot, family=binomial)
    
    ### Bootstrap performance
    pref_boot <- ext_validation(model=fit_boot, target=data_boot)

    ### Original performance
    pref_orig <- ext_validation(model=fit_boot, target=data)

    optimism <- unlist(pref_boot) - unlist(pref_orig)
    return(optimism)
  })
  optimism <- colMeans(do.call(rbind, optimism_list))
  optimism_adj_pref <- pref_appr - optimism
  
  return(optimism_adj_pref)
}


boot_validation <- function(boot_n, data, f1) {
  ## Bootstrap validation
  ## boot_n = 200, # bootstrap iteratons
  
  # develop model in bootstrap sample
  boot <- lapply(1:boot_n, function(i) {
      boot_index <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
      data_dev <- data[boot_index,]
      fit_boot <- glm(as.formula(f1), data=data_dev, family=binomial)
      # validate model in original sample
      boot_val <- unlist(ext_validation(model=fit_boot, target=data))
  })
  boot <- data.frame(do.call(rbind, boot))
  
  result <- list(auc=mean(boot$auc), brier=mean(boot$brier), citl=mean(boot$citl), cslope=mean(boot$cslope))
  return(result)
}