library(dplyr)
library(tidyr)
library(pROC)
## Density Ratio Weighting

drw <- function(dpop, tpop, columns, f1) {
  dpop$group <- factor(1, levels = c(0, 1)) # 1 = source (development)
  tpop$group <- factor(0, levels = c(0, 1)) # 0 = target
  
  membership <- rbind(dpop, tpop)
  
  predictors <- setdiff(columns, 'y')
  f_drw <- as.formula(paste("group ~", paste(predictors, collapse = "+")))
  
  fit_drw <- glm(f_drw, data = membership, family = binomial)
  
  # Log odds
  logodds_source <- predict(fit_drw, newdata = dpop, type = 'link')
  dpop$inv_ps_source <- exp(-logodds_source)
  # fit DR weight LR
  fit <- glm(f1, data=dpop, weights=inv_ps_source, family=binomial)
  return(fit)
}
## Graph surgery estimator
gse <- function(N, pop_csv, gse_argv, type, dpop, tpop,save_image,tau,tau_beta) {
    python_gse <- paste("gse-", type, ".py", sep="")
    results <- system2(
      paste(Sys.getenv("HOME"),"/miniconda3/envs/causal/bin/python", sep=""), 
      paste(python_gse, gse_argv, save_image, tau, tau_beta), 
      input = pop_csv, 
      stdout = TRUE, 
      stderr = FALSE
      ) 
    return(read.csv(text=results))
    #return(list(predicted_gse = predicted_gse, sum_vars = res$sum_vars, f_gse = f_gse))
}

hsictest <- function(x, y) {
  ht <- dHSIC::dhsic.test(X = matrix(x, ncol = 1), Y = matrix(y, ncol = 1), method = "gamma", kernel = c("gaussian", "discrete"))
  return(ht$p.value)
}

gen_all_combinations <- function(x_cols) {
    all_combinations <- unlist(
    lapply(1:length(x_cols), function(r) {
      combn(x_cols, r, simplify = FALSE)
    }), 
    recursive = FALSE
  )
  # add empty set y ~ 1
  all_combinations <- c(list("1"), all_combinations)
  return(all_combinations)
}


## Traditional ICP (on residuals)
## Using Wilcoxon test
icp_res <- function(dpop, columns, delta, valid_split, method='wilcoxon') {
  x_cols <- setdiff(columns, 'y')
  

  # Train-valid split
  train_size <- floor((1 - valid_split) * nrow(dpop))
  train_index <- sample(seq_len(nrow(dpop)), size = train_size)
  dpop_train <- dpop[train_index, ]
  dpop_valid  <- dpop[-train_index, ]

  # All Possible Combinations 
  all_combinations <- gen_all_combinations(x_cols)

  results <- data.frame(
    f = character(),
    pval = numeric(),
    stringsAsFactors = FALSE
  )

    # test for each combination
  for (i in seq_along(all_combinations)) {
    current_x <- all_combinations[[i]]
    f <- paste("y ~", paste( current_x, collapse=" + "))
    # Y ~ X model
    fit <- glm(as.formula(f), data = dpop_train, family = binomial)
    residual <- dpop_valid$y - predict(fit, newdata = dpop_valid, type = 'response')

    if (method == 'wilcoxon') {  
        # ดึงรายชื่อ environment/cluster ทั้งหมดจากข้อมูล validation
        environments <- unique(dpop_valid$cluster)
        k <- length(environments)

        p_values <- c()
        # one vs all test
        for (e in environments) {
          group_e <- residual[dpop_valid$cluster == e]
          group_others <- residual[dpop_valid$cluster != e]
          

          w_res <- wilcox.test(x = group_e, y = group_others)
          p_values <- c(p_values, w_res$p.value)
        }
        # Bonferroni correction 
        
        pval <- min(min(p_values) * k, 1.0)
    } else if (method == 'hsic') {  
      
      pval <- hsictest(residual, dpop_valid$cluster)
      
    }
    # store p values
    results <- rbind(results, data.frame(
      f = paste("y ~", paste(current_x, collapse = " + ")),
      pval = pval
    ))
  }
  subsets <- results[results$pval > delta,]
  if (nrow(subsets) > 0) {
    #finding the best mse on dpop_valid
    subset_valid <- data.frame(
      f = character(),
      pval = numeric(),
      mse = numeric()
    )

    for (i in 1:nrow(subsets)) {
      f_subset <- subsets[i, "f"]
      pval <- subsets[i, "pval"]
      fit <- glm(as.formula(f_subset), data = dpop_train, family = binomial)
      pred <- predict(fit, newdata = dpop_valid, type = "response")
      mse <- mean((pred - dpop_valid$y) ^ 2)
      subset_valid <- rbind(subset_valid, data.frame(f = f_subset, pval = pval, mse = mse))
    }
    best_subsets <- subset_valid[order(subset_valid$mse), ]
    rownames(best_subsets) <- NULL
    # return the formula of the best subset
    return(best_subsets[1,"f"])
  } else {
    return("no accepted set")
  }
}



icp_ICP <- function(dpop, columns, delta, valid_split, test = "normal") {
  ## use InvariantCausalPrediction package, which only takes the whole data, not splitted train-test
  x_cols <- setdiff(columns, 'y')


  train_size <- floor((1 - valid_split) * nrow(dpop))
  train_index <- sample(seq_len(nrow(dpop)), size = train_size)


  dpop_train <- dpop[train_index, ]
  dpop_valid  <- dpop[-train_index, ]

  X_mat <- as.matrix(dpop_train[, x_cols])
  Y_vec <- factor(dpop_train$y)
  ExpInd_vec <- factor(factor(dpop_train$cluster)) # แปลง cluster เป็นตัวเลขระบุ Environment
  
  # ICP test
  fit_icp <- InvariantCausalPrediction::ICP(X = X_mat, Y = Y_vec, ExpInd = ExpInd_vec, alpha = delta, test = test, showAcceptedSets = FALSE, showCompletion = FALSE)

  accepted <- fit_icp$acceptedSets

  if (length(accepted) > 0) {
       
      
      accepted_valid <- data.frame(
        f = character(),
        mse = numeric(),
        stringsAsFactors = FALSE
      )

      for (i in seq_along(accepted)) {
        indices <- accepted[[i]]
        
        if (length(indices) == 0) {
          f_subset <- "y ~ 1"
        } else {
          best_vars <- x_cols[indices]
          f_subset <- paste("y ~", paste(best_vars, collapse = " + "))
        }
        
        # MSE
        fit <- glm(as.formula(f_subset), data = dpop_train, family = binomial)
        pred <- predict(fit, newdata = dpop_valid, type = "response")
        mse <- mean((pred - dpop_valid$y) ^ 2)
        
  
        accepted_valid <- rbind(accepted_valid, data.frame(
          f = f_subset, 
          mse = mse
        ))
      }
      
      best_subsets <- accepted_valid[order(accepted_valid$mse), ]
      rownames(best_subsets) <- NULL
      
      return(best_subsets[1, "f"])
      
  } else {
      return("no accepted set")
  }
}



## Invariant Target Prediction
## Using AUC and MSE
itp <- function(dpop, columns, delta, valid_split, method='auc') {
  x_cols <- setdiff(columns, 'y')

  train_size <- floor((1 - valid_split) * nrow(dpop))
  train_index <- sample(seq_len(nrow(dpop)), size = train_size)

  dpop_train <- dpop[train_index, ]
  dpop_valid  <- dpop[-train_index, ]

  all_combinations <- gen_all_combinations(x_cols)

  results <- data.frame(
    f = character(),
    pval = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(all_combinations)) {
    current_x <- all_combinations[[i]]
    f <- paste("y ~", paste( current_x, collapse=" + "))
    # Y ~ X model
    fit_x <- glm(as.formula(f), data = dpop_train, family = binomial)
    # Y ~ X + E model
    fit_xe <- glm(as.formula(paste(f, "as.factor(cluster)", sep = " + ")), data = dpop_train, family = binomial)
    if (method == 'auc') {  
      roc_x <- suppressMessages(roc(factor(dpop_valid$y), predict(fit_x, newdata = dpop_valid, type = 'response')))
      roc_xe <- suppressMessages(roc(factor(dpop_valid$y), predict(fit_xe, newdata = dpop_valid, type = 'response')))
      roc_x_xe_test <- roc.test(roc_x, roc_xe, method = "delong")
      pval <- roc_x_xe_test$p.value
    } else if (method == 'mse') {
      se_xe <- (predict(fit_xe, newdata = dpop_valid, type = 'response') - dpop_valid$y)^2
      se_x <- (predict(fit_x, newdata = dpop_valid, type = 'response') - dpop_valid$y)^2
      wilcox.res <- wilcox.test(x = se_xe, y = se_x, alternative = "less")
      pval <- wilcox.res$p.value
    } else if (method == 'lrt') {
      lrt_result <- anova(fit_x, fit_xe, test = "Chisq")
      pval <- lrt_result$`Pr(>Chi)`[2]
    }

    results <- rbind(results, data.frame(
      f = paste("y ~", paste(current_x, collapse = " + ")),
      pval = pval
    ))
  }
  subsets <- results[results$pval > delta,]
  if (nrow(subsets) > 0) {
    #finding the best mse on dpop_valid
    subset_valid <- data.frame(
      f = character(),
      pval = numeric(),
      mse = numeric()
    )

    for (i in 1:nrow(subsets)) {
      f_subset <- subsets[i, "f"]
      pval <- subsets[i, "pval"]
      fit <- glm(as.formula(f_subset), data = dpop_train, family = binomial)
      pred <- predict(fit, newdata = dpop_valid, type = "response")
      mse <- mean((pred - dpop_valid$y) ^ 2)
      subset_valid <- rbind(subset_valid, data.frame(f = f_subset, pval = pval, mse = mse))
    }
    best_subsets <- subset_valid[order(subset_valid$mse), ]
    rownames(best_subsets) <- NULL
    # return the formula of the best subset
    return(best_subsets[1,"f"])
  } else {
    return("no accepted set")
  }
}

