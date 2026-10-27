suppressPackageStartupMessages(library(pROC, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(metamisc, quietly = TRUE))
#suppressPackageStartupMessages(library(parallel, quietly = TRUE))
#suppressPackageStartupMessages(library(lme4, quietly = TRUE))
suppressPackageStartupMessages(library(MASS, quietly = TRUE))
source("function_generate_beta.R")
source("function_generate_clusters.R")
source("function_ext_validation.R")
source("function_models.R")
source("function_utils.R")

iteration_fnc <- function(
                    iter,
                    n_iter,
                    N,
                    K,
                    structure,
                    tau,
                    tau_beta,
                    tau_sd,
                    p,
                    perfFUN='mse',
                    f0="y ~ 1",
                    f1="y ~ x1 + x2 + x3 + x4 + x5",
                    type='c',
                    target=1, # kth subpopulation as target population
                    development=-1, #kth subpopulation as development population, -1 for all others not in target
                    seed=68250340
) { 
  #  target=10; development=-1; iter=1; type='c'; tau=1; tau_beta=0; tau_sd=0; p=0.2; K=10;  N=1000; structure='unobserved';seed=68250340;f0="y ~ 1";f1="y ~  x2 + x3 + x4 + x5"
  if (development == -1) {
    development_type <- "multiple"
  } else if (development %in% 1:K) {
    development_type <- "single"
  }
  
  cat(sprintf(paste("\r", as.character(Sys.time()), " : Iteration: %d", sep=""), iter))
  flush.console()
  
  # assign columns
  if (structure=='observed') {
    columns <- c("y", "x1", "x2", "x3", "x4", "x5")
  } else if (structure=='unobserved') {
    columns <- c("y", "x2", "x3", "x4", "x5")
  } else if (structure=='proxy') {
    columns <- c("y", "l1", "x2", "x3", "x4", "x5")
  }
  
  # Load data
  iter_filename <- paste("iter", iter, "N", N, "tau", tau, "tau_beta", tau_beta, "tau_sd", tau_sd, "p", p, "K", K, "seed", seed)
  iter_fullpath <- paste("~/scratch/simulated_data/",type,"/",structure, "/", iter_filename, ".rds", sep="")
  print(iter_fullpath)
  clusters <- readRDS(iter_fullpath)
  # set target pop as the first ranked cluster (among 1st, 3rd, 5th cluster)
  print(paste(development_type,"Development population =",development,"Target population = ", target))

  if (development_type == "multiple") {
    tpop <- clusters[clusters$cluster == target, c(columns, "cluster")] # target pop
    dpop <- clusters[clusters$cluster != target , c(columns, "cluster")] # development pop
  } else if (development_type == "single") {
    tpop <- clusters[clusters$cluster == target, c(columns, "cluster")] # target pop
    dpop <- clusters[clusters$cluster == development , c(columns, "cluster")] # development pop
  }

  tpop_dev_index <- readRDS(paste("~/simulation/simulated_data/tpop_dev_index_N_", N,".rds", sep=""))
  tpop_dev <- tpop[tpop_dev_index, ]
  tpop_val <- tpop[-tpop_dev_index, ]
  clusters$iter <- iter
  clusters$cluster <- factor(clusters$cluster)

  #print("Saved")
  f1 <- as.formula(f1)
  f0 <- as.formula(f0)

  # Fit models 0.0001"
  if (structure=='observed') {

    if (development_type == "multiple") {
      # SIECV
      # (lambda = 0)
      iecv0 <- metapred(data = dpop, strata = "cluster", formula = f1, scope = f0, retest = TRUE, estFUN = "glm", family = binomial, perfFUN = "mse",genFUN = "rema.tau")
      # (lambda = 1)
      iecv1 <- metapred(data = dpop, strata = "cluster", formula = f1, scope = f0, retest = TRUE, estFUN = "glm", family = binomial, perfFUN = "mse",genFUN = "rema.beta")

      ## Invariant causal prediction
        # Invariant target prediction (MSE)
      f_itp1 <- itp(dpop, columns, delta=0.05, valid_split=0.2, method = "mse")
        # Chi-square likelihood ratio test
      f_itp2 <- itp(dpop, columns, delta=0.05, valid_split=0.2, method = "lrt")
        # Residual test (ttest)
      f_icp1 <- icp_ICP(dpop, columns, delta=0.05, valid_split=0.2,  test='normal')
        # Residual test (hsic)
      f_icp2 <- icp_res(dpop, columns, delta=0.05, valid_split=0.2,  method="hsic")

      if (f_itp1 != "no accepted set") {
        itp1 <- glm(as.formula(f_itp1), data=dpop[, columns], family='binomial')
      }
      if (f_itp2 != "no accepted set") {
        itp2 <- glm(as.formula(f_itp2), data=dpop[, columns], family='binomial')
      }
      if (f_icp1 != "no accepted set") {
        icp1 <- glm(as.formula(f_icp1), data=dpop[, columns], family='binomial')
      }
      if (f_icp2 != "no accepted set") {
        icp2 <- glm(as.formula(f_icp2), data=dpop[, columns], family='binomial')
      }

      f_iecv0 <- deparse(iecv0$formula)
      f_iecv1 <- deparse(iecv1$formula)
    }

    # Standard
    null_model <- glm(as.formula(f0), data=dpop[,columns], family='binomial')
    m1 <- glm(as.formula(f1), data=dpop[, columns], family='binomial')
  
    if (development_type == 'single' & tau_beta > 0) {
      # fit logistic regression with external knowledge for predictor outcome shift
      b1 <-readRDS(paste("~/simulation/simulated_data/b1_tau_beta",tau_beta,".rds", sep="_"))
      tau_extknowledge <- 0.05
      b1_extknowledge <- rnorm(1, b1[target], tau_extknowledge)
      offset_x1 <- b1_extknowledge * dpop$x1
      m2 <- glm(y ~ x2 + x3 + x4 + x5, data = dpop, family = binomial, offset = offset_x1)
    }

    ### Graph Surgery Estimator
    # Pop for GSE
    pop <- rbind(mutate(dpop, pop='dpop'), mutate(tpop, pop='tpop'))
    pop_csv <- paste(capture.output(write.csv(pop,row.names = FALSE)), collapse = "\n")

    # GSE
    if (iter == 1) {
      gse1 <- gse(N = N, pop_csv = pop_csv, gse_argv = "observed GSE1", type = type, save_image = "T", tau = tau, tau_beta = tau_beta)
      gse2 <- gse(N = N, pop_csv = pop_csv, gse_argv = "observed GSE2", type = type, save_image = "T", tau = tau, tau_beta = tau_beta)
    } else {
      gse1 <- gse(N = N, pop_csv = pop_csv, gse_argv = "observed GSE1", type = type, save_image = "F", tau = tau, tau_beta = tau_beta)
      gse2 <- gse(N = N, pop_csv = pop_csv, gse_argv = "observed GSE2", type = type, save_image = "F", tau = tau, tau_beta = tau_beta)
    }
    f_gse1 <- gse1$expr[1]
    f_gse2 <- gse2$expr[1]
    # Density ratio
    
    dr <- drw(dpop, tpop, columns, f1)
    f_dr <- deparse(dr$formula)

    # Baseline (Bootstrap in tpop)
    results_base <- boot_validation(boot_n = 200, data = tpop, f1 = f1)
    #m_base <- glm(f1, data=tpop_dev[, columns], family='binomial')

    # get formula
    
    f_m1 <- deparse(m1$formula)

    
  } else if (structure=='unobserved') {
    if (development_type == "multiple") {
      # SIECV
      # (lambda = 0)
      iecv0 <- metapred(data = dpop, strata = "cluster", formula = f1, scope = f0, retest = TRUE, estFUN = "glm", family = binomial, perfFUN = "mse",genFUN = "rema.tau")
      # (lambda = 1)
      iecv1 <- metapred(data = dpop, strata = "cluster", formula = f1, scope = f0, retest = TRUE, estFUN = "glm", family = binomial, perfFUN = "mse",genFUN = "rema.beta")

      ## Invariant causal prediction
        # Invariant target prediction (MSE)
      f_itp1 <- itp(dpop, columns, delta=0.05, valid_split=0.2, method = "mse")
        # Chi-square likelihood ratio test
      f_itp2 <- itp(dpop, columns, delta=0.05, valid_split=0.2, method = "lrt")
        # Residual test (ttest)
      f_icp1 <- icp_ICP(dpop, columns, delta=0.05, valid_split=0.2,  test='normal')
        # Residual test (hsic)
      f_icp2 <- icp_res(dpop, columns, delta=0.05, valid_split=0.2,  method="hsic")

      if (f_itp1 != "no accepted set") {
        itp1 <- glm(as.formula(f_itp1), data=dpop[, columns], family='binomial')
      }
      if (f_itp2 != "no accepted set") {
        itp2 <- glm(as.formula(f_itp2), data=dpop[, columns], family='binomial')
      }
      if (f_icp1 != "no accepted set") {
        icp1 <- glm(as.formula(f_icp1), data=dpop[, columns], family='binomial')
      }
      if (f_icp2 != "no accepted set") {
        icp2 <- glm(as.formula(f_icp2), data=dpop[, columns], family='binomial')
      }
    
      f_iecv0 <- deparse(iecv0$formula)
      f_iecv1 <- deparse(iecv1$formula)

    }
    # Standard
    
    #full_model <- glm(as.formula(f1), data=dpop[, columns], family='binomial')
    null_model <- glm(as.formula(f0), data=dpop[, columns], family='binomial')
    #m1 <- stepAIC(full_model, direction='both', scope=list(lower=null_model, upper=full_model), trace=FALSE)
    m1 <- glm(as.formula(f1), data=dpop[, columns], family='binomial')

    ### Graph Surgery Estimator
    # Pop for GSE
    pop <- rbind(mutate(dpop, pop='dpop'), mutate(tpop, pop='tpop'))
    pop_csv <- paste(capture.output(write.csv(pop,row.names = FALSE)), collapse = "\n")

    # GSE
    if (iter == 1) {
      gse1 <- gse(N = N, pop_csv = pop_csv, gse_argv = "unobserved GSE1", type = type, save_image = "T", tau = tau, tau_beta = tau_beta)
      gse2 <- gse(N = N, pop_csv = pop_csv, gse_argv = "unobserved GSE2", type = type, save_image = "T", tau = tau, tau_beta = tau_beta)
    } else {
      gse1 <- gse(N = N, pop_csv = pop_csv, gse_argv = "unobserved GSE1", type = type, save_image = "F", tau = tau, tau_beta = tau_beta)
      gse2 <- gse(N = N, pop_csv = pop_csv, gse_argv = "unobserved GSE2", type = type, save_image = "F", tau = tau, tau_beta = tau_beta)
    }
    f_gse1 <- gse1$expr[1]
    f_gse2 <- gse2$expr[1]
    # Density ratio
    dr <- drw(dpop, tpop, columns, f1)
    f_dr <- deparse(dr$formula)

    # Baseline (Bootstrap in tpop)
    results_base <- boot_validation(boot_n = 200, data = tpop, f1 = f1)
    #m_base <- glm(f1, data=tpop_dev[, columns], family='binomial')

    # get formula
    
    f_m1 <- deparse(m1$formula)



  } else if (structure=='proxy') {

    if (development_type == "multiple") {
      # SIECV
      # (lambda = 0)
      iecv0 <- metapred(data = dpop, strata = "cluster", formula = f1, scope = f0, retest = TRUE, estFUN = "glm", family = binomial, perfFUN = "mse",genFUN = "rema.tau")
      # (lambda = 1)
      iecv1 <- metapred(data = dpop, strata = "cluster", formula = f1, scope = f0, retest = TRUE, estFUN = "glm", family = binomial, perfFUN = "mse",genFUN = "rema.beta")

      ## Invariant causal prediction
        # Invariant target prediction (MSE)
      f_itp1 <- itp(dpop, columns, delta=0.05, valid_split=0.2, method = "mse")
        # Chi-square likelihood ratio test
      f_itp2 <- itp(dpop, columns, delta=0.05, valid_split=0.2, method = "lrt")
        # Residual test (ttest)
      f_icp1 <- icp_ICP(dpop, columns, delta=0.05, valid_split=0.2,  test='normal')
        # Residual test (hsic)
      f_icp2 <- icp_res(dpop, columns, delta=0.05, valid_split=0.2,  method="hsic")

      if (f_itp1 != "no accepted set") {
        itp1 <- glm(as.formula(f_itp1), data=dpop[, columns], family='binomial')
      }
      if (f_itp2 != "no accepted set") {
        itp2 <- glm(as.formula(f_itp2), data=dpop[, columns], family='binomial')
      }
      if (f_icp1 != "no accepted set") {
        icp1 <- glm(as.formula(f_icp1), data=dpop[, columns], family='binomial')
      }
      if (f_icp2 != "no accepted set") {
        icp2 <- glm(as.formula(f_icp2), data=dpop[, columns], family='binomial')
      }
    
      f_iecv0 <- deparse(iecv0$formula)
      f_iecv1 <- deparse(iecv1$formula)
    }
    # Standard
    #full_model <- glm(as.formula(f1), data=dpop[, columns], family='binomial')
    null_model <- glm(as.formula(f0), data=dpop[, columns], family='binomial')
    #m1 <- stepAIC(full_model, direction='both', scope=list(lower=null_model, upper=full_model), trace=FALSE)
    m1 <- glm(as.formula(f1), data=dpop[, columns], family='binomial')


    ### Graph Surgery Estimator
    # Pop for GSE
    pop <- rbind(mutate(dpop, pop='dpop'), mutate(tpop, pop='tpop'))
    pop_csv <- paste(capture.output(write.csv(pop,row.names = FALSE)), collapse = "\n")

    # GSE
    if (iter == 1) {
      gse1 <- gse(N = N, pop_csv = pop_csv, gse_argv = "proxy GSE1", type = type, save_image = "T", tau = tau, tau_beta = tau_beta)
      gse2 <- gse(N = N, pop_csv = pop_csv, gse_argv = "proxy GSE2", type = type, save_image = "T", tau = tau, tau_beta = tau_beta)
    } else {
      gse1 <- gse(N = N, pop_csv = pop_csv, gse_argv = "proxy GSE1", type = type, save_image = "F", tau = tau, tau_beta = tau_beta)
      gse2 <- gse(N = N, pop_csv = pop_csv, gse_argv = "proxy GSE2", type = type, save_image = "F", tau = tau, tau_beta = tau_beta)
    }
    f_gse1 <- gse1$expr[1]
    f_gse2 <- gse2$expr[1]
    # Density ratio
    dr <- drw(dpop, tpop, columns, f1)
    f_dr <- deparse(dr$formula)

    # Baseline (Bootstrap in tpop)
    results_base <- boot_validation(boot_n = 200, data = tpop, f1 = f1)
    #m_base <- glm(f1, data=tpop_dev[, columns], family='binomial')
    
    # get formula
    
    f_m1 <- deparse(m1$formula)
  }

  # external validation in tpop
  if (development_type == "multiple") {
    results_iecv0 <- ext_validation(model=iecv0, target=tpop)
    results_iecv1 <- ext_validation(model=iecv1, target=tpop)

    if (f_itp1 != "no accepted set") {
      results_itp1 <- ext_validation(model=itp1, target=tpop)
    } else {
      results_itp1 <- list(auc=NA, brier=NA, citl=NA, cslope=NA)
    }
    if (f_itp2 != "no accepted set") {
      results_itp2 <- ext_validation(model=itp2, target=tpop)
    } else {
      results_itp2 <- list(auc=NA, brier=NA, citl=NA, cslope=NA)
    }
    if (f_icp1 != "no accepted set") {
      results_icp1 <- ext_validation(model=icp1, target=tpop)
    } else {
      results_icp1 <- list(auc=NA, brier=NA, citl=NA, cslope=NA)
    }
    if (f_icp2 != "no accepted set") {
      results_icp2 <- ext_validation(model=icp2, target=tpop)
    } else {
      results_icp2 <- list(auc=NA, brier=NA, citl=NA, cslope=NA)
    }
  }  

  results_null <- ext_validation(model=null_model, target=tpop)
  results_m1 <- ext_validation(model=m1, target=tpop)

  # Graph surgery
  results_gse1 <- ext_validation_pred(y=tpop$y, y_pred=gse1$predicted_gse)
  results_gse2 <- ext_validation_pred(y=tpop$y, y_pred=gse2$predicted_gse)

  # DRW
  results_dr <- ext_validation(model=dr, target=tpop)

  if (development_type == 'single' & structure=='observed' & tau_beta > 0) {
    results_m2 <- ext_validation(model=m2, target=tpop)
  }

# save predicted files
  paramname <- paste("iter", iter,"N", N, "tau", tau, "tau_beta", tau_beta, "tau_sd", tau_sd,  "p", p, "target", target)

#  data.frame(
#    iter = iter,
#    observed = tpop$y,
#    predicted_m1 = predict(m1, newdata=tpop, type='response'),
#    predicted_prunning = predict(prunning, newdata=tpop, type='response'),
#    predicted_iecv0 = predict(iecv0, newdata=tpop, type='response'),
#    predicted_iecv1 = predict(iecv1, newdata=tpop, type='response'),
#    predicted_dr = predict(dr, newdata=tpop, type='response'),
#    predicted_gse1 = gse1$predicted_gse,
#    predicted_gse2 = gse2$predicted_gse
#    #predicted_gse1 = predict(gse1, newdata=tpop, type='response')
#    #predicted_gse2 = predict(gse2, newdata=tpop, type='response')
#  ) %>%
#  saveRDS(paste("~/scratch/simulated_predicted/", type, "/", structure,"/", paramname,".rds", sep=""))
#
  # collect results into a list
  if (development_type == "multiple") {
    result <- data.frame(
      iter = iter,
      data_path = iter_fullpath,
      auc_null = results_null$auc,
      brier_null = results_null$brier,
      citl_null = results_null$citl,
      cslope_null = results_null$cslope,
      auc_m1 = results_m1$auc,
      brier_m1 = results_m1$brier,
      citl_m1 = results_m1$citl,
      cslope_m1 = results_m1$cslope,
      # IECV
      auc_iecv0 = results_iecv0$auc,
      brier_iecv0 = results_iecv0$brier,
      citl_iecv0 = results_iecv0$citl,
      cslope_iecv0 = results_iecv0$cslope,

      auc_iecv1 = results_iecv1$auc,
      brier_iecv1 = results_iecv1$brier,
      citl_iecv1 = results_iecv1$citl,
      cslope_iecv1 = results_iecv1$cslope,

      # GSE
      auc_gse1 = results_gse1$auc,
      brier_gse1 = results_gse1$brier,
      citl_gse1 = results_gse1$citl,
      cslope_gse1 = results_gse1$cslope,  
      auc_gse2 = results_gse2$auc,
      brier_gse2 = results_gse2$brier,
      citl_gse2 = results_gse2$citl,
      cslope_gse2 = results_gse2$cslope,


      auc_base = results_base$auc,
      brier_base = results_base$brier,
      citl_base = results_base$citl,
      cslope_base = results_base$cslope,

      # ICP
      
      auc_itp1 = results_itp1$auc,
      brier_itp1 = results_itp1$brier,
      citl_itp1 = results_itp1$citl,
      cslope_itp1 = results_itp1$cslope,

      auc_itp2 = results_itp2$auc,
      brier_itp2 = results_itp2$brier,
      citl_itp2 = results_itp2$citl,
      cslope_itp2 = results_itp2$cslope,

      auc_icp1 = results_icp1$auc,
      brier_icp1 = results_icp1$brier,
      citl_icp1 = results_icp1$citl,
      cslope_icp1 = results_icp1$cslope,

      auc_icp2 = results_icp2$auc,
      brier_icp2 = results_icp2$brier,
      citl_icp2 = results_icp2$citl,
      cslope_icp2 = results_icp2$cslope,

      # DRW
      auc_dr = results_dr$auc,
      brier_dr = results_dr$brier,
      citl_dr = results_dr$citl,
      cslope_dr = results_dr$cslope,
      

      f_m1 = as.character(f_m1), 
      f_iecv0 = f_iecv0, 
      f_iecv1 = f_iecv1, 
      f_itp1 = f_itp1,
      f_itp2 = f_itp2,
      f_icp1 = f_icp1,
      f_icp2 = f_icp2,
      f_gse1 = f_gse1,
      f_gse2 = f_gse2,
      f_dr = f_dr
    )
  } else if (development_type == 'single') {
    result <- data.frame(
      iter = iter,
      data_path = iter_fullpath,
      auc_null = results_null$auc,
      brier_null = results_null$brier,
      citl_null = results_null$citl,
      cslope_null = results_null$cslope,
      auc_m1 = results_m1$auc,
      brier_m1 = results_m1$brier,
      citl_m1 = results_m1$citl,
      cslope_m1 = results_m1$cslope,
      # GSE
      auc_gse1 = results_gse1$auc,
      brier_gse1 = results_gse1$brier,
      citl_gse1 = results_gse1$citl,
      cslope_gse1 = results_gse1$cslope,  
      auc_gse2 = results_gse2$auc,
      brier_gse2 = results_gse2$brier,
      citl_gse2 = results_gse2$citl,
      cslope_gse2 = results_gse2$cslope,

      auc_base = results_base$auc,
      brier_base = results_base$brier,
      citl_base = results_base$citl,
      cslope_base = results_base$cslope,

      # DRW
      auc_dr = results_dr$auc,
      brier_dr = results_dr$brier,
      citl_dr = results_dr$citl,
      cslope_dr = results_dr$cslope,
  
      f_m1 = as.character(f_m1), 
      f_gse1 = f_gse1,
      f_gse2 = f_gse2,
      f_dr = f_dr
    )
  }


  if (development_type == 'single' & structure=='observed' & tau_beta > 0) {
    result <- cbind(
      result,
      auc_m2 = results_m2$auc,
      brier_m2 = results_m2$brier,
      citl_m2 = results_m2$citl,
      cslope_m2 = results_m2$cslope
    )

  }
#fname_string <- paste("n_iter", n_iter, "N", N, "tau", tau, "tau_beta", tau_beta, "tau_sd", tau_sd, "p", p, "K", K, "seed", seed)
#folder <- paste("output-multiple/", type, "/", structure,"/", fname_string, sep = "")
#if (!dir.exists(folder)) {
#  dir.create(folder)
#}
#saveRDS(result, paste(folder, "/", iter,".rds", sep=""))

print(paste("\nFinish Iter", iter,"\n"))
  #saveRDS(result, paste("debug/iter_", iter, ".rds", sep=""))
#  fin_msg <- paste(Sys.time(),"Finish iter", iter)
#  writeRDS(fin_msg, paste("logs/",type,"/",structure, "/", iter_filename, ".rds", sep=""))
  return(result)
}

simulation_nrun_fnc <- function(n_iter,
                                N,
                                K,
                                structure,
                                tau,
                                tau_beta,
                                tau_sd,
                                p,
                                perfFUN='mse',
                                f0,
                                f1,
                                save_data=FALSE,
                                msg='',
                                type='c',
                                seed=1,
                                target=1,
                                development=-1,
                                iters_run=1:n_iter,
                                cores=cores) {
  #perfFUN = 'mse'
  # n_iter = number of iterations to repeat simulations over
  # N = number of observations each run
  # K = number of populations (heterogeneity)
  # structure = causal structure (observed cause, unobserved cause, proxy)
  # tau = tau of intercept
  # tau_beta = tau of slope
  # p = event_Rate
  # perfFUN = performance metric for IECV held out evaluation
  # f0 = fitting formula (model specification)
  # f1 = scope formula (initial stepwise specification)
  # save_data = whether the simulated data will be saved on disk

  ## Repeat the simulation across iter number of iterations


# register randomed simulations id

timestamp <- as.character(as.integer(Sys.time()))
print(paste(as.character(Sys.time()), "Timestamp:", timestamp ))

## set seed
set.seed(seed)
f1 <- as.formula(f1)
f0 <- as.formula(f0)

#results_iter <- parallel_sapply(1:n_iter, iteration_fnc, mc.cores = cores)
#results_iter <- parallel_lapply(1:n_iter, iteration_fnc, mc.cores = cores)

results_iter <- mclapply(iters_run, iteration_fnc,  
                                    mc.cores = cores,
                                    n_iter = n_iter,
                                    N = N,
                                    K = K,
                                    structure = structure,
                                    tau = tau,
                                    tau_beta = tau_beta,
                                    tau_sd = tau_sd,
                                    p = p,
                                    perfFUN='mse',
                                    f0 = f0,
                                    f1 = f1,
                                    type='c',
                                    target = target,
                                    development = development,
                                    seed = seed)

    
print("Finish parallel_sapply")
#saveRDS(results_iter, "~/simulation/debug.rds")

# combine results
perfResults <- do.call(rbind, results_iter)

  #perfResults <- data.frame(t(data.frame(results_iter)))
  #rownames(perfResults) <- 1:length(perfResults$auc_m1)

  return(perfResults)
  print("Finished")
  print(msg)
  

  
  
}



## siecv (h) 
#iecv1 <- metapred(data = dpop, strata = "cluster", scope = f0,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "mse", genFUN = "rema.tau")
#iecv1.1 <- suppressWarnings(suppressMessages(metapred(data = dpop, strata = "cluster", scope = y ~ x1,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "auc", genFUN = "rema.tau")))
#iecv1.2 <- suppressWarnings(suppressMessages(metapred(data = dpop, strata = "cluster", scope = y ~ x2,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "auc", genFUN = "rema.tau")))
#iecv1.3 <- suppressWarnings(suppressMessages(metapred(data = dpop, strata = "cluster", scope = y ~ x3,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "auc", genFUN = "rema.tau")))
#iecv1.4 <- suppressWarnings(suppressMessages(metapred(data = dpop, strata = "cluster", scope = y ~ x4,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "auc", genFUN = "rema.tau")))
#iecv1.5 <- suppressWarnings(suppressMessages(metapred(data = dpop, strata = "cluster", scope = y ~ x5,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "auc", genFUN = "rema.tau")))
#gen_values <- c(gen(iecv1.1), gen(iecv1.2), gen(iecv1.3), gen(iecv1.4),gen(iecv1.5))
#iecv1 <- list(iecv1.1, iecv1.2, iecv1.3, iecv1.4, iecv1.5)[[which.min(gen_values)]]
#iecv3 <- suppressWarnings(suppressMessages(metapred(data = dpop, strata = "cluster", scope = f0,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "cal.int", genFUN = "rema.tau")))
## siecv (h) 
#iecv4.1 <- metapred(data = dpop, strata = "cluster", scope = y ~ x1,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "cal.slope", genFUN = "rema.tau")
#iecv4.2 <- metapred(data = dpop, strata = "cluster", scope = y ~ x2,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "cal.slope", genFUN = "rema.tau")
#iecv4.3 <- metapred(data = dpop, strata = "cluster", scope = y ~ x3,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "cal.slope", genFUN = "rema.tau")
#iecv4.4 <- metapred(data = dpop, strata = "cluster", scope = y ~ x4,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "cal.slope", genFUN = "rema.tau")
#iecv4.5 <- metapred(data = dpop, strata = "cluster", scope = y ~ x5,formula = f1,retest = TRUE,estFUN = "glm",family = binomial,perfFUN = "cal.slope", genFUN = "rema.tau")
#gen_values <- c(gen(iecv4.1), gen(iecv4.2), gen(iecv4.3), gen(iecv4.4),gen(iecv4.5))
#iecv4 <- list(iecv4.1, iecv4.2, iecv4.3, iecv4.4, iecv4.5)[[which.min(gen_values)]]