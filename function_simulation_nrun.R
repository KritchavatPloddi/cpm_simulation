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
                    seed=1
) { 
  #iter <- 8
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
  tpop <- clusters[clusters$cluster == 1, c(columns, "cluster")] # target pop
  dpop <- clusters[clusters$cluster != 1 , c(columns, "cluster")] # development pop
  

  tpop_dev_index <- readRDS(paste("simulated_data/tpop_dev_index_N_", N,".rds", sep=""))
  tpop_dev <- tpop[tpop_dev_index, ]
  tpop_val <- tpop[-tpop_dev_index, ]
  clusters$iter <- iter
  clusters$cluster <- factor(clusters$cluster)

  #print("Saved")
  f1 <- as.formula(f1)
  f0 <- as.formula(f0)


  # Fit models
  if (structure=='observed') {
    
    # Standard
    full_model <- glm(as.formula(f1), data=dpop[,columns], family='binomial')
    null_model <- glm(as.formula(f0), data=dpop[,columns], family='binomial')
    m1 <- stepAIC(full_model, direction='both', scope=list(lower=null_model, upper=full_model), trace=FALSE)
  

    # SIECV
    # Base SIECV model, not used
    
    iecv0 <- metapred(data = dpop,
                      strata = "cluster",
                      formula = f1, 
                      scope = f1,
                      retest = TRUE,
                      estFUN = "glm",
                      family = binomial,
                      perfFUN = "mse",
                      genFUN = "abs.mean")
                      
    iecv1 <- update(iecv0, scope = f0)
    iecv2 <- update(iecv1, estFUN = "logistfirth", genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 1 )
    iecv3 <- update(iecv2, genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 1/2 )
    iecv4 <- update(iecv3, genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 0 )


    # Baseline (Bootstrap in tpop)
    #results_base <- boot_validation(boot_n = 200, data = tpop, f1 = f1)
    m_base <- glm(f1, data=tpop_dev[, columns], family='binomial')

    # Prunning
    dpop_csv <- paste(capture.output(write.csv(dpop, row.names = FALSE)), collapse = "\n")
    f_prunning <- system2(paste(Sys.getenv("HOME"),"/miniconda3/envs/ctl/bin/python", sep=""), 
                          paste("causal_transfer_learning/prunning.py","x1+x2+x3+x4+x5 y"), 
                          input=dpop_csv, stdout = TRUE, stderr = NULL)
    prunning <- glm(f_prunning, data=dpop[, columns], family='binomial')
    
    
    ### Graph Surgery Estimator
    # Pop for GSE
    pop <- rbind(mutate(dpop, pop='dpop'), mutate(tpop, pop='tpop'))
    pop_csv <- paste(capture.output(write.csv(pop,row.names = FALSE)), collapse = "\n")

    # GSE
    gse1 <- gse(pop_csv = pop_csv, gse_argv = "observed GSE1", type = type)
    gse2 <- gse(pop_csv = pop_csv, gse_argv = "observed GSE2", type = type)
    f_gse1 <- gse1$expr[1]
    f_gse2 <- gse2$expr[1]
    # Density ratio
    dr <- drw(dpop, tpop, columns, f1)
    f_dr <- deparse(dr$formula)

    # get formula
    
    f_m1 <- deparse(m1$formula)
    f_iecv1 <- deparse(iecv1$formula)
    f_iecv2 <- deparse(iecv2$formula)
    f_iecv3 <- deparse(iecv3$formula)
    f_iecv4 <- deparse(iecv4$formula)
    #f_base
    
    print(paste("Formula m1:", f_m1, sep=" "))
    print(paste("Formula iecv1:", f_iecv1, sep=" "))
    print(paste("Formula iecv2:", f_iecv2, sep=" "))
    print(paste("Formula iecv3:", f_iecv3, sep=" "))
    print(paste("Formula iecv4:", f_iecv4, sep=" "))
    #print(paste("Formula base:", f_base, sep=" "))
    print(paste("Formula prunning:", f_prunning, sep=" "))
    #print(paste("Formula GSE1:", gse1$intervention_conditioning_set[1], sep=" "))
    #print(paste("Formula GSE2:", gse2$intervention_conditioning_set[1], sep=" "))
    print(paste("Formula GSE1:",f_gse1, sep=" "))
    print(paste("Formula GSE2:",f_gse2, sep=" "))
    print(paste("Formula DRW:", f_dr, sep=" "))
  } else if (structure=='unobserved') {
    # Standard
    full_model <- glm(as.formula(f1), data=dpop[, columns], family='binomial')
    null_model <- glm(as.formula(f0), data=dpop[, columns], family='binomial')
    m1 <- stepAIC(full_model, direction='both', scope=list(lower=null_model, upper=full_model), trace=FALSE)

    # SIECV
    # Base SIECV model, not used
    
    iecv0 <- metapred(data = dpop,
                      strata = "cluster",
                      formula = f1, 
                      scope = f1,
                      retest = TRUE,
                      estFUN = "glm",
                      family = binomial,
                      perfFUN = "mse",
                      genFUN = "abs.mean")
                      
    iecv1 <- update(iecv0, scope = f0)
    iecv2 <- update(iecv1, estFUN = "logistfirth", genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 1 )
    iecv3 <- update(iecv2, genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 1/2 )
    iecv4 <- update(iecv3, genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 0 )


    # Baseline (Bootstrap in tpop)
    #results_base <- boot_validation(boot_n = 200, data = tpop, f1 = f1)
    m_base <- glm(f1, data=tpop_dev[, columns], family='binomial')

    # Prunning
    dpop_csv <- paste(capture.output(write.csv(dpop, row.names = FALSE)), collapse = "\n")
    f_prunning <- system2(paste(Sys.getenv("HOME"),"/miniconda3/envs/ctl/bin/python", sep=""), 
                          paste("causal_transfer_learning/prunning.py","x2+x3+x4+x5 y"), 
                          input=dpop_csv, stdout = TRUE, stderr = NULL)
    prunning <- glm(f_prunning, data=dpop[, columns], family='binomial')
    
    
    ### Graph Surgery Estimator
    # Pop for GSE
    pop <- rbind(mutate(dpop, pop='dpop'), mutate(tpop, pop='tpop'))
    pop_csv <- paste(capture.output(write.csv(pop,row.names = FALSE)), collapse = "\n")

    # GSE
    gse1 <- gse(pop_csv = pop_csv, gse_argv = "unobserved GSE1", type = type)
    gse2 <- gse(pop_csv = pop_csv, gse_argv = "unobserved GSE2", type = type)
    f_gse1 <- gse1$expr[1]
    f_gse2 <- gse2$expr[1]
    # Density ratio
    dr <- drw(dpop, tpop, columns, f1)
    f_dr <- deparse(dr$formula)

    # get formula
    
    f_m1 <- deparse(m1$formula)
    f_iecv1 <- deparse(iecv1$formula)
    f_iecv2 <- deparse(iecv2$formula)
    f_iecv3 <- deparse(iecv3$formula)
    f_iecv4 <- deparse(iecv4$formula)
    #f_base
    
    print(paste("Formula m1:", f_m1, sep=" "))
    print(paste("Formula iecv1:", f_iecv1, sep=" "))
    print(paste("Formula iecv2:", f_iecv2, sep=" "))
    print(paste("Formula iecv3:", f_iecv3, sep=" "))
    print(paste("Formula iecv4:", f_iecv4, sep=" "))
    #print(paste("Formula base:", f_base, sep=" "))
    print(paste("Formula prunning:", f_prunning, sep=" "))
    #print(paste("Formula GSE1:", gse1$intervention_conditioning_set[1], sep=" "))
    #print(paste("Formula GSE2:", gse2$intervention_conditioning_set[1], sep=" "))
    print(paste("Formula GSE1:",f_gse1, sep=" "))
    print(paste("Formula GSE2:",f_gse2, sep=" "))
    print(paste("Formula DRW:", f_dr, sep=" "))
  } else if (structure=='proxy') {
    # Standard
    full_model <- glm(as.formula(f1), data=dpop[, columns], family='binomial')
    null_model <- glm(as.formula(f0), data=dpop[, columns], family='binomial')
    m1 <- stepAIC(full_model, direction='both', scope=list(lower=null_model, upper=full_model), trace=FALSE)

    # SIECV
    # Base SIECV model, not used
    
    iecv0 <- metapred(data = dpop,
                      strata = "cluster",
                      formula = f1, 
                      scope = f1,
                      retest = TRUE,
                      estFUN = "glm",
                      family = binomial,
                      perfFUN = "mse",
                      genFUN = "abs.mean")
                      
    iecv1 <- update(iecv0, scope = f0)
    iecv2 <- update(iecv1, estFUN = "logistfirth", genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 1 )
    iecv3 <- update(iecv2, genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 1/2 )
    iecv4 <- update(iecv3, genFUN = list("rema", "rema.beta", "rema.tau"), lambda = 0 )


    # Baseline (Bootstrap in tpop)
    #results_base <- boot_validation(boot_n = 200, data = tpop, f1 = f1)
    m_base <- glm(f1, data=tpop_dev[, columns], family='binomial')

    # Prunning
    dpop_csv <- paste(capture.output(write.csv(dpop, row.names = FALSE)), collapse = "\n")
    f_prunning <- system2(paste(Sys.getenv("HOME"),"/miniconda3/envs/ctl/bin/python", sep=""), 
                          paste("causal_transfer_learning/prunning.py","l1+x2+x3+x4+x5 y"), 
                          input=dpop_csv, stdout = TRUE, stderr = NULL)
    prunning <- glm(f_prunning, data=dpop[, columns], family='binomial')
    
    
    ### Graph Surgery Estimator
    # Pop for GSE
    pop <- rbind(mutate(dpop, pop='dpop'), mutate(tpop, pop='tpop'))
    pop_csv <- paste(capture.output(write.csv(pop,row.names = FALSE)), collapse = "\n")

    # GSE
    gse1 <- gse(pop_csv = pop_csv, gse_argv = "proxy GSE1", type = type)
    gse2 <- gse(pop_csv = pop_csv, gse_argv = "proxy GSE2", type = type)
    f_gse1 <- gse1$expr[1]
    f_gse2 <- gse2$expr[1]
    # Density ratio
    dr <- drw(dpop, tpop, columns, f1)
    f_dr <- deparse(dr$formula)

    # get formula
    
    f_m1 <- deparse(m1$formula)
    f_iecv1 <- deparse(iecv1$formula)
    f_iecv2 <- deparse(iecv2$formula)
    f_iecv3 <- deparse(iecv3$formula)
    f_iecv4 <- deparse(iecv4$formula)
    #f_base
    
    print(paste("Formula m1:", f_m1, sep=" "))
    print(paste("Formula iecv1:", f_iecv1, sep=" "))
    print(paste("Formula iecv2:", f_iecv2, sep=" "))
    print(paste("Formula iecv3:", f_iecv3, sep=" "))
    print(paste("Formula iecv4:", f_iecv4, sep=" "))
    #print(paste("Formula base:", f_base, sep=" "))
    print(paste("Formula prunning:", f_prunning, sep=" "))
    #print(paste("Formula GSE1:", gse1$intervention_conditioning_set[1], sep=" "))
    #print(paste("Formula GSE2:", gse2$intervention_conditioning_set[1], sep=" "))
    print(paste("Formula GSE1:",f_gse1, sep=" "))
    print(paste("Formula GSE2:",f_gse2, sep=" "))
    print(paste("Formula DRW:", f_dr, sep=" "))
  }
  
  # evaluate in tpop (external validation)
  results_null <- ext_validation(model=null_model, target=tpop)
  results_m1 <- ext_validation(model=m1, target=tpop)
  #results_iecv <- ext_validation(model=iecv, target=tpop)
  results_iecv1 <- ext_validation(model=iecv1, target=tpop)
  results_iecv2 <- ext_validation(model=iecv2, target=tpop)
  results_iecv3 <- ext_validation(model=iecv3, target=tpop)
  results_iecv4 <- ext_validation(model=iecv4, target=tpop)
  results_base <- ext_validation(model=m_base, target=tpop_val)
  
  # Graph surgery
  results_gse1 <- ext_validation_pred(y=tpop$y, y_pred=gse1$predicted_gse)
  results_gse2 <- ext_validation_pred(y=tpop$y, y_pred=gse2$predicted_gse)
  #results_gse1 <- ext_validation(model=gse1, target=tpop)
  #results_gse2 <- ext_validation(model=gse2, target=tpop)
  # Prunning
  results_prunning <- ext_validation(model=prunning, target=tpop)
  # DRW
  results_dr <- ext_validation(model=dr, target=tpop)
  
  print(paste("Model_m1, iter", iter, "cslope", results_m1$cslope, paste(round(coef(m1), 3), collapse=","   )))

# save predicted files
  paramname <- paste("iter", iter,"N", N, "tau", tau, "tau_beta", tau_beta, "tau_sd", tau_sd,  "p", p)

  data.frame(
    iter = iter,
    observed = tpop$y,
    predicted_m1 = predict(m1, newdata=tpop, type='response'),
    predicted_prunning = predict(prunning, newdata=tpop, type='response'),
    #predicted_iecv = predict(iecv, newdata=tpop, type='response'),
    predicted_iecv1 = predict(iecv1, newdata=tpop, type='response'),
    predicted_iecv2 = predict(iecv2, newdata=tpop, type='response'),
    predicted_iecv3 = predict(iecv3, newdata=tpop, type='response'),
    predicted_iecv4 = predict(iecv4, newdata=tpop, type='response'),
    predicted_dr = predict(dr, newdata=tpop, type='response'),
    predicted_gse1 = gse1$predicted_gse,
    predicted_gse2 = gse2$predicted_gse
    #predicted_gse1 = predict(gse1, newdata=tpop, type='response')
    #predicted_gse2 = predict(gse2, newdata=tpop, type='response')
  ) %>%
  saveRDS(paste("~/scratch/simulated_predicted/", type, "/", structure,"/", paramname,".rds", sep=""))

  # collect results into a list
  
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
    #auc_iecv = results_iecv$auc,
    #brier_iecv = results_iecv$brier,
    #citl_iecv = results_iecv$citl,
    #cslope_iecv = results_iecv$cslope,
    auc_iecv1 = results_iecv1$auc,
    brier_iecv1 = results_iecv1$brier,
    citl_iecv1 = results_iecv1$citl,
    cslope_iecv1 = results_iecv1$cslope,
    auc_iecv2 = results_iecv2$auc,
    brier_iecv2 = results_iecv2$brier,
    citl_iecv2 = results_iecv2$citl,
    cslope_iecv2 = results_iecv2$cslope,
    auc_iecv3 = results_iecv3$auc,
    brier_iecv3 = results_iecv3$brier,
    citl_iecv3 = results_iecv3$citl,
    cslope_iecv3 = results_iecv3$cslope,
    auc_iecv4 = results_iecv4$auc,
    brier_iecv4 = results_iecv4$brier,
    citl_iecv4 = results_iecv4$citl,
    cslope_iecv4 = results_iecv4$cslope,

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

    
    # prunning
    auc_prunning = results_prunning$auc,
    brier_prunning = results_prunning$brier,
    citl_prunning = results_prunning$citl,
    cslope_prunning = results_prunning$cslope,

    # DRW
    auc_dr = results_dr$auc,
    brier_dr = results_dr$brier,
    citl_dr = results_dr$citl,
    cslope_dr = results_dr$cslope,
    

    f_m1 = as.character(f_m1), 
    f_iecv1 = f_iecv1, 
    f_iecv2 = f_iecv2, 
    f_iecv3 = f_iecv3, 
    f_iecv4 = f_iecv4, 
    #f_base = f_base, 
    f_prunning = f_prunning,
    #f_gse1 = gse1$intervention_conditioning_set[1],
    #f_gse2 = gse2$intervention_conditioning_set[1],
    f_gse1 = f_gse1,
    f_gse2 = f_gse2,
    f_dr = f_dr

    #dpop_mean = mean(dpop$x1),
    #dpop_sd = sd(dpop$x1),
    #tpop_mean = mean(tpop$x1),
    #tpop_sd = sd(tpop$x1)
  )
  
  # combine results and all information and return this each iteration
  #result <- c(iter = iter, 
  #            data_path = iter_fullpath,
  #            result, 
  #            f_m1 = as.character(f_m1), 
  #            f_iecv1 = f_iecv1, 
  #            f_iecv2 = f_iecv2, 
  #            f_iecv3 = f_iecv3, 
  #            f_iecv4 = f_iecv4, 
  #            #f_base = f_base, 
  #            f_prunning = f_prunning,
  #            #f_gse1 = gse1$intervention_conditioning_set[1],
  #            #f_gse2 = gse2$intervention_conditioning_set[1],
  #            f_gse1 = f_gse1,
  #            f_gse2 = f_gse2,
  #            f_dr = f_dr
  #            )
fname_string <- paste("n_iter", n_iter, "N", N, "tau", tau, "tau_beta", tau_beta, "tau_sd", tau_sd, "p", p, "K", K, "seed", seed)
folder <- paste("output-multiple/", type, "/", structure,"/", fname_string, sep = "")
if (!dir.exists(folder)) {
  dir.create(folder)
}
saveRDS(result, paste(folder, "/", iter,".rds", sep=""))

print(paste("Finish Iter", iter))
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



