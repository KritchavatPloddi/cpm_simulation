Simulation steps

1. Generate intercept and betas (generate_betas.R)
2. Generate subpopulations by taking generated betas and applydataset shift settings and other parameters (simulate_datasets.R/function_generate_clusters.R, done by submitting jobs at jobs/simulate_data.sh)
3. Apply models and store results (simulation_cmd.R/function_simulation_nrun.R, done by submitting jobs at jobs/run_multiple.sh)
