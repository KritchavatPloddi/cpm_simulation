Simulation steps

1. Generate overarching population without dataset shift (generate_overarching.R).
2. Generate subpopulations by modifying the overarching population with dataset shift settings and other parameters (generate_clusters.R, done by submitting jobs at jobs/simulate_data.sh)
3. Apply models and store results (simulation_cmd.R, done by submitting jobs at jobs/run_multiple.sh)
