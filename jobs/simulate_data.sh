#!/bin/bash --login
#SBATCH -p serial 
#SBATCH -t 0-1
#SBATCH -a 0-24
#SBATCH --output=logs/simulate_data_%a.out
conda activate r
CMD_LIST=(
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="observed""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0.5 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="observed""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="observed""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="unobserved""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0.5 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="unobserved""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="unobserved""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0.5 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0 --tau_sd=1 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0.5 --tau_beta=0 --tau_sd=1 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=1 --tau_beta=0 --tau_sd=1 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0 --tau_sd=2 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0.5 --tau_beta=0 --tau_sd=2 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=1 --tau_beta=0 --tau_sd=2 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""

"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0.5 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="observed""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="observed""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0.5 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="unobserved""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="unobserved""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0.5 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0.5 --tau_sd=1 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=1 --tau_sd=1 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=0.5 --tau_sd=2 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
"Rscript ~/simulation/simulate_datasets.R --type="c" --tau=0 --tau_beta=1 --tau_sd=2 --N=1000 --n_iter=1000 --p=0.2 --structure="proxy""
)

# Run the command for this array index
eval ${CMD_LIST[$SLURM_ARRAY_TASK_ID]}