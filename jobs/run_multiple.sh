#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -n 30
#SBATCH -t 0-00:40
#SBATCH -a 0-24
#SBATCH --output=logs/multiple_%a.out 

source ~/.bashrc 
conda activate r
    
# Define the list of Rscript commands
CMD_LIST=(

"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="observed" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0.5 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="observed" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="observed" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="unobserved" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0.5 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="unobserved" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="unobserved" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0.5 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0 --tau_sd=1 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0.5 --tau_beta=0 --tau_sd=1 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=1 --tau_beta=0 --tau_sd=1 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0 --tau_sd=2 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0.5 --tau_beta=0 --tau_sd=2 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=1 --tau_beta=0 --tau_sd=2 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"

"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0.5 --tau_sd=0 --N=1000 --n_iter=1000 --structure="observed" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --N=1000 --n_iter=1000 --structure="observed" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0.5 --tau_sd=0 --N=1000 --n_iter=1000 --structure="unobserved" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --N=1000 --n_iter=1000 --structure="unobserved" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0.5 --tau_sd=0 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0.5 --tau_sd=1 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=1 --tau_sd=1 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=0.5 --tau_sd=2 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"
"Rscript ~/simulation/simulation_cmd.R --development=-1 --target=1 --type="c" --tau=0 --tau_beta=1 --tau_sd=2 --N=1000 --n_iter=1000 --structure="proxy" --cores=30"

)





# Run the command for this array index
eval ${CMD_LIST[$SLURM_ARRAY_TASK_ID]}
