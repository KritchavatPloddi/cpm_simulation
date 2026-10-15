#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -n 24
#SBATCH -t 0-10
#SBATCH -a 0-24
#SBATCH --output=logs/multiple_%a.out 

source ~/.bashrc 
conda activate r
    
# Define the list of Rscript commands
CMD_LIST=(
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="observed" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="observed" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=2 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="observed" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="observed" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=2 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="observed" --cores=24"

"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="unobserved" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="unobserved" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=2 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="unobserved" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="unobserved" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=2 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="unobserved" --cores=24"
    
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=1 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=2 --tau_beta=0 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=1 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=2 --tau_sd=0 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"

"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=0 --tau_sd=1 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=1 --tau_beta=0 --tau_sd=1 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=2 --tau_beta=0 --tau_sd=1 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=1 --tau_sd=1 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=2 --tau_sd=1 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"

"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=0 --tau_sd=2 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=1 --tau_beta=0 --tau_sd=2 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=2 --tau_beta=0 --tau_sd=2 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=1 --tau_sd=2 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"
"Rscript ~/simulation/simulation_cmd.R --type="c" --tau=0 --tau_beta=2 --tau_sd=2 --dpop_num=9 --N=1000 --n_iter=1000 --structure="proxy" --cores=24"

)





# Run the command for this array index
eval ${CMD_LIST[$SLURM_ARRAY_TASK_ID]}
