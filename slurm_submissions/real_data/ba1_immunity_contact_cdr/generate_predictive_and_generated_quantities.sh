#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 1:00:00   ## 2 hr run time limit
#SBATCH --mem=4G
#SBATCH -o generate_predictive_and_generated_quantities-%A_%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=71-126

module purge
module load julia/1.8.5
cd //pub/bayerd/immunity_semi_parametric_model/

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/real_data/ba1_immunity_contact_cdr/tidy_results.sh
fi

julia --project scripts/real_data/ba1_immunity_contact_cdr/generate_predictive_and_generated_quantities.jl $SLURM_ARRAY_TASK_ID
julia --project scripts/real_data/ba1_immunity_contact_cdr/generate_peak_quantities.jl $SLURM_ARRAY_TASK_ID
