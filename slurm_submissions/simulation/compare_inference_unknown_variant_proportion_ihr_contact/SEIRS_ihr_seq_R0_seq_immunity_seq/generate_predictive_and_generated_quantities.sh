#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 00:10:00   ## 2 hr run time limit
#SBATCH --mem=4G
#SBATCH -o generate_predictive_and_generated_quantities-%A_%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=28-37

module purge
module load julia/1.8.5
cd //pub/bayerd/immunity_semi_parametric_model/

if [ $SLURM_ARRAY_TASK_ID == 28 ]; then
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/simulation/compare_inference_unknown_variant_proportion_ihr_contact/SEIRS_ihr_seq_R0_seq_immunity_seq/tidy_results.sh
fi

julia --project --threads 1 scripts/simulation/compare_inference_unknown_variant_proportion_ihr_contact/SEIRS_ihr_seq_R0_seq_immunity_seq/generate_predictive_and_generated_quantities.jl $SLURM_ARRAY_TASK_ID
