#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 6:00:00   ## 2 hr run time limit
#SBATCH --mem=4G
#SBATCH -o fit_SEIRS_variant_models-%A_%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=28-48

module purge
module load julia/1.8.5
cd //pub/bayerd/immunity_semi_parametric_model/

julia --project --threads 4 scripts/simulation/compare_inference_known_variant_proportion/SEIRS_variant/fit_SEIRS_variant_model.jl $SLURM_ARRAY_TASK_ID