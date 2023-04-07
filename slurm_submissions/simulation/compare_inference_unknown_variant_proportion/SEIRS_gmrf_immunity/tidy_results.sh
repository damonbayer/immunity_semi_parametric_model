#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 1 task (1 CPUs)
#SBATCH -t 01:00:00   ## 1 hr run time limit
#SBATCH --mem=3G
#SBATCH -o tidy_results-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=20-35

module purge
module load R
cd //pub/bayerd/immunity_semi_parametric_model

Rscript scripts/simulation/compare_inference_unknown_variant_proportion/SEIRS_GMRF_immunity/tidy_results.R $SLURM_ARRAY_TASK_ID
