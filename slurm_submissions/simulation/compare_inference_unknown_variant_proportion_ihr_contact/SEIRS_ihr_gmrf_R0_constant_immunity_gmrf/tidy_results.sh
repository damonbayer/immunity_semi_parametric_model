#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 1 task (1 CPUs)
#SBATCH -t 01:00:00   ## 1 hr run time limit
#SBATCH --mem=3G
#SBATCH -o tidy_results-%A_%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=28-37

module purge
module load R
cd //pub/bayerd/immunity_semi_parametric_model

Rscript scripts/simulation/compare_inference_unknown_variant_proportion_ihr_contact/SEIRS_ihr_gmrf_R0_constant_immunity_gmrf/tidy_results.R $SLURM_ARRAY_TASK_ID
