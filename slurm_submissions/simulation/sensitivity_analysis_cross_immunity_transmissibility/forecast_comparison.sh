#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 20          ## request 4 tasks (4 CPUs)
#SBATCH -t 1:00:00   ## 2 hr run time limit
#SBATCH --mem=16G
#SBATCH -o forecast_comparison-%A_%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=1-3

module purge
module load R
cd //pub/bayerd/immunity_semi_parametric_model/


Rscript scripts/simulation/sensitivity_analysis_cross_immunity_transmissibility/forecast_comparison.R $SLURM_ARRAY_TASK_ID
