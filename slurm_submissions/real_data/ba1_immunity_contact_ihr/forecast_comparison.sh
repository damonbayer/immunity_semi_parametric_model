#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 20          ## request 20 tasks (4 CPUs)
#SBATCH -t 01:00:00   ## 1 hr run time limit
#SBATCH --mem=16G
#SBATCH -o forecast_comparison-%A_%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu

module purge
module load R
cd //pub/bayerd/immunity_semi_parametric_model

Rscript scripts/real_data/ba1_immunity_contact_ihr/forecast_comparison.R
