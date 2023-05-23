#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 1:00:00   ## 2 hr run time limit
#SBATCH --mem=4G
#SBATCH -o tidy_results-%A_%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=1-56

module purge
module load R
cd //pub/bayerd/immunity_semi_parametric_model/


Rscript scripts/real_data/ba1_immunity_contact_cdr/tidy_results.R $SLURM_ARRAY_TASK_ID
