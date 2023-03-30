#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 24:00:00   ## 2 hr run time limit
#SBATCH --mem=16G
#SBATCH -o fit_SEIRS_gmrf_immunity_models-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=34-48
#SBATCH --exclude=hpc3-15-29,hpc3-21-30,hpc3-21-23

module purge
module load julia/1.8.5
cd //pub/bayerd/immunity_semi_parametric_model/

julia --project --threads 4 scripts/simulation/compare_inference_known_variant_proportion/SEIRS_GMRF_immunity/fit_SEIRS_GMRF_immunity_model.jl $SLURM_ARRAY_TASK_ID
