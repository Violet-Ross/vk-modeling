#!/bin/bash
#SBATCH --job-name=prior_info
#SBATCH --time=4:00:00                           ## Per-task time limit (not total)
#SBATCH --partition=long-cascadelake 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4                        ## 4 cores for Stan's 4 chains
#SBATCH --mem=8G                                 ## Per-task memory
#SBATCH --array=1-750                            ## One task per (prior_sd, l_0) pair
#SBATCH --output=/scratch/isilon/viro3403/slurm_%A_%a_%x.out   ## %A = array job ID, %a = task ID
#SBATCH --error=/scratch/isilon/viro3403/slurm_%A_%a_%x.err

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rstan-clean
cd ~/vk-modeling

Rscript scripts/experiments/comparing_prior_onephase/run_condition.R
