#!/bin/bash
#SBATCH --job-name=restricted_prior
#SBATCH --time=4:00:00
#SBATCH --partition=long-cascadelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4                        ## 4 cores for Stan's 4 chains
#SBATCH --mem=8G
#SBATCH --array=1-150                            ## 50 prior_sd values x 3 l_0 values
#SBATCH --output=/scratch/isilon/viro3403/slurm_%A_%a_%x.out
#SBATCH --error=/scratch/isilon/viro3403/slurm_%A_%a_%x.err

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rstan-clean
cd ~/vk-modeling

Rscript scripts/experiments/parallel_comparing_prior/restricted/run_condition.R