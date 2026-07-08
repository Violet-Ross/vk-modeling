#!/bin/bash
#SBATCH --job-name=acquisition
#SBATCH --time=8:00:00
#SBATCH --partition=long-cascadelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4                        ## 4 cores for Stan's 4 chains
#SBATCH --mem=8G
#SBATCH --array=1-300                            ## Update limit after running setup script
#SBATCH --output=/scratch/isilon/viro3403/slurm_%A_%a_%x.out
#SBATCH --error=/scratch/isilon/viro3403/slurm_%A_%a_%x.err

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rstan-clean
cd ~/vk-modeling

Rscript scripts/experiments/parallel_comparing_acquisition/run_acquisition_condition.R