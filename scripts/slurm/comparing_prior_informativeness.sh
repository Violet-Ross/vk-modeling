#!/bin/bash
#SBATCH --job-name=prior_info                      
#SBATCH --time=9-00:00:00                       
#SBATCH --partition=long-cascadelake       
#SBATCH --nodes=1                            
#SBATCH --ntasks=8                               ## RESERVE TWO CPUS/CORES; CHANGE IF NEEDED
#SBATCH --mem=12G                                ## RESERVE 12GB MEMORY (RAM); CHANGE IF NEEDED
#SBATCH --output=/scratch/isilon/viro3403/slurm_%j_%x.out 
#SBATCH --error=/scratch/isilon/viro3403/slurm_%j_%x.err
        

##  RUN YOUR CODE BELOW
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rstan-clean
cd ~/vk-modeling
Rscript scripts/experiments/comparing_prior_informativeness.R