#!/bin/bash
#SBATCH -A introtogds
#SBATCH -p normal_q
#SBATCH --cpus-per-task=8
#SBATCH -J Integration
#SBATCH --mail-user=arazan@vt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o Integration_%j.out


echo "Start date: $(date)"
echo "Running on host: $HOSTNAME"
scontrol show job --details $SLURM_JOB_ID

# Load R module
module load R/4.4.2-gfbf-2024a

# Set custom R library path
export R_LIBS_USER="/home/arazan/R/owl-genoa/4.4.2"

# Confirm path and Rscript
echo "Using Rscript at: $(which Rscript)"
echo "R_LIBS_USER is set to: $R_LIBS_USER"

# Run R script
Rscript Integration_Pathogen_NonPathogen.R
