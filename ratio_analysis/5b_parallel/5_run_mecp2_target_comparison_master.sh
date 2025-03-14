#!/bin/bash
#SBATCH --job-name=5_mecp2_master
#SBATCH --output=logs/5_mecp2_master.out
#SBATCH --error=logs/5_mecp2_master.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=10:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Master script for the parallelized Mecp2 target comparison heatmap generation.
#
# This script:
# 1. Creates necessary directories
# 2. Submits the array job for parallel metaprofile computation
# 3. The array job will automatically submit the plot generation job when all tasks complete

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

base_dir="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"

# Set variables
OUTPUT_DIR="${base_dir}/ratio_analysis/results/mecp2_target_comparison_heatmaps_parallel"

# Create output directory
mkdir -p ${OUTPUT_DIR}/temp

echo "Starting parallelized Mecp2 target comparison heatmap generation..."

# Submit the array job
sbatch 5_run_mecp2_target_comparison_heatmap_parallel.sh

echo "Submitted array job for parallel metaprofile computation."
echo "The plot generation job will be automatically submitted when all array tasks complete."
echo "Check the status of the jobs with: squeue -u \$USER"
echo "Final results will be available in: ${OUTPUT_DIR}" 