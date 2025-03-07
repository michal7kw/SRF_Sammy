#!/bin/bash
#SBATCH --job-name=plots_GFP_M2
#SBATCH --output=logs/plots_GFP_M2.out
#SBATCH --error=logs/plots_GFP_M2.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs
mkdir -p /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

python plots_GFP_M2.py \
  --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation \
  --output-file /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation/chromatin_transitions.png