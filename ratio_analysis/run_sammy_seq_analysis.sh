#!/bin/bash
#SBATCH --job-name=run_sammy_seq_analysis
#SBATCH --output=logs/run_sammy_seq_analysis.out
#SBATCH --error=logs/run_sammy_seq_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
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

# Run the Python script with all necessary arguments
python sammy_seq_analysis.py \
  --workdir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis \
  --outdir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation \
  --bin-size 50000 \
  --eu-threshold 0.1 \
  --hetero-threshold -0.1 \
  --threads 16 \