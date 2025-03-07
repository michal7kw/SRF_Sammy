#!/bin/bash
#SBATCH --job-name=plots_shift_GFP_M2
#SBATCH --output=logs/plots_shift_GFP_M2.out
#SBATCH --error=logs/plots_shift_GFP_M2.err
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

# python plots_shift_GPF_M2_binary.py \
#   --genome-size 2652783500 \
#   --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation \
#   --output-file /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation/chromatin_changes.png

# python plots_shift_GPF_M2_normal.py \
#   --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation \
#   --method ratio_comparison \
#   --debug \
#   --genome-size 2652783500 \
#   --output-file /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation/chromatin_changes.png

# python plots_shift_GPF_M2_loosen.py \
#   --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation \
#   --ratio-change-threshold 0.2 \
#   --debug \
#   --genome-size 2652783500 \
#   --output-file /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation/chromatin_changes.png
  

python plots_shift_GPF_M2_loosen_with_window.py \
  --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation \
  --ratio-change-threshold 0.3 \
  --window-size 150000 \
  --debug \
  --genome-size 2652783500 \
  --output-file /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation/chromatin_changes_150kb.png