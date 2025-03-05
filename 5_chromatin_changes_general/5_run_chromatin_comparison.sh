#!/bin/bash
#SBATCH --job-name=chromatin_comparison
#SBATCH --output=logs/chromatin_comparison.out
#SBATCH --error=logs/chromatin_comparison.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Error handling
set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Create logs directory if it doesn't exist
mkdir -p logs

# Start logging
log "Starting chromatin comparison analysis"

# Activate conda environment
log "Activating conda environment"
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/5_chromatin_changes_general"
log "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log "ERROR: Failed to change to working directory"; exit 1; }

# Check if the Python script exists
if [ ! -f "5_compare_chromatin_states.py" ]; then
    log "ERROR: 5_compare_chromatin_states.py not found in $WORKDIR"
    exit 1
fi

# Create output directory
mkdir -p ../results/comparison

# Run the comparison script
log "Running chromatin comparison analysis"
python 5_compare_chromatin_states.py --analysis_dir ../results/analysis --output_dir ../results/comparison

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    log "Chromatin comparison analysis completed successfully"
else
    log "ERROR: Chromatin comparison analysis failed"
    exit 1
fi

log "Results saved to ../results/comparison/" 