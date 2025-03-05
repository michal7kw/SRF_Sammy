#!/bin/bash
#SBATCH --job-name=2_merge_endogenous_target_results
#SBATCH --output=logs/2_merge_endogenous_target_results.out
#SBATCH --error=logs/2_merge_endogenous_target_results.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
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
log "Starting merge of endogenous Mecp2 target gene analysis results"

# Activate conda environment
log "Activating conda environment"
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/2_chromatin_changes_endogenous_targets"
log "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log "ERROR: Failed to change to working directory"; exit 1; }

# Check if the Python script exists
if [ ! -f "2_merge_endogenous_target_results.py" ]; then
    log "ERROR: 2_merge_endogenous_target_results.py not found in $WORKDIR"
    exit 1
fi

# Define results directory
RESULTS_DIR="../results"
if [ ! -d "$RESULTS_DIR" ]; then
    log "ERROR: Results directory not found: $RESULTS_DIR"
    exit 1
fi

# Check if the endogenous_target_analysis directory exists
TARGET_DIR="${RESULTS_DIR}/endogenous_target_analysis"
if [ ! -d "$TARGET_DIR" ]; then
    log "ERROR: Target analysis directory not found: $TARGET_DIR"
    exit 1
fi

# Run the merging script
log "Starting Python merging script"
python 2_merge_endogenous_target_results.py --results_dir "$RESULTS_DIR"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    log "Merging completed successfully"
else
    log "ERROR: Merging script failed"
    exit 1
fi 