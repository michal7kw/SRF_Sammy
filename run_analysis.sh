#!/bin/bash
#SBATCH --job-name=analyze_chromatin
#SBATCH --output=logs/analyze_chromatin.out
#SBATCH --error=logs/analyze_chromatin.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
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
log "Starting chromatin analysis job"

# Activate conda environment
log "Activating conda environment"
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Sammy"
log "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log "ERROR: Failed to change to working directory"; exit 1; }

# Check if the Python script exists
if [ ! -f "analyze_chromatin_changes.py" ]; then
    log "ERROR: analyze_chromatin_changes.py not found in $WORKDIR"
    exit 1
fi

# Run the analysis
log "Starting Python analysis script"
python analyze_chromatin_changes.py

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    log "Analysis completed successfully"
else
    log "ERROR: Analysis script failed"
    exit 1
fi
