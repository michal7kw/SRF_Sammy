#!/bin/bash
#SBATCH --job-name=1_analyze_chromatin_changes_GFP_M2
#SBATCH --output=logs/1_analyze_chromatin_changes_GFP_M2_%a.out
#SBATCH --error=logs/1_analyze_chromatin_changes_GFP_M2_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3:00:00  # Reduced time since we're skipping p-value calculations
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --array=1-21  # 19 autosomes + X + Y for mouse (mm10)

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
log "Starting chromatin analysis job for array task ${SLURM_ARRAY_TASK_ID}"

# Map array task ID to chromosome
if [ ${SLURM_ARRAY_TASK_ID} -le 19 ]; then
    CHROM="${SLURM_ARRAY_TASK_ID}"
elif [ ${SLURM_ARRAY_TASK_ID} -eq 20 ]; then
    CHROM="X"
else
    CHROM="Y"
fi

log "Processing chromosome: ${CHROM}"

# Activate conda environment
log "Activating conda environment"
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/1_chromatin_changes_GFP_M2"
log "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log "ERROR: Failed to change to working directory"; exit 1; }

# Check if the Python script exists
if [ ! -f "1_analyze_chromatin_changes_GFP_M2.py" ]; then
    log "ERROR: 1_analyze_chromatin_changes_GFP_M2.py not found in $WORKDIR"
    exit 1
fi

# Run the analysis for this chromosome
log "Starting Python analysis script for chromosome ${CHROM}"
python 1_analyze_chromatin_changes_GFP_M2.py --chrom "${CHROM}"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    log "Analysis completed successfully for chromosome ${CHROM}"
else
    log "ERROR: Analysis script failed for chromosome ${CHROM}"
    exit 1
fi
