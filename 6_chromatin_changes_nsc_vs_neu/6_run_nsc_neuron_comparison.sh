#!/bin/bash
#SBATCH --job-name=nsc_neuron_comparison
#SBATCH --output=logs/nsc_neuron_comparison.out
#SBATCH --error=logs/nsc_neuron_comparison.err
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
log "Starting NSC vs Neuron comparison analysis"

# Activate conda environment
log "Activating conda environment"
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"
log "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log "ERROR: Failed to change to working directory"; exit 1; }

# Check if the Python script exists
if [ ! -f "compare_nsc_neuron.py" ]; then
    log "ERROR: compare_nsc_neuron.py not found in $WORKDIR"
    exit 1
fi

# Create output directory
mkdir -p results/nsc_neuron_comparison

# Run the comparison script
log "Running NSC vs Neuron comparison analysis"
python compare_nsc_neuron.py --analysis_dir results/analysis --output_dir results/nsc_neuron_comparison

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    log "NSC vs Neuron comparison analysis completed successfully"
else
    log "ERROR: NSC vs Neuron comparison analysis failed"
    exit 1
fi

log "Results saved to results/nsc_neuron_comparison/" 