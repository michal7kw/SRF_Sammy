#!/bin/bash
#SBATCH --job-name=full_chromatin_analysis
#SBATCH --output=logs/full_chromatin_analysis.out
#SBATCH --error=logs/full_chromatin_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
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
log "Starting full chromatin analysis"

# Activate conda environment
log "Activating conda environment"
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"
log "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log "ERROR: Failed to change to working directory"; exit 1; }

# Check if the Python scripts exist
if [ ! -f "1_analyze_chromatin_changes_GFP_M2.py" ]; then
    log "ERROR: 1_analyze_chromatin_changes_GFP_M2.py not found in $WORKDIR"
    exit 1
fi

if [ ! -f "compare_chromatin_states.py" ]; then
    log "ERROR: compare_chromatin_states.py not found in $WORKDIR"
    exit 1
fi

# Create output directories
mkdir -p results/analysis
mkdir -p results/comparison
mkdir -p results/nsc_neuron_comparison

# Step 1: Run the original analysis for all chromosomes
log "Running original chromatin analysis for all chromosomes"

for CHROM in {1..19} X Y; do
    log "Processing chromosome: ${CHROM}"
    python 1_analyze_chromatin_changes_GFP_M2.py --chrom "${CHROM}"
    
    # Check if the script completed successfully
    if [ $? -ne 0 ]; then
        log "ERROR: Chromatin analysis failed for chromosome ${CHROM}"
        exit 1
    fi
done

log "Original chromatin analysis completed successfully"

# Step 2: Run the comparison script
log "Running chromatin comparison analysis"
python compare_chromatin_states.py --analysis_dir results/analysis --output_dir results/comparison

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    log "Chromatin comparison analysis completed successfully"
else
    log "ERROR: Chromatin comparison analysis failed"
    exit 1
fi

# Step 3: Run the NSC vs Neuron comparison script
log "Running NSC vs Neuron comparison analysis"
python compare_nsc_neuron.py --analysis_dir results/analysis --output_dir results/nsc_neuron_comparison

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    log "NSC vs Neuron comparison analysis completed successfully"
else
    log "ERROR: NSC vs Neuron comparison analysis failed"
    exit 1
fi

log "Full analysis completed. Results saved to:"
log "- Original analysis: results/analysis/"
log "- Comparison results: results/comparison/"
log "- NSC vs Neuron comparison: results/nsc_neuron_comparison/" 