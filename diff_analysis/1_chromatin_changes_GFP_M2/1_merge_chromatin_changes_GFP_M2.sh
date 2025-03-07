#!/bin/bash
#SBATCH --job-name=1_merge_chromatin_changes_GFP_M2
#SBATCH --output=logs/1_merge_chromatin_changes_GFP_M2.out
#SBATCH --error=logs/1_merge_chromatin_changes_GFP_M2.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:30:00  # Reduced time since we're skipping p-value calculations
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
log "Starting results merge job"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/1_chromatin_changes_GFP_M2"
log "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log "ERROR: Failed to change to working directory"; exit 1; }

# Check if the Python script exists
if [ ! -f "1_merge_chromatin_changes_GFP_M2.py" ]; then
    log "ERROR: 1_merge_chromatin_changes_GFP_M2.py not found in $WORKDIR"
    exit 1
fi

# Debug: Check directory structure before merging
log "Checking directory structure before merging"
RESULTS_DIR="$WORKDIR/../results/analysis"
log "Results directory: $RESULTS_DIR"

# Check for chromosome directories
log "Checking for chromosome directories:"
for CHROM in {1..19} X Y; do
    CHROM_DIR="$RESULTS_DIR/chr$CHROM/chromosomes/chr$CHROM"
    if [ -d "$CHROM_DIR" ]; then
        log "  Found directory: $CHROM_DIR"
        log "  Files in directory:"
        ls -la "$CHROM_DIR" | grep -v "^total" | awk '{print "    "$9}'
    else
        log "  Directory not found: $CHROM_DIR"
    fi
done

# Run the merge operation
log "Starting merge operation"
python 1_merge_chromatin_changes_GFP_M2.py --results_dir results

# Check if the script completed successfully
MERGE_STATUS=$?
if [ $MERGE_STATUS -eq 0 ]; then
    log "Merge completed successfully"
    
    # Check merged results
    log "Checking merged results:"
    
    # Check Neu results
    if [ -f "$RESULTS_DIR/Neu_chromatin_changes_all.csv" ]; then
        NEU_LINES=$(wc -l < "$RESULTS_DIR/Neu_chromatin_changes_all.csv")
        log "  Neu results found: $NEU_LINES lines"
        
        if [ -f "$RESULTS_DIR/Neu_chromatin_changes_significant.csv" ]; then
            NEU_SIG_LINES=$(wc -l < "$RESULTS_DIR/Neu_chromatin_changes_significant.csv")
            log "  Neu significant results: $NEU_SIG_LINES lines"
        else
            log "  WARNING: Neu significant results file not found"
        fi
    else
        log "  ERROR: Neu results not found"
    fi
    
    # Check NSC results
    if [ -f "$RESULTS_DIR/NSC_chromatin_changes_all.csv" ]; then
        NSC_LINES=$(wc -l < "$RESULTS_DIR/NSC_chromatin_changes_all.csv")
        log "  NSC results found: $NSC_LINES lines"
        
        if [ -f "$RESULTS_DIR/NSC_chromatin_changes_significant.csv" ]; then
            NSC_SIG_LINES=$(wc -l < "$RESULTS_DIR/NSC_chromatin_changes_significant.csv")
            log "  NSC significant results: $NSC_SIG_LINES lines"
        else
            log "  WARNING: NSC significant results file not found"
        fi
    else
        log "  ERROR: NSC results not found"
    fi
    
    # List all CSV files in the results directory
    log "All CSV files in results directory:"
    find "$RESULTS_DIR" -name "*.csv" -maxdepth 1 -type f | sort | while read -r file; do
        log "  $(basename "$file"): $(wc -l < "$file") lines"
    done
else
    log "ERROR: Merge operation failed with status $MERGE_STATUS"
    exit 1
fi 