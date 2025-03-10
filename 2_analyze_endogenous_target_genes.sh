#!/bin/bash
#SBATCH --job-name=2_analyze_endogenous_target_genes
#SBATCH --output=logs/2_analyze_endogenous_target_genes_%a.out
#SBATCH --error=logs/2_analyze_endogenous_target_genes_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3:00:00
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
log "Starting endogenous Mecp2 target gene analysis for array task ${SLURM_ARRAY_TASK_ID}"

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
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"
log "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log "ERROR: Failed to change to working directory"; exit 1; }

# Check if the Python script exists
if [ ! -f "2_analyze_endogenous_target_genes.py" ]; then
    log "ERROR: 2_analyze_endogenous_target_genes.py not found in $WORKDIR"
    exit 1
fi

# Check if gene list files exist
if [ ! -f "gene_lists/endo_gene_list_neu.txt" ] || [ ! -f "gene_lists/endo_gene_list_nsc.txt" ]; then
    log "ERROR: Gene list files not found in results directory"
    exit 1
fi

# Check if GTF file exists
GTF_FILE="./gencode.vM25.basic.annotation.gtf"
if [ ! -f "$GTF_FILE" ]; then
    log "ERROR: GTF file not found: $GTF_FILE"
    exit 1
fi

# Run the analysis for this chromosome
log "Starting Python analysis script for chromosome ${CHROM}"
python 2_analyze_endogenous_target_genes.py --chrom "${CHROM}" --gtf_file "$GTF_FILE"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    log "Analysis completed successfully for chromosome ${CHROM}"
else
    log "ERROR: Analysis script failed for chromosome ${CHROM}"
    exit 1
fi 