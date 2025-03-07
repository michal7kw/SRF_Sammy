#!/bin/bash
#SBATCH --job-name=run_gene_chromatin_state_comparison_neu
#SBATCH --output=logs/run_gene_chromatin_state_comparison_neu.out
#SBATCH --error=logs/run_gene_chromatin_state_comparison_neu.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

# Set variables
GENE_LIST="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/gene_lists/neu_enriched_gene_list_exo_vs_endo.txt"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results"
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/neu_enriched_comparison"
GENOME_GTF="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/DATA/gencode.vM25.basic.annotation.gtf"
PROMOTER_SIZE=2000

mkdir -p ${OUTPUT_DIR}

# Run the analysis
python gene_chromatin_state_comparison.py \
    --gene-list ${GENE_LIST} \
    --results-dir ${RESULTS_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --genome-gtf ${GENOME_GTF} \
    --promoter-size ${PROMOTER_SIZE}

echo "Analysis completed."
