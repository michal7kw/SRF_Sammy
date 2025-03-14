#!/bin/bash
#SBATCH --job-name=4_run_mecp2_repressed_genes_heatmap
#SBATCH --output=logs/4_run_mecp2_repressed_genes_heatmap.out
#SBATCH --error=logs/4_run_mecp2_repressed_genes_heatmap.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

"""
This script generates a metaprofile heatmap specifically for Mecp2 repressed genes,
showing the log2 ratio of S2S/S3 across gene bodies (from -5kb to TSS to TES to +5kb).

The heatmap uses a blue-white-red color scheme where:
- Blue represents heterochromatin (log2(S2S/S3) < 0)
- White represents undefined regions (log2(S2S/S3) = 0)
- Red represents euchromatin (log2(S2S/S3) > 0)

The script creates a heatmap similar to the one in the provided image, with rows for
different conditions (NSC2_GFP, NSC2_M2, NSC3_GFP, NSC3_M2) and columns representing
positions from -5kb through TSS to TES and +5kb.
"""

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

base_dir="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"

# Set variables
RESULTS_DIR="${base_dir}/ratio_analysis/results/shift_GFP_M2_with_blacklist_with_mirko_gene_lists"
OUTPUT_DIR="${base_dir}/ratio_analysis/results/mecp2_repressed_genes_heatmap"
GENOME_GTF="${base_dir}/DATA/gencode.vM25.basic.annotation.gtf"
GENE_LISTS_DIR="${base_dir}/figures/genelists"
UPSTREAM_DISTANCE=5000
DOWNSTREAM_DISTANCE=5000
BIN_SIZE=100
GENE_BODY_BINS=100

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Run the analysis
python 4_mecp2_binding_genes_heatmap.py \
    --results-dir ${RESULTS_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --genome-gtf ${GENOME_GTF} \
    --gene-lists-dir ${GENE_LISTS_DIR} \
    --upstream-distance ${UPSTREAM_DISTANCE} \
    --downstream-distance ${DOWNSTREAM_DISTANCE} \
    --bin-size ${BIN_SIZE} \
    --gene-body-bins ${GENE_BODY_BINS}

echo "Mecp2 repressed genes heatmap generation completed." 