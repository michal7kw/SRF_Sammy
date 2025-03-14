#!/bin/bash
#SBATCH --job-name=3_run_mecp2_metaprofile_heatmap
#SBATCH --output=logs/3_run_mecp2_metaprofile_heatmap.out
#SBATCH --error=logs/3_run_mecp2_metaprofile_heatmap.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

"""
This script generates metaprofile heatmaps for Mecp2 target and non-target genes in NSC and NEU conditions.

It creates heatmaps showing the log2 ratio of S2S/S3 across gene bodies (from -5kb to TSS to TES to +5kb)
using a blue-white-red color scheme where:
- Blue represents heterochromatin (log2(S2S/S3) < -0.1)
- White represents undefined regions (-0.1 <= log2(S2S/S3) <= 0.1)
- Red represents euchromatin (log2(S2S/S3) > 0.1)

The script processes the following gene lists:
- NEU_exo.txt: Mecp2 target genes in NEU cells
- no_NEU_exo.txt: Non-target genes in NEU cells
- NSC_exo.txt: Mecp2 target genes in NSC cells
- no_NSC_exo.txt: Non-target genes in NSC cells

For each gene list, it creates a metaprofile heatmap showing the chromatin state across all conditions.
It also creates a combined heatmap for each condition showing the comparison between different gene lists.
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
OUTPUT_DIR="${base_dir}/ratio_analysis/results/mecp2_metaprofile_heatmaps"
GENOME_GTF="${base_dir}/DATA/gencode.vM25.basic.annotation.gtf"
GENE_LISTS_DIR="${base_dir}/figures/genelists"
UPSTREAM_DISTANCE=5000
DOWNSTREAM_DISTANCE=5000
BIN_SIZE=100
GENE_BODY_BINS=100

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Run the analysis
python 3_create_mecp2_metaprofile_heatmap.py \
    --results-dir ${RESULTS_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --genome-gtf ${GENOME_GTF} \
    --gene-lists-dir ${GENE_LISTS_DIR} \
    --upstream-distance ${UPSTREAM_DISTANCE} \
    --downstream-distance ${DOWNSTREAM_DISTANCE} \
    --bin-size ${BIN_SIZE} \
    --gene-body-bins ${GENE_BODY_BINS}

echo "Metaprofile heatmap generation completed." 