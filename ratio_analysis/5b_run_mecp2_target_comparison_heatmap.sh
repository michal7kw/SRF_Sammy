#!/bin/bash
#SBATCH --job-name=5_run_mecp2_target_comparison
#SBATCH --output=logs/5_run_mecp2_target_comparison.out
#SBATCH --error=logs/5_run_mecp2_target_comparison.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

"""
This script generates a comprehensive comparison of metaprofile heatmaps for Mecp2 target and non-target genes
in NSC and NEU conditions, showing the log2 ratio of S2S/S3 across gene bodies.

It creates a 2x2 grid of heatmaps for:
1. NEU Mecp2 target genes (NEU_exo.txt)
2. NEU non-target genes (not_NEU_exo.txt)
3. NSC Mecp2 target genes (NSC_exo.txt)
4. NSC non-target genes (not_NSC_exo.txt)

Each heatmap shows the chromatin state across all conditions (NSC_GFP, NSC_M2, Neu_GFP, Neu_M2)
using a blue-white-red color scheme where:
- Blue represents heterochromatin (log2(S2S/S3) < 0)
- White represents undefined regions (log2(S2S/S3) = 0)
- Red represents euchromatin (log2(S2S/S3) > 0)

The script also generates line plots to better visualize differences between conditions.
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
OUTPUT_DIR="${base_dir}/ratio_analysis/results/mecp2_target_comparison_heatmaps"
GENOME_GTF="${base_dir}/DATA/gencode.vM25.basic.annotation.gtf"
GENE_LISTS_DIR="${base_dir}/figures/genelists"
UPSTREAM_DISTANCE=5000
DOWNSTREAM_DISTANCE=5000
BIN_SIZE=100
GENE_BODY_BINS=100

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Install required packages if needed
pip install scipy matplotlib==3.5.3 seaborn==0.12.2

# Run the analysis
python 5_create_mecp2_target_comparison_heatmap.py \
    --results-dir ${RESULTS_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --genome-gtf ${GENOME_GTF} \
    --gene-lists-dir ${GENE_LISTS_DIR} \
    --upstream-distance ${UPSTREAM_DISTANCE} \
    --downstream-distance ${DOWNSTREAM_DISTANCE} \
    --bin-size ${BIN_SIZE} \
    --gene-body-bins ${GENE_BODY_BINS}

echo "Mecp2 target comparison heatmap generation completed."

# Create a directory for additional visualizations
ENHANCED_OUTPUT_DIR="${OUTPUT_DIR}/enhanced_visualizations"
mkdir -p ${ENHANCED_OUTPUT_DIR}

# Copy the generated files to the enhanced directory
cp ${OUTPUT_DIR}/mecp2_target_comparison_heatmap.png ${ENHANCED_OUTPUT_DIR}/
cp ${OUTPUT_DIR}/mecp2_target_comparison_heatmap_line_plots.png ${ENHANCED_OUTPUT_DIR}/

echo "Enhanced visualizations are available in ${ENHANCED_OUTPUT_DIR}"