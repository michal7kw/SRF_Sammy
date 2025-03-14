#!/bin/bash
#SBATCH --job-name=5_mecp2_parallel
#SBATCH --output=logs/5_mecp2_parallel_%A_%a.out
#SBATCH --error=logs/5_mecp2_parallel_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --array=0-15

# This script parallelizes the computation of metaprofiles for Mecp2 target and non-target genes
# across multiple SLURM jobs using array jobs.
#
# The array job is split into 16 tasks (0-15), each handling a specific combination of:
# - Gene set (4 options): NEU_target, NEU_non_target, NSC_target, NSC_non_target
# - Condition (4 options per gene set): NSC_GFP, NSC_M2, Neu_GFP, Neu_M2
#
# After all array jobs complete, a separate script combines the results and generates the final heatmaps.

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

base_dir="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"

# Set variables
RESULTS_DIR="${base_dir}/ratio_analysis/results/shift_GFP_M2_with_blacklist_with_mirko_gene_lists"
OUTPUT_DIR="${base_dir}/ratio_analysis/results/mecp2_target_comparison_heatmaps_parallel"
GENOME_GTF="${base_dir}/DATA/gencode.vM25.basic.annotation.gtf"
GENE_LISTS_DIR="${base_dir}/figures/genelists"
UPSTREAM_DISTANCE=5000
DOWNSTREAM_DISTANCE=5000
BIN_SIZE=100
GENE_BODY_BINS=100

# Create output directory
mkdir -p ${OUTPUT_DIR}/temp

# Define gene lists and conditions
GENE_LISTS=("NEU_exo.txt" "not_NEU_exo.txt" "NSC_exo.txt" "not_NSC_exo.txt")
GENE_SET_NAMES=("NEU_target" "NEU_non_target" "NSC_target" "NSC_non_target")
CONDITIONS=("NSC_GFP" "NSC_M2" "Neu_GFP" "Neu_M2")

# Calculate which gene set and condition to process based on array task ID
# Each gene set has 4 conditions, so we have 16 total combinations
GENE_SET_INDEX=$((SLURM_ARRAY_TASK_ID / 4))
CONDITION_INDEX=$((SLURM_ARRAY_TASK_ID % 4))

GENE_LIST=${GENE_LISTS[$GENE_SET_INDEX]}
GENE_SET=${GENE_SET_NAMES[$GENE_SET_INDEX]}
CONDITION=${CONDITIONS[$CONDITION_INDEX]}

echo "Processing gene set: ${GENE_SET} (${GENE_LIST}) with condition: ${CONDITION}"

# Install required packages if needed
pip install scipy matplotlib==3.5.3 seaborn==0.12.2 pyBigWig pybedtools tqdm

# Run the parallel metaprofile computation
python 5a_compute_metaprofile_parallel.py \
    --results-dir ${RESULTS_DIR} \
    --output-dir ${OUTPUT_DIR}/temp \
    --genome-gtf ${GENOME_GTF} \
    --gene-list-file ${GENE_LISTS_DIR}/${GENE_LIST} \
    --gene-set-name ${GENE_SET} \
    --condition ${CONDITION} \
    --upstream-distance ${UPSTREAM_DISTANCE} \
    --downstream-distance ${DOWNSTREAM_DISTANCE} \
    --bin-size ${BIN_SIZE} \
    --gene-body-bins ${GENE_BODY_BINS}

echo "Metaprofile computation completed for gene set ${GENE_SET} and condition ${CONDITION}"

# If this is the last array job, trigger the plot generation
if [ $SLURM_ARRAY_TASK_ID -eq 15 ]; then
    echo "All metaprofile computations completed. Submitting job to create plots..."
    
    # Submit the job to create the final plots
    sbatch 5b_run_create_plots.sh
fi 