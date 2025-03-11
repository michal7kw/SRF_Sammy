#!/bin/bash
#SBATCH --job-name=2_run_gene_chromatin_state_comparison_nsc
#SBATCH --output=logs/2_run_gene_chromatin_state_comparison_nsc.out
#SBATCH --error=logs/2_run_gene_chromatin_state_comparison_nsc.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq


"""
This code analyzes chromatin states for specific genes, focusing solely on the comparison between NSC (Neural Stem Cell) and Neu (Neuron) samples in GFP condition. 
It creates visualizations of how chromatin states differ between these cell types.

**Inputs:**
1. Gene list file (--gene-list): A text file containing gene names, one per line
2. Results directory (--results-dir): Contains chromatin state BED files in a "chromatin_states" subdirectory
3. Genome GTF file (--genome-gtf): Gene annotation file for determining promoter coordinates

**Outputs (saved to --output-dir):**
1. gene_chromatin_states.csv: Table with chromatin state data for each gene
2. Visualization plots:
   - chromatin_states_distribution.png: Pie charts showing proportion of genes in each state
   - chromatin_state_transitions.png: Bar chart showing state changes from NSC to Neu
   - key_chromatin_state_changes.png: Focused bar chart on key transitions
   - chromatin_state_change_proportion.png: Pie chart showing percentage of genes that changed state

Unlike the previous code, this version:
1. Only performs NSC vs Neu comparison (no M2 condition analysis)
2. Uses a simpler approach to determine chromatin states
3. Does not include the numerical signal strength analysis
4. Has a more straightforward file path structure for loading chromatin states
5. Creates fewer visualization plots

The script classifies each gene promoter as being in euchromatin (open), heterochromatin (closed), or unknown state based on overlap with chromatin regions from BED files.
"""

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

base_dir="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"

# # Set variables
# GENE_LIST="${base_dir}/gene_lists/nsc_enriched_gene_list_exo_vs_endo.txt"
# RESULTS_DIR="${base_dir}/ratio_analysis/results"
# OUTPUT_DIR="${base_dir}/ratio_analysis/results/nsc_enriched_comparison"
# GENOME_GTF="${base_dir}/DATA/gencode.vM25.basic.annotation.gtf"
# PROMOTER_SIZE=2000

# mkdir -p ${OUTPUT_DIR}

# # Run the analysis
# python 2_gene_chromatin_state_comparison.py \
#     --gene-list ${GENE_LIST} \
#     --results-dir ${RESULTS_DIR} \
#     --output-dir ${OUTPUT_DIR} \
#     --genome-gtf ${GENOME_GTF} \
#     --promoter-size ${PROMOTER_SIZE}

##########################################################################################################

# Set variables
GENE_LIST="${base_dir}/figures/genelists/NSC_exo.txt"
RESULTS_DIR="${base_dir}/ratio_analysis/results/shift_GFP_M2_with_blacklist_with_mirko_gene_lists"
OUTPUT_DIR="${base_dir}/ratio_analysis/results/nsc_enriched_comparison_mirko_gene_lists"
GENOME_GTF="${base_dir}/DATA/gencode.vM25.basic.annotation.gtf"
PROMOTER_SIZE=2000

# mkdir -p ${OUTPUT_DIR}

# Run the analysis
python 2_gene_chromatin_state_comparison.py \
    --gene-list ${GENE_LIST} \
    --results-dir ${RESULTS_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --genome-gtf ${GENOME_GTF} \
    --promoter-size ${PROMOTER_SIZE}

echo "Analysis completed."
