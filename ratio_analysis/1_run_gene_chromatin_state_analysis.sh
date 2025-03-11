#!/bin/bash
#SBATCH --job-name=1_run_gene_chromatin_state_analysis
#SBATCH --output=logs/1_run_gene_chromatin_state_analysis.out
#SBATCH --error=logs/1_run_gene_chromatin_state_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

"""
This code performs chromatin state analysis for specific genes in NSC (Neural Stem Cell) and Neu (Neuron) samples under GFP (or M2 conditions). 
It examines how chromatin states (euchromatin and heterochromatin) differ across these conditions, focusing on gene promoter regions.

**Inputs:**
1. Gene list file (--gene-list): A text file with gene names, one per line
2. Results directory (--results-dir): Contains chromatin state BED files for different conditions
3. Genome GTF file (--genome-gtf): Gene annotation file used to get promoter coordinates

**Outputs (saved to --output-dir):**
1. CSV file with all chromatin state data
2. Visualization plots in two subdirectories (gfp_vs_m2/ and nsc_vs_neu/):
   - Distribution pie charts showing proportion of euchromatin vs heterochromatin
   - Transition bar charts showing how states change between conditions
   - Boxplots of signal strengths
   - Scatter plots comparing signals between conditions
   - Histograms of euchromatin/heterochromatin ratios

The script analyzes two key comparisons:
1. GFP vs M2 in NSCs (to study effect of M2 treatment)
2. NSC vs Neu in GFP condition (to study differentiation effects)

For each gene, it determines if the promoter is in euchromatin (open) or heterochromatin (closed) state and tracks changes across conditions.
"""

work_dir="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"

# Create necessary directories
cd ${work_dir}/ratio_analysis
mkdir -p logs

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

##########################################################################################################

# genelist_path="${work_dir}/gene_lists/"
# genelist_file="endo_enriched_gene_list_nsc_vs_neu.txt"
# genelist_path="${genelist_path}${genelist_file}"

# results_dir="${work_dir}/ratio_analysis/results/shift_GFP_M2"
# output_dir="${work_dir}/ratio_analysis/results/gene_promoter_chromatin_state_analysis"
# genome_gtf="${work_dir}/DATA/gencode.vM25.basic.annotation.gtf"
# promoter_size=2000

# python 1_gene_chromatin_state_analysis.py \
#   --gene-list ${genelist_path} \
#   --results-dir ${results_dir} \
#   --output-dir ${output_dir} \
#   --genome-gtf ${genome_gtf} \
#   --promoter-size ${promoter_size}

##########################################################################################################

genelist_path="${work_dir}/figures/genelists/"
genelist_file="endo.txt"
genelist_path="${genelist_path}${genelist_file}"

results_dir="${work_dir}/ratio_analysis/results/shift_GFP_M2_with_blacklist_with_mirko_gene_lists"
output_dir="${work_dir}/ratio_analysis/results/gene_promoter_chromatin_state_analysis_mirko_gene_lists"
genome_gtf="${work_dir}/DATA/gencode.vM25.basic.annotation.gtf"
promoter_size=2000

mkdir -p ${output_dir}

python 1_gene_chromatin_state_analysis.py \
  --gene-list ${genelist_path} \
  --results-dir ${results_dir} \
  --output-dir ${output_dir} \
  --genome-gtf ${genome_gtf} \
  --promoter-size ${promoter_size}