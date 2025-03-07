#!/bin/bash
#SBATCH --job-name=run_gene_chromatin_state_analysis
#SBATCH --output=logs/run_gene_chromatin_state_analysis.out
#SBATCH --error=logs/run_gene_chromatin_state_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs
mkdir -p /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/gene_chromatin_state_analysis

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

# Run the Python script with all necessary arguments
python gene_chromatin_state_analysis.py \
  --gene-list /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/gene_lists/endo_enriched_gene_list_nsc_vs_neu.txt \
  --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/py_implementation \
  --output-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/gene_promoter_chromatin_state_analysis \
  --genome-gtf /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/DATA/gencode.vM25.basic.annotation.gtf \
  --promoter-size 2000