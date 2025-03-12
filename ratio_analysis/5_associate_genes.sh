#!/bin/bash
#SBATCH --job-name=5_associate_genes
#SBATCH --output=logs/5_associate_genes.out
#SBATCH --error=logs/5_associate_genes.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

work_dir="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"

# Create necessary directories
cd ${work_dir}/ratio_analysis

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

##########################################################################################################
python associate_regions_with_genes.py \
  --results-dir ${work_dir}/ratio_analysis/results/shift_GFP_M2_with_blacklist_with_mirko_gene_lists/shifted_regions \
  --gtf-file ${work_dir}/DATA/gencode.vM25.basic.annotation.gtf \
  --upstream 5000 \
  --downstream 5000