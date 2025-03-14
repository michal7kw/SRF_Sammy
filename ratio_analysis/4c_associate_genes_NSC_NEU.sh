#!/bin/bash
#SBATCH --job-name=4c_associate_promoters_NSC_NEU
#SBATCH --output=logs/4c_associate_promoters_NSC_NEU.out
#SBATCH --error=logs/4c_associate_promoters_NSC_NEU.err
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
mkdir -p logs

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

##########################################################################################################
# Run the gene association script for NSC_NEU comparison, focusing only on promoters
python 4c_associate_regions_with_genes_NSC_NEU.py \
  --results-dir ${work_dir}/ratio_analysis/results/shift_GFP_M2_with_blacklist_with_mirko_gene_lists/shifted_regions_NSC_NEU \
  --gtf-file ${work_dir}/DATA/gencode.vM25.basic.annotation.gtf \
  --promoter-size 5000 \
  --debug

echo "Promoter association analysis completed for NSC_GFP vs Neu_GFP comparison"

# Optional: Run a second analysis that includes gene bodies
# python 4c_associate_regions_with_genes_NSC_NEU.py \
#   --results-dir ${work_dir}/ratio_analysis/results/shift_NSC_NEU/shifted_regions_NSC_NEU \
#   --gtf-file ${work_dir}/DATA/gencode.vM25.basic.annotation.gtf \
#   --promoter-size 5000 \
#   --include-gene-body \
#   --debug