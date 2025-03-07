#!/bin/bash
#SBATCH --job-name=build_mm10_index
#SBATCH --output=logs/build_mm10_index.out
#SBATCH --error=logs/build_mm10_index.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis 

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

# Create directory for mm10 genome
mkdir -p logs
mkdir -p mm10
cd mm10

# Download mm10 genome from UCSC
echo "Downloading mm10 genome from UCSC..."
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz

# Decompress the genome
echo "Decompressing mm10 genome..."
gunzip mm10.fa.gz

# Create BWA index
echo "Creating BWA index for mm10 genome..."
bwa index -a bwtsw mm10.fa

echo "BWA index creation complete. The following index files have been created:"
ls -lh mm10.fa.*

echo "BWA index files for mm10 are ready to use."