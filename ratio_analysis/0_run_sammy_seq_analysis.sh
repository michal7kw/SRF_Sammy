#!/bin/bash
#SBATCH --job-name=run_sammy_seq_analysis
#SBATCH --output=logs/run_sammy_seq_analysis.out
#SBATCH --error=logs/run_sammy_seq_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

"""
This code implements a complete SAMMY-seq analysis pipeline for identifying heterochromatin and euchromatin regions based on the comparison of S2S (accessible) and S3 (inaccessible) chromatin fractions. 

**Inputs:**
1. Working directory (--workdir): Contains processed BAM files in a "results/bam" subdirectory
2. Optional blacklist (--blacklist): BED file with regions to exclude from analysis
3. Configuration parameters:
   - Bin size (--bin-size): Genomic window size (default: 50,000bp)
   - Euchromatin threshold (--eu-threshold): Log2 ratio threshold (default: 0.1)
   - Heterochromatin threshold (--hetero-threshold): Log2 ratio threshold (default: -0.1)

**Outputs (saved to --outdir or workdir/chromatin_states_results/):**
1. Merged BAM files for each condition and fraction ("merged_bams/")
2. Coverage files in BigWig format ("coverage/")
3. Log2 ratio files comparing S2S to S3 ("ratio/")
4. BedGraph versions of ratio files ("bedgraph/")
5. BED files of identified euchromatin and heterochromatin regions ("chromatin_states/")
6. Summary TSV files:
   - chromatin_state_summary.tsv: Statistics on identified regions
   - condition_comparison_summary.tsv: Comparisons between conditions
7. Visualization plots ("plots/"):
   - Bar charts showing proportions of chromatin states
   - Pie charts comparing conditions
   - Chromosome heatmaps for visualization

The pipeline processes multiple conditions (NSC_GFP, NSC_M2, Neu_GFP, Neu_M2) and performs key comparisons to identify:
1. Effects of M2 treatment (GFP vs M2) in both NSC and Neu cells
2. Differences between cell types (NSC vs Neu) 

This script provides comprehensive chromatin state analysis based on SAMMY-seq methodology, where euchromatin regions are defined by log2(S2S/S3) > 0.1 and heterochromatin regions by log2(S2S/S3) < -0.1.
"""

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs
mkdir -p /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/shift_GFP_M2

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

# Run the Python script with all necessary arguments
python sammy_seq_analysis.py \
  --workdir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis \
  --outdir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/shift_GFP_M2 \
  --bin-size 50000 \
  --eu-threshold 0.1 \
  --hetero-threshold -0.1 \
  --threads 16 \