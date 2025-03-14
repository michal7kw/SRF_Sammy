#!/bin/bash
#SBATCH --job-name=4_plots_GFP_M2
#SBATCH --output=logs/4_plots_GFP_M2.out
#SBATCH --error=logs/4_plots_GFP_M2.err
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
mkdir -p ${work_dir}/ratio_analysis/results/shift_GFP_M2

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

##########################################################################################################

# results_dir="${work_dir}/ratio_analysis/results/shift_GFP_M2"
# output_file="${results_dir}/chromatin_transitions.png"

# python plots_GFP_M2.py \
#   --results-dir ${results_dir} \
#   --output-file ${output_file}


##########################################################################################################

results_dir="${work_dir}/ratio_analysis/results/shift_GFP_M2_with_blacklist_with_mirko_gene_lists"
output_file="${results_dir}/chromatin_transitions.png"

python 4_plots_GFP_M2.py \
  --results-dir ${results_dir} \
  --output-file ${output_file}