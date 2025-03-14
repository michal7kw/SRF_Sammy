#!/bin/bash
#SBATCH --job-name=4_plots_shift_GFP_M2
#SBATCH --output=logs/4_plots_shift_GFP_M2.out
#SBATCH --error=logs/4_plots_shift_GFP_M2.err
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

# python 4_plots_shift_GPF_M2_binary.py \
#   --genome-size 2652783500 \
#   --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/shift_GFP_M2 \
#   --output-file /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/shift_GFP_M2/chromatin_changes.png

# python 4_plots_shift_GPF_M2_normal.py \
#   --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/shift_GFP_M2 \
#   --method ratio_comparison \
#   --debug \
#   --genome-size 2652783500 \
#   --output-file /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/shift_GFP_M2/chromatin_changes.png

# python 4_plots_shift_GPF_M2_loosen.py \
#   --results-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/shift_GFP_M2 \
#   --ratio-change-threshold 0.2 \
#   --debug \
#   --genome-size 2652783500 \
#   --output-file /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/results/shift_GFP_M2/chromatin_changes.png
  
##########################################################################################################
# results_dir="${work_dir}/ratio_analysis/results/shift_GFP_M2"
# output_file="${results_dir}/chromatin_changes_150kb.png"

# python 4_plots_shift_GPF_M2_loosen_with_window.py \
#   --results-dir ${results_dir} \
#   --ratio-change-threshold 0.1 \
#   --window-size 150000 \
#   --debug \
#   --genome-size 2652783500 \
#   --output-file ${output_file}

##########################################################################################################
results_dir="${work_dir}/ratio_analysis/results/shift_GFP_M2_with_blacklist_with_mirko_gene_lists"
output_file="${results_dir}/chromatin_changes_150kb.png"

python 4_plots_shift_GPF_M2_loosen_with_window.py \
  --results-dir ${results_dir} \
  --ratio-change-threshold 0.1 \
  --window-size 150000 \
  --debug \
  --genome-size 2652783500 \
  --output-file ${output_file}
