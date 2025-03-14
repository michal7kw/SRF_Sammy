#!/bin/bash
#SBATCH --job-name=4_plots_shift_NSC_NEU
#SBATCH --output=logs/4_plots_shift_NSC_NEU.out
#SBATCH --error=logs/4_plots_shift_NSC_NEU.err
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
mkdir -p ${work_dir}/ratio_analysis/results/shift_NSC_NEU
mkdir -p ${work_dir}/ratio_analysis/results/shift_NSC_NEU/bedgraph

# Check if input files exist and copy them if needed
if [ ! -f "${work_dir}/ratio_analysis/results/shift_NSC_NEU/bedgraph/NSC_GFP_S2S_vs_S3_ratio.bedgraph" ]; then
    echo "Copying NSC_GFP bedgraph file..."
    # You may need to adjust the source path based on where your original files are located
    cp ${work_dir}/ratio_analysis/results/shift_GFP_M2/bedgraph/NSC_GFP_S2S_vs_S3_ratio.bedgraph \
       ${work_dir}/ratio_analysis/results/shift_NSC_NEU/bedgraph/ || echo "Warning: Could not copy NSC_GFP file"
fi

if [ ! -f "${work_dir}/ratio_analysis/results/shift_NSC_NEU/bedgraph/Neu_GFP_S2S_vs_S3_ratio.bedgraph" ]; then
    echo "Copying Neu_GFP bedgraph file..."
    # You may need to adjust the source path based on where your original files are located
    cp ${work_dir}/ratio_analysis/results/shift_GFP_M2/bedgraph/Neu_GFP_S2S_vs_S3_ratio.bedgraph \
       ${work_dir}/ratio_analysis/results/shift_NSC_NEU/bedgraph/ || echo "Warning: Could not copy Neu_GFP file"
fi

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
output_file="${results_dir}/NSC_Neu/chromatin_changes_150kb.png"

echo "Running analysis with results directory: ${results_dir}"
echo "Output will be saved to: ${output_file}"

# Make sure the output directory exists
mkdir -p $(dirname "${output_file}")

python 4b_plots_shift_NSC_NEU_loosen_with_window.py \
  --results-dir ${results_dir} \
  --ratio-change-threshold 0.1 \
  --window-size 150000 \
  --debug \
  --genome-size 2652783500 \
  --output-file ${output_file}

echo "Analysis completed"

# Run with blacklist and gene lists if needed
# results_dir="${work_dir}/ratio_analysis/results/shift_NSC_NEU_with_blacklist_with_gene_lists"
# output_file="${results_dir}/chromatin_changes_150kb.png"

# python 4b_plots_shift_NSC_NEU_loosen_with_window.py \
#   --results-dir ${results_dir} \
#   --ratio-change-threshold 0.1 \
#   --window-size 150000 \
#   --debug \
#   --genome-size 2652783500 \
#   --output-file ${output_file}
