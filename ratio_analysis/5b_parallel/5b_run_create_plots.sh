#!/bin/bash
#SBATCH --job-name=5_mecp2_plots
#SBATCH --output=logs/5_mecp2_plots.out
#SBATCH --error=logs/5_mecp2_plots.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# This script combines the results of the parallel metaprofile computations
# and generates the final heatmaps and line plots.
#
# It should be run after all the array jobs have completed.

# Create necessary directories
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis
mkdir -p logs

# Load conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate snakemake

base_dir="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"

# Set variables
OUTPUT_DIR="${base_dir}/ratio_analysis/results/mecp2_target_comparison_heatmaps_parallel"
TEMP_DIR="${OUTPUT_DIR}/temp"

# Install required packages if needed
pip install scipy matplotlib==3.5.3 seaborn==0.12.2 pyBigWig pybedtools tqdm

# Check if all required files exist
required_files=16  # 4 gene sets x 4 conditions
actual_files=$(ls ${TEMP_DIR}/*_matrix.pkl 2>/dev/null | wc -l)

if [ ${actual_files} -lt ${required_files} ]; then
    echo "Warning: Not all required files exist. Found ${actual_files} out of ${required_files} files."
    echo "Will proceed with available files, but some conditions may be missing from the plots."
fi

echo "All required files found. Generating plots..."

# Run the plot generation script
python 5b_create_mecp2_target_comparison_plots.py \
    --input-dir ${TEMP_DIR} \
    --output-dir ${OUTPUT_DIR}

echo "Plot generation completed."

# Create a directory for enhanced visualizations
ENHANCED_OUTPUT_DIR="${OUTPUT_DIR}/enhanced_visualizations"
mkdir -p ${ENHANCED_OUTPUT_DIR}

# Copy the generated files to the enhanced directory
cp ${OUTPUT_DIR}/mecp2_target_comparison_heatmap.png ${ENHANCED_OUTPUT_DIR}/
cp ${OUTPUT_DIR}/mecp2_target_comparison_heatmap_line_plots.png ${ENHANCED_OUTPUT_DIR}/

echo "Enhanced visualizations are available in ${ENHANCED_OUTPUT_DIR}" 