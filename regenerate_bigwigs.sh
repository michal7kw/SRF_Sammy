#!/bin/bash
#SBATCH --job-name=regenerate_bigwigs
#SBATCH --output=logs/regenerate_bigwigs_%a.out
#SBATCH --error=logs/regenerate_bigwigs_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --array=1-63

# Error handling
set -e
set -u

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define directories
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Sammy"
OUTPUT_DIR="$WORKDIR/results"
GENOME_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome"
THREADS=16

cd $WORKDIR

# Create output directory for new bigWigs
mkdir -p $OUTPUT_DIR/bigwig_new

# Get the BAM file for this array job
BAM_FILE=$(ls $OUTPUT_DIR/bam/*.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file $BAM_FILE not found"
    exit 1
fi

# Get base name for output files
BASE=$(basename $BAM_FILE .bam)
echo "Processing $BASE..."

# Create bigWig with improved parameters
echo "Creating bigWig for $BASE..."
bamCoverage \
    --bam $BAM_FILE \
    --outFileName $OUTPUT_DIR/bigwig_new/${BASE}.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --effectiveGenomeSize 2913022398 \
    --ignoreForNormalization chrX chrY chrM \
    --minMappingQuality 30 \
    --extendReads \
    --centerReads \
    --smoothLength 50 \
    --numberOfProcessors $THREADS

echo "Completed processing $BASE"
