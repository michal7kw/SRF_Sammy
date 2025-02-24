#!/bin/bash
#SBATCH --job-name=process_samples
#SBATCH --output=logs/process_samples_%a.out
#SBATCH --error=logs/process_samples_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --array=1-63

# Error handling
set -e  # Exit on error
set -u  # Exit on undefined variable

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Sammy"
cd $WORKDIR

# Directory setup
FASTQ_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Sammy/DATA/Sammy_Seq_fastq"
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/Sammy/results"
GENOME_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome"
THREADS=16  # Using all CPUs allocated to the job

# Create output directories
mkdir -p $OUTPUT_DIR/bam
mkdir -p $OUTPUT_DIR/bigwig

# Get the FASTQ file for this array job
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_list.txt)
FASTQ="$FASTQ_DIR/$FASTQ"

if [ ! -f "$FASTQ" ]; then
    echo "Error: FASTQ file $FASTQ not found"
    exit 1
fi

# Get base name for output files
BASE=$(basename $FASTQ _R1_001.fastq.gz)
echo "Processing $BASE..."

# Align with bowtie2
echo "Aligning with bowtie2..."
bowtie2 -p $THREADS \
    -x $GENOME_DIR/GRCh38 \
    -U $FASTQ \
    -S $OUTPUT_DIR/bam/${BASE}.sam

# Convert SAM to BAM
echo "Converting SAM to BAM..."
samtools view -bS $OUTPUT_DIR/bam/${BASE}.sam > $OUTPUT_DIR/bam/${BASE}.bam
rm $OUTPUT_DIR/bam/${BASE}.sam

# Sort BAM
echo "Sorting BAM..."
samtools sort -@ $THREADS $OUTPUT_DIR/bam/${BASE}.bam -o $OUTPUT_DIR/bam/${BASE}.sorted.bam
mv $OUTPUT_DIR/bam/${BASE}.sorted.bam $OUTPUT_DIR/bam/${BASE}.bam

# Index BAM
echo "Indexing BAM..."
samtools index $OUTPUT_DIR/bam/${BASE}.bam

# Create bigWig
echo "Creating bigWig..."
bamCoverage -b $OUTPUT_DIR/bam/${BASE}.bam \
    -o $OUTPUT_DIR/bigwig/${BASE}.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --smoothLength 50 \
    -p $THREADS

echo "Processing complete for $BASE"

echo "Processing complete. BigWig files are in $OUTPUT_DIR/bigwig"
