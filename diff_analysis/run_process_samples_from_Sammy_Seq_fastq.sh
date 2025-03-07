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
#SBATCH --array=1-34%10  # This will be updated based on the number of files found

# Error handling
set -e  # Exit on error
set -u  # Exit on undefined variable

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy"
cd $WORKDIR

# Directory setup
FASTQ_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/DATA/Sammy_Seq_fastq"
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/results"
GENOME_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/mm10_bowtie2_index"
THREADS=16  # Using all CPUs allocated to the job

# Check if the FASTQ directory exists
if [ ! -d "$FASTQ_DIR" ]; then
  echo "Error: FASTQ directory $FASTQ_DIR does not exist"
  exit 1
fi
echo "Using FASTQ directory: $FASTQ_DIR"

# Create output directories
mkdir -p $OUTPUT_DIR/bam
mkdir -p logs

# Check if the filtered_sample_list.txt exists
if [ ! -f "filtered_sample_list.txt" ]; then
  echo "Error: filtered_sample_list.txt not found"
  echo "Please run check_number_of_files.sh first to generate the file list"
  exit 1
fi

# Count how many files we have
FILE_COUNT=$(wc -l < filtered_sample_list.txt)
echo "Found $FILE_COUNT files in the filtered sample list"

# Display the first few files to verify
echo -e "\nFirst 10 files to process:"
head -n 10 filtered_sample_list.txt

# Check if the array size matches the number of files
if [[ "$SLURM_ARRAY_TASK_MAX" ]]; then
  ARRAY_SIZE=$((SLURM_ARRAY_TASK_MAX))
  if [ "$ARRAY_SIZE" != "$FILE_COUNT" ]; then
    echo "Warning: Array size ($ARRAY_SIZE) does not match the number of files ($FILE_COUNT)"
    echo "For optimal processing, you should cancel this job and resubmit with:"
    echo "sbatch --array=1-${FILE_COUNT}%10 $0"
  fi
else
  echo "Note: Running in non-array mode or array information not available"
fi

# Get the FASTQ file for this array job
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filtered_sample_list.txt)
if [ -z "$FASTQ" ]; then
  echo "Error: No file found for task ID ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi
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
    -x $GENOME_DIR/mm10 \
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

# Create output directory for bigWigs
mkdir -p $OUTPUT_DIR/bigwig

# Create bigWig file
echo "Creating bigWig for $BASE..."
bamCoverage \
    --bam $OUTPUT_DIR/bam/${BASE}.bam \
    --outFileName $OUTPUT_DIR/bigwig/${BASE}.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --effectiveGenomeSize 2652783500 \
    --ignoreForNormalization chrX chrY chrM \
    --minMappingQuality 30 \
    --extendReads 200 \
    --centerReads \
    --smoothLength 50 \
    --numberOfProcessors $THREADS

echo "Processing complete for $BASE"