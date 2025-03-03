#!/bin/bash

# Navigate to your working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy

# Define the correct FASTQ directory
FASTQ_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/DATA/Sammy_Seq_fastq"

# Check if the directory exists
if [ ! -d "$FASTQ_DIR" ]; then
  echo "Error: Directory $FASTQ_DIR does not exist"
  echo "Please check the path and try again"
  exit 1
fi

echo "Using FASTQ directory: $FASTQ_DIR"

# Create a temporary file to store the list
> temp_filtered_list.txt

# Find all files matching the pattern
for cell_type in Neu NSC; do
  for condition in GFP M2; do
    for state in S2S S3; do
      echo "Checking pattern: ${cell_type}*_${condition}_${state}_*4F*_*.fastq.gz"
      find "$FASTQ_DIR" -name "${cell_type}*_${condition}_${state}_*4F*_*.fastq.gz" -type f >> temp_filtered_list.txt
    done
  done
done

# Count and display the files
FILE_COUNT=$(wc -l < temp_filtered_list.txt)
echo "Found $FILE_COUNT files matching the pattern"

# Display the first few files to verify
echo -e "\nFirst 10 matching files:"
head -n 10 temp_filtered_list.txt

# Display the total count for each combination
echo -e "\nCount by combination:"
for cell_type in Neu NSC; do
  for condition in GFP M2; do
    for state in S2S S3; do
      count=$(grep -c "${cell_type}.*_${condition}_${state}.*4F" temp_filtered_list.txt || echo 0)
      echo "${cell_type}_${condition}_${state}: $count files"
    done
  done
done

# Create the final filtered list for processing
echo -e "\nCreating filtered_sample_list.txt for processing script..."
cat temp_filtered_list.txt | sed "s|$FASTQ_DIR/||" > filtered_sample_list.txt
echo "Saved $FILE_COUNT files to filtered_sample_list.txt"
echo "To process these files, run: sbatch --array=1-${FILE_COUNT}%10 run_process_samples_from_Sammy_Seq_fastq.sh"

# Show some examples of files that don't match the pattern
echo -e "\nSample of files in directory that don't match the pattern:"
find "$FASTQ_DIR" -name "*.fastq.gz" | grep -v -E "(Neu|NSC).*_(GFP|M2)_(S2S|S3).*4F" | head -n 10

# Check for files with 4F but different naming pattern
echo -e "\nFiles containing '4F' but not matching our pattern:"
find "$FASTQ_DIR" -name "*4F_*.fastq.gz" | grep -v -E "(Neu|NSC).*_(GFP|M2)_(S2S|S3).*4F" | head -n 10

# Suggest alternative patterns if needed
echo -e "\nPossible alternative patterns to consider:"
echo "1. Files with 4F as prefix: $(find "$FASTQ_DIR" -name "4F_*.fastq.gz" | wc -l) files"
echo "2. Files with S2S but different format: $(find "$FASTQ_DIR" -name "*S2S_*.fastq.gz" | grep -v -E "(Neu|NSC).*_(GFP|M2)_S2S" | wc -l) files"
echo "3. Files with S3 but different format: $(find "$FASTQ_DIR" -name "*S3_*.fastq.gz" | grep -v -E "(Neu|NSC).*_(GFP|M2)_S3" | wc -l) files"