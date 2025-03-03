import csv
import os
import sys
from pathlib import Path

def set_working_directory():
    """Set the working directory to the specified path."""
    work_dir = Path("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/gene_lists")
    
    if not work_dir.exists():
        print(f"Error: Working directory does not exist: {work_dir}")
        sys.exit(1)
    
    try:
        os.chdir(work_dir)
        print(f"Working directory set to: {work_dir}")
    except Exception as e:
        print(f"Error setting working directory: {e}")
        sys.exit(1)

def process_file(input_file, output_file):
    """Process a CSV file and extract SYMBOL column to a text file."""
    print(f"Processing {input_file} -> {output_file}")
    try:
        # Read the CSV file and extract SYMBOL column
        with open(input_file, 'r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            symbols = [row['SYMBOL'] for row in csv_reader]

        # Write symbols to text file
        with open(output_file, 'w') as txt_file:
            txt_file.write('\n'.join(symbols))
        
        print(f"Successfully processed {len(symbols)} symbols")
    except Exception as e:
        print(f"Error processing {input_file}: {e}")

# Set the working directory at the beginning of script execution
set_working_directory()

# Define all input-output file pairs
file_pairs = [
    ("endo_list_neu_tss_peaks.csv", "endo_gene_list_neu.txt"),
    ("endo_list_nsc_tss_peaks.csv", "endo_gene_list_nsc.txt"),
    ("neu_nsc_enriched_signal_2_endo.csv", "endo_enriched_gene_list_nsc_vs_neu.txt"),
    ("exo_endo_enriched_signal_2_neu.csv", "neu_enriched_gene_list_exo_vs_endo.txt"),
    ("exo_endo_enriched_signal_2_nsc.csv", "nsc_enriched_gene_list_exo_vs_endo.txt")
]

# Process each file pair
for input_file, output_file in file_pairs:
    process_file(input_file, output_file)