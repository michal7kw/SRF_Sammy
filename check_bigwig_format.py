#!/usr/bin/env python
import os
import sys
import pyBigWig
import glob

def inspect_bigwig_files(bigwig_dir):
    """Inspect bigWig files to determine chromosome naming convention."""
    print(f"Inspecting bigWig files in: {bigwig_dir}")
    
    # Find all bigWig files in the directory
    bigwig_files = glob.glob(os.path.join(bigwig_dir, "**/*.bw"), recursive=True)
    
    if not bigwig_files:
        print("No bigWig files found!")
        return
    
    print(f"Found {len(bigwig_files)} bigWig files")
    
    # Sample the first few files
    sample_size = min(3, len(bigwig_files))
    for i in range(sample_size):
        bw_file = bigwig_files[i]
        print(f"\nInspecting file: {os.path.basename(bw_file)}")
        
        try:
            bw = pyBigWig.open(bw_file)
            chroms = bw.chroms()
            
            print(f"Chromosome naming examples:")
            for j, (chrom, length) in enumerate(chroms.items()):
                print(f"  {chrom}: {length} bp")
                if j >= 4:  # Show only first 5 chromosomes
                    print(f"  ... and {len(chroms) - 5} more")
                    break
                    
            bw.close()
        except Exception as e:
            print(f"Error inspecting file: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        bigwig_dir = sys.argv[1]
    else:
        bigwig_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/results/bigwig"
    
    inspect_bigwig_files(bigwig_dir)