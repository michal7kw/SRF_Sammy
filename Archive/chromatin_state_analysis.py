#!/usr/bin/env python3

import os
import glob
import pandas as pd
import subprocess
from pathlib import Path

class ChromatinStateAnalysis:
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        self.output_dir = self.data_dir.parent / "results"
        self.output_dir.mkdir(exist_ok=True)
        
    def get_sample_files(self, condition, cell_type, state):
        """Get FASTQ files for a specific condition, cell type and chromatin state."""
        pattern = f"{cell_type}*_{condition}_{state}_*_R1_001.fastq.gz"
        return sorted(glob.glob(str(self.data_dir / pattern)))
    
    def create_bigwig(self, bam_file, output_file):
        """Convert BAM to bigWig format."""
        # Using deepTools bamCoverage for bigWig conversion
        cmd = [
            "bamCoverage",
            "-b", bam_file,
            "-o", output_file,
            "--binSize", "10",
            "--normalizeUsing", "RPKM",
            "--smoothLength", "50",
            "-p", "8"  # number of processors
        ]
        subprocess.run(cmd, check=True)
    
    def analyze_state_changes(self):
        """Analyze chromatin state changes between S2S and S3."""
        cell_types = ["Neu", "NSC"]  # Neurons and NPCs
        conditions = ["GFP", "M2"]    # Control and Mecp2
        states = ["S2S", "S3"]        # Euchromatin and heterochromatin
        
        for cell_type in cell_types:
            print(f"\nProcessing {cell_type}...")
            
            # Create output directories for this cell type
            cell_output_dir = self.output_dir / cell_type
            cell_output_dir.mkdir(exist_ok=True)
            
            for condition in conditions:
                print(f"  Analyzing {condition}...")
                
                for state in states:
                    files = self.get_sample_files(condition, cell_type, state)
                    print(f"    Found {len(files)} files for {state}")
                    
                    # Process each replicate
                    for i, fastq in enumerate(files, 1):
                        base_name = f"{cell_type}_{condition}_{state}_rep{i}"
                        print(f"    Processing {base_name}")
                        
                        # Output bigWig file path
                        bw_file = cell_output_dir / f"{base_name}.bw"
                        
                        # Here:
                        # 1. Align FASTQ to genome (e.g., using bowtie2)
                        # 2. Convert SAM to BAM
                        # 3. Sort and index BAM
                        # 4. Create bigWig
                        print(f"      Would create bigWig: {bw_file}")

def main():
    # Initialize analysis
    data_dir = Path("/beegfs/scratch/ric.broccoli/kubacki.michal/Sammy/DATA/Sammy_Seq_fastq")
    analysis = ChromatinStateAnalysis(data_dir)
    
    # Run analysis
    analysis.analyze_state_changes()

if __name__ == "__main__":
    main()
