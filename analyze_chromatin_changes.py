#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
from pathlib import Path
import pyBigWig
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import sys
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('chromatin_analysis.log')
    ]
)

class ChromatinStateAnalyzer:
    def __init__(self, results_dir):
        self.results_dir = Path(results_dir)
        self.bigwig_dir = self.results_dir / "bigwig"
        self.output_dir = self.results_dir / "analysis"
        
        # Verify directories exist
        if not self.results_dir.exists():
            raise FileNotFoundError(f"Results directory not found: {self.results_dir}")
        if not self.bigwig_dir.exists():
            raise FileNotFoundError(f"BigWig directory not found: {self.bigwig_dir}")
            
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(exist_ok=True)
        logging.info(f"Using results directory: {self.results_dir}")
        logging.info(f"Using bigwig directory: {self.bigwig_dir}")
        logging.info(f"Output will be saved to: {self.output_dir}")
        
        # Parameters for analysis
        self.window_size = 1000  # Size of windows to analyze
        self.min_fold_change = 2  # Minimum fold change to consider significant
        self.pvalue_threshold = 0.05
        
        logging.info(f"Analysis parameters: window_size={self.window_size}, min_fold_change={self.min_fold_change}, pvalue_threshold={self.pvalue_threshold}")
        
    def get_bigwig_files(self, cell_type, condition, state):
        """Get bigWig files for a specific condition."""
        pattern = f"{cell_type}*_{condition}_{state}_*.bw"
        files = sorted(glob.glob(str(self.bigwig_dir / pattern)))
        if not files:
            logging.warning(f"No bigWig files found for pattern: {pattern}")
        else:
            logging.info(f"Found {len(files)} bigWig files for {cell_type} {condition} {state}")
            for f in files:
                logging.debug(f"Found file: {f}")
        return files
    
    def compute_average_signal(self, bigwig_files, chrom, start, end):
        """Compute average signal across replicates for a region."""
        signals = []
        for bw_file in bigwig_files:
            with pyBigWig.open(bw_file) as bw:
                try:
                    values = bw.values(chrom, start, end)
                    if values is not None:
                        # Replace negative values with 0
                        values = np.array(values)
                        values[values < 0] = 0
                        signal = np.nanmean(values)
                        if not np.isnan(signal) and signal >= 0:
                            signals.append(signal)
                except Exception as e:
                    print(f"Warning: Error processing {bw_file}: {str(e)}")
                    continue
        return np.nanmean(signals) if signals else 0.0  # Return 0 instead of nan for regions with no signal
    
    def identify_state_changes(self, cell_type):
        """Identify regions that change state between S2S and S3."""
        print(f"\nAnalyzing {cell_type}...")
        
        # Get chromosome sizes from one of the bigWig files
        sample_bw = next(self.bigwig_dir.glob("*.bw"))
        with pyBigWig.open(str(sample_bw)) as bw:
            chroms = [(chrom, length) for chrom, length in bw.chroms().items()]
        
        results = []
        
        for chrom, length in chroms:
            print(f"Processing {chrom}...")
            
            # Analyze windows across the chromosome
            for start in range(0, length, self.window_size):
                end = min(start + self.window_size, length)
                
                # Get signals for each condition
                gfp_s2s = self.get_bigwig_files(cell_type, "GFP", "S2S")
                gfp_s3 = self.get_bigwig_files(cell_type, "GFP", "S3")
                m2_s2s = self.get_bigwig_files(cell_type, "M2", "S2S")
                m2_s3 = self.get_bigwig_files(cell_type, "M2", "S3")
                
                # Compute average signals
                gfp_s2s_signal = self.compute_average_signal(gfp_s2s, chrom, start, end)
                gfp_s3_signal = self.compute_average_signal(gfp_s3, chrom, start, end)
                m2_s2s_signal = self.compute_average_signal(m2_s2s, chrom, start, end)
                m2_s3_signal = self.compute_average_signal(m2_s3, chrom, start, end)
                
                # Add small epsilon to avoid division by zero
                epsilon = 1e-10
                
                # Calculate state changes with safety checks
                if gfp_s2s_signal > epsilon and m2_s2s_signal > epsilon:
                    gfp_ratio = gfp_s3_signal / gfp_s2s_signal
                    m2_ratio = m2_s3_signal / m2_s2s_signal
                    
                    # Calculate log2 fold change safely
                    if gfp_ratio > 0 and m2_ratio > 0:
                        log2fc = np.log2(m2_ratio/gfp_ratio)
                        
                        # Identify significant changes
                        if abs(log2fc) >= np.log2(self.min_fold_change):
                            results.append({
                                'chrom': chrom,
                                'start': start,
                                'end': end,
                                'GFP_S2S': gfp_s2s_signal,
                                'GFP_S3': gfp_s3_signal,
                                'M2_S2S': m2_s2s_signal,
                                'M2_S3': m2_s3_signal,
                                'GFP_ratio': gfp_ratio,
                                'M2_ratio': m2_ratio,
                                'log2FC': log2fc,
                                'direction': 'S2S→S3' if log2fc > 0 else 'S3→S2S'
                            })
        
        # Convert results to DataFrame
        df = pd.DataFrame(results)
        
        # Save results
        output_file = self.output_dir / f"{cell_type}_chromatin_changes.csv"
        df.to_csv(output_file, index=False)
        
        # Create summary plots
        self.plot_results(df, cell_type)
        
        return df
    
    def plot_results(self, df, cell_type):
        """Create visualization of chromatin state changes."""
        if len(df) == 0:
            print(f"\nNo significant changes found for {cell_type}")
            return
            
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot 1: Distribution of log2 fold changes
        sns.histplot(data=df, x='log2FC', bins=50, ax=ax1)
        ax1.axvline(x=0, color='r', linestyle='--')
        ax1.set_title(f'Distribution of Chromatin State Changes\nin {cell_type}')
        ax1.set_xlabel('log2(M2_ratio/GFP_ratio)')
        ax1.set_ylabel('Count')
        
        # Plot 2: Direction of changes
        direction_counts = df['direction'].value_counts()
        sns.barplot(x=direction_counts.index, y=direction_counts.values, ax=ax2)
        ax2.set_title(f'Direction of Chromatin State Changes\nin {cell_type}')
        ax2.set_ylabel('Number of regions')
        plt.xticks(rotation=45)
        
        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{cell_type}_chromatin_changes_dist.pdf")
        plt.close()
        
        # Print summary statistics
        total_regions = len(df)
        s2s_to_s3 = sum(df['direction'] == 'S2S→S3')
        s3_to_s2s = sum(df['direction'] == 'S3→S2S')
        
        print(f"\nSummary for {cell_type}:")
        print(f"Total regions with significant changes: {total_regions}")
        print(f"Regions changing S2S → S3: {s2s_to_s3} ({s2s_to_s3/total_regions*100:.1f}%)")
        print(f"Regions changing S3 → S2S: {s3_to_s2s} ({s3_to_s2s/total_regions*100:.1f}%)")
        
        # Calculate and print average fold changes
        s2s_to_s3_fc = df[df['direction'] == 'S2S→S3']['log2FC'].mean()
        s3_to_s2s_fc = df[df['direction'] == 'S3→S2S']['log2FC'].mean()
        print(f"Average log2 fold change for S2S → S3: {s2s_to_s3_fc:.2f}")
        print(f"Average log2 fold change for S3 → S2S: {s3_to_s2s_fc:.2f}")

def main():
    try:
        logging.info("Starting chromatin state analysis")
        start_time = datetime.now()
        
        # Initialize analysis
        results_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/Sammy/results"
        logging.info(f"Initializing ChromatinStateAnalyzer with results_dir: {results_dir}")
        analyzer = ChromatinStateAnalyzer(results_dir)
        
        # Analyze both cell types
        cell_types = ["Neu", "NSC"]
        for cell_type in cell_types:
            logging.info(f"\n{'='*50}\nProcessing cell type: {cell_type}\n{'='*50}")
            analyzer.identify_state_changes(cell_type)
            
        end_time = datetime.now()
        duration = end_time - start_time
        logging.info(f"\nAnalysis completed successfully in {duration}\n")
        
    except Exception as e:
        logging.error(f"An error occurred during analysis: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
