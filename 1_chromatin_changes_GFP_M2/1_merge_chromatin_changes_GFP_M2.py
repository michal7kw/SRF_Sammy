#!/usr/bin/env python3
import os
import glob
import argparse
import logging
import pandas as pd
from pathlib import Path

def merge_results(results_dir='results'):
    """
    Merge chromosome-specific results into genome-wide results.
    
    Args:
        results_dir (str): Path to the results directory
    
    Returns:
        bool: True if merge was successful, False otherwise
    """
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    logging.info("Starting merge operation")
    
    # Convert to Path object
    results_dir = Path(results_dir)
    analysis_dir = results_dir / "analysis"
    
    # Check if analysis directory exists
    if not analysis_dir.exists():
        logging.error(f"Analysis directory not found: {analysis_dir}")
        return False
    
    # Process each cell type
    for cell_type in ['Neu', 'NSC']:
        logging.info(f"Processing {cell_type} results")
        
        # Find all chromosome-specific result files for this cell type
        result_files = []
        for chrom in list(range(1, 20)) + ['X', 'Y']:
            chrom_dir = analysis_dir / f"chr{chrom}/chromosomes/chr{chrom}"
            result_file = chrom_dir / f"{cell_type}_chromatin_changes.csv"
            
            if result_file.exists():
                result_files.append(result_file)
                logging.info(f"  Found results for chromosome {chrom}")
            else:
                logging.warning(f"  No results found for chromosome {chrom}")
        
        # Merge results if any were found
        if result_files:
            logging.info(f"  Merging {len(result_files)} chromosome files")
            
            # Read and concatenate all dataframes
            dfs = []
            for file in result_files:
                try:
                    df = pd.read_csv(file)
                    dfs.append(df)
                    logging.info(f"  Read {len(df)} rows from {file}")
                except Exception as e:
                    logging.error(f"  Error reading {file}: {str(e)}")
            
            if dfs:
                # Concatenate all dataframes
                merged_df = pd.concat(dfs, ignore_index=True)
                logging.info(f"  Merged dataframe has {len(merged_df)} rows")
                
                # Save all results
                all_output_file = analysis_dir / f"{cell_type}_chromatin_changes_all.csv"
                merged_df.to_csv(all_output_file, index=False)
                logging.info(f"  Saved all results to {all_output_file}")
                
                # Save significant results
                significant_df = merged_df[merged_df['is_significant'] == True]
                sig_output_file = analysis_dir / f"{cell_type}_chromatin_changes_significant.csv"
                significant_df.to_csv(sig_output_file, index=False)
                logging.info(f"  Saved {len(significant_df)} significant results to {sig_output_file}")
                
                # Save results by change type
                for change_type in ['S2S to S3', 'S3 to S2S']:
                    change_df = merged_df[merged_df['change_type'] == change_type]
                    if len(change_df) > 0:
                        change_file = analysis_dir / f"{cell_type}_chromatin_changes_{change_type.replace(' to ', '_to_')}.csv"
                        change_df.to_csv(change_file, index=False)
                        logging.info(f"  Saved {len(change_df)} {change_type} changes to {change_file}")
            else:
                logging.error(f"  No valid dataframes found for {cell_type}")
        else:
            logging.error(f"  No results found for {cell_type}")
    
    logging.info("Merge operation completed")
    return True

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Merge chromosome-specific results')
    parser.add_argument('--results_dir', default='results', help='Directory containing results')
    args = parser.parse_args()
    
    # Merge results
    success = merge_results(args.results_dir)
    
    # Return appropriate exit code
    return 0 if success else 1

if __name__ == "__main__":
    exit(main())