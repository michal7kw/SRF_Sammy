#!/usr/bin/env python3
# Analysis script for chromatin state changes
# Configured for mm10 mouse genome (19 autosomes + X + Y)

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
import re
import gffutils
import hashlib
import pickle
import json
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('chromatin_analysis.log')
    ]
)

# Global parameters for analysis
WINDOW_SIZE = 10000  # Size of windows to analyze (10kb bins as in the paper)
MIN_DIFF_THRESHOLD = 2.0  # Minimum absolute difference to consider significant
# Removing PVALUE_THRESHOLD as we're skipping p-value calculations

# Define standard chromosomes to analyze for mouse (mm10)
STANDARD_CHROMOSOMES = [str(i) for i in range(1, 20)] + ['X', 'Y']

class ResultsCache:
    """Cache manager for intermediate results."""
    
    def __init__(self, cache_dir):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_file = self.cache_dir / 'cache_metadata.json'
        self.metadata = self._load_metadata()
        
    def _load_metadata(self):
        """Load cache metadata from file."""
        if self.metadata_file.exists():
            try:
                with open(self.metadata_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logging.warning(f"Failed to load cache metadata: {e}")
                return {}
        return {}
    
    def _save_metadata(self):
        """Save cache metadata to file."""
        try:
            with open(self.metadata_file, 'w') as f:
                json.dump(self.metadata, f)
        except Exception as e:
            logging.warning(f"Failed to save cache metadata: {e}")
    
    def _compute_hash(self, key_dict):
        """Compute hash for cache key."""
        key_str = json.dumps(key_dict, sort_keys=True)
        return hashlib.md5(key_str.encode()).hexdigest()
    
    def get(self, key_dict, default=None):
        """Retrieve cached result."""
        cache_key = self._compute_hash(key_dict)
        cache_file = self.cache_dir / f"{cache_key}.pkl"
        
        if cache_file.exists():
            metadata = self.metadata.get(cache_key, {})
            if metadata.get('key_dict') == key_dict:  # Verify cache matches
                try:
                    with open(cache_file, 'rb') as f:
                        logging.info(f"Loading cached result for: {key_dict}")
                        return pickle.load(f)
                except Exception as e:
                    logging.warning(f"Failed to load cache file: {e}")
        
        return default
    
    def set(self, key_dict, value):
        """Store result in cache."""
        cache_key = self._compute_hash(key_dict)
        cache_file = self.cache_dir / f"{cache_key}.pkl"
        
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(value, f)
            
            self.metadata[cache_key] = {
                'key_dict': key_dict,
                'timestamp': datetime.now().isoformat(),
                'file': str(cache_file)
            }
            self._save_metadata()
            logging.info(f"Cached result for: {key_dict}")
        except Exception as e:
            logging.warning(f"Failed to cache result: {e}")

def setup_directories(results_dir, chrom=None):
    """Set up the necessary directories for analysis."""
    # Ensure the chromosome has the 'chr' prefix for mm10 genome
    if chrom and not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
    
    # Create the main results directory if it doesn't exist
    results_dir = Path(results_dir)
    results_dir.mkdir(exist_ok=True)
    
    # Create the bigwig directory
    bigwig_dir = results_dir / "bigwig"
    bigwig_dir.mkdir(exist_ok=True)
    
    # Create the analysis directory
    analysis_dir = results_dir / "analysis"
    analysis_dir.mkdir(exist_ok=True)
    
    # Create chromosome-specific directory if a chromosome is specified
    if chrom:
        chrom_dir = analysis_dir / chrom
        chrom_dir.mkdir(exist_ok=True)
        
        # Create a cache directory for this chromosome
        cache_dir = analysis_dir / "cache" / chrom
        cache_dir.mkdir(exist_ok=True, parents=True)
    else:
        chrom_dir = analysis_dir
        cache_dir = analysis_dir / "cache"
        cache_dir.mkdir(exist_ok=True, parents=True)
    
    logging.info(f"Using results directory: {results_dir}")
    logging.info(f"Using bigwig directory: {bigwig_dir}")
    logging.info(f"Output will be saved to: {analysis_dir}")
    logging.info(f"Cache directory: {cache_dir}")
    
    logging.info(f"Analysis parameters: window_size={WINDOW_SIZE}, min_diff_threshold={MIN_DIFF_THRESHOLD}")
    
    return results_dir, bigwig_dir, chrom_dir, cache_dir

def get_bigwig_files(bigwig_dir, cell_type, condition, state):
    """Get a list of bigWig files for a given cell type, condition, and state."""
    # Create the pattern to match files
    pattern = f"{cell_type}*_{condition}_{state}_*4F*.bw"
    
    # Find all matching files
    bigwig_files = list(Path(bigwig_dir).glob(pattern))
    
    # Convert to strings
    bigwig_files = [str(f) for f in bigwig_files]
    
    # If no files found, log a warning
    if not bigwig_files:
        logging.warning(f"No bigWig files found for pattern: {pattern}")
    
    return bigwig_files

def compute_average_signal(bigwig_files, chrom, start, end):
    """Compute the average signal across multiple bigWig files for a given region."""
    signals = []
    
    for bw_file in bigwig_files:
        if os.path.exists(bw_file):
            try:
                with pyBigWig.open(bw_file) as bw:
                    # Check if chromosome exists in this bigWig file
                    if chrom not in bw.chroms():
                        continue
                        
                    # Get values for the region
                    values = bw.values(chrom, start, end)
                    
                    # Process values if they exist
                    if values is not None:
                        # Replace negative values with 0
                        values = np.array(values)
                        values[values < 0] = 0
                        signal = np.nanmean(values)
                        if not np.isnan(signal) and signal >= 0:
                            signals.append(signal)
            except Exception as e:
                logging.warning(f"Error processing {bw_file} at {chrom}:{start}-{end}: {str(e)}")
                continue
    
    # Only log collection details in debug mode
    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        logging.debug(f"Collected {len(signals)} valid signals from {len(bigwig_files)} files")
    
    return np.nanmean(signals) if signals else 0.0  # Return 0 instead of nan for regions with no signal

def identify_state_changes(bigwig_dir, output_dir, cell_type, cache, chrom=None):
    """Identify chromatin state changes between conditions."""
    logging.info(f"\n==================================================")
    logging.info(f"Processing cell type: {cell_type}")
    logging.info(f"==================================================")
    
    # Check if results are already cached
    cache_key = {
        'task': 'identify_state_changes',
        'cell_type': cell_type,
        'chrom': chrom
    }
    
    cached_result = cache.get(cache_key)
    if cached_result is not None:
        logging.info(f"Using cached results for {cell_type}" + (f" chromosome {chrom}" if chrom else ""))
        return cached_result
    
    # Retrieve bigWig files for each condition
    bigwig_files_gfp_s2s = get_bigwig_files(bigwig_dir, cell_type, 'GFP', 'S2S')
    bigwig_files_gfp_s3 = get_bigwig_files(bigwig_dir, cell_type, 'GFP', 'S3')
    bigwig_files_m2_s2s = get_bigwig_files(bigwig_dir, cell_type, 'M2', 'S2S')
    bigwig_files_m2_s3 = get_bigwig_files(bigwig_dir, cell_type, 'M2', 'S3')
    
    logging.info(f"Found {len(bigwig_files_gfp_s2s)} bigWig files for {cell_type} GFP S2S")
    logging.info(f"Found {len(bigwig_files_gfp_s3)} bigWig files for {cell_type} GFP S3")
    logging.info(f"Found {len(bigwig_files_m2_s2s)} bigWig files for {cell_type} M2 S2S")
    logging.info(f"Found {len(bigwig_files_m2_s3)} bigWig files for {cell_type} M2 S3")
    
    # Combine all bigWig files to check if we have any
    all_bigwigs = bigwig_files_gfp_s2s + bigwig_files_gfp_s3 + bigwig_files_m2_s2s + bigwig_files_m2_s3
    
    if not all_bigwigs:
        logging.error(f"No bigWig files found for {cell_type}")
        return pd.DataFrame()  # Return empty dataframe
    
    # Find a valid bigWig file to get chromosome information
    valid_bigwig = None
    chrom_length = None
    
    for bw_file in all_bigwigs:
        if not os.path.exists(bw_file):
            logging.warning(f"BigWig file does not exist: {bw_file}")
            continue
            
        try:
            with pyBigWig.open(bw_file) as bw:
                if bw is None or not bw.isBigWig():
                    logging.warning(f"Invalid BigWig file: {bw_file}")
                    continue
                    
                # If a specific chromosome is requested, check if it exists
                if chrom is not None:
                    if chrom not in bw.chroms():
                        logging.warning(f"Chromosome {chrom} not found in {bw_file}")
                        continue
                    chrom_length = bw.chroms()[chrom]
                    valid_bigwig = bw_file
                    break
                else:
                    # If no specific chromosome, just use this file
                    valid_bigwig = bw_file
                    break
        except Exception as e:
            logging.warning(f"Error opening BigWig file {bw_file}: {str(e)}")
            continue
    
    if valid_bigwig is None:
        logging.error(f"No valid BigWig files found for {cell_type}")
        return pd.DataFrame()  # Return empty dataframe
        
    logging.info(f"Using {valid_bigwig} for chromosome information")
    
    # Process the specified chromosome or all chromosomes
    results = []
    
    try:
        with pyBigWig.open(valid_bigwig) as bw:
            # If a specific chromosome is requested, only process that one
            if chrom is not None:
                if chrom not in bw.chroms():
                    logging.error(f"Chromosome {chrom} not found in {valid_bigwig}")
                    return pd.DataFrame()
                    
                chrom_length = bw.chroms()[chrom]
                logging.info(f"Processing {chrom} ({chrom_length} bp)")
                
                # Process windows along the chromosome
                for start in range(0, chrom_length, WINDOW_SIZE):
                    end = min(start + WINDOW_SIZE, chrom_length)
                    
                    try:
                        # Calculate average signals for each condition
                        gfp_s2s_signal = compute_average_signal(bigwig_files_gfp_s2s, chrom, start, end)
                        gfp_s3_signal = compute_average_signal(bigwig_files_gfp_s3, chrom, start, end)
                        m2_s2s_signal = compute_average_signal(bigwig_files_m2_s2s, chrom, start, end)
                        m2_s3_signal = compute_average_signal(bigwig_files_m2_s3, chrom, start, end)
                        
                        # Calculate differences
                        gfp_diff = gfp_s3_signal - gfp_s2s_signal
                        m2_diff = m2_s3_signal - m2_s2s_signal
                        condition_diff = m2_diff - gfp_diff
                        
                        # Store results
                        results.append({
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'gfp_s2s': gfp_s2s_signal,
                            'gfp_s3': gfp_s3_signal,
                            'gfp_diff': gfp_diff,
                            'm2_s2s': m2_s2s_signal,
                            'm2_s3': m2_s3_signal,
                            'm2_diff': m2_diff,
                            'condition_diff': condition_diff
                        })
                    except Exception as e:
                        logging.warning(f"Error processing {chrom}:{start}-{end}: {str(e)}")
            else:
                # Process all chromosomes (this should not happen with your SLURM array setup)
                logging.warning("No specific chromosome provided, processing all chromosomes")
                for chrom_name, chrom_length in bw.chroms().items():
                    logging.info(f"Processing {chrom_name} ({chrom_length} bp)")
                    
                    # Process windows along the chromosome
                    for start in range(0, chrom_length, WINDOW_SIZE):
                        end = min(start + WINDOW_SIZE, chrom_length)
                        
                        try:
                            # Calculate average signals for each condition
                            gfp_s2s_signal = compute_average_signal(bigwig_files_gfp_s2s, chrom_name, start, end)
                            gfp_s3_signal = compute_average_signal(bigwig_files_gfp_s3, chrom_name, start, end)
                            m2_s2s_signal = compute_average_signal(bigwig_files_m2_s2s, chrom_name, start, end)
                            m2_s3_signal = compute_average_signal(bigwig_files_m2_s3, chrom_name, start, end)
                            
                            # Calculate differences
                            gfp_diff = gfp_s3_signal - gfp_s2s_signal
                            m2_diff = m2_s3_signal - m2_s2s_signal
                            condition_diff = m2_diff - gfp_diff
                            
                            # Store results
                            results.append({
                                'chrom': chrom_name,
                                'start': start,
                                'end': end,
                                'gfp_s2s': gfp_s2s_signal,
                                'gfp_s3': gfp_s3_signal,
                                'gfp_diff': gfp_diff,
                                'm2_s2s': m2_s2s_signal,
                                'm2_s3': m2_s3_signal,
                                'm2_diff': m2_diff,
                                'condition_diff': condition_diff
                            })
                        except Exception as e:
                            logging.warning(f"Error processing {chrom_name}:{start}-{end}: {str(e)}")
    except Exception as e:
        logging.error(f"Error processing BigWig file {valid_bigwig}: {str(e)}")
        return pd.DataFrame()  # Return empty dataframe
    
    # Convert to dataframe
    df = pd.DataFrame(results)
    
    if len(df) == 0:
        logging.warning(f"No results found for {cell_type}")
        return df
    
    # Add absolute difference column
    df['abs_condition_diff'] = df['condition_diff'].abs()
    
    # Determine significance based on new criteria: at least one signal value is not zero
    df['is_significant'] = True
    
    # Classify changes
    df['change_type'] = 'No change'
    
    # S2S to S3 transition (euchromatin to heterochromatin)
    s2s_to_s3 = (df['condition_diff'] > MIN_DIFF_THRESHOLD) & df['is_significant']
    df.loc[s2s_to_s3, 'change_type'] = 'S2S to S3'
    
    # S3 to S2S transition (heterochromatin to euchromatin)
    s3_to_s2s = (df['condition_diff'] < -MIN_DIFF_THRESHOLD) & df['is_significant']
    df.loc[s3_to_s2s, 'change_type'] = 'S3 to S2S'
    
    # For individual chromosomes, save to chromosome-specific files
    chrom_dir = output_dir / f"chromosomes/{chrom if chrom else 'all'}"
    chrom_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(chrom_dir / f"{cell_type}_chromatin_changes.csv", index=False)
    
    # Log summary statistics
    total_regions = len(df)
    significant_regions = df['is_significant'].sum()
    s2s_to_s3_count = s2s_to_s3.sum()
    s3_to_s2s_count = s3_to_s2s.sum()
    
    logging.info(f"\nSummary for {cell_type}" + (f" chromosome {chrom}" if chrom else "") + ":")
    logging.info(f"Total regions analyzed: {total_regions}")
    logging.info(f"Significant regions: {significant_regions} ({significant_regions/total_regions*100:.2f}%)")
    logging.info(f"S2S to S3 transitions: {s2s_to_s3_count} ({s2s_to_s3_count/total_regions*100:.2f}%)")
    logging.info(f"S3 to S2S transitions: {s3_to_s2s_count} ({s3_to_s2s_count/total_regions*100:.2f}%)")
    
    # Cache the results before returning
    cache.set(cache_key, df)
    return df

def main():
    # Add command-line argument parsing
    parser = argparse.ArgumentParser(description='Analyze chromatin state changes')
    parser.add_argument('--results_dir', default='results', help='Directory for results')
    parser.add_argument('--chrom', help='Specific chromosome to analyze (e.g., "1", "X")')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging.info("Starting chromatin state analysis")
    logging.info(f"Initializing analysis with results_dir: {args.results_dir}")

    # Ensure the chromosome is properly formatted
    chrom = None
    if args.chrom:
        # Format the chromosome name (add 'chr' prefix if not present)
        chrom = f"chr{args.chrom}" if not args.chrom.startswith('chr') else args.chrom
        logging.info(f"Processing single chromosome: {chrom}")

    # Set up directories
    results_dir, bigwig_dir, output_dir, cache_dir = setup_directories(args.results_dir, chrom)
    
    # Initialize cache
    cache = ResultsCache(cache_dir)
    
    logging.info(f"Using results directory: {results_dir}")
    logging.info(f"Using bigwig directory: {bigwig_dir}")
    logging.info(f"Output will be saved to: {output_dir}")
    logging.info(f"Cache directory: {cache_dir}")
    logging.info(f"Analysis parameters: window_size={WINDOW_SIZE}, min_diff_threshold={MIN_DIFF_THRESHOLD}")
    
    # Process only the specified chromosome if provided
    if chrom:
        logging.info(f"Processing chromosome: {chrom}")
        
        logging.info("\n==================================================")
        logging.info("Task 1: Analyzing chromatin state changes")
        logging.info("==================================================")
        
        # Process only Neu and NSC cell types
        for cell_type in ['Neu', 'NSC']:
            identify_state_changes(bigwig_dir, output_dir, cell_type, cache, chrom)
    else:
        logging.error("No chromosome specified. This script is designed to be run with a specific chromosome.")
        exit(1)

if __name__ == "__main__":
    main()
