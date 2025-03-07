#!/usr/bin/env python3
"""
This script analyzes chromatin state changes using bigWig files. It is configured for the mm10 mouse genome (19 autosomes + X + Y).

The script performs the following steps:
1.  Identifies bigWig files for different conditions (e.g., GFP S2S, GFP S3, M2 S2S, M2 S3).
2.  Calculates the average signal intensity within defined genomic windows (default: 10kb).
3.  Computes the difference in signal between conditions to identify chromatin state changes.
4.  Classifies these changes as either S2S to S3 (euchromatin to heterochromatin) or S3 to S2S (heterochromatin to euchromatin).
5.  Stores the results in a pandas DataFrame and saves it to a CSV file.
6.  Caches intermediate results to speed up subsequent runs.

The script uses a minimum difference threshold to filter out noise and focuses on identifying significant changes in chromatin state.
"""

"""
no input gene list
"""

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

# Set up logging: Configure the logging system to output messages to both the console and a file.
logging.basicConfig(
    level=logging.INFO,  # Set the logging level to INFO (other levels: DEBUG, WARNING, ERROR, CRITICAL)
    format='%(asctime)s - %(levelname)s - %(message)s',  # Define the format of log messages
    handlers=[
        logging.StreamHandler(sys.stdout),  # Output log messages to the console (standard output)
        logging.FileHandler('chromatin_analysis.log')  # Output log messages to a file named 'chromatin_analysis.log'
    ]
)

# Global parameters for analysis: Define global variables that control the analysis.
WINDOW_SIZE = 10000  # Size of windows to analyze (10kb bins as in the paper). This determines the resolution of the analysis.
MIN_DIFF_THRESHOLD = 0.0  # Minimum absolute difference in signal to consider a change significant.  This helps filter out noise.
# Removing PVALUE_THRESHOLD as we're skipping p-value calculations. P-value calculations were originally intended for statistical significance testing.

# Define standard chromosomes to analyze for mouse (mm10):  List of standard chromosomes to process.  This ensures consistency and avoids processing unnecessary chromosomes.
STANDARD_CHROMOSOMES = [str(i) for i in range(1, 20)] + ['X', 'Y']

class ResultsCache:
    """Cache manager for intermediate results.  This class handles caching of analysis results to speed up subsequent runs."""
    
    def __init__(self, cache_dir):
        """Initialize the ResultsCache with a directory to store cache files."""
        self.cache_dir = Path(cache_dir)  # Store the cache directory as a Path object
        self.cache_dir.mkdir(parents=True, exist_ok=True)  # Create the cache directory if it doesn't exist
        self.metadata_file = self.cache_dir / 'cache_metadata.json'  # Define the path to the metadata file
        self.metadata = self._load_metadata()  # Load existing metadata from file
        
    def _load_metadata(self):
        """Load cache metadata from file.  This metadata tracks the cached results."""
        if self.metadata_file.exists():  # Check if the metadata file exists
            try:
                with open(self.metadata_file, 'r') as f:  # Open the metadata file in read mode
                    return json.load(f)  # Load the JSON data from the file
            except Exception as e:
                logging.warning(f"Failed to load cache metadata: {e}")  # Log a warning if loading fails
                return {}  # Return an empty dictionary if loading fails
        return {}  # Return an empty dictionary if the metadata file doesn't exist
    
    def _save_metadata(self):
        """Save cache metadata to file.  This ensures the cache is persistent."""
        try:
            with open(self.metadata_file, 'w') as f:  # Open the metadata file in write mode
                json.dump(self.metadata, f)  # Save the metadata to the file as JSON
        except Exception as e:
            logging.warning(f"Failed to save cache metadata: {e}")  # Log a warning if saving fails
    
    def _compute_hash(self, key_dict):
        """Compute hash for cache key.  This creates a unique identifier for each set of parameters."""
        # Add a script-specific prefix to ensure uniqueness across scripts
        key_dict_with_prefix = {"script": "task1_chromatin_changes_GFP_M2", **key_dict}
        key_str = json.dumps(key_dict_with_prefix, sort_keys=True)  # Convert the key dictionary to a JSON string (sorted for consistency)
        return hashlib.md5(key_str.encode()).hexdigest()  # Compute the MD5 hash of the JSON string
    
    def get(self, key_dict, default=None):
        """Retrieve cached result.  If the result is not cached, return the default value."""
        cache_key = self._compute_hash(key_dict)  # Compute the cache key based on the input dictionary
        cache_file = self.cache_dir / f"{cache_key}.pkl"  # Define the path to the cache file
        
        if cache_file.exists():  # Check if the cache file exists
            metadata = self.metadata.get(cache_key, {})  # Retrieve the metadata for this cache key
            if metadata.get('key_dict') == key_dict:  # Verify that the cached data matches the input key
                try:
                    with open(cache_file, 'rb') as f:  # Open the cache file in read binary mode
                        logging.info(f"Loading cached result for: {key_dict}")  # Log that we're loading from the cache
                        return pickle.load(f)  # Load the data from the cache file using pickle
                except Exception as e:
                    logging.warning(f"Failed to load cache file: {e}")  # Log a warning if loading fails
        
        return default  # Return the default value if the result is not cached or loading fails
    
    def set(self, key_dict, value):
        """Store result in cache.  This saves the result for future use."""
        cache_key = self._compute_hash(key_dict)  # Compute the cache key based on the input dictionary
        cache_file = self.cache_dir / f"{cache_key}.pkl"  # Define the path to the cache file
        
        try:
            with open(cache_file, 'wb') as f:  # Open the cache file in write binary mode
                pickle.dump(value, f)  # Save the data to the cache file using pickle
            
            self.metadata[cache_key] = {  # Update the metadata with information about the cached file
                'key_dict': key_dict,  # Store the key dictionary
                'timestamp': datetime.now().isoformat(),  # Store the current timestamp
                'file': str(cache_file)  # Store the path to the cache file
            }
            self._save_metadata()  # Save the updated metadata to file
            logging.info(f"Cached result for: {key_dict}")  # Log that we've cached the result
        except Exception as e:
            logging.warning(f"Failed to cache result: {e}")  # Log a warning if caching fails

def setup_directories(results_dir, chrom=None):
    """Set up the necessary directories for analysis.  This function creates the directory structure for storing results."""
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
        
        # Create a script-specific cache directory for this chromosome
        cache_dir = analysis_dir / "cache" / "task1_cache" / chrom
        cache_dir.mkdir(exist_ok=True, parents=True)
    else:
        chrom_dir = analysis_dir
        cache_dir = analysis_dir / "cache" / "task1_cache"
        cache_dir.mkdir(exist_ok=True, parents=True)
    
    logging.info(f"Using results directory: {results_dir}")
    logging.info(f"Using bigwig directory: {bigwig_dir}")
    logging.info(f"Output will be saved to: {analysis_dir}")
    logging.info(f"Cache directory: {cache_dir}")
    
    logging.info(f"Analysis parameters: window_size={WINDOW_SIZE}, min_diff_threshold={MIN_DIFF_THRESHOLD}")
    
    return results_dir, bigwig_dir, chrom_dir, cache_dir

def get_bigwig_files(bigwig_dir, cell_type, condition, state):
    """Get a list of bigWig files for a given cell type, condition, and state.  This function finds the relevant bigWig files for analysis."""
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
    """Compute the average signal across multiple bigWig files for a given region.  This function calculates the average signal intensity within a specified genomic region."""
    signals = []  # Initialize an empty list to store signal values from each bigWig file
    
    for bw_file in bigwig_files:  # Iterate through each bigWig file
        if os.path.exists(bw_file):  # Check if the bigWig file exists
            try:
                with pyBigWig.open(bw_file) as bw:  # Open the bigWig file using pyBigWig
                    # Check if chromosome exists in this bigWig file
                    if chrom not in bw.chroms():
                        continue  # Skip this file if the chromosome is not present
                        
                    # Get values for the region
                    values = bw.values(chrom, start, end)  # Extract signal values for the specified chromosome, start, and end positions
                    
                    # Process values if they exist
                    if values is not None:
                        # Replace negative values with 0
                        values = np.array(values)  # Convert the values to a NumPy array
                        values[values < 0] = 0  # Set any negative values to 0 (important for signal processing)
                        signal = np.nanmean(values)  # Calculate the mean of the signal values, ignoring NaN values
                        if not np.isnan(signal) and signal >= 0:  # Check if the signal is a valid number and non-negative
                            signals.append(signal)  # Add the signal to the list of signals
            except Exception as e:
                logging.warning(f"Error processing {bw_file} at {chrom}:{start}-{end}: {str(e)}")  # Log a warning if an error occurs during processing
                continue  # Continue to the next bigWig file
    
    # Only log collection details in debug mode
    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        logging.debug(f"Collected {len(signals)} valid signals from {len(bigwig_files)} files")  # Log the number of valid signals collected
    
    return np.nanmean(signals) if signals else 0.0  # Return the average signal across all files, or 0.0 if no signals were collected
def identify_state_changes(bigwig_dir, output_dir, cell_type, cache, chrom=None):
    """
    Identify chromatin state changes between conditions (GFP S2S/S3 vs M2 S2S/S3).
    This is the core analysis function that orchestrates the entire process:
    1. Retrieves relevant bigWig files for each condition.
    2. Processes the bigWig files to calculate signal differences.
    3. Identifies significant chromatin state changes based on defined thresholds.
    4. Stores the results in a pandas DataFrame and saves it to a file.
    5. Caches the results for faster retrieval in subsequent runs.

    Args:
        bigwig_dir (str): Directory containing the bigWig files.
        output_dir (str): Directory to save the analysis results.
        cell_type (str): The cell type to analyze (e.g., 'Neu', 'NSC').
        cache (ResultsCache): A cache object to store and retrieve results.
        chrom (str, optional): Specific chromosome to analyze (e.g., 'chr1'). If None, all chromosomes are processed. Defaults to None.

    Returns:
        pandas.DataFrame: A DataFrame containing the analysis results, including signal differences and identified state changes.
                          Returns an empty DataFrame if no data is found or an error occurs.
    """
    logging.info(f"\n==================================================")
    logging.info(f"Processing cell type: {cell_type}")
    logging.info(f"==================================================")
    
    # Define a unique cache key based on the function's input parameters
    cache_key = {
        'task': 'identify_state_changes',  # Task identifier
        'cell_type': cell_type,            # Cell type being analyzed
        'chrom': chrom                     # Chromosome being analyzed (or None if all)
    }
    
    # Check if the results for this specific analysis are already stored in the cache
    cached_result = cache.get(cache_key)
    if cached_result is not None:
        logging.info(f"Using cached results for {cell_type}" + (f" chromosome {chrom}" if chrom else ""))
        return cached_result  # Return the cached DataFrame if available
    
    # If the results are not cached, retrieve the necessary bigWig files for each condition
    bigwig_files_gfp_s2s = get_bigwig_files(bigwig_dir, cell_type, 'GFP', 'S2S')  # GFP S2S condition
    bigwig_files_gfp_s3 = get_bigwig_files(bigwig_dir, cell_type, 'GFP', 'S3')    # GFP S3 condition
    bigwig_files_m2_s2s = get_bigwig_files(bigwig_dir, cell_type, 'M2', 'S2S')    # M2 S2S condition
    bigwig_files_m2_s3 = get_bigwig_files(bigwig_dir, cell_type, 'M2', 'S3')      # M2 S3 condition
    
    # Log the number of bigWig files found for each condition
    logging.info(f"Found {len(bigwig_files_gfp_s2s)} bigWig files for {cell_type} GFP S2S")
    logging.info(f"Found {len(bigwig_files_gfp_s3)} bigWig files for {cell_type} GFP S3")
    logging.info(f"Found {len(bigwig_files_m2_s2s)} bigWig files for {cell_type} M2 S2S")
    logging.info(f"Found {len(bigwig_files_m2_s3)} bigWig files for {cell_type} M2 S3")
    
    # Combine all bigWig files into a single list to check for the existence of any files
    all_bigwigs = bigwig_files_gfp_s2s + bigwig_files_gfp_s3 + bigwig_files_m2_s2s + bigwig_files_m2_s3
    
    # If no bigWig files were found for any condition, log an error and return an empty DataFrame
    if not all_bigwigs:
        logging.error(f"No bigWig files found for {cell_type}")
        return pd.DataFrame()  # Return empty dataframe
    
    # Find a valid bigWig file to extract chromosome information (length)
    valid_bigwig = None  # Initialize the variable to store the path to a valid bigWig file
    chrom_length = None  # Initialize the variable to store the chromosome length
    
    # Iterate through all bigWig files to find a valid one
    for bw_file in all_bigwigs:
        if not os.path.exists(bw_file):  # Check if the file exists
            logging.warning(f"BigWig file does not exist: {bw_file}")
            continue  # Skip to the next file if it doesn't exist
            
        try:
            with pyBigWig.open(bw_file) as bw:  # Open the bigWig file using pyBigWig
                if bw is None or not bw.isBigWig():  # Check if the file is a valid bigWig file
                    logging.warning(f"Invalid BigWig file: {bw_file}")
                    continue  # Skip to the next file if it's invalid
                    
                # If a specific chromosome is requested, check if it exists in the current bigWig file
                if chrom is not None:
                    if chrom not in bw.chroms():  # Check if the chromosome exists in the bigWig file
                        logging.warning(f"Chromosome {chrom} not found in {bw_file}")
                        continue  # Skip to the next file if the chromosome is not found
                    chrom_length = bw.chroms()[chrom]  # Get the length of the chromosome
                    valid_bigwig = bw_file  # Set the current file as the valid bigWig file
                    break  # Exit the loop since we found a valid file and chromosome
                else:
                    # If no specific chromosome is requested, use the first valid bigWig file found
                    valid_bigwig = bw_file  # Set the current file as the valid bigWig file
                    break  # Exit the loop since we found a valid file
        except Exception as e:
            logging.warning(f"Error opening BigWig file {bw_file}: {str(e)}")
            continue  # Skip to the next file if an error occurred while opening it
    
    # If no valid bigWig file was found, log an error and return an empty DataFrame
    if valid_bigwig is None:
        logging.error(f"No valid BigWig files found for {cell_type}")
        return pd.DataFrame()  # Return empty dataframe
        
    logging.info(f"Using {valid_bigwig} for chromosome information")
    
    # Initialize an empty list to store the analysis results for each genomic region
    results = []
    
    try:
        with pyBigWig.open(valid_bigwig) as bw:  # Open the valid bigWig file
            # If a specific chromosome is requested, process only that chromosome
            if chrom is not None:
                if chrom not in bw.chroms():  # Verify that the chromosome exists in the bigWig file
                    logging.error(f"Chromosome {chrom} not found in {valid_bigwig}")
                    return pd.DataFrame()  # Return an empty DataFrame if the chromosome is not found
                    
                chrom_length = bw.chroms()[chrom]  # Get the length of the chromosome from the bigWig file
                logging.info(f"Processing {chrom} ({chrom_length} bp)")
                
                # Iterate over the chromosome in windows of size WINDOW_SIZE
                for start in range(0, chrom_length, WINDOW_SIZE):
                    end = min(start + WINDOW_SIZE, chrom_length)  # Calculate the end position of the window
                    
                    try:
                        # Calculate the average signal for each condition within the current window
                        gfp_s2s_signal = compute_average_signal(bigwig_files_gfp_s2s, chrom, start, end)
                        gfp_s3_signal = compute_average_signal(bigwig_files_gfp_s3, chrom, start, end)
                        m2_s2s_signal = compute_average_signal(bigwig_files_m2_s2s, chrom, start, end)
                        m2_s3_signal = compute_average_signal(bigwig_files_m2_s3, chrom, start, end)
                        
                        # Calculate the difference in signal between S3 and S2S for each condition
                        gfp_diff = gfp_s3_signal - gfp_s2s_signal
                        m2_diff = m2_s3_signal - m2_s2s_signal
                        
                        # Calculate the difference between the M2 and GFP signal differences
                        condition_diff = m2_diff - gfp_diff
                        
                        # Store the results for the current window in a dictionary
                        results.append({
                            'chrom': chrom,              # Chromosome name
                            'start': start,              # Start position of the window
                            'end': end,                  # End position of the window
                            'gfp_s2s': gfp_s2s_signal,    # Average signal for GFP S2S
                            'gfp_s3': gfp_s3_signal,      # Average signal for GFP S3
                            'gfp_diff': gfp_diff,        # Difference between GFP S3 and S2S
                            'm2_s2s': m2_s2s_signal,      # Average signal for M2 S2S
                            'm2_s3': m2_s3_signal,        # Average signal for M2 S3
                            'm2_diff': m2_diff,          # Difference between M2 S3 and S2S
                            'condition_diff': condition_diff  # Difference between M2 and GFP signal differences
                        })
                    except Exception as e:
                        logging.warning(f"Error processing {chrom}:{start}-{end}: {str(e)}")  # Log any errors that occur during processing
            else:
                # If no specific chromosome is provided, process all chromosomes in the bigWig file
                logging.warning("No specific chromosome provided, processing all chromosomes")
                for chrom_name, chrom_length in bw.chroms().items():  # Iterate over each chromosome in the bigWig file
                    logging.info(f"Processing {chrom_name} ({chrom_length} bp)")
                    
                    # Iterate over the chromosome in windows of size WINDOW_SIZE
                    for start in range(0, chrom_length, WINDOW_SIZE):
                        end = min(start + WINDOW_SIZE, chrom_length)  # Calculate the end position of the window
                        
                        try:
                            # Calculate the average signal for each condition within the current window
                            gfp_s2s_signal = compute_average_signal(bigwig_files_gfp_s2s, chrom_name, start, end)
                            gfp_s3_signal = compute_average_signal(bigwig_files_gfp_s3, chrom_name, start, end)
                            m2_s2s_signal = compute_average_signal(bigwig_files_m2_s2s, chrom_name, start, end)
                            m2_s3_signal = compute_average_signal(bigwig_files_m2_s3, chrom_name, start, end)
                            
                            # Calculate the difference in signal between S3 and S2S for each condition
                            gfp_diff = gfp_s3_signal - gfp_s2s_signal
                            m2_diff = m2_s3_signal - m2_s2s_signal
                            
                            # Calculate the difference between the M2 and GFP signal differences
                            condition_diff = m2_diff - gfp_diff
                            
                            # Store the results for the current window in a dictionary
                            results.append({
                                'chrom': chrom_name,          # Chromosome name
                                'start': start,              # Start position of the window
                                'end': end,                  # End position of the window
                                'gfp_s2s': gfp_s2s_signal,    # Average signal for GFP S2S
                                'gfp_s3': gfp_s3_signal,      # Average signal for GFP S3
                                'gfp_diff': gfp_diff,        # Difference between GFP S3 and S2S
                                'm2_s2s': m2_s2s_signal,      # Average signal for M2 S2S
                                'm2_s3': m2_s3_signal,        # Average signal for M2 S3
                                'm2_diff': m2_diff,          # Difference between M2 S3 and S2S
                                'condition_diff': condition_diff  # Difference between M2 and GFP signal differences
                            })
                        except Exception as e:
                            logging.warning(f"Error processing {chrom_name}:{start}-{end}: {str(e)}")  # Log any errors that occur during processing
    except Exception as e:
        logging.error(f"Error processing BigWig file {valid_bigwig}: {str(e)}")  # Log any errors that occur while processing the bigWig file
        return pd.DataFrame()  # Return empty dataframe
    
    # Convert the list of results into a pandas DataFrame
    df = pd.DataFrame(results)
    
    # If no results were found, log a warning and return the empty DataFrame
    if len(df) == 0:
        logging.warning(f"No results found for {cell_type}")
        return df
    
    # Add a column for the absolute value of the condition difference
    df['abs_condition_diff'] = df['condition_diff'].abs()
    
    # Determine significance based on new criteria: at least one signal value is not zero
    df['is_significant'] = True
    
    # Classify changes
    df['change_type'] = 'No change'
    
    # Identify regions where the chromatin state transitions from S2S to S3 (euchromatin to heterochromatin)
    s2s_to_s3 = (df['condition_diff'] > MIN_DIFF_THRESHOLD) & df['is_significant']
    df.loc[s2s_to_s3, 'change_type'] = 'S2S to S3'
    
    # Identify regions where the chromatin state transitions from S3 to S2S (heterochromatin to euchromatin)
    s3_to_s2s = (df['condition_diff'] < -MIN_DIFF_THRESHOLD) & df['is_significant']
    df.loc[s3_to_s2s, 'change_type'] = 'S3 to S2S'
    
    # Create a directory to store the results for the current chromosome
    chrom_dir = output_dir / f"chromosomes/{chrom if chrom else 'all'}"
    chrom_dir.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist
    df.to_csv(chrom_dir / f"{cell_type}_chromatin_changes.csv", index=False)  # Save the DataFrame to a CSV file
    
    # Calculate and log summary statistics
    total_regions = len(df)  # Total number of regions analyzed
    significant_regions = df['is_significant'].sum()  # Number of significant regions
    s2s_to_s3_count = s2s_to_s3.sum()  # Number of S2S to S3 transitions
    s3_to_s2s_count = s3_to_s2s.sum()  # Number of S3 to S2S transitions
    
    logging.info(f"\nSummary for {cell_type}" + (f" chromosome {chrom}" if chrom else "") + ":")
    logging.info(f"Total regions analyzed: {total_regions}")
    logging.info(f"Significant regions: {significant_regions} ({significant_regions/total_regions*100:.2f}%)")
    logging.info(f"S2S to S3 transitions: {s2s_to_s3_count} ({s2s_to_s3_count/total_regions*100:.2f}%)")
    logging.info(f"S3 to S2S transitions: {s3_to_s2s_count} ({s3_to_s2s_count/total_regions*100:.2f}%)")
    
    # Store the results in the cache before returning
    cache.set(cache_key, df)
    return df  # Return the DataFrame containing the analysis results

def main():
    """Main function to execute the chromatin state analysis."""
    # Add command-line argument parsing
    parser = argparse.ArgumentParser(description='Analyze chromatin state changes')
    parser.add_argument('--results_dir', default='../results', help='Directory for results')
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
