#!/usr/bin/env python3
# 3_analyze_exo_enriched_genes.py
# Script to analyze chromatin state changes in genes where exogenous Mecp2 is enriched (FC>2, exo vs endo)
# compared to the GFP control in both NSCs and Neurons

import os  # For interacting with the operating system
import sys  # For system-specific parameters and functions
import logging  # For logging messages
import argparse  # For parsing command-line arguments
import hashlib  # For creating hash values
import json  # For working with JSON data
import numpy as np  # For numerical opediffns
import pandas as pd  # For data manipulation and analysis
import matplotlib.pyplot as plt  # For plotting
import seaborn as sns  # For statistical data visualization
from pathlib import Path  # For working with file paths
import pyBigWig  # For accessing BigWig files
import pickle  # For serializing and de-serializing Python objects
from collections import defaultdict  # For creating dictionaries with default values

# Constants
WINDOW_SIZE = 10000  # Window size for gene promoter analysis (10kb).  Defines the region around the TSS to analyze.
MIN_DIFF_THRESHOLD = 2.0  # Minimum difference threshold for state changes.  Used to classify chromatin states based on S3/S2S diffs.

class ResultsCache:
    """Cache for storing computation results to avoid redundant calculations.
    This class uses a directory to store pickled results, indexed by a hash of the input parameters.
    """
    
    def __init__(self, cache_dir):
        """Initialize the cache with the specified directory.
        
        Args:
            cache_dir (str): The directory where cached results will be stored.
        """
        self.cache_dir = Path(cache_dir)  # Convert cache directory to a Path object
        self.cache_dir.mkdir(parents=True, exist_ok=True)  # Create the cache directory if it doesn't exist
        self.metadata_file = self.cache_dir / "metadata.json"  # Define the metadata file path
        self.metadata = self._load_metadata()  # Load existing metadata or initialize an empty dictionary
        
    def _load_metadata(self):
        """Load metadata from file or initialize if it doesn't exist.
        
        Returns:
            dict: A dictionary containing the cache metadata.
        """
        if self.metadata_file.exists():  # Check if the metadata file exists
            try:
                with open(self.metadata_file, 'r') as f:  # Open the metadata file for reading
                    return json.load(f)  # Load the metadata from the JSON file
            except (json.JSONDecodeError, IOError) as e:  # Handle potential errors during loading
                logging.warning(f"Failed to load cache metadata: {e}")  # Log a warning message
                return {}  # Return an empty dictionary if loading fails
        return {}  # Return an empty dictionary if the metadata file doesn't exist
    
    def _save_metadata(self):
        """Save metadata to file.
        This function saves the current metadata to the metadata file in JSON format.
        """
        try:
            with open(self.metadata_file, 'w') as f:  # Open the metadata file for writing
                json.dump(self.metadata, f)  # Save the metadata to the JSON file
        except IOError as e:  # Handle potential errors during saving
            logging.warning(f"Failed to save cache metadata: {e}")  # Log a warning message
    
    def _compute_hash(self, key_dict):
        """Compute a hash for the cache key.
        
        Args:
            key_dict (dict): The dictionary to hash.
        
        Returns:
            str: The hexadecimal hash of the key dictionary.
        """
        # Add a script-specific prefix to ensure uniqueness across scripts
        key_dict_with_prefix = {"script": "task3_exo_enriched", **key_dict}
        key_str = json.dumps(key_dict_with_prefix, sort_keys=True)  # Convert the key dictionary to a sorted JSON string
        return hashlib.md5(key_str.encode()).hexdigest()  # Compute the MD5 hash of the JSON string
    
    def get(self, key_dict, default=None):
        """Get a value from the cache using the key dictionary.
        
        Args:
            key_dict (dict): A dictionary containing the key parameters.
            default (any, optional): The default value to return if the key is not found in the cache. Defaults to None.
        
        Returns:
            any: The cached value if found, otherwise the default value.
        """
        key_hash = self._compute_hash(key_dict)  # Compute the hash of the key dictionary
        if key_hash in self.metadata:  # Check if the hash exists in the metadata
            cache_file = self.cache_dir / f"{key_hash}.pkl"  # Define the cache file path
            if cache_file.exists():  # Check if the cache file exists
                try:
                    with open(cache_file, 'rb') as f:  # Open the cache file for reading in binary mode
                        logging.debug(f"Cache hit for {key_hash}")  # Log a debug message
                        return pickle.load(f)  # Load and return the cached value
                except (pickle.PickleError, IOError) as e:  # Handle potential errors during loading
                    logging.warning(f"Failed to load cached data: {e}")  # Log a warning message
        return default  # Return the default value if the key is not found
    
    def set(self, key_dict, value):
        """Set a value in the cache using the key dictionary.
        
        Args:
            key_dict (dict): A dictionary containing the key parameters.
            value (any): The value to cache.
        
        Returns:
            bool: True if the value was successfully cached, False otherwise.
        """
        key_hash = self._compute_hash(key_dict)  # Compute the hash of the key dictionary
        cache_file = self.cache_dir / f"{key_hash}.pkl"  # Define the cache file path
        try:
            with open(cache_file, 'wb') as f:  # Open the cache file for writing in binary mode
                pickle.dump(value, f)  # Serialize and save the value to the cache file
            self.metadata[key_hash] = {  # Update the metadata
                'key': key_dict,
                'file': str(cache_file)
            }
            self._save_metadata()  # Save the updated metadata
            logging.debug(f"Cached data for {key_hash}")  # Log a debug message
            return True  # Return True if caching was successful
        except (pickle.PickleError, IOError) as e:  # Handle potential errors during saving
            logging.warning(f"Failed to cache data: {e}")  # Log a warning message
            return False  # Return False if caching failed

def setup_directories(results_dir, chrom=None):
    """Set up the necessary directories for analysis.
    
    Args:
        results_dir (str): The main directory for storing results.
        chrom (str, optional): The specific chromosome being analyzed. Defaults to None.
    
    Returns:
        tuple: A tuple containing the results directory, bigwig directory, output directory, and cache directory.
    """
    # Create main results directory
    results_dir = Path(results_dir)  # Convert results directory to a Path object
    results_dir.mkdir(parents=True, exist_ok=True)  # Create the results directory if it doesn't exist
    
    # Create subdirectories
    bigwig_dir = results_dir / "bigwig"  # Define the bigwig directory path
    
    # Create output directory with chromosome-specific subfolder if needed
    if chrom:  # If a specific chromosome is being analyzed
        output_dir = results_dir / "exo_enriched_analysis" / chrom  # Create a chromosome-specific output directory
    else:  # If all chromosomes are being analyzed
        output_dir = results_dir / "exo_enriched_analysis"  # Create a general output directory
    
    output_dir.mkdir(parents=True, exist_ok=True)  # Create the output directory if it doesn't exist
    
    # Create script-specific cache directory
    cache_dir = results_dir / "cache" / "task3_cache"  # Define the script-specific cache directory path
    cache_dir.mkdir(parents=True, exist_ok=True)  # Create the cache directory if it doesn't exist
    
    return results_dir, bigwig_dir, output_dir, cache_dir  # Return the directory paths

def get_bigwig_files(bigwig_dir, cell_type, condition, state):
    """Get a list of BigWig files for the specified parameters.
    
    Args:
        bigwig_dir (str): The directory containing the BigWig files.
        cell_type (str): The cell type (e.g., "Neu", "NSC").
        condition (str): The experimental condition (e.g., "M2", "GFP").
        state (str): The chromatin state (e.g., "S2S", "S3").
    
    Returns:
        list: A list of Path objects representing the BigWig files that match the specified parameters.
    """
    bigwig_dir = Path(bigwig_dir)  # Convert bigwig directory to a Path object
    pattern = f"*{cell_type}*{condition}*{state}*.bw"  # Define the file pattern to search for
    files = list(bigwig_dir.glob(pattern))  # Find all files that match the pattern
    
    if not files:  # If no files were found
        logging.warning(f"No BigWig files found for {cell_type} {condition} {state} in {bigwig_dir}")  # Log a warning message
    else:  # If files were found
        logging.info(f"Found {len(files)} BigWig files for {cell_type} {condition} {state}")  # Log an info message
        for f in files:  # Iterate over the found files
            logging.debug(f"  - {f.name}")  # Log the name of each file
    
    return files  # Return the list of files

def compute_average_signal(bigwig_files, chrom, start, end):
    """Compute the average signal across multiple BigWig files for a genomic region.
    
    Args:
        bigwig_files (list): A list of Path objects representing the BigWig files to analyze.
        chrom (str): The chromosome.
        start (int): The start coordinate of the region.
        end (int): The end coordinate of the region.
    
    Returns:
        float: The average signal across all BigWig files for the specified region, or None if no signal was found.
    """
    if not bigwig_files:  # If no BigWig files were provided
        return None  # Return None
    
    signals = []  # Initialize an empty list to store the signals from each BigWig file
    for bw_file in bigwig_files:  # Iterate over the BigWig files
        try:
            bw = pyBigWig.open(str(bw_file))  # Open the BigWig file
            if chrom in bw.chroms():  # Check if the chromosome exists in the BigWig file
                # Ensure the region is within chromosome bounds
                chrom_size = bw.chroms()[chrom]  # Get the size of the chromosome
                valid_start = max(0, start)  # Ensure the start coordinate is not negative
                valid_end = min(end, chrom_size)  # Ensure the end coordinate is not greater than the chromosome size
                
                if valid_end > valid_start:  # If the region is valid
                    # Get values and compute mean, ignoring NaN values
                    values = bw.values(chrom, valid_start, valid_end)  # Get the signal values for the region
                    values = np.array(values)  # Convert the values to a NumPy array
                    values = values[~np.isnan(values)]  # Remove NaN values
                    if len(values) > 0:  # If there are any valid values
                        signals.append(np.mean(values))  # Compute the mean of the values and append it to the list of signals
            bw.close()  # Close the BigWig file
        except Exception as e:  # Handle potential errors during processing
            logging.error(f"Error processing {bw_file}: {e}")  # Log an error message
    
    if signals:  # If any signals were found
        return np.mean(signals)  # Return the average of the signals
    return None  # Return None if no signals were found

def load_gene_list(gene_list_file):
    """Load a list of genes from a file.
    
    Args:
        gene_list_file (str): The path to the file containing the list of genes.
    
    Returns:
        list: A list of gene names.
    """
    genes = []  # Initialize an empty list to store the genes
    try:
        with open(gene_list_file, 'r') as f:  # Open the gene list file for reading
            for line in f:  # Iterate over the lines in the file
                gene = line.strip()  # Remove leading/trailing whitespace from the line
                if gene:  # Skip empty lines
                    genes.append(gene)  # Append the gene to the list
        logging.info(f"Loaded {len(genes)} genes from {gene_list_file}")  # Log an info message
    except Exception as e:  # Handle potential errors during loading
        logging.error(f"Error loading gene list from {gene_list_file}: {e}")  # Log an error message
    
    return genes  # Return the list of genes

def load_gene_coordinates(genes, gtf_file):
    """Load gene coordinates from a GTF file for the specified genes.
    
    Args:
        genes (list): A list of gene names to load coordinates for.
        gtf_file (str): The path to the GTF file.
    
    Returns:
        dict: A dictionary mapping gene names to their coordinates (chromosome, start, end, strand, gene_id, gene_name).
    """
    gene_coords = {}  # Initialize an empty dictionary to store the gene coordinates
    
    try:
        # Read GTF file and extract gene coordinates
        with open(gtf_file, 'r') as f:  # Open the GTF file for reading
            for line in f:  # Iterate over the lines in the file
                if line.startswith('#'):  # Skip header lines
                    continue
                
                fields = line.strip().split('\t')  # Split the line into fields
                if len(fields) < 9 or fields[2] != 'gene':  # Skip lines that are not gene entries
                    continue
                
                # Extract gene ID from attributes
                attributes = fields[8]  # Get the attributes field
                gene_id = None  # Initialize gene_id
                gene_name = None  # Initialize gene_name
                
                # Parse attributes to find gene_id and gene_name
                for attr in attributes.split(';'):  # Iterate over the attributes
                    attr = attr.strip()  # Remove leading/trailing whitespace from the attribute
                    if attr.startswith('gene_id'):  # If the attribute is gene_id
                        gene_id = attr.split('"')[1] if '"' in attr else attr.split('=')[1]  # Extract the gene ID
                    elif attr.startswith('gene_name'):  # If the attribute is gene_name
                        gene_name = attr.split('"')[1] if '"' in attr else attr.split('=')[1]  # Extract the gene name
                
                # Check if this gene is in our list
                if gene_name in genes or gene_id in genes:  # If the gene is in the list of genes to load
                    chrom = fields[0]  # Get the chromosome
                    start = int(fields[3])  # Get the start coordinate
                    end = int(fields[4])  # Get the end coordinate
                    strand = fields[6]  # Get the strand
                    
                    # Use gene_name if available, otherwise use gene_id
                    gene_key = gene_name if gene_name in genes else gene_id  # Determine the gene key
                    
                    # Store gene coordinates
                    gene_coords[gene_key] = {  # Store the gene coordinates in the dictionary
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'gene_id': gene_id,
                        'gene_name': gene_name
                    }
        
        logging.info(f"Mapped coordinates for {len(gene_coords)} out of {len(genes)} genes")  # Log an info message
    except Exception as e:  # Handle potential errors during loading
        logging.error(f"Error loading gene coordinates from {gtf_file}: {e}")  # Log an error message
    
    return gene_coords  # Return the dictionary of gene coordinates

def analyze_gene_chromatin_state(gene_coords, s2s_files, s3_files, cache):
    """Analyze the chromatin state of genes by comparing S2S and S3 signals.
    
    Args:
        gene_coords (dict): A dictionary mapping gene names to their coordinates.
        s2s_files (list): A list of Path objects representing the BigWig files for the S2S condition.
        s3_files (list): A list of Path objects representing the BigWig files for the S3 condition.
        cache (ResultsCache): The cache object to use for caching results.
    
    Returns:
        pd.DataFrame: A DataFrame containing the analysis results for each gene.
    """
    results = []  # Initialize an empty list to store the results
    
    for gene, coords in gene_coords.items():  # Iterate over the genes and their coordinates
        chrom = coords['chrom']  # Get the chromosome
        
        # Define promoter region (5kb upstream of TSS)
        if coords['strand'] == '+':  # If the gene is on the positive strand
            promoter_start = max(0, coords['start'] - WINDOW_SIZE)  # Define the promoter start coordinate
            promoter_end = coords['start']  # Define the promoter end coordinate
        else:  # If the gene is on the negative strand
            promoter_start = coords['end']  # Define the promoter start coordinate
            promoter_end = coords['end'] + WINDOW_SIZE  # Define the promoter end coordinate
        
        # Create cache key for this computation
        cache_key = {  # Define the cache key
            'function': 'analyze_gene_chromatin_state',
            'gene': gene,
            'chrom': chrom,
            'promoter_start': promoter_start,
            'promoter_end': promoter_end,
            's2s_files': [str(f) for f in s2s_files],
            's3_files': [str(f) for f in s3_files]
        }
        
        # Try to get from cache
        cached_result = cache.get(cache_key)  # Try to get the result from the cache
        if cached_result is not None:  # If the result was found in the cache
            results.append(cached_result)  # Append the cached result to the list of results
            continue  # Skip to the next gene
        
        # Compute average signals
        s2s_signal = compute_average_signal(s2s_files, chrom, promoter_start, promoter_end)  # Compute the average signal for the S2S condition
        s3_signal = compute_average_signal(s3_files, chrom, promoter_start, promoter_end)  # Compute the average signal for the S3 condition
        
        # Determine chromatin state
        if s2s_signal is not None and s3_signal is not None:  # If both signals were computed
            diff = s3_signal - s2s_signal  # Compute the diff of S3 to S2S signals
            
            # Classify based on diff
            if diff > 1 + MIN_DIFF_THRESHOLD:  # If the diff is greater than 1 + the minimum difference threshold
                state = "Heterochromatin"  # S3 dominant. Classify as heterochromatin
            elif diff < 1 - MIN_DIFF_THRESHOLD:  # If the diff is less than 1 - the minimum difference threshold
                state = "Euchromatin"  # S2S dominant. Classify as euchromatin
            else:  # Otherwise
                state = "Mixed"  # No clear dominance. Classify as mixed
            
            result = {  # Create a dictionary to store the results for this gene
                'gene': gene,
                'chrom': chrom,
                'start': coords['start'],
                'end': coords['end'],
                'strand': coords['strand'],
                'promoter_start': promoter_start,
                'promoter_end': promoter_end,
                's2s_signal': s2s_signal,
                's3_signal': s3_signal,
                'diff': diff,
                'state': state
            }
            
            # Cache and append result
            cache.set(cache_key, result)  # Cache the result
            results.append(result)  # Append the result to the list of results
        else:  # If either signal could not be computed
            logging.warning(f"Could not compute signals for gene {gene} at {chrom}:{promoter_start}-{promoter_end}")  # Log a warning message
    
    # Convert results to DataFrame
    df = pd.DataFrame(results)  # Convert the list of results to a DataFrame
    return df  # Return the DataFrame

def create_state_change_plot(mecp2_results, gfp_results, cell_type, output_file):
    """Create a plot showing chromatin state changes between Mecp2 and GFP conditions.
    
    Args:
        mecp2_results (pd.DataFrame): A DataFrame containing the analysis results for the Mecp2 condition.
        gfp_results (pd.DataFrame): A DataFrame containing the analysis results for the GFP condition.
        cell_type (str): The cell type being analyzed.
        output_file (str): The path to save the plot to.
    
    Returns:
        pd.DataFrame: A DataFrame containing the merged results for Mecp2 and GFP conditions.
    """
    # Merge results on gene
    merged = pd.merge(  # Merge the Mecp2 and GFP results on the gene name
        mecp2_results[['gene', 'state', 'diff']].rename(columns={'state': 'Mecp2_state', 'diff': 'Mecp2_diff'}),  # Select the gene, state, and diff columns from the Mecp2 results and rename the columns
        gfp_results[['gene', 'state', 'diff']].rename(columns={'state': 'GFP_state', 'diff': 'GFP_diff'}),  # Select the gene, state, and diff columns from the GFP results and rename the columns
        on='gene', how='inner'  # Merge the DataFrames on the gene column using an inner join
    )
    
    if merged.empty:  # If no common genes were found
        logging.warning(f"No common genes found between Mecp2 and GFP for {cell_type}")  # Log a warning message
        return  # Return None
    
    # Count state transitions
    state_transitions = merged.groupby(['GFP_state', 'Mecp2_state']).size().reset_index(name='count')  # Group the merged DataFrame by the GFP state and Mecp2 state and count the number of occurrences of each combination
    
    # Create a pivot table for better visualization
    pivot_table = state_transitions.pivot(index='GFP_state', columns='Mecp2_state', values='count').fillna(0)  # Create a pivot table from the state transitions DataFrame
    
    # Create a heatmap of state transitions
    plt.figure(figsize=(10, 8))  # Create a new figure
    sns.heatmap(pivot_table, annot=True, fmt='g', cmap='viridis')  # Create a heatmap of the pivot table
    plt.title(f"Chromatin State Changes in {cell_type}: GFP vs Mecp2")  # Set the title of the plot
    plt.tight_layout()  # Adjust the plot layout
    plt.savefig(output_file)  # Save the plot to a file
    plt.close()  # Close the plot
    
    # Calculate percentages for each transition
    total = state_transitions['count'].sum()  # Calculate the total number of transitions
    state_transitions['percentage'] = state_transitions['count'] / total * 100  # Calculate the percentage of each transition
    
    # Save transition data
    transition_file = output_file.with_suffix('.csv')  # Create the path to the transition file
    state_transitions.to_csv(transition_file, index=False)  # Save the state transitions DataFrame to a CSV file
    
    logging.info(f"Created state change plot for {cell_type}: {output_file}")  # Log an info message
    return merged  # Return the merged DataFrame

def create_diff_comparison_plot(mecp2_results, gfp_results, cell_type, output_file):
    """Create a scatter plot comparing S3-S2S diffs between Mecp2 and GFP conditions.
    
    Args:
        mecp2_results (pd.DataFrame): A DataFrame containing the analysis results for the Mecp2 condition.
        gfp_results (pd.DataFrame): A DataFrame containing the analysis results for the GFP condition.
        cell_type (str): The cell type being analyzed.
        output_file (str): The path to save the plot to.
    
    Returns:
        pd.DataFrame: A DataFrame containing the merged results for Mecp2 and GFP conditions.
    """
    # Merge results on gene
    merged = pd.merge(  # Merge the Mecp2 and GFP results on the gene name
        mecp2_results[['gene', 'diff']].rename(columns={'diff': 'Mecp2_diff'}),  # Select the gene and diff columns from the Mecp2 results and rename the diff column
        gfp_results[['gene', 'diff']].rename(columns={'diff': 'GFP_diff'}),  # Select the gene and diff columns from the GFP results and rename the diff column
        on='gene', how='inner'  # Merge the DataFrames on the gene column using an inner join
    )
    
    if merged.empty:  # If no common genes were found
        logging.warning(f"No common genes found between Mecp2 and GFP for {cell_type}")  # Log a warning message
        return  # Return None
    
    # Log transform diffs for better visualization
    # merged['Mecp2_log_diff'] = np.log2(merged['Mecp2_diff'])  # Calculate the log2 of the Mecp2 diff
    # merged['GFP_log_diff'] = np.log2(merged['GFP_diff'])  # Calculate the log2 of the GFP diff
    merged['Mecp2_log_diff'] = merged['Mecp2_diff']  # do not calculate log2 of the Mecp2 diff
    merged['GFP_log_diff'] = merged['GFP_diff']  # do not calculate log2 of the GFP diff
    
    # Calculate diff change
    merged['log_diff_change'] = merged['Mecp2_log_diff'] - merged['GFP_log_diff']  # Calculate the difference between the log2 diffs
    
    # Create scatter plot
    plt.figure(figsize=(10, 8))  # Create a new figure
    
    # Plot diagonal line (no change)
    max_val = max(merged['Mecp2_log_diff'].max(), merged['GFP_log_diff'].max())  # Find the maximum log2 diff
    min_val = min(merged['Mecp2_log_diff'].min(), merged['GFP_log_diff'].min())  # Find the minimum log2 diff
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)  # Plot a diagonal line
    
    # Color points by diff change
    scatter = plt.scatter(  # Create a scatter plot
        merged['GFP_log_diff'],  # X-axis: GFP log2 diff
        merged['Mecp2_log_diff'],  # Y-axis: Mecp2 log2 diff
        c=merged['log_diff_change'],  # Color: log2 diff change
        cmap='coolwarm',  # Colormap: coolwarm
        alpha=0.7,  # Alpha: 0.7
        s=50  # Size: 50
    )
    
    # Add colorbar
    cbar = plt.colorbar(scatter)  # Add a colorbar
    cbar.set_label('Mecp2 - GFP diff')  # Set the colorbar label
    
    # Add labels and title
    plt.xlabel('S3 - S2S in GFP')  # Set the x-axis label
    plt.ylabel('S3 - S2S in Mecp2')  # Set the y-axis label
    plt.title(f"Chromatin State diff Changes in {cell_type}")  # Set the title
    
    # Add grid
    plt.grid(True, alpha=0.3)  # Add a grid
    
    # Save the figure
    plt.tight_layout()  # Adjust the plot layout
    plt.savefig(output_file)  # Save the plot to a file
    plt.close()  # Close the plot
    
    # Save data for further analysis
    diff_file = output_file.with_suffix('.csv')  # Create the path to the diff file
    merged.to_csv(diff_file, index=False)  # Save the merged DataFrame to a CSV file
    
    logging.info(f"Created diff comparison plot for {cell_type}: {output_file}")  # Log an info message
    return merged  # Return the merged DataFrame

def task3_analyze_exo_enriched_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache, chrom=None):
    """Task 3: Analyze chromatin state changes in genes where exogenous Mecp2 is enriched.
    
    Args:
        gene_list_files (list): A list of paths to the gene list files.
        gtf_file (str): The path to the GTF file.
        bigwig_dir (str): The path to the directory containing the BigWig files.
        output_dir (str): The path to the output directory.
        cache (ResultsCache): The cache object to use for caching results.
        chrom (str, optional): The specific chromosome to analyze. Defaults to None.
    
    Returns:
        pd.DataFrame: A DataFrame containing the combined results for all cell types.
    """
    # Create cache key
    cache_key = {  # Define the cache key
        'function': 'task3_analyze_exo_enriched_genes',
        'gene_list_files': [str(f) for f in gene_list_files],
        'gtf_file': str(gtf_file),
        'chrom': chrom
    }
    
    # Try to load from cache
    cached_result = cache.get(cache_key)  # Try to get the result from the cache
    if cached_result is not None:  # If the result was found in the cache
        return cached_result  # Return the cached result
    
    logging.info(f"\nPerforming Task 3: Analyzing chromatin state changes in exogenous Mecp2 enriched genes")  # Log an info message
    
    # Process each cell type with its corresponding gene list
    cell_types = ['Neu', 'NSC']  # Define the cell types to analyze
    results = {}  # Initialize an empty dictionary to store the results for each cell type
    
    for i, cell_type in enumerate(cell_types):  # Iterate over the cell types
        if i < len(gene_list_files):  # If there is a gene list file for this cell type
            gene_list_file = gene_list_files[i]  # Get the path to the gene list file
            logging.info(f"Processing {cell_type} with gene list: {gene_list_file}")  # Log an info message
            
            # Load gene list and coordinates
            genes = load_gene_list(gene_list_file)  # Load the list of genes from the gene list file
            if not genes:  # If no genes were loaded
                logging.error(f"No genes loaded for {cell_type} in Task 3, skipping")  # Log an error message
                continue  # Skip to the next cell type
                
            gene_coords = load_gene_coordinates(genes, gtf_file)  # Load the coordinates for the genes from the GTF file
            if not gene_coords:  # If no gene coordinates were loaded
                logging.error(f"No gene coordinates mapped for {cell_type} in Task 3, skipping")  # Log an error message
                continue  # Skip to the next cell type
            
            # Filter genes by chromosome if specified
            if chrom:  # If a specific chromosome was specified
                gene_coords = {gene: coords for gene, coords in gene_coords.items()  # Filter the gene coordinates to only include genes on the specified chromosome
                              if coords['chrom'] == chrom}
                logging.info(f"Filtered to {len(gene_coords)} genes on chromosome {chrom}")  # Log an info message
            
            # Get bigWig files for this cell type
            mecp2_s2s_files = get_bigwig_files(bigwig_dir, cell_type, "M2", "S2S")  # Get the BigWig files for the Mecp2 S2S condition
            mecp2_s3_files = get_bigwig_files(bigwig_dir, cell_type, "M2", "S3")  # Get the BigWig files for the Mecp2 S3 condition
            gfp_s2s_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S2S")  # Get the BigWig files for the GFP S2S condition
            gfp_s3_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S3")  # Get the BigWig files for the GFP S3 condition
            
            # Analyze gene phase in both conditions
            logging.info(f"Analyzing chromatin state in {cell_type} with Mecp2 overexpression...")  # Log an info message
            mecp2_results = analyze_gene_chromatin_state(gene_coords, mecp2_s2s_files, mecp2_s3_files, cache)  # Analyze the chromatin state for the Mecp2 condition
            
            # Analyze chromatin state in the current cell type with GFP control.
            logging.info(f"Analyzing chromatin state in {cell_type} with GFP control...")
            gfp_results = analyze_gene_chromatin_state(gene_coords, gfp_s2s_files, gfp_s3_files, cache)
            
            # Check if either Mecp2 or GFP results are empty. If so, skip this cell type.
            if mecp2_results.empty or gfp_results.empty:
                logging.warning(f"No results for {cell_type}, skipping")
                continue
            
            # Add cell type and condition information to the results DataFrames.
            mecp2_results['cell_type'] = cell_type
            mecp2_results['condition'] = 'Mecp2'
            gfp_results['cell_type'] = cell_type
            gfp_results['condition'] = 'GFP'
            
            # Save the individual results for Mecp2 and GFP to CSV files.
            mecp2_results.to_csv(output_dir / f"task3_{cell_type}_Mecp2_chromatin_state.csv", index=False)
            gfp_results.to_csv(output_dir / f"task3_{cell_type}_GFP_chromatin_state.csv", index=False)
            
            # Create comparison visualizations: state change plot and diff comparison plot.
            state_changes = create_state_change_plot(
                mecp2_results, gfp_results, 
                cell_type, 
                output_dir / f"task3_{cell_type}_state_changes.pdf"
            )
            
            diff_comparison = create_diff_comparison_plot(
                mecp2_results, gfp_results, 
                cell_type, 
                output_dir / f"task3_{cell_type}_diff_comparison.pdf"
            )
            
            # Store the results for this cell type in the 'results' dictionary.
            # This includes Mecp2 results, GFP results, and the state change comparison.
            results[cell_type] = {
                'mecp2': mecp2_results,
                'gfp': gfp_results,
                'comparison': state_changes if state_changes is not None else pd.DataFrame()
            }
    
    # After processing all cell types, check if we have results for both 'Neu' and 'NSC'.
    # If so, create cross-cell-type comparisons.
    if 'Neu' in results and 'NSC' in results:
        # Compare state changes between Neurons and NSCs using the previously generated comparison DataFrames.
        compare_state_changes_between_cell_types(
            results['Neu']['comparison'], 
            results['NSC']['comparison'],
            output_dir
        )
    
    # Combine all the results from different cell types and conditions into a single DataFrame.
    all_results = []
    for cell_type, cell_results in results.items():
        all_results.append(cell_results['mecp2'])
        all_results.append(cell_results['gfp'])
    
    combined_results = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    
    # If there are combined results, save them to a CSV file.
    if not combined_results.empty:
        combined_results.to_csv(output_dir / "task3_all_chromatin_states.csv", index=False)
    
    # Cache the combined results before returning them. This speeds up subsequent runs with the same parameters.
    # Store the combined results in the cache using the generated cache key.
    cache.set(cache_key, combined_results)
    # Return the combined results DataFrame for further analysis or reporting.
    return combined_results

def compare_state_changes_between_cell_types(neu_comparison, nsc_comparison, output_dir):
    """
    Compare chromatin state changes between Neurons (Neu) and Neural Stem Cells (NSCs).
    This function analyzes the differences in chromatin states between Mecp2 and GFP conditions
    in both cell types and generates a summary of these changes.

    Args:
        neu_comparison (pd.DataFrame): DataFrame containing the chromatin state comparison for Neurons.
        nsc_comparison (pd.DataFrame): DataFrame containing the chromatin state comparison for NSCs.
        output_dir (Path): Directory to save the comparison results and visualizations.
    """
    # Check if either of the input DataFrames is empty. If so, log a warning and exit the function.
    if neu_comparison.empty or nsc_comparison.empty:
        logging.warning("Cannot compare state changes between cell types: missing data")
        return
    
    # Calculate the percentage of genes that change state in Neurons.
    # 'GFP_state' and 'Mecp2_state' columns are compared to identify genes with state changes.
    neu_changed = neu_comparison[neu_comparison['GFP_state'] != neu_comparison['Mecp2_state']].shape[0]
    neu_total = neu_comparison.shape[0]
    neu_pct_changed = (neu_changed / neu_total * 100) if neu_total > 0 else 0
    
    # Calculate the percentage of genes that change state in NSCs.
    # Similar to Neurons, 'GFP_state' and 'Mecp2_state' columns are compared.
    nsc_changed = nsc_comparison[nsc_comparison['GFP_state'] != nsc_comparison['Mecp2_state']].shape[0]
    nsc_total = nsc_comparison.shape[0]
    nsc_pct_changed = (nsc_changed / nsc_total * 100) if nsc_total > 0 else 0
    
    # Create a bar chart comparing the percentage of genes that change state in Neurons and NSCs.
    plt.figure(figsize=(8, 6))
    bars = plt.bar(['Neurons', 'NSCs'], [neu_pct_changed, nsc_pct_changed])
    
    # Add value labels on top of the bars, displaying the percentage of genes changed.
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{height:.1f}%', ha='center', va='bottom')
    
    # Customize the plot with a title, y-axis label, y-axis limits, and gridlines.
    plt.title('Percentage of Genes with Chromatin State Changes')
    plt.ylabel('Percentage of Genes')
    plt.ylim(0, max(neu_pct_changed, nsc_pct_changed) * 1.2)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    # Save the generated figure to a PDF file in the specified output directory.
    plt.savefig(output_dir / "task3_cell_type_comparison.pdf")
    plt.close()
    
    # Create a summary table containing the total number of genes, number of changed genes,
    # and percentage of changed genes for both Neurons and NSCs.
    summary = pd.DataFrame({
        'Cell Type': ['Neurons', 'NSCs'],
        'Total Genes': [neu_total, nsc_total],
        'Changed Genes': [neu_changed, nsc_changed],
        'Percentage Changed': [neu_pct_changed, nsc_pct_changed]
    })
    
    # Save the summary table to a CSV file in the specified output directory.
    summary.to_csv(output_dir / "task3_cell_type_comparison.csv", index=False)
    
    # Log an informative message summarizing the comparison results.
    logging.info(f"Compared state changes between cell types: Neurons ({neu_pct_changed:.1f}%) vs NSCs ({nsc_pct_changed:.1f}%)")

def main():
    """
    Main function to execute the chromatin state change analysis for exogenous Mecp2 enriched genes.
    This function parses command-line arguments, sets up the logging, initializes directories and cache,
    loads gene lists and GTF file, and runs the main analysis function.
    """
    # Add command-line argument parsing using argparse.
    parser = argparse.ArgumentParser(description='Analyze chromatin state changes in exogenous Mecp2 enriched genes')
    parser.add_argument('--results_dir', default='../results', help='Directory for results')
    parser.add_argument('--chrom', help='Specific chromosome to analyze (e.g., "1", "X")')
    parser.add_argument('--gtf_file', default='../gencode.vM25.basic.annotation.gtf', help='GTF file with gene annotations')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()

    # Set up logging based on the debug flag.
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Log the start of the analysis and the specified results directory.
    logging.info("Starting analysis of exogenous Mecp2 enriched genes")
    logging.info(f"Initializing analysis with results_dir: {args.results_dir}")

    # Ensure the chromosome is properly formatted.
    chrom = None
    if args.chrom:
        # Format the chromosome name (add 'chr' prefix if not present).
        chrom = f"chr{args.chrom}" if not args.chrom.startswith('chr') else args.chrom
        logging.info(f"Processing single chromosome: {chrom}")

    # Set up directories for results, BigWig files, output, and cache.
    results_dir, bigwig_dir, output_dir, cache_dir = setup_directories(args.results_dir, chrom)
    
    # Initialize the ResultsCache to store intermediate results.
    cache = ResultsCache(cache_dir)
    
    # Log the directories being used for the analysis.
    logging.info(f"Using results directory: {results_dir}")
    logging.info(f"Using bigwig directory: {bigwig_dir}")
    logging.info(f"Output will be saved to: {output_dir}")
    logging.info(f"Cache directory: {cache_dir}")
    
    # Define the gene list files for Neurons and NSCs.
    gene_list_files = [
        Path("gene_lists") / "neu_enriched_gene_list_exo_vs_endo.txt",
        Path("gene_lists") / "nsc_enriched_gene_list_exo_vs_endo.txt"
    ]
    
    # Check if the gene list files exist. Exit if any file is not found.
    for file in gene_list_files:
        if not file.exists():
            logging.error(f"Gene list file not found: {file}")
            sys.exit(1)
    
    # Check if the GTF file exists. Exit if the file is not found.
    gtf_file = Path(args.gtf_file)
    if not gtf_file.exists():
        logging.error(f"GTF file not found: {gtf_file}")
        sys.exit(1)
    
    # Run the main analysis function: task3_analyze_exo_enriched_genes.
    task3_analyze_exo_enriched_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache, chrom)
    
    # Log a success message upon completion of the analysis.
    logging.info("Analysis completed successfully")

# Entry point of the script.
if __name__ == "__main__":
    main()