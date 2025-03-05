#!/usr/bin/env python3
# 2_analyze_endogenous_target_genes.py
"""
This script analyzes the chromatin phase of target genes of endogenous Mecp2 by comparing neurons (Neu) and neural stem cells (NSCs) in the GFP group.
It performs the following steps:
1.  Loads gene lists for Neurons and NSCs.
2.  Loads gene coordinates from a GTF file.
3.  Filters genes by chromosome if specified.
4.  Retrieves relevant BigWig files for S2S and S3 chromatin states.
5.  Analyzes the chromatin state of each gene by comparing S2S and S3 signals at the promoter region.
6.  Classifies genes into "Euchromatin", "Heterochromatin", or "Mixed" states based on the S3/S2S diff.
7.  Generates a phase comparison plot showing the distribution of chromatin states in Neurons and NSCs.
8.  Creates a heatmap comparing S3/S2S diffs between Neurons and NSCs.
9.  Compares gene states between cell types and generates a heatmap of state transitions.
10. Uses a caching mechanism to store and retrieve intermediate results, improving performance.
"""

"""
input gene list:

outputs: /results/endogenous_target_analysis/chr<n>

"""

import os  # For interacting with the operating system
import sys  # For system-specific parameters and functions
import logging  # For logging messages and debugging
import argparse  # For parsing command-line arguments
import hashlib  # For creating hash values for caching
import json  # For working with JSON data
import numpy as np  # For numerical opediffns
import pandas as pd  # For data manipulation and analysis
import matplotlib.pyplot as plt  # For plotting
import seaborn as sns  # For statistical data visualization
from pathlib import Path  # For working with file paths
import pyBigWig  # For accessing and manipulating BigWig files
import pickle  # For serializing and de-serializing Python objects
from collections import defaultdict  # For creating dictionaries with default values

# Constants
WINDOW_SIZE = 10000  # Window size for gene promoter analysis (5kb).  Defines the region around the transcription start site (TSS) to be analyzed.
MIN_DIFF_THRESHOLD = 2.0  # Minimum difference threshold for state changes.  Used to determine if the diff of S3/S2S signal is significantly different to classify chromatin state.

class ResultsCache:
    """
    Cache for storing computation results to avoid redundant calculations.
    This class uses a directory to store pickled results and a metadata file (JSON) to track the cached items.
    """
    
    def __init__(self, cache_dir):
        """
        Initialize the cache with the specified directory.
        
        Args:
            cache_dir (str): The directory where cached results will be stored.
        """
        self.cache_dir = Path(cache_dir)  # Convert cache directory to a Path object
        self.cache_dir.mkdir(parents=True, exist_ok=True)  # Create the cache directory if it doesn't exist
        self.metadata_file = self.cache_dir / "metadata.json"  # Define the metadata file path
        self.metadata = self._load_metadata()  # Load existing metadata or initialize an empty dictionary
        
    def _load_metadata(self):
        """
        Load metadata from file or initialize if it doesn't exist.
        
        Returns:
            dict: A dictionary containing the cache metadata.
        """
        if self.metadata_file.exists():  # Check if the metadata file exists
            try:
                with open(self.metadata_file, 'r') as f:  # Open the metadata file in read mode
                    return json.load(f)  # Load the JSON data from the file
            except (json.JSONDecodeError, IOError) as e:  # Handle potential errors during loading
                logging.warning(f"Failed to load cache metadata: {e}")  # Log a warning message
                return {}  # Return an empty dictionary if loading fails
        return {}  # Return an empty dictionary if the metadata file doesn't exist
    
    def _save_metadata(self):
        """Save metadata to file."""
        try:
            with open(self.metadata_file, 'w') as f:  # Open the metadata file in write mode
                json.dump(self.metadata, f)  # Save the metadata to the file as JSON
        except IOError as e:  # Handle potential errors during saving
            logging.warning(f"Failed to save cache metadata: {e}")  # Log a warning message
    
    def _compute_hash(self, key_dict):
        """
        Compute a hash for the cache key.
        
        Args:
            key_dict (dict): The dictionary to hash.
        
        Returns:
            str: The hexadecimal hash of the key dictionary.
        """
        # Add a script-specific prefix to ensure uniqueness across scripts
        key_dict_with_prefix = {"script": "task2_endogenous_targets", **key_dict}
        key_str = json.dumps(key_dict_with_prefix, sort_keys=True)  # Convert the key dictionary to a sorted JSON string
        return hashlib.md5(key_str.encode()).hexdigest()  # Compute the MD5 hash of the JSON string
    
    def get(self, key_dict, default=None):
        """
        Get a value from the cache using the key dictionary.
        
        Args:
            key_dict (dict): A dictionary representing the key for the cache entry.
            default (any, optional): The value to return if the key is not found in the cache. Defaults to None.
        
        Returns:
            any: The cached value if found, otherwise the default value.
        """
        key_hash = self._compute_hash(key_dict)  # Compute the hash of the key dictionary
        if key_hash in self.metadata:  # Check if the hash exists in the metadata
            cache_file = self.cache_dir / f"{key_hash}.pkl"  # Construct the file path for the cached data
            if cache_file.exists():  # Check if the cache file exists
                try:
                    with open(cache_file, 'rb') as f:  # Open the cache file in read binary mode
                        logging.debug(f"Cache hit for {key_hash}")  # Log a debug message indicating a cache hit
                        return pickle.load(f)  # Load the data from the file using pickle and return it
                except (pickle.PickleError, IOError) as e:  # Handle potential errors during loading
                    logging.warning(f"Failed to load cached data: {e}")  # Log a warning message
        return default  # Return the default value if the key is not found or loading fails
    
    def set(self, key_dict, value):
        """
        Set a value in the cache using the key dictionary.
        
        Args:
            key_dict (dict): A dictionary representing the key for the cache entry.
            value (any): The value to store in the cache.
        
        Returns:
            bool: True if the value was successfully cached, False otherwise.
        """
        key_hash = self._compute_hash(key_dict)  # Compute the hash of the key dictionary
        cache_file = self.cache_dir / f"{key_hash}.pkl"  # Construct the file path for the cached data
        try:
            with open(cache_file, 'wb') as f:  # Open the cache file in write binary mode
                pickle.dump(value, f)  # Dump the value to the file using pickle
            self.metadata[key_hash] = {  # Update the metadata with the new entry
                'key': key_dict,
                'file': str(cache_file)
            }
            self._save_metadata()  # Save the updated metadata
            logging.debug(f"Cached data for {key_hash}")  # Log a debug message indicating successful caching
            return True  # Return True to indicate success
        except (pickle.PickleError, IOError) as e:  # Handle potential errors during saving
            logging.warning(f"Failed to cache data: {e}")  # Log a warning message
            return False  # Return False to indicate failure

def setup_directories(results_dir, chrom=None):
    """
    Set up the necessary directories for analysis.
    
    Args:
        results_dir (str): The main directory for storing results.
        chrom (str, optional): The specific chromosome being analyzed. Defaults to None.
    
    Returns:
        tuple: A tuple containing the paths to the results directory, BigWig directory, output directory, and cache directory.
    """
    # Create main results directory
    results_dir = Path(results_dir)  # Convert results directory to a Path object
    results_dir.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist
    
    # Create subdirectories
    bigwig_dir = results_dir / "bigwig"  # Define the BigWig directory path
    
    # Create output directory with chromosome-specific subfolder if needed
    if chrom:  # If a specific chromosome is specified
        output_dir = results_dir / "endogenous_target_analysis" / chrom  # Create a chromosome-specific output directory
    else:  # If no specific chromosome is specified
        output_dir = results_dir / "endogenous_target_analysis"  # Create a general output directory
    
    output_dir.mkdir(parents=True, exist_ok=True)  # Create the output directory if it doesn't exist
    
    # Create script-specific cache directory
    cache_dir = results_dir / "cache" / "task2_cache"  # Define the script-specific cache directory path
    cache_dir.mkdir(parents=True, exist_ok=True)  # Create the cache directory if it doesn't exist
    
    return results_dir, bigwig_dir, output_dir, cache_dir  # Return the paths to the created directories

def get_bigwig_files(bigwig_dir, cell_type, condition, state):
    """
    Get a list of BigWig files for the specified parameters.
    
    Args:
        bigwig_dir (str): The directory containing the BigWig files.
        cell_type (str): The cell type to filter for (e.g., "Neu", "NSC").
        condition (str): The condition to filter for (e.g., "GFP").
        state (str): The state to filter for (e.g., "S2S", "S3").
    
    Returns:
        list: A list of Path objects representing the BigWig files that match the specified parameters.
    """
    bigwig_dir = Path(bigwig_dir)  # Convert BigWig directory to a Path object
    pattern = f"*{cell_type}*{condition}*{state}*.bw"  # Create a glob pattern to match the desired files
    files = list(bigwig_dir.glob(pattern))  # Find all files that match the pattern
    
    if not files:  # If no files were found
        logging.warning(f"No BigWig files found for {cell_type} {condition} {state} in {bigwig_dir}")  # Log a warning message
    else:  # If files were found
        logging.info(f"Found {len(files)} BigWig files for {cell_type} {condition} {state}")  # Log an info message
        for f in files:  # Iterate over the found files
            logging.debug(f"  - {f.name}")  # Log the name of each file as a debug message
    
    return files  # Return the list of found files

def compute_average_signal(bigwig_files, chrom, start, end):
    """
    Compute the average signal across multiple BigWig files for a genomic region.
    
    Args:
        bigwig_files (list): A list of Path objects representing the BigWig files to analyze.
        chrom (str): The chromosome of the genomic region.
        start (int): The start coordinate of the genomic region.
        end (int): The end coordinate of the genomic region.
    
    Returns:
        float: The average signal across all BigWig files for the specified genomic region, or None if no signal could be computed.
    """
    if not bigwig_files:  # If no BigWig files were provided
        return None  # Return None
    
    signals = []  # Initialize an empty list to store the average signal from each BigWig file
    for bw_file in bigwig_files:  # Iterate over the BigWig files
        try:
            bw = pyBigWig.open(str(bw_file))  # Open the BigWig file
            if chrom in bw.chroms():  # Check if the chromosome exists in the BigWig file
                # Ensure the region is within chromosome bounds
                chrom_size = bw.chroms()[chrom]  # Get the size of the chromosome
                valid_start = max(0, start)  # Ensure the start coordinate is not negative
                valid_end = min(end, chrom_size)  # Ensure the end coordinate is not greater than the chromosome size
                
                if valid_end > valid_start:  # If the valid region is not empty
                    # Get values and compute mean, ignoring NaN values
                    values = bw.values(chrom, valid_start, valid_end)  # Get the signal values for the region
                    values = np.array(values)  # Convert the values to a NumPy array
                    values = values[~np.isnan(values)]  # Remove NaN values
                    if len(values) > 0:  # If there are any valid values
                        signals.append(np.mean(values))  # Compute the mean of the values and append it to the list of signals
            bw.close()  # Close the BigWig file
        except Exception as e:  # Handle potential errors during processing
            logging.error(f"Error processing {bw_file}: {e}")  # Log an error message
    
    if signals:  # If any signals were collected
        return np.mean(signals)  # Return the average of the signals
    return None  # Return None if no signals were collected

def load_gene_list(gene_list_file):
    """
    Load a list of genes from a file.
    
    Args:
        gene_list_file (str): The path to the file containing the list of genes.
    
    Returns:
        list: A list of gene names.
    """
    genes = []  # Initialize an empty list to store the genes
    try:
        with open(gene_list_file, 'r') as f:  # Open the gene list file in read mode
            for line in f:  # Iterate over the lines in the file
                gene = line.strip()  # Remove leading/trailing whitespace from the line
                if gene:  # Skip empty lines
                    genes.append(gene)  # Add the gene to the list
        logging.info(f"Loaded {len(genes)} genes from {gene_list_file}")  # Log an info message
    except Exception as e:  # Handle potential errors during loading
        logging.error(f"Error loading gene list from {gene_list_file}: {e}")  # Log an error message
    
    return genes  # Return the list of genes

def load_gene_coordinates(genes, gtf_file):
    """
    Load gene coordinates from a GTF file for the specified genes.
    
    Args:
        genes (list): A list of gene names to load coordinates for.
        gtf_file (str): The path to the GTF file.
    
    Returns:
        dict: A dictionary mapping gene names to their coordinates (chromosome, start, end, strand, gene_id, gene_name).
    """
    gene_coords = {}  # Initialize an empty dictionary to store the gene coordinates
    
    try:
        # Read GTF file and extract gene coordinates
        with open(gtf_file, 'r') as f:  # Open the GTF file in read mode
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
                    attr = attr.strip()  # Remove leading/trailing whitespace
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
    """
    Analyze the chromatin state of genes by comparing S2S and S3 signals.
    
    Args:
        gene_coords (dict): A dictionary mapping gene names to their coordinates.
        s2s_files (list): A list of Path objects representing the S2S BigWig files.
        s3_files (list): A list of Path objects representing the S3 BigWig files.
        cache (ResultsCache): The cache object to use for caching results.
    
    Returns:
        pandas.DataFrame: A DataFrame containing the analysis results for each gene.
    """
    results = []  # Initialize an empty list to store the results
    
    for gene, coords in gene_coords.items():  # Iterate over the gene coordinates
        chrom = coords['chrom']  # Get the chromosome
        
        # Define promoter region (5kb upstream of TSS)
        if coords['strand'] == '+':  # If the gene is on the positive strand
            promoter_start = max(0, coords['start'] - WINDOW_SIZE)  # Calculate the promoter start coordinate
            promoter_end = coords['start']  # The promoter end coordinate is the TSS
        else:  # '-' strand  # If the gene is on the negative strand
            promoter_start = coords['end']  # The promoter start coordinate is the end of the gene
            promoter_end = coords['end'] + WINDOW_SIZE  # Calculate the promoter end coordinate
        
        # Create cache key for this computation
        cache_key = {  # Create a dictionary representing the cache key
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
            results.append(cached_result)  # Add the cached result to the list of results
            continue  # Skip to the next gene
        
        # Compute average signals
        s2s_signal = compute_average_signal(s2s_files, chrom, promoter_start, promoter_end)  # Compute the average S2S signal
        s3_signal = compute_average_signal(s3_files, chrom, promoter_start, promoter_end)  # Compute the average S3 signal
        
        # Determine chromatin state
        if s2s_signal is not None and s3_signal is not None:  # If both signals were computed
            diff = s3_signal - s2s_signal  # Calculate the diff of S3 to S2S signal
            
            # Classify based on diff
            if diff > MIN_DIFF_THRESHOLD:  # If the diff is significantly greater than 1
                state = "Heterochromatin"  # S3 dominant  # Classify as heterochromatin
            elif diff < -MIN_DIFF_THRESHOLD:  # If the diff is significantly less than 1
                state = "Euchromatin"  # S2S dominant  # Classify as euchromatin
            else:  # If the diff is not significantly different from 1
                state = "Mixed"  # No clear dominance  # Classify as mixed
            
            result = {  # Create a dictionary representing the result
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
            results.append(result)  # Add the result to the list of results
        else:  # If either signal could not be computed
            logging.warning(f"Could not compute signals for gene {gene} at {chrom}:{promoter_start}-{promoter_end}")  # Log a warning message
    
    # Convert results to DataFrame
    df = pd.DataFrame(results)  # Convert the list of results to a DataFrame
    return df  # Return the DataFrame

def create_phase_comparison_plot(neu_results, nsc_results, title, output_file):
    """
    Create a plot comparing chromatin phase distribution between Neurons and NSCs.
    
    Args:
        neu_results (pandas.DataFrame): DataFrame containing results for Neurons.
        nsc_results (pandas.DataFrame): DataFrame containing results for NSCs.
        title (str): Title of the plot.
        output_file (str): Path to save the plot.
    """
    # Count states in each cell type
    neu_counts = neu_results['state'].value_counts().to_dict()  # Count the occurrences of each state in Neurons
    nsc_counts = nsc_results['state'].value_counts().to_dict()  # Count the occurrences of each state in NSCs
    
    # Ensure all states are represented
    states = ['Euchromatin', 'Heterochromatin', 'Mixed']  # Define the possible states
    for state in states:  # Iterate over the states
        if state not in neu_counts:  # If the state is not in the Neuron counts
            neu_counts[state] = 0  # Set the count to 0
        if state not in nsc_counts:  # If the state is not in the NSC counts
            nsc_counts[state] = 0  # Set the count to 0
    
    # Create DataFrame for plotting
    plot_data = pd.DataFrame({  # Create a DataFrame for plotting
        'State': states * 2,
        'Count': [neu_counts[s] for s in states] + [nsc_counts[s] for s in states],
        'Cell Type': ['Neurons'] * len(states) + ['NSCs'] * len(states)
    })
    
    # Calculate percentages
    for cell_type in ['Neurons', 'NSCs']:  # Iterate over the cell types
        mask = plot_data['Cell Type'] == cell_type  # Create a mask for the current cell type
        total = plot_data.loc[mask, 'Count'].sum()  # Calculate the total count for the current cell type
        if total > 0:  # If the total count is greater than 0
            plot_data.loc[mask, 'Percentage'] = plot_data.loc[mask, 'Count'] / total * 100  # Calculate the percentage for each state
        else:  # If the total count is 0
            plot_data.loc[mask, 'Percentage'] = 0  # Set the percentage to 0
    
    # Create the plot
    plt.figure(figsize=(10, 6))  # Create a new figure
    
    # Bar plot
    ax = sns.barplot(x='State', y='Percentage', hue='Cell Type', data=plot_data)  # Create a bar plot
    
    # Add value labels on top of bars
    for p in ax.patches:  # Iterate over the bars
        ax.annotate(f'{p.get_height():.1f}%', 
                    (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='bottom', fontsize=10)  # Add the percentage label to the bar
    
    plt.title(f"Chromatin State Distribution: {title}")  # Set the title of the plot
    plt.ylabel("Percentage of Genes")  # Set the y-axis label
    plt.ylim(0, 100)  # Set the y-axis limits
    plt.tight_layout()  # Adjust the layout
    
    # Save the figure
    plt.savefig(output_file)  # Save the figure to a file
    plt.close()  # Close the figure
    
    logging.info(f"Created phase comparison plot: {output_file}")  # Log an info message

def create_diff_heatmap(neu_results, nsc_results, title, output_file):
    """
    Create a heatmap comparing S3/S2S diffs between Neurons and NSCs.
    
    Args:
        neu_results (pandas.DataFrame): DataFrame containing results for Neurons.
        nsc_results (pandas.DataFrame): DataFrame containing results for NSCs.
        title (str): Title of the heatmap.
        output_file (str): Path to save the heatmap.
    """
    # Merge the dataframes on gene
    merged = pd.merge(  # Merge the two DataFrames on the 'gene' column
        neu_results[['gene', 'diff']].rename(columns={'diff': 'Neurons_diff'}),  # Select the 'gene' and 'diff' columns from the Neuron results and rename the 'diff' column
        nsc_results[['gene', 'diff']].rename(columns={'diff': 'NSCs_diff'}),  # Select the 'gene' and 'diff' columns from the NSC results and rename the 'diff' column
        on='gene', how='inner'  # Perform an inner merge on the 'gene' column
    )
    
    if merged.empty:  # If the merged DataFrame is empty
        logging.warning("No common genes found between Neurons and NSCs for heatmap")  # Log a warning message
        return  # Return from the function
    
    # Log transform diffs for better visualization
    # merged['Neurons_log_diff'] = np.log2(merged['Neurons_diff'])  # Calculate the log2 of the Neuron diffs
    # merged['NSCs_log_diff'] = np.log2(merged['NSCs_diff'])  # Calculate the log2 of the NSC diffs
    merged['Neurons_log_diff'] = merged['Neurons_diff']  # do not use log2 of the Neuron diffs
    merged['NSCs_log_diff'] = merged['NSCs_diff']  # do not use log2 of the NSC diffs
    
    # Create a pivot table for the heatmap
    diff_data = merged[['gene', 'Neurons_log_diff', 'NSCs_log_diff']]  # Select the 'gene', 'Neurons_log_diff', and 'NSCs_log_diff' columns
    
    # Sort by the difference in log diffs
    diff_data['diff'] = diff_data['Neurons_log_diff'] - diff_data['NSCs_log_diff']  # Calculate the difference between the Neuron and NSC log diffs
    diff_data = diff_data.sort_values('diff')  # Sort the DataFrame by the difference
    
    # Select top and bottom genes for visualization (to avoid overcrowding)
    max_genes = 50  # Define the maximum number of genes to display
    if len(diff_data) > max_genes:  # If there are more genes than the maximum
        top_half = diff_data.iloc[-max_genes//2:]  # Select the top half of the genes
        bottom_half = diff_data.iloc[:max_genes//2]  # Select the bottom half of the genes
        diff_data = pd.concat([bottom_half, top_half])  # Concatenate the top and bottom halves
    
    # Prepare data for heatmap
    heatmap_data = diff_data.drop('diff', axis=1).set_index('gene')  # Drop the 'diff' column and set the 'gene' column as the index
    
    # Create the heatmap
    plt.figure(figsize=(12, max(8, len(heatmap_data) * 0.2)))  # Create a new figure
    
    # Define a diverging colormap centered at 0
    cmap = sns.diverging_palette(240, 10, as_cmap=True)  # Create a diverging colormap
    
    # Create heatmap
    sns.heatmap(heatmap_data, cmap=cmap, center=0, 
                yticklabels=True, xticklabels=True, 
                linewidths=0.5, cbar_kws={'label': 'S3-S2S diff'})  # Create the heatmap
    
    plt.title(f"Chromatin State diff Comparison: {title}")  # Set the title of the heatmap
    plt.tight_layout()  # Adjust the layout
    
    # Save the figure
    plt.savefig(output_file)  # Save the heatmap to a file
    plt.close()  # Close the figure
    
    logging.info(f"Created diff heatmap: {output_file}")  # Log an info message

def task2_analyze_endogenous_target_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache, chrom=None):
    """
    Task 2: Check phase of target genes of endogenous Mecp2.
    
    Args:
        gene_list_files (list): A list of paths to the gene list files for each cell type.
        gtf_file (str): The path to the GTF file.
        bigwig_dir (str): The directory containing the BigWig files.
        output_dir (str): The directory to save the output files.
        cache (ResultsCache): The cache object to use for caching results.
        chrom (str, optional): The specific chromosome to analyze. Defaults to None.
    
    Returns:
        pandas.DataFrame: A DataFrame containing the combined results for all cell types.
    """
    # Create cache key
    cache_key = {  # Create a dictionary representing the cache key
        'function': 'task2_analyze_endogenous_target_genes',
        'gene_list_files': [str(f) for f in gene_list_files],
        'gtf_file': str(gtf_file),
        'chrom': chrom
    }
    
    # Try to load from cache
    cached_result = cache.get(cache_key)  # Try to get the result from the cache
    if cached_result is not None:  # If the result was found in the cache
        return cached_result  # Return the cached result
    
    logging.info(f"\nPerforming Task 2: Analyzing phase of endogenous Mecp2 target genes")  # Log an info message
    
    # Process each cell type with its corresponding gene list
    cell_types = ['Neu', 'NSC']  # Define the cell types
    results = {}  # Initialize an empty dictionary to store the results for each cell type
    
    for i, cell_type in enumerate(cell_types):  # Iterate over the cell types
        if i < len(gene_list_files):  # If there is a gene list file for the current cell type
            gene_list_file = gene_list_files[i]  # Get the gene list file
            logging.info(f"Processing {cell_type} with gene list: {gene_list_file}")  # Log an info message
            
            # Load gene list and coordinates
            genes = load_gene_list(gene_list_file)  # Load the gene list
            if not genes:  # If no genes were loaded
                logging.error(f"No genes loaded for {cell_type} in Task 2, skipping")  # Log an error message
                continue  # Skip to the next cell type
                
            gene_coords = load_gene_coordinates(genes, gtf_file)  # Load the gene coordinates
            if not gene_coords:  # If no gene coordinates were loaded
                logging.error(f"No gene coordinates mapped for {cell_type} in Task 2, skipping")  # Log an error message
                continue  # Skip to the next cell type
            
            # Filter genes by chromosome if specified
            if chrom:  # If a specific chromosome was specified
                gene_coords = {gene: coords for gene, coords in gene_coords.items() 
                              if coords['chrom'] == chrom}  # Filter the gene coordinates to only include genes on the specified chromosome
                logging.info(f"Filtered to {len(gene_coords)} genes on chromosome {chrom}")  # Log an info message
            
            # Get GFP bigWig files for this cell type
            s2s_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S2S")  # Get the S2S BigWig files
            s3_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S3")  # Get the S3 BigWig files
            
            # Analyze gene phase
            logging.info(f"Analyzing chromatin state in {cell_type}...")  # Log an info message
            cell_results = analyze_gene_chromatin_state(gene_coords, s2s_files, s3_files, cache)  # Analyze the chromatin state of the genes
            
            if cell_results.empty:  # If no results were obtained
                logging.warning(f"No results for {cell_type}, skipping")  # Log a warning message
                continue  # Skip to the next cell type
                
            cell_results['cell_type'] = cell_type  # Add a column to the results indicating the cell type
            results[cell_type] = cell_results  # Store the results for the current cell type
            
            # Save individual results
            output_file = output_dir / f"task2_{cell_type}_target_genes_phase.csv"
            cell_results.to_csv(output_file, index=False)
            logging.info(f"Saved {cell_type} results to {output_file}")
    
    # If we have results for both cell types, create comparison visualizations
    if 'Neu' in results and 'NSC' in results:
        # Visualize phase distribution
        create_phase_comparison_plot(
            results['Neu'], results['NSC'], 
            "Endogenous Mecp2 Target Genes", 
            output_dir / "task2_phase_comparison.pdf"
        )
        
        # Create a heatmap of S3-S2S diffs
        create_diff_heatmap(
            results['Neu'], results['NSC'], 
            "Endogenous Mecp2 Target Genes",
            output_dir / "task2_diff_heatmap.pdf"
        )
        
        # Compare gene states between cell types
        compare_gene_states_between_cell_types(results['Neu'], results['NSC'], output_dir)
    
    # Combine all results
    combined_results = pd.concat(list(results.values()), ignore_index=True) if results else pd.DataFrame()
    if not combined_results.empty:
        combined_results.to_csv(output_dir / "task2_target_genes_phase.csv", index=False)
    
    # Cache the results before returning
    cache.set(cache_key, combined_results)
    return combined_results

def compare_gene_states_between_cell_types(neu_results, nsc_results, output_dir):
    """
    Compare chromatin states of the same genes between Neurons and NSCs.
    This function analyzes how the chromatin state of a gene changes (or stays the same)
    when comparing neurons and neural stem cells (NSCs).  It generates a heatmap and
    saves the state transition counts and percentages to files.

    Args:
        neu_results (pd.DataFrame): DataFrame containing chromatin state analysis results for Neurons.
                                     Must include 'gene', 'state', and 'diff' columns.
        nsc_results (pd.DataFrame): DataFrame containing chromatin state analysis results for NSCs.
                                     Must include 'gene', 'state', and 'diff' columns.
        output_dir (Path): Directory to save the output files (state transitions, percentages, and heatmap).
    """
    # Merge results on gene:  Find common genes between the two cell types based on the 'gene' column.
    # We rename the 'state' and 'diff' columns to distinguish between Neurons and NSCs.
    merged = pd.merge(
        neu_results[['gene', 'state', 'diff']].rename(columns={'state': 'Neurons_state', 'diff': 'Neurons_diff'}),
        nsc_results[['gene', 'state', 'diff']].rename(columns={'state': 'NSCs_state', 'diff': 'NSCs_diff'}),
        on='gene',  # Merge based on the 'gene' column
        how='inner'  # Perform an inner merge to keep only common genes
    )
    
    # Check if any common genes were found
    if merged.empty:
        logging.warning("No common genes found between Neurons and NSCs for comparison")
        return  # Exit the function if no common genes are found
    
    # Count state transitions:  Determine how many genes transition from one chromatin state in Neurons
    # to another chromatin state in NSCs.  For example, how many genes are in an "active" state in
    # Neurons but a "repressed" state in NSCs?
    state_transitions = merged.groupby(['Neurons_state', 'NSCs_state']).size().reset_index(name='count')
    
    # Create a pivot table for better visualization:  Reshape the state transition counts into a pivot table
    # for easier interpretation and visualization as a heatmap.  The rows represent Neuron states,
    # the columns represent NSC states, and the values represent the number of genes transitioning
    # between those states.
    pivot_table = state_transitions.pivot(index='Neurons_state', columns='NSCs_state', values='count').fillna(0)
    
    # Save the transition counts to a CSV file
    pivot_table.to_csv(output_dir / "task2_state_transitions.csv")
    
    # Create a heatmap of state transitions:  Visualize the state transitions using a heatmap.
    # The heatmap provides a color-coded representation of the number of genes transitioning
    # between different chromatin states.
    plt.figure(figsize=(10, 8))  # Set the figure size
    sns.heatmap(pivot_table, annot=True, fmt='g', cmap='viridis')  # Create the heatmap with annotations, integer formatting, and a colormap
    plt.title("Chromatin State Transitions: Neurons vs NSCs")  # Set the title of the heatmap
    plt.tight_layout()  # Adjust the layout to prevent labels from overlapping
    plt.savefig(output_dir / "task2_state_transitions.pdf")  # Save the heatmap to a PDF file
    plt.close()  # Close the figure to free up memory
    
    # Calculate percentages of genes in each state for each cell type:  Determine the percentage of genes
    # in each chromatin state for both Neurons and NSCs.  This provides a normalized view of the
    # distribution of chromatin states in each cell type.
    neu_state_counts = merged['Neurons_state'].value_counts(normalize=True) * 100  # Calculate percentages for Neurons
    nsc_state_counts = merged['NSCs_state'].value_counts(normalize=True) * 100  # Calculate percentages for NSCs
    
    # Combine into a DataFrame:  Create a DataFrame to store the state percentages for both cell types.
    state_percentages = pd.DataFrame({
        'Neurons': neu_state_counts,
        'NSCs': nsc_state_counts
    }).fillna(0).reset_index().rename(columns={'index': 'State'})
    
    # Save state percentages to a CSV file
    state_percentages.to_csv(output_dir / "task2_state_percentages.csv", index=False)
    
    logging.info(f"Compared gene states between cell types, found {len(merged)} common genes")

def main():
    """
    Main function to execute the endogenous Mecp2 target gene analysis.
    This function parses command-line arguments, sets up the logging,
    initializes directories and cache, and runs the main analysis task.
    """
    # Add command-line argument parsing
    parser = argparse.ArgumentParser(description='Analyze endogenous Mecp2 target genes')
    parser.add_argument('--results_dir', default='../results', help='Directory for results')
    parser.add_argument('--chrom', help='Specific chromosome to analyze (e.g., "1", "X")')
    parser.add_argument('--gtf_file', default='../gencode.vM25.basic.annotation.gtf', help='GTF file with gene annotations')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging.info("Starting endogenous Mecp2 target gene analysis")
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
    logging.info(f"Using output directory: {output_dir}")
    logging.info(f"Using cache directory: {cache_dir}")
    
    # Define gene list files
    gene_list_files = [
        Path("../gene_lists") / "endo_gene_list_neu.txt",
        Path("../gene_lists") / "endo_gene_list_nsc.txt"
    ]
    # Check if gene list files exist
    for file in gene_list_files:
        if not file.exists():
            logging.error(f"Gene list file not found: {file}")
            sys.exit(1)
    
    # Check if GTF file exists
    gtf_file = Path(args.gtf_file)
    if not gtf_file.exists():
        logging.error(f"GTF file not found: {gtf_file}")
        sys.exit(1)
    
    # Run the analysis
    task2_analyze_endogenous_target_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache, chrom)
    
    logging.info("Analysis completed successfully")

if __name__ == "__main__":
    main()