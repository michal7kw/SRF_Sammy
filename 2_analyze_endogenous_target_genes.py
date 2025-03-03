#!/usr/bin/env python3
# 2_analyze_endogenous_target_genes.py
# Script to analyze the phase of target genes of endogenous Mecp2 by comparing neurons and NSCs in GFP group

import os
import sys
import logging
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pyBigWig
import pickle
from collections import defaultdict

# Constants
WINDOW_SIZE = 10000  # Window size for gene promoter analysis (5kb)
MIN_DIFF_THRESHOLD = 2.0  # Minimum difference threshold for state changes

class ResultsCache:
    """Cache for storing computation results to avoid redundant calculations."""
    
    def __init__(self, cache_dir):
        """Initialize the cache with the specified directory."""
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_file = self.cache_dir / "metadata.json"
        self.metadata = self._load_metadata()
        
    def _load_metadata(self):
        """Load metadata from file or initialize if it doesn't exist."""
        if self.metadata_file.exists():
            try:
                with open(self.metadata_file, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                logging.warning(f"Failed to load cache metadata: {e}")
                return {}
        return {}
    
    def _save_metadata(self):
        """Save metadata to file."""
        try:
            with open(self.metadata_file, 'w') as f:
                json.dump(self.metadata, f)
        except IOError as e:
            logging.warning(f"Failed to save cache metadata: {e}")
    
    def _compute_hash(self, key_dict):
        """Compute a hash for the given key dictionary."""
        key_str = json.dumps(key_dict, sort_keys=True)
        return hashlib.md5(key_str.encode()).hexdigest()
    
    def get(self, key_dict, default=None):
        """Get a value from the cache using the key dictionary."""
        key_hash = self._compute_hash(key_dict)
        if key_hash in self.metadata:
            cache_file = self.cache_dir / f"{key_hash}.pkl"
            if cache_file.exists():
                try:
                    with open(cache_file, 'rb') as f:
                        logging.debug(f"Cache hit for {key_hash}")
                        return pickle.load(f)
                except (pickle.PickleError, IOError) as e:
                    logging.warning(f"Failed to load cached data: {e}")
        return default
    
    def set(self, key_dict, value):
        """Set a value in the cache using the key dictionary."""
        key_hash = self._compute_hash(key_dict)
        cache_file = self.cache_dir / f"{key_hash}.pkl"
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(value, f)
            self.metadata[key_hash] = {
                'key': key_dict,
                'file': str(cache_file)
            }
            self._save_metadata()
            logging.debug(f"Cached data for {key_hash}")
            return True
        except (pickle.PickleError, IOError) as e:
            logging.warning(f"Failed to cache data: {e}")
            return False

def setup_directories(results_dir, chrom=None):
    """Set up the necessary directories for analysis."""
    # Create main results directory
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    bigwig_dir = results_dir / "bigwig"
    
    # Create output directory with chromosome-specific subfolder if needed
    if chrom:
        output_dir = results_dir / "endogenous_target_analysis" / chrom
    else:
        output_dir = results_dir / "endogenous_target_analysis"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create cache directory
    cache_dir = results_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    return results_dir, bigwig_dir, output_dir, cache_dir

def get_bigwig_files(bigwig_dir, cell_type, condition, state):
    """Get a list of BigWig files for the specified parameters."""
    bigwig_dir = Path(bigwig_dir)
    pattern = f"*{cell_type}*{condition}*{state}*.bw"
    files = list(bigwig_dir.glob(pattern))
    
    if not files:
        logging.warning(f"No BigWig files found for {cell_type} {condition} {state} in {bigwig_dir}")
    else:
        logging.info(f"Found {len(files)} BigWig files for {cell_type} {condition} {state}")
        for f in files:
            logging.debug(f"  - {f.name}")
    
    return files

def compute_average_signal(bigwig_files, chrom, start, end):
    """Compute the average signal across multiple BigWig files for a genomic region."""
    if not bigwig_files:
        return None
    
    signals = []
    for bw_file in bigwig_files:
        try:
            bw = pyBigWig.open(str(bw_file))
            if chrom in bw.chroms():
                # Ensure the region is within chromosome bounds
                chrom_size = bw.chroms()[chrom]
                valid_start = max(0, start)
                valid_end = min(end, chrom_size)
                
                if valid_end > valid_start:
                    # Get values and compute mean, ignoring NaN values
                    values = bw.values(chrom, valid_start, valid_end)
                    values = np.array(values)
                    values = values[~np.isnan(values)]  # Remove NaN values
                    if len(values) > 0:
                        signals.append(np.mean(values))
            bw.close()
        except Exception as e:
            logging.error(f"Error processing {bw_file}: {e}")
    
    if signals:
        return np.mean(signals)
    return None

def load_gene_list(gene_list_file):
    """Load a list of genes from a file."""
    genes = []
    try:
        with open(gene_list_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:  # Skip empty lines
                    genes.append(gene)
        logging.info(f"Loaded {len(genes)} genes from {gene_list_file}")
    except Exception as e:
        logging.error(f"Error loading gene list from {gene_list_file}: {e}")
    
    return genes

def load_gene_coordinates(genes, gtf_file):
    """Load gene coordinates from a GTF file for the specified genes."""
    gene_coords = {}
    
    try:
        # Read GTF file and extract gene coordinates
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):  # Skip header lines
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                
                # Extract gene ID from attributes
                attributes = fields[8]
                gene_id = None
                gene_name = None
                
                # Parse attributes to find gene_id and gene_name
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1] if '"' in attr else attr.split('=')[1]
                    elif attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1] if '"' in attr else attr.split('=')[1]
                
                # Check if this gene is in our list
                if gene_name in genes or gene_id in genes:
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    
                    # Use gene_name if available, otherwise use gene_id
                    gene_key = gene_name if gene_name in genes else gene_id
                    
                    # Store gene coordinates
                    gene_coords[gene_key] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'gene_id': gene_id,
                        'gene_name': gene_name
                    }
        
        logging.info(f"Mapped coordinates for {len(gene_coords)} out of {len(genes)} genes")
    except Exception as e:
        logging.error(f"Error loading gene coordinates from {gtf_file}: {e}")
    
    return gene_coords

def analyze_gene_chromatin_state(gene_coords, s2s_files, s3_files, cache):
    """Analyze the chromatin state of genes by comparing S2S and S3 signals."""
    results = []
    
    for gene, coords in gene_coords.items():
        chrom = coords['chrom']
        
        # Define promoter region (5kb upstream of TSS)
        if coords['strand'] == '+':
            promoter_start = max(0, coords['start'] - WINDOW_SIZE)
            promoter_end = coords['start']
        else:  # '-' strand
            promoter_start = coords['end']
            promoter_end = coords['end'] + WINDOW_SIZE
        
        # Create cache key for this computation
        cache_key = {
            'function': 'analyze_gene_chromatin_state',
            'gene': gene,
            'chrom': chrom,
            'promoter_start': promoter_start,
            'promoter_end': promoter_end,
            's2s_files': [str(f) for f in s2s_files],
            's3_files': [str(f) for f in s3_files]
        }
        
        # Try to get from cache
        cached_result = cache.get(cache_key)
        if cached_result is not None:
            results.append(cached_result)
            continue
        
        # Compute average signals
        s2s_signal = compute_average_signal(s2s_files, chrom, promoter_start, promoter_end)
        s3_signal = compute_average_signal(s3_files, chrom, promoter_start, promoter_end)
        
        # Determine chromatin state
        if s2s_signal is not None and s3_signal is not None:
            ratio = s3_signal / s2s_signal if s2s_signal > 0 else float('inf')
            
            # Classify based on ratio
            if ratio > 1 + MIN_DIFF_THRESHOLD:
                state = "Heterochromatin"  # S3 dominant
            elif ratio < 1 - MIN_DIFF_THRESHOLD:
                state = "Euchromatin"  # S2S dominant
            else:
                state = "Mixed"  # No clear dominance
            
            result = {
                'gene': gene,
                'chrom': chrom,
                'start': coords['start'],
                'end': coords['end'],
                'strand': coords['strand'],
                'promoter_start': promoter_start,
                'promoter_end': promoter_end,
                's2s_signal': s2s_signal,
                's3_signal': s3_signal,
                'ratio': ratio,
                'state': state
            }
            
            # Cache and append result
            cache.set(cache_key, result)
            results.append(result)
        else:
            logging.warning(f"Could not compute signals for gene {gene} at {chrom}:{promoter_start}-{promoter_end}")
    
    # Convert results to DataFrame
    df = pd.DataFrame(results)
    return df

def create_phase_comparison_plot(neu_results, nsc_results, title, output_file):
    """Create a plot comparing chromatin phase distribution between Neurons and NSCs."""
    # Count states in each cell type
    neu_counts = neu_results['state'].value_counts().to_dict()
    nsc_counts = nsc_results['state'].value_counts().to_dict()
    
    # Ensure all states are represented
    states = ['Euchromatin', 'Heterochromatin', 'Mixed']
    for state in states:
        if state not in neu_counts:
            neu_counts[state] = 0
        if state not in nsc_counts:
            nsc_counts[state] = 0
    
    # Create DataFrame for plotting
    plot_data = pd.DataFrame({
        'State': states * 2,
        'Count': [neu_counts[s] for s in states] + [nsc_counts[s] for s in states],
        'Cell Type': ['Neurons'] * len(states) + ['NSCs'] * len(states)
    })
    
    # Calculate percentages
    for cell_type in ['Neurons', 'NSCs']:
        mask = plot_data['Cell Type'] == cell_type
        total = plot_data.loc[mask, 'Count'].sum()
        if total > 0:
            plot_data.loc[mask, 'Percentage'] = plot_data.loc[mask, 'Count'] / total * 100
        else:
            plot_data.loc[mask, 'Percentage'] = 0
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    
    # Bar plot
    ax = sns.barplot(x='State', y='Percentage', hue='Cell Type', data=plot_data)
    
    # Add value labels on top of bars
    for p in ax.patches:
        ax.annotate(f'{p.get_height():.1f}%', 
                    (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='bottom', fontsize=10)
    
    plt.title(f"Chromatin State Distribution: {title}")
    plt.ylabel("Percentage of Genes")
    plt.ylim(0, 100)
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file)
    plt.close()
    
    logging.info(f"Created phase comparison plot: {output_file}")

def create_ratio_heatmap(neu_results, nsc_results, title, output_file):
    """Create a heatmap comparing S3/S2S ratios between Neurons and NSCs."""
    # Merge the dataframes on gene
    merged = pd.merge(
        neu_results[['gene', 'ratio']].rename(columns={'ratio': 'Neurons_ratio'}),
        nsc_results[['gene', 'ratio']].rename(columns={'ratio': 'NSCs_ratio'}),
        on='gene', how='inner'
    )
    
    if merged.empty:
        logging.warning("No common genes found between Neurons and NSCs for heatmap")
        return
    
    # Log transform ratios for better visualization
    merged['Neurons_log_ratio'] = np.log2(merged['Neurons_ratio'])
    merged['NSCs_log_ratio'] = np.log2(merged['NSCs_ratio'])
    
    # Create a pivot table for the heatmap
    ratio_data = merged[['gene', 'Neurons_log_ratio', 'NSCs_log_ratio']]
    
    # Sort by the difference in log ratios
    ratio_data['diff'] = ratio_data['Neurons_log_ratio'] - ratio_data['NSCs_log_ratio']
    ratio_data = ratio_data.sort_values('diff')
    
    # Select top and bottom genes for visualization (to avoid overcrowding)
    max_genes = 50
    if len(ratio_data) > max_genes:
        top_half = ratio_data.iloc[-max_genes//2:]
        bottom_half = ratio_data.iloc[:max_genes//2]
        ratio_data = pd.concat([bottom_half, top_half])
    
    # Prepare data for heatmap
    heatmap_data = ratio_data.drop('diff', axis=1).set_index('gene')
    
    # Create the heatmap
    plt.figure(figsize=(12, max(8, len(heatmap_data) * 0.2)))
    
    # Define a diverging colormap centered at 0
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    
    # Create heatmap
    sns.heatmap(heatmap_data, cmap=cmap, center=0, 
                yticklabels=True, xticklabels=True, 
                linewidths=0.5, cbar_kws={'label': 'log2(S3/S2S Ratio)'})
    
    plt.title(f"Chromatin State Ratio Comparison: {title}")
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file)
    plt.close()
    
    logging.info(f"Created ratio heatmap: {output_file}")

def task2_analyze_endogenous_target_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache, chrom=None):
    """Task 2: Check phase of target genes of endogenous Mecp2."""
    # Create cache key
    cache_key = {
        'function': 'task2_analyze_endogenous_target_genes',
        'gene_list_files': [str(f) for f in gene_list_files],
        'gtf_file': str(gtf_file),
        'chrom': chrom
    }
    
    # Try to load from cache
    cached_result = cache.get(cache_key)
    if cached_result is not None:
        return cached_result
    
    logging.info(f"\nPerforming Task 2: Analyzing phase of endogenous Mecp2 target genes")
    
    # Process each cell type with its corresponding gene list
    cell_types = ['Neu', 'NSC']
    results = {}
    
    for i, cell_type in enumerate(cell_types):
        if i < len(gene_list_files):
            gene_list_file = gene_list_files[i]
            logging.info(f"Processing {cell_type} with gene list: {gene_list_file}")
            
            # Load gene list and coordinates
            genes = load_gene_list(gene_list_file)
            if not genes:
                logging.error(f"No genes loaded for {cell_type} in Task 2, skipping")
                continue
                
            gene_coords = load_gene_coordinates(genes, gtf_file)
            if not gene_coords:
                logging.error(f"No gene coordinates mapped for {cell_type} in Task 2, skipping")
                continue
            
            # Filter genes by chromosome if specified
            if chrom:
                gene_coords = {gene: coords for gene, coords in gene_coords.items() 
                              if coords['chrom'] == chrom}
                logging.info(f"Filtered to {len(gene_coords)} genes on chromosome {chrom}")
            
            # Get GFP bigWig files for this cell type
            s2s_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S2S")
            s3_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S3")
            
            # Analyze gene phase
            logging.info(f"Analyzing chromatin state in {cell_type}...")
            cell_results = analyze_gene_chromatin_state(gene_coords, s2s_files, s3_files, cache)
            
            if cell_results.empty:
                logging.warning(f"No results for {cell_type}, skipping")
                continue
                
            cell_results['cell_type'] = cell_type
            results[cell_type] = cell_results
            
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
        
        # Create a heatmap of S3/S2S ratios
        create_ratio_heatmap(
            results['Neu'], results['NSC'], 
            "Endogenous Mecp2 Target Genes",
            output_dir / "task2_ratio_heatmap.pdf"
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
    """Compare chromatin states of the same genes between Neurons and NSCs."""
    # Merge results on gene
    merged = pd.merge(
        neu_results[['gene', 'state', 'ratio']].rename(columns={'state': 'Neurons_state', 'ratio': 'Neurons_ratio'}),
        nsc_results[['gene', 'state', 'ratio']].rename(columns={'state': 'NSCs_state', 'ratio': 'NSCs_ratio'}),
        on='gene', how='inner'
    )
    
    if merged.empty:
        logging.warning("No common genes found between Neurons and NSCs for comparison")
        return
    
    # Count state transitions
    state_transitions = merged.groupby(['Neurons_state', 'NSCs_state']).size().reset_index(name='count')
    
    # Create a pivot table for better visualization
    pivot_table = state_transitions.pivot(index='Neurons_state', columns='NSCs_state', values='count').fillna(0)
    
    # Save the transition counts
    pivot_table.to_csv(output_dir / "task2_state_transitions.csv")
    
    # Create a heatmap of state transitions
    plt.figure(figsize=(10, 8))
    sns.heatmap(pivot_table, annot=True, fmt='g', cmap='viridis')
    plt.title("Chromatin State Transitions: Neurons vs NSCs")
    plt.tight_layout()
    plt.savefig(output_dir / "task2_state_transitions.pdf")
    plt.close()
    
    # Calculate percentages of genes in each state for each cell type
    neu_state_counts = merged['Neurons_state'].value_counts(normalize=True) * 100
    nsc_state_counts = merged['NSCs_state'].value_counts(normalize=True) * 100
    
    # Combine into a DataFrame
    state_percentages = pd.DataFrame({
        'Neurons': neu_state_counts,
        'NSCs': nsc_state_counts
    }).fillna(0).reset_index().rename(columns={'index': 'State'})
    
    # Save state percentages
    state_percentages.to_csv(output_dir / "task2_state_percentages.csv", index=False)
    
    logging.info(f"Compared gene states between cell types, found {len(merged)} common genes")

def main():
    # Add command-line argument parsing
    parser = argparse.ArgumentParser(description='Analyze endogenous Mecp2 target genes')
    parser.add_argument('--results_dir', default='results', help='Directory for results')
    parser.add_argument('--chrom', help='Specific chromosome to analyze (e.g., "1", "X")')
    parser.add_argument('--gtf_file', default='reference/mm10.refGene.gtf', help='GTF file with gene annotations')
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
    logging.info(f"Output will be saved to: {output_dir}")
    logging.info(f"Cache directory: {cache_dir}")
    
    # Define gene list files
    gene_list_files = [
        Path("gene_lists") / "endo_gene_list_neu.txt",
        Path("gene_lists") / "endo_gene_list_nsc.txt"
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