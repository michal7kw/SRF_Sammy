#!/usr/bin/env python3
# 3_analyze_exo_enriched_genes.py
# Script to analyze chromatin state changes in genes where exogenous Mecp2 is enriched (FC>2, exo vs endo)
# compared to the GFP control in both NSCs and Neurons

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
WINDOW_SIZE = 10000  # Window size for gene promoter analysis (10kb)
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
        output_dir = results_dir / "exo_enriched_analysis" / chrom
    else:
        output_dir = results_dir / "exo_enriched_analysis"
    
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

def create_state_change_plot(mecp2_results, gfp_results, cell_type, output_file):
    """Create a plot showing chromatin state changes between Mecp2 and GFP conditions."""
    # Merge results on gene
    merged = pd.merge(
        mecp2_results[['gene', 'state', 'ratio']].rename(columns={'state': 'Mecp2_state', 'ratio': 'Mecp2_ratio'}),
        gfp_results[['gene', 'state', 'ratio']].rename(columns={'state': 'GFP_state', 'ratio': 'GFP_ratio'}),
        on='gene', how='inner'
    )
    
    if merged.empty:
        logging.warning(f"No common genes found between Mecp2 and GFP for {cell_type}")
        return
    
    # Count state transitions
    state_transitions = merged.groupby(['GFP_state', 'Mecp2_state']).size().reset_index(name='count')
    
    # Create a pivot table for better visualization
    pivot_table = state_transitions.pivot(index='GFP_state', columns='Mecp2_state', values='count').fillna(0)
    
    # Create a heatmap of state transitions
    plt.figure(figsize=(10, 8))
    sns.heatmap(pivot_table, annot=True, fmt='g', cmap='viridis')
    plt.title(f"Chromatin State Changes in {cell_type}: GFP vs Mecp2")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    
    # Calculate percentages for each transition
    total = state_transitions['count'].sum()
    state_transitions['percentage'] = state_transitions['count'] / total * 100
    
    # Save transition data
    transition_file = output_file.with_suffix('.csv')
    state_transitions.to_csv(transition_file, index=False)
    
    logging.info(f"Created state change plot for {cell_type}: {output_file}")
    return merged

def create_ratio_comparison_plot(mecp2_results, gfp_results, cell_type, output_file):
    """Create a scatter plot comparing S3/S2S ratios between Mecp2 and GFP conditions."""
    # Merge results on gene
    merged = pd.merge(
        mecp2_results[['gene', 'ratio']].rename(columns={'ratio': 'Mecp2_ratio'}),
        gfp_results[['gene', 'ratio']].rename(columns={'ratio': 'GFP_ratio'}),
        on='gene', how='inner'
    )
    
    if merged.empty:
        logging.warning(f"No common genes found between Mecp2 and GFP for {cell_type}")
        return
    
    # Log transform ratios for better visualization
    merged['Mecp2_log_ratio'] = np.log2(merged['Mecp2_ratio'])
    merged['GFP_log_ratio'] = np.log2(merged['GFP_ratio'])
    
    # Calculate ratio change
    merged['log_ratio_change'] = merged['Mecp2_log_ratio'] - merged['GFP_log_ratio']
    
    # Create scatter plot
    plt.figure(figsize=(10, 8))
    
    # Plot diagonal line (no change)
    max_val = max(merged['Mecp2_log_ratio'].max(), merged['GFP_log_ratio'].max())
    min_val = min(merged['Mecp2_log_ratio'].min(), merged['GFP_log_ratio'].min())
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)
    
    # Color points by ratio change
    scatter = plt.scatter(
        merged['GFP_log_ratio'], 
        merged['Mecp2_log_ratio'],
        c=merged['log_ratio_change'], 
        cmap='coolwarm', 
        alpha=0.7,
        s=50
    )
    
    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('log2(Mecp2/GFP Ratio)')
    
    # Add labels and title
    plt.xlabel('log2(S3/S2S) in GFP')
    plt.ylabel('log2(S3/S2S) in Mecp2')
    plt.title(f"Chromatin State Ratio Changes in {cell_type}")
    
    # Add grid
    plt.grid(True, alpha=0.3)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    
    # Save data for further analysis
    ratio_file = output_file.with_suffix('.csv')
    merged.to_csv(ratio_file, index=False)
    
    logging.info(f"Created ratio comparison plot for {cell_type}: {output_file}")
    return merged

def task3_analyze_exo_enriched_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache, chrom=None):
    """Task 3: Analyze chromatin state changes in genes where exogenous Mecp2 is enriched."""
    # Create cache key
    cache_key = {
        'function': 'task3_analyze_exo_enriched_genes',
        'gene_list_files': [str(f) for f in gene_list_files],
        'gtf_file': str(gtf_file),
        'chrom': chrom
    }
    
    # Try to load from cache
    cached_result = cache.get(cache_key)
    if cached_result is not None:
        return cached_result
    
    logging.info(f"\nPerforming Task 3: Analyzing chromatin state changes in exogenous Mecp2 enriched genes")
    
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
                logging.error(f"No genes loaded for {cell_type} in Task 3, skipping")
                continue
                
            gene_coords = load_gene_coordinates(genes, gtf_file)
            if not gene_coords:
                logging.error(f"No gene coordinates mapped for {cell_type} in Task 3, skipping")
                continue
            
            # Filter genes by chromosome if specified
            if chrom:
                gene_coords = {gene: coords for gene, coords in gene_coords.items() 
                              if coords['chrom'] == chrom}
                logging.info(f"Filtered to {len(gene_coords)} genes on chromosome {chrom}")
            
            # Get bigWig files for this cell type
            mecp2_s2s_files = get_bigwig_files(bigwig_dir, cell_type, "M2", "S2S")
            mecp2_s3_files = get_bigwig_files(bigwig_dir, cell_type, "M2", "S3")
            gfp_s2s_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S2S")
            gfp_s3_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S3")
            
            # Analyze gene phase in both conditions
            logging.info(f"Analyzing chromatin state in {cell_type} with Mecp2 overexpression...")
            mecp2_results = analyze_gene_chromatin_state(gene_coords, mecp2_s2s_files, mecp2_s3_files, cache)
            
            logging.info(f"Analyzing chromatin state in {cell_type} with GFP control...")
            gfp_results = analyze_gene_chromatin_state(gene_coords, gfp_s2s_files, gfp_s3_files, cache)
            
            if mecp2_results.empty or gfp_results.empty:
                logging.warning(f"No results for {cell_type}, skipping")
                continue
            
            # Add cell type information
            mecp2_results['cell_type'] = cell_type
            mecp2_results['condition'] = 'Mecp2'
            gfp_results['cell_type'] = cell_type
            gfp_results['condition'] = 'GFP'
            
            # Save individual results
            mecp2_results.to_csv(output_dir / f"task3_{cell_type}_Mecp2_chromatin_state.csv", index=False)
            gfp_results.to_csv(output_dir / f"task3_{cell_type}_GFP_chromatin_state.csv", index=False)
            
            # Create comparison visualizations
            state_changes = create_state_change_plot(
                mecp2_results, gfp_results, 
                cell_type, 
                output_dir / f"task3_{cell_type}_state_changes.pdf"
            )
            
            ratio_comparison = create_ratio_comparison_plot(
                mecp2_results, gfp_results, 
                cell_type, 
                output_dir / f"task3_{cell_type}_ratio_comparison.pdf"
            )
            
            # Store results for this cell type
            results[cell_type] = {
                'mecp2': mecp2_results,
                'gfp': gfp_results,
                'comparison': state_changes if state_changes is not None else pd.DataFrame()
            }
    
    # If we have results for both cell types, create cross-cell-type comparisons
    if 'Neu' in results and 'NSC' in results:
        # Compare state changes between cell types
        compare_state_changes_between_cell_types(
            results['Neu']['comparison'], 
            results['NSC']['comparison'],
            output_dir
        )
    
    # Combine all results
    all_results = []
    for cell_type, cell_results in results.items():
        all_results.append(cell_results['mecp2'])
        all_results.append(cell_results['gfp'])
    
    combined_results = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    if not combined_results.empty:
        combined_results.to_csv(output_dir / "task3_all_chromatin_states.csv", index=False)
    
    # Cache the results before returning
    cache.set(cache_key, combined_results)
    return combined_results

def compare_state_changes_between_cell_types(neu_comparison, nsc_comparison, output_dir):
    """Compare chromatin state changes between Neurons and NSCs."""
    if neu_comparison.empty or nsc_comparison.empty:
        logging.warning("Cannot compare state changes between cell types: missing data")
        return
    
    # Calculate the percentage of genes that change state in each cell type
    neu_changed = neu_comparison[neu_comparison['GFP_state'] != neu_comparison['Mecp2_state']].shape[0]
    neu_total = neu_comparison.shape[0]
    neu_pct_changed = (neu_changed / neu_total * 100) if neu_total > 0 else 0
    
    nsc_changed = nsc_comparison[nsc_comparison['GFP_state'] != nsc_comparison['Mecp2_state']].shape[0]
    nsc_total = nsc_comparison.shape[0]
    nsc_pct_changed = (nsc_changed / nsc_total * 100) if nsc_total > 0 else 0
    
    # Create a bar chart comparing the percentage of genes that change state
    plt.figure(figsize=(8, 6))
    bars = plt.bar(['Neurons', 'NSCs'], [neu_pct_changed, nsc_pct_changed])
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{height:.1f}%', ha='center', va='bottom')
    
    plt.title('Percentage of Genes with Chromatin State Changes')
    plt.ylabel('Percentage of Genes')
    plt.ylim(0, max(neu_pct_changed, nsc_pct_changed) * 1.2)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_dir / "task3_cell_type_comparison.pdf")
    plt.close()
    
    # Create a summary table
    summary = pd.DataFrame({
        'Cell Type': ['Neurons', 'NSCs'],
        'Total Genes': [neu_total, nsc_total],
        'Changed Genes': [neu_changed, nsc_changed],
        'Percentage Changed': [neu_pct_changed, nsc_pct_changed]
    })
    
    summary.to_csv(output_dir / "task3_cell_type_comparison.csv", index=False)
    
    logging.info(f"Compared state changes between cell types: Neurons ({neu_pct_changed:.1f}%) vs NSCs ({nsc_pct_changed:.1f}%)")

def main():
    # Add command-line argument parsing
    parser = argparse.ArgumentParser(description='Analyze chromatin state changes in exogenous Mecp2 enriched genes')
    parser.add_argument('--results_dir', default='results', help='Directory for results')
    parser.add_argument('--chrom', help='Specific chromosome to analyze (e.g., "1", "X")')
    parser.add_argument('--gtf_file', default='gencode.vM25.basic.annotation.gtf', help='GTF file with gene annotations')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging.info("Starting analysis of exogenous Mecp2 enriched genes")
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
        Path("gene_lists") / "neu_enriched_gene_list_exo_vs_endo.txt",
        Path("gene_lists") / "nsc_enriched_gene_list_exo_vs_endo.txt"
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
    task3_analyze_exo_enriched_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache, chrom)
    
    logging.info("Analysis completed successfully")

if __name__ == "__main__":
    main() 