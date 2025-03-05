#!/usr/bin/env python3
# 4_analyze_endogenous_enriched_genes.py
# Script to analyze the phase of genes where endogenous Mecp2 is enriched (FC>2, Neurons vs NSCs) in GFP-treated samples

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
FOLD_CHANGE_THRESHOLD = 2.0  # Fold change threshold for enrichment

class ResultsCache:
    def __init__(self, cache_dir):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_file = self.cache_dir / "metadata.json"
        self.metadata = self._load_metadata()
        
    def _load_metadata(self):
        if self.metadata_file.exists():
            try:
                with open(self.metadata_file, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                logging.warning(f"Failed to load cache metadata: {e}")
                return {}
        return {}
    
    def _save_metadata(self):
        try:
            with open(self.metadata_file, 'w') as f:
                json.dump(self.metadata, f)
        except IOError as e:
            logging.warning(f"Failed to save cache metadata: {e}")
    
    def _compute_hash(self, key_dict):
        """Compute a hash for the cache key."""
        # Add a script-specific prefix to ensure uniqueness across scripts
        key_dict_with_prefix = {"script": "task4_endo_enriched_nsc_vs_neu", **key_dict}
        key_str = json.dumps(key_dict_with_prefix, sort_keys=True)
        return hashlib.md5(key_str.encode()).hexdigest()
    
    def get(self, key_dict, default=None):
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
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    
    bigwig_dir = Path("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_1b/aligned/bigwig_files")
    
    if chrom:
        output_dir = results_dir / "endogenous_enriched_analysis" / chrom
    else:
        output_dir = results_dir / "endogenous_enriched_analysis"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create script-specific cache directory
    cache_dir = results_dir / "cache" / "task4_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    return results_dir, bigwig_dir, output_dir, cache_dir

def get_bigwig_files(bigwig_dir, cell_type, condition, state):
    """Get a list of BigWig files for the specified parameters.
    
    Args:
        bigwig_dir (str): The directory containing the BigWig files.
        cell_type (str): The cell type (e.g., "Neu", "NSC", "GFP").
        condition (str): The experimental condition (e.g., "M2", "GFP").
        state (str): The chromatin state (e.g., "S2S", "S3").
    
    Returns:
        list: A list of Path objects representing the BigWig files that match the specified parameters.
    """
    bigwig_dir = Path(bigwig_dir)
    
    # The BigWig files follow the pattern NeuV[1-3].bw or NSCV[1-3].bw
    if cell_type == 'Neu':
        pattern = "NeuV[1-3].bw"
    elif cell_type == 'NSC':
        pattern = "NSCv[1-3].bw"
    else:
        pattern = "*.bw"
    
    # List all matching files in the directory
    files = list(bigwig_dir.glob(pattern))
    
    if not files:
        logging.warning(f"No BigWig files found for {cell_type} {condition} {state} in {bigwig_dir}")
    else:
        logging.info(f"Found {len(files)} BigWig files for {cell_type} {condition} {state}")
        for f in files:
            logging.debug(f"  - {f.name}")
    
    return files

def compute_average_signal(bigwig_files, chrom, start, end):
    if not bigwig_files:
        return None
    
    signals = []
    for bw_file in bigwig_files:
        try:
            bw = pyBigWig.open(str(bw_file))
            if not bw.chroms().get(chrom):
                logging.warning(f"Chromosome {chrom} not found in {bw_file}")
                continue
                
            values = bw.values(chrom, start, end)
            if values:
                # Replace None values with 0 and compute mean
                values = [v if v is not None else 0 for v in values]
                signals.append(np.mean(values))
            bw.close()
        except Exception as e:
            logging.error(f"Error processing {bw_file}: {e}")
            continue
    
    return np.mean(signals) if signals else None

def load_enriched_genes(gene_list_file):
    """Load the list of enriched genes."""
    try:
        with open(gene_list_file, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        logging.info(f"Loaded {len(genes)} enriched genes from {gene_list_file}")
        return genes
    except Exception as e:
        logging.error(f"Error loading gene list from {gene_list_file}: {e}")
        return []

def load_gene_coordinates(genes, gtf_file):
    gene_coords = {}
    try:
        for line in open(gtf_file):
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            
            # Parse attributes
            attrs = dict(x.strip().split(' ', 1) for x in fields[8].rstrip(';').split('; '))
            gene_name = attrs.get('gene_name', '').strip('"')
            gene_id = attrs.get('gene_id', '').strip('"')
            
            if gene_name in genes:
                gene_coords[gene_name] = {
                    'chrom': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'strand': fields[6],
                    'gene_id': gene_id,
                    'gene_name': gene_name
                }
        
        logging.info(f"Loaded coordinates for {len(gene_coords)} genes")
        return gene_coords
    except Exception as e:
        logging.error(f"Error loading gene coordinates from {gtf_file}: {e}")
        return {}

def analyze_gene_chromatin_state(gene_coords, s2s_files, s3_files, cache):
    results = []
    
    for gene_name, coords in gene_coords.items():
        chrom = coords['chrom']
        if not chrom.startswith('chr'):
            chrom = f"chr{chrom}"
            
        # Define promoter region
        if coords['strand'] == '+':
            start = coords['start'] - WINDOW_SIZE // 2
            end = coords['start'] + WINDOW_SIZE // 2
        else:
            start = coords['end'] - WINDOW_SIZE // 2
            end = coords['end'] + WINDOW_SIZE // 2
            
        # Ensure start is not negative
        start = max(0, start)
        
        # Create cache key
        cache_key = {
            'gene': gene_name,
            'chrom': chrom,
            'start': start,
            'end': end,
            's2s_files': [str(f) for f in s2s_files],
            's3_files': [str(f) for f in s3_files]
        }
        
        # Try to get from cache
        cached_result = cache.get(cache_key)
        if cached_result is not None:
            results.append(cached_result)
            continue
            
        # Compute signals
        s2s_signal = compute_average_signal(s2s_files, chrom, start, end)
        s3_signal = compute_average_signal(s3_files, chrom, start, end)
        
        if s2s_signal is None or s3_signal is None:
            logging.warning(f"Could not compute signal for {gene_name}")
            continue
            
        # Compute ratio and determine state
        ratio = s3_signal / s2s_signal if s2s_signal > 0 else float('inf')
        if ratio > MIN_DIFF_THRESHOLD:
            state = 'condensed'
        elif 1/ratio > MIN_DIFF_THRESHOLD:
            state = 'decondensed'
        else:
            state = 'intermediate'
            
        result = {
            'gene': gene_name,
            'chrom': chrom,
            'start': start,
            'end': end,
            's2s_signal': s2s_signal,
            's3_signal': s3_signal,
            'ratio': ratio,
            'state': state
        }
        
        # Cache the result
        cache.set(cache_key, result)
        results.append(result)
    
    return pd.DataFrame(results)

def create_phase_comparison_plot(gfp_results, title, output_file):
    """Create a plot showing chromatin phase distribution for GFP-treated samples."""
    plt.figure(figsize=(10, 6))
    
    # Count states
    state_counts = gfp_results['state'].value_counts()
    total = len(gfp_results)
    
    # Create bar plot
    colors = {'condensed': 'red', 'decondensed': 'blue', 'intermediate': 'gray'}
    bars = plt.bar(state_counts.index, state_counts.values / total * 100)
    
    # Color bars
    for bar, state in zip(bars, state_counts.index):
        bar.set_color(colors[state])
    
    plt.title(title)
    plt.ylabel('Percentage of Genes')
    plt.ylim(0, 100)
    
    # Add percentage labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}%', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def create_ratio_heatmap(gfp_results, title, output_file):
    """Create a heatmap of S3/S2S ratios for GFP-treated samples."""
    plt.figure(figsize=(12, 8))
    
    # Prepare data for heatmap
    ratios = gfp_results.set_index('gene')['ratio']
    ratios = np.log2(ratios)  # Log transform for better visualization
    ratios = ratios.sort_values()
    
    # Create heatmap
    sns.heatmap(ratios.values.reshape(1, -1),
                cmap='RdBu_r',
                center=0,
                cbar_kws={'label': 'log2(S3/S2S ratio)'},
                yticklabels=False,
                xticklabels=False)
    
    plt.title(title)
    plt.xlabel('Genes (sorted by ratio)')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def analyze_enriched_genes(gene_list_file, gtf_file, bigwig_dir, output_dir, cache, chrom=None):
    """Analyze chromatin state of Mecp2-enriched genes in GFP-treated samples."""
    # Load enriched genes
    genes = load_enriched_genes(gene_list_file)
    if not genes:
        logging.error("No enriched genes found")
        return None
        
    # Load gene coordinates
    gene_coords = load_gene_coordinates(genes, gtf_file)
    if not gene_coords:
        logging.error("No gene coordinates found")
        return None
        
    # Filter genes by chromosome if specified
    if chrom:
        gene_coords = {name: coords for name, coords in gene_coords.items()
                      if coords['chrom'] == chrom or coords['chrom'] == f"chr{chrom}"}
        
    if not gene_coords:
        logging.warning(f"No genes found for chromosome {chrom}")
        return None
        
    # Get BigWig files for Neurons and NSCs
    neu_files = get_bigwig_files(bigwig_dir, "Neu", None, None)
    nsc_files = get_bigwig_files(bigwig_dir, "NSC", None, None)
    
    # Get BigWig files for GFP-treated samples (for phase analysis)
    gfp_s2s_files = get_bigwig_files(bigwig_dir, "GFP", "GFP", "S2S")
    gfp_s3_files = get_bigwig_files(bigwig_dir, "GFP", "GFP", "S3")
    
    if not gfp_s2s_files or not gfp_s3_files:
        logging.warning("Missing required GFP BigWig files for phase analysis")
    
    # Log the search path and files found
    logging.info(f"Searching for BigWig files in: {bigwig_dir}")
    logging.info(f"Found {len(neu_files)} Neuron files: {[f.name for f in neu_files]}")
    logging.info(f"Found {len(nsc_files)} NSC files: {[f.name for f in nsc_files]}")
    
    if not neu_files or not nsc_files:
        logging.error("Missing required BigWig files")
        return None
        
    # Calculate average signal for each cell type
    results = []
    for gene_name, coords in gene_coords.items():
        chrom = coords['chrom']
        if not chrom.startswith('chr'):
            chrom = f"chr{chrom}"
            
        # Define promoter region
        if coords['strand'] == '+':
            start = coords['start'] - WINDOW_SIZE // 2
            end = coords['start'] + WINDOW_SIZE // 2
        else:
            start = coords['end'] - WINDOW_SIZE // 2
            end = coords['end'] + WINDOW_SIZE // 2
            
        # Ensure start is not negative
        start = max(0, start)
        
        # Compute average signal for each cell type
        neu_signal = compute_average_signal(neu_files, chrom, start, end)
        nsc_signal = compute_average_signal(nsc_files, chrom, start, end)
        
        if neu_signal is not None and nsc_signal is not None:
            # Calculate fold change (Neurons vs NSCs)
            fold_change = neu_signal / nsc_signal if nsc_signal > 0 else float('inf')
            
            result = {
                'gene': gene_name,
                'chrom': chrom,
                'start': start,
                'end': end,
                'neu_signal': neu_signal,
                'nsc_signal': nsc_signal,
                'fold_change': fold_change,
                'enriched_in': 'Neurons' if fold_change > FOLD_CHANGE_THRESHOLD else ('NSCs' if fold_change < 1/FOLD_CHANGE_THRESHOLD else 'Neither')
            }
            results.append(result)
        
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    if results_df.empty:
        logging.warning("No results generated from analysis")
        return None
    
    # Analyze chromatin phase for enriched genes in GFP-treated samples
    if gfp_s2s_files and gfp_s3_files:
        logging.info("Analyzing chromatin phase in GFP-treated samples for enriched genes")
        # Filter gene coordinates to only include enriched genes
        enriched_genes = results_df[results_df['enriched_in'] == 'Neurons']['gene'].tolist()
        enriched_gene_coords = {name: coords for name, coords in gene_coords.items() if name in enriched_genes}
        
        # Analyze chromatin state for enriched genes in GFP-treated samples
        gfp_phase_results = analyze_gene_chromatin_state(enriched_gene_coords, gfp_s2s_files, gfp_s3_files, cache)
        
        if not gfp_phase_results.empty:
            # Save phase results to CSV
            phase_results_file = output_dir / "gfp_phase_analysis.csv"
            gfp_phase_results.to_csv(phase_results_file, index=False)
            logging.info(f"Saved GFP phase analysis results to {phase_results_file}")
            
            # Create phase comparison plot
            create_phase_comparison_plot(
                gfp_phase_results,
                "Chromatin Phase Distribution of Mecp2-Enriched Genes in GFP-Treated Samples",
                output_dir / "gfp_phase_distribution.png"
            )
            
            # Create ratio heatmap
            create_ratio_heatmap(
                gfp_phase_results,
                "S3/S2S Ratio of Mecp2-Enriched Genes in GFP-Treated Samples",
                output_dir / "gfp_ratio_heatmap.png"
            )
            
            # Merge phase information with results_df
            phase_info = gfp_phase_results[['gene', 'state', 'ratio']].rename(
                columns={'state': 'gfp_phase', 'ratio': 'gfp_s3_s2s_ratio'}
            )
            results_df = pd.merge(results_df, phase_info, on='gene', how='left')
        else:
            logging.warning("No phase analysis results generated for GFP-treated samples")
    
    # Create output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save results to CSV
    results_file = output_dir / "enriched_genes_signal_comparison.csv"
    results_df.to_csv(results_file, index=False)
    logging.info(f"Saved results to {results_file}")
    
    # Create summary plots
    plt.figure(figsize=(10, 6))
    enrichment_counts = results_df['enriched_in'].value_counts()
    plt.bar(enrichment_counts.index, enrichment_counts.values)
    plt.title('Distribution of Mecp2 Enrichment')
    plt.ylabel('Number of Genes')
    plt.tight_layout()
    plt.savefig(output_dir / "enrichment_distribution.png")
    plt.close()
    
    # Create fold change heatmap
    plt.figure(figsize=(12, 8))
    fold_changes = results_df.set_index('gene')['fold_change']
    fold_changes = np.log2(fold_changes)  # Log transform for better visualization
    fold_changes = fold_changes.sort_values()
    
    sns.heatmap(fold_changes.values.reshape(1, -1),
                cmap='RdBu_r',
                center=0,
                cbar_kws={'label': 'log2(Neurons/NSCs signal ratio)'},
                yticklabels=False,
                xticklabels=False)
    
    plt.title('Mecp2 Signal Fold Change (Neurons vs NSCs)')
    plt.xlabel('Genes (sorted by fold change)')
    plt.tight_layout()
    plt.savefig(output_dir / "fold_change_heatmap.png")
    plt.close()
    
    return results_df

def main():
    parser = argparse.ArgumentParser(description="Analyze Mecp2 signal enrichment and chromatin phase in GFP-treated samples")
    parser.add_argument('--chrom', type=str, help='Chromosome to analyze')
    parser.add_argument('--gtf_file', type=str, required=True, help='Path to GTF file')
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Setup directories
    results_dir = Path("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/results")
    results_dir, bigwig_dir, output_dir, cache_dir = setup_directories(results_dir, args.chrom)
    
    # Initialize cache
    cache = ResultsCache(cache_dir)
    
    # Define input files
    gene_list_file = Path("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/gene_lists/endo_enriched_gene_list_nsc_vs_neu.txt")
    
    # Run analysis
    results = analyze_enriched_genes(
        gene_list_file=gene_list_file,
        gtf_file=args.gtf_file,
        bigwig_dir=bigwig_dir,
        output_dir=output_dir,
        cache=cache,
        chrom=args.chrom
    )
    
    if results is None:
        logging.error("Analysis failed")
        sys.exit(1)
    
    # Log summary statistics
    enrichment_counts = results['enriched_in'].value_counts()
    logging.info("\nEnrichment Summary:")
    for category, count in enrichment_counts.items():
        logging.info(f"{category}: {count} genes ({count/len(results)*100:.1f}%)")
    
    # Log phase analysis summary if available
    if 'gfp_phase' in results.columns:
        phase_counts = results['gfp_phase'].value_counts()
        logging.info("\nGFP Phase Analysis Summary:")
        for phase, count in phase_counts.items():
            if pd.notna(phase):  # Only log non-NA values
                logging.info(f"{phase}: {count} genes ({count/len(results.dropna(subset=['gfp_phase']))*100:.1f}%)")
    
    logging.info("\nAnalysis completed successfully")

if __name__ == "__main__":
    main()
