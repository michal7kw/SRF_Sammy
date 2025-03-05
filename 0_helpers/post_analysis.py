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

def create_browser_view(df, cell_type, output_dir):
    """Create a genome browser-like visualization for a representative chromosome."""
    # Select a chromosome with significant changes
    if len(df) == 0:
        return
    
    # Find chromosome with most changes
    chrom_counts = df['chrom'].value_counts()
    if len(chrom_counts) == 0:
        return
    
    top_chrom = chrom_counts.idxmax()
    chrom_df = df[df['chrom'] == top_chrom].sort_values('start')
    
    if len(chrom_df) < 5:  # Need enough points to make a meaningful plot
        return
    
    # Create a browser-like view
    fig, axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True)
    
    # Get chromosome length
    chrom_length = chrom_df['end'].max()
    
    # Plot GFP S3-S2S differential
    axes[0].fill_between(chrom_df['start'], chrom_df['GFP_diff'], 
                        where=chrom_df['GFP_diff'] > 0, color='green', alpha=0.7)
    axes[0].fill_between(chrom_df['start'], chrom_df['GFP_diff'], 
                        where=chrom_df['GFP_diff'] < 0, color='yellow', alpha=0.7)
    axes[0].set_ylabel('GFP\nS3 vs S2S')
    axes[0].axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # Plot M2 S3-S2S differential
    axes[1].fill_between(chrom_df['start'], chrom_df['M2_diff'], 
                        where=chrom_df['M2_diff'] > 0, color='green', alpha=0.7)
    axes[1].fill_between(chrom_df['start'], chrom_df['M2_diff'], 
                        where=chrom_df['M2_diff'] < 0, color='yellow', alpha=0.7)
    axes[1].set_ylabel('M2\nS3 vs S2S')
    axes[1].axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # Plot differential between M2 and GFP
    axes[2].fill_between(chrom_df['start'], chrom_df['diff_signal'], 
                        where=chrom_df['diff_signal'] > 0, color='blue', alpha=0.7)
    axes[2].fill_between(chrom_df['start'], chrom_df['diff_signal'], 
                        where=chrom_df['diff_signal'] < 0, color='yellow', alpha=0.7)
    axes[2].set_ylabel('M2 vs GFP\nDifferential')
    axes[2].axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # Add domains as gray boxes
    significant_regions = chrom_df[abs(chrom_df['diff_signal']) >= MIN_DIFF_THRESHOLD]
    for _, region in significant_regions.iterrows():
        rect = plt.Rectangle((region['start'], axes[2].get_ylim()[0]), 
                            region['end'] - region['start'], 
                            0.1 * (axes[2].get_ylim()[1] - axes[2].get_ylim()[0]),
                            facecolor='gray', alpha=0.5)
        axes[2].add_patch(rect)
    
    # Set x-axis limits and labels
    axes[2].set_xlim(0, chrom_length)
    axes[2].set_xlabel(f'Chromosome {top_chrom} Position (bp)')
    
    # Add title
    plt.suptitle(f'Chromatin State Changes in {cell_type} - Chromosome {top_chrom}', fontsize=16)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.savefig(output_dir / f"{cell_type}_chr{top_chrom}_browser_view.pdf")
    plt.close()

def create_bed_files(df, cell_type, output_dir):
    """Create BED files for visualization in genome browsers."""
    try:
        # Check if required columns exist
        required_columns = ['chrom', 'start', 'end', 'change_type', 'is_significant']
        if not all(col in df.columns for col in required_columns):
            logging.warning(f"Cannot create BED files: missing required columns. Available columns: {df.columns.tolist()}")
            return
        
        # Filter for significant changes
        significant_df = df[df['is_significant']]
        
        # Create BED files for each change type
        s2s_to_s3 = significant_df[significant_df['change_type'] == 'S2S to S3']
        s3_to_s2s = significant_df[significant_df['change_type'] == 'S3 to S2S']
        
        # Save BED files
        s2s_to_s3_file = output_dir / f"{cell_type}_S2S_to_S3.bed"
        s3_to_s2s_file = output_dir / f"{cell_type}_S3_to_S2S.bed"
        
        # Write S2S to S3 changes
        with open(s2s_to_s3_file, 'w') as f:
            f.write('track name="{0}_S2S_to_S3" description="Regions changing from S2S to S3 in {0}" visibility=2 itemRgb="On"\n'.format(cell_type))
            for _, row in s2s_to_s3.iterrows():
                # BED format: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    row['chrom'],
                    row['start'],
                    row['end'],
                    f"{row['chrom']}:{row['start']}-{row['end']}",
                    int(min(1000, abs(row['condition_diff']) * 100)),  # Score based on difference
                    '.',  # No strand
                    row['start'],
                    row['end'],
                    '255,0,0'  # Red for S2S to S3
                ))
        
        # Write S3 to S2S changes
        with open(s3_to_s2s_file, 'w') as f:
            f.write('track name="{0}_S3_to_S2S" description="Regions changing from S3 to S2S in {0}" visibility=2 itemRgb="On"\n'.format(cell_type))
            for _, row in s3_to_s2s.iterrows():
                # BED format: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    row['chrom'],
                    row['start'],
                    row['end'],
                    f"{row['chrom']}:{row['start']}-{row['end']}",
                    int(min(1000, abs(row['condition_diff']) * 100)),  # Score based on difference
                    '.',  # No strand
                    row['start'],
                    row['end'],
                    '0,0,255'  # Blue for S3 to S2S
                ))
        
        logging.info(f"Created BED files for {cell_type}: {s2s_to_s3_file} and {s3_to_s2s_file}")
    except Exception as e:
        logging.error(f"Error creating BED files for {cell_type}: {str(e)}")
        # Don't raise the exception, just log it and continue

def load_gene_list(gene_list_file):
    """Load a list of genes from a text file."""
    try:
        with open(gene_list_file, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        logging.info(f"Loaded {len(genes)} genes from {gene_list_file}")
        return genes
    except FileNotFoundError:
        logging.error(f"Gene list file not found: {gene_list_file}")
        return []
    except Exception as e:
        logging.error(f"Error loading gene list: {str(e)}")
        return []

def load_gene_coordinates(genes, gtf_file):
    """Map genes to their genomic coordinates using a GTF file."""
    gene_coords = {}
    
    # Check if database exists, if not create it
    db_file = str(gtf_file) + '.db'
    try:
        if not os.path.exists(db_file):
            logging.info(f"Creating GTF database: {db_file}")
            db = gffutils.create_db(gtf_file, db_file, force=True, 
                                   disable_infer_genes=False, disable_infer_transcripts=False)
        else:
            db = gffutils.FeatureDB(db_file)
        
        # Map genes to coordinates
        for gene_id in genes:
            try:
                # Try direct lookup by ID
                gene = db[gene_id]
                if gene.seqid in STANDARD_CHROMOSOMES:  # Only include standard chromosomes
                    gene_coords[gene_id] = (gene.seqid, gene.start, gene.end, gene.strand)
            except Exception:
                try:
                    # Try lookup by gene name
                    matches = list(db.features_of_type('gene', {'gene_name': gene_id}))
                    if matches:
                        gene = matches[0]
                        if gene.seqid in STANDARD_CHROMOSOMES:  # Only include standard chromosomes
                            gene_coords[gene_id] = (gene.seqid, gene.start, gene.end, gene.strand)
                    else:
                        logging.warning(f"Gene {gene_id} not found in annotation")
                except Exception:
                    logging.warning(f"Gene {gene_id} not found in annotation")
        
        logging.info(f"Successfully mapped {len(gene_coords)} genes to coordinates")
        return gene_coords
    
    except Exception as e:
        logging.error(f"Error loading gene coordinates: {str(e)}")
        return {}

def analyze_gene_chromatin_state(gene_coords, bigwig_files_s2s, bigwig_files_s3, cache):
    """Analyze chromatin state for a set of genes."""
    # Create cache key
    cache_key = {
        'function': 'analyze_gene_chromatin_state',
        'gene_coords': sorted([(g, *coords) for g, coords in gene_coords.items()]),
        'bigwig_files_s2s': sorted([str(f) for f in bigwig_files_s2s]),
        'bigwig_files_s3': sorted([str(f) for f in bigwig_files_s3])
    }
    
    # Try to load from cache
    cached_result = cache.get(cache_key)
    if cached_result is not None:
        return cached_result
    
    results = []
    
    total_genes = len(gene_coords)
    processed = 0
    
    for gene_id, (chrom, start, end, strand) in gene_coords.items():
        processed += 1
        if processed % 100 == 0:
            logging.info(f"Processed {processed}/{total_genes} genes")
            
        # Compute average signals
        s2s_signal = compute_average_signal(bigwig_files_s2s, chrom, start, end)
        s3_signal = compute_average_signal(bigwig_files_s3, chrom, start, end)
        
        # Determine dominant phase
        if s2s_signal > 0:
            ratio = s3_signal / s2s_signal
        else:
            ratio = float('inf') if s3_signal > 0 else 0
            
        phase = "Heterochromatin (S3)" if ratio > 1 else "Euchromatin (S2S)"
        
        results.append({
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'end': end,
            'S2S_signal': s2s_signal,
            'S3_signal': s3_signal,
            'S3_to_S2S_ratio': ratio,
            'dominant_phase': phase
        })
    
    results_df = pd.DataFrame(results)
    
    # Cache the results before returning
    cache.set(cache_key, results_df)
    return results_df

def create_phase_comparison_plot(df1, df2, title, output_file):
    """Create a plot comparing chromatin phases between two conditions."""
    # Count phases in each dataset
    df1_counts = df1['dominant_phase'].value_counts().reset_index()
    df1_counts.columns = ['phase', 'count']
    df1_counts['condition'] = df1['cell_type'].iloc[0] if 'cell_type' in df1.columns else 'Condition 1'
    
    df2_counts = df2['dominant_phase'].value_counts().reset_index()
    df2_counts.columns = ['phase', 'count']
    df2_counts['condition'] = df2['cell_type'].iloc[0] if 'cell_type' in df2.columns else 'Condition 2'
    
    # Combine data
    combined = pd.concat([df1_counts, df2_counts])
    
    # Create plot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=combined, x='phase', y='count', hue='condition')
    plt.title(f"Chromatin Phase Distribution: {title}")
    plt.ylabel("Number of Genes")
    plt.xlabel("Dominant Chromatin Phase")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    
    # Log results
    logging.info(f"\nChromatin phase distribution for {title}:")
    for condition, df in [(df1_counts['condition'].iloc[0], df1), 
                         (df2_counts['condition'].iloc[0], df2)]:
        s3_count = sum(df['dominant_phase'] == 'Heterochromatin (S3)')
        s2s_count = sum(df['dominant_phase'] == 'Euchromatin (S2S)')
        total = len(df)
        logging.info(f"{condition}: S3 (Heterochromatin): {s3_count}/{total} ({s3_count/total*100:.1f}%), "
                    f"S2S (Euchromatin): {s2s_count}/{total} ({s2s_count/total*100:.1f}%)")




def task2_analyze_endogenous_target_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache):
    """Task 2: Check phase of target genes of endogenous Mecp2."""
    # Create cache key
    cache_key = {
        'function': 'task2_analyze_endogenous_target_genes',
        'gene_list_files': [str(f) for f in gene_list_files],
        'gtf_file': str(gtf_file)
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
            
            # Get GFP bigWig files for this cell type
            s2s_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S2S")
            s3_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S3")
            
            # Analyze gene phase
            logging.info(f"Analyzing chromatin state in {cell_type}...")
            cell_results = analyze_gene_chromatin_state(gene_coords, s2s_files, s3_files, cache)
            cell_results['cell_type'] = cell_type
            results[cell_type] = cell_results
            
            # Save individual results
            cell_results.to_csv(output_dir / f"task2_{cell_type}_target_genes_phase.csv", index=False)
    
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
    
    # Combine all results
    combined_results = pd.concat(list(results.values()), ignore_index=True) if results else pd.DataFrame()
    if not combined_results.empty:
        combined_results.to_csv(output_dir / "task2_target_genes_phase.csv", index=False)
    
    # Cache the results before returning
    cache.set(cache_key, combined_results)
    return combined_results


def analyze_state_changes_for_gene_set(gfp_df, m2_df, cell_type, output_dir):
    """Analyze chromatin state changes between GFP and M2 for a gene set."""
    # Merge dataframes
    gfp_subset = gfp_df[['gene_id', 'S3_to_S2S_ratio', 'dominant_phase']].copy()
    gfp_subset.columns = ['gene_id', 'GFP_ratio', 'GFP_phase']
    
    m2_subset = m2_df[['gene_id', 'S3_to_S2S_ratio', 'dominant_phase']].copy()
    m2_subset.columns = ['gene_id', 'M2_ratio', 'M2_phase']
    
    merged = pd.merge(gfp_subset, m2_subset, on='gene_id')
    
    # Determine if state changed
    merged['state_changed'] = merged['GFP_phase'] != merged['M2_phase']
    merged['change_direction'] = 'No Change'
    
    # For genes that changed state, determine direction
    changed_mask = merged['state_changed']
    s2s_to_s3_mask = (merged['GFP_phase'] == 'Euchromatin (S2S)') & (merged['M2_phase'] == 'Heterochromatin (S3)')
    s3_to_s2s_mask = (merged['GFP_phase'] == 'Heterochromatin (S3)') & (merged['M2_phase'] == 'Euchromatin (S2S)')
    
    merged.loc[s2s_to_s3_mask, 'change_direction'] = 'S2S→S3'
    merged.loc[s3_to_s2s_mask, 'change_direction'] = 'S3→S2S'
    
    # Calculate fold change in ratio
    merged['ratio_fold_change'] = merged['M2_ratio'] / merged['GFP_ratio']
    
    # Save results
    merged.to_csv(output_dir / f"task4_{cell_type}_state_changes.csv", index=False)
    
    # Count changes
    total_genes = len(merged)
    changed_genes = sum(merged['state_changed'])
    s2s_to_s3 = sum(s2s_to_s3_mask)
    s3_to_s2s = sum(s3_to_s2s_mask)
    
    logging.info(f"\nChromatin state changes in {cell_type} for exogenous Mecp2 enriched genes:")
    logging.info(f"Total genes analyzed: {total_genes}")
    logging.info(f"Genes with changed state: {changed_genes} ({changed_genes/total_genes*100:.1f}%)")
    logging.info(f"  S2S → S3 transitions: {s2s_to_s3} ({s2s_to_s3/total_genes*100:.1f}%)")
    logging.info(f"  S3 → S2S transitions: {s3_to_s2s} ({s3_to_s2s/total_genes*100:.1f}%)")
    
    # Create visualization of state changes
    plt.figure(figsize=(10, 8))
    
    # Create pie chart of state changes
    plt.subplot(2, 2, 1)
    labels = ['No Change', 'S2S → S3', 'S3 → S2S']
    sizes = [total_genes - changed_genes, s2s_to_s3, s3_to_s2s]
    colors = ['lightgray', 'lightgreen', 'lightblue']
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%')
    plt.title(f'Chromatin State Changes in {cell_type}')
    
    # Create scatter plot of ratios
    plt.subplot(2, 2, 2)
    plt.scatter(np.log2(merged['GFP_ratio']), np.log2(merged['M2_ratio']), 
               alpha=0.6, c=merged['state_changed'].map({True: 'red', False: 'blue'}))
    plt.axhline(y=0, color='black', linestyle='--', alpha=0.3)
    plt.axvline(x=0, color='black', linestyle='--', alpha=0.3)
    plt.xlabel('Log2(S3/S2S) in GFP')
    plt.ylabel('Log2(S3/S2S) in M2')
    plt.title('Chromatin State Ratio Comparison')
    
    # Create histogram of fold changes
    plt.subplot(2, 2, 3)
    plt.hist(np.log2(merged['ratio_fold_change'].clip(0.1, 10)), bins=20)
    plt.axvline(x=0, color='red', linestyle='--')
    plt.xlabel('Log2(M2/GFP Ratio)')
    plt.ylabel('Number of Genes')
    plt.title('Distribution of Ratio Changes')
    
    # Create bar plot of phase distribution
    plt.subplot(2, 2, 4)
    phase_counts = pd.DataFrame({
        'GFP': gfp_df['dominant_phase'].value_counts(),
        'M2': m2_df['dominant_phase'].value_counts()
    }).fillna(0)
    phase_counts.plot(kind='bar')
    plt.xlabel('Chromatin Phase')
    plt.ylabel('Number of Genes')
    plt.title('Phase Distribution Comparison')
    
    plt.tight_layout()
    plt.savefig(output_dir / f"task4_{cell_type}_state_change_analysis.pdf")
    plt.close()

def compare_cell_types(results_dict, output_dir):
    """Compare chromatin state changes between cell types."""
    cell_types = list(results_dict.keys())
    if len(cell_types) < 2:
        return
    
    # Create a Venn diagram of overlapping regions
    try:
        from matplotlib_venn import venn2
        
        # Convert regions to sets for comparison
        region_sets = {}
        for cell_type, df in results_dict.items():
            regions = set()
            for _, row in df.iterrows():
                regions.add(f"{row['chrom']}:{row['start']}-{row['end']}")
            region_sets[cell_type] = regions
        
        # Create Venn diagram
        plt.figure(figsize=(8, 8))
        venn2([region_sets[ct] for ct in cell_types], set_labels=cell_types)
        plt.title('Overlap of Chromatin State Changes Between Cell Types')
        plt.savefig(output_dir / "cell_type_comparison_venn.pdf")
        plt.close()
        
        # Find common regions and analyze direction consistency
        common_regions = region_sets[cell_types[0]].intersection(region_sets[cell_types[1]])
        if common_regions:
            logging.info(f"Found {len(common_regions)} regions that change in both cell types")
            
            # Create lookup dictionaries for each cell type
            lookups = {}
            for cell_type, df in results_dict.items():
                lookup = {}
                for _, row in df.iterrows():
                    region_key = f"{row['chrom']}:{row['start']}-{row['end']}"
                    lookup[region_key] = row['change_direction']
                lookups[cell_type] = lookup
            
            # Count consistent and inconsistent directions
            consistent = 0
            for region in common_regions:
                if lookups[cell_types[0]][region] == lookups[cell_types[1]][region]:
                    consistent += 1
            
            logging.info(f"Direction consistency in common regions: {consistent}/{len(common_regions)} ({consistent/len(common_regions)*100:.1f}%)")
            
    except ImportError:
        logging.warning("matplotlib_venn not installed, skipping Venn diagram")

def create_ratio_heatmap(df, output_file, title='Chromatin State Ratio Heatmap', column='condition_diff'):
    """Create a heatmap visualization of chromatin state ratio data.
    
    Args:
        df: DataFrame containing chromosome and position data with ratio values
        output_file: Path to save the output heatmap image
        title: Title for the heatmap plot
        column: Column name containing the ratio values to plot
    """
    try:
        if len(df) == 0:
            logging.warning(f"No data available for heatmap visualization")
            return
            
        # Create a pivot table for the heatmap
        if 'chrom' in df.columns and 'start' in df.columns and column in df.columns:
            # Select only chromosomes with standard names (chr1-19, chrX, chrY for mouse)
            std_chroms = [f"chr{c}" for c in STANDARD_CHROMOSOMES]
            df_filtered = df[df['chrom'].isin(std_chroms)]
            
            if len(df_filtered) == 0:
                logging.warning(f"No standard chromosomes found in data for heatmap visualization")
                return
                
            # Sort by chromosome and start position
            df_filtered = df_filtered.sort_values(['chrom', 'start'])
            
            # Create the figure
            plt.figure(figsize=(12, 8))
            
            # Create the heatmap using seaborn
            heatmap_data = df_filtered.pivot_table(
                index='chrom', 
                columns='start', 
                values=column,
                aggfunc='mean'
            )
            
            # Use a diverging colormap centered at 0
            sns.heatmap(
                heatmap_data,
                cmap='RdBu_r',
                center=0,
                robust=True
            )
            
            plt.title(title)
            plt.tight_layout()
            plt.savefig(output_file, dpi=300)
            plt.close()
            
            logging.info(f"Saved ratio heatmap to {output_file}")
        else:
            logging.warning(f"Required columns for heatmap not found in data")
    except Exception as e:
        logging.error(f"Error creating ratio heatmap: {e}")
        
    return

def main():
    # Add command-line argument parsing
    parser = argparse.ArgumentParser(description='Analyze chromatin state changes.')
    parser.add_argument('--results_dir', type=str, default='results', help='Directory for results')
    parser.add_argument('--chrom', type=str, help='Chromosome to analyze (e.g., "1", "X")')
    parser.add_argument('--skip_pvalues', action='store_true', help='Skip p-value calculations to speed up analysis')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging (DEBUG level)')
    args = parser.parse_args()
    
    # Configure logging - set level to INFO by default, DEBUG if verbose flag is used
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    start_time = datetime.now()
    logging.info("Starting chromatin state analysis")
    
    # Initialize with the provided results directory
    logging.info(f"Initializing analysis with results_dir: {args.results_dir}")
    
    # Set up directories
    results_dir, bigwig_dir, output_dir, cache_dir = setup_directories(args.results_dir, args.chrom)
    
    # Ensure chromosome has the correct format for mm10 genome
    chrom = args.chrom
    if chrom and not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
    
    # Only log directory info once
    logging.info(f"Using results directory: {results_dir}")
    logging.info(f"Using bigwig directory: {bigwig_dir}")
    logging.info(f"Output will be saved to: {output_dir}")
    logging.info(f"Cache directory: {cache_dir}")
    
    # Analysis parameters
    window_size = 10000  # 10kb windows
    
    if chrom:
        logging.info(f"Processing chromosome: {chrom}")
    
    # Initialize cache
    cache = ResultsCache(cache_dir)
    
    # Task 1: Identify chromatin state changes
    logging.info("\n==================================================")
    logging.info("Task 1: Analyzing chromatin state changes")
    logging.info("==================================================")
    
    # Analyze each cell type
    cell_types = ["Neu", "NSC"]
    results_dict = {}
    
    # Compare results between cell types
    if all(df is not None for df in results_dict.values()):
        compare_cell_types(results_dict, output_dir)
    
    
    # Task 2: Analyze endogenous target genes
    if not chrom:
        gene_list_files = {
            "Neu": "gene_lists/Neu_target_genes.txt",
            "NSC": "gene_lists/NSC_target_genes.txt"
        }
        gtf_file = "references/gencode.vM25.annotation.gtf"
        task2_analyze_endogenous_target_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache)
    
    # Task 3: Analyze endogenous enriched genes
    if not chrom:
        gene_list_file = "gene_lists/endogenous_enriched_genes.txt"
        task3_analyze_endogenous_enriched_genes(gene_list_file, gtf_file, bigwig_dir, output_dir, cache)
    
    # Task 4: Analyze exogenous enriched genes
    if not chrom:
        gene_list_files = {
            "Neu": "gene_lists/Neu_exogenous_enriched_genes.txt",
            "NSC": "gene_lists/NSC_exogenous_enriched_genes.txt"
        }
        task4_analyze_exogenous_enriched_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache)
    
    end_time = datetime.now()
    logging.info(f"\nAnalysis completed successfully in {end_time - start_time}")

if __name__ == "__main__":
    main()










def task3_analyze_endogenous_enriched_genes(gene_list_file, gtf_file, bigwig_dir, output_dir, cache):
    """Task 3: Check phase of genes where endogenous Mecp2 is enriched."""
    # Create cache key
    cache_key = {
        'function': 'task3_analyze_endogenous_enriched_genes',
        'gene_list_file': str(gene_list_file),
        'gtf_file': str(gtf_file)
    }
    
    # Try to load from cache
    cached_result = cache.get(cache_key)
    if cached_result is not None:
        return cached_result
    
    logging.info(f"\nPerforming Task 3: Analyzing phase of genes where endogenous Mecp2 is enriched")
    
    # Load gene list and coordinates
    genes = load_gene_list(gene_list_file)
    if not genes:
        logging.error("No genes loaded for Task 3, skipping analysis")
        return None
        
    gene_coords = load_gene_coordinates(genes, gtf_file)
    if not gene_coords:
        logging.error("No gene coordinates mapped for Task 3, skipping analysis")
        return None
    
    # Get GFP bigWig files for both cell types
    neu_s2s_files = get_bigwig_files(bigwig_dir, "Neu", "GFP", "S2S")
    neu_s3_files = get_bigwig_files(bigwig_dir, "Neu", "GFP", "S3")
    nsc_s2s_files = get_bigwig_files(bigwig_dir, "NSC", "GFP", "S2S")
    nsc_s3_files = get_bigwig_files(bigwig_dir, "NSC", "GFP", "S3")
    
    # Analyze gene phase in each cell type
    logging.info("Analyzing chromatin state in Neurons...")
    neu_results = analyze_gene_chromatin_state(gene_coords, neu_s2s_files, neu_s3_files, cache)
    logging.info("Analyzing chromatin state in NSCs...")
    nsc_results = analyze_gene_chromatin_state(gene_coords, nsc_s2s_files, nsc_s3_files, cache)
    
    # Add cell type info
    neu_results['cell_type'] = 'Neurons'
    nsc_results['cell_type'] = 'NSC'
    
    # Combine results
    combined_results = pd.concat([neu_results, nsc_results], ignore_index=True)
    combined_results.to_csv(output_dir / "task3_endogenous_enriched_genes_phase.csv", index=False)
    
    # Visualize phase distribution
    create_phase_comparison_plot(neu_results, nsc_results, "Endogenous Mecp2 Enriched Genes", 
                               output_dir / "task3_phase_comparison.pdf")
    
    # Create a heatmap of S3/S2S ratios
    create_ratio_heatmap(neu_results, nsc_results, "Endogenous Mecp2 Enriched Genes",
                       output_dir / "task3_ratio_heatmap.pdf")
    
    # Cache the results before returning
    cache.set(cache_key, combined_results)
    return combined_results

def task4_analyze_exogenous_enriched_genes(gene_list_files, gtf_file, bigwig_dir, output_dir, cache):
    """Task 4: Check phase of genes where exogenous Mecp2 is enriched."""
    # Create cache key
    cache_key = {
        'function': 'task4_analyze_exogenous_enriched_genes',
        'gene_list_files': [str(f) for f in gene_list_files],
        'gtf_file': str(gtf_file)
    }
    
    # Try to load from cache
    cached_result = cache.get(cache_key)
    if cached_result is not None:
        return cached_result
    
    logging.info(f"\nPerforming Task 4: Analyzing phase of genes where exogenous Mecp2 is enriched")
    
    # Process each cell type with its corresponding gene list
    cell_types = ['Neu', 'NSC']
    all_results = {}
    
    for i, cell_type in enumerate(cell_types):
        if i < len(gene_list_files):
            gene_list_file = gene_list_files[i]
            logging.info(f"Processing {cell_type} with gene list: {gene_list_file}")
            
            # Load gene list and coordinates
            genes = load_gene_list(gene_list_file)
            if not genes:
                logging.error(f"No genes loaded for {cell_type} in Task 4, skipping")
                continue
                
            gene_coords = load_gene_coordinates(genes, gtf_file)
            if not gene_coords:
                logging.error(f"No gene coordinates mapped for {cell_type} in Task 4, skipping")
                continue
            
            # Get bigWig files for both conditions for this cell type
            gfp_s2s_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S2S")
            gfp_s3_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S3")
            m2_s2s_files = get_bigwig_files(bigwig_dir, cell_type, "M2", "S2S")
            m2_s3_files = get_bigwig_files(bigwig_dir, cell_type, "M2", "S3")
            
            # Analyze chromatin state for each condition
            logging.info(f"Analyzing chromatin state in {cell_type} with GFP...")
            gfp_results = analyze_gene_chromatin_state(gene_coords, gfp_s2s_files, gfp_s3_files, cache)
            gfp_results['condition'] = 'GFP'
            gfp_results['cell_type'] = cell_type
            
            logging.info(f"Analyzing chromatin state in {cell_type} with M2...")
            m2_results = analyze_gene_chromatin_state(gene_coords, m2_s2s_files, m2_s3_files, cache)
            m2_results['condition'] = 'M2'
            m2_results['cell_type'] = cell_type
            
            # Save individual results
            gfp_results.to_csv(output_dir / f"task4_{cell_type}_GFP_chromatin_state.csv", index=False)
            m2_results.to_csv(output_dir / f"task4_{cell_type}_M2_chromatin_state.csv", index=False)
            
            # Create comparison plot
            create_phase_comparison_plot(
                gfp_results, m2_results, 
                f"Exogenous Mecp2 Enriched Genes - {cell_type}", 
                output_dir / f"task4_{cell_type}_phase_comparison.pdf"
            )
            
            # Create ratio heatmap
            create_ratio_heatmap(
                gfp_results, m2_results,
                f"Exogenous Mecp2 Enriched Genes - {cell_type}",
                output_dir / f"task4_{cell_type}_ratio_heatmap.pdf"
            )
            
            # Analyze state changes
            analyze_state_changes_for_gene_set(gfp_results, m2_results, cell_type, output_dir)
            
            # Store results
            all_results[f"{cell_type}_GFP"] = gfp_results
            all_results[f"{cell_type}_M2"] = m2_results
    
    # Cache the results before returning
    cache.set(cache_key, all_results)
    return all_results
