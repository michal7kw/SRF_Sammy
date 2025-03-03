#!/usr/bin/env python3
"""
Script to merge chromosome-specific results from endogenous Mecp2 target gene analysis.
This script combines results from individual chromosome analyses into genome-wide results.
"""

import argparse
import logging
import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
import re

def setup_logging(debug=False):
    """Set up logging configuration."""
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def find_chromosome_results(results_dir):
    """Find all chromosome-specific result directories."""
    # Look for directories named like 'chr1', 'chr2', etc. in the endogenous_target_analysis subdirectory
    target_dir = os.path.join(results_dir, "endogenous_target_analysis")
    if not os.path.exists(target_dir):
        logging.warning(f"Directory {target_dir} does not exist")
        return []
    
    chrom_dirs = []
    for item in os.listdir(target_dir):
        item_path = os.path.join(target_dir, item)
        if os.path.isdir(item_path) and re.match(r'^chr[0-9XY]+$', item):
            chrom_dirs.append(item_path)
    
    logging.info(f"Found {len(chrom_dirs)} chromosome result directories")
    return chrom_dirs

def merge_gene_results(chrom_dirs, output_dir):
    """Merge gene-level results from all chromosomes."""
    logging.info("Merging gene-level results")
    
    # Initialize empty DataFrames for each cell type
    neu_results_all = pd.DataFrame()
    nsc_results_all = pd.DataFrame()
    
    # Collect results from each chromosome
    for chrom_dir in chrom_dirs:
        chrom = os.path.basename(chrom_dir)
        logging.info(f"Processing results from {chrom}")
        
        # Paths to gene results for this chromosome
        neu_file = os.path.join(chrom_dir, "task2_Neu_target_genes_phase.csv")
        nsc_file = os.path.join(chrom_dir, "task2_NSC_target_genes_phase.csv")
        
        # Check if files exist
        if os.path.exists(neu_file) and os.path.exists(nsc_file):
            # Load data
            neu_data = pd.read_csv(neu_file)
            nsc_data = pd.read_csv(nsc_file)
            
            # Add chromosome information if not already present
            if 'chromosome' not in neu_data.columns:
                neu_data['chromosome'] = chrom
            if 'chromosome' not in nsc_data.columns:
                nsc_data['chromosome'] = chrom
            
            # Append to combined results
            neu_results_all = pd.concat([neu_results_all, neu_data], ignore_index=True)
            nsc_results_all = pd.concat([nsc_results_all, nsc_data], ignore_index=True)
        else:
            logging.warning(f"Results files not found for {chrom}")
    
    # Save merged results
    neu_output = os.path.join(output_dir, "task2_Neu_target_genes_phase_all.csv")
    nsc_output = os.path.join(output_dir, "task2_NSC_target_genes_phase_all.csv")
    
    neu_results_all.to_csv(neu_output, index=False)
    nsc_results_all.to_csv(nsc_output, index=False)
    
    logging.info(f"Merged results saved to {neu_output} and {nsc_output}")
    
    return neu_results_all, nsc_results_all

def create_genome_wide_plots(neu_results, nsc_results, output_dir):
    """Create genome-wide summary plots from merged results."""
    logging.info("Creating genome-wide summary plots")
    
    # Create output directory for plots if it doesn't exist
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # 1. Phase comparison plot (S2S vs S3)
    create_phase_comparison_plot(neu_results, nsc_results, 
                                "Genome-wide Phase Comparison", 
                                os.path.join(plots_dir, "genome_wide_phase_comparison.png"))
    
    # 2. Ratio heatmap
    create_ratio_heatmap(neu_results, nsc_results,
                        "Genome-wide Signal Ratio Heatmap",
                        os.path.join(plots_dir, "genome_wide_ratio_heatmap.png"))
    
    # 3. Chromosome distribution plot
    create_chromosome_distribution_plot(neu_results, nsc_results,
                                       os.path.join(plots_dir, "chromosome_distribution.png"))
    
    logging.info(f"Genome-wide plots saved to {plots_dir}")

def create_phase_comparison_plot(neu_results, nsc_results, title, output_file):
    """Create a scatter plot comparing S2S and S3 signals across cell types."""
    plt.figure(figsize=(12, 10))
    
    # Neurons
    plt.subplot(2, 2, 1)
    plt.scatter(neu_results['s2s_signal'], neu_results['s3_signal'], alpha=0.5)
    plt.title('Neurons')
    plt.xlabel('S2S Signal')
    plt.ylabel('S3 Signal')
    plt.grid(True, alpha=0.3)
    
    # NSCs
    plt.subplot(2, 2, 2)
    plt.scatter(nsc_results['s2s_signal'], nsc_results['s3_signal'], alpha=0.5)
    plt.title('NSCs')
    plt.xlabel('S2S Signal')
    plt.ylabel('S3 Signal')
    plt.grid(True, alpha=0.3)
    
    # Combined with different colors
    plt.subplot(2, 1, 2)
    plt.scatter(neu_results['s2s_signal'], neu_results['s3_signal'], alpha=0.5, label='Neurons', color='blue')
    plt.scatter(nsc_results['s2s_signal'], nsc_results['s3_signal'], alpha=0.5, label='NSCs', color='red')
    plt.title('Combined Comparison')
    plt.xlabel('S2S Signal')
    plt.ylabel('S3 Signal')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def create_ratio_heatmap(neu_results, nsc_results, title, output_file):
    """Create a heatmap of signal ratios between conditions."""
    # Calculate ratios
    neu_ratio = neu_results['s3_signal'] / neu_results['s2s_signal']
    nsc_ratio = nsc_results['s3_signal'] / nsc_results['s2s_signal']
    
    # Replace infinities and NaNs
    neu_ratio = neu_ratio.replace([np.inf, -np.inf], np.nan).fillna(0)
    nsc_ratio = nsc_ratio.replace([np.inf, -np.inf], np.nan).fillna(0)
    
    # Create a DataFrame for the heatmap
    ratio_df = pd.DataFrame({
        'Neurons': neu_ratio,
        'NSCs': nsc_ratio
    })
    
    # Sort by neuron ratio for better visualization
    ratio_df = ratio_df.sort_values('Neurons')
    
    # Create heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(ratio_df.T, cmap='viridis', cbar_kws={'label': 'S3/S2S Ratio'})
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def create_chromosome_distribution_plot(neu_results, nsc_results, output_file):
    """Create a plot showing the distribution of genes across chromosomes."""
    # Count genes per chromosome
    neu_chrom_counts = neu_results['chromosome'].value_counts().sort_index()
    nsc_chrom_counts = nsc_results['chromosome'].value_counts().sort_index()
    
    # Combine into a DataFrame
    chrom_df = pd.DataFrame({
        'Neurons': neu_chrom_counts,
        'NSCs': nsc_chrom_counts
    })
    
    # Fill NaN with zeros
    chrom_df = chrom_df.fillna(0)
    
    # Create plot
    plt.figure(figsize=(14, 8))
    chrom_df.plot(kind='bar')
    plt.title('Gene Distribution Across Chromosomes')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Genes')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def compare_cell_types(neu_results, nsc_results, output_dir):
    """Compare results between cell types and generate comparison metrics."""
    logging.info("Comparing results between cell types")
    
    # Create a merged dataset with genes present in both cell types
    merged_df = pd.merge(
        neu_results, 
        nsc_results,
        on='gene_id',
        suffixes=('_neu', '_nsc')
    )
    
    # Calculate differential metrics
    merged_df['s2s_diff'] = merged_df['s2s_signal_neu'] - merged_df['s2s_signal_nsc']
    merged_df['s3_diff'] = merged_df['s3_signal_neu'] - merged_df['s3_signal_nsc']
    merged_df['ratio_diff'] = (merged_df['s3_signal_neu'] / merged_df['s2s_signal_neu']) - (merged_df['s3_signal_nsc'] / merged_df['s2s_signal_nsc'])
    
    # Save comparison results
    comparison_file = os.path.join(output_dir, "cell_type_comparison.csv")
    merged_df.to_csv(comparison_file, index=False)
    
    logging.info(f"Cell type comparison saved to {comparison_file}")
    
    # Create comparison plots
    plots_dir = os.path.join(output_dir, "plots")
    
    # Differential signal plot
    plt.figure(figsize=(12, 10))
    plt.subplot(2, 1, 1)
    plt.scatter(merged_df['s2s_signal_neu'], merged_df['s2s_signal_nsc'], alpha=0.5)
    plt.title('S2S Signal: Neurons vs NSCs')
    plt.xlabel('Neurons S2S Signal')
    plt.ylabel('NSCs S2S Signal')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 1, 2)
    plt.scatter(merged_df['s3_signal_neu'], merged_df['s3_signal_nsc'], alpha=0.5)
    plt.title('S3 Signal: Neurons vs NSCs')
    plt.xlabel('Neurons S3 Signal')
    plt.ylabel('NSCs S3 Signal')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "cell_type_signal_comparison.png"), dpi=300)
    plt.close()
    
    return merged_df

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Merge chromosome-specific results from endogenous Mecp2 target gene analysis')
    parser.add_argument('--results_dir', default='results', help='Directory containing chromosome-specific results')
    parser.add_argument('--output_dir', default=None, help='Directory for merged results (defaults to results_dir/endogenous_target_analysis/merged)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.debug)
    
    # Set output directory
    if args.output_dir is None:
        args.output_dir = os.path.join(args.results_dir, "endogenous_target_analysis", "merged")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    logging.info(f"Starting merge of chromosome-specific results")
    logging.info(f"Results directory: {args.results_dir}")
    logging.info(f"Output directory: {args.output_dir}")
    
    # Find chromosome-specific result directories
    chrom_dirs = find_chromosome_results(args.results_dir)
    
    if not chrom_dirs:
        logging.error("No chromosome result directories found. Exiting.")
        return
    
    # Merge gene-level results
    neu_results, nsc_results = merge_gene_results(chrom_dirs, args.output_dir)
    
    # Create genome-wide summary plots
    create_genome_wide_plots(neu_results, nsc_results, args.output_dir)
    
    # Compare cell types
    compare_cell_types(neu_results, nsc_results, args.output_dir)
    
    logging.info("Merge completed successfully")

if __name__ == "__main__":
    main() 