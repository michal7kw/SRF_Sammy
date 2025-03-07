#!/usr/bin/env python3
"""
Script to merge chromosome-specific results from endogenous Mecp2 enriched genes analysis.
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
import re

def setup_logging(debug=False):
    """Set up logging configuration.

    Args:
        debug (bool, optional): If True, set logging level to DEBUG; otherwise, set to INFO. Defaults to False.
    """
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def find_chromosome_results(results_dir):
    """Find all chromosome-specific result directories within the main results directory.

    Args:
        results_dir (str): The main directory containing the results.

    Returns:
        list: A list of paths to chromosome-specific result directories.
    """
    # Construct the path to the 'endogenous_enriched_analysis' subdirectory
    endo_dir = os.path.join(results_dir, "endogenous_enriched_analysis")
    
    # Check if the target directory exists
    if not os.path.exists(endo_dir):
        logging.warning(f"Directory {endo_dir} does not exist")
        return []
    
    chrom_dirs = []
    
    # Iterate through the items in the target directory
    for item in os.listdir(endo_dir):
        item_path = os.path.join(endo_dir, item)
        
        # Check if the item is a directory and its name matches the chromosome pattern
        if os.path.isdir(item_path) and re.match(r'^chr[0-9XY]+$', item):
            chrom_dirs.append(item_path)
    
    logging.info(f"Found {len(chrom_dirs)} chromosome result directories")
    return chrom_dirs

def merge_chromatin_state_results(chrom_dirs, output_dir):
    """Merge chromatin state results from all chromosomes into a combined dataframe.

    Args:
        chrom_dirs (list): A list of paths to chromosome-specific result directories.
        output_dir (str): The directory where the merged results will be saved.

    Returns:
        pd.DataFrame: A DataFrame containing the merged results.
    """
    logging.info("Merging chromatin state results")
    
    # Initialize an empty DataFrame for the combined results
    all_results = pd.DataFrame()
    
    # Iterate through each chromosome directory
    for chrom_dir in chrom_dirs:
        chrom = os.path.basename(chrom_dir)
        logging.info(f"Processing results from {chrom}")
        
        # Construct the path to the result file
        result_file = os.path.join(chrom_dir, "enriched_genes_chromatin_state.csv")
        
        if os.path.exists(result_file):
            # Load the data from the CSV file
            data = pd.read_csv(result_file)
            
            # Append the data to the combined results DataFrame
            all_results = pd.concat([all_results, data], ignore_index=True)
        else:
            logging.warning(f"Results file not found: {result_file}")
    
    # Save the merged results to a CSV file
    if not all_results.empty:
        output_file = os.path.join(output_dir, "enriched_genes_chromatin_state_all.csv")
        all_results.to_csv(output_file, index=False)
        logging.info(f"Saved merged results to {output_file}")
    else:
        logging.warning("No results to merge")
    
    return all_results

def create_genome_wide_plots(results, output_dir):
    """Create genome-wide summary plots from merged results.

    Args:
        results (pd.DataFrame): DataFrame containing merged chromatin state results.
        output_dir (str): The directory where the plots will be saved.
    """
    if results.empty:
        logging.warning("No results available for creating plots")
        return
    
    logging.info("Creating genome-wide summary plots")
    
    # Create output directory for plots if it doesn't exist
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # 1. Create phase distribution plot
    create_phase_distribution_plot(results, plots_dir)
    
    # 2. Create ratio distribution plot
    create_ratio_distribution_plot(results, plots_dir)
    
    # 3. Create chromosome distribution plot
    create_chromosome_distribution_plot(results, plots_dir)
    
    logging.info(f"Genome-wide plots saved to {plots_dir}")

def create_phase_distribution_plot(results, plots_dir):
    """Create a plot showing the distribution of chromatin phases.

    Args:
        results (pd.DataFrame): DataFrame containing chromatin state results.
        plots_dir (str): Directory to save the plot.
    """
    plt.figure(figsize=(10, 6))
    
    # Count states
    state_counts = results['state'].value_counts()
    total = len(results)
    
    # Create bar plot
    colors = {'condensed': 'red', 'decondensed': 'blue', 'intermediate': 'gray'}
    bars = plt.bar(state_counts.index, state_counts.values / total * 100)
    
    # Color bars
    for bar, state in zip(bars, state_counts.index):
        if state in colors:
            bar.set_color(colors[state])
    
    plt.title('Genome-wide Chromatin Phase Distribution in Mecp2-Enriched Genes')
    plt.ylabel('Percentage of Genes')
    plt.ylim(0, 100)
    
    # Add percentage labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}%', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "genome_wide_phase_distribution.png"), dpi=300)
    plt.close()
    
    # Save the data
    state_df = pd.DataFrame({
        'State': state_counts.index,
        'Count': state_counts.values,
        'Percentage': state_counts.values / total * 100
    })
    state_df.to_csv(os.path.join(plots_dir, "genome_wide_phase_distribution.csv"), index=False)

def create_ratio_distribution_plot(results, plots_dir):
    """Create plots showing the distribution of S3/S2S ratios.

    Args:
        results (pd.DataFrame): DataFrame containing chromatin state results.
        plots_dir (str): Directory to save the plots.
    """
    # 1. Create histogram of ratios
    plt.figure(figsize=(12, 6))
    
    # Log transform ratios for better visualization
    log_ratios = np.log2(results['ratio'])
    
    # Create histogram
    plt.hist(log_ratios, bins=50, alpha=0.7)
    plt.title('Distribution of S3/S2S Ratios in Mecp2-Enriched Genes')
    plt.xlabel('log2(S3/S2S ratio)')
    plt.ylabel('Number of Genes')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "ratio_distribution_histogram.png"), dpi=300)
    plt.close()
    
    # 2. Create boxplot of ratios by state
    plt.figure(figsize=(10, 6))
    
    # Group by state and create boxplot
    sns.boxplot(x='state', y=log_ratios, data=results)
    plt.title('S3/S2S Ratio by Chromatin State')
    plt.xlabel('Chromatin State')
    plt.ylabel('log2(S3/S2S ratio)')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "ratio_by_state_boxplot.png"), dpi=300)
    plt.close()
    
    # 3. Create heatmap of ratios
    plt.figure(figsize=(14, 4))
    
    # Sort ratios
    sorted_ratios = log_ratios.sort_values()
    
    # Create heatmap
    sns.heatmap(sorted_ratios.values.reshape(1, -1),
                cmap='RdBu_r',
                center=0,
                cbar_kws={'label': 'log2(S3/S2S ratio)'},
                yticklabels=False,
                xticklabels=False)
    
    plt.title('Genome-wide S3/S2S Ratio Distribution in Mecp2-Enriched Genes')
    plt.xlabel('Genes (sorted by ratio)')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "genome_wide_ratio_heatmap.png"), dpi=300)
    plt.close()

def create_chromosome_distribution_plot(results, plots_dir):
    """Create plots showing the distribution of genes and chromatin states across chromosomes.

    Args:
        results (pd.DataFrame): DataFrame containing chromatin state results.
        plots_dir (str): Directory to save the plots.
    """
    # 1. Count genes per chromosome
    chrom_counts = results['chrom'].value_counts().sort_index()
    
    # Create bar plot
    plt.figure(figsize=(14, 8))
    chrom_counts.plot(kind='bar')
    plt.title('Gene Distribution Across Chromosomes')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Genes')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "chromosome_distribution.png"), dpi=300)
    plt.close()
    
    # Save the data
    chrom_df = pd.DataFrame({
        'Chromosome': chrom_counts.index,
        'Gene Count': chrom_counts.values
    })
    chrom_df.to_csv(os.path.join(plots_dir, "chromosome_distribution.csv"), index=False)
    
    # 2. Create a heatmap of states by chromosome
    # Create a pivot table of states by chromosome
    state_pivot = pd.pivot_table(
        results, 
        values='gene', 
        index='chrom',
        columns='state', 
        aggfunc='count',
        fill_value=0
    )
    
    # Create heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(state_pivot, cmap='viridis', annot=True, fmt='g')
    plt.title('Chromatin State Distribution by Chromosome')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "state_by_chromosome.png"), dpi=300)
    plt.close()
    
    # Save the pivot table
    state_pivot.to_csv(os.path.join(plots_dir, "state_by_chromosome.csv"))
    
    # 3. Create a stacked bar chart of states by chromosome
    plt.figure(figsize=(14, 8))
    state_pivot.plot(kind='bar', stacked=True)
    plt.title('Chromatin State Distribution by Chromosome')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Genes')
    plt.grid(True, alpha=0.3)
    plt.legend(title='Chromatin State')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "state_by_chromosome_stacked.png"), dpi=300)
    plt.close()

def analyze_signal_distribution(results, output_dir):
    """Analyze the distribution of S2S and S3 signals.

    Args:
        results (pd.DataFrame): DataFrame containing chromatin state results.
        output_dir (str): Directory to save the analysis results.
    """
    if results.empty:
        logging.warning("No results available for signal distribution analysis")
        return
    
    logging.info("Analyzing signal distribution")
    
    plots_dir = os.path.join(output_dir, "plots")
    
    # 1. Create scatter plot of S2S vs S3 signals
    plt.figure(figsize=(10, 8))
    plt.scatter(results['s2s_signal'], results['s3_signal'], alpha=0.5)
    plt.title('S2S vs S3 Signal in Mecp2-Enriched Genes')
    plt.xlabel('S2S Signal')
    plt.ylabel('S3 Signal')
    plt.grid(True, alpha=0.3)
    
    # Add diagonal line
    max_val = max(results['s2s_signal'].max(), results['s3_signal'].max())
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "s2s_vs_s3_scatter.png"), dpi=300)
    plt.close()
    
    # 2. Create histograms of signals
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # S2S signal histogram
    ax1.hist(results['s2s_signal'], bins=50, alpha=0.7)
    ax1.set_title('Distribution of S2S Signal')
    ax1.set_xlabel('S2S Signal')
    ax1.set_ylabel('Number of Genes')
    ax1.grid(True, alpha=0.3)
    
    # S3 signal histogram
    ax2.hist(results['s3_signal'], bins=50, alpha=0.7)
    ax2.set_title('Distribution of S3 Signal')
    ax2.set_xlabel('S3 Signal')
    ax2.set_ylabel('Number of Genes')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "signal_histograms.png"), dpi=300)
    plt.close()
    
    # 3. Create boxplots of signals by state
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # S2S signal boxplot
    sns.boxplot(x='state', y='s2s_signal', data=results, ax=ax1)
    ax1.set_title('S2S Signal by Chromatin State')
    ax1.set_xlabel('Chromatin State')
    ax1.set_ylabel('S2S Signal')
    ax1.grid(True, alpha=0.3)
    
    # S3 signal boxplot
    sns.boxplot(x='state', y='s3_signal', data=results, ax=ax2)
    ax2.set_title('S3 Signal by Chromatin State')
    ax2.set_xlabel('Chromatin State')
    ax2.set_ylabel('S3 Signal')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "signal_by_state_boxplots.png"), dpi=300)
    plt.close()
    
    # 4. Calculate and save summary statistics
    summary = pd.DataFrame({
        'Statistic': ['Mean', 'Median', 'Min', 'Max', 'Std Dev'],
        'S2S Signal': [
            results['s2s_signal'].mean(),
            results['s2s_signal'].median(),
            results['s2s_signal'].min(),
            results['s2s_signal'].max(),
            results['s2s_signal'].std()
        ],
        'S3 Signal': [
            results['s3_signal'].mean(),
            results['s3_signal'].median(),
            results['s3_signal'].min(),
            results['s3_signal'].max(),
            results['s3_signal'].std()
        ],
        'Ratio': [
            results['ratio'].mean(),
            results['ratio'].median(),
            results['ratio'].min(),
            results['ratio'].max(),
            results['ratio'].std()
        ]
    })
    
    summary.to_csv(os.path.join(plots_dir, "signal_summary_statistics.csv"), index=False)

def main():
    """Main function to execute the script."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Merge chromosome-specific results from endogenous Mecp2 enriched genes analysis')
    parser.add_argument('--results_dir', default='results', help='Directory containing chromosome-specific results')
    parser.add_argument('--output_dir', default=None, help='Directory for merged results (defaults to results_dir/endogenous_enriched_analysis/merged)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.debug)
    
    # Set output directory
    if args.output_dir is None:
        args.output_dir = os.path.join(args.results_dir, "endogenous_enriched_analysis", "merged")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    logging.info(f"Starting merge of chromosome-specific results")
    logging.info(f"Results directory: {args.results_dir}")
    logging.info(f"Output directory: {args.output_dir}")
    
    # Find chromosome-specific result directories
    chrom_dirs = find_chromosome_results(args.results_dir)
    
    # Check if any chromosome directories were found
    if not chrom_dirs:
        logging.error("No chromosome result directories found. Exiting.")
        return
    
    # Merge chromatin state results
    results = merge_chromatin_state_results(chrom_dirs, args.output_dir)
    
    # Create genome-wide summary plots
    create_genome_wide_plots(results, args.output_dir)
    
    # Analyze signal distribution
    analyze_signal_distribution(results, args.output_dir)
    
    logging.info("Merge completed successfully")

if __name__ == "__main__":
    main() 