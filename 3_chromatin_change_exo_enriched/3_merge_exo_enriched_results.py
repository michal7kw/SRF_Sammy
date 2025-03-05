#!/usr/bin/env python3
"""
Script to merge chromosome-specific results from exogenous Mecp2 enriched genes analysis.
This script combines results from individual chromosome analyses into genome-wide results.
"""

import argparse  # For parsing command-line arguments
import logging  # For logging information and errors
import os  # For interacting with the operating system
import json  # For working with JSON data
import pandas as pd  # For data manipulation and analysis
import numpy as np  # For numerical operations
import matplotlib.pyplot as plt  # For creating plots
import seaborn as sns  # For enhanced statistical data visualization
from pathlib import Path  # For working with file paths
import re  # For regular expressions

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
    # Construct the path to the 'exo_enriched_analysis' subdirectory
    exo_dir = os.path.join(results_dir, "exo_enriched_analysis")
    
    # Check if the target directory exists
    if not os.path.exists(exo_dir):
        logging.warning(f"Directory {exo_dir} does not exist")
        return []
    
    chrom_dirs = []
    
    # Iterate through the items in the target directory
    for item in os.listdir(exo_dir):
        item_path = os.path.join(exo_dir, item)
        
        # Check if the item is a directory and its name matches the chromosome pattern
        if os.path.isdir(item_path) and re.match(r'^chr[0-9XY]+$', item):
            chrom_dirs.append(item_path)
    
    logging.info(f"Found {len(chrom_dirs)} chromosome result directories")
    return chrom_dirs

def merge_chromatin_state_results(chrom_dirs, output_dir):
    """Merge chromatin state results from all chromosomes into combined dataframes.

    Args:
        chrom_dirs (list): A list of paths to chromosome-specific result directories.
        output_dir (str): The directory where the merged results will be saved.

    Returns:
        dict: A dictionary containing DataFrames for each cell type and condition.
    """
    logging.info("Merging chromatin state results")
    
    # Initialize empty DataFrames for each cell type and condition
    results = {
        'Neu': {
            'Mecp2': pd.DataFrame(),
            'GFP': pd.DataFrame()
        },
        'NSC': {
            'Mecp2': pd.DataFrame(),
            'GFP': pd.DataFrame()
        }
    }
    
    # Iterate through each chromosome directory
    for chrom_dir in chrom_dirs:
        chrom = os.path.basename(chrom_dir)
        logging.info(f"Processing results from {chrom}")
        
        # Process each cell type and condition
        for cell_type in ['Neu', 'NSC']:
            for condition in ['Mecp2', 'GFP']:
                # Construct the path to the result file
                result_file = os.path.join(chrom_dir, f"task3_{cell_type}_{condition}_chromatin_state.csv")
                
                if os.path.exists(result_file):
                    # Load the data from the CSV file
                    data = pd.read_csv(result_file)
                    
                    # Add chromosome information if it's not already present
                    if 'chrom' not in data.columns:
                        data['chrom'] = chrom
                    
                    # Append the data to the combined results DataFrame
                    results[cell_type][condition] = pd.concat([results[cell_type][condition], data], ignore_index=True)
                else:
                    logging.warning(f"Results file not found: {result_file}")
    
    # Save the merged results to CSV files
    for cell_type in ['Neu', 'NSC']:
        for condition in ['Mecp2', 'GFP']:
            if not results[cell_type][condition].empty:
                output_file = os.path.join(output_dir, f"task3_{cell_type}_{condition}_chromatin_state_all.csv")
                results[cell_type][condition].to_csv(output_file, index=False)
                logging.info(f"Saved merged results to {output_file}")
    
    return results

def merge_comparison_results(chrom_dirs, output_dir):
    """Merge comparison results from all chromosomes.

    Args:
        chrom_dirs (list): A list of paths to chromosome-specific result directories.
        output_dir (str): The directory where the merged results will be saved.

    Returns:
        dict: A dictionary containing DataFrames for each comparison type.
    """
    logging.info("Merging comparison results")
    
    # Initialize empty DataFrames for each comparison type
    comparisons = {
        'state_changes': {
            'Neu': pd.DataFrame(),
            'NSC': pd.DataFrame()
        },
        'ratio_comparison': {
            'Neu': pd.DataFrame(),
            'NSC': pd.DataFrame()
        }
    }
    
    # Iterate through each chromosome directory
    for chrom_dir in chrom_dirs:
        chrom = os.path.basename(chrom_dir)
        logging.info(f"Processing comparison results from {chrom}")
        
        # Process each cell type and comparison type
        for cell_type in ['Neu', 'NSC']:
            # State changes CSV
            state_file = os.path.join(chrom_dir, f"task3_{cell_type}_state_changes.csv")
            if os.path.exists(state_file):
                data = pd.read_csv(state_file)
                data['chromosome'] = chrom
                comparisons['state_changes'][cell_type] = pd.concat(
                    [comparisons['state_changes'][cell_type], data], 
                    ignore_index=True
                )
            
            # Ratio comparison CSV
            ratio_file = os.path.join(chrom_dir, f"task3_{cell_type}_ratio_comparison.csv")
            if os.path.exists(ratio_file):
                data = pd.read_csv(ratio_file)
                data['chromosome'] = chrom
                comparisons['ratio_comparison'][cell_type] = pd.concat(
                    [comparisons['ratio_comparison'][cell_type], data], 
                    ignore_index=True
                )
    
    # Save the merged comparison results
    for comp_type in comparisons:
        for cell_type in ['Neu', 'NSC']:
            if not comparisons[comp_type][cell_type].empty:
                output_file = os.path.join(output_dir, f"task3_{cell_type}_{comp_type}_all.csv")
                comparisons[comp_type][cell_type].to_csv(output_file, index=False)
                logging.info(f"Saved merged {comp_type} for {cell_type} to {output_file}")
    
    return comparisons

def create_genome_wide_plots(results, comparisons, output_dir):
    """Create genome-wide summary plots from merged results.

    Args:
        results (dict): Dictionary containing merged chromatin state results.
        comparisons (dict): Dictionary containing merged comparison results.
        output_dir (str): The directory where the plots will be saved.
    """
    logging.info("Creating genome-wide summary plots")
    
    # Create output directory for plots if it doesn't exist
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # 1. Create chromatin state distribution plots
    create_state_distribution_plots(results, plots_dir)
    
    # 2. Create ratio comparison plots
    create_ratio_comparison_plots(comparisons, plots_dir)
    
    # 3. Create cell type comparison plots
    create_cell_type_comparison_plots(results, comparisons, plots_dir)
    
    logging.info(f"Genome-wide plots saved to {plots_dir}")

def create_state_distribution_plots(results, plots_dir):
    """Create plots showing the distribution of chromatin states.

    Args:
        results (dict): Dictionary containing merged chromatin state results.
        plots_dir (str): Directory to save the plots.
    """
    # Create state distribution plots for each cell type
    for cell_type in ['Neu', 'NSC']:
        if not results[cell_type]['Mecp2'].empty and not results[cell_type]['GFP'].empty:
            # Count states for each condition
            mecp2_states = results[cell_type]['Mecp2']['state'].value_counts()
            gfp_states = results[cell_type]['GFP']['state'].value_counts()
            
            # Combine into a DataFrame
            state_df = pd.DataFrame({
                'Mecp2': mecp2_states,
                'GFP': gfp_states
            }).fillna(0)
            
            # Create the plot
            plt.figure(figsize=(10, 6))
            state_df.plot(kind='bar')
            plt.title(f'Chromatin State Distribution in {cell_type}')
            plt.xlabel('Chromatin State')
            plt.ylabel('Number of Genes')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{cell_type}_state_distribution.png"), dpi=300)
            plt.close()
            
            # Save the data
            state_df.to_csv(os.path.join(plots_dir, f"{cell_type}_state_distribution.csv"))

def create_ratio_comparison_plots(comparisons, plots_dir):
    """Create plots comparing S3/S2S ratios between conditions.

    Args:
        comparisons (dict): Dictionary containing merged comparison results.
        plots_dir (str): Directory to save the plots.
    """
    # Create ratio comparison plots for each cell type
    for cell_type in ['Neu', 'NSC']:
        ratio_data = comparisons['ratio_comparison'][cell_type]
        if not ratio_data.empty:
            # Create scatter plot of log ratios
            plt.figure(figsize=(10, 8))
            plt.scatter(ratio_data['GFP_log_ratio'], ratio_data['Mecp2_log_ratio'], alpha=0.5)
            
            # Add diagonal line
            max_val = max(ratio_data['Mecp2_log_ratio'].max(), ratio_data['GFP_log_ratio'].max())
            min_val = min(ratio_data['Mecp2_log_ratio'].min(), ratio_data['GFP_log_ratio'].min())
            plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)
            
            plt.title(f'S3/S2S Ratio Comparison in {cell_type}')
            plt.xlabel('log2(S3/S2S) in GFP')
            plt.ylabel('log2(S3/S2S) in Mecp2')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{cell_type}_ratio_comparison.png"), dpi=300)
            plt.close()
            
            # Create histogram of ratio changes
            plt.figure(figsize=(10, 6))
            plt.hist(ratio_data['log_ratio_change'], bins=50, alpha=0.7)
            plt.title(f'Distribution of S3/S2S Ratio Changes in {cell_type}')
            plt.xlabel('log2(Mecp2/GFP Ratio)')
            plt.ylabel('Number of Genes')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{cell_type}_ratio_change_distribution.png"), dpi=300)
            plt.close()

def create_cell_type_comparison_plots(results, comparisons, plots_dir):
    """Create plots comparing results between cell types.

    Args:
        results (dict): Dictionary containing merged chromatin state results.
        comparisons (dict): Dictionary containing merged comparison results.
        plots_dir (str): Directory to save the plots.
    """
    # 1. Compare state changes between cell types
    neu_changes = comparisons['state_changes']['Neu']
    nsc_changes = comparisons['state_changes']['NSC']
    
    if not neu_changes.empty and not nsc_changes.empty:
        # Calculate percentage of genes that change state
        neu_changed = neu_changes[neu_changes['GFP_state'] != neu_changes['Mecp2_state']].shape[0]
        neu_total = neu_changes.shape[0]
        neu_pct_changed = (neu_changed / neu_total * 100) if neu_total > 0 else 0
        
        nsc_changed = nsc_changes[nsc_changes['GFP_state'] != nsc_changes['Mecp2_state']].shape[0]
        nsc_total = nsc_changes.shape[0]
        nsc_pct_changed = (nsc_changed / nsc_total * 100) if nsc_total > 0 else 0
        
        # Create bar chart
        plt.figure(figsize=(8, 6))
        bars = plt.bar(['Neurons', 'NSCs'], [neu_pct_changed, nsc_pct_changed])
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{height:.1f}%', ha='center', va='bottom')
        
        plt.title('Percentage of Genes with Chromatin State Changes')
        plt.ylabel('Percentage of Genes')
        plt.ylim(0, max(neu_pct_changed, nsc_pct_changed) * 1.2)
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "cell_type_state_change_comparison.png"), dpi=300)
        plt.close()
        
        # Create summary table
        summary = pd.DataFrame({
            'Cell Type': ['Neurons', 'NSCs'],
            'Total Genes': [neu_total, nsc_total],
            'Changed Genes': [neu_changed, nsc_changed],
            'Percentage Changed': [neu_pct_changed, nsc_pct_changed]
        })
        
        summary.to_csv(os.path.join(plots_dir, "cell_type_state_change_summary.csv"), index=False)
    
    # 2. Compare ratio changes between cell types
    neu_ratios = comparisons['ratio_comparison']['Neu']
    nsc_ratios = comparisons['ratio_comparison']['NSC']
    
    if not neu_ratios.empty and not nsc_ratios.empty:
        # Create boxplot of ratio changes
        plt.figure(figsize=(10, 6))
        data = [neu_ratios['log_ratio_change'], nsc_ratios['log_ratio_change']]
        plt.boxplot(data, labels=['Neurons', 'NSCs'])
        plt.title('Comparison of S3/S2S Ratio Changes Between Cell Types')
        plt.ylabel('log2(Mecp2/GFP Ratio)')
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "cell_type_ratio_change_comparison.png"), dpi=300)
        plt.close()
        
        # Calculate summary statistics
        neu_mean = neu_ratios['log_ratio_change'].mean()
        neu_median = neu_ratios['log_ratio_change'].median()
        nsc_mean = nsc_ratios['log_ratio_change'].mean()
        nsc_median = nsc_ratios['log_ratio_change'].median()
        
        # Create summary table
        summary = pd.DataFrame({
            'Cell Type': ['Neurons', 'NSCs'],
            'Mean Ratio Change': [neu_mean, nsc_mean],
            'Median Ratio Change': [neu_median, nsc_median]
        })
        
        summary.to_csv(os.path.join(plots_dir, "cell_type_ratio_change_summary.csv"), index=False)

def analyze_chromosome_distribution(results, output_dir):
    """Analyze the distribution of genes and chromatin states across chromosomes.

    Args:
        results (dict): Dictionary containing merged chromatin state results.
        output_dir (str): Directory to save the analysis results.
    """
    logging.info("Analyzing chromosome distribution")
    
    plots_dir = os.path.join(output_dir, "plots")
    
    # For each cell type
    for cell_type in ['Neu', 'NSC']:
        mecp2_data = results[cell_type]['Mecp2']
        gfp_data = results[cell_type]['GFP']
        
        if not mecp2_data.empty and not gfp_data.empty:
            # Count genes per chromosome
            mecp2_chrom_counts = mecp2_data['chrom'].value_counts().sort_index()
            gfp_chrom_counts = gfp_data['chrom'].value_counts().sort_index()
            
            # Combine into DataFrame
            chrom_df = pd.DataFrame({
                'Mecp2': mecp2_chrom_counts,
                'GFP': gfp_chrom_counts
            }).fillna(0)
            
            # Create chromosome distribution plot
            plt.figure(figsize=(14, 8))
            chrom_df.plot(kind='bar')
            plt.title(f'Gene Distribution Across Chromosomes in {cell_type}')
            plt.xlabel('Chromosome')
            plt.ylabel('Number of Genes')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{cell_type}_chromosome_distribution.png"), dpi=300)
            plt.close()
            
            # Save the data
            chrom_df.to_csv(os.path.join(plots_dir, f"{cell_type}_chromosome_distribution.csv"))
            
            # Analyze state distribution per chromosome
            # Create a pivot table of states by chromosome
            mecp2_pivot = pd.pivot_table(
                mecp2_data, 
                values='gene', 
                index='chrom',
                columns='state', 
                aggfunc='count',
                fill_value=0
            )
            
            gfp_pivot = pd.pivot_table(
                gfp_data, 
                values='gene', 
                index='chrom',
                columns='state', 
                aggfunc='count',
                fill_value=0
            )
            
            # Save the pivot tables
            mecp2_pivot.to_csv(os.path.join(plots_dir, f"{cell_type}_Mecp2_state_by_chromosome.csv"))
            gfp_pivot.to_csv(os.path.join(plots_dir, f"{cell_type}_GFP_state_by_chromosome.csv"))
            
            # Create heatmaps
            plt.figure(figsize=(12, 8))
            sns.heatmap(mecp2_pivot, cmap='viridis', annot=True, fmt='g')
            plt.title(f'Chromatin State Distribution by Chromosome in {cell_type} (Mecp2)')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{cell_type}_Mecp2_state_by_chromosome.png"), dpi=300)
            plt.close()
            
            plt.figure(figsize=(12, 8))
            sns.heatmap(gfp_pivot, cmap='viridis', annot=True, fmt='g')
            plt.title(f'Chromatin State Distribution by Chromosome in {cell_type} (GFP)')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{cell_type}_GFP_state_by_chromosome.png"), dpi=300)
            plt.close()

def main():
    """Main function to execute the script."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Merge chromosome-specific results from exogenous Mecp2 enriched genes analysis')
    parser.add_argument('--results_dir', default='results', help='Directory containing chromosome-specific results')
    parser.add_argument('--output_dir', default=None, help='Directory for merged results (defaults to results_dir/exo_enriched_analysis/merged)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.debug)
    
    # Set output directory
    if args.output_dir is None:
        args.output_dir = os.path.join(args.results_dir, "exo_enriched_analysis", "merged")
    
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
    
    # Merge comparison results
    comparisons = merge_comparison_results(chrom_dirs, args.output_dir)
    
    # Create genome-wide summary plots
    create_genome_wide_plots(results, comparisons, args.output_dir)
    
    # Analyze chromosome distribution
    analyze_chromosome_distribution(results, args.output_dir)
    
    logging.info("Merge completed successfully")

if __name__ == "__main__":
    main() 