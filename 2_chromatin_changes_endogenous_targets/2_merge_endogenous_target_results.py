#!/usr/bin/env python3
"""
Script to merge chromosome-specific results from endogenous Mecp2 target gene analysis.
This script combines results from individual chromosome analyses into genome-wide results.
"""

import argparse  # For parsing command-line arguments
import logging  # For logging information and errors
import os  # For interacting with the operating system (e.g., file paths)
import json  # For working with JSON data (not currently used, but might be useful in the future)
import pandas as pd  # For data manipulation and analysis using DataFrames
import numpy as np  # For numerical opediffns
import matplotlib.pyplot as plt  # For creating plots
import seaborn as sns  # For enhanced statistical data visualization
from pathlib import Path  # For working with file paths in an object-oriented way
import glob  # For finding files matching a pattern (not currently used)
import re  # For regular expressions

def setup_logging(debug=False):
    """Set up logging configudiffn.

    Args:
        debug (bool, optional): If True, set logging level to DEBUG; otherwise, set to INFO. Defaults to False.
    """
    log_level = logging.DEBUG if debug else logging.INFO  # Determine logging level based on debug flag
    logging.basicConfig(  # Configure basic logging settings
        level=log_level,  # Set the logging level
        format='%(asctime)s - %(levelname)s - %(message)s',  # Define the log message format
        datefmt='%Y-%m-%d %H:%M:%S'  # Define the date format in log messages
    )

def find_chromosome_results(results_dir):
    """Find all chromosome-specific result directories within the main results directory.

    Args:
        results_dir (str): The main directory containing the results.

    Returns:
        list: A list of paths to chromosome-specific result directories.
    """
    # Construct the path to the 'endogenous_target_analysis' subdirectory
    target_dir = os.path.join(results_dir, "endogenous_target_analysis")
    
    # Check if the target directory exists
    if not os.path.exists(target_dir):
        logging.warning(f"Directory {target_dir} does not exist")  # Log a warning if the directory is not found
        return []  # Return an empty list if the target directory doesn't exist
    
    chrom_dirs = []  # Initialize an empty list to store chromosome directory paths
    
    # Iterate through the items in the target directory
    for item in os.listdir(target_dir):
        item_path = os.path.join(target_dir, item)  # Construct the full path to the item
        
        # Check if the item is a directory and its name matches the chromosome pattern (e.g., 'chr1', 'chrX')
        if os.path.isdir(item_path) and re.match(r'^chr[0-9XY]+$', item):
            chrom_dirs.append(item_path)  # Add the directory path to the list
    
    logging.info(f"Found {len(chrom_dirs)} chromosome result directories")  # Log the number of chromosome directories found
    return chrom_dirs  # Return the list of chromosome directory paths

def merge_gene_results(chrom_dirs, output_dir):
    """Merge gene-level results from all chromosomes into combined dataframes and save to CSV files.

    Args:
        chrom_dirs (list): A list of paths to chromosome-specific result directories.
        output_dir (str): The directory where the merged results will be saved.

    Returns:
        tuple: A tuple containing two pandas DataFrames: neu_results_all and nsc_results_all.
    """
    logging.info("Merging gene-level results")
    
    # Initialize empty DataFrames for each cell type to store the combined results
    neu_results_all = pd.DataFrame()
    nsc_results_all = pd.DataFrame()
    
    # Iterate through each chromosome directory
    for chrom_dir in chrom_dirs:
        chrom = os.path.basename(chrom_dir)  # Extract the chromosome name from the directory path
        logging.info(f"Processing results from {chrom}")
        
        # Construct the paths to the gene result files for this chromosome
        neu_file = os.path.join(chrom_dir, "task2_Neu_target_genes_phase.csv")
        nsc_file = os.path.join(chrom_dir, "task2_NSC_target_genes_phase.csv")
        
        # Check if both result files exist for this chromosome
        if os.path.exists(neu_file) and os.path.exists(nsc_file):
            # Load the data from the CSV files into pandas DataFrames
            neu_data = pd.read_csv(neu_file)
            nsc_data = pd.read_csv(nsc_file)
            
            # Add chromosome information to the DataFrames if it's not already present
            if 'chromosome' not in neu_data.columns:
                neu_data['chromosome'] = chrom
            if 'chromosome' not in nsc_data.columns:
                nsc_data['chromosome'] = chrom
            
            # Append the data to the combined results DataFrames
            neu_results_all = pd.concat([neu_results_all, neu_data], ignore_index=True)
            nsc_results_all = pd.concat([nsc_results_all, nsc_data], ignore_index=True)
        else:
            logging.warning(f"Results files not found for {chrom}")  # Log a warning if the result files are missing
    
    # Construct the output file paths for the merged results
    neu_output = os.path.join(output_dir, "task2_Neu_target_genes_phase_all.csv")
    nsc_output = os.path.join(output_dir, "task2_NSC_target_genes_phase_all.csv")
    
    # Save the merged results DataFrames to CSV files
    neu_results_all.to_csv(neu_output, index=False)
    nsc_results_all.to_csv(nsc_output, index=False)
    
    logging.info(f"Merged results saved to {neu_output} and {nsc_output}")  # Log the successful saving of merged results
    
    return neu_results_all, nsc_results_all  # Return the merged DataFrames

def create_genome_wide_plots(neu_results, nsc_results, output_dir):
    """Create genome-wide summary plots from merged results.

    Args:
        neu_results (pd.DataFrame): DataFrame containing merged results for neurons.
        nsc_results (pd.DataFrame): DataFrame containing merged results for NSCs.
        output_dir (str): The directory where the plots will be saved.
    """
    logging.info("Creating genome-wide summary plots")
    
    # Create output directory for plots if it doesn't exist
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # 1. Phase comparison plot (S2S vs S3)
    create_phase_comparison_plot(neu_results, nsc_results, 
                                "Genome-wide Phase Comparison", 
                                os.path.join(plots_dir, "genome_wide_phase_comparison.png"))
    
    # 2. diff heatmap
    create_diff_heatmap(neu_results, nsc_results,
                        "Genome-wide Signal diff Heatmap",
                        os.path.join(plots_dir, "genome_wide_diff_heatmap.png"))
    
    # 3. Chromosome distribution plot
    create_chromosome_distribution_plot(neu_results, nsc_results,
                                       os.path.join(plots_dir, "chromosome_distribution.png"))
    
    logging.info(f"Genome-wide plots saved to {plots_dir}")

def create_phase_comparison_plot(neu_results, nsc_results, title, output_file):
    """Create a scatter plot comparing S2S and S3 signals across cell types.

    Args:
        neu_results (pd.DataFrame): DataFrame containing results for neurons.
        nsc_results (pd.DataFrame): DataFrame containing results for NSCs.
        title (str): The title of the plot.
        output_file (str): The path to save the plot.
    """
    plt.figure(figsize=(12, 10))
    
    # Neurons
    plt.subplot(2, 2, 1)  # Create a subplot in a 2x2 grid, this is the first subplot
    plt.scatter(neu_results['s2s_signal'], neu_results['s3_signal'], alpha=0.5)  # Scatter plot of S2S vs S3 signal
    plt.title('Neurons')  # Set the title of the subplot
    plt.xlabel('S2S Signal')  # Set the x-axis label
    plt.ylabel('S3 Signal')  # Set the y-axis label
    plt.grid(True, alpha=0.3)  # Add a grid for better readability
    
    # NSCs
    plt.subplot(2, 2, 2)  # Create a subplot in a 2x2 grid, this is the second subplot
    plt.scatter(nsc_results['s2s_signal'], nsc_results['s3_signal'], alpha=0.5)  # Scatter plot of S2S vs S3 signal
    plt.title('NSCs')  # Set the title of the subplot
    plt.xlabel('S2S Signal')  # Set the x-axis label
    plt.ylabel('S3 Signal')  # Set the y-axis label
    plt.grid(True, alpha=0.3)  # Add a grid for better readability
    
    # Combined with different colors
    plt.subplot(2, 1, 2)  # Create a subplot in a 2x1 grid, this is the second subplot (occupies the bottom half)
    plt.scatter(neu_results['s2s_signal'], neu_results['s3_signal'], alpha=0.5, label='Neurons', color='blue')  # Scatter plot for Neurons in blue
    plt.scatter(nsc_results['s2s_signal'], nsc_results['s3_signal'], alpha=0.5, label='NSCs', color='red')  # Scatter plot for NSCs in red
    plt.title('Combined Comparison')  # Set the title of the subplot
    plt.xlabel('S2S Signal')  # Set the x-axis label
    plt.ylabel('S3 Signal')  # Set the y-axis label
    plt.legend()  # Show the legend to distinguish between cell types
    plt.grid(True, alpha=0.3)  # Add a grid for better readability
    
    plt.suptitle(title)  # Set the overall title for the figure
    plt.tight_layout()  # Adjust subplot parameters for a tight layout
    plt.savefig(output_file, dpi=300)  # Save the figure to a file with a specified DPI
    plt.close()  # Close the figure to free memory

def create_diff_heatmap(neu_results, nsc_results, title, output_file):
    """Create a heatmap of signal diffs between conditions.

    Args:
        neu_results (pd.DataFrame): DataFrame containing results for neurons.
        nsc_results (pd.DataFrame): DataFrame containing results for NSCs.
        title (str): The title of the heatmap.
        output_file (str): The path to save the heatmap.
    """
    # Calculate diffs of S3 signal to S2S signal for each cell type
    neu_diff = neu_results['s3_signal'] - neu_results['s2s_signal']
    nsc_diff = nsc_results['s3_signal'] - nsc_results['s2s_signal']
    
    # Replace infinities and NaNs with 0 to handle division by zero or missing data
    # neu_diff = neu_diff.replace([np.inf, -np.inf], np.nan).fillna(0)
    # nsc_diff = nsc_diff.replace([np.inf, -np.inf], np.nan).fillna(0)
    
    # Create a DataFrame for the heatmap, with diffs for Neurons and NSCs
    diff_df = pd.DataFrame({
        'Neurons': neu_diff,
        'NSCs': nsc_diff
    })
    
    # Sort the DataFrame by the neuron diff for better visualization
    diff_df = diff_df.sort_values('Neurons')
    
    # Create the heatmap
    plt.figure(figsize=(10, 8))  # Set the figure size
    sns.heatmap(diff_df.T, cmap='viridis', cbar_kws={'label': 'S3-S2S diff'})  # Create the heatmap using seaborn
    plt.title(title)  # Set the title of the heatmap
    plt.tight_layout()  # Adjust subplot parameters for a tight layout
    plt.savefig(output_file, dpi=300)  # Save the heatmap to a file
    plt.close()  # Close the figure

def create_chromosome_distribution_plot(neu_results, nsc_results, output_file):
    """Create a plot showing the distribution of genes across chromosomes.

    Args:
        neu_results (pd.DataFrame): DataFrame containing results for neurons.
        nsc_results (pd.DataFrame): DataFrame containing results for NSCs.
        output_file (str): The path to save the plot.
    """
    # Count the number of genes per chromosome for each cell type
    neu_chrom_counts = neu_results['chromosome'].value_counts().sort_index()
    nsc_chrom_counts = nsc_results['chromosome'].value_counts().sort_index()
    
    # Combine the counts into a single DataFrame
    chrom_df = pd.DataFrame({
        'Neurons': neu_chrom_counts,
        'NSCs': nsc_chrom_counts
    })
    
    # Fill NaN values with zeros (in case a chromosome has no genes in one cell type)
    chrom_df = chrom_df.fillna(0)
    
    # Create the plot
    plt.figure(figsize=(14, 8))  # Set the figure size
    chrom_df.plot(kind='bar')  # Create a bar plot of the chromosome distribution
    plt.title('Gene Distribution Across Chromosomes')  # Set the title of the plot
    plt.xlabel('Chromosome')  # Set the x-axis label
    plt.ylabel('Number of Genes')  # Set the y-axis label
    plt.grid(True, alpha=0.3)  # Add a grid for better readability
    plt.legend()  # Show the legend
    plt.tight_layout()  # Adjust subplot parameters for a tight layout
    plt.savefig(output_file, dpi=300)  # Save the plot to a file
    plt.close()  # Close the figure

def compare_cell_types(neu_results, nsc_results, output_dir):
    """Compare results between cell types and generate comparison metrics.

    Args:
        neu_results (pd.DataFrame): DataFrame containing results for neurons.
        nsc_results (pd.DataFrame): DataFrame containing results for NSCs.
        output_dir (str): The directory to save the comparison results.

    Returns:
        pd.DataFrame: A DataFrame containing the merged results and comparison metrics.
    """
    logging.info("Comparing results between cell types")
    
    # Create a merged dataset with genes present in both cell types
    # The merge is based on the 'gene_id' column, and suffixes are added to distinguish columns from each cell type
    merged_df = pd.merge(
        neu_results, 
        nsc_results,
        on='gene_id',
        suffixes=('_neu', '_nsc')
    )
    
    # Calculate differential metrics between cell types
    merged_df['s2s_diff'] = merged_df['s2s_signal_neu'] - merged_df['s2s_signal_nsc']  # Difference in S2S signal
    merged_df['s3_diff'] = merged_df['s3_signal_neu'] - merged_df['s3_signal_nsc']  # Difference in S3 signal
    merged_df['diff_diff'] = (merged_df['s3_signal_neu'] - merged_df['s2s_signal_neu']) - (merged_df['s3_signal_nsc'] - merged_df['s2s_signal_nsc'])  # Difference in S3-S2S diff
    
    # Save the comparison results to a CSV file
    comparison_file = os.path.join(output_dir, "cell_type_comparison.csv")
    merged_df.to_csv(comparison_file, index=False)
    
    logging.info(f"Cell type comparison saved to {comparison_file}")
    
    # Create comparison plots
    plots_dir = os.path.join(output_dir, "plots")
    
    # Differential signal plot: scatter plots of S2S and S3 signals for Neurons vs NSCs
    plt.figure(figsize=(12, 10))
    plt.subplot(2, 1, 1)  # Create the first subplot (top)
    plt.scatter(merged_df['s2s_signal_neu'], merged_df['s2s_signal_nsc'], alpha=0.5)  # Scatter plot of S2S signals
    plt.title('S2S Signal: Neurons vs NSCs')  # Set the title
    plt.xlabel('Neurons S2S Signal')  # Set the x-axis label
    plt.ylabel('NSCs S2S Signal')  # Set the y-axis label
    plt.grid(True, alpha=0.3)  # Add a grid
    
    plt.subplot(2, 1, 2)  # Create the second subplot (bottom)
    plt.scatter(merged_df['s3_signal_neu'], merged_df['s3_signal_nsc'], alpha=0.5)  # Scatter plot of S3 signals
    plt.title('S3 Signal: Neurons vs NSCs')  # Set the title
    plt.xlabel('Neurons S3 Signal')  # Set the x-axis label
    plt.ylabel('NSCs S3 Signal')  # Set the y-axis label
    plt.grid(True, alpha=0.3)  # Add a grid
    
    plt.tight_layout()  # Adjust subplot parameters for a tight layout
    plt.savefig(os.path.join(plots_dir, "cell_type_signal_comparison.png"), dpi=300)  # Save the plot
    plt.close()  # Close the figure
    
    return merged_df  # Return the merged DataFrame

def main():
    """Main function to execute the script."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Merge chromosome-specific results from endogenous Mecp2 target gene analysis')
    parser.add_argument('--results_dir', default='../results/endogenous_target_analysis', help='Directory containing chromosome-specific results')
    parser.add_argument('--output_dir', default=None, help='Directory for merged results (defaults to ../results_dir/endogenous_target_analysis/merged)')
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
    
    # Check if any chromosome directories were found
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
    main()  # Execute the main function when the script is run