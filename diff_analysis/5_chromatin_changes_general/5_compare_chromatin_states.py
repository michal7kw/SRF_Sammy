#!/usr/bin/env python3
"""
This script compares chromatin states between different cell types and conditions:
1. NSCs (GFP) vs Neurons (GFP) - to show how many regions switch from A (euchromatin - S2S) to B (heterochromatin - S3)
2. NSCs (GFP) vs NSCs (M2) - to demonstrate that only a few regions modify their chromatin state
3. Neurons (GFP) vs Neurons (M2) - to demonstrate that only a few regions modify their chromatin state

The script generates stacked bar charts showing the proportion of regions in different chromatin states.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import sys
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('chromatin_comparison.log')
    ]
)

# Constants
MIN_DIFF_THRESHOLD = 2.0  # Same threshold as in the original analysis

def load_data(analysis_dir):
    """
    Load the chromatin state data for NSCs and Neurons from the analysis directory.
    
    Args:
        analysis_dir (str): Path to the directory containing analysis results
        
    Returns:
        tuple: DataFrames for NSC and Neuron chromatin states
    """
    analysis_dir = Path(analysis_dir)
    
    # Check if the analysis directory exists
    if not analysis_dir.exists():
        logging.error(f"Analysis directory {analysis_dir} does not exist")
        return None, None
    
    # Load NSC data
    nsc_file = analysis_dir / "NSC_chromatin_changes_all.csv"
    if not nsc_file.exists():
        logging.error(f"NSC data file {nsc_file} does not exist")
        return None, None
    
    # Load Neuron data
    neu_file = analysis_dir / "Neu_chromatin_changes_all.csv"
    if not neu_file.exists():
        logging.error(f"Neuron data file {neu_file} does not exist")
        return None, None
    
    logging.info(f"Loading NSC data from {nsc_file}")
    nsc_df = pd.read_csv(nsc_file)
    
    logging.info(f"Loading Neuron data from {neu_file}")
    neu_df = pd.read_csv(neu_file)
    
    return nsc_df, neu_df

def classify_regions(df, condition='gfp'):
    """
    Classify genomic regions based on their chromatin state.
    
    Args:
        df (pandas.DataFrame): DataFrame containing chromatin state data
        condition (str): Which condition to analyze ('gfp' or 'm2')
        
    Returns:
        pandas.DataFrame: DataFrame with added classification column
    """
    # Create a copy of the DataFrame to avoid modifying the original
    result_df = df.copy()
    
    # Get the column names for the condition
    s2s_col = f'{condition}_s2s'
    s3_col = f'{condition}_s3'
    diff_col = f'{condition}_diff'
    
    # Classify regions based on the difference between S3 and S2S signals
    # A: Euchromatin (S2S > S3)
    # B: Heterochromatin (S3 > S2S)
    # A->B: Transition from euchromatin to heterochromatin
    # B->A: Transition from heterochromatin to euchromatin
    
    # Initialize the classification column
    result_df[f'{condition}_class'] = 'Unknown'
    
    # Classify as A (euchromatin) if S2S signal is higher than S3
    result_df.loc[result_df[diff_col] < -MIN_DIFF_THRESHOLD, f'{condition}_class'] = 'A'
    
    # Classify as B (heterochromatin) if S3 signal is higher than S2S
    result_df.loc[result_df[diff_col] > MIN_DIFF_THRESHOLD, f'{condition}_class'] = 'B'
    
    # Classify as intermediate if the difference is small
    result_df.loc[(result_df[diff_col] >= -MIN_DIFF_THRESHOLD) & 
                 (result_df[diff_col] <= MIN_DIFF_THRESHOLD), f'{condition}_class'] = 'Intermediate'
    
    return result_df

def compare_cell_types(nsc_df, neu_df):
    """
    Compare chromatin states between NSCs and Neurons in GFP condition.
    
    Args:
        nsc_df (pandas.DataFrame): DataFrame containing NSC chromatin state data
        neu_df (pandas.DataFrame): DataFrame containing Neuron chromatin state data
        
    Returns:
        pandas.DataFrame: DataFrame with comparison results
    """
    # Classify regions for both cell types
    nsc_classified = classify_regions(nsc_df, 'gfp')
    neu_classified = classify_regions(neu_df, 'gfp')
    
    # Merge the DataFrames on chromosome, start, and end
    merged_df = pd.merge(
        nsc_classified[['chrom', 'start', 'end', 'gfp_class']],
        neu_classified[['chrom', 'start', 'end', 'gfp_class']],
        on=['chrom', 'start', 'end'],
        suffixes=('_nsc', '_neu')
    )
    
    # Classify the transitions
    merged_df['transition'] = 'No change'
    
    # A to B transition (euchromatin to heterochromatin)
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'A') & 
                 (merged_df['gfp_class_neu'] == 'B'), 'transition'] = 'A->B'
    
    # B to A transition (heterochromatin to euchromatin)
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'B') & 
                 (merged_df['gfp_class_neu'] == 'A'), 'transition'] = 'B->A'
    
    # A to A (no change in euchromatin)
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'A') & 
                 (merged_df['gfp_class_neu'] == 'A'), 'transition'] = 'A'
    
    # B to B (no change in heterochromatin)
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'B') & 
                 (merged_df['gfp_class_neu'] == 'B'), 'transition'] = 'B'
    
    # Handle intermediate states
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'Intermediate') & 
                 (merged_df['gfp_class_neu'] == 'A'), 'transition'] = 'I->A'
    
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'Intermediate') & 
                 (merged_df['gfp_class_neu'] == 'B'), 'transition'] = 'I->B'
    
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'A') & 
                 (merged_df['gfp_class_neu'] == 'Intermediate'), 'transition'] = 'A->I'
    
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'B') & 
                 (merged_df['gfp_class_neu'] == 'Intermediate'), 'transition'] = 'B->I'
    
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'Intermediate') & 
                 (merged_df['gfp_class_neu'] == 'Intermediate'), 'transition'] = 'I'
    
    return merged_df

def compare_conditions(df, cell_type):
    """
    Compare chromatin states between GFP and M2 conditions for a given cell type.
    
    Args:
        df (pandas.DataFrame): DataFrame containing chromatin state data
        cell_type (str): Cell type name for labeling
        
    Returns:
        pandas.DataFrame: DataFrame with comparison results
    """
    # Classify regions for both conditions
    classified_df = df.copy()
    classified_df = classify_regions(classified_df, 'gfp')
    classified_df = classify_regions(classified_df, 'm2')
    
    # Classify the transitions
    classified_df['transition'] = 'No change'
    
    # A to B transition (euchromatin to heterochromatin)
    classified_df.loc[(classified_df['gfp_class'] == 'A') & 
                     (classified_df['m2_class'] == 'B'), 'transition'] = 'A->B'
    
    # B to A transition (heterochromatin to euchromatin)
    classified_df.loc[(classified_df['gfp_class'] == 'B') & 
                     (classified_df['m2_class'] == 'A'), 'transition'] = 'B->A'
    
    # A to A (no change in euchromatin)
    classified_df.loc[(classified_df['gfp_class'] == 'A') & 
                     (classified_df['m2_class'] == 'A'), 'transition'] = 'A'
    
    # B to B (no change in heterochromatin)
    classified_df.loc[(classified_df['gfp_class'] == 'B') & 
                     (classified_df['m2_class'] == 'B'), 'transition'] = 'B'
    
    # Handle intermediate states
    classified_df.loc[(classified_df['gfp_class'] == 'Intermediate') & 
                     (classified_df['m2_class'] == 'A'), 'transition'] = 'I->A'
    
    classified_df.loc[(classified_df['gfp_class'] == 'Intermediate') & 
                     (classified_df['m2_class'] == 'B'), 'transition'] = 'I->B'
    
    classified_df.loc[(classified_df['gfp_class'] == 'A') & 
                     (classified_df['m2_class'] == 'Intermediate'), 'transition'] = 'A->I'
    
    classified_df.loc[(classified_df['gfp_class'] == 'B') & 
                     (classified_df['m2_class'] == 'Intermediate'), 'transition'] = 'B->I'
    
    classified_df.loc[(classified_df['gfp_class'] == 'Intermediate') & 
                     (classified_df['m2_class'] == 'Intermediate'), 'transition'] = 'I'
    
    # Add cell type information
    classified_df['cell_type'] = cell_type
    
    return classified_df

def plot_chromatin_transitions(comparison_dfs, output_dir):
    """
    Create stacked bar charts showing chromatin state transitions.
    
    Args:
        comparison_dfs (list): List of DataFrames with comparison results
        output_dir (str): Directory to save the plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Combine all comparison DataFrames
    combined_df = pd.concat(comparison_dfs)
    
    # Create a pivot table to count the number of regions in each transition category
    pivot_df = combined_df.pivot_table(
        index='cell_type',
        columns='transition',
        values='chrom',
        aggfunc='count',
        fill_value=0
    )
    
    # Calculate the percentage of regions in each category
    pivot_pct = pivot_df.div(pivot_df.sum(axis=1), axis=0) * 100
    
    # Define the order of transitions for the plot (to match the example image)
    # We want to prioritize the main categories: A, B, A->B, B->A
    transition_order = ['B', 'B->A', 'A->B', 'A']
    
    # Filter columns to include only the transitions in the order
    # and ensure all columns exist (some might be missing in the data)
    available_columns = [col for col in transition_order if col in pivot_pct.columns]
    pivot_pct = pivot_pct[available_columns]
    
    # Define colors for each transition category to match the example image
    colors = {
        'A': '#2a9d8f',       # Teal for euchromatin (A)
        'B': '#e9c46a',       # Yellow/gold for heterochromatin (B)
        'A->B': '#f4a261',    # Orange for euchromatin to heterochromatin (A->B)
        'B->A': '#a8dadc',    # Light blue for heterochromatin to euchromatin (B->A)
    }
    
    # Create the stacked bar chart
    plt.figure(figsize=(12, 8))
    
    # Plot the stacked bars
    ax = pivot_pct.plot(
        kind='barh',
        stacked=True,
        color=[colors.get(x, '#333333') for x in pivot_pct.columns],
        figsize=(12, 8),
        width=0.7
    )
    
    # Customize the plot to match the example
    plt.title('Neu 10K', fontsize=16)
    plt.xlabel('Percentage of genome', fontsize=14)
    plt.ylabel('', fontsize=14)
    
    # Calculate the maximum percentage value to set appropriate x-axis limits
    max_total_pct = pivot_pct.sum(axis=1).max()
    # Add a 10% buffer to the maximum value for better visualization
    x_limit = min(max_total_pct * 1.1, 100)
    # Ensure the x-axis limit is at least 20% for better visualization
    x_limit = max(x_limit, 20)
    # Round up to the nearest 5 for cleaner tick marks
    x_limit = np.ceil(x_limit / 5) * 5
    
    plt.xlim(0, x_limit)
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    
    # Add x-axis ticks at regular intervals
    tick_interval = 5 if x_limit <= 30 else 10
    plt.xticks(np.arange(0, x_limit + 1, tick_interval))
    
    # Remove the frame on the right and top
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Add a legend on the right side
    # Create a custom legend with the correct symbols and labels
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor=colors['B'], label='B'),
        plt.Rectangle((0, 0), 1, 1, facecolor=colors['B->A'], label='B→A'),
        plt.Rectangle((0, 0), 1, 1, facecolor=colors['A->B'], label='A→B'),
        plt.Rectangle((0, 0), 1, 1, facecolor=colors['A'], label='A')
    ]
    
    # Place the legend to the right of the plot
    plt.legend(handles=legend_elements, loc='center right', bbox_to_anchor=(1.15, 0.5), fontsize=12)
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_dir / 'chromatin_transitions.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'chromatin_transitions.pdf', bbox_inches='tight')
    
    logging.info(f"Plot saved to {output_dir / 'chromatin_transitions.png'}")
    
    # Also save the data used for the plot
    pivot_pct.to_csv(output_dir / 'chromatin_transitions_percentages.csv')
    logging.info(f"Data saved to {output_dir / 'chromatin_transitions_percentages.csv'}")
    
    # Print the actual percentages for debugging
    logging.info("Percentage of genome in each category:")
    for cell_type, row in pivot_pct.iterrows():
        logging.info(f"{cell_type}: {dict(row)}")
        logging.info(f"Total: {row.sum():.2f}%")

def main():
    """Main function to execute the chromatin state comparison."""
    # Add command-line argument parsing
    parser = argparse.ArgumentParser(description='Compare chromatin states between cell types and conditions')
    parser.add_argument('--analysis_dir', default='../results/analysis', help='Directory containing analysis results')
    parser.add_argument('--output_dir', default='../results/comparison', help='Directory to save comparison results')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging.info("Starting chromatin state comparison")
    
    # Load the data
    nsc_df, neu_df = load_data(args.analysis_dir)
    
    if nsc_df is None or neu_df is None:
        logging.error("Failed to load data. Exiting.")
        return
    
    # Compare NSCs vs Neurons (GFP condition)
    logging.info("Comparing NSCs vs Neurons (GFP condition)")
    nsc_vs_neu = compare_cell_types(nsc_df, neu_df)
    # Set a descriptive label for the NSC vs Neuron comparison
    nsc_vs_neu['cell_type'] = 'NSC_vs_Neu'
    
    # Create comparisons for GFP vs M2 in both cell types
    comparisons = []
    
    # NSC comparisons
    logging.info("Comparing NSC: GFP vs M2")
    nsc_gfp_vs_m2 = compare_conditions(nsc_df, 'NSC_GFP')
    comparisons.append(nsc_gfp_vs_m2)
    
    # Neuron comparisons
    logging.info("Comparing Neurons: GFP vs M2")
    neu_gfp_vs_m2 = compare_conditions(neu_df, 'Neu_GFP')
    comparisons.append(neu_gfp_vs_m2)
    
    # Combine all comparisons for plotting
    all_comparisons = [nsc_vs_neu] + comparisons
    
    # Plot the results
    logging.info("Creating plots")
    plot_chromatin_transitions(all_comparisons, args.output_dir)
    
    logging.info("Chromatin state comparison completed")

if __name__ == "__main__":
    main() 