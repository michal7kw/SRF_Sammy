#!/usr/bin/env python3
"""
This script specifically compares chromatin states between NSCs and Neurons in GFP condition.
It generates a detailed visualization showing how many regions switch from A (euchromatin - S2S) 
to B (heterochromatin - S3) during neuronal differentiation.

The script uses 10K windows for Sammy-seq signal aggregation.
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
        logging.FileHandler('nsc_neuron_comparison.log')
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

def compare_nsc_neuron(nsc_df, neu_df, output_dir):
    """
    Compare chromatin states between NSCs and Neurons in GFP condition.
    
    Args:
        nsc_df (pandas.DataFrame): DataFrame containing NSC chromatin state data
        neu_df (pandas.DataFrame): DataFrame containing Neuron chromatin state data
        output_dir (str): Directory to save the results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Classify regions for both cell types
    nsc_classified = classify_regions(nsc_df, 'gfp')
    neu_classified = classify_regions(neu_df, 'gfp')
    
    # Log the distribution of chromatin states in each cell type
    logging.info("Chromatin state distribution in NSCs:")
    nsc_class_counts = nsc_classified['gfp_class'].value_counts()
    nsc_class_percentages = (nsc_class_counts / len(nsc_classified) * 100).round(2)
    for state, percentage in nsc_class_percentages.items():
        logging.info(f"  {state}: {percentage:.2f}%")
    
    logging.info("Chromatin state distribution in Neurons:")
    neu_class_counts = neu_classified['gfp_class'].value_counts()
    neu_class_percentages = (neu_class_counts / len(neu_classified) * 100).round(2)
    for state, percentage in neu_class_percentages.items():
        logging.info(f"  {state}: {percentage:.2f}%")
    
    # Merge the DataFrames on chromosome, start, and end
    merged_df = pd.merge(
        nsc_classified[['chrom', 'start', 'end', 'gfp_class', 'gfp_s2s', 'gfp_s3']],
        neu_classified[['chrom', 'start', 'end', 'gfp_class', 'gfp_s2s', 'gfp_s3']],
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
    merged_df.loc[(merged_df['gfp_class_nsc'] == 'Intermediate') | 
                 (merged_df['gfp_class_neu'] == 'Intermediate'), 'transition'] = 'Intermediate'
    
    # Save the merged DataFrame
    merged_df.to_csv(output_dir / 'nsc_neuron_comparison.csv', index=False)
    
    # Count the number of regions in each transition category
    transition_counts = merged_df['transition'].value_counts()
    
    # Calculate percentages
    transition_percentages = (transition_counts / len(merged_df) * 100).round(2)
    
    # Log the transition percentages
    logging.info("Transition percentages between NSCs and Neurons:")
    for transition, percentage in transition_percentages.items():
        logging.info(f"  {transition}: {percentage:.2f}%")
    
    # Save the transition counts and percentages
    transition_df = pd.DataFrame({
        'count': transition_counts,
        'percentage': transition_percentages
    })
    transition_df.to_csv(output_dir / 'nsc_neuron_transition_counts.csv')
    
    # Create a pie chart of the transitions
    plt.figure(figsize=(10, 8))
    
    # Define colors for each transition category
    colors = {
        'A': '#2a9d8f',       # Teal for euchromatin (A)
        'B': '#e9c46a',       # Yellow/gold for heterochromatin (B)
        'A->B': '#f4a261',    # Orange for euchromatin to heterochromatin (A->B)
        'B->A': '#a8dadc',    # Light blue for heterochromatin to euchromatin (B->A)
        'Intermediate': '#e5e5e5'  # Light gray for intermediate
    }
    
    # Filter out transitions with very small percentages
    min_percentage = 1.0  # Minimum percentage to include in the pie chart
    filtered_transitions = transition_percentages[transition_percentages >= min_percentage]
    
    # Create the pie chart
    plt.pie(
        filtered_transitions,
        labels=filtered_transitions.index,
        autopct='%1.1f%%',
        colors=[colors.get(x, '#333333') for x in filtered_transitions.index],
        startangle=90,
        wedgeprops={'edgecolor': 'w', 'linewidth': 1}
    )
    
    plt.title('Chromatin State Transitions: NSC to Neuron (10K windows)', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_dir / 'nsc_neuron_transitions_pie.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'nsc_neuron_transitions_pie.pdf', bbox_inches='tight')
    
    # Create a bar chart showing the signal distribution for each transition type
    plt.figure(figsize=(12, 8))
    
    # Calculate the mean signal for each transition type
    signal_by_transition = merged_df.groupby('transition').agg({
        'gfp_s2s_nsc': 'mean',
        'gfp_s3_nsc': 'mean',
        'gfp_s2s_neu': 'mean',
        'gfp_s3_neu': 'mean'
    }).reset_index()
    
    # Log the mean signal values
    logging.info("Mean signal values by transition type:")
    for _, row in signal_by_transition.iterrows():
        logging.info(f"  {row['transition']}:")
        logging.info(f"    NSC S2S: {row['gfp_s2s_nsc']:.2f}")
        logging.info(f"    NSC S3: {row['gfp_s3_nsc']:.2f}")
        logging.info(f"    Neuron S2S: {row['gfp_s2s_neu']:.2f}")
        logging.info(f"    Neuron S3: {row['gfp_s3_neu']:.2f}")
    
    # Melt the DataFrame for easier plotting
    signal_melted = pd.melt(
        signal_by_transition,
        id_vars=['transition'],
        value_vars=['gfp_s2s_nsc', 'gfp_s3_nsc', 'gfp_s2s_neu', 'gfp_s3_neu'],
        var_name='signal_type',
        value_name='signal'
    )
    
    # Create labels for the signal types
    signal_melted['cell_type'] = signal_melted['signal_type'].apply(
        lambda x: 'NSC' if 'nsc' in x else 'Neuron'
    )
    signal_melted['mark'] = signal_melted['signal_type'].apply(
        lambda x: 'S2S (Euchromatin)' if 's2s' in x else 'S3 (Heterochromatin)'
    )
    
    # Create the grouped bar chart
    sns.barplot(
        data=signal_melted,
        x='transition',
        y='signal',
        hue='mark',
        palette={'S2S (Euchromatin)': '#2a9d8f', 'S3 (Heterochromatin)': '#e9c46a'},
        errorbar=None,
        alpha=0.7
    )
    
    plt.title('Average Signal by Transition Type: NSC to Neuron', fontsize=16)
    plt.xlabel('Transition Type', fontsize=14)
    plt.ylabel('Average Signal', fontsize=14)
    plt.xticks(rotation=45)
    plt.legend(title='Mark')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / 'nsc_neuron_signal_by_transition.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'nsc_neuron_signal_by_transition.pdf', bbox_inches='tight')
    
    # Create a stacked bar chart showing the percentage of regions in each transition category
    plt.figure(figsize=(12, 8))
    
    # Create a DataFrame with the transition percentages
    transition_df = pd.DataFrame({
        'transition': transition_percentages.index,
        'percentage': transition_percentages.values
    })
    
    # Define the order of transitions for the plot
    transition_order = ['B', 'B->A', 'A->B', 'A', 'Intermediate']
    
    # Filter and order the transitions
    available_transitions = [t for t in transition_order if t in transition_df['transition'].values]
    transition_df = transition_df[transition_df['transition'].isin(available_transitions)]
    transition_df['transition'] = pd.Categorical(
        transition_df['transition'], 
        categories=available_transitions, 
        ordered=True
    )
    transition_df = transition_df.sort_values('transition')
    
    # Create the bar chart
    ax = sns.barplot(
        data=transition_df,
        x='percentage',
        y='transition',
        palette={
            'A': '#2a9d8f',
            'B': '#e9c46a',
            'A->B': '#f4a261',
            'B->A': '#a8dadc',
            'Intermediate': '#e5e5e5'
        },
        orient='h'
    )
    
    # Calculate the maximum percentage value to set appropriate x-axis limits
    max_pct = transition_df['percentage'].max()
    # Add a 10% buffer to the maximum value for better visualization
    x_limit = min(max_pct * 1.1, 100)
    # Ensure the x-axis limit is at least 20% for better visualization
    x_limit = max(x_limit, 20)
    # Round up to the nearest 5 for cleaner tick marks
    x_limit = np.ceil(x_limit / 5) * 5
    
    plt.xlim(0, x_limit)
    
    # Add x-axis ticks at regular intervals
    tick_interval = 5 if x_limit <= 30 else 10
    plt.xticks(np.arange(0, x_limit + 1, tick_interval))
    
    plt.title('Chromatin State Transitions: NSC to Neuron (10K windows)', fontsize=16)
    plt.xlabel('Percentage of Regions', fontsize=14)
    plt.ylabel('Transition Type', fontsize=14)
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    
    # Add percentage labels to the bars
    for i, p in enumerate(ax.patches):
        width = p.get_width()
        plt.text(
            width + 0.5,
            p.get_y() + p.get_height() / 2,
            f'{width:.1f}%',
            ha='left',
            va='center'
        )
    
    plt.tight_layout()
    plt.savefig(output_dir / 'nsc_neuron_transitions_bar.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'nsc_neuron_transitions_bar.pdf', bbox_inches='tight')
    
    # Create a heatmap showing the distribution of transitions across chromosomes
    plt.figure(figsize=(14, 10))
    
    # Count transitions by chromosome
    chrom_transitions = pd.crosstab(
        merged_df['chrom'],
        merged_df['transition'],
        normalize='index'
    ) * 100  # Convert to percentages
    
    # Sort chromosomes in a more natural order
    def chrom_sort_key(chrom):
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        try:
            return int(chrom)
        except ValueError:
            return float('inf') if chrom == 'Y' else float('inf') - 1  # X before Y
    
    sorted_chroms = sorted(chrom_transitions.index, key=chrom_sort_key)
    chrom_transitions = chrom_transitions.loc[sorted_chroms]
    
    # Create the heatmap
    sns.heatmap(
        chrom_transitions,
        cmap='YlGnBu',
        annot=True,
        fmt='.1f',
        cbar_kws={'label': 'Percentage of Regions'}
    )
    
    plt.title('Chromatin State Transitions by Chromosome: NSC to Neuron', fontsize=16)
    plt.xlabel('Transition Type', fontsize=14)
    plt.ylabel('Chromosome', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_dir / 'nsc_neuron_transitions_by_chromosome.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'nsc_neuron_transitions_by_chromosome.pdf', bbox_inches='tight')
    
    logging.info(f"NSC vs Neuron comparison results saved to {output_dir}")
    
    return merged_df

def main():
    """Main function to execute the NSC vs Neuron comparison."""
    # Add command-line argument parsing
    parser = argparse.ArgumentParser(description='Compare chromatin states between NSCs and Neurons')
    parser.add_argument('--analysis_dir', default='results/analysis', help='Directory containing analysis results')
    parser.add_argument('--output_dir', default='results/nsc_neuron_comparison', help='Directory to save comparison results')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging.info("Starting NSC vs Neuron comparison")
    
    # Load the data
    nsc_df, neu_df = load_data(args.analysis_dir)
    
    if nsc_df is None or neu_df is None:
        logging.error("Failed to load data. Exiting.")
        return
    
    # Compare NSCs vs Neurons
    compare_nsc_neuron(nsc_df, neu_df, args.output_dir)
    
    logging.info("NSC vs Neuron comparison completed")

if __name__ == "__main__":
    main() 