#!/usr/bin/env python3
"""
SAMMY-seq Chromatin Transition States Visualization

This script creates a stacked bar chart showing the distribution of chromatin states,
including the transition zones between heterochromatin and euchromatin.

Usage:
    python plot_transition_states.py --results-dir /path/to/results --output-file chromatin_states.png
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='SAMMY-seq Chromatin Transition States Visualization')
    parser.add_argument('--results-dir', required=True, help='Directory containing SAMMY-seq analysis results')
    parser.add_argument('--output-file', default='chromatin_transition_states.png', help='Output file name')
    parser.add_argument('--hetero-threshold', type=float, default=-0.1, help='Log2 ratio threshold for heterochromatin (default: -0.1)')
    parser.add_argument('--eu-threshold', type=float, default=0.1, help='Log2 ratio threshold for euchromatin (default: 0.1)')
    return parser.parse_args()

def calculate_transition_states(results_dir, hetero_threshold, eu_threshold):
    """
    Calculate the proportion of the genome in each chromatin state,
    including transition states
    """
    # Find all bedgraph files
    bedgraph_files = glob.glob(os.path.join(results_dir, "bedgraph", "*_S2S_vs_S3_ratio.bedgraph"))
    
    if not bedgraph_files:
        print(f"No bedgraph files found in {os.path.join(results_dir, 'bedgraph')}")
        return None
    
    results = {}
    
    for bedgraph_file in bedgraph_files:
        # Extract condition name from filename
        filename = os.path.basename(bedgraph_file)
        condition_match = re.match(r'(.+)_S2S_vs_S3_ratio\.bedgraph', filename)
        if not condition_match:
            continue
        
        condition = condition_match.group(1)
        print(f"Processing {condition}...")
        
        # Read bedgraph file
        try:
            df = pd.read_csv(bedgraph_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
        except Exception as e:
            print(f"Error reading {bedgraph_file}: {e}")
            continue
        
        # Calculate region sizes
        df['size'] = df['end'] - df['start']
        total_size = df['size'].sum()
        
        # Classify regions
        state_A = df[df['score'] <= hetero_threshold]  # Heterochromatin
        state_B = df[df['score'] >= eu_threshold]      # Euchromatin
        state_A_to_B = df[(df['score'] > hetero_threshold) & (df['score'] < 0)]  # Transition to euchromatin
        state_B_to_A = df[(df['score'] >= 0) & (df['score'] < eu_threshold)]     # Transition to heterochromatin
        
        # Calculate percentages
        results[condition] = {
            'A': state_A['size'].sum() / total_size * 100,  # Heterochromatin
            'B': state_B['size'].sum() / total_size * 100,  # Euchromatin
            'A->B': state_A_to_B['size'].sum() / total_size * 100,  # Transition to euchromatin
            'B->A': state_B_to_A['size'].sum() / total_size * 100   # Transition to heterochromatin
        }
        
        # Also store sample group information
        if '_GFP' in condition:
            results[condition]['group'] = 'GFP'
        elif '_M2' in condition:
            results[condition]['group'] = 'M2'
        
        # Get sample number
        if re.search(r'[1-3]', condition):
            sample_num = re.search(r'[1-3]', condition).group(0)
            results[condition]['sample_num'] = sample_num
        else:
            results[condition]['sample_num'] = ''
            
        # Get cell type
        if condition.startswith('NSC'):
            results[condition]['cell_type'] = 'NSC'
        elif condition.startswith('Neu'):
            results[condition]['cell_type'] = 'Neu'
        else:
            results[condition]['cell_type'] = ''
            
    return results

def plot_transition_states(results, output_file):
    """Create a stacked bar chart showing chromatin transition states"""
    if not results:
        print("No results to plot")
        return
    
    # Convert results to DataFrame
    df = pd.DataFrame.from_dict(results, orient='index')
    
    # Sort by cell type, group, and sample number
    df['sort_key'] = df.apply(lambda x: f"{x['cell_type']}_{x['group']}_{x['sample_num']}", axis=1)
    df = df.sort_values('sort_key')
    
    # Drop extra columns for plotting
    plot_df = df[['A', 'A->B', 'B->A', 'B']]
    
    # Create figure
    plt.figure(figsize=(10, 8))
    
    # Create stacked bar plot
    ax = plot_df.plot(kind='barh', stacked=True, figsize=(10, 8), 
                     color=['#4daf4a', '#b3de69', '#fccde5', '#dfc27d'], 
                     width=0.8)
    
    # Customize plot
    plt.xlim(0, 100)
    plt.xlabel('Percentage of genome')
    plt.ylabel('')
    plt.title('Chromatin State Distribution', fontsize=14)
    
    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#4daf4a', label='A (Heterochromatin)'),
        Patch(facecolor='#b3de69', label='A→B (Transition to Euchromatin)'),
        Patch(facecolor='#fccde5', label='B→A (Transition to Heterochromatin)'),
        Patch(facecolor='#dfc27d', label='B (Euchromatin)')
    ]
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Tight layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    
    # Also save data to TSV
    tsv_file = os.path.splitext(output_file)[0] + '.tsv'
    plot_df.to_csv(tsv_file, sep='\t')
    print(f"Data saved to {tsv_file}")
    
    plt.close()

def main():
    """Main function"""
    args = parse_arguments()
    
    # Calculate chromatin transition states
    results = calculate_transition_states(
        args.results_dir, 
        args.hetero_threshold, 
        args.eu_threshold
    )
    
    # Plot results
    if results:
        plot_transition_states(results, args.output_file)
    else:
        print("No results to plot")

if __name__ == "__main__":
    main()