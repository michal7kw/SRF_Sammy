#!/usr/bin/env python3
"""
Visualize Chromatin State Changes Between Conditions

This script analyzes chromatin state bed files to determine how chromatin structure
changes between GFP and M2 conditions, showing the percentage of the genome that
transitions between heterochromatin and euchromatin states.

Usage:
    python chromatin_state_changes.py --results-dir /path/to/results --output-file changes.png
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import tempfile
import glob
import re

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Chromatin State Changes Between Conditions')
    parser.add_argument('--results-dir', required=True, help='Directory containing SAMMY-seq analysis results')
    parser.add_argument('--output-file', default='chromatin_state_changes.png', help='Output file name')
    parser.add_argument('--genome-size', type=int, default=2800000000, help='Genome size in bp (default: 2.8Gb for mouse)')
    return parser.parse_args()

def analyze_chromatin_changes(results_dir, genome_size):
    """
    Analyze changes in chromatin state between GFP and M2 conditions
    for both Neu and NSC cell types
    """
    # Find euchromatin and heterochromatin bed files
    chromatin_dir = os.path.join(results_dir, "chromatin_states")
    
    # Define the comparisons we want to make
    comparisons = [
        {
            'name': 'Neu_GFP vs Neu_M2',
            'gfp_hetero': os.path.join(chromatin_dir, 'Neu_GFP_heterochromatin.bed'),
            'gfp_eu': os.path.join(chromatin_dir, 'Neu_GFP_euchromatin.bed'),
            'm2_hetero': os.path.join(chromatin_dir, 'Neu_M2_heterochromatin.bed'),
            'm2_eu': os.path.join(chromatin_dir, 'Neu_M2_euchromatin.bed'),
        },
        {
            'name': 'NSC_GFP vs NSC_M2',
            'gfp_hetero': os.path.join(chromatin_dir, 'NSC_GFP_heterochromatin.bed'),
            'gfp_eu': os.path.join(chromatin_dir, 'NSC_GFP_euchromatin.bed'),
            'm2_hetero': os.path.join(chromatin_dir, 'NSC_M2_heterochromatin.bed'),
            'm2_eu': os.path.join(chromatin_dir, 'NSC_M2_euchromatin.bed'),
        }
    ]
    
    results = {}
    
    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmp_dir:
        for comp in comparisons:
            # Check if all required files exist
            files_exist = all(os.path.exists(comp[key]) for key in ['gfp_hetero', 'gfp_eu', 'm2_hetero', 'm2_eu'])
            if not files_exist:
                print(f"Skipping {comp['name']} - missing required files")
                continue
            
            # Initialize results for this comparison
            results[comp['name']] = {'A_to_B': 0, 'B_to_A': 0}
            
            print(f"Processing {comp['name']}...")
            
            # 1. Identify A_to_B: Regions that were heterochromatin (A) in GFP and became euchromatin (B) in M2
            a_to_b_file = os.path.join(tmp_dir, f"{comp['name'].replace(' ', '_')}_A_to_B.bed")
            cmd = f"bedtools intersect -a {comp['gfp_hetero']} -b {comp['m2_eu']} -wa > {a_to_b_file}"
            subprocess.run(cmd, shell=True, check=True)
            
            # 2. Identify B_to_A: Regions that were euchromatin (B) in GFP and became heterochromatin (A) in M2
            b_to_a_file = os.path.join(tmp_dir, f"{comp['name'].replace(' ', '_')}_B_to_A.bed")
            cmd = f"bedtools intersect -a {comp['gfp_eu']} -b {comp['m2_hetero']} -wa > {b_to_a_file}"
            subprocess.run(cmd, shell=True, check=True)
            
            # Calculate total base pairs for each transition
            if os.path.exists(a_to_b_file) and os.stat(a_to_b_file).st_size > 0:
                a_to_b_df = pd.read_csv(a_to_b_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
                a_to_b_bp = int(a_to_b_df['end'].sum() - a_to_b_df['start'].sum())
                results[comp['name']]['A_to_B'] = (a_to_b_bp / genome_size) * 100
            
            if os.path.exists(b_to_a_file) and os.stat(b_to_a_file).st_size > 0:
                b_to_a_df = pd.read_csv(b_to_a_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
                b_to_a_bp = int(b_to_a_df['end'].sum() - b_to_a_df['start'].sum())
                results[comp['name']]['B_to_A'] = (b_to_a_bp / genome_size) * 100
            
            print(f"  A_to_B: {results[comp['name']]['A_to_B']:.2f}%")
            print(f"  B_to_A: {results[comp['name']]['B_to_A']:.2f}%")
    
    return results

def plot_chromatin_changes(results, output_file):
    """Create a bar chart showing chromatin state changes between conditions"""
    if not results:
        print("No results to plot")
        return
    
    # Convert results to a format suitable for plotting
    plot_data = []
    for comp_name, values in results.items():
        plot_data.append({
            'comparison': comp_name,
            'transition': 'A_to_B',
            'percentage': values['A_to_B']
        })
        plot_data.append({
            'comparison': comp_name,
            'transition': 'B_to_A',
            'percentage': values['B_to_A']
        })
    
    df = pd.DataFrame(plot_data)
    
    # Reshape data for easier plotting
    pivot_df = df.pivot(index='comparison', columns='transition', values='percentage')
    
    # Create figure
    fig, axes = plt.subplots(1, len(pivot_df), figsize=(10, 4), sharey=True)
    
    # Check if we have multiple subplots
    if len(pivot_df) == 1:
        axes = [axes]
    
    # Colors
    colors = {
        'A_to_B': '#8B0A50',  # Dark purple
        'B_to_A': '#F4A460'   # Sandy brown
    }
    
    # Create individual plots for each comparison
    for i, (idx, row) in enumerate(pivot_df.iterrows()):
        ax = axes[i]
        
        # Plot A_to_B (heterochromatin to euchromatin)
        ax.bar(0, row['A_to_B'], width=0.8, color=colors['A_to_B'], align='center')
        ax.text(0, row['A_to_B']/2, f"{row['A_to_B']:.2f}", ha='center', va='center', color='white', fontweight='bold')
        
        # Plot B_to_A (euchromatin to heterochromatin)
        ax.bar(1, row['B_to_A'], width=0.8, color=colors['B_to_A'], align='center')
        ax.text(1, row['B_to_A']/2, f"{row['B_to_A']:.2f}", ha='center', va='center', color='black', fontweight='bold')
        
        # Set title and remove axes
        ax.set_title(idx)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
    
    # Create a custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors['A_to_B'], label='A_to_B'),
        Patch(facecolor=colors['B_to_A'], label='B_to_A')
    ]
    
    # Add the legend below the plots
    fig.legend(
        handles=legend_elements,
        loc='lower center',
        ncol=2,
        bbox_to_anchor=(0.5, 0),
        frameon=False
    )
    
    plt.tight_layout(rect=[0, 0.1, 1, 1])  # Adjust layout to make room for the legend
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    
    # Also save data to TSV
    tsv_file = os.path.splitext(output_file)[0] + '.tsv'
    pivot_df.to_csv(tsv_file, sep='\t')
    print(f"Data saved to {tsv_file}")
    
    plt.close()

def main():
    """Main function"""
    args = parse_arguments()
    
    # Analyze chromatin state changes
    results = analyze_chromatin_changes(args.results_dir, args.genome_size)
    
    # Plot results
    if results:
        plot_chromatin_changes(results, args.output_file)
    else:
        print("No results to plot")

if __name__ == "__main__":
    main()