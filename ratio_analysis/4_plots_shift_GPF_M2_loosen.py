#!/usr/bin/env python3
"""
Subtle Chromatin State Changes Analysis

This script analyzes the magnitude of changes in chromatin accessibility between
GFP and M2 conditions, detecting even subtle shifts that don't necessarily cross
the strict heterochromatin/euchromatin thresholds.

Usage:
    python subtle_chromatin_changes.py --results-dir /path/to/results --output-file changes.png
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger()

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Subtle Chromatin State Changes Analysis')
    parser.add_argument('--results-dir', required=True, help='Directory containing SAMMY-seq analysis results')
    parser.add_argument('--output-file', default='subtle_chromatin_changes.png', help='Output file name')
    parser.add_argument('--genome-size', type=int, default=2800000000, help='Genome size in bp (default: 2.8Gb for mouse)')
    parser.add_argument('--ratio-change-threshold', type=float, default=0.2, 
                        help='Minimum absolute change in log2 ratio to consider a region changed (default: 0.2)')
    parser.add_argument('--hetero-threshold', type=float, default=-0.1, 
                        help='Log2 ratio threshold for heterochromatin (default: -0.1)')
    parser.add_argument('--eu-threshold', type=float, default=0.1, 
                        help='Log2 ratio threshold for euchromatin (default: 0.1)')
    parser.add_argument('--debug', action='store_true', help='Print additional debugging information')
    return parser.parse_args()

def analyze_subtle_chromatin_changes(results_dir, genome_size, ratio_change_threshold=0.2, 
                                  hetero_threshold=-0.1, eu_threshold=0.1, debug=False):
    """
    Analyze subtle changes in chromatin state between GFP and M2 conditions
    looking for significant shifts in log2 ratio even if they don't cross thresholds
    """
    bedgraph_dir = os.path.join(results_dir, "bedgraph")
    
    # Define comparisons
    comparisons = [
        {
            'name': 'Neu_GFP vs Neu_M2',
            'gfp_ratio': os.path.join(bedgraph_dir, 'Neu_GFP_S2S_vs_S3_ratio.bedgraph'),
            'm2_ratio': os.path.join(bedgraph_dir, 'Neu_M2_S2S_vs_S3_ratio.bedgraph'),
        },
        {
            'name': 'NSC_GFP vs NSC_M2',
            'gfp_ratio': os.path.join(bedgraph_dir, 'NSC_GFP_S2S_vs_S3_ratio.bedgraph'),
            'm2_ratio': os.path.join(bedgraph_dir, 'NSC_M2_S2S_vs_S3_ratio.bedgraph'),
        }
    ]
    
    results = {}
    histograms = {}
    
    for comp in comparisons:
        # Check if required files exist
        if not os.path.exists(comp['gfp_ratio']) or not os.path.exists(comp['m2_ratio']):
            logger.warning(f"Skipping {comp['name']} - missing ratio files")
            continue
        
        logger.info(f"Processing {comp['name']} using subtle changes method...")
        
        # Read ratio files
        try:
            gfp_df = pd.read_csv(comp['gfp_ratio'], sep='\t', header=None, 
                                 names=['chrom', 'start', 'end', 'score'])
            m2_df = pd.read_csv(comp['m2_ratio'], sep='\t', header=None, 
                               names=['chrom', 'start', 'end', 'score'])
            
            if debug:
                logger.info(f"  GFP ratio file: {len(gfp_df)} regions")
                logger.info(f"  M2 ratio file: {len(m2_df)} regions")
                logger.info(f"  GFP ratio range: {gfp_df['score'].min()} to {gfp_df['score'].max()}")
                logger.info(f"  M2 ratio range: {m2_df['score'].min()} to {m2_df['score'].max()}")
        except Exception as e:
            logger.error(f"Error reading ratio files for {comp['name']}: {e}")
            continue
        
        # Merge the dataframes on genomic coordinates
        merged_df = pd.merge(
            gfp_df, 
            m2_df, 
            on=['chrom', 'start', 'end'],
            suffixes=('_gfp', '_m2')
        )
        
        if debug:
            logger.info(f"  Merged dataframe: {len(merged_df)} regions")
        
        # Calculate the change in log2 ratio (M2 - GFP)
        merged_df['ratio_change'] = merged_df['score_m2'] - merged_df['score_gfp']
        
        # Save histogram data for later visualization
        histograms[comp['name']] = merged_df['ratio_change'].copy()
        
        # Find regions where the ratio changed significantly towards euchromatin (more positive)
        a_to_b_regions = merged_df[
            (merged_df['ratio_change'] >= ratio_change_threshold)
        ]
        
        # Find regions where the ratio changed significantly towards heterochromatin (more negative)
        b_to_a_regions = merged_df[
            (merged_df['ratio_change'] <= -ratio_change_threshold)
        ]
        
        # Add a filter for strict A_to_B and B_to_A transitions
        strict_a_to_b = merged_df[
            (merged_df['score_gfp'] <= hetero_threshold) & 
            (merged_df['score_m2'] >= eu_threshold)
        ]
        
        strict_b_to_a = merged_df[
            (merged_df['score_gfp'] >= eu_threshold) & 
            (merged_df['score_m2'] <= hetero_threshold)
        ]
        
        # Calculate total base pairs for each transition
        a_to_b_bp = int(a_to_b_regions['end'].sum() - a_to_b_regions['start'].sum())
        b_to_a_bp = int(b_to_a_regions['end'].sum() - b_to_a_regions['start'].sum())
        
        # Calculate percentages
        a_to_b_percent = (a_to_b_bp / genome_size) * 100
        b_to_a_percent = (b_to_a_bp / genome_size) * 100
        
        if debug:
            logger.info(f"  Significant A→B regions: {len(a_to_b_regions)}, total bp: {a_to_b_bp}")
            logger.info(f"  Significant B→A regions: {len(b_to_a_regions)}, total bp: {b_to_a_bp}")
            logger.info(f"  Strict A→B transitions: {len(strict_a_to_b)}")
            logger.info(f"  Strict B→A transitions: {len(strict_b_to_a)}")
            
            # Calculate mean change
            mean_change = merged_df['ratio_change'].mean()
            median_change = merged_df['ratio_change'].median()
            logger.info(f"  Mean ratio change: {mean_change:.4f}")
            logger.info(f"  Median ratio change: {median_change:.4f}")
        
        results[comp['name']] = {
            'A_to_B': a_to_b_percent,
            'B_to_A': b_to_a_percent
        }
        
        logger.info(f"  A_to_B: {a_to_b_percent:.2f}%")
        logger.info(f"  B_to_A: {b_to_a_percent:.2f}%")
    
    return results, histograms

def plot_chromatin_changes(results, histograms, output_file, ratio_change_threshold=0.2):
    """Create visualizations of chromatin state changes"""
    if not results:
        logger.error("No results to plot")
        return
    
    # Create a figure with multiple subplots
    fig = plt.figure(figsize=(15, 10))
    
    # 1. Bar charts showing percentage of genome with significant changes
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2)
    
    # Convert results to a dataframe for plotting
    comparisons = list(results.keys())
    a_to_b_values = [results[comp]['A_to_B'] for comp in comparisons]
    b_to_a_values = [results[comp]['B_to_A'] for comp in comparisons]
    
    # Bar positions
    bar_width = 0.35
    x = np.arange(len(comparisons))
    
    # Create bars
    bars1 = ax1.bar(x - bar_width/2, a_to_b_values, bar_width, label='A→B', color='#8B0A50')
    bars2 = ax1.bar(x + bar_width/2, b_to_a_values, bar_width, label='B→A', color='#F4A460')
    
    # Add labels and titles
    ax1.set_xlabel('Comparison')
    ax1.set_ylabel('Percentage of Genome')
    ax1.set_title(f'Chromatin Accessibility Changes\n(Changes ≥ {ratio_change_threshold} log2 ratio)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(comparisons)
    ax1.legend()
    
    # Add value labels on bars
    def add_labels(bars):
        for bar in bars:
            height = bar.get_height()
            ax1.annotate(f'{height:.2f}%',
                         xy=(bar.get_x() + bar.get_width() / 2, height),
                         xytext=(0, 3),
                         textcoords="offset points",
                         ha='center', va='bottom')
    
    add_labels(bars1)
    add_labels(bars2)
    
    # 2. Individual visualization for each comparison
    for i, comp_name in enumerate(comparisons):
        # Calculate position
        ax = plt.subplot2grid((2, len(comparisons)), (1, i))
        
        # Plot histogram of log2 ratio changes
        if comp_name in histograms:
            hist_data = histograms[comp_name]
            n, bins, patches = ax.hist(hist_data, bins=50, density=True, alpha=0.75)
            
            # Color the histogram based on whether changes are towards A or B
            for j, patch in enumerate(patches):
                if bins[j] > 0:
                    patch.set_facecolor('#8B0A50')  # A to B (purple)
                else:
                    patch.set_facecolor('#F4A460')  # B to A (orange)
            
            # Add vertical lines for thresholds
            ax.axvline(x=ratio_change_threshold, color='#8B0A50', linestyle='--', alpha=0.7)
            ax.axvline(x=-ratio_change_threshold, color='#F4A460', linestyle='--', alpha=0.7)
            ax.axvline(x=0, color='black', linestyle='-', alpha=0.5)
            
            # Calculate and display the percentage of regions with significant changes
            a_to_b_count = sum(hist_data >= ratio_change_threshold)
            b_to_a_count = sum(hist_data <= -ratio_change_threshold)
            total_count = len(hist_data)
            
            a_to_b_pct = (a_to_b_count / total_count) * 100
            b_to_a_pct = (b_to_a_count / total_count) * 100
            
            # Add text annotations
            ax.text(0.05, 0.95, f"A→B: {a_to_b_pct:.1f}% of regions", 
                    transform=ax.transAxes, verticalalignment='top', color='#8B0A50')
            ax.text(0.05, 0.87, f"B→A: {b_to_a_pct:.1f}% of regions", 
                    transform=ax.transAxes, verticalalignment='top', color='#F4A460')
            
            # Add labels and title
            ax.set_xlabel('Log2 Ratio Change (M2 - GFP)')
            ax.set_ylabel('Density')
            ax.set_title(comp_name)
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logger.info(f"Plot saved to {output_file}")
    
    # Generate the side-by-side bar chart like the example
    create_side_by_side_chart(results, os.path.splitext(output_file)[0] + '_simple.png')
    
    # Save data to TSV
    tsv_file = os.path.splitext(output_file)[0] + '.tsv'
    df = pd.DataFrame(results).T
    df.to_csv(tsv_file, sep='\t')
    logger.info(f"Data saved to {tsv_file}")
    
    plt.close()

def create_side_by_side_chart(results, output_file):
    """Create a simple side-by-side chart like in the example"""
    if not results:
        return
    
    fig, axes = plt.subplots(1, len(results), figsize=(10, 4), sharey=True)
    
    # Check if we have multiple subplots
    if len(results) == 1:
        axes = [axes]
    
    # Colors
    colors = {
        'A_to_B': '#8B0A50',  # Dark purple for A→B
        'B_to_A': '#F4A460'   # Sandy brown for B→A
    }
    
    # Create individual plots for each comparison
    for i, (comp_name, values) in enumerate(results.items()):
        ax = axes[i]
        
        # Plot the bars side by side
        a_to_b_height = values['A_to_B']
        b_to_a_height = values['B_to_A']
        
        # Plot A_to_B (heterochromatin to euchromatin)
        ax.bar(0, a_to_b_height, width=0.8, color=colors['A_to_B'], align='center')
        ax.text(0, a_to_b_height/2, f"{a_to_b_height:.2f}", ha='center', va='center', 
                color='white', fontweight='bold')
        
        # Plot B_to_A (euchromatin to heterochromatin)
        ax.bar(1, b_to_a_height, width=0.8, color=colors['B_to_A'], align='center')
        ax.text(1, b_to_a_height/2, f"{b_to_a_height:.2f}", ha='center', va='center', 
                color='black', fontweight='bold')
        
        # Set title and remove axes
        ax.set_title(comp_name)
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
    logger.info(f"Simple chart saved to {output_file}")
    
    plt.close()

def main():
    """Main function"""
    args = parse_arguments()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    # Analyze subtle chromatin changes
    results, histograms = analyze_subtle_chromatin_changes(
        args.results_dir, 
        args.genome_size, 
        args.ratio_change_threshold,
        args.hetero_threshold,
        args.eu_threshold,
        args.debug
    )
    
    # Plot results
    if results:
        plot_chromatin_changes(results, histograms, args.output_file, args.ratio_change_threshold)
    else:
        logger.error("No results to plot")

if __name__ == "__main__":
    main()