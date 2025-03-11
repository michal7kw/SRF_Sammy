#!/usr/bin/env python3
"""
Broad Region Chromatin State Changes Analysis

This script analyzes chromatin state changes between conditions over larger genomic windows,
allowing for detection of broader patterns of chromatin remodeling.

Usage:
    python broad_chromatin_changes.py --results-dir /path/to/results --window-size 150000 --output-file changes.png
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from scipy.stats import ttest_ind

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger()

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Broad Region Chromatin State Changes Analysis')
    parser.add_argument('--results-dir', required=True, help='Directory containing SAMMY-seq analysis results')
    parser.add_argument('--output-file', default='broad_chromatin_changes.png', help='Output file name')
    parser.add_argument('--genome-size', type=int, default=2800000000, help='Genome size in bp (default: 2.8Gb for mouse)')
    parser.add_argument('--window-size', type=int, default=150000, 
                        help='Size of genomic windows to analyze in bp (default: 150000)')
    parser.add_argument('--ratio-change-threshold', type=float, default=0.2, 
                        help='Minimum absolute change in log2 ratio to consider a region changed (default: 0.2)')
    parser.add_argument('--pvalue-threshold', type=float, default=0.05, 
                        help='P-value threshold for significant changes (default: 0.05)')
    parser.add_argument('--hetero-threshold', type=float, default=-0.1, 
                        help='Log2 ratio threshold for heterochromatin (default: -0.1)')
    parser.add_argument('--eu-threshold', type=float, default=0.1, 
                        help='Log2 ratio threshold for euchromatin (default: 0.1)')
    parser.add_argument('--debug', action='store_true', help='Print additional debugging information')
    return parser.parse_args()

def aggregate_into_windows(df, window_size, value_column='score'):
    """
    Aggregate bedgraph data into larger windows
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Bedgraph data with chrom, start, end, and score columns
    window_size : int
        Size of windows to aggregate into
    value_column : str
        Name of the column containing values to aggregate
        
    Returns:
    --------
    pandas.DataFrame
        Data aggregated into larger windows
    """
    # Create a new DataFrame to store aggregated data
    aggregated = []
    
    # Process each chromosome separately
    for chrom in df['chrom'].unique():
        chrom_df = df[df['chrom'] == chrom].copy()
        
        # Calculate the range of genomic positions for this chromosome
        chrom_min = chrom_df['start'].min()
        chrom_max = chrom_df['end'].max()
        
        # Create windows
        window_starts = np.arange(chrom_min, chrom_max, window_size)
        
        for win_start in window_starts:
            win_end = win_start + window_size
            
            # Find bins that overlap with this window
            overlap_df = chrom_df[(chrom_df['end'] > win_start) & (chrom_df['start'] < win_end)]
            
            if len(overlap_df) > 0:
                # Calculate weighted average based on the amount of overlap
                total_overlap = 0
                weighted_sum = 0
                
                for _, row in overlap_df.iterrows():
                    # Calculate overlap between bin and window
                    overlap_start = max(row['start'], win_start)
                    overlap_end = min(row['end'], win_end)
                    overlap_size = overlap_end - overlap_start
                    
                    # Add to weighted sum
                    weighted_sum += row[value_column] * overlap_size
                    total_overlap += overlap_size
                
                # Calculate weighted average if there is any overlap
                if total_overlap > 0:
                    avg_value = weighted_sum / total_overlap
                    
                    # Add to results
                    aggregated.append({
                        'chrom': chrom,
                        'start': win_start,
                        'end': win_end,
                        value_column: avg_value,
                        'n_bins': len(overlap_df),
                        'coverage': total_overlap / window_size
                    })
    
    # Convert to DataFrame
    return pd.DataFrame(aggregated)

def analyze_broad_chromatin_changes(results_dir, window_size, genome_size, ratio_change_threshold=0.2, 
                                   pvalue_threshold=0.05, hetero_threshold=-0.1, eu_threshold=0.1, 
                                   debug=False):
    """
    Analyze chromatin state changes between GFP and M2 conditions over broad genomic regions
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
    window_results = {}
    
    for comp in comparisons:
        # Check if required files exist
        if not os.path.exists(comp['gfp_ratio']) or not os.path.exists(comp['m2_ratio']):
            logger.warning(f"Skipping {comp['name']} - missing ratio files")
            continue
        
        logger.info(f"Processing {comp['name']} with {window_size}bp windows...")
        
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
        
        # Aggregate data into larger windows
        logger.info(f"  Aggregating GFP data into {window_size}bp windows...")
        gfp_windows = aggregate_into_windows(gfp_df, window_size)
        
        logger.info(f"  Aggregating M2 data into {window_size}bp windows...")
        m2_windows = aggregate_into_windows(m2_df, window_size)
        
        if debug:
            logger.info(f"  GFP windows: {len(gfp_windows)}")
            logger.info(f"  M2 windows: {len(m2_windows)}")
        
        # Merge the window dataframes on genomic coordinates
        logger.info("  Merging window data...")
        merged_windows = pd.merge(
            gfp_windows, 
            m2_windows, 
            on=['chrom', 'start', 'end'],
            suffixes=('_gfp', '_m2')
        )
        
        if debug:
            logger.info(f"  Merged windows: {len(merged_windows)}")
        
        # Filter for windows with good coverage
        min_coverage = 0.5  # Require at least 50% coverage
        merged_windows = merged_windows[
            (merged_windows['coverage_gfp'] >= min_coverage) & 
            (merged_windows['coverage_m2'] >= min_coverage)
        ]
        
        if debug:
            logger.info(f"  Windows with good coverage: {len(merged_windows)}")
        
        # Calculate the change in log2 ratio (M2 - GFP)
        merged_windows['ratio_change'] = merged_windows['score_m2'] - merged_windows['score_gfp']
        
        # Find windows where the ratio changed significantly towards euchromatin (more positive)
        a_to_b_windows = merged_windows[
            (merged_windows['ratio_change'] >= ratio_change_threshold)
        ]
        
        # Find windows where the ratio changed significantly towards heterochromatin (more negative)
        b_to_a_windows = merged_windows[
            (merged_windows['ratio_change'] <= -ratio_change_threshold)
        ]
        
        # Calculate total base pairs for each transition
        a_to_b_bp = a_to_b_windows['end'].sum() - a_to_b_windows['start'].sum()
        b_to_a_bp = b_to_a_windows['end'].sum() - b_to_a_windows['start'].sum()
        
        # Calculate percentages
        a_to_b_percent = (a_to_b_bp / genome_size) * 100
        b_to_a_percent = (b_to_a_bp / genome_size) * 100
        
        if debug:
            logger.info(f"  A→B windows: {len(a_to_b_windows)}, total bp: {a_to_b_bp}")
            logger.info(f"  B→A windows: {len(b_to_a_windows)}, total bp: {b_to_a_bp}")
            
            # Calculate statistics on changes
            mean_change = merged_windows['ratio_change'].mean()
            median_change = merged_windows['ratio_change'].median()
            logger.info(f"  Mean ratio change: {mean_change:.4f}")
            logger.info(f"  Median ratio change: {median_change:.4f}")
            
            # Statistical test - are the GFP and M2 distributions significantly different?
            t_stat, p_value = ttest_ind(
                merged_windows['score_gfp'], 
                merged_windows['score_m2'],
                equal_var=False
            )
            logger.info(f"  t-test: t={t_stat:.4f}, p={p_value:.6f}")
            
            # What percentage of windows show significant changes?
            total_windows = len(merged_windows)
            pct_a_to_b = (len(a_to_b_windows) / total_windows) * 100
            pct_b_to_a = (len(b_to_a_windows) / total_windows) * 100
            logger.info(f"  {pct_a_to_b:.2f}% of windows show A→B changes")
            logger.info(f"  {pct_b_to_a:.2f}% of windows show B→A changes")
        
        results[comp['name']] = {
            'A_to_B': a_to_b_percent,
            'B_to_A': b_to_a_percent,
            'p_value': p_value if 'p_value' in locals() else 1.0
        }
        
        # Store window data for plotting
        window_results[comp['name']] = merged_windows.copy()
        
        logger.info(f"  A_to_B: {a_to_b_percent:.2f}%")
        logger.info(f"  B_to_A: {b_to_a_percent:.2f}%")
    
    return results, window_results

def plot_broad_chromatin_changes(results, window_results, output_file, window_size, ratio_change_threshold=0.2):
    """Create visualizations of chromatin state changes across broad regions"""
    if not results:
        logger.error("No results to plot")
        return
    
    # Set style parameters
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14
    })
    
    # Create a figure with multiple subplots with adjusted height
    fig = plt.figure(figsize=(15, 15))  # Increased height from 12 to 15
    
    # Create GridSpec for better control over subplot layout
    gs = plt.GridSpec(3, len(results), height_ratios=[1, 1, 1.2])  # Increased ratio for bottom row
    gs.update(hspace=0.4, wspace=0.3)  # Increased spacing between plots
    
    # 1. Bar charts showing percentage of genome with significant changes
    ax1 = plt.subplot(gs[0, :])
    
    # Convert results to a dataframe for plotting
    comparisons = list(results.keys())
    a_to_b_values = [results[comp]['A_to_B'] for comp in comparisons]
    b_to_a_values = [results[comp]['B_to_A'] for comp in comparisons]
    
    # Bar positions
    bar_width = 0.35
    x = np.arange(len(comparisons))
    
    # Create bars with enhanced colors
    bars1 = ax1.bar(x - bar_width/2, a_to_b_values, bar_width, label='A→B', color='#8B0A50', alpha=0.8)
    bars2 = ax1.bar(x + bar_width/2, b_to_a_values, bar_width, label='B→A', color='#F4A460', alpha=0.8)
    
    # Add statistical significance indicators with better positioning
    for i, comp in enumerate(comparisons):
        p_value = results[comp].get('p_value', 1.0)
        if p_value < 0.001:
            significance = '***'
        elif p_value < 0.01:
            significance = '**'
        elif p_value < 0.05:
            significance = '*'
        else:
            significance = 'ns'
        
        max_height = max(a_to_b_values[i], b_to_a_values[i])
        ax1.text(i, max_height + max_height*0.05, significance, 
                 ha='center', va='bottom', fontweight='bold', fontsize=12)
    
    # Add labels and titles with improved formatting
    ax1.set_xlabel('Comparison', fontsize=11, labelpad=10)
    ax1.set_ylabel('Percentage of Genome', fontsize=11, labelpad=10)
    ax1.set_title(f'Chromatin Accessibility Changes in {window_size/1000:.0f}kb Windows\n(Changes ≥ {ratio_change_threshold} log2 ratio)',
                  pad=20, fontsize=14)
    ax1.set_xticks(x)
    ax1.set_xticklabels(comparisons, rotation=0)
    ax1.legend(frameon=True, fancybox=True, shadow=True)
    
    # Add value labels on bars with improved positioning
    def add_labels(bars):
        for bar in bars:
            height = bar.get_height()
            ax1.annotate(f'{height:.2f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom',
                        fontsize=10,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
    
    add_labels(bars1)
    add_labels(bars2)
    
    # 2. Histograms for each comparison
    for i, comp_name in enumerate(comparisons):
        ax2 = plt.subplot(gs[1, i])
        
        if comp_name in window_results:
            data = window_results[comp_name]
            
            # Plot histogram with enhanced appearance
            n, bins, patches = ax2.hist(data['ratio_change'], bins=40, density=True, alpha=0.75)
            
            # Color the histogram with improved colors
            for j, patch in enumerate(patches):
                if bins[j] > 0:
                    patch.set_facecolor('#8B0A50')  # A to B (purple)
                else:
                    patch.set_facecolor('#F4A460')  # B to A (orange)
            
            # Add vertical lines with better styling
            ax2.axvline(x=ratio_change_threshold, color='#8B0A50', linestyle='--', alpha=0.7, linewidth=1.5)
            ax2.axvline(x=-ratio_change_threshold, color='#F4A460', linestyle='--', alpha=0.7, linewidth=1.5)
            ax2.axvline(x=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
            
            # Add labels with improved formatting
            ax2.set_xlabel('Log2 Ratio Change (M2 - GFP)', labelpad=8)
            ax2.set_ylabel('Density', labelpad=8)
            ax2.set_title(f'{comp_name} - Ratio Changes', pad=10)
            ax2.grid(True, alpha=0.3)
    
    # 3. Scatter plots comparing GFP vs M2 values
    for i, comp_name in enumerate(comparisons):
        ax3 = plt.subplot(gs[2, i])
        
        if comp_name in window_results:
            data = window_results[comp_name]
            
            # Create scatter plot with enhanced appearance
            sc = ax3.scatter(data['score_gfp'], data['score_m2'], 
                           c=data['ratio_change'], cmap='coolwarm', 
                           alpha=0.6, s=5)
            
            # Add colorbar with better formatting
            cbar = plt.colorbar(sc, ax=ax3)
            cbar.set_label('Log2 Ratio Change', rotation=270, labelpad=15)
            
            # Add diagonal line with improved style
            min_val = min(data['score_gfp'].min(), data['score_m2'].min())
            max_val = max(data['score_gfp'].max(), data['score_m2'].max())
            ax3.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.7, linewidth=1.5)
            
            # Add thresholds with better styling
            ax3.axhline(y=0.1, color='green', linestyle=':', alpha=0.5, linewidth=1)
            ax3.axhline(y=-0.1, color='blue', linestyle=':', alpha=0.5, linewidth=1)
            ax3.axvline(x=0.1, color='green', linestyle=':', alpha=0.5, linewidth=1)
            ax3.axvline(x=-0.1, color='blue', linestyle=':', alpha=0.5, linewidth=1)
            
            # Add labels with improved formatting
            ax3.set_xlabel('GFP Log2 Ratio', labelpad=8)
            ax3.set_ylabel('M2 Log2 Ratio', labelpad=8)
            ax3.set_title(f'{comp_name} - Value Comparison', pad=10)
            ax3.grid(True, alpha=0.3)
    
    # Save figure with high quality
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    logger.info(f"Plot saved to {output_file}")
    
    # Generate the side-by-side bar chart
    create_side_by_side_chart(results, os.path.splitext(output_file)[0] + '_simple.png')
    
    # Save data to TSV
    tsv_file = os.path.splitext(output_file)[0] + '.tsv'
    df = pd.DataFrame(results).T
    df.to_csv(tsv_file, sep='\t')
    logger.info(f"Data saved to {tsv_file}")
    
    plt.close()

def create_side_by_side_chart(results, output_file):
    """Create a simple side-by-side chart with improved formatting"""
    if not results:
        return
    
    # Set style parameters
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14
    })
    
    # Create figure with better dimensions
    fig, axes = plt.subplots(1, len(results), figsize=(12, 6), sharey=True)
    
    # Check if we have multiple subplots
    if len(results) == 1:
        axes = [axes]
    
    # Enhanced colors with better contrast
    colors = {
        'A_to_B': '#8B0A50',  # Dark purple for A→B
        'B_to_A': '#F4A460'   # Sandy brown for B→A
    }
    
    # Create individual plots for each comparison
    max_height = 0
    for i, (comp_name, values) in enumerate(results.items()):
        ax = axes[i]
        
        # Plot the bars side by side with enhanced styling
        a_to_b_height = values['A_to_B']
        b_to_a_height = values['B_to_A']
        max_height = max(max_height, a_to_b_height, b_to_a_height)
        
        # Plot A_to_B with improved appearance
        bar1 = ax.bar(0, a_to_b_height, width=0.8, color=colors['A_to_B'], 
                     align='center', alpha=0.8)
        ax.text(0, a_to_b_height/2, f"{a_to_b_height:.2f}%", 
                ha='center', va='center', color='white', 
                fontweight='bold', fontsize=11)
        
        # Plot B_to_A with improved appearance
        bar2 = ax.bar(1, b_to_a_height, width=0.8, color=colors['B_to_A'], 
                     align='center', alpha=0.8)
        ax.text(1, b_to_a_height/2, f"{b_to_a_height:.2f}%", 
                ha='center', va='center', color='black', 
                fontweight='bold', fontsize=11)
        
        # Set title and clean up axes with better formatting
        ax.set_title(comp_name, pad=15, fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Remove spines with slight padding
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        # Add subtle grid
        ax.grid(False)
        
        # Add statistical significance if available
        p_value = values.get('p_value', 1.0)
        if p_value < 0.001:
            significance = '***'
        elif p_value < 0.01:
            significance = '**'
        elif p_value < 0.05:
            significance = '*'
        else:
            significance = 'ns'
        
        ax.text(0.5, max(a_to_b_height, b_to_a_height) * 1.1, 
                significance, ha='center', va='bottom', 
                fontweight='bold', fontsize=12)
    
    # Create a custom legend with better styling
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors['A_to_B'], label='A→B', alpha=0.8),
        Patch(facecolor=colors['B_to_A'], label='B→A', alpha=0.8)
    ]
    
    # Add the legend below the plots with improved positioning
    fig.legend(
        handles=legend_elements,
        loc='lower center',
        ncol=2,
        bbox_to_anchor=(0.5, -0.1),
        frameon=True,
        fancybox=True,
        shadow=True,
        fontsize=10
    )
    
    # Add main title
    plt.suptitle('Chromatin Accessibility Changes', 
                 y=1.02, fontsize=14, fontweight='bold')
    
    # Adjust layout with better spacing
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    
    # Save figure with high quality
    plt.savefig(output_file, dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    logger.info(f"Simple chart saved to {output_file}")
    
    plt.close()

def main():
    """Main function"""
    args = parse_arguments()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    # Analyze broad region chromatin changes
    results, window_results = analyze_broad_chromatin_changes(
        args.results_dir, 
        args.window_size,
        args.genome_size, 
        args.ratio_change_threshold,
        args.pvalue_threshold,
        args.hetero_threshold,
        args.eu_threshold,
        args.debug
    )
    
    # Plot results
    if results:
        plot_broad_chromatin_changes(
            results, 
            window_results, 
            args.output_file, 
            args.window_size,
            args.ratio_change_threshold
        )
    else:
        logger.error("No results to plot")

if __name__ == "__main__":
    main()