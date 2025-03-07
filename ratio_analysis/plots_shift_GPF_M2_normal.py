#!/usr/bin/env python3
"""
Improved Visualization of Chromatin State Changes Between Conditions

This script analyzes chromatin state bed files and bedgraph files to determine how 
chromatin structure changes between GFP and M2 conditions, with multiple analysis methods
and flexible overlap criteria.

Usage:
    python improved_chromatin_changes.py --results-dir /path/to/results --output-file changes.png
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
    parser = argparse.ArgumentParser(description='Improved Chromatin State Changes Between Conditions')
    parser.add_argument('--results-dir', required=True, help='Directory containing SAMMY-seq analysis results')
    parser.add_argument('--output-file', default='chromatin_state_changes.png', help='Output file name')
    parser.add_argument('--genome-size', type=int, default=2800000000, help='Genome size in bp (default: 2.8Gb for mouse)')
    parser.add_argument('--method', choices=['intersect', 'ratio_comparison', 'combined'], default='combined', 
                        help='Method to detect chromatin state changes (default: combined)')
    parser.add_argument('--min-overlap', type=float, default=0.1, 
                        help='Minimum overlap required as fraction of region (default: 0.1)')
    parser.add_argument('--debug', action='store_true', help='Print additional debugging information')
    return parser.parse_args()

def analyze_chromatin_changes_by_intersection(results_dir, genome_size, min_overlap=0.1, debug=False):
    """
    Analyze changes in chromatin state between GFP and M2 conditions
    using bedtools intersect with flexible overlap criteria
    """
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
                logger.warning(f"Skipping {comp['name']} - missing required files")
                continue
            
            # Log file sizes for debugging
            if debug:
                for key in ['gfp_hetero', 'gfp_eu', 'm2_hetero', 'm2_eu']:
                    file_size = os.path.getsize(comp[key])
                    count = sum(1 for _ in open(comp[key]))
                    logger.info(f"  {key}: {file_size} bytes, {count} lines")
            
            # Initialize results for this comparison
            results[comp['name']] = {'A_to_B': 0, 'B_to_A': 0}
            
            logger.info(f"Processing {comp['name']}...")
            
            # 1. Identify A_to_B: Regions that were heterochromatin (A) in GFP and became euchromatin (B) in M2
            a_to_b_file = os.path.join(tmp_dir, f"{comp['name'].replace(' ', '_')}_A_to_B.bed")
            cmd = f"bedtools intersect -a {comp['gfp_hetero']} -b {comp['m2_eu']} -wa -f {min_overlap} > {a_to_b_file}"
            subprocess.run(cmd, shell=True, check=True)
            
            # 2. Identify B_to_A: Regions that were euchromatin (B) in GFP and became heterochromatin (A) in M2
            b_to_a_file = os.path.join(tmp_dir, f"{comp['name'].replace(' ', '_')}_B_to_A.bed")
            cmd = f"bedtools intersect -a {comp['gfp_eu']} -b {comp['m2_hetero']} -wa -f {min_overlap} > {b_to_a_file}"
            subprocess.run(cmd, shell=True, check=True)
            
            # Calculate total base pairs for each transition
            if os.path.exists(a_to_b_file) and os.stat(a_to_b_file).st_size > 0:
                a_to_b_df = pd.read_csv(a_to_b_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
                a_to_b_bp = int(a_to_b_df['end'].sum() - a_to_b_df['start'].sum())
                results[comp['name']]['A_to_B'] = (a_to_b_bp / genome_size) * 100
                
                if debug:
                    logger.info(f"  A_to_B regions: {len(a_to_b_df)}, total bp: {a_to_b_bp}")
            else:
                logger.warning(f"  No A_to_B regions found for {comp['name']}")
            
            if os.path.exists(b_to_a_file) and os.stat(b_to_a_file).st_size > 0:
                b_to_a_df = pd.read_csv(b_to_a_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
                b_to_a_bp = int(b_to_a_df['end'].sum() - b_to_a_df['start'].sum())
                results[comp['name']]['B_to_A'] = (b_to_a_bp / genome_size) * 100
                
                if debug:
                    logger.info(f"  B_to_A regions: {len(b_to_a_df)}, total bp: {b_to_a_bp}")
            else:
                logger.warning(f"  No B_to_A regions found for {comp['name']}")
            
            logger.info(f"  A_to_B: {results[comp['name']]['A_to_B']:.2f}%")
            logger.info(f"  B_to_A: {results[comp['name']]['B_to_A']:.2f}%")
    
    return results

def analyze_chromatin_changes_by_ratio_comparison(results_dir, genome_size, debug=False):
    """
    Alternative method: Compare the ratio bedgraph files directly to identify 
    regions where the ratio changes from above threshold to below threshold or vice versa
    """
    bedgraph_dir = os.path.join(results_dir, "bedgraph")
    
    # Define comparisons and thresholds
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
    
    eu_threshold = 0.1
    hetero_threshold = -0.1
    
    results = {}
    
    for comp in comparisons:
        # Check if required files exist
        if not os.path.exists(comp['gfp_ratio']) or not os.path.exists(comp['m2_ratio']):
            logger.warning(f"Skipping {comp['name']} - missing ratio files")
            continue
        
        logger.info(f"Processing {comp['name']} using ratio comparison method...")
        
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
        # We'll use an inner merge to only keep regions that appear in both datasets
        merged_df = pd.merge(
            gfp_df, 
            m2_df, 
            on=['chrom', 'start', 'end'],
            suffixes=('_gfp', '_m2')
        )
        
        if debug:
            logger.info(f"  Merged dataframe: {len(merged_df)} regions")
        
        # Identify regions with state changes
        a_to_b_regions = merged_df[
            (merged_df['score_gfp'] <= hetero_threshold) & # Heterochromatin in GFP
            (merged_df['score_m2'] >= eu_threshold)        # Euchromatin in M2
        ]
        
        b_to_a_regions = merged_df[
            (merged_df['score_gfp'] >= eu_threshold) &     # Euchromatin in GFP
            (merged_df['score_m2'] <= hetero_threshold)    # Heterochromatin in M2
        ]
        
        # Calculate total base pairs for each transition
        a_to_b_bp = int(a_to_b_regions['end'].sum() - a_to_b_regions['start'].sum())
        b_to_a_bp = int(b_to_a_regions['end'].sum() - b_to_a_regions['start'].sum())
        
        # Calculate percentages
        a_to_b_percent = (a_to_b_bp / genome_size) * 100
        b_to_a_percent = (b_to_a_bp / genome_size) * 100
        
        if debug:
            logger.info(f"  A_to_B regions: {len(a_to_b_regions)}, total bp: {a_to_b_bp}")
            logger.info(f"  B_to_A regions: {len(b_to_a_regions)}, total bp: {b_to_a_bp}")
        
        results[comp['name']] = {
            'A_to_B': a_to_b_percent,
            'B_to_A': b_to_a_percent
        }
        
        logger.info(f"  A_to_B: {a_to_b_percent:.2f}%")
        logger.info(f"  B_to_A: {b_to_a_percent:.2f}%")
    
    return results

def plot_chromatin_changes(results, output_file):
    """Create a bar chart showing chromatin state changes between conditions"""
    if not results:
        logger.error("No results to plot")
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
    fig, axes = plt.subplots(1, len(pivot_df), figsize=(12, 4), sharey=True)
    
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
        ax.text(0, row['A_to_B']/2, f"{row['A_to_B']:.2f}", ha='center', va='center', 
                color='white', fontweight='bold')
        
        # Plot B_to_A (euchromatin to heterochromatin)
        ax.bar(1, row['B_to_A'], width=0.8, color=colors['B_to_A'], align='center')
        ax.text(1, row['B_to_A']/2, f"{row['B_to_A']:.2f}", ha='center', va='center', 
                color='black', fontweight='bold')
        
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
    logger.info(f"Plot saved to {output_file}")
    
    # Also save data to TSV
    tsv_file = os.path.splitext(output_file)[0] + '.tsv'
    pivot_df.to_csv(tsv_file, sep='\t')
    logger.info(f"Data saved to {tsv_file}")
    
    plt.close()

def main():
    """Main function"""
    args = parse_arguments()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    # Use the selected method for analysis
    if args.method == 'intersect':
        logger.info("Using bedtools intersect method for analysis")
        results = analyze_chromatin_changes_by_intersection(
            args.results_dir, args.genome_size, args.min_overlap, args.debug
        )
    elif args.method == 'ratio_comparison':
        logger.info("Using ratio comparison method for analysis")
        results = analyze_chromatin_changes_by_ratio_comparison(
            args.results_dir, args.genome_size, args.debug
        )
    elif args.method == 'combined':
        logger.info("Using both methods and combining results")
        results_intersect = analyze_chromatin_changes_by_intersection(
            args.results_dir, args.genome_size, args.min_overlap, args.debug
        )
        results_ratio = analyze_chromatin_changes_by_ratio_comparison(
            args.results_dir, args.genome_size, args.debug
        )
        
        # Combine results, taking the maximum value from either method
        results = {}
        for comp_name in set(list(results_intersect.keys()) + list(results_ratio.keys())):
            results[comp_name] = {
                'A_to_B': max(
                    results_intersect.get(comp_name, {'A_to_B': 0})['A_to_B'],
                    results_ratio.get(comp_name, {'A_to_B': 0})['A_to_B']
                ),
                'B_to_A': max(
                    results_intersect.get(comp_name, {'B_to_A': 0})['B_to_A'],
                    results_ratio.get(comp_name, {'B_to_A': 0})['B_to_A']
                )
            }
    
    # Plot results
    if results:
        plot_chromatin_changes(results, args.output_file)
    else:
        logger.error("No results to plot")

if __name__ == "__main__":
    main()