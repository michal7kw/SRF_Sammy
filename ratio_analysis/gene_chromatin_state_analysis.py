#!/usr/bin/env python3
"""
Analysis of chromatin states for specific genes in NSC and Neu samples (GFP condition)
This script visualizes the chromatin states of genes from a provided list
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import logging
import subprocess
import tempfile

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("gene_chromatin_state_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('Gene-Chromatin-Analysis')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Gene Chromatin State Analysis')
    parser.add_argument('--gene-list', required=True, help='Path to gene list file')
    parser.add_argument('--results-dir', required=True, help='Directory containing chromatin state results')
    parser.add_argument('--output-dir', required=True, help='Output directory for analysis results')
    parser.add_argument('--genome-gtf', default='/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/ratio_analysis/reference/Mus_musculus.GRCm38.102.gtf', 
                        help='Path to genome GTF file for gene coordinates')
    parser.add_argument('--promoter-size', type=int, default=2000, 
                        help='Size of promoter region in base pairs upstream of transcription start site (default: 2000)')
    return parser.parse_args()

def load_gene_list(file_path):
    """Load gene list from file"""
    try:
        with open(file_path, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(genes)} genes from {file_path}")
        return genes
    except Exception as e:
        logger.error(f"Error loading gene list: {e}")
        sys.exit(1)

def load_chromatin_states(results_dir, condition, state_type):
    """Load chromatin state data for a specific condition and state type"""
    state_file = os.path.join(results_dir, "chromatin_states", f"{condition}_{state_type}.bed")
    
    if not os.path.exists(state_file):
        logger.error(f"Chromatin state file not found: {state_file}")
        return pd.DataFrame()
    
    try:
        # The 4th column is the score, not a state name
        states = pd.read_csv(state_file, sep='\t',
                            names=['chromosome', 'start', 'end', 'score'])
        # Add a state type column
        states['state_type'] = state_type
        logger.info(f"Loaded {len(states)} {state_type} regions for {condition}")
        return states
    except Exception as e:
        logger.error(f"Error loading chromatin states: {e}")
        return pd.DataFrame()

def get_gene_promoter_coordinates(genes, gtf_file, promoter_size=2000):
    """Extract promoter coordinates for genes from GTF file
    
    Promoter is defined as a region upstream of the transcription start site (TSS)
    Default promoter size is 2kb upstream of TSS
    """
    logger.info(f"Extracting promoter coordinates for {len(genes)} genes from GTF file")
    
    try:
        # Use grep to extract gene entries from GTF file
        gene_promoters = {}
        for gene in genes:
            # Extract gene entries for this gene using the gene_name attribute
            cmd = f"grep 'gene_name \"{gene}\"' {gtf_file} | grep -w 'gene' | head -n 1"
            try:
                result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                if result.stdout.strip():
                    # Parse the GTF line
                    fields = result.stdout.strip().split('\t')
                    if len(fields) >= 8:  # Need at least 8 fields to get strand information
                        chromosome = fields[0]
                        feature_type = fields[2]
                        start = int(fields[3])
                        end = int(fields[4])
                        strand = fields[6]  # '+' or '-'
                        
                        if feature_type == 'gene':
                            # Define promoter region based on strand
                            if strand == '+':
                                # For genes on + strand, promoter is upstream of start
                                promoter_start = max(1, start - promoter_size)  # Ensure start is not negative
                                promoter_end = start
                            else:  # strand == '-'
                                # For genes on - strand, promoter is upstream of end
                                promoter_start = end
                                promoter_end = end + promoter_size
                            
                            gene_promoters[gene] = {
                                'chromosome': chromosome,
                                'start': promoter_start,
                                'end': promoter_end,
                                'strand': strand
                            }
                            logger.debug(f"Found promoter coordinates for gene {gene}: {chromosome}:{promoter_start}-{promoter_end} ({strand})")
            except subprocess.CalledProcessError:
                # Try an alternative approach - some gene names might have special characters
                # that need to be escaped in the grep pattern
                try:
                    # Try with a more flexible pattern
                    cmd = f"grep 'gene_name' {gtf_file} | grep '{gene}' | grep -w 'gene' | head -n 1"
                    result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                    if result.stdout.strip():
                        # Parse the GTF line
                        fields = result.stdout.strip().split('\t')
                        if len(fields) >= 8:
                            chromosome = fields[0]
                            feature_type = fields[2]
                            start = int(fields[3])
                            end = int(fields[4])
                            strand = fields[6]  # '+' or '-'
                            
                            if feature_type == 'gene':
                                # Define promoter region based on strand
                                if strand == '+':
                                    promoter_start = max(1, start - promoter_size)
                                    promoter_end = start
                                else:  # strand == '-'
                                    promoter_start = end
                                    promoter_end = end + promoter_size
                                
                                gene_promoters[gene] = {
                                    'chromosome': chromosome,
                                    'start': promoter_start,
                                    'end': promoter_end,
                                    'strand': strand
                                }
                                logger.debug(f"Found promoter coordinates for gene {gene} (alt method): {chromosome}:{promoter_start}-{promoter_end} ({strand})")
                except subprocess.CalledProcessError:
                    logger.warning(f"Could not find coordinates for gene {gene}")
        
        logger.info(f"Found promoter coordinates for {len(gene_promoters)} out of {len(genes)} genes")
        return gene_promoters
    
    except Exception as e:
        logger.error(f"Error extracting gene promoter coordinates: {e}")
        return {}

def assign_chromatin_states(gene_coords, eu_states, hetero_states):
    """Assign chromatin states to genes based on overlap with chromatin regions"""
    gene_states = {}
    
    for gene, coords in gene_coords.items():
        chrom = coords['chromosome']
        gene_start = coords['start']
        gene_end = coords['end']
        
        # Find overlapping euchromatin regions
        eu_overlaps = eu_states[
            (eu_states['chromosome'] == chrom) & 
            (eu_states['end'] >= gene_start) & 
            (eu_states['start'] <= gene_end)
        ]
        
        # Find overlapping heterochromatin regions
        hetero_overlaps = hetero_states[
            (hetero_states['chromosome'] == chrom) & 
            (hetero_states['end'] >= gene_start) & 
            (hetero_states['start'] <= gene_end)
        ]
        
        # Determine the predominant state
        if not eu_overlaps.empty and not hetero_overlaps.empty:
            # Calculate the total overlap length for each state
            eu_overlap_len = 0
            for _, region in eu_overlaps.iterrows():
                overlap_start = max(region['start'], gene_start)
                overlap_end = min(region['end'], gene_end)
                eu_overlap_len += (overlap_end - overlap_start)
            
            hetero_overlap_len = 0
            for _, region in hetero_overlaps.iterrows():
                overlap_start = max(region['start'], gene_start)
                overlap_end = min(region['end'], gene_end)
                hetero_overlap_len += (overlap_end - overlap_start)
            
            # Assign state based on which has more overlap
            if eu_overlap_len > hetero_overlap_len:
                gene_states[gene] = "Euchromatin"
            else:
                gene_states[gene] = "Heterochromatin"
        elif not eu_overlaps.empty:
            gene_states[gene] = "Euchromatin"
        elif not hetero_overlaps.empty:
            gene_states[gene] = "Heterochromatin"
        else:
            gene_states[gene] = "Unknown"
    
    return gene_states

def analyze_gene_chromatin_states(genes, results_dir, gtf_file, promoter_size=2000):
    """Analyze chromatin states for gene promoters in both NSC and Neu (GFP condition)
    
    This function focuses on the promoter regions of genes rather than the entire gene body.
    Promoters are defined as regions upstream of the transcription start site (TSS).
    """
    # Get gene promoter coordinates
    gene_promoters = get_gene_promoter_coordinates(genes, gtf_file, promoter_size)
    
    # Load chromatin states for NSC_GFP
    nsc_eu = load_chromatin_states(results_dir, "NSC_GFP", "euchromatin")
    nsc_hetero = load_chromatin_states(results_dir, "NSC_GFP", "heterochromatin")
    
    # Load chromatin states for Neu_GFP
    neu_eu = load_chromatin_states(results_dir, "Neu_GFP", "euchromatin")
    neu_hetero = load_chromatin_states(results_dir, "Neu_GFP", "heterochromatin")
    
    # Assign chromatin states for gene promoters
    nsc_states = assign_chromatin_states(gene_promoters, nsc_eu, nsc_hetero)
    neu_states = assign_chromatin_states(gene_promoters, neu_eu, neu_hetero)
    
    # Combine results
    combined_states = {}
    for gene in genes:
        if gene in gene_promoters:
            combined_states[gene] = {
                "NSC": nsc_states.get(gene, "Unknown"),
                "Neu": neu_states.get(gene, "Unknown")
            }
        else:
            combined_states[gene] = {
                "NSC": "Unknown",
                "Neu": "Unknown"
            }
    
    logger.info(f"Completed chromatin state analysis for promoters of {len(combined_states)} genes")
    return combined_states

def plot_chromatin_states(gene_states, output_dir, gene_list_name):
    """Create visualization of chromatin states for the genes"""
    # Prepare data for plotting
    data = []
    for gene, states in gene_states.items():
        data.append({
            'Gene': gene,
            'NSC': states['NSC'],
            'Neu': states['Neu']
        })
    df = pd.DataFrame(data)
    
    # Save the data to CSV
    df.to_csv(os.path.join(output_dir, 'gene_chromatin_states.csv'), index=False)
    logger.info(f"Saved gene chromatin states to CSV")
    
    # Create state distribution plot
    plt.figure(figsize=(14, 7))
    
    # Define consistent color mapping
    color_map = {
        'Unknown': '#66b3ff',      # Blue for Unknown
        'Euchromatin': '#99ff99',  # Green for Euchromatin
        'Heterochromatin': '#ff9999', # Red for Heterochromatin
        'Mixed': '#ffcc99'         # Orange for Mixed
    }
    
    # Plot for NSC
    nsc_counts = df['NSC'].value_counts()
    plt.subplot(1, 2, 1)
    
    # Map colors to the actual states present in the data
    nsc_colors = [color_map[state] for state in nsc_counts.index]
    plt.pie(nsc_counts.values, labels=nsc_counts.index, autopct='%1.1f%%', colors=nsc_colors)
    plt.title('Chromatin States Distribution in NSC (GFP)')
    
    # Plot for Neu
    neu_counts = df['Neu'].value_counts()
    plt.subplot(1, 2, 2)
    
    # Map colors to the actual states present in the data
    neu_colors = [color_map[state] for state in neu_counts.index]
    plt.pie(neu_counts.values, labels=neu_counts.index, autopct='%1.1f%%', colors=neu_colors)
    plt.title('Chromatin States Distribution in Neu (GFP)')
    
    plt.suptitle(f'Chromatin State Distribution for {gene_list_name} Genes', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromatin_states_distribution.png'), dpi=300)
    plt.close()
    logger.info(f"Created chromatin states distribution plot")
    
    # Create state transition plot
    transitions = {}
    for gene, states in gene_states.items():
        transition = f"{states['NSC']} → {states['Neu']}"
        transitions[transition] = transitions.get(transition, 0) + 1
    
    # Sort transitions to group them meaningfully
    sorted_transitions = {}
    # First add transitions that maintain the same state
    for state in ['Euchromatin', 'Heterochromatin', 'Mixed', 'Unknown']:
        key = f"{state} → {state}"
        if key in transitions:
            sorted_transitions[key] = transitions[key]
    
    # Then add transitions from open to closed chromatin
    key = "Euchromatin → Heterochromatin"
    if key in transitions:
        sorted_transitions[key] = transitions[key]
    
    # Then add transitions from closed to open chromatin
    key = "Heterochromatin → Euchromatin"
    if key in transitions:
        sorted_transitions[key] = transitions[key]
    
    # Add all other transitions
    for key, value in transitions.items():
        if key not in sorted_transitions:
            sorted_transitions[key] = value
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(sorted_transitions.keys(), sorted_transitions.values(), color='skyblue')
    plt.xlabel('Chromatin State Transition (NSC → Neu)')
    plt.ylabel('Number of Genes')
    plt.title(f'Chromatin State Transitions for {gene_list_name} Genes')
    plt.xticks(rotation=45, ha='right')
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                 f'{height}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromatin_state_transitions.png'), dpi=300)
    plt.close()
    logger.info(f"Created chromatin state transitions plot")
    
    # Create a focused plot on chromatin state changes (open to closed and vice versa)
    plt.figure(figsize=(10, 6))
    
    # Define the key transitions we're interested in
    key_transitions = {
        'Euchromatin → Heterochromatin': 'Open to Closed',
        'Heterochromatin → Euchromatin': 'Closed to Open',
        'Euchromatin → Euchromatin': 'Remained Open',
        'Heterochromatin → Heterochromatin': 'Remained Closed'
    }
    
    # Define consistent colors for transitions
    transition_color_map = {
        'Open to Closed': '#ff9999',      # Red for closing (like heterochromatin)
        'Closed to Open': '#99ff99',      # Green for opening (like euchromatin)
        'Remained Open': '#99cc99',       # Darker green for remained open
        'Remained Closed': '#cc9999'      # Darker red for remained closed
    }
    
    # Extract counts for these transitions
    transition_counts = []
    transition_labels = []
    transition_colors = []
    
    for key, label in key_transitions.items():
        if key in transitions:
            transition_counts.append(transitions[key])
            transition_labels.append(label)
            transition_colors.append(transition_color_map[label])
        else:
            transition_counts.append(0)
            transition_labels.append(label)
            transition_colors.append(transition_color_map[label])
    
    # Create the bar plot
    bars = plt.bar(transition_labels, transition_counts, color=transition_colors)
    plt.xlabel('Chromatin State Change')
    plt.ylabel('Number of Genes')
    plt.title(f'Key Chromatin State Changes from NSC to Neu')
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:  # Only add text if the bar has height
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{height}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'key_chromatin_state_changes.png'), dpi=300)
    plt.close()
    logger.info(f"Created key chromatin state changes plot")
    
    # Create a pie chart showing the proportion of genes that changed state vs. remained the same
    changed_state = 0
    same_state = 0
    unknown_state = 0
    
    for gene, states in gene_states.items():
        if states['NSC'] == 'Unknown' or states['Neu'] == 'Unknown':
            unknown_state += 1
        elif states['NSC'] == states['Neu']:
            same_state += 1
        else:
            changed_state += 1
    
    # Use consistent colors - blue for unknown, red for changed (like heterochromatin), green for maintained (like euchromatin)
    change_color_map = {
        'Changed State': '#ff9999',    # Red for changes
        'Maintained State': '#99ff99', # Green for maintained
        'Unknown': '#66b3ff'           # Blue for unknown
    }
    
    plt.figure(figsize=(8, 8))
    plt.pie([changed_state, same_state, unknown_state], 
            labels=['Changed State', 'Maintained State', 'Unknown'],
            autopct='%1.1f%%', colors=[change_color_map['Changed State'], 
                                      change_color_map['Maintained State'], 
                                      change_color_map['Unknown']])
    plt.title(f'Proportion of Genes Changing Chromatin State from NSC to Neu')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromatin_state_change_proportion.png'), dpi=300)
    plt.close()
    logger.info(f"Created chromatin state change proportion plot")

def main():
    # Parse arguments
    args = parse_arguments()
    
    # Extract gene list name from file path
    gene_list_name = os.path.basename(args.gene_list).split('.')[0]
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    logger.info(f"Analysis output will be saved to {args.output_dir}")
    logger.info(f"Analyzing promoter regions with size: {args.promoter_size} bp")
    
    # Load gene list
    genes = load_gene_list(args.gene_list)
    
    # Analyze chromatin states in promoter regions
    gene_states = analyze_gene_chromatin_states(genes, args.results_dir, args.genome_gtf, args.promoter_size)
    
    # Create visualizations
    plot_chromatin_states(gene_states, args.output_dir, gene_list_name)
    
    logger.info(f"Analysis complete. Results saved in {args.output_dir}")

if __name__ == "__main__":
    main()
