#!/usr/bin/env python3
"""
SAMMY-seq Mecp2 Target Comparison - Parallel Metaprofile Computation

This script computes metaprofiles for a specific gene set and condition in parallel,
saving the results to a file for later combination.

It is designed to be run as part of a SLURM array job, with each job processing
a specific combination of gene set and condition.
"""

import os
import sys
import argparse
import logging
import numpy as np
import pyBigWig
import pybedtools
from tqdm import tqdm
import pickle

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("metaprofile_parallel.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('SAMMY-seq-Parallel-Metaprofile')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='SAMMY-seq Parallel Metaprofile Computation')
    parser.add_argument('--results-dir', required=True, help='Directory containing SAMMY-seq analysis results')
    parser.add_argument('--output-dir', required=True, help='Output directory for metaprofiles')
    parser.add_argument('--genome-gtf', required=True, help='GTF file with gene annotations')
    parser.add_argument('--gene-list-file', required=True, help='File containing gene list')
    parser.add_argument('--gene-set-name', required=True, help='Name of the gene set (e.g., NEU_target)')
    parser.add_argument('--condition', required=True, help='Condition to process (e.g., NSC_GFP)')
    parser.add_argument('--upstream-distance', type=int, default=5000, help='Distance upstream of TSS (default: 5000)')
    parser.add_argument('--downstream-distance', type=int, default=5000, help='Distance downstream of TES (default: 5000)')
    parser.add_argument('--bin-size', type=int, default=100, help='Bin size for metaprofile (default: 100)')
    parser.add_argument('--gene-body-bins', type=int, default=100, help='Number of bins for gene body (default: 100)')
    return parser.parse_args()

def load_gene_list(file_path):
    """Load gene list from file"""
    logger.info(f"Loading gene list from {file_path}")
    try:
        with open(file_path, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(genes)} genes from {file_path}")
        return genes
    except Exception as e:
        logger.error(f"Error loading gene list from {file_path}: {e}")
        return []

def get_gene_coordinates(genes, gtf_file):
    """
    Get gene coordinates (TSS, TES, strand) from GTF file.
    """
    logger.info(f"Getting coordinates for {len(genes)} genes from GTF file")
    
    gene_coords = {}
    
    try:
        gtf = pybedtools.BedTool(gtf_file)
        gene_features = gtf.filter(lambda x: x[2] == 'gene')
        
        for feature in gene_features:
            attributes = dict(item.strip().split(' ') for item in feature.fields[8].split(';') if item.strip())
            
            gene_name = None
            for attr in ['gene_name', 'gene_id']:
                if attr in attributes:
                    gene_name = attributes[attr].strip('"\'')
                    break
            
            if gene_name and gene_name in genes:
                chrom = feature.chrom
                start = int(feature.start)
                end = int(feature.end)
                strand = feature.strand
                
                if strand == '+':
                    tss = start
                    tes = end
                else:
                    tss = end
                    tes = start
                
                gene_coords[gene_name] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'tss': tss,
                    'tes': tes
                }
        
        logger.info(f"Found coordinates for {len(gene_coords)} out of {len(genes)} genes")
        return gene_coords
    
    except Exception as e:
        logger.error(f"Error getting gene coordinates: {e}")
        return {}

def compute_metaprofile(gene_coords, bigwig_file, upstream_distance, downstream_distance, bin_size, gene_body_bins):
    """
    Compute metaprofile for a set of genes from a BigWig file.
    """
    logger.info(f"Computing metaprofile from {bigwig_file}")
    
    try:
        bw = pyBigWig.open(bigwig_file)
        
        upstream_bins = int(upstream_distance / bin_size)
        downstream_bins = int(downstream_distance / bin_size)
        
        total_bins = upstream_bins + gene_body_bins + downstream_bins
        
        matrix = np.zeros((len(gene_coords), total_bins))
        matrix.fill(np.nan)
        
        for i, (gene_name, coords) in enumerate(tqdm(gene_coords.items(), desc="Processing genes")):
            chrom = coords['chrom']
            tss = coords['tss']
            tes = coords['tes']
            strand = coords['strand']
            
            if chrom not in bw.chroms():
                continue
            
            if strand == '+':
                upstream_start = max(0, tss - upstream_distance)
                upstream_end = tss
                gene_body_start = tss
                gene_body_end = tes
                downstream_start = tes
                downstream_end = min(bw.chroms()[chrom], tes + downstream_distance)
            else:
                upstream_start = max(0, tes)
                upstream_end = min(bw.chroms()[chrom], tes + upstream_distance)
                gene_body_start = tes
                gene_body_end = tss
                downstream_start = max(0, tss - downstream_distance)
                downstream_end = tss
            
            if upstream_end > upstream_start:
                try:
                    upstream_values = bw.stats(chrom, upstream_start, upstream_end, 
                                              nBins=upstream_bins, type="mean")
                    
                    if strand == '-':
                        upstream_values = upstream_values[::-1] if upstream_values else upstream_values
                    
                    for j, val in enumerate(upstream_values):
                        if val is not None:
                            matrix[i, j] = val
                except Exception as e:
                    logger.debug(f"Error processing upstream region for {gene_name}: {e}")
            
            if gene_body_end > gene_body_start:
                try:
                    gene_body_values = bw.stats(chrom, gene_body_start, gene_body_end, 
                                               nBins=gene_body_bins, type="mean")
                    
                    if strand == '-':
                        gene_body_values = gene_body_values[::-1] if gene_body_values else gene_body_values
                    
                    for j, val in enumerate(gene_body_values):
                        if val is not None:
                            matrix[i, upstream_bins + j] = val
                except Exception as e:
                    logger.debug(f"Error processing gene body for {gene_name}: {e}")
            
            if downstream_end > downstream_start:
                try:
                    downstream_values = bw.stats(chrom, downstream_start, downstream_end, 
                                                nBins=downstream_bins, type="mean")
                    
                    if strand == '-':
                        downstream_values = downstream_values[::-1] if downstream_values else downstream_values
                    
                    for j, val in enumerate(downstream_values):
                        if val is not None:
                            matrix[i, upstream_bins + gene_body_bins + j] = val
                except Exception as e:
                    logger.debug(f"Error processing downstream region for {gene_name}: {e}")
        
        bw.close()
        
        return matrix
    
    except Exception as e:
        logger.error(f"Error computing metaprofile: {e}")
        return None

def main():
    """Main function"""
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load gene list
    genes = load_gene_list(args.gene_list_file)
    
    if not genes:
        logger.error(f"No genes found in {args.gene_list_file}, exiting")
        sys.exit(1)
    
    # Get gene coordinates
    gene_coords = get_gene_coordinates(genes, args.genome_gtf)
    
    if not gene_coords:
        logger.error(f"No gene coordinates found for {args.gene_list_file}, exiting")
        sys.exit(1)
    
    # Find the ratio BigWig file
    ratio_file = os.path.join(args.results_dir, "ratio", f"{args.condition}_S2S_vs_S3_ratio.bw")
    
    if not os.path.exists(ratio_file):
        logger.error(f"Ratio file not found: {ratio_file}, exiting")
        sys.exit(1)
    
    # Compute metaprofile
    matrix = compute_metaprofile(
        gene_coords, 
        ratio_file, 
        args.upstream_distance, 
        args.downstream_distance, 
        args.bin_size, 
        args.gene_body_bins
    )
    
    if matrix is None:
        logger.error("Failed to compute metaprofile, exiting")
        sys.exit(1)
    
    # Apply data normalization to enhance differences
    # Center the data around zero
    matrix = matrix - np.nanmean(matrix)
    
    # Save the matrix to a file
    output_file = os.path.join(args.output_dir, f"{args.gene_set_name}_{args.condition}_matrix.pkl")
    
    with open(output_file, 'wb') as f:
        pickle.dump({
            'matrix': matrix,
            'gene_set': args.gene_set_name,
            'condition': args.condition,
            'upstream_distance': args.upstream_distance,
            'downstream_distance': args.downstream_distance,
            'bin_size': args.bin_size,
            'gene_body_bins': args.gene_body_bins,
            'gene_coords': gene_coords
        }, f)
    
    logger.info(f"Saved metaprofile matrix to {output_file}")

if __name__ == "__main__":
    main() 