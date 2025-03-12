#!/usr/bin/env python3
"""
Associate Shifted Regions with Genes

This script takes the BED files containing shifted chromatin regions and associates them
with genes from the mm10 genome using the provided GTF annotation file.
"""

import os
import sys
import pandas as pd
import logging
from collections import defaultdict
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger()

def parse_gtf(gtf_file):
    """
    Parse GTF file and extract gene information.
    
    Returns a dictionary with gene coordinates and information.
    """
    logger.info(f"Parsing GTF file: {gtf_file}")
    genes = {}
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                
                # Parse attributes
                attr_dict = {}
                for attr in fields[8].split(';'):
                    if not attr.strip():
                        continue
                    key, value = attr.strip().split(' ', 1)
                    attr_dict[key] = value.strip('"')
                
                gene_id = attr_dict.get('gene_id', '')
                gene_name = attr_dict.get('gene_name', gene_id)
                gene_type = attr_dict.get('gene_type', '')
                
                genes[gene_id] = {
                    'chrom': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'strand': fields[6],
                    'gene_name': gene_name,
                    'gene_type': gene_type
                }
        
        logger.info(f"Successfully parsed {len(genes)} genes from GTF file")
        return genes
    
    except Exception as e:
        logger.error(f"Error parsing GTF file: {e}")
        sys.exit(1)

def find_overlapping_genes(region, genes, upstream=5000, downstream=5000):
    """
    Find genes that overlap with a given region, including promoter regions.
    
    Parameters:
    -----------
    region : dict
        Dictionary containing region information (chrom, start, end)
    genes : dict
        Dictionary of gene information
    upstream : int
        Number of base pairs upstream to consider as promoter region
    downstream : int
        Number of base pairs downstream to consider as promoter region
    
    Returns:
    --------
    list
        List of overlapping gene IDs and their overlap type
    """
    overlapping = []
    
    for gene_id, gene in genes.items():
        if gene['chrom'] != region['chrom']:
            continue
        
        # Define gene boundaries including promoter regions
        if gene['strand'] == '+':
            gene_start = gene['start'] - upstream
            gene_end = gene['end'] + downstream
        else:
            gene_start = gene['start'] - downstream
            gene_end = gene['end'] + upstream
        
        # Check for overlap
        if not (region['end'] < gene_start or region['start'] > gene_end):
            # Determine overlap type
            if region['start'] >= gene['start'] and region['end'] <= gene['end']:
                overlap_type = 'gene_body'
            elif region['end'] >= gene_start and region['end'] <= gene['start']:
                overlap_type = 'promoter_upstream'
            elif region['start'] <= gene_end and region['start'] >= gene['end']:
                overlap_type = 'promoter_downstream'
            else:
                overlap_type = 'partial'
            
            overlapping.append({
                'gene_id': gene_id,
                'gene_name': gene['gene_name'],
                'gene_type': gene['gene_type'],
                'overlap_type': overlap_type,
                'strand': gene['strand']
            })
    
    return overlapping

def process_bed_file(bed_file, genes, output_file):
    """
    Process a BED file and find overlapping genes for each region.
    """
    logger.info(f"Processing BED file: {bed_file}")
    
    try:
        # Read BED file
        regions_df = pd.read_csv(bed_file, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'ratio_change'])
        
        if len(regions_df) == 0:
            logger.warning(f"No regions found in {bed_file}")
            return
        
        # Process each region
        results = []
        for _, region in regions_df.iterrows():
            overlapping = find_overlapping_genes(region, genes)
            
            for gene in overlapping:
                results.append({
                    'chrom': region['chrom'],
                    'start': region['start'],
                    'end': region['end'],
                    'ratio_change': region['ratio_change'],
                    'gene_id': gene['gene_id'],
                    'gene_name': gene['gene_name'],
                    'gene_type': gene['gene_type'],
                    'overlap_type': gene['overlap_type'],
                    'strand': gene['strand']
                })
        
        # Create results DataFrame and save
        if results:
            results_df = pd.DataFrame(results)
            results_df.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Found {len(results)} gene associations")
            logger.info(f"Results saved to: {output_file}")
            
            # Create summary statistics
            summary_file = output_file.replace('.tsv', '_summary.tsv')
            gene_type_summary = results_df.groupby('gene_type').size().reset_index(name='count')
            gene_type_summary.to_csv(summary_file, sep='\t', index=False)
            logger.info(f"Summary statistics saved to: {summary_file}")
        else:
            logger.warning("No gene associations found")
            
    except Exception as e:
        logger.error(f"Error processing BED file: {e}")
        return

def main():
    parser = argparse.ArgumentParser(description='Associate shifted regions with genes')
    parser.add_argument('--results-dir', required=True, help='Directory containing shifted regions BED files')
    parser.add_argument('--gtf-file', required=True, 
                       default='/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/DATA/gencode.vM25.basic.annotation.gtf',
                       help='Path to GTF annotation file')
    parser.add_argument('--upstream', type=int, default=5000,
                       help='Number of base pairs upstream to consider as promoter region')
    parser.add_argument('--downstream', type=int, default=5000,
                       help='Number of base pairs downstream to consider as promoter region')
    args = parser.parse_args()
    
    # Parse GTF file
    genes = parse_gtf(args.gtf_file)
    
    # Create output directory
    output_dir = os.path.join(args.results_dir, 'gene_associations')
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each BED file in the shifted regions directory
    for filename in os.listdir(args.results_dir):
        if filename.endswith('.bed'):
            bed_file = os.path.join(args.results_dir, filename)
            output_file = os.path.join(output_dir, filename.replace('.bed', '_genes.tsv'))
            process_bed_file(bed_file, genes, output_file)
    
    logger.info("Analysis complete")

if __name__ == "__main__":
    main() 