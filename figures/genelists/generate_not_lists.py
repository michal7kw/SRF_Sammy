#!/usr/bin/env python3

import os
import re

# Define paths
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ENDO_FILE = os.path.join(CURRENT_DIR, "endo.txt")
NOT_ENDO_FILE = os.path.join(CURRENT_DIR, "not_endo.txt")
GTF_FILE = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Sammy/DATA/gencode.vM25.basic.annotation.gtf"

def read_gene_list(file_path):
    """Read gene list from a file, one gene per line."""
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def is_valid_gene_name(gene_name):
    """
    Check if a gene name is valid based on common naming conventions.
    
    Valid gene names typically:
    - Start with a letter (not a number)
    - Don't contain special characters except hyphens and dots
    - Are not too short (at least 2 characters)
    - Are not temporary/predicted IDs (like ENSMUSG...)
    """
    # Skip very short names
    if len(gene_name) < 2:
        return False
    
    # Skip names that start with numbers (like 0610005C13Rik)
    if gene_name[0].isdigit():
        return False
    
    # Skip Ensembl IDs
    if gene_name.startswith(('ENSMUSG', 'ENSMUST')):
        return False
    
    # Skip names with unusual characters
    if not re.match(r'^[A-Za-z0-9\.\-_]+$', gene_name):
        return False
    
    return True

def extract_genes_from_gtf(gtf_file):
    """Extract canonical gene names from GTF file."""
    genes = set()
    gene_name_pattern = re.compile(r'gene_name "([^"]+)"')
    gene_type_pattern = re.compile(r'gene_type "([^"]+)"')
    
    # Gene types to include (protein-coding genes and well-characterized non-coding genes)
    valid_gene_types = {
        'protein_coding', 
        'lincRNA', 
        'antisense', 
        'processed_transcript'
    }
    
    print(f"Parsing GTF file: {gtf_file}")
    with open(gtf_file, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            
            # Only process gene entries, not transcripts/exons
            fields = line.split('\t')
            if len(fields) >= 3 and fields[2] == 'gene':
                
                # Check gene type
                gene_type_match = gene_type_pattern.search(line)
                if gene_type_match:
                    gene_type = gene_type_match.group(1)
                    if gene_type not in valid_gene_types:
                        continue
                
                # Extract gene name
                match = gene_name_pattern.search(line)
                if match:
                    gene_name = match.group(1)
                    if is_valid_gene_name(gene_name):
                        genes.add(gene_name)
            
            if i % 100000 == 0:
                print(f"Processed {i} lines...")
    
    return genes

def create_not_list(all_genes, genes_in_list):
    """Create a list of genes not in the specified list."""
    return all_genes - genes_in_list

def write_gene_list(genes, output_file):
    """Write gene list to a file, one gene per line."""
    with open(output_file, 'w') as f:
        for gene in sorted(genes):
            f.write(f"{gene}\n")

def main():
    # Read endo gene list
    print(f"Reading gene list from {ENDO_FILE}")
    endo_genes = read_gene_list(ENDO_FILE)
    print(f"Found {len(endo_genes)} genes in endo list")
    
    # Extract all genes from GTF
    all_genes = extract_genes_from_gtf(GTF_FILE)
    print(f"Found {len(all_genes)} valid canonical genes in GTF file")
    
    # Create not_endo list
    not_endo_genes = create_not_list(all_genes, endo_genes)
    print(f"Generated {len(not_endo_genes)} genes for not_endo list")
    
    # Write not_endo list to file
    write_gene_list(not_endo_genes, NOT_ENDO_FILE)
    print(f"Wrote not_endo list to {NOT_ENDO_FILE}")

if __name__ == "__main__":
    main()
