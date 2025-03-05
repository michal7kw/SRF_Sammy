#!/bin/bash

# First, extract a few sequences from your FASTQ file
zcat DATA/Sammy_Seq_fastq/NSC2_M2_S4_4F_B3_S18_L001_R1_001.fastq.gz | head -n 4000 > sample_reads.fastq

# Convert the first few FASTQ entries to FASTA format
seqtk seq -A sample_reads.fastq > sample_reads.fasta

# Now use these sequences for BLAST through NCBI web interface
# Or use command line BLAST