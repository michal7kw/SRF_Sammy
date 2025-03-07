# Input

- ("endo_list_neu_tss_peaks.csv", "endo_gene_list_neu.txt")
- ("endo_list_nsc_tss_peaks.csv", "endo_gene_list_nsc.txt")

# Analysis format

Get GFP bigWig files for this cell type
```
s2s_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S2S")  # Get the S2S BigWig files
s3_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S3")  # Get the S3 BigWig files
```

Define promoter region (5kb upstream of TSS)
```
if coords['strand'] == '+':  # If the gene is on the positive strand
    promoter_start = max(0, coords['start'] - WINDOW_SIZE)  # Calculate the promoter start coordinate
    promoter_end = coords['start']  # The promoter end coordinate is the TSS
else:  # '-' strand  # If the gene is on the negative strand
    promoter_start = coords['end']  # The promoter start coordinate is the end of the gene
    promoter_end = coords['end'] + WINDOW_SIZE  # Calculate the promoter end coordinate
```

Compute average signals
```
s2s_signal = compute_average_signal(s2s_files, chrom, promoter_start, promoter_end)  # Compute the average S2S signal
s3_signal = compute_average_signal(s3_files, chrom, promoter_start, promoter_end)  # Compute the average S3 signal
```

Determine chromatin state
```
if s2s_signal is not None and s3_signal is not None:  # If both signals were computed
    diff = s3_signal - s2s_signal  # Calculate the diff of S3 to S2S signal
    
    # Classify based on diff
    if diff > MIN_DIFF_THRESHOLD:  # If the diff is significantly greater than 1
        state = "Heterochromatin"  # S3 dominant  # Classify as heterochromatin
    elif diff < -MIN_DIFF_THRESHOLD:  # If the diff is significantly less than 1
        state = "Euchromatin"  # S2S dominant  # Classify as euchromatin
    else:  # If the diff is not significantly different from 1
        state = "Mixed"  # No clear dominance  # Classify as mixed
    
    result = {  # Create a dictionary representing the result
        'gene': gene,
        'chrom': chrom,
        'start': coords['start'],
        'end': coords['end'],
        'strand': coords['strand'],
        'promoter_start': promoter_start,
        'promoter_end': promoter_end,
        's2s_signal': s2s_signal,
        's3_signal': s3_signal,
        'diff': diff,
        'state': state
    }

```