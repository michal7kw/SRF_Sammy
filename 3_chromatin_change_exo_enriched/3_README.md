
# Analyze chromatin state changes in genes where exogenous Mecp2 is enriched




Define promoter region (5kb upstream of TSS)
```
if coords['strand'] == '+':  # If the gene is on the positive strand
    promoter_start = max(0, coords['start'] - WINDOW_SIZE)  # Define the promoter start coordinate
    promoter_end = coords['start']  # Define the promoter end coordinate
else:  # If the gene is on the negative strand
    promoter_start = coords['end']  # Define the promoter start coordinate
    promoter_end = coords['end'] + WINDOW_SIZE  # Define the promoter end coordinate
```


Get bigWig files for this cell type
```
mecp2_s2s_files = get_bigwig_files(bigwig_dir, cell_type, "M2", "S2S")  # Get the BigWig files for the Mecp2 S2S condition
mecp2_s3_files = get_bigwig_files(bigwig_dir, cell_type, "M2", "S3")  # Get the BigWig files for the Mecp2 S3 condition
gfp_s2s_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S2S")  # Get the BigWig files for the GFP S2S condition
gfp_s3_files = get_bigwig_files(bigwig_dir, cell_type, "GFP", "S3")  # Get the BigWig files for the GFP S3 condition
```    


Compute average signals
```
s2s_signal = compute_average_signal(s2s_files, chrom, promoter_start, promoter_end)  # Compute the average signal for the S2S condition
s3_signal = compute_average_signal(s3_files, chrom, promoter_start, promoter_end)  # Compute the average signal for the S3 condition
```

Determine chromatin state
```
diff = s3_signal - s2s_signal  # Compute the diff of S3 to S2S signals

# Classify based on diff
if diff > 1 + MIN_DIFF_THRESHOLD:  # If the diff is greater than 1 + the minimum difference threshold
    state = "Heterochromatin"  # S3 dominant. Classify as heterochromatin
elif diff < 1 - MIN_DIFF_THRESHOLD:  # If the diff is less than 1 - the minimum difference threshold
    state = "Euchromatin"  # S2S dominant. Classify as euchromatin
else:  # Otherwise
    state = "Mixed"  # No clear dominance. Classify as mixed
```