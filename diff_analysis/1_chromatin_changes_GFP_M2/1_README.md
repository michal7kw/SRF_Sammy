# Input

- *No* gene list

# Analysis format

```
gfp_s2s_signal = compute_average_signal(bigwig_files_gfp_s2s, chrom, start, end)
gfp_s3_signal = compute_average_signal(bigwig_files_gfp_s3, chrom, start, end)
m2_s2s_signal = compute_average_signal(bigwig_files_m2_s2s, chrom, start, end)
m2_s3_signal = compute_average_signal(bigwig_files_m2_s3, chrom, start, end)
```

Calculate the difference in signal between S3 and S2S for each condition
```
gfp_diff = gfp_s3_signal - gfp_s2s_signal
m2_diff = m2_s3_signal - m2_s2s_signal
```

Calculate the difference between the M2 and GFP signal differences
```
condition_diff = m2_diff - gfp_diff
```

Store the results for the current window in a dictionary
```
results.append({
    'chrom': chrom,              # Chromosome name
    'start': start,              # Start position of the window
    'end': end,                  # End position of the window
    'gfp_s2s': gfp_s2s_signal,    # Average signal for GFP S2S
    'gfp_s3': gfp_s3_signal,      # Average signal for GFP S3
    'gfp_diff': gfp_diff,        # Difference between GFP S3 and S2S
    'm2_s2s': m2_s2s_signal,      # Average signal for M2 S2S
    'm2_s3': m2_s3_signal,        # Average signal for M2 S3
    'm2_diff': m2_diff,          # Difference between M2 S3 and S2S
    'condition_diff': condition_diff  # Difference between M2 and GFP signal differences
})
```
