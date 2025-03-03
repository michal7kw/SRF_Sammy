def benjamini_hochberg_correction(pvalues):
    """
    Implements Benjamini-Hochberg FDR correction for multiple hypothesis testing.
    
    Args:
        pvalues: numpy array or pandas Series of p-values
        
    Returns:
        numpy array of adjusted p-values (q-values)
    """
    if isinstance(pvalues, pd.Series):
        pvalues = pvalues.values
    
    n = len(pvalues)
    if n == 0:
        return np.array([])
    
    # Create ranks (1-based)
    ranks = np.argsort(pvalues)
    sorted_pvalues = pvalues[ranks]
    ranks_rev = np.argsort(ranks)
    
    # Calculate adjusted p-values
    adjusted = np.minimum(1, sorted_pvalues * n / np.arange(1, n + 1))
    
    # Ensure monotonicity
    for i in range(n-2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i+1])
    
    # Return to original order
    return adjusted[ranks_rev]


def plot_results(df, cell_type, output_dir):
    """Create plots to visualize the results."""
    try:
        # Check which columns are available in the dataframe
        available_columns = df.columns.tolist()
        logging.info(f"Available columns for plotting: {available_columns}")
        
        # Create a figure with multiple subplots
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'Chromatin State Changes in {cell_type}', fontsize=16)
        
        # Plot 1: Distribution of condition differences
        if 'condition_diff' in available_columns:
            sns.histplot(data=df, x='condition_diff', bins=50, ax=axes[0, 0])
            axes[0, 0].set_title('Distribution of Condition Differences')
            axes[0, 0].set_xlabel('M2 - GFP Difference')
            axes[0, 0].set_ylabel('Count')
        else:
            axes[0, 0].text(0.5, 0.5, 'Column "condition_diff" not available', 
                          horizontalalignment='center', verticalalignment='center')
        
        # Plot 2: Scatter plot of GFP vs M2 differences
        if 'GFP_diff' in available_columns and 'M2_diff' in available_columns:
            axes[0, 1].scatter(df['GFP_diff'], df['M2_diff'], alpha=0.1)
            axes[0, 1].set_title('GFP vs M2 Differences')
            axes[0, 1].set_xlabel('GFP Difference (S3 - S2S)')
            axes[0, 1].set_ylabel('M2 Difference (S3 - S2S)')
            # Add diagonal line
            lims = [
                min(axes[0, 1].get_xlim()[0], axes[0, 1].get_ylim()[0]),
                max(axes[0, 1].get_xlim()[1], axes[0, 1].get_ylim()[1])
            ]
            axes[0, 1].plot(lims, lims, 'k--', alpha=0.5, zorder=0)
        else:
            axes[0, 1].text(0.5, 0.5, 'Columns "GFP_diff" or "M2_diff" not available', 
                          horizontalalignment='center', verticalalignment='center')
        
        # Plot 3: Change type distribution
        if 'change_type' in available_columns:
            change_counts = df['change_type'].value_counts()
            axes[1, 0].bar(change_counts.index, change_counts.values)
            axes[1, 0].set_title('Distribution of Change Types')
            axes[1, 0].set_xlabel('Change Type')
            axes[1, 0].set_ylabel('Count')
            # Rotate x-axis labels if needed
            plt.setp(axes[1, 0].get_xticklabels(), rotation=45, ha='right')
        else:
            axes[1, 0].text(0.5, 0.5, 'Column "change_type" not available', 
                          horizontalalignment='center', verticalalignment='center')
        
        # Plot 4: Significant changes by chromosome
        if 'chrom' in available_columns and 'is_significant' in available_columns:
            # Group by chromosome and count significant changes
            chrom_counts = df[df['is_significant']].groupby('chrom').size().reset_index(name='count')
            
            # Sort chromosomes naturally
            def chrom_to_order(chrom):
                if isinstance(chrom, int):
                    return chrom
                if str(chrom).isdigit():
                    return int(chrom)
                if chrom == 'X':
                    return 23
                if chrom == 'Y':
                    return 24
                return 25  # For any other chromosomes
                
            chrom_counts['chrom_order'] = chrom_counts['chrom'].apply(chrom_to_order)
            chrom_counts = chrom_counts.sort_values('chrom_order')
            
            axes[1, 1].bar(chrom_counts['chrom'], chrom_counts['count'])
            axes[1, 1].set_title('Significant Changes by Chromosome')
            axes[1, 1].set_xlabel('Chromosome')
            axes[1, 1].set_ylabel('Count')
            # Rotate x-axis labels if needed
            plt.setp(axes[1, 1].get_xticklabels(), rotation=45, ha='right')
        else:
            axes[1, 1].text(0.5, 0.5, 'Columns "chrom" or "is_significant" not available', 
                          horizontalalignment='center', verticalalignment='center')
        
        # Adjust layout and save
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.savefig(output_dir / f"{cell_type}_chromatin_changes_plot.png", dpi=300)
        plt.close()
        
        logging.info(f"Plots created for {cell_type}")
    except Exception as e:
        logging.error(f"Error creating plots for {cell_type}: {str(e)}")