# Chromatin State Analysis

This repository contains scripts for analyzing chromatin state changes between different cell types (NSCs and Neurons) and conditions (GFP and M2).

## Overview

The analysis consists of two main steps:

1. **Chromatin State Analysis**: Analyzes bigWig files to identify regions where chromatin state changes between conditions.
2. **Chromatin State Comparison**: Compares chromatin states between different cell types and conditions, and generates visualizations.

## Scripts

### 1. Chromatin State Analysis

- `1_analyze_chromatin_changes_GFP_M2.py`: Analyzes chromatin state changes between GFP and M2 conditions for NSCs and Neurons.
- `1_analyze_chromatin_changes_GFP_M2.sh`: SLURM job script to run the analysis for a specific chromosome.

### 2. Chromatin State Comparison

- `compare_chromatin_states.py`: Compares chromatin states between different cell types and conditions, and generates visualizations.
- `run_chromatin_comparison.sh`: SLURM job script to run the comparison analysis.

### 3. NSC vs Neuron Comparison

- `compare_nsc_neuron.py`: Specifically compares chromatin states between NSCs and Neurons in GFP condition, with detailed visualizations.
- `run_nsc_neuron_comparison.sh`: SLURM job script to run the NSC vs Neuron comparison analysis.

### 4. Full Analysis

- `run_full_analysis.sh`: SLURM job script to run both the chromatin state analysis and comparison.

## Running the Analysis

### Option 1: Run the full analysis

```bash
sbatch run_full_analysis.sh
```

This will:
1. Run the chromatin state analysis for all chromosomes
2. Run the chromatin state comparison
3. Generate visualizations

### Option 2: Run the analysis for a specific chromosome

```bash
sbatch --array=1-21 1_analyze_chromatin_changes_GFP_M2.sh
```

This will run the chromatin state analysis for all chromosomes in parallel (1-19, X, Y).

### Option 3: Run only the comparison

```bash
sbatch run_chromatin_comparison.sh
```

This will run the chromatin state comparison and generate visualizations.

### Option 4: Run only the NSC vs Neuron comparison

```bash
sbatch run_nsc_neuron_comparison.sh
```

This will run the detailed comparison between NSCs and Neurons in GFP condition and generate visualizations.

## Output

The analysis generates the following output:

- `results/analysis/`: Contains the results of the chromatin state analysis.
  - `NSC_chromatin_changes_all.csv`: All chromatin changes for NSCs.
  - `Neu_chromatin_changes_all.csv`: All chromatin changes for Neurons.
  - `NSC_chromatin_changes_S2S_to_S3.csv`: Regions where chromatin state changes from S2S to S3 in NSCs.
  - `NSC_chromatin_changes_S3_to_S2S.csv`: Regions where chromatin state changes from S3 to S2S in NSCs.
  - `Neu_chromatin_changes_S2S_to_S3.csv`: Regions where chromatin state changes from S2S to S3 in Neurons.
  - `Neu_chromatin_changes_S3_to_S2S.csv`: Regions where chromatin state changes from S3 to S2S in Neurons.

- `results/comparison/`: Contains the results of the chromatin state comparison.
  - `chromatin_transitions.png`: Visualization of chromatin state transitions.
  - `chromatin_transitions.pdf`: PDF version of the visualization.
  - `chromatin_transitions_percentages.csv`: Data used to generate the visualization.

- `results/nsc_neuron_comparison/`: Contains the detailed results of the NSC vs Neuron comparison.
  - `nsc_neuron_comparison.csv`: Detailed data on chromatin state transitions between NSCs and Neurons.
  - `nsc_neuron_transition_counts.csv`: Counts and percentages of each transition type.
  - `nsc_neuron_transitions_pie.png`: Pie chart showing the distribution of transition types.
  - `nsc_neuron_transitions_bar.png`: Horizontal bar chart showing the percentage of regions in each transition type.
  - `nsc_neuron_signal_by_transition.png`: Bar chart showing the average signal by transition type.
  - `nsc_neuron_transitions_by_chromosome.png`: Heatmap showing the distribution of transitions across chromosomes.

## Interpretation

The visualizations show the percentage of genomic regions in different chromatin states:

- **A**: Euchromatin (S2S > S3)
- **B**: Heterochromatin (S3 > S2S)
- **A→B**: Transition from euchromatin to heterochromatin
- **B→A**: Transition from heterochromatin to euchromatin

The comparison between NSCs and Neurons in GFP condition shows how many regions switch from A (euchromatin) to B (heterochromatin) during neuronal differentiation.

The comparisons between GFP and M2 conditions within each cell type demonstrate that only a few regions modify their chromatin state when comparing these conditions.

## Visualization Notes

The visualizations are designed to show the data with appropriate scaling:

1. The x-axis limits in the bar charts are automatically adjusted based on the actual data range, ensuring that the data is clearly visible and not compressed into a small portion of the plot.

2. The stacked bar charts show the percentage of genomic regions in each chromatin state category, with the total adding up to 100% for each cell type or comparison.

3. The pie charts show the relative proportions of different transition types, filtering out very small percentages for clarity.

4. The heatmaps show the distribution of transitions across chromosomes, with color intensity indicating the percentage of regions in each category. 