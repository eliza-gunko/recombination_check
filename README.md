# Recombination check

This script processes a tab-delimited genotype file.
It assumes the input file is prepared with columns:
- chr, position, parent1, parent2, sample1, sample2, ... etc.

The script performs the following steps:
1. Reads the input file.
2. Calculates chi-square p-values for genotype frequencies (0/0, 0/1, 1/1) using only the genotype columns (when header starts with "sample").
3. Filters rows with chi-square p-values > 0.05. This step helps to get rid of the artefacts (sequencing and calling mistakes)
4. Converts genotype strings to numeric values:
    - '0/0' → -1, '0/1' or '1/0' → 0, '1/1' → 1.
5. Applies a rolling median (window = window_size, centered) on each genotype column grouped by chromosome.
6. Extracts recombination coordinate changes when genotype values change.
7. Writes two output tables:
    - recombination_coordinates.tsv (all detected coordinates), neighboring rows are dedicated to the same 'switching event'
    - filtered_recombination_coordinates.tsv (only those with genotype = -0.5 or 0.5), because these coordinates are supposed to be closer to the place of the actual genotype switching, then the other rows in unfiltered recombination_coordinates.tsv
8. Generates plots:
    - A separate plot per chromosome (recombination coordinates; blue dots) saved in the "chromosome_plots" directory. 
        Here the samples will be plotted together (one plot - many samples)
    - Genotype-median plots for one or more user-specified chromosomes. For each chromosome, a separate output directory (e.g. "genotype_plots_chromosome_1") is created and the corresponding plots are saved there. Here the samples are analysed one by one (one plot - one sample)
            
Usage:
```bash
python recombination_check.py --input input_file.txt --output_dir output_directory --window window_size [--plot_chromosomes "chromosome_1,chromosome_2,..."]
```