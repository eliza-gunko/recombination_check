import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import chisquare
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

tqdm.pandas()

def calc_chisq(x):
    counts = [
        np.sum(x == '0/0'),
        np.sum((x == '0/1') | (x == '1/0')),
        np.sum(x == '1/1')
    ]
    total_counts = sum(counts)
    if total_counts == 0:
        return np.nan
    expected_counts = [0.25 * total_counts, 0.5 * total_counts, 0.25 * total_counts]
    _, p_value = chisquare(counts, f_exp=expected_counts)
    return p_value

def geno_to_number(x):
    x = np.where(x == '0/0', -1, x)
    x = np.where((x == '0/1') | (x == '1/0'), 0, x)
    x = np.where(x == '1/1', 1, x)
    return np.where(pd.isna(x), np.nan, x).astype(np.float64)

def plot_genotype_median(df, out_dir):
    target_columns = [col for col in df.columns if col.lower().startswith('sample')]
    discrete_y_values = [-1, -0.5, 0, 0.5, 1]
    
    for column in target_columns:
        genotype_median = df[column]
        fig, ax = plt.subplots(figsize=(14, 8))
        
        ax.scatter(df['position'] / 1e6, genotype_median, marker='o', color='b', s=2)
        ax.set_ylim(-1.5, 1.5)
        ax.set_yticks(discrete_y_values)
        ax.set_xlabel('Position (Mb)', fontsize=20)
        ax.set_ylabel('Genotype', fontsize=20)
        ax.set_title(f'{column}', fontsize=25)
        ax.grid(True, linestyle='--', alpha=0.5)
        
        formatter = mticker.ScalarFormatter(useOffset=False, useMathText=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        
        output_path = os.path.join(out_dir, f"{column}_median_plot.png")
        fig.savefig(output_path, facecolor='white', dpi=600)
        plt.close(fig)

def main(args):
    window_size = args.window

    try:
        df = pd.read_csv(args.input, sep="\t")
    except Exception as e:
        sys.exit(f"Error reading input file: {e}")

    print("Columns in the input file:")
    print(list(df.columns))
    
    df = df.replace('./.', np.nan)

    try:
        df = df.sort_values(by=["chr", "position"]).reset_index(drop=True)
    except Exception as e:
        sys.exit(f"Error sorting DataFrame: {e}")

    genotype_columns = [col for col in df.columns if col.lower().startswith('sample')]
    parent_columns = [col for col in df.columns if col.lower().startswith('parent')]

    if not genotype_columns:
        sys.exit("No genotype columns (starting with 'sample') found. Please check your input table.")

    df['chisqp'] = df.progress_apply(lambda x: calc_chisq(x[genotype_columns].values.flatten()), axis=1)
    
    df_filtered = df[df['chisqp'] > 0.05].copy()
    print(f"Rows after chi-square filter: {df_filtered.shape[0]}")

    for col in genotype_columns:
        df_filtered[col] = geno_to_number(df_filtered[col])
    for col in parent_columns:
        df_filtered[col] = geno_to_number(df_filtered[col])

    sample_columns_chr_pos = ["chr", "position"] + genotype_columns
    samples_df = df_filtered[sample_columns_chr_pos].copy()

    for col in genotype_columns:
        samples_df[col] = samples_df.groupby("chr")[col].apply(lambda x: x.rolling(window_size, center=True).median())
    
    nan_columns = samples_df.columns.difference(["chr", "position"])
    cleaned_df = samples_df.dropna(subset=nan_columns, how='all').reset_index(drop=True)

    list_of_samples = []
    for chr_name, chr_data in tqdm(cleaned_df.groupby("chr"), desc="Processing chromosomes"):
        for col in tqdm(chr_data.columns[2:], desc=f"Processing samples in {chr_name}", leave=False):
            list_of_positions = []
            list_of_genotypes = []
            for i in range(1, len(chr_data)):
                if chr_data[col].iloc[i] != chr_data[col].iloc[i - 1]:
                    list_of_positions.append(chr_data["position"].iloc[i - 1])
                    list_of_genotypes.append(chr_data[col].iloc[i - 1])
                    list_of_positions.append(chr_data["position"].iloc[i])
                    list_of_genotypes.append(chr_data[col].iloc[i])
            list_of_samples.append((chr_name, col, list_of_positions, list_of_genotypes))
    
    rows = []
    for chr_name, col, positions, genotypes in tqdm(list_of_samples, desc="Compiling rows"):
        for pos, geno in zip(positions, genotypes):
            rows.append((chr_name, col, pos, geno))
    
    recomb_df = pd.DataFrame(rows, columns=['chr', 'Sample', 'position', 'Genotype'])
    recomb_df = recomb_df.drop_duplicates()
    recomb_df['Genotype'] = recomb_df['Genotype'].round(1)
    
    recomb_df_filtered = recomb_df[(recomb_df['Genotype'] == -0.5) | (recomb_df['Genotype'] == 0.5)]
    
    os.makedirs(args.output_dir, exist_ok=True)
    recomb_out = os.path.join(args.output_dir, "recombination_coordinates.tsv")
    filtered_recomb_out = os.path.join(args.output_dir, "filtered_recombination_coordinates.tsv")
    
    recomb_df.to_csv(recomb_out, sep="\t", index=False)
    recomb_df_filtered.to_csv(filtered_recomb_out, sep="\t", index=False)
    
    print(f"Recombination coordinates table saved to: {recomb_out}")
    print(f"Filtered recombination coordinates table saved to: {filtered_recomb_out}")

    chromosome_plots_dir = os.path.join(args.output_dir, "chromosome_plots")
    os.makedirs(chromosome_plots_dir, exist_ok=True)

    for chr_name in tqdm(recomb_df['chr'].unique(), desc="Generating per-chromosome recombination plots"):
        chr_data = recomb_df[recomb_df['chr'] == chr_name]
        unique_samples = sorted(chr_data['Sample'].unique())
        sample_y_map = {sample: i for i, sample in enumerate(unique_samples)}
        y_positions = chr_data['Sample'].map(sample_y_map)
        
        plt.figure(figsize=(20, 5))
        min_pos = recomb_df['position'].min() / 1e6
        max_pos = recomb_df['position'].max() / 1e6
        plt.xlim(min_pos - 1, max_pos + 1)
        
        plt.scatter(
            chr_data['position'] / 1e6,
            y_positions,
            color='blue',
            s=20,
            alpha=0.7
        )
        plt.title(f"{chr_name}", fontsize=20)
        plt.xlabel("Position (Mb)", fontsize=18)
        plt.ylabel("Summary for individual samples", fontsize=18)
        plt.yticks([])
        plt.grid(axis='x', linestyle='--', alpha=0.5)
        
        chr_plot_file = os.path.join(chromosome_plots_dir, f"{chr_name}_blue_dots.png")
        plt.savefig(chr_plot_file, dpi=600, facecolor='white')
        plt.close()

    plt.figure(figsize=(20, 5))
    min_pos = recomb_df['position'].min() / 1e6
    max_pos = recomb_df['position'].max() / 1e6
    plt.xlim(min_pos - 1, max_pos + 1)
    
    for chr_name in tqdm(recomb_df['chr'].unique(), desc="Plotting all chromosomes combined"):
        chr_data = recomb_df[recomb_df['chr'] == chr_name]
        unique_samples = sorted(chr_data['Sample'].unique())
        sample_y_map = {sample: i for i, sample in enumerate(unique_samples)}
        y_positions = chr_data['Sample'].map(sample_y_map)
        plt.scatter(
            chr_data['position'] / 1e6,
            y_positions,
            color='blue',
            s=20,
            alpha=0.7
        )
    
    plt.title("Recombination Positions Across All Chromosomes", fontsize=20)
    plt.xlabel("Position (Mb)", fontsize=18)
    plt.ylabel("Summary for individual samples", fontsize=18)
    plt.yticks([])
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    
    combined_plot_file = os.path.join(chromosome_plots_dir, "all_chromosomes_combined.png")
    plt.savefig(combined_plot_file, dpi=600, facecolor='white')
    plt.close()
    
    print(f"Recombination plots saved to: {chromosome_plots_dir}")

    if args.plot_chromosomes:
        chromosomes_to_plot = [chrom.strip() for chrom in args.plot_chromosomes.split(",")]
        genotype_df = df_filtered.copy()
        genotype_columns = [col for col in genotype_df.columns if col.lower().startswith('sample')]
        for col in genotype_columns:
            genotype_df[col] = geno_to_number(genotype_df[col])
            genotype_df[col] = genotype_df.groupby("chr")[col].apply(lambda x: x.rolling(window_size, center=True).median())
        
        for chrom in chromosomes_to_plot:
            df_chrom = genotype_df[genotype_df["chr"] == chrom].copy()
            if df_chrom.empty:
                print(f"No data for chromosome: {chrom}")
                continue
            df_chrom = df_chrom.sort_values(by="position").reset_index(drop=True)
            chrom_genotype_dir = os.path.join(args.output_dir, f"genotype_plots_{chrom}")
            os.makedirs(chrom_genotype_dir, exist_ok=True)
            plot_genotype_median(df_chrom, chrom_genotype_dir)
            print(f"Genotype median plots saved to: {chrom_genotype_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process a genotype table (with columns chr, position, parent1, parent2, sample1, sample2, ...) and create recombination coordinate tables and plots."
    )
    parser.add_argument("--input", required=True, help="Path to the input tab-delimited genotype file.")
    parser.add_argument("--output_dir", required=True, help="Directory where the output tables and plots will be saved.")
    parser.add_argument("--window", type=int, required=True, help="Window size for the rolling median.")
    parser.add_argument("--plot_chromosomes", required=False,
                        help="Comma-separated list of chromosomes for genotype median plots (e.g., 'chromosome_1,chromosome_4').",
                        default="")
    args = parser.parse_args()
    main(args)
