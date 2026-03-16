import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt


def build_density(df, genome_length, bin_size):
    counts = np.zeros(genome_length, dtype=int)
    for _, row in df.iterrows():
        ref_pos = row.get("ref_pos")
        if pd.isna(ref_pos):
            continue
        ref_pos = int(ref_pos)
        vtype = row.get("type")
        length = row.get("length")
        if pd.isna(length):
            length = 1
        else:
            length = int(length)
        if vtype == "DEL":
            start = max(ref_pos, 0)
            end = min(ref_pos + length, genome_length)
            if start < end:
                counts[start:end] += 1
        else:
            if 0 <= ref_pos < genome_length:
                counts[ref_pos] += 1

    bins = np.arange(0, genome_length + bin_size, bin_size)
    binned = np.add.reduceat(counts, bins[:-1])
    centers = bins[:-1] + bin_size / 2
    return centers, binned


def main():
    parser = argparse.ArgumentParser(description="Summarize mutation CSV and generate plots")
    parser.add_argument("--csv", required=True, help="mutation_report.csv")
    parser.add_argument("--reference", required=True, help="reference.fasta")
    parser.add_argument("--gene-out", required=True, help="Output gene summary CSV")
    parser.add_argument("--density-plot", required=True, help="Output mutation density plot")
    parser.add_argument("--gene-plot", required=True, help="Output gene bar chart")
    parser.add_argument("--bin-size", type=int, default=500, help="Bin size for density plot")
    parser.add_argument("--top-genes", type=int, default=12, help="Number of genes in bar chart")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    df["gene"] = df["gene"].fillna("intergenic")

    # Gene-level summary
    summary = (
        df.groupby(["gene", "type"])\
        .size()\
        .unstack(fill_value=0)\
        .reset_index()
    )
    summary["total"] = summary.drop(columns=["gene"]).sum(axis=1)
    summary = summary.sort_values("total", ascending=False)
    summary.to_csv(args.gene_out, index=False)

    # Mutation density plot
    ref_len = len(SeqIO.read(args.reference, "fasta").seq)
    centers, binned = build_density(df, ref_len, args.bin_size)

    plt.figure(figsize=(10, 3.5))
    plt.plot(centers, binned, color="#0b6aa2", linewidth=1.5)
    plt.fill_between(centers, binned, color="#0b6aa2", alpha=0.15)
    plt.title("Mutation Density Across Genome")
    plt.xlabel("Genome position (bp)")
    plt.ylabel(f"Mutations per {args.bin_size} bp")
    plt.tight_layout()
    plt.savefig(args.density_plot, dpi=200)
    plt.close()

    # Gene bar chart
    top = summary.head(args.top_genes)
    plt.figure(figsize=(7, 4))
    plt.barh(top["gene"], top["total"], color="#2a9d8f")
    plt.gca().invert_yaxis()
    plt.title("Top Genes by Mutation Count")
    plt.xlabel("Mutation count")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig(args.gene_plot, dpi=200)
    plt.close()


if __name__ == "__main__":
    main()
