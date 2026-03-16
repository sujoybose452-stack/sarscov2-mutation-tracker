import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO


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


def top_aa_mutations(df, top_n):
    aa = df.copy()
    aa = aa[aa["ref_aa"].notna() & aa["alt_aa"].notna()]
    aa = aa[aa["ref_aa"].astype(str).str.len() == 1]
    aa = aa[aa["alt_aa"].astype(str).str.len() == 1]
    aa = aa[~aa["ref_aa"].isin(["X", "*"])]
    aa = aa[~aa["alt_aa"].isin(["X", "*"])]
    aa = aa[aa["ref_aa"] != aa["alt_aa"]]
    aa = aa[aa["gene"].notna()]
    aa["aa_pos"] = aa["aa_pos"].fillna(0).astype(int)
    aa["label"] = (
        aa["gene"].astype(str)
        + ":"
        + aa["ref_aa"].astype(str)
        + aa["aa_pos"].astype(str)
        + aa["alt_aa"].astype(str)
    )
    counts = aa.groupby("label").size().sort_values(ascending=False)
    return counts.head(top_n)


def main():
    parser = argparse.ArgumentParser(description="Build a Results figure panel")
    parser.add_argument("--variants", required=True, help="cohort_variants.csv")
    parser.add_argument(
        "--gene-summary",
        help="cohort_gene_summary.csv (optional; computed from variants if omitted)",
    )
    parser.add_argument("--reference", required=True, help="reference.fasta")
    parser.add_argument("--out", required=True, help="Output PNG")
    parser.add_argument("--bin-size", type=int, default=500)
    parser.add_argument("--top-genes", type=int, default=10)
    parser.add_argument("--top-aa", type=int, default=10)
    args = parser.parse_args()

    variants = pd.read_csv(args.variants)
    if args.gene_summary:
        gene_summary = pd.read_csv(args.gene_summary)
    else:
        tmp = variants.copy()
        tmp["gene"] = tmp["gene"].fillna("intergenic")
        gene_summary = (
            tmp.groupby(["gene", "type"])
            .size()
            .unstack(fill_value=0)
            .reset_index()
        )
        gene_summary["total"] = gene_summary.drop(columns=["gene"]).sum(axis=1)

    ref_len = len(SeqIO.read(args.reference, "fasta").seq)
    centers, binned = build_density(variants, ref_len, args.bin_size)

    top_genes = gene_summary.sort_values("total", ascending=False).head(args.top_genes)
    top_aa = top_aa_mutations(variants, args.top_aa)

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), gridspec_kw={"height_ratios": [1.1, 1, 1.1]})

    # Panel A: mutation density
    ax = axes[0]
    ax.plot(centers, binned, color="#0b6aa2", linewidth=1.6)
    ax.fill_between(centers, binned, color="#0b6aa2", alpha=0.2)
    ax.set_title("A. Mutation Density Across Genome")
    ax.set_ylabel(f"Mutations / {args.bin_size} bp")
    ax.set_xlabel("Genome position (bp)")

    # Panel B: gene bar chart
    ax = axes[1]
    ax.barh(top_genes["gene"], top_genes["total"], color="#2a9d8f")
    ax.invert_yaxis()
    ax.set_title("B. Top Genes by Mutation Count")
    ax.set_xlabel("Mutation count")
    ax.set_ylabel("Gene")

    # Panel C: top AA mutations
    ax = axes[2]
    ax.barh(top_aa.index[::-1], top_aa.values[::-1], color="#e76f51")
    ax.set_title("C. Top Amino-Acid Changes")
    ax.set_xlabel("Occurrences in cohort")
    ax.set_ylabel("AA change")

    plt.tight_layout()
    plt.savefig(args.out, dpi=240)
    plt.close(fig)


if __name__ == "__main__":
    main()
