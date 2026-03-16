import argparse
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

import mutation_annotator as ma


def qc_sequence(seq, min_length, max_n_fraction):
    seq = ma.normalize_seq(seq)
    n_fraction = seq.count("N") / len(seq) if seq else 1.0
    if len(seq) < min_length:
        return False, n_fraction
    if n_fraction > max_n_fraction:
        return False, n_fraction
    return True, n_fraction


def summarize_rows(rows):
    summary = defaultdict(int)
    for row in rows:
        summary[row["type"]] += 1
        effect = row.get("effect")
        if effect:
            summary[effect] += 1
    frameshift = sum(1 for row in rows if row.get("frameshift"))
    summary["frameshift"] = frameshift
    return summary


def main():
    parser = argparse.ArgumentParser(description="Batch-run mutation annotation on a cohort FASTA")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--genbank", required=True)
    parser.add_argument("--cohort", required=True, help="Multi-FASTA of query genomes")
    parser.add_argument("--out-variants", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-qc", required=True)
    parser.add_argument("--min-length", type=int, default=29000)
    parser.add_argument("--max-n", type=float, default=0.01)
    parser.add_argument("--max-samples", type=int, default=0, help="Limit number of sequences (0 = all)")
    args = parser.parse_args()

    ref_record = SeqIO.read(args.reference, "fasta")
    ref_seq = ma.normalize_seq(ref_record.seq)

    gb_record = SeqIO.read(args.genbank, "genbank")
    cds_info, pos_to_cds = ma.build_cds_maps(gb_record)

    variants_rows = []
    summary_rows = []
    qc_rows = []

    processed = 0
    for record in SeqIO.parse(args.cohort, "fasta"):
        sample_id = record.id
        query_seq = ma.normalize_seq(record.seq)
        passed, n_fraction = qc_sequence(query_seq, args.min_length, args.max_n)

        if not passed:
            qc_rows.append(
                {
                    "sample_id": sample_id,
                    "length": len(query_seq),
                    "n_fraction": n_fraction,
                    "status": "filtered",
                }
            )
            continue

        aligned = ma.align_parasail(ref_seq, query_seq)
        if aligned is None:
            raise SystemExit("parasail alignment not available")

        aligned_ref, aligned_query = aligned
        warnings = []
        ma.validate_alignment(aligned_ref, aligned_query, warnings)

        ref_to_query, events = ma.build_ref_to_query_map(aligned_ref, aligned_query)
        rows = ma.annotate_variants(events, cds_info, pos_to_cds)
        ma.validate_query_cds(cds_info, ref_to_query, warnings)

        for row in rows:
            row["sample_id"] = sample_id
            variants_rows.append(row)

        summary = summarize_rows(rows)
        summary_rows.append(
            {
                "sample_id": sample_id,
                "length": len(query_seq),
                "n_fraction": n_fraction,
                "snp": summary.get("SNP", 0),
                "ins": summary.get("INS", 0),
                "del": summary.get("DEL", 0),
                "synonymous": summary.get("synonymous", 0),
                "nonsynonymous": summary.get("nonsynonymous", 0),
                "stop_gained": summary.get("stop_gained", 0),
                "stop_lost": summary.get("stop_lost", 0),
                "frameshift": summary.get("frameshift", 0),
                "warnings": "; ".join(warnings) if warnings else "",
            }
        )
        qc_rows.append(
            {
                "sample_id": sample_id,
                "length": len(query_seq),
                "n_fraction": n_fraction,
                "status": "kept",
            }
        )
        processed += 1
        if args.max_samples and processed >= args.max_samples:
            break

    if variants_rows:
        pd.DataFrame(variants_rows).to_csv(args.out_variants, index=False)
    pd.DataFrame(summary_rows).to_csv(args.out_summary, index=False)
    pd.DataFrame(qc_rows).to_csv(args.out_qc, index=False)


if __name__ == "__main__":
    main()
