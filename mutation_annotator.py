import argparse
import csv
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

VALID_BASES = set("ACGTN")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
CODON_TABLE = CodonTable.unambiguous_dna_by_id[1]


def normalize_seq(seq):
    return str(seq).upper().replace("U", "T")


def complement_base(base):
    return base.translate(COMPLEMENT)


def safe_translate_codon(codon):
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "X"
    if codon in CODON_TABLE.stop_codons:
        return "*"
    return CODON_TABLE.forward_table.get(codon, "X")


def load_fasta(path):
    return SeqIO.read(path, "fasta")


def load_genbank(path):
    return SeqIO.read(path, "genbank")


def get_feature_label(feature):
    qualifiers = feature.qualifiers
    gene = (qualifiers.get("gene") or [None])[0]
    locus = (qualifiers.get("locus_tag") or [None])[0]
    product = (qualifiers.get("product") or [None])[0]
    protein_id = (qualifiers.get("protein_id") or [None])[0]
    label = gene or locus or protein_id or "CDS"
    return label, product, protein_id


def iter_coding_positions(feature):
    location = feature.location
    parts = location.parts if hasattr(location, "parts") else [location]
    strand = location.strand or 1
    if strand == 1:
        for part in parts:
            for pos in range(int(part.start), int(part.end)):
                yield pos
    else:
        for part in reversed(parts):
            for pos in range(int(part.end) - 1, int(part.start) - 1, -1):
                yield pos


def build_cds_maps(ref_record):
    cds_info = {}
    pos_to_cds = defaultdict(list)

    for feature in ref_record.features:
        if feature.type != "CDS":
            continue
        label, product, protein_id = get_feature_label(feature)
        cds_id = f"{label}|{protein_id or 'NA'}"
        coding_positions = list(iter_coding_positions(feature))
        cds_seq = normalize_seq(feature.extract(ref_record.seq))
        strand = feature.location.strand or 1

        cds_info[cds_id] = {
            "label": label,
            "product": product,
            "protein_id": protein_id,
            "strand": strand,
            "seq": cds_seq,
            "coding_positions": coding_positions,
        }

        for idx, pos in enumerate(coding_positions):
            pos_to_cds[int(pos)].append((cds_id, idx))

    return cds_info, pos_to_cds


def validate_sequence(seq, name, warnings, n_threshold=0.05):
    seq = normalize_seq(seq)
    invalid = sorted(set(seq) - VALID_BASES)
    if invalid:
        warnings.append(f"{name}: contains invalid bases: {''.join(invalid)}")
    n_count = seq.count("N")
    if len(seq) > 0 and (n_count / len(seq)) > n_threshold:
        warnings.append(
            f"{name}: N content is {n_count/len(seq):.2%} (> {n_threshold:.0%})"
        )


def validate_cds(cds_info, warnings):
    for cds_id, cds in cds_info.items():
        if len(cds["seq"]) % 3 != 0:
            warnings.append(
                f"CDS {cds['label']} length {len(cds['seq'])} not multiple of 3"
            )


def validate_alignment(aligned_ref, aligned_query, warnings, gap_threshold=0.1):
    if len(aligned_ref) != len(aligned_query):
        warnings.append("Alignment length mismatch between reference and query")
        return
    total = len(aligned_ref)
    if total == 0:
        warnings.append("Empty alignment")
        return
    gaps = aligned_ref.count("-") + aligned_query.count("-")
    gap_ratio = gaps / (2 * total)
    if gap_ratio > gap_threshold:
        warnings.append(f"High gap ratio in alignment: {gap_ratio:.2%}")


def align_pairwise(ref_seq, query_seq):
    # Slow for long genomes; intended for short segments or testing only.
    from Bio import pairwise2

    alignments = pairwise2.align.globalms(ref_seq, query_seq, 2, -1, -5, -1)
    best = alignments[0]
    return best[0], best[1]


def align_parasail(ref_seq, query_seq, match=2, mismatch=-1, gap_open=5, gap_extend=1):
    try:
        import parasail
    except Exception:
        return None

    def sanitize(seq):
        return "".join(b if b in "ACGTN" else "N" for b in seq.upper().replace("U", "T"))

    ref_seq = sanitize(ref_seq)
    query_seq = sanitize(query_seq)

    matrix = parasail.matrix_create("ACGTN", match, mismatch)
    result = parasail.nw_trace_scan_32(ref_seq, query_seq, gap_open, gap_extend, matrix)
    tb = result.traceback
    return tb.query, tb.ref


def load_aligned_fasta(path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 2:
        raise ValueError("Aligned FASTA must contain exactly two sequences")
    return str(records[0].seq), str(records[1].seq)


def build_ref_to_query_map(aligned_ref, aligned_query):
    ref_pos = -1
    query_pos = -1
    ref_to_query = {}
    events = []

    i = 0
    while i < len(aligned_ref):
        r = aligned_ref[i]
        q = aligned_query[i]

        if r != "-" and q != "-":
            ref_pos += 1
            query_pos += 1
            ref_to_query[ref_pos] = q
            if r != q:
                events.append(
                    {
                        "type": "SNP",
                        "ref_pos": ref_pos,
                        "ref_base": r,
                        "alt_base": q,
                        "length": 1,
                    }
                )
            i += 1
            continue

        if r != "-" and q == "-":
            start_pos = ref_pos + 1
            del_bases = []
            while i < len(aligned_ref) and aligned_ref[i] != "-" and aligned_query[i] == "-":
                ref_pos += 1
                del_bases.append(aligned_ref[i])
                ref_to_query[ref_pos] = "-"
                i += 1
            events.append(
                {
                    "type": "DEL",
                    "ref_pos": start_pos,
                    "ref_base": "".join(del_bases),
                    "alt_base": "-",
                    "length": len(del_bases),
                }
            )
            continue

        if r == "-" and q != "-":
            ins_bases = []
            while i < len(aligned_ref) and aligned_ref[i] == "-" and aligned_query[i] != "-":
                query_pos += 1
                ins_bases.append(aligned_query[i])
                i += 1
            events.append(
                {
                    "type": "INS",
                    "ref_pos": ref_pos,
                    "ref_base": "-",
                    "alt_base": "".join(ins_bases),
                    "length": len(ins_bases),
                }
            )
            continue

        i += 1

    return ref_to_query, events


def annotate_snp(event, cds_info, pos_to_cds):
    rows = []
    pos = event["ref_pos"]
    for cds_id, idx in pos_to_cds.get(pos, []):
        cds = cds_info[cds_id]
        codon_index = idx // 3
        codon_offset = idx % 3
        codon_start = codon_index * 3
        codon_end = codon_start + 3
        if codon_end > len(cds["seq"]):
            continue
        ref_codon = cds["seq"][codon_start:codon_end]
        alt_base = event["alt_base"].upper()
        if cds["strand"] == -1:
            alt_base = complement_base(alt_base)
        alt_codon_list = list(ref_codon)
        alt_codon_list[codon_offset] = alt_base
        alt_codon = "".join(alt_codon_list)
        ref_aa = safe_translate_codon(ref_codon)
        alt_aa = safe_translate_codon(alt_codon)

        if ref_aa == alt_aa:
            effect = "synonymous"
        elif alt_aa == "*":
            effect = "stop_gained"
        elif ref_aa == "*":
            effect = "stop_lost"
        else:
            effect = "nonsynonymous"

        rows.append(
            {
                "gene": cds["label"],
                "product": cds["product"],
                "protein_id": cds["protein_id"],
                "codon_index": codon_index + 1,
                "codon_pos": codon_offset + 1,
                "ref_codon": ref_codon,
                "alt_codon": alt_codon,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "aa_pos": codon_index + 1,
                "effect": effect,
            }
        )

    if not rows:
        rows.append({"gene": None})
    return rows


def annotate_indel(event, cds_info, pos_to_cds):
    rows = []
    if event["type"] == "DEL":
        start = event["ref_pos"]
        length = event["length"]
        affected = defaultdict(list)
        for offset in range(length):
            pos = start + offset
            for cds_id, idx in pos_to_cds.get(pos, []):
                affected[cds_id].append(idx)

        for cds_id, idxs in affected.items():
            cds = cds_info[cds_id]
            del_len = len(idxs)
            frameshift = del_len % 3 != 0
            codon_index = min(idxs) // 3 + 1
            effect = "frameshift_deletion" if frameshift else "inframe_deletion"
            rows.append(
                {
                    "gene": cds["label"],
                    "product": cds["product"],
                    "protein_id": cds["protein_id"],
                    "codon_index": codon_index,
                    "effect": effect,
                    "frameshift": frameshift,
                }
            )

        if not rows:
            rows.append({"gene": None})
        return rows

    if event["type"] == "INS":
        pos = event["ref_pos"]
        candidate = set()
        for cds_id, _ in pos_to_cds.get(pos, []):
            candidate.add(cds_id)
        for cds_id, _ in pos_to_cds.get(pos + 1, []):
            candidate.add(cds_id)

        for cds_id in candidate:
            cds = cds_info[cds_id]
            frameshift = event["length"] % 3 != 0
            effect = "frameshift_insertion" if frameshift else "inframe_insertion"
            rows.append(
                {
                    "gene": cds["label"],
                    "product": cds["product"],
                    "protein_id": cds["protein_id"],
                    "codon_index": None,
                    "effect": effect,
                    "frameshift": frameshift,
                }
            )

        if not rows:
            rows.append({"gene": None})
        return rows

    return rows


def build_query_cds_sequence(cds, ref_to_query):
    seq = []
    for pos in cds["coding_positions"]:
        base = ref_to_query.get(pos, "N")
        if base == "-":
            base = "N"
        if cds["strand"] == -1:
            base = complement_base(base)
        seq.append(base)
    return "".join(seq)


def validate_query_cds(cds_info, ref_to_query, warnings):
    for cds in cds_info.values():
        query_cds = build_query_cds_sequence(cds, ref_to_query)
        if len(query_cds) < 3:
            continue
        translation = str(Seq(query_cds).translate(to_stop=False))
        if "*" in translation[:-1]:
            warnings.append(
                f"CDS {cds['label']}: premature stop codon detected in query"
            )


def annotate_variants(events, cds_info, pos_to_cds):
    rows = []
    for event in events:
        if event["type"] == "SNP":
            annotations = annotate_snp(event, cds_info, pos_to_cds)
        else:
            annotations = annotate_indel(event, cds_info, pos_to_cds)
        for ann in annotations:
            row = {
                "type": event["type"],
                "ref_pos": event["ref_pos"],
                "ref_base": event["ref_base"],
                "alt_base": event["alt_base"],
                "length": event["length"],
                "gene": ann.get("gene"),
                "product": ann.get("product"),
                "protein_id": ann.get("protein_id"),
                "codon_index": ann.get("codon_index"),
                "codon_pos": ann.get("codon_pos"),
                "ref_codon": ann.get("ref_codon"),
                "alt_codon": ann.get("alt_codon"),
                "ref_aa": ann.get("ref_aa"),
                "alt_aa": ann.get("alt_aa"),
                "aa_pos": ann.get("aa_pos"),
                "effect": ann.get("effect"),
                "frameshift": ann.get("frameshift"),
            }
            rows.append(row)
    return rows


def write_csv(rows, path):
    if not rows:
        return
    fieldnames = [
        "type",
        "ref_pos",
        "ref_base",
        "alt_base",
        "length",
        "gene",
        "product",
        "protein_id",
        "codon_index",
        "codon_pos",
        "ref_codon",
        "alt_codon",
        "ref_aa",
        "alt_aa",
        "aa_pos",
        "effect",
        "frameshift",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_events(events, rows):
    summary = defaultdict(int)
    for event in events:
        summary[event["type"]] += 1
    frameshift = sum(1 for row in rows if row.get("frameshift"))
    return summary, frameshift


def main():
    parser = argparse.ArgumentParser(
        description="Codon-aware mutation annotation for aligned SARS-CoV-2 genomes"
    )
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--query", required=True, help="Query FASTA")
    parser.add_argument("--genbank", required=True, help="Reference GenBank for CDS")
    parser.add_argument("--aligned-fasta", help="Aligned FASTA with exactly two sequences")
    parser.add_argument(
        "--out-csv", default="mutation_report.csv", help="Output CSV for variants"
    )
    parser.add_argument(
        "--max-pairwise-length",
        type=int,
        default=5000,
        help="Max length for in-script alignment (otherwise require --aligned-fasta)",
    )

    args = parser.parse_args()

    ref_record = load_fasta(args.reference)
    query_record = load_fasta(args.query)
    ref_seq = normalize_seq(ref_record.seq)
    query_seq = normalize_seq(query_record.seq)

    warnings = []
    validate_sequence(ref_seq, "Reference", warnings)
    validate_sequence(query_seq, "Query", warnings)

    gb_record = load_genbank(args.genbank)
    cds_info, pos_to_cds = build_cds_maps(gb_record)
    validate_cds(cds_info, warnings)

    if args.aligned_fasta:
        aligned_ref, aligned_query = load_aligned_fasta(args.aligned_fasta)
    else:
        if max(len(ref_seq), len(query_seq)) > args.max_pairwise_length:
            aligned = align_parasail(ref_seq, query_seq)
            if aligned is None:
                raise SystemExit(
                    "Sequences too long for in-script alignment. Provide --aligned-fasta "
                    "from MAFFT/minimap2, or install parasail for fast alignment."
                )
            aligned_ref, aligned_query = aligned
        else:
            aligned_ref, aligned_query = align_pairwise(ref_seq, query_seq)

    validate_alignment(aligned_ref, aligned_query, warnings)

    ref_to_query, events = build_ref_to_query_map(aligned_ref, aligned_query)
    rows = annotate_variants(events, cds_info, pos_to_cds)
    validate_query_cds(cds_info, ref_to_query, warnings)
    write_csv(rows, args.out_csv)

    summary, frameshift = summarize_events(events, rows)

    print(f"Reference: {ref_record.id}")
    print(f"Query: {query_record.id}")
    print(f"Variants: SNP={summary['SNP']} INS={summary['INS']} DEL={summary['DEL']}")
    print(f"Frameshift annotations: {frameshift}")
    print(f"CSV report: {args.out_csv}")

    if warnings:
        print("\nValidation warnings:")
        for w in warnings:
            print(f"- {w}")


if __name__ == "__main__":
    main()
