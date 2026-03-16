# Viral Mutation Tracker (SARS-CoV-2)

Python pipeline to compare SARS-CoV-2 genomes against the NCBI reference (NC_045512.2), call SNPs/indels from global alignment, and annotate codon-aware amino-acid effects using GenBank CDS features. Includes cohort batch processing and figure generation (mutation density, gene-level mutation burden, and top recurrent AA changes).

## Setup

```bash
pip install -r requirements.txt
```

## Data

Reference genome (FASTA): `NC_045512.2`  
Reference annotations (GenBank): `NC_045512.2` (GenBank format)

You can download from NCBI manually, or via NCBI E-utilities (example URLs):

- FASTA: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text`
- GenBank: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=gbwithparts&retmode=text`

## Single Genome Comparison

```bash
python mutation_annotator.py --reference reference.fasta --query query.fasta --genbank reference.gb --out-csv mutation_report.csv
```

Outputs a per-variant CSV with nucleotide calls plus codon/AA effects when variants fall in CDS features.

## Cohort Run (Multi-FASTA)

Fetch a cohort FASTA from NCBI:

```bash
python fetch_ncbi_cohort.py --query "SARS-CoV-2 Omicron[All Fields] AND complete genome[Title]" --retmax 50 --out data/cohort.fasta
```

Run batch annotation:

```bash
python batch_run.py --reference reference.fasta --genbank reference.gb --cohort data/cohort.fasta --out-variants reports/cohort_variants.csv --out-summary reports/cohort_summary.csv --out-qc reports/cohort_qc.csv --max-samples 200
```

## Manuscript-Style Figure Panel

```bash
python make_results_panel.py --variants reports/cohort_variants.csv --reference reference.fasta --out figures/results_panel.png
```

## Notes

- This is a research/education pipeline (not for clinical use).
- If your environment has SSL certificate issues reaching NCBI, `fetch_ncbi_cohort.py` supports `--insecure` to disable SSL verification.
