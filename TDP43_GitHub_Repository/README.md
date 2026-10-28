# Donor-aware TDP-43 neuronal-nuclei RNA-seq reanalysis

Reproducible analysis code and exact inputs for a donor-aware secondary analysis of **GSE126542** (BioProject **PRJNA522295**).

The cohort comprises 12 frontal-cortex, FACS-isolated NeuN-positive neuronal-nuclei libraries from seven donors, including five complete matched pairs. Source “TDP-43 positive” means nuclear TDP-43 retained; “TDP-43 negative” means nuclear TDP-43 loss. Every contrast is **loss minus retained**.

## What was and was not run

- The reported differential-expression results were produced from the supplied gene-level count matrix using **PyDESeq2 0.5.4**.
- `code/02_run_differential_expression.py` implements input validation and the three reported models.
- `code/03_annotate_results.py` annotates results through MyGene.info and preserves raw responses and unmapped identifiers.
- `code/04_generate_figures.py` reproduces eight figures and their derived tables.
- `code/01_bulk_rnaseq_pipeline.sh` is a configurable historical/reference FASTQ-to-counts workflow. **It was not executed for the reported reanalysis.**

No cryptic-exon inclusion, PSI, isoform, or protein-level measurements were generated. Whole-gene UNC13A expression is not a surrogate for UNC13A cryptic splicing or functional protein.

## Repository structure

```text
.
├── code/        # FASTQ workflow, models, annotation, and figures
├── data/        # Exact count matrix and featureCounts summary
├── metadata/    # Verified analyzed-sample and full-series metadata
├── tables/      # Complete results and derived main/supplementary tables
├── references/  # Verified bibliographic records
├── requirements.txt
├── environment.yml
└── CITATION.cff
```

## Installation

Python 3.12 is recommended.

```bash
python3 -m pip install -r requirements.txt
```

Or use the supplied `environment.yml` with a Conda-compatible environment manager.

## 1. Validate the exact inputs

Run from the repository root:

```bash
python3 code/02_run_differential_expression.py \
  --counts data/counts_matrix.txt \
  --summary data/raw_counts.txt.summary \
  --metadata metadata/metadata_verified.csv \
  --outdir regenerated/validation \
  --validate-only
```

Expected validation landmarks:

- 62,754 genes and 12 libraries
- 10,564 all-zero genes
- 44,225 genes with at least 10 total counts
- seven donors and five complete matched pairs
- exact matrix-sum agreement with featureCounts Assigned counts
- assigned-read range 18,584,661–29,152,965; median 24,343,171

## 2. Rerun differential expression

```bash
python3 code/02_run_differential_expression.py \
  --counts data/counts_matrix.txt \
  --summary data/raw_counts.txt.summary \
  --metadata metadata/metadata_verified.csv \
  --outdir regenerated/models \
  --n-cpus 4
```

Models:

- primary: `~ donor + nuclear_status`
- sensitivity: `~ nuclear_status`
- sensitivity: `~ diagnosis + nuclear_status`

The script uses Cook’s-distance refitting/filtering, independent filtering, and Benjamini–Hochberg-adjusted P values. A fresh run writes TSV results; the supplied historical results are CSV. The figure script accepts either format.

## 3. Annotate results

Network access is required for new MyGene.info queries.

```bash
python3 code/03_annotate_results.py \
  --results tables/donor_aware_results.csv \
  --outdir regenerated/annotation \
  --top 50
```

Annotation is used for labeling, not statistical inference. External annotation may change over time.

## 4. Regenerate figures

To reproduce the supplied figures from the supplied tables:

```bash
python3 code/04_generate_figures.py --package-dir .
```

Outputs are written to `figures/` as vector PDF, 400-dpi PNG, and 400-dpi TIFF. Derived tables are written to `tables/`.

Each figure also has a separate runnable script under `code/figure_scripts/`, for example:

```bash
python3 code/figure_scripts/Figure_4_marker_expression.py --package-dir .
```

See `code/figure_scripts/README.md` for all eight entry points.

For a fully regenerated workflow, place the fresh model results and annotation outputs in a separate package copy with the same directory layout, then run the figure script against that copy. Preserve the original results for comparison.

## FASTQ-to-counts workflow

`code/01_bulk_rnaseq_pipeline.sh` requires SRA Toolkit, FastQC, Trimmomatic, STAR, samtools, featureCounts, a reference FASTA, GTF, adapter FASTA, and a verified featureCounts strandedness setting. Configure all paths and record every tool version before a new run. Do not claim that its default parameters generated the supplied matrix.

## Key results for consistency checks

- Donor-aware model: 3,442 genes at FDR < 0.05; 1,865 lower and 1,577 higher in loss nuclei.
- STMN2: LFC −1.812666; FDR 6.223865 × 10⁻¹⁵; approximately 71.5% lower.
- UNC13A: LFC +0.372086; FDR 0.006220552; approximately 29.4% higher whole-gene expression in the donor-aware model.
- TARDBP: LFC −0.606290; FDR 4.979965 × 10⁻⁶.
- Donor-aware versus diagnosis-adjusted effect sizes: Pearson r = 0.831; Spearman ρ = 0.829.

## Input SHA-256 checksums

- `counts_matrix.txt`: `b21f38d7ceb9be1c290b8e809af18be2fee4cc28919f1af6064c0eeb3b1decdc`
- `raw_counts.txt.summary`: `5c1f87eb3591305313c800ef792f695a65070909629e50060687abc5b27e9649`
- `metadata_verified.csv`: `801eb454979465121d507e2fd270b4d3f33f019c491d7a186ebeddb7e3caa596`

## Limitations

Small cohort; only five complete pairs; two eligible sorted libraries absent from the supplied matrix for undocumented reasons; no healthy controls; postmortem frontal-cortex neuronal nuclei only; source metadata annotate a C9orf72 repeat-expansion background; original FASTQ processing not rerun; gene-level rather than splicing-level evidence; no independent external validation.

## Attribution

Please cite both the accompanying reanalysis and the original Liu et al. GSE126542 study. The underlying project and analysis ownership remain with the authors. No GitHub upload is claimed by this local repository package.
