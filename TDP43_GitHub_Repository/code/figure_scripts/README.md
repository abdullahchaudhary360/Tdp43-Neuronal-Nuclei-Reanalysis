# Single-figure reproduction scripts

Each main or supplementary figure has a separate runnable entry point. All scripts reuse the verified functions in `code/04_generate_figures.py`; they do not duplicate or alter the statistical analysis.

Run from the package root:

```bash
python3 code/figure_scripts/Figure_1_study_design.py --package-dir .
python3 code/figure_scripts/Figure_2_QC_and_PCA.py --package-dir .
python3 code/figure_scripts/Figure_3_differential_expression.py --package-dir .
python3 code/figure_scripts/Figure_4_marker_expression.py --package-dir .
python3 code/figure_scripts/Figure_5_leading_gene_heatmap.py --package-dir .
python3 code/figure_scripts/Figure_S1_featureCounts_assignment.py --package-dir .
python3 code/figure_scripts/Figure_S2_sensitivity_volcanoes.py --package-dir .
python3 code/figure_scripts/Figure_S3_cryptic_exon_evidence_boundary.py --package-dir .
```

Each command writes the corresponding vector PDF, 400-dpi PNG, and 400-dpi TIFF to `figures/`. Scripts for Figures 2–5 also refresh their derived table files. `_common.py` is a shared loader, not a figure deliverable.

To regenerate all eight figures in one command:

```bash
python3 code/04_generate_figures.py --package-dir .
```

No script estimates cryptic-exon inclusion, PSI, isoform abundance, or protein-level consequences.
