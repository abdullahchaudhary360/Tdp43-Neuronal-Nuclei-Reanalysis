# Table files

All tab-delimited tables are editable and use plain text or numeric fields.

## Main tables

- `Table_1_cohort_structure.tsv` — donor-level cohort structure and library accessions.
- `Table_2_model_comparison.tsv` — global differential-expression summaries for all models.
- `Table_3_marker_effects.tsv` — primary marker estimates, intervals, and sensitivity-model FDRs.

## Complete model and annotation outputs

- `donor_aware_results.csv`, `diagnosis_adjusted_results.csv`, `status_only_results.csv` — complete filtered model results.
- `donor_aware_results_annotated.tsv` — donor-aware results merged with MyGene.info annotation.
- `mygene_raw_annotation.tsv` — preserved raw annotation responses.
- `unmapped_gene_ids.txt` — identifiers without a returned gene symbol.

## Supplementary tables

- **Table S1** — `Table_S1_sample_metadata.tsv`: Verified sample-level metadata for the 12 analyzed libraries.
- **Table S2** — `Table_S2_donor_aware_significant_genes.tsv`: Annotated donor-aware genes with FDR < 0.05.
- **Table S3** — `Table_S3_top_genes.tsv`: Fifty leading donor-aware genes by adjusted P value.
- **Table S4** — `Table_S4_marker_model_statistics.tsv`: Marker estimates, standard errors, Wald intervals, and FDRs for all three models.
- **Table S5** — `Table_S5_model_comparison.tsv`: Global and marker-level summaries across models.
- **Table S6** — `Table_S6_PCA_coordinates.tsv`: PCA coordinates and sample annotations.
- **Table S7** — `Table_S7_effect_size_concordance.tsv`: Gene-wise donor-aware and diagnosis-adjusted effect sizes.
- **Table S8** — `Table_S8_marker_normalized_counts.tsv`: Normalized STMN2, UNC13A, and TARDBP counts by library.
- **Table S9** — `Table_S9_heatmap_gene_selection.tsv`: Selection statistics for the 30 heatmap genes.
- **Table S10** — `Table_S10_size_factors.tsv`: Descriptive median-of-ratios size factors.
