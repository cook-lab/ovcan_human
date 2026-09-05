# Integration diff — regenerated `output/` vs `reports/_snapshot_postfix/output/`

Generated 2026-07-24 16:38 EDT by `reports/integration_logs/compare_outputs.R`.

**Baseline:** `reports/_snapshot_postfix/output/` — 122 files, the post-repair state.
**Regenerated:** `output/` after the certified pass-2 clean run (23 script invocations,
all exit 0, from a genuinely empty `output/` with `metadata/line_family_map.csv` removed).

## Method

Every baseline file is matched by name. CSVs are parsed and compared on row count,
column count, column-name set, and cell-by-cell content — numerics with a relative
tolerance of 1e-8, everything else as text. Non-CSV files are compared byte-for-byte.
Columns and rows whose only purpose is to record *when* a run happened
(`recorded_utc`, `timestamp`, `*_utc`) are excluded from the verdict and reported
separately, because they change on every run by design.

## Summary

| Status | n | Meaning |
|---|---|---|
| OK — identical | 109 | Byte- or value-identical to the baseline |
| OK* — provenance only | 12 | `session_info*.txt`; differ in run timestamp / session introspection |
| CHANGED — explained | 1 | `wes_signature_refit_bootstrap.csv` (see below) |
| MISSING | 0 | — |
| **Total** | **122** | |

**No file differs in a numeric value, a column name, or a row count except the one
documented below. Nothing is missing.**

## The one changed file

`wes_signature_refit_bootstrap.csv` — 1,250 rows (baseline) -> **1,804 rows** (regenerated).

This is the bootstrap *cache*, not a result. 1,804 = 22 lines x 60 signatures (1,320,
`full` reference set) + 22 lines x 22 signatures (484, `restricted`). The regenerated
file is the complete grid; the baseline held only the ever-selected subset
(916 + 334), i.e. it was written by a superseded version of `boot_summary()`.

Evidence that this is **not** non-determinism:

* Pass 1 and pass 2 — two independent from-scratch 200-replicate bootstraps —
  produced **byte-identical** cache files (1,804 rows each).
* A third recompute, forced by planting the short baseline cache, again produced a
  **byte-identical** cache and byte-identical `wes_signature_refit_exposures.csv`
  and `wes_signature_refit_summary.csv`.
* Both derived result tables are IDENTICAL to the baseline, so no reported statistic
  depends on the extra rows.
* The cache-key fields agree across all three files:
  `n_boots` 200, `cache_max_delta` 0.004, `cache_n_ref_full` 60, `cache_n_ref_restricted` 22.

### A real latent defect this surfaced (fixed)

The cache-validity test at `22_wes_signature_refit.R:194-201` checked only cache
*keys*. The short 1,250-row baseline cache satisfied every one of them, so it would
have been silently reused — and the join at section 5 uses
`replace_na(boot_selected_frac, 0)` in its final `filter()`, so a **missing** cache row
is indistinguishable from a *never-selected* signature and the row is dropped from
`wes_signature_refit_exposures.csv` without a warning.

Two fixes:

* `22:202-213` — cache validity now also requires the **complete grid**
  (`n_lines x n_signatures` per reference set), so an incomplete cache is rejected.
* `22:258-265` — a `stopifnot` asserts the bootstrap summary covers every
  (line, signature, reference_set) triple, and has no duplicates, before the join.

Both paths were then tested end to end: with the short cache planted the script logs
`Cached bootstrap does not match the current parameters, recomputing` and rebuilds the
full cache; with the complete cache present it logs `Reusing cached bootstrap summary`
and reuses it.

### Footnote: cache reuse changes float *printing*, not values

`wes_signature_refit_exposures.csv` written via the reuse path is numerically identical
to the recompute path (`all.equal` TRUE, same row order) but not byte-identical: the CSV
round-trip drops a final significant digit on a few values, e.g. `96.72750068668061`
-> `96.7275006866806`. Harmless, and invisible to the tolerance-based comparison. The
delivered `output/` holds the recompute-path (pass-2) files.

## Provenance files that differ by design (12)

`00_setup.R` rewrites `output/session_info.txt` on every `source()`, and 11 scripts
write their own `session_info_<script>.txt`. These record *when* a script ran and what
its R session had loaded, so they cannot be identical across runs:

* `session_info.txt` — differs in `caller` and timestamp only. Baseline recorded
  `03_rna_de_signatures.R`; the regenerated file records the last script sourced.
* `session_info_01_rna_load_qc.txt` — timestamp, plus a different
  "loaded via a namespace (and not attached)" list. That list depends on which
  namespaces R lazy-loaded during the run, so it is session introspection rather than
  an analysis output.
* The remaining 10 — identical size, differ only in the `run_at` line.

`output/package_versions.csv` compared **IDENTICAL** once its `recorded_utc` row is
excluded, which is why it is the better single citation for the environment (as
`fix_report_wes.md` §4.1 argued).

## Full file table

| File | Status | rows (base -> new) | cols | Notes |
|---|---|---|---|---|
| `adc_expression.csv` | OK |     527 | 8 |  |
| `adc_folr1_bimodality.csv` | OK |      21 | 6 |  |
| `adc_modality_line_sets.csv` | OK |       2 | 3 |  |
| `adc_subtype_summary.csv` | OK |      17 | 15 |  |
| `auth_mucinous.csv` | OK |       3 | 22 |  |
| `auth_mucinous_marker_ranks.csv` | OK |     496 | 17 |  |
| `auth_mucinous_sensitivity.csv` | OK |     180 | 11 |  |
| `auth_perline_table.csv` | OK |      42 | 27 |  |
| `auth_swisnf_long.csv` | OK |     248 | 16 |  |
| `auth_swisnf_panel.csv` | OK |      42 | 43 |  |
| `cellosaurus_str_status.csv` | OK |      42 | 12 |  |
| `consensusov_calls.csv` | OK |      42 | 22 |  |
| `external_ccle_concordance.csv` | OK |     155 | 11 |  |
| `external_depmap_burden.csv` | OK |       5 | 6 |  |
| `external_depmap_driver_crosscheck.csv` | OK |      55 | 12 |  |
| `external_depmap_spearman_all.csv` | OK |   2,077 | 9 |  |
| `external_selfmatch_margin.csv` | OK |       5 | 28 |  |
| `hgs_hallmark_cluster_means.csv` | OK |      50 | 4 |  |
| `hgs_heterogeneity.csv` | OK |      15 | 22 |  |
| `integ_rnaprot_cor.csv` | OK |   7,924 | 6 |  |
| `integ_rnaprot_cor_summary.csv` | OK |       6 | 13 |  |
| `integ_rnaprot_n_distribution.csv` | OK |      14 | 4 |  |
| `integ_rnaprot_n_thresholds.csv` | OK |       5 | 4 |  |
| `package_versions.csv` | OK |      45 | 3 |  |
| `prot_abundance_matrix.csv` | OK |   8,427 | 32 |  |
| `prot_block_missingness.csv` | OK |   8,427 | 8 |  |
| `prot_bridge_agreement.csv` | OK |       4 | 15 |  |
| `prot_bridge_cor.csv` | OK |       4 | 7 |  |
| `prot_compression_floor_check.csv` | OK |       2 | 7 |  |
| `prot_cv_by_abundance.csv` | OK |      10 | 9 |  |
| `prot_dynamic_range.csv` | OK |   7,896 | 21 |  |
| `prot_feature_accounting.csv` | OK |      14 | 3 |  |
| `prot_matrix.rds` | OK | — | — |  |
| `prot_pc_confounder.csv` | OK |       5 | 8 |  |
| `prot_pca.rds` | OK | — | — |  |
| `prot_qc.csv` | OK |   8,427 | 15 |  |
| `prot_sample_annotation.csv` | OK |      31 | 7 |  |
| `prot_sample_qc.csv` | OK |      31 | 7 |  |
| `prot_silhouette.csv` | OK |       6 | 3 |  |
| `prot_variancepartition.csv` | OK |       5 | 12 |  |
| `prot_zero_plex_proteins.csv` | OK |      70 | 12 |  |
| `rna_counts.csv` | OK |  22,544 | 32 |  |
| `rna_dds.rds` | OK | — | — |  |
| `rna_de_CC.csv` | OK |  22,544 | 11 |  |
| `rna_de_EC.csv` | OK |  22,544 | 11 |  |
| `rna_de_HGS.csv` | OK |  22,544 | 11 |  |
| `rna_de_MC.csv` | OK |  22,544 | 11 |  |
| `rna_de_MMMT.csv` | OK |  22,544 | 11 |  |
| `rna_de_SCCOHT.csv` | OK |  22,544 | 11 |  |
| `rna_de_all.csv` | OK | 135,264 | 11 |  |
| `rna_de_gsea_go.csv` | OK |  36,594 | 13 |  |
| `rna_de_gsea_recovery.csv` | OK |       5 | 9 |  |
| `rna_de_gsea_recovery_seeds.csv` | OK |      15 | 10 |  |
| `rna_de_gsea_recovery_stability.csv` | OK |       5 | 10 |  |
| `rna_ec_markers.csv` | OK |     124 | 7 |  |
| `rna_geneset_scores.csv` | OK |      31 | 52 |  |
| `rna_marker_effectsizes.csv` | OK |      26 | 23 |  |
| `rna_markers_summary.csv` | OK |      26 | 15 |  |
| `rna_passage_check.csv` | OK |      18 | 10 |  |
| `rna_passage_discordance.csv` | OK |      22 | 10 |  |
| `rna_pc_confounder.csv` | OK |       5 | 4 |  |
| `rna_pc_confounder_joint.csv` | OK |       5 | 24 |  |
| `rna_pc_confounder_permutation.csv` | OK |       5 | 12 |  |
| `rna_pca.rds` | OK | — | — |  |
| `rna_qc_metrics.csv` | OK |      31 | 7 |  |
| `rna_qc_site_comparison.csv` | OK |      11 | 10 |  |
| `rna_reference_reconciliation.csv` | OK |      31 | 8 |  |
| `rna_sample_annotation.csv` | OK |      31 | 24 |  |
| `rna_signatures_CC.csv` | OK |     365 | 7 |  |
| `rna_signatures_EC.csv` | OK |     270 | 7 |  |
| `rna_signatures_HGS.csv` | OK |   1,858 | 7 |  |
| `rna_signatures_MC.csv` | OK |     153 | 7 |  |
| `rna_signatures_MMMT.csv` | OK |     191 | 7 |  |
| `rna_signatures_SCCOHT.csv` | OK |     278 | 7 |  |
| `rna_signatures_all.csv` | OK |   3,115 | 7 |  |
| `rna_silhouette.csv` | OK |       7 | 3 |  |
| `rna_swisnf.csv` | OK |      31 | 14 |  |
| `rna_tpm.csv` | OK |  39,568 | 32 |  |
| `rna_txi.rds` | OK | — | — |  |
| `rna_variancepartition.csv` | OK |       4 | 12 |  |
| `rna_variancepartition_dropped_genes.txt` | OK | — | — |  |
| `rna_vst.rds` | OK | — | — |  |
| `rna_within_cc_site.csv` | OK |       6 | 10 |  |
| `sensitivity_patient_reps.csv` | OK |      32 | 8 |  |
| `sensitivity_patient_reps_de.csv` | OK | 135,264 | 10 |  |
| `sensitivity_patient_reps_go_comparison.csv` | OK |       5 | 16 |  |
| `sensitivity_patient_reps_gsea_go.csv` | OK |  36,594 | 11 |  |
| `sensitivity_patient_reps_pca.csv` | OK |       9 | 10 |  |
| `sensitivity_patient_reps_silhouette.csv` | OK |      21 | 4 |  |
| `session_info.txt` | OK* | — | — | timestamped provenance; 1000 -> 1002 bytes |
| `session_info_01_rna_load_qc.txt` | OK* | — | — | timestamped provenance; 5434 -> 6778 bytes |
| `session_info_02_rna_separation.txt` | OK* | — | — | timestamped provenance; 6312 -> 6312 bytes |
| `session_info_03_rna_de_signatures.txt` | OK* | — | — | timestamped provenance; 6306 -> 6306 bytes |
| `session_info_04_rna_markers_genesets.txt` | OK* | — | — | timestamped provenance; 7517 -> 7517 bytes |
| `session_info_05_proteomics_load_qc.txt` | OK* | — | — | timestamped provenance; 3553 -> 3553 bytes |
| `session_info_06_proteomics_separation.txt` | OK* | — | — | timestamped provenance; 4833 -> 4833 bytes |
| `session_info_12_rna_protein_concordance.txt` | OK* | — | — | timestamped provenance; 3901 -> 3901 bytes |
| `session_info_13_adc_atlas.txt` | OK* | — | — | timestamped provenance; 4824 -> 4824 bytes |
| `session_info_17_variance_confounders.txt` | OK* | — | — | timestamped provenance; 5820 -> 5820 bytes |
| `session_info_19_proteomics_dynamic_range.txt` | OK* | — | — | timestamped provenance; 3557 -> 3557 bytes |
| `session_info_21_rna_sensitivity.txt` | OK* | — | — | timestamped provenance; 6001 -> 6001 bytes |
| `silhouette_by_modality.csv` | OK |      14 | 5 |  |
| `supplement_per_line.csv` | OK |      42 | 55 |  |
| `tx2gene_ensembl_rel105.csv` | OK | 266,615 | 4 |  |
| `variancepartition_sensitivity.csv` | OK |      16 | 12 |  |
| `wes_cnv_arm_freq_patient.csv` | OK |      44 | 15 |  |
| `wes_cnv_arm_freq_sensitivity.csv` | OK |      60 | 14 |  |
| `wes_cnv_fga.csv` | OK |      23 | 20 |  |
| `wes_cnv_segments.csv` | OK |   5,428 | 10 |  |
| `wes_driver_freq_by_subtype.csv` | OK |      95 | 10 |  |
| `wes_driver_freq_patient.csv` | OK |      95 | 10 |  |
| `wes_driver_tiers.csv` | OK |      50 | 13 |  |
| `wes_hrd_feasibility.md` | OK | — | — |  |
| `wes_msi_mmr.csv` | OK |      22 | 38 |  |
| `wes_mutation_load.csv` | OK |      22 | 21 |  |
| `wes_mutations_filtered.csv` | OK |   6,036 | 23 |  |
| `wes_pipeline_parameters.csv` | OK |      15 | 6 |  |
| `wes_sbs_context.csv` | OK |      96 | 23 |  |
| `wes_sbs_cosine.csv` | OK |   1,320 | 7 |  |
| `wes_signature_refit_bootstrap.csv` | CHANGED |   1,250 ->   1,804 | 14 | complete grid vs baseline's selected-only subset; see above |
| `wes_signature_refit_exposures.csv` | OK |   1,250 | 22 |  |
| `wes_signature_refit_summary.csv` | OK |      22 | 34 |  |
