# Fix report — RNA / proteomics / integration arm

**Scope:** `scripts/01`, `02`, `03`, `04`, `05`, `06`, `12`, `13`, `17`, `19` + new `scripts/21_rna_sensitivity.R`.
**Environment:** R 4.5.2, `Rscript`, `OVCAN_PROJ` default. Every script below was executed to completion after editing, and the whole set was then re-executed in dependency order (`01 → 05 → 02 → 06 → 03 → 04 → 12 → 13 → 19 → 17 → 21`), all 11 completing without error.
**Not touched:** `00_setup.R`, `07`–`11`, `15`, `16`, `18`, `20`, all `3*` figure scripts, `docs/`, `reports/_snapshot_pre_revision/`, `output/external/`.

**Concurrent edit by another agent — please read.** Part-way through this work an "[integration revision]" pass by a parallel agent edited six of the files in my scope (`06`, `12`, `13`, `17`, `19`, `21`), replacing my inline protein-matrix loads and family-map reads with two new shared helpers it added to `00_setup.R`: `read_prot_matrix()` (which carries the zero-plex row names as an attribute) and `ensure_family_map()`. I left those edits in place and re-ran all 11 scripts afterwards; all 11 completed without error and 55 of 55 headline values in this report reproduced.

**Caveat on that verification, stated plainly.** The same agent then took ownership of all of `scripts/` and began a clean-room reproducibility run, emptying `output/` into `output/_preclean/` at 15:15. My re-verification ran at ~15:30, i.e. against an `output/` that two processes had been writing to. The values are reproducible *from the code* — every script is deterministic (`fgsea` is seeded per call at `nPermSimple = 50000`; the only stochastic steps are seeded permutation and bootstrap loops) — so I expect them to hold, but **the authoritative check is the integration agent's clean-room run, not mine.** Any value in this report that its run does not reproduce should be treated as a real discrepancy and investigated, not reconciled away.

Two things I confirmed read-only afterwards, both relevant to that run: `output/tx2gene_ensembl_rel105.csv` was *not* preserved through the wipe — `01` regenerated it via a **live biomaRt query** at 15:17 — but it came back **byte-identical** to the cached copy (md5 `119dfe0ab4f856e6b81efc6ce78f4ba7`, 266,615 transcripts), so the pinned Ensembl 105 archive reproduces exactly. And `metadata/line_family_map.csv`, recreated by `ensure_family_map()`, is byte-identical to `reports/_snapshot_pre_revision/`.

Note also that `06`/`12`/`13`/`17`/`19`/`21` now have a hard dependency on `read_prot_matrix()` and `ensure_family_map()` existing in `00_setup.R`, and that `scripts/` is no longer mine to edit.

---

## 1. Fix table

| # | File:line | What changed | Before → after |
|---|---|---|---|
| 1 | `01:36-49` (comments), `01:75-124` (new §3b) | Reference provenance corrected in the comments (index = **Ensembl-104-era**, `index_version 10`, `n_targets 185299`, quantification started 2021-10-20; `tx2gene` = 105 as nearest live archive). New §3b quantifies the loss per sample, writes `output/rna_reference_reconciliation.csv`, and guards it with 5 `stopifnot`s. | undocumented → **3,529 / 185,299 targets (1.90%)** dropped, carrying **2.22% of TPM** (range 1.60–3.41%) and **1.53% of estimated counts** (1.38–1.77%) |
| 2 | `01:141-182` | `lib_size` → `assigned_gene_counts` (it is `colSums` of post-filter gene-level estimated counts, not a read count). New `n_processed_fragments` parsed from each `run_info.json`, with a comment stating these are **fragments** (paired-end), not reads. `jsonlite` added to the library block. | `lib_size` median 56,497,657 → same value, renamed `assigned_gene_counts`. New `n_processed_fragments` **median 64.3 M, range 45.1–97.0 M** (matches the expected value exactly) |
| 3 | `01:184-243` (new §6b) | Per-site QC comparison table `output/rna_qc_site_comparison.csv`: per-site n + medians for pseudoalignment / genes detected / both depth measures, Wilcoxon rank-sum tests, and three Pearson correlations across all 31 lines. Narrative comment rewritten. | untested absence claim → Mes-Masson (n=19) 92.2% / 19,817 genes; Huntsman (n=11) 88.1% / 20,765 genes. **Wilcoxon genes detected p = 0.00819** (exact) / **0.00982** (normal approx). **Pearson(pseudoalign, detected) r = −0.635, p = 1.26e−4** (expected ≈ −0.63) |
| 4 | `05:109-196` | Feature hygiene rewritten. The 3 no-symbol rows are dropped with a logged count and their accessions printed (A6NIZ1, A6NNZ2, Q6ZSR9) instead of entering as `NA` / `NA.1` / `NA.2`. Duplicate symbols keep all rows but exactly one is the documented deterministic **symbol representative** (most non-missing lines → most `Npeptides_quant` → lowest q-value → first accession); non-representatives are renamed `SYMBOL\|UNIPROT`. New `output/prot_feature_accounting.csv`. | matrix 8,430 → **8,427** rows; complete-case 6,856 → **6,855**; absent-from-≥1-plex 1,573 → **1,572**; **6,855 + 1,572 = 8,427** (the off-by-one is gone); ≥50%-presence 7,734 → **7,733**; distinct symbols **8,396**; non-representative duplicate rows **31** |
| 5 | `05:198-262`, `06:159-166` | Zero-plex proteins characterised to `output/prot_zero_plex_proteins.csv` and flagged `zero_plex = TRUE` in `prot_qc.csv`. `n_any` is used (it feeds the accounting table and the console summary — note it was already used at the old `05:316`, so the audit's "never used" is incorrect). `06`'s `match(mk$gene, feat$symbol)` now resolves through the representative instead of taking the first row. | zero-plex 71 → **70** (one of the dropped no-symbol rows was zero-plex); `n_any` 8,359 → **8,357** |
| 6 | `17:214-290` (new §A2) | Within-clear-cell site control refitted **within clear cell only** (top-2000 HVG re-selected on the 7 CC lines). Both statistics written to `output/rna_within_cc_site.csv` with explicit `approach` labels, plus an ANOVA p, a within-group permutation p (2,000 draws) and the null expectation. | OLD (global PCs, wrong axes) PC1–3 R² = 0.044 / 0.058 / 0.057 → **CORRECT (within-CC PCA) 0.157 / 0.126 / 0.273**, ANOVA p **0.380 / 0.435 / 0.228**, permutation p **0.490 / 0.541 / 0.246**, E[R²] = 1/6 = **0.167**. All expected values reproduced exactly |
| 7 | `17:110-212` | Adjusted commonality components added (`adj_unique_subtype`, `adj_unique_site`, `adj_shared`) alongside the raw ones; `adjr2_site_2lvl` added; a `headline_site_variable` column states the 3-level variable is the headline; permutation nulls (1,000 draws each) written to `output/rna_pc_confounder_permutation.csv`; PC1–PC3 unique-site printed against the claimed ≤0.2%. | PC1 raw unique-site 0.0019 → **adjusted −0.0251**; PC1 raw unique-subtype 0.4238 → **adjusted 0.3928**; raw joint R² 0.7369 vs **adjusted 0.6569**. **PC3 unique-site = 0.24%, which exceeds the claimed ≤0.2%** (confirmed) |
| 8 | `17:333-368`, `17:370-448` | `vp_lmer` now returns the dropped features as an attribute; `vp_summary` gains `model`, `mean_pct`, `q25_pct`, `q75_pct` (already present), `n_features_input`, `n_features_dropped`, `n_samples` and a `design` string recording patient-level and site-level counts. Two RNA sensitivity models + one protein sensitivity model written to `output/variancepartition_sensitivity.csv`. Dropped genes listed in `output/rna_variancepartition_dropped_genes.txt`. Figure facet labels now use the fitted feature count. | dropped genes silent → **2 logged** (`ENSG00000223672`, `ENSG00000233377`); `n_features` = **22,542** (not 22,544). **`patient` dropped**: subtype 5.95 → **7.09%**, site 3.53 → **3.95%**, residual 32.8 → **72.3%**. **`family` restricted to the 3 replicated families**: subtype **5.95%**, site **3.35%**, family median **0.76%** (mean 11.4%) — i.e. the 27.3% "patient identity" term collapses to 0.76% once it is estimated from the replication that actually exists |
| 9 | `17:120-137`, `17:333-368`, `17:527-556` | Passage added to the joint PC model (`r2_passage`, `r2_joint3_with_passage`, `unique_passage_beyond_subtype_site`) and to the per-gene decomposition as a scaled **fixed** covariate (variance = `var(predict(fit, re.form = NA))`, the same quantity `variancePartition` attributes to fixed terms). RNA-vs-WES discordance written to `output/rna_passage_discordance.csv` with a summary block. | PC1 R²(passage) **0.078**; PC2 **0.159**. Passage beyond subtype **and** site: PC1 **0.0041**, PC2 0.0201, PC4 0.0204. Per gene: passage median **4.62%** / mean 9.15%, and adding it drops the **site median from 3.53% to 0.00%**. Discordance: 13 lines, median \|diff\| **4.0**, mean 5.8, **range −17 to +20 passages**, 2 lines ≥10, 1 line ≥20 |
| 10 | `17:630-712` | BH adjustment added (`wilcox_p_bh`); the smallest attainable exact p per group size added (`wilcox_p_min_attainable`) with a `wilcox_p_floor_limited` flag; stratified bootstrap CIs (2,000 replicates, resampling within the intended group and within the rest) for **every** Cohen's *d* and oriented AUC, plus `n_boot_valid`. | 22 unadjusted p → **25 with BH**. **4/25 survive BH at 0.05** (HNF1B, SPP1, KRT20, SMARCA2, all BH = 0.0437). Min attainable p: **5.29e−3** for n=2 (matches the expected ~0.005), 6.11e−4 for n=3, 1.69e−6 for n=7. **1/25 floor-limited (SMARCA2)**. Example CIs: SMARCA2 *d* = −2.08 [−3.36, −1.53]; KRT20 *d* = 2.84 [1.72, **18.97**]; CDX2 *d* = 3.08 [0.96, 9.55] |
| 11 | `17:591-628` | Marker quantification harmonised to **script 04's rule** — collapse `gene_id → symbol` by summing TPM, then `log2(TPM+1)` — replacing the previous VST / max-mean-ENSG rule, so `lands_right` and the effect sizes are now the same quantity. Rule stated in a comment. | **Same 22 markers:** median \|*d*\| 0.720 → **0.698**; median oriented AUC 0.690 → **0.706**; AUC ≥ 0.80 **8 → 7**; \|*d*\| ≥ 0.8 10 → 10. Per marker: KRT20 3.082 → **2.842**, CDX2 2.774 → **3.075**, SPP1 1.461 → **1.504**, HNF1B 1.364 → **1.490**, PAX8 0.529 → **0.408**, MUC2 −0.294 → **−0.107**, MECOM −0.064 → **−0.274**, SMARCA2 −2.776 → **−2.082**. `lands_right` unchanged for all 22 |
| 12 | `04:88-90`, `04:110-125`, `04:145-148`, `04:190-232` | Endometrioid marker set added (**ESR1** up, **PGR** up, **ARID1A** down, per Hollis et al. 2020), scored with the same top-2-of-6 rule and the same effect-size machinery in `17`. New `output/rna_ec_markers.csv` reports per-line values and z-scores, at n=2 **and** for VOA4395 alone. A binomial null for the top-2-of-6 rule was also added. | **No EC marker lands.** ESR1 EC mean **0.03** log2 TPM (rank 6/6), PGR **0.00** (rank 5/6), ARID1A rank 2 (needs ≥5). Effect sizes: ESR1 oriented AUC **0.173** [0.00–0.42], PGR **0.519**, ARID1A **0.231**; none significant (BH 0.35 / 0.96 / 0.46). Landing count **16/22 → 16/25 graded**; binomial vs p=1/3: **p = 0.0016** (all graded), **0.0035** ('up' only) |
| 13 | `12:96-104`, `12:193-263` | Full per-gene n distribution written (`integ_rnaprot_n_distribution.csv`, `integ_rnaprot_n_thresholds.csv`); complete-case subset added to the summary as a headline-comparable statistic; negative fraction recorded; the two conflated counts written into `integ_rnaprot_cor_summary.csv` as columns. | shared symbols **8,212** and per-gene correlations **7,893 → 7,894** now both recorded and printed as separate numbers. Complete case 6,685 → **6,688 genes, median 0.408**; overall median **0.397**; **10.3% of genes negative**; n per gene **10–30**, 353 at n<15, 623 at n<20 |
| 14 | `19:91-94`, `19:132-152`, `19:186-243`, `19:360-437` | IQR ratio and SD ratio made the **primary** compression estimate; `range_ratio_floor_sensitive` added and the range ratio labelled secondary; `rna_n_zero_lines` / `rna_has_zero` added and a floor-split diagnostic written (`prot_compression_floor_check.csv`); Bland–Altman/MA bridge agreement written (`prot_bridge_agreement.csv`) with per-link n; CV and bridge noise stratified by abundance decile (`prot_cv_by_abundance.csv`). | **Compression: IQR ratio 0.30, SD ratio 0.32, range ratio 0.34 — all agree at ~3.3×.** MSLN: range ratio 0.189 but **IQR ratio 0.104, SD ratio 0.141**. Bridge: per-link n **7,276–7,436**, bias −0.028 to +0.010, repeatability SD **0.211–0.268 log2 (15.7–20.4% CV)**, 95% LoA span **0.83–1.05 log2** vs a median cross-line protein IQR of **0.34**. Vendor CV falls **11.0% → 3.0%** and bridge SD **0.430 → 0.106** from the lowest to the highest abundance decile |
| 15 | `03:195-224` and `03:267-303` | `set.seed(SEED)` moved inside `run_gsea()`. This exposed a second defect: at fgsea's default `nPermSimple = 1000` the p-values are not resolved near the 0.05 boundary, so `nPermSimple` was raised to **50,000** and a 3-seed stability check added (`rna_de_gsea_recovery_stability.csv`, `rna_de_gsea_recovery_seeds.csv`). | See §2 — **two GO recovery grades change.** All five grades are now stable across 3 seeds |
| 16 | `13:86-89`, `13:131-166`, `13:215-290` | Per-subtype n emitted for each modality (`n_measured`, `n_lines_in_modality` in `adc_subtype_summary.csv`); the modality-specific lines written to `adc_modality_line_sets.csv`. FOLR1 bimodality formally tested (`adc_folr1_bimodality.csv`). | Modality asymmetry recorded: **VOA6861 RNA-only, VOA14993 protein-only** (both CC). **FOLR1 is NOT bimodal within HGSC at n=15**: bimodality coefficient **0.440** (threshold 0.555, not flagged); mclust ΔBIC(G2−G1) = **−1.06** (equal-variance) / **−0.95** (variable-variance) → **one component preferred**. The 2-component fit does split 4 low lines (mean 0.65) from 11 high (mean 6.49), but the data do not support it |
| 17 | new `scripts/21_rna_sensitivity.R` | Patient-representative sensitivity: PCA + silhouette + PC confounder + commonality on the 28 representatives (two HVG variants), and the previously uncomputed one-vs-rest DESeq2 + GO recovery for all 6 subtypes on the 28 reps. | See §2. Every expected value reproduced exactly: HGS silhouette **0.162 → 0.142**, overall **0.220 → 0.208**, PC1 **20.75% → 19.93%**, PC1 R²(subtype) **0.735 → 0.792** |
| 18 | `02:56-66`, `06:98-115` | `02` now writes the overall RNA silhouette as an `ALL` row. `06` writes `output/silhouette_by_modality.csv` with an explicit `modality` column, an `ALL` row per modality and a `unit` column. | Both modalities in one labelled table. Protein: HGS **0.178**, CC **−0.003**, EC **0.128**, MC **0.141**, MMMT **0.409**, SCCOHT **0.028**, overall **0.135**. RNA overall **0.220** now on disk |
| 19 | end of `01`, `02`, `03`, `04`, `05`, `06`, `12`, `13`, `17`, `19`, `21` | `write_session_info()` writes `output/session_info_<script>.txt` containing script name, run timestamp, seed, full `sessionInfo()` and a sorted attached+loaded package-version table. | 11 session records written |

---

## 2. Every number that changed, for the manuscript

### 2.1 New / corrected numbers

| Claim as written | Corrected value | Source |
|---|---|---|
| "index built from Ensembl release 105" | index is **Ensembl-104-era** (built before 2021-10-20); `tx2gene` is release 105, the nearest live archive. **1.90% of index targets (3,529/185,299)** are silently dropped, carrying **2.22% of TPM** (1.60–3.41%) and **1.53% of estimated counts** (1.38–1.77%) | `rna_reference_reconciliation.csv` |
| "median library size of 56.5 M reads" | **assigned gene counts** median 56,497,657 (not a read count). True sequenced depth = **64.3 M fragments** (45.1–97.0 M); paired-end, so ~128.6 M reads if reads are meant | `rna_qc_metrics.csv` |
| "does not translate into a difference in genes detected" | It does: Huntsman detects **948 more** genes (20,765 vs 19,817), **Wilcoxon p = 0.0082**; across 31 lines pseudoalignment and gene detection are **negatively** correlated (r = −0.635, p = 1.3e−4) | `rna_qc_site_comparison.csv` |
| "8,430 proteins … 6,856 (81.3%) … 1,573 (18.7%)" | **8,427 = 6,855 (81.3%) + 1,572 (18.7%)**. Search output is 8,430 rows; 3 have no gene symbol and are dropped | `prot_feature_accounting.csv` |
| ≥50%-presence filter retains 7,734 | **7,733** | `prot_qc.csv` |
| 71 proteins present in 0 plexes | **70** | `prot_zero_plex_proteins.csv` |
| "within clear cell, source lab explains only 4–6% of the leading components" | Wrong axes. Within-CC PCA gives **15.7 / 12.6 / 27.3%**, ANOVA p **0.38 / 0.44 / 0.23**, permutation p **0.49 / 0.54 / 0.25**, against a null expectation of **16.7%**. Correct conclusion: **underpowered, indistinguishable from chance** | `rna_within_cc_site.csv` |
| "source site adds ≤0.2% of PC1 variance beyond histotype (unique-subtype 42%, unique-site 0.2%, shared 31%)" | Raw values correct. **Adjusted: unique-subtype 39.3%, unique-site −2.5%, shared 28.9%; adjusted joint R² 0.657 vs raw 0.737.** Permutation p: unique-subtype **0.002**, unique-site **0.562** (observed below its own null mean) | `rna_pc_confounder_joint.csv`, `rna_pc_confounder_permutation.csv` |
| "the same holds on PC2–PC3 (unique-site ≤0.2%)" | **PC3 unique-site = 0.24%**, above 0.2%. PC2 = 0.05%. Also note unique-subtype is **not** significant on PC4 (perm p = 0.369) or PC5 (0.075) | same |
| "subtype ≥ site in RNA (median 5.9% vs 3.5%)" | Medians confirmed (5.949 / 3.525) but **means are 14.94% vs 14.50%** and IQRs are 0–24.79% vs 0–25.19%. All four now in the output table | `rna_variancepartition.csv` |
| "line/patient identity the largest structured term (RNA 27.3%)" | Design artefact. Restricting the term to the 3 replicated families gives **family median 0.76%**; dropping `patient` moves subtype 5.95→7.09% and site 3.53→3.95%, preserving the ordering | `variancepartition_sensitivity.csv` |
| Fig 3D "22,544 genes" | **22,542** (2 genes dropped, now listed) | `rna_variancepartition.csv`, `rna_variancepartition_dropped_genes.txt` |
| Fig 3D "6,856 proteins" | **6,855** | `prot_variancepartition.csv` |
| passage "not an independent structure driver" | R²(passage) on PC1 = **0.078**, PC2 = **0.159**; partial R² after site alone = 0.148, but **after subtype *and* site it is only 0.004 on PC1**. Per gene, passage takes **4.62% median** and **absorbs the entire site median (3.53% → 0.00%)** | `rna_pc_confounder_joint.csv`, `variancepartition_sensitivity.csv` |
| — (absent) | RNA-vs-WES passage discordance: 13 lines, median \|diff\| **4**, **range −17 to +20 passages**, 1 line ≥20 | `rna_passage_discordance.csv` |
| "on 30 lines and 8,212 shared genes, per-gene median Spearman is 0.40" | **8,212 shared symbols** but only **7,894 receive a correlation** (n≥10). Median **0.397**; **complete-case (n=30 for every gene, 6,688 genes) median 0.408**; **10.3% negative**; per-gene n **10–30** | `integ_rnaprot_cor_summary.csv` |
| "TMT compression narrows target ranges 3–5× (MSLN protein 1.8 vs RNA 9.3)" | Use the floor-insensitive statistics: **IQR ratio 0.30, SD ratio 0.32 → ~3.3×**. For MSLN specifically the floor-insensitive ratios show **more** compression (IQR 0.104 = 9.6×, SD 0.141 = 7.1×), not less | `prot_dynamic_range.csv`, `prot_compression_floor_check.csv` |
| bridge "Pearson 0.991–0.994" | Report agreement instead: per-link n **7,276–7,436**, bias **−0.028 to +0.010** log2, repeatability SD **0.211–0.268** log2 (**15.7–20.4% CV**), 95% LoA span **0.83–1.05** log2 — **wider than the median cross-line protein IQR of 0.34** | `prot_bridge_agreement.csv` |
| "FOLR1 is strongly bimodal within HGSC" | **Not supported.** Bimodality coefficient **0.440** (< 0.555); mclust **ΔBIC = −1.06 / −0.95** → one component preferred. Describe it as a wide range with a near-zero cluster, not as bimodal | `adc_folr1_bimodality.csv` |
| silhouettes "HGSC 0.16, CC 0.12, MC 0.15, MMMT 0.74, SCCOHT 0.82" | These are **RNA**. Protein: HGS **0.178**, CC **−0.003**, EC **0.128**, MC **0.141**, MMMT **0.409**, SCCOHT **0.028**; overall RNA **0.220** vs protein **0.135** | `silhouette_by_modality.csv` |
| "16 of 22 markers land" | **16 of 25 graded** after the EC set is added (the 3 EC markers all fail). Binomial vs the 1/3 null: **p = 0.0016** | `rna_markers_summary.csv` |
| marker median \|*d*\| 0.72, median AUC 0.69, 8/22 at AUC ≥ 0.80 | On the harmonised quantification, same 22 markers: **0.698 / 0.706 / 7 of 22**. Across all 25: **0.710 / 0.667 / 7 of 25** | `rna_marker_effectsizes.csv` |
| — (absent) | Endometrioid identity is **not marker-validated**: ESR1 and PGR are effectively unexpressed in all 31 lines (EC means 0.03 and 0.00 log2 TPM) and ARID1A is not low in EC. VOA4395, the sole EC model, has ESR1 = 0.000 and PGR = 0.000 log2 TPM | `rna_ec_markers.csv` |

### 2.2 GO recovery — two grades change from the `fgsea` precision fix (item 15)

At the archived setting (`nPermSimple = 1000`, seed set once) the grades were partly RNG artefacts. At `nPermSimple = 50,000` all five are stable across 3 seeds.

| Program | Published grade | Corrected grade | padj |
|---|---|---|---|
| HGS DNA repair (HR/DSB) | recovered (padj = 0.0498) | **suggestive (nominal p<0.05)** | 0.0498 → **0.0748** (3-seed range 0.068–0.075) |
| CC oxidative/glutathione | recovered | recovered (unchanged) | 0.0045 → 0.0072 |
| MC glycan/oligosaccharide | suggestive | suggestive (unchanged) | 0.117 → 0.120 |
| MMMT EMT / migration / ECM | **not recovered** | **suggestive (nominal p = 0.046)** | 0.458 → 0.338 |
| SCCOHT cell cycle / mitosis | recovered | recovered (unchanged) | 6.0e−5 → 8.9e−5 |

The HGS "recovered" claim was a coin flip on the RNG (0.0498 on one draw, 0.0500 on the next) and does not survive a precise estimate. `n_pos_sig` for SCCOHT also moves 24 → 26.

### 2.3 Patient-representative sensitivity (item 17) — the previously open question

31 lines / 28 patients; HGS 15 lines / 12 patients. Full table in `output/sensitivity_patient_reps.csv`.

**Structure (all four expected values reproduced exactly):**

| Metric | 31 lines | 28 reps | Δ |
|---|---|---|---|
| PC1 variance explained | 20.75% | **19.93%** | −0.81 |
| PC1 R²(subtype) | 0.735 | **0.792** | +0.057 |
| PC1 R²(site, 3-level) | 0.313 | 0.272 | −0.041 |
| PC1 unique-subtype | 0.424 | 0.523 | +0.099 |
| PC1 unique-site | 0.0020 | 0.0026 | +0.0006 |
| PC1 shared | 0.311 | 0.270 | −0.042 |
| PC1 adjusted joint R² | 0.657 | 0.723 | +0.066 |
| Silhouette HGS | 0.162 | **0.142** | −0.020 |
| Silhouette CC | 0.121 | 0.104 | −0.017 |
| Silhouette EC | −0.014 | −0.048 | −0.034 |
| Silhouette MC | 0.154 | 0.146 | −0.008 |
| Silhouette MMMT | 0.738 | 0.724 | −0.014 |
| Silhouette SCCOHT | 0.818 | 0.798 | −0.021 |
| Silhouette overall | 0.220 | **0.208** | −0.012 |

The separation conclusion is robust and the direction is favourable: collapsing to patients *raises* PC1 R²(subtype) and unique-subtype and *lowers* the site and shared components.

**DE / GSEA — not computed by anyone before. This is the part that moves.**

| Contrast | n_group | DE genes padj<0.05 | Up-signature (padj<0.05 & log2FC>1) |
|---|---|---|---|
| HGS | 15 → 12 | 7,913 → **6,372** (−19.5%) | 2,088 → **1,646** (−21.2%) |
| CC | 7 → 7 | 1,256 → **821** (−34.6%) | 400 → **286** (−28.5%) |
| EC | 2 → 2 | 1,880 → 1,780 | 277 → 266 |
| MC | 3 → 3 | 565 → 580 | 168 → 175 |
| MMMT | 2 → 2 | 1,387 → 1,287 | 199 → 186 |
| SCCOHT | 2 → 2 | 2,232 → 2,152 | 290 → 287 |

Read carefully: the CC/EC/MMMT/SCCOHT groups do **not** lose members, yet their DE counts move, because the "rest" group shrinks from 24 to 21 lines. So part of every one-vs-rest contrast's power came from the duplicate genomes, not only the HGS contrast's.

**GO recovery under patient collapse** (`sensitivity_patient_reps_go_comparison.csv`):

| Program | n | padj (31) | padj (28) | Grade |
|---|---|---|---|---|
| HGS DNA repair | 15 → 12 | 0.0748 | 0.132 | suggestive → suggestive (best term shifts from "regulation of" to "positive regulation of" DSB repair via HR) |
| CC oxidative/glutathione | 7 → 7 | 0.00723 | 0.00773 | recovered → recovered |
| MC glycan | 3 → 3 | 0.120 | 0.153 | suggestive → suggestive |
| **MMMT EMT/migration/ECM** | 2 → 2 | 0.338 | 0.460 | **suggestive → not recovered** (nominal p 0.046 → 0.083) |
| SCCOHT cell cycle | 2 → 2 | 8.9e−5 | 7.4e−5 | recovered → recovered |

**Summary for the manuscript:** the GO recovery of expected biology is qualitatively robust to patient collapse for HGS, CC, MC and SCCOHT. It is **not** robust for MMMT, which loses its nominal significance. Combined with §2.2, the honest statement is: **CC and SCCOHT recover their expected programs robustly; HGS DNA repair and MC glycan metabolism are suggestive only; MMMT EMT is not recovered.**

---

## 3. New output files

| File | Contents |
|---|---|
| `output/rna_reference_reconciliation.csv` | Per line: index target count, unmapped-target count, and the unmapped fraction of targets / TPM / estimated counts, plus index and tx2gene provenance strings |
| `output/rna_qc_site_comparison.csv` | Three labelled blocks: per-site medians (pseudoalignment, genes detected, both depth measures); Wilcoxon rank-sum tests (exact and normal-approximation for genes detected); three Pearson correlations across all 31 lines |
| `output/prot_feature_accounting.csv` | 14 rows reconciling every protein count from the 8,430 search rows down to the 7,733 presence-filtered set, each with a definition string |
| `output/prot_zero_plex_proteins.csv` | The 70 identified-but-unquantified proteins: row name, symbol, accession, peptide counts, q-value, CV, isoDoping flag |
| `output/silhouette_by_modality.csv` | RNA and protein silhouettes in one table with `modality`, `subtype` (incl. an `ALL` row), `n` and a `unit` column |
| `output/rna_ec_markers.csv` | Per-line log2 TPM and z-across-31-lines for ESR1, PGR, ARID1A and VIM, with the group each is scored under |
| `output/rna_de_gsea_recovery_seeds.csv` | The recovery grading repeated at 3 fgsea seeds |
| `output/rna_de_gsea_recovery_stability.csv` | Per program: padj and p ranges across seeds, distinct statuses, and a `stable` flag |
| `output/integ_rnaprot_n_distribution.csv` | Genes per per-gene n value, cumulative and percentage |
| `output/integ_rnaprot_n_thresholds.csv` | Genes below/at-or-above n thresholds 10/15/20/25/30 with the median Spearman at each |
| `output/adc_modality_line_sets.csv` | The lines present in only one modality (VOA6861 RNA-only, VOA14993 protein-only) |
| `output/adc_folr1_bimodality.csv` | FOLR1-within-HGS: n, min/max/range/SD, bimodality coefficient, mclust BIC for G=1 and G=2 under equal- and variable-variance models, ΔBIC, fitted component means and proportions, and an explicit "dip test not computed" row |
| `output/prot_compression_floor_check.csv` | Genes split by whether the RNA zero-floor is engaged, with median range/IQR/SD ratios and median RNA range/IQR each side |
| `output/prot_bridge_agreement.csv` | Per bridge link: n proteins, bias, SD of difference, 95% limits of agreement, median absolute difference, repeatability CV%, MA-plot slope and its p, plus Pearson/Spearman for comparison |
| `output/prot_cv_by_abundance.csv` | Mean-abundance decile: n proteins, abundance bounds, vendor CV median/IQR, bridge SD of difference, bridge median absolute difference |
| `output/rna_within_cc_site.csv` | Within-clear-cell site control, both the old (global-PC subset) and correct (within-CC PCA) approaches, with ANOVA p, permutation p and the null expectation |
| `output/rna_pc_confounder_permutation.csv` | Observed unique-subtype and unique-site with their permutation null mean/median/95th percentile and empirical p, for PC1–PC5 |
| `output/variancepartition_sensitivity.csv` | The three RNA and one protein sensitivity models, with full mean/median/IQR, feature counts, dropped counts and a `design` string |
| `output/rna_variancepartition_dropped_genes.txt` | The 2 gene IDs the lmer decomposition drops |
| `output/rna_passage_discordance.csv` | Per-line RNA vs WES passage with signed and absolute difference, plus a summary block (n, median/mean \|diff\|, min/max, counts ≥10 and ≥20, Spearman) |
| `output/sensitivity_patient_reps.csv` | 31 metrics as `analysis`, `metric`, `value_31_lines`, `value_28_reps`, `delta`, `note`, `unit_31`, `unit_28` |
| `output/sensitivity_patient_reps_pca.csv` | PC1–PC3 variance, R², commonality components for 31 lines and both 28-rep HVG variants |
| `output/sensitivity_patient_reps_silhouette.csv` | Silhouettes per subtype for the same three sets |
| `output/sensitivity_patient_reps_de.csv` | Full one-vs-rest DESeq2 results for all 6 subtypes on the 28 reps |
| `output/sensitivity_patient_reps_gsea_go.csv` | Full fgsea GO-BP results on the 28 reps |
| `output/sensitivity_patient_reps_go_comparison.csv` | Side-by-side recovery grading, 31 lines vs 28 reps, with `status_changed` and `same_best_term` flags |
| `output/session_info_*.txt` | 11 files: script, timestamp, seed, `sessionInfo()`, sorted package-version table |

**Changed in shape (columns added, existing columns preserved):** `rna_qc_metrics.csv` (`lib_size` → `assigned_gene_counts`, + `n_processed_fragments`), `prot_qc.csv` (+ `row`, `symbol_representative`, `zero_plex`), `prot_abundance_matrix.csv` / `prot_matrix.rds` / `prot_block_missingness.csv` (8,430 → 8,427 rows; `SYMBOL.N` → `SYMBOL|UNIPROT`), `rna_silhouette.csv` (+ `ALL` row), `rna_pc_confounder_joint.csv` (+ 12 columns), `rna_variancepartition.csv` / `prot_variancepartition.csv` (+ `model`, `n_features_input`, `n_features_dropped`, `n_samples`, `design`), `rna_marker_effectsizes.csv` (26 rows, + 8 columns), `rna_markers_summary.csv` (23 → 26 rows), `prot_dynamic_range.csv` (+ `range_ratio_floor_sensitive`, `rna_n_zero_lines`, `rna_has_zero`), `integ_rnaprot_cor_summary.csv` (+ 2 rows, + 5 columns), `adc_subtype_summary.csv` (+ 3 n columns).

---

## 4. Items I could not complete

Nothing in the brief was skipped. Two items were delivered with a substitution, both flagged in the code:

1. **Item 16, Hartigan dip test.** `diptest` is not installed and I did not add a dependency. The bimodality question is answered instead by a 2-component-vs-1-component `mclust` mixture fit (BIC, both equal- and variable-variance) and by the moment-based bimodality coefficient. Both agree that one component is preferred. `adc_folr1_bimodality.csv` carries an explicit `hartigan_dip_test = NA` row stating why. Installing `diptest` and adding `dip.test(x)` at `13:242` would complete it; the conclusion is unlikely to change.
2. **Item 5, "flag them in the deposited matrix".** `prot_abundance_matrix.csv` is a plain numeric matrix read as `as.matrix(prot[, -1])` by `12`, `13`, `19` **and by `10:95` and `11:62`, which I do not own**, so adding a flag column would break those two scripts. The 70 zero-plex rows are therefore retained in the matrix (not dropped) and flagged in `prot_qc.csv` (`zero_plex`) plus listed in `prot_zero_plex_proteins.csv`. If you want the flag inside the matrix file itself, `10` and `11` need a matching change.

---

## 5. Needed from files I do not own

| File | What is needed | Why |
|---|---|---|
| `35_fig3_biology.R:116-117` | Change the Fig 3D facet labels `"RNA (22,544 genes)"` → **22,542** and `"Protein (6,856 proteins)"` → **6,855**. Better: read `n_features` from `rna_variancepartition.csv` / `prot_variancepartition.csv`, which now carry it. | The hard-coded counts are wrong; `recode()` matches on the literal string so a stale label silently mislabels the panel |
| `35_fig3_biology.R` panel E | Fig 3E will now render **26 marker rows** (23 before) with a new `EC` group band. It will not error — `expected_cell` already has `EC` in its factor levels and `gord` at `:143` is dead code — but the panel geometry and height need a look. Consider annotating the EC block as *not* recovered. | Item 12 added ESR1/PGR/ARID1A to `rna_markers_summary.csv` |
| Manuscript text, §2 (not a script fix) | `35_fig3_biology.R` prints no *n* for panel C, so nothing breaks — but the **text**'s "8,212 shared genes" must become **7,894** (the number that receives a correlation), and 10.3% of genes are negatively correlated. Both are now columns in `integ_rnaprot_cor_summary.csv` (`n_pergene_reported`, `frac_pergene_negative`) so the figure caption can carry them. | Item 13 |
| `37_supp_rnaprot.R:186` | Add `"EC"` to `grp_lvl <- c("HGS","CC","MC","MMMT","SCCOHT")`. Without it the 3 EC markers get `grp = NA` and drop out of the Fig S7 grouping/colouring. Fig S7 is now **25 markers**, not 22. | Item 12 |
| `37_supp_rnaprot.R` Fig S7 | `rna_marker_effectsizes.csv` now carries `cohens_d_lo/hi` and `auc_oriented_lo/hi` — bootstrap CIs are available as error bars, which is what referee M-minor 9 asked for. `n_boot_valid` marks the unstable ones. | Item 10 (opportunity, not a break) |
| `34_fig2_qc.R` | No change required — it uses `pseudoalign`, `detected`, `dr$iqr_ratio`, `rna_iqr`, `prot_iqr` and `bm$present_n_plex`, all preserved. But Fig 2C now sums to **8,427** (not 8,429/8,430) and the "0 plexes" bar is **70**. Panel D's `n_total` will read 8,427 automatically. | Items 4–5 |
| `20_supplement_table.R` | **Already handled** by the parallel agent (`:98-111` reads `assigned_gene_counts` with a `lib_size` fallback and adds `rna_sequenced_fragments_M`). No action needed; noting it so it is not undone. | Item 2 |
| Whoever wipes `output/` for a clean run | **`output/tx2gene_ensembl_rel105.csv` must be preserved.** It is a hand-cached input, not a reproducible intermediate: `01:49-65` falls back to a **live biomaRt query against the Ensembl 105 archive** on a cache miss. biomaRt 2.64.0 is installed so it will not hard-fail, but it needs network access to `dec2021.archive.ensembl.org` — which is exactly the non-reproducibility the pinned cache exists to avoid — and any archive drift silently changes the 39,568-gene matrix and every number downstream of it (`01`, `03`, `04`, `12`, `13`, `17`, `19`, `21`). Copy it out and back in, or exclude it from the wipe. | Pre-existing; surfaced by item 1 |
| `00_setup.R` | Nothing required. It now writes `output/package_versions.csv`, which the per-script `session_info_*.txt` files complement rather than duplicate (they record what each script actually loaded, at the moment it ran). If you would rather have one shared helper than the same 12-line `write_session_info()` in 11 scripts, moving it into `00_setup.R` would remove the duplication — I did not do this because I do not own that file. | Item 19 |

---

## 6. Missing packages

| Package | Status | Impact |
|---|---|---|
| `diptest` | **not installed** | Hartigan dip test for item 16 not computed; substituted with `mclust` BIC + bimodality coefficient (see §4.1) |
| `mclust` 6.1.2 | installed | Used for the FOLR1 mixture fit. Note it must be **attached**, not namespace-loaded: `mclust::Mclust()` calls `mclustBIC()` unqualified. Handled at `13:236`. |
| `jsonlite` | installed | Newly used in `01` to parse `run_info.json` |
| `boot`, `pROC`, `effsize` | not needed | Bootstrap CIs are computed directly (stratified resampling), so no new dependency was added |
| `variancePartition` 1.38.1 | installed but broken in this R (calls `lme4::findbars`, which moved to `reformulas`) | Unchanged from the audit; the documented `lme4` fallback is used and is now labelled with the exact model formulas and the fixed-effect variance definition |
