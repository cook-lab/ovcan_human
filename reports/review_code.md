# Independent code & statistics audit — OvCAN multi-omic Data Descriptor

**Auditor:** external code/reproducibility reviewer (no prior knowledge of the authors or project).
**Scope reviewed:** all 33 analysis/figure/report scripts in `scripts/` (~7,950 lines, read in full), `metadata/*.csv`, `judy_archive/data/` inputs (RNA `abundance.tsv` + `run_info.json`, 22 Mutect2 MAFs, CNVkit `.cns`, `proteomics/*.xlsx`), all 74 files in `output/` and `output/external/`, and the manuscript at `…/scratchpad/manuscript.md`.
**Not read** (by instruction, to preserve independence): `PROJECT_SPEC.md`, `ANALYSIS_PLAN.md`, `ANALYSIS_LOG.md`, `CLAUDE.md`, anything else in `reports/`.
**Read-only:** no project file was modified; scratch work went to the session scratchpad. No pipeline script was executed.

---

## 1. Overall assessment

**The analysis is unusually careful and, on the numbers I could check, unusually accurate — but it is not reproducible as it stands, and three specific defects reach published tables.** I independently recomputed the WES filtering cascade from the 22 raw MAFs (557,392 records → 15,692 PASS → 6,036 coding non-synonymous, matching `wes_mutations_filtered.csv` row-for-row), the OV2295 waterfall (25,914 rows → 493 PASS), the shared-gene count for RNA–protein concordance (8,212), the proteomic block-missingness tiers, the TMT layout/bridge topology, the DepMap overlap and reciprocal-best matching, and ~55 further manuscript numbers; essentially all of them are exactly as reported. I also stress-tested the two load-bearing statistical arguments (subtype-vs-site, patient collapse) and both survive, though one of them is quantified with the wrong statistic.

The blocking problems are: (i) `output/wes_mutations_filtered.csv` — the input to five downstream scripts — is written by no script (the only `write_csv` is commented out), so a clean re-run breaks; (ii) the authentication script consumes the chrX-contaminated FGA that the paper itself declares an artifact, which flips one line's reported call; (iii) Table S1 lists Tier-3 driver calls unflagged, directly contradicting the paper's "defensible somatic BRCA1/2 = 0". Separately, the Methods misstate the transcriptome reference (the kallisto index predates Ensembl 105) and mislabel a gene-count sum as a read count. I found **no** matrix/annotation misalignment, no join row-explosion, no wrong-margin `apply`, and no multiple-testing abuse — the classes of bug that usually dominate this kind of pipeline are absent.

**Verdict: trustworthy in its numbers, not yet trustworthy as a reproducible artifact.** Fix C1–C3 and M1–M2 before submission; the rest are disclosure and framing.

---

## 2. Pipeline map

| Script | Key inputs | Key outputs | Computes |
|---|---|---|---|
| `00_setup.R` | `metadata/samples.csv` | — | Paths, package **check** (not versions), `SEED <- 1234` |
| `00b_figure_theme.R` | — | — | Manuscript theme, locked subtype/site/tier palettes, `save_fig()` |
| `01_rna_load_qc.R` | 31 × `rna_seq/*/abundance.tsv`, cached `tx2gene` | `rna_txi/dds/vst.rds`, `rna_tpm/counts.csv`, `rna_sample_annotation.csv`, `rna_qc_metrics.csv` | tximport → 39,568 genes; filter ≥10 counts in ≥2 lines → 22,544; VST; QC |
| `02_rna_separation.R` | `rna_vst.rds`, `rna_dds.rds` | `rna_pca.rds`, `rna_silhouette.csv`, `rna_pc_confounder.csv` | PCA (top-2000 HVG), t-SNE, silhouette, univariate PC~subtype/site R² |
| `03_rna_de_signatures.R` | `rna_dds.rds`, `tx2gene` | `rna_de_*.csv`, `rna_signatures_*.csv`, `rna_de_gsea_go.csv`, `rna_de_gsea_recovery.csv` | DESeq2 one-vs-rest ×6; up-signatures; fgsea GO-BP; targeted recovery grading |
| `04_rna_markers_genesets.R` | `rna_tpm.csv`, `tx2gene`, `h.all.v7.4.gmt` | `rna_markers_summary.csv`, `rna_swisnf.csv`, `rna_geneset_scores.csv` | 23-marker panel, top-2-of-6 rule; SWI/SNF ranks; Hallmark singscore |
| `05_proteomics_load_qc.R` | `proteomics/{protein_relative_abundance,peptide_ratio,tmt.layout}.xlsx`, samples | `prot_abundance_matrix.csv`, `prot_matrix.rds`, `prot_qc.csv`, `prot_sample_{qc,annotation}.csv`, `prot_bridge_cor.csv` | Channel→line map (45 data cols), 8,430 × 31 matrix, per-plex missingness, 4 bridge correlations |
| `06_proteomics_separation.R` | `prot_matrix.rds` | `prot_pca.rds`, `prot_silhouette.csv`, `prot_pc_confounder.csv` | Complete-case PCA/t-SNE/silhouette; subtype vs plex vs site R²; QC gate |
| `07_wes_mutations.R` | 22 MAFs + VCFs, `line_family_map.csv` | `wes_driver_tiers.csv`, `wes_driver_freq_{patient,by_subtype}.csv`; **reads** `wes_mutations_filtered.csv` | FILTER==PASS + pop-AF filter, VAF join, 3-tier somatic confidence, patient collapse, oncoprint |
| `08_wes_cnv.R` | 23 CNVkit `.cns`, family map | `wes_cnv_segments.csv`, `wes_cnv_fga.csv`, `wes_cnv_arm_freq_patient.csv` | Autosome median-centring, chrX diagnostic, FGA, arm calls per patient, CNV heatmap |
| `09_wes_hrd.R` | CNVkit headers, `commands.txt` | `wes_hrd_feasibility.md` | Feasibility gate: no allele-specific CN → HRD not computable (no score fabricated) |
| `10_authentication.R` | `rna_tpm`, `prot_abundance_matrix`, `wes_mutations_filtered`, `wes_cnv_fga` | `auth_swisnf_panel.csv`, `auth_perline_table.csv` | Multi-omic SWI/SNF panel; per-line expression/genomics consistency; flags |
| `11_mucinous_authenticity.R` | same three layers | `auth_mucinous.csv` | Ovarian-vs-GI marker index and verdicts for 3 MC lines |
| `12_rna_protein_concordance.R` | `rna_tpm`, `prot_abundance_matrix` | `integ_rnaprot_cor{,_summary}.csv` | Per-gene (n≥10 lines) and per-line Spearman/Pearson + threshold sensitivity |
| `13_adc_atlas.R` | `rna_tpm`, `prot_abundance_matrix` | `adc_expression.csv`, `adc_subtype_summary.csv` | 9-target ADC atlas, subtype means, expected-association check |
| `14_hgs_heterogeneity.R` | `rna_geneset_scores.csv`, `rna_tpm` | `hgs_heterogeneity.csv`, `hgs_hallmark_cluster_means.csv` | ward.D2 k=3 on Hallmark z (15 HGS); PROGENy; theme scores |
| `15_patient_family_map.R` | `metadata/samples.csv` **only** | `metadata/line_family_map.csv` | 42 lines → 34 patients; 5 families; deterministic representative |
| `16_wes_signatures_msi.R` | `wes_mutations_filtered.csv`, 22 raw MAFs, hg38 BSgenome | `wes_mutation_load.csv`, `wes_msi_mmr.csv`, `wes_sbs_context.csv` | Load/indel/Ts-Tv, MMR & POLE scan, SBS-96, COSMIC v3.2 cosine, MSI call |
| `17_variance_confounders.R` | `rna_vst/pca/dds`, `prot_matrix.rds`, family map | `rna_pc_confounder_joint.csv`, `rna/prot_variancepartition.csv`, `rna_passage_check.csv`, `rna_marker_effectsizes.csv` | Commonality analysis; per-feature lme4 REML decomposition; passage check; Cohen's d + AUC on 28 patient reps |
| `18_external_benchmarking.R` | `output/external/{Model,depmap_expr_ovarian,…}.csv`, Cellosaurus JSON, `rna_tpm` | `external_ccle_concordance.csv`, `consensusov_calls.csv`, `cellosaurus_str_status.csv` | DepMap overlap + expression identity (2000 HVG), driver cross-check, ConsensusOV **parse**, STR status |
| `19_proteomics_dynamic_range.R` | `prot_abundance_matrix`, `rna_tpm`, `integ_rnaprot_cor`, `adc_expression`, layout xlsx | `prot_dynamic_range.csv`, `prot_block_missingness.csv` | n-matched per-gene RNA vs protein spread; block missingness; bridge noise floor |
| `20_supplement_table.R` | family map + 8 `output/*.csv` | `supplement_per_line.csv` (Table S1) | Pure join, 42 × 41 |
| `30_fig1_overview.R` | `line_family_map.csv` | `fig1.{pdf,png}` | Workflow schematic, assay coverage matrix, counts |
| `31_fig4_genomics.R` | `wes_driver_tiers`, `wes_mutation_load`, `wes_cnv_arm_freq_patient`, `external_ccle_concordance` | `fig4.*`, `assets/f_wes_oncoplot.png`, `assets/f_external_concordance.png` | Waterfall (hard-coded), ggplot oncoprint, arm frequency, DepMap heatmap |
| `32_fig5_rare.R` | `wes_mutation_load`, `wes_msi_mmr`, `wes_sbs_context`, `auth_mucinous` | `fig5.*`, `assets/f_wes_hypermutation.png`, `assets/f_auth_mucinous.png` | Load bar, SBS-96, COSMIC cosine, mucinous markers |
| `33_supp_genomics.R` | `wes_cnv_segments`, `consensusov_calls`, `hgs_heterogeneity` | `figs5.*`, `figs6.*`, `assets/f_wes_cnv.png` | Per-line CNV heatmap; ConsensusOV vs intrinsic strata |
| `34_fig2_qc.R` | `rna_qc_metrics`, `prot_bridge_cor`, `prot_dynamic_range`, `prot_block_missingness`, abundance xlsx | `fig2.*`, 3 assets | QC scatter, bridge scatter (with assert vs stored r), compression, missingness |
| `35_fig3_biology.R` | `rna_pca.rds`, `integ_rnaprot_cor*`, `*_variancepartition`, `rna_markers_summary`, `rna_pc_confounder_joint` | `fig3.*`, 5 assets | PCA ×2, concordance density, variance decomposition, marker heatmap |
| `36_fig6_adc.R` | `adc_expression`, `adc_subtype_summary` | `fig6.*`, `assets/f_adc.png` | Rebranded RNA-over-protein ADC heatmap |
| `37_supp_rnaprot.R` | `prot_pca.rds`, `rna_passage_check`, `rna_vst.rds`, `prot_block_missingness`, `rna_marker_effectsizes`, `hgs_heterogeneity` | `figs{1,2,3,4,7,8}.*`, 4 assets | Protein PCA, passage, sample correlation, presence patterns, effect sizes, HGS strata |
| `build_report.py` | `reports/01_…md`, **`reports/01_…html`** (own prior output) | `reports/01_…html` | HTML presentation build |
| `sandbox/{explore_hypermutation,probe_sbs}.R` | filtered MAF / raw MAF | console only | Exploration (consistent with `16`) |

**Dependency-order finding:** `15_patient_family_map.R` writes `metadata/line_family_map.csv`, which `07` and `08` consume. A numeric top-to-bottom run works only because that file is committed.

---

## 3. Critical issues

### C1 — `output/wes_mutations_filtered.csv` is generated by no script; a clean re-run breaks 5 scripts
**`scripts/07_wes_mutations.R:424-433`**
```r
# readr::write_csv(as_tibble(filtered_out), file.path(OUT, "wes_mutations_filtered.csv"))
stopifnot("recomputed retained-variant count drifted from the on-disk table" =
            nrow(filtered_out) == nrow(data.table::fread(file.path(OUT, "wes_mutations_filtered.csv"))))
```
*What it does:* the only write of this file is commented out; the script instead **reads** the file to assert its own row count.
*Should do:* write the file it produces.
*Consequence:* from a clean checkout with an empty `output/`, script 07 aborts at line 430, and `10_authentication.R:101`, `11_mucinous_authenticity.R:100`, `16_wes_signatures_msi.R:75`, `18_external_benchmarking.R:240` all fail. Every mutation-derived number in the manuscript (tiers, TP53, CDK12, BRCA2, MSI/TOV21G, SWI/SNF, mucinous, DepMap driver cross-check, Table S1 drivers) becomes non-reproducible. The on-disk file is dated 23 Jul 14:50, before the current script (22:43) — it is the product of an older version.
*Verification:* **CONFIRMED.** I recomputed the cascade independently from the 22 raw MAFs (PASS ∧ pop-AF ≤ 0.001 ∧ maftools non-synonymous set): 557,392 raw → 15,692 PASS → 15,609 → **6,036** coding non-synonymous, exactly the 6,036 rows on disk; per-line median 206.5, range 133–1,416. So the archived table *is* currently correct and the guard would pass — the defect is purely that nothing regenerates it.
*Fix:* uncomment line 429. Keep the assertion, but compare a content hash, not just `nrow` (see Minor).

### C2 — `10_authentication.R` uses the chrX-contaminated FGA the paper calls an artifact, and it changes a reported count
**`scripts/10_authentication.R:309-311`** (thresholds at `:386-388`, `:408`)
```r
fga_l <- fga %>% transmute(cell_line, FGA = round(fga_0.2, 3),
  cnv_instability = cut(fga_0.2, c(-Inf, 0.15, 0.30, Inf), …))
```
*What it does:* reads `fga_0.2` — the **legacy, chrX-inclusive** column — from `wes_cnv_fga.csv`. `08_wes_cnv.R` was explicitly revised to make `fga_auto_0.2` the headline metric because "chrX is a pooled-normal sex-composition artifact" (`08:32-36`), and the manuscript reports autosome-restricted FGA throughout.
*Should do:* use `fga_auto_0.2`.
*Consequence:* **OV90 is the one line that crosses the 0.30 HGS threshold** — `fga_0.2 = 0.3084` vs `fga_auto_0.2 = 0.2689`. With the corrected metric, OV90's `genomics_consistent` goes `consistent` → `partial` and `cnv_instability` `high` → `intermediate`. The manuscript's "**20 genomics-consistent**" becomes 19 (and one more "partial"). OV90 is one of the specifically discussed flagged lines. Table S1 is already internally inconsistent because `20_supplement_table.R:66` correctly uses `fga_auto_0.2` for `fga_autosome` while inheriting `genomics_consistent` from the with-chrX rule: OV90's row reads `fga_autosome = 0.269` next to `genomics_consistent = "consistent"` for a rule stated as FGA > 0.30. Also every `FGA` value in `auth_perline_table.csv` is the artifact-inclusive one.
*Verification:* **CONFIRMED** — I compared both columns across all 23 CNV lines; OV90 is the only 0.30-threshold flip and the only `cnv_instability` category change.
*Fix:* switch to `fga_auto_0.2` in `10:309-311`; re-run `10` then `20`; update the "20 genomics-consistent" figure.

### C3 — Table S1 lists Tier-3 driver calls unflagged, contradicting "defensible somatic BRCA1/2 is zero"
**`scripts/10_authentication.R:306-308`** → **`scripts/20_supplement_table.R:69-73`**
```r
drivers <- wm %>% filter(is_driver) %>% group_by(cell_line) %>%
  summarise(key_drivers = paste(sort(unique(Hugo_Symbol)), collapse = ", "))
```
*What it does:* `key_drivers` is built from all retained driver calls with **no tier filter**, and `supplement_per_line.csv` (Table S1) carries it with **no tier column at all**.
*Should do:* restrict to Tier 1–2, or add per-gene tier annotation.
*Consequence:* Technical Validation §3 states "Defensible somatic BRCA1/2 is zero: both candidate BRCA2 calls are Tier 3 and excluded", and Usage Notes tell reusers Table S1 "is the recommended entry point for reuse" providing "tiered driver calls". But Table S1 shows `key_drivers = "BRCA2"` as TOV81D's **only** driver, and `"BRCA2, CDK12, TP53"` for TOV3133D — exactly the two calls the paper says are not defensible. Nine (line, gene) pairs are Tier-3-only and appear unqualified: OV2295 *SMARCA2*, OV2295-R2 *ARID1A*, TOV1369 *NF1*, TOV21G *CTNNB1* + *SMARCA2*, TOV2929D *CDK12* + *SMARCA2*, TOV3133D *BRCA2*, TOV81D *BRCA2*. The hypermutator passenger-risk `context` column that exists in `wes_driver_tiers.csv` is also not carried into Table S1.
*Verification:* **CONFIRMED** by joining `wes_driver_tiers.csv` tier minima to `supplement_per_line.csv`.
*Fix:* in `20_supplement_table.R`, join `wes_driver_tiers.csv` and emit `drivers_tier12` and `drivers_tier3` (or `key_drivers` annotated as `GENE(T3)`); Table S1 is the artifact reusers will actually use.

---

## 4. Major issues

### M1 — Methods misstate the transcriptome reference: the kallisto index is **not** Ensembl 105
Manuscript Methods: "pseudoaligned with kallisto v0.46 against a transcriptome index built from **Ensembl release 105**".
`judy_archive/data/rna_seq/*/run_info.json` records `kallisto_version 0.46.0`, `index_version 10`, `n_targets 185299`, `start_time Wed Oct 20 23:13:59 2021`, index `…/kallisto_index/Homo_sapiens_GRCh38` — i.e. the index was built **before** release 105 existed (Dec 2021). `01_rna_load_qc.R:36-40` says so: "Index built Oct-2021 (run_info): current Ensembl release then was 104 … release 105 (Dec 2021) is the nearest available pinned archive".
*Consequence:* the deposited quantification's actual annotation is misstated, and the release mismatch has a measurable, undocumented cost: **3,529 of 185,299 index targets (1.9%) are absent from the release-105 `tx2gene` and are silently dropped by tximport, carrying 2.07% of TPM and 1.64% of estimated counts.** Nothing logs or asserts this. A reuser rebuilding against 105 will not reproduce the matrices.
*Verification:* **CONFIRMED** (parsed one `abundance.tsv` against `output/tx2gene_ensembl_rel105.csv`; index version/date consistent across all 31 samples).
*Fix:* state that the index derives from Ensembl 104 and that the transcript→gene map is pinned to 105 as the nearest live archive, with the 1.9% / 2.1% loss quantified; add a `stopifnot` on the unmapped fraction in `01`.

### M2 — "median library size of 56.5M reads" is a sum of gene-level estimated counts, not a read count
**`scripts/01_rna_load_qc.R:92`** — `lib_size = colSums(counts)`, where `counts` are tximport gene-level estimated counts **after** the ≥10-in-≥2-lines filter. Median = 56,497,657.
*Consequence:* the true median sequenced fragment count (`n_processed`, `kallisto_qc.tsv`, 31 generated lines) is **64.3M** (range 45.1–97.0M); for paired-end that is 64.3M pairs ≈ 128.6M reads. A Data Records / Technical Validation number readers will take as sequencing depth is wrong by ~12% (or ~2× if "reads" is read literally).
*Verification:* **CONFIRMED.**
*Fix:* report `n_processed` and rename the column `assigned_gene_counts`.

### M3 — The "clean cross-site control" (site explains 4–6% within clear cell) is computed on the wrong axes and is underpowered
**`scripts/17_variance_confounders.R:143-151`** subsets 7 clear-cell rows out of the PCA fitted to **all 31 lines**, then regresses those global PC scores on site. Those PCs are defined by *between-histotype* variance, so the within-CC slice is a residual and site's R² on it is structurally small. Output: PC1 0.044, PC2 0.058, PC3 0.057 → the manuscript's "4–6%".
*Should do:* decompose within-group variance within the group.
*Consequence:* recomputing the PCA **within clear cell only** (n=7: 2 CHUM, 5 BC) gives site R² = **0.157 / 0.126 / 0.273** on PC1–PC3 (p = 0.38 / 0.44 / 0.23); per-gene within-CC R²(site) has median **0.194**, mean 0.267 over 22,387 genes. Note the null expectation for a 1-df/5-df split is E[R²] ≈ 1/6 = 0.167, so 0.194 genuinely indicates no detectable site effect — but the correct conclusion is "**underpowered: with 2 vs 5 lines a within-clear-cell site effect is indistinguishable from chance**", not "site explains only 4–6%". As written the number reads 3–5× more reassuring than the data support.
*Verification:* **CONFIRMED** (both statistics recomputed from `rna_vst.rds`).
*Fix:* reword; the commonality result (unique-site 0.19% of PC1) is the strong leg and stands on its own.

### M4 — Undisclosed pseudoreplication throughout the expression arm
Patient collapse is applied to genomic frequencies (`07`, `08`) and to marker effect sizes (`17D`, 28 reps) — but **not** to PCA/t-SNE/silhouette (`02`), one-vs-rest DE and GSEA (`03`), Hallmark scoring (`04`), RNA–protein concordance (`12`), ADC subtype means (`13`), or HGS clustering (`14`). The RNA set is 31 lines / 28 patients; HGS is **15 lines / 12 patients** (pairs 1369, 2295, 3133). `03_rna_de_signatures.R:70-82` gives DESeq2 15 HGS "replicates", 3 of which are duplicate genomes.
*Consequence:* the HGS Wald/BH statistics — and therefore the GO recovery in Technical Validation §2 — are anti-conservative by an undisclosed amount. Methods correctly scope collapse to "all genomic frequency statistics" but the expression arm is never flagged.
*Verification:* **CONFIRMED, and quantified.** Recomputing PCA + silhouette on the 28 patient representatives: HGS silhouette 0.162 → **0.142**; overall 0.220 → 0.208; PC1 20.7% → 19.9%; PC1 R²(subtype) 0.735 → **0.792**. So the separation conclusion is robust and the direction is if anything favourable — the issue is disclosure, and the DE/GSEA p-values, which I did not re-derive.
*Fix:* state the limitation; run the HGS contrast on the 28 reps as a one-paragraph sensitivity analysis.

### M5 — ConsensusOV calls are parsed from free-text notes, not computed
Manuscript Methods: "ConsensusOV HGSC subtype calls **were generated** for the 15 HGSC RNA lines … (script `10`)".
**`scripts/18_external_benchmarking.R:254`**
```r
consensusov_call = str_match(notes, "ConsensusOV:([A-Za-z]+)")[, 2]
```
No script loads or runs `consensusOV`; `10_authentication.R` never mentions it. The calls are string-extracted from `metadata/samples.csv` `notes`, inherited from the archived analysis with no recorded classifier version or settings.
*Consequence:* Table S5 / Fig S6 / the §2 caveat paragraph rest on a result no code in this repository can reproduce; the script attribution (10) is also wrong (18).
*Verification:* **CONFIRMED** (counts DIF 6 / MES 4 / IMR 3 / PRO 2 = 15 match the manuscript, but they are transcriptions).
*Fix:* either run consensusOV in the pipeline, or say explicitly that the calls are inherited and unreproducible here; correct the attribution.

### M6 — TOV21G's EPCAM truncating variant is dropped from the MSI narrative
Manuscript Usage Notes: "**No MMR-enzyme coding mutation is present**, consistent with the epigenetic (MLH1-methylation) mechanism".
**`scripts/16_wes_signatures_msi.R:59`** defines `MMR_INDIRECT <- c("EPCAM")  # 3' deletions silence MSH2 (indirect)`, and `output/wes_msi_mmr.csv` records for TOV21G: `epcam_hit = TRUE`, `mmr_genes_hit = "EPCAM(truncating)"` — the **only** MMR-panel hit anywhere in the panel. `32_fig5_rare.R`'s panel-B subtitle even says "TOV21G: only EPCAM".
*Consequence:* the pipeline's own MMR panel surfaces a truncating hit in the one gene it flags as an indirect MMR mechanism, and the text says none is present. Whether the call is real (tumour-only, so germline cannot be excluded) or artefactual, dropping it silently is not defensible in a paper that otherwise reports its negatives.
*Verification:* **CONFIRMED.**
*Fix:* report the EPCAM call with its tumour-only caveat.

### M7 — Mutational-signature evidence is cosine screening, not signature fitting, on a germline-contaminated substrate
**`scripts/16_wes_signatures_msi.R:265-298`**: `mut_matrix()` → `cos_sim_matrix()` against COSMIC v3.2, then `sbs_call` thresholds a *group maximum* at an unjustified 0.75. There is no `fit_to_signatures()` refit, no bootstrap/uncertainty, and no attempt to remove the germline component. Substrate = exome-wide tumour-only PASS ∧ popAF ≤ 0.001 SNVs: **303–2,417 per line** (TOV21G 2,417); the paper itself puts residual germline at ~206 coding variants/line, so a comparable germline fraction sits in the exome-wide set.
*Consequence:* "SBS6 cosine 0.88" is a similarity to a single reference spectrum, not an exposure estimate, and Methods describe it as comparison "to COSMIC SBS reference signatures by cosine similarity" (accurate) while the Usage Notes read as attribution.
*Verification:* **CONFIRMED, and the direction survives:** TOV21G cos(MMR-d) 0.877 / cos(POLE) 0.315 / cos(clock) 0.651; all 21 other lines sit at cos(MMR-d) 0.51–0.70 with cos(clock) > cos(MMR-d). TOV21G is a genuine, well-separated outlier. The calibrated "candidate … pending MMR IHC/MSI-PCR" framing largely covers this.
*Fix:* say "cosine-similarity screening" in the Usage Notes too; add a refit (`fit_to_signatures_strict`) or state why not; justify the 0.75 cut or drop it in favour of the reported margin.

### M8 — The variance decomposition's `patient` term is near-unidentifiable, so "patient is the largest structured term" is structural
**`scripts/17_variance_confounders.R:196-198, 207-218`**: `y ~ (1|subtype) + (1|source_site) + (1|patient)` on 31 RNA lines with **28 patient levels** (only 3 patients replicated) and **3 site levels, one with n=1** (BIN67). The patient variance component and the residual are separable only from those 3 replicate pairs; the site component is estimated from 3 groups.
*Consequence:* the manuscript's "with line/patient identity the largest structured term, as expected for distinct cell lines" (RNA 27.3%, protein 19.8%) restates the model's design rather than reporting a finding, and the site components (RNA 3.5%, protein 0.0%) carry very wide implicit uncertainty. The subtype ≥ site ordering is still directionally supportable. Additionally, **2 of 22,544 genes are silently dropped** (`vp_lmer`'s `keep <- !vapply(rows, is.null, …)`; `n_features = 22,542`) with no warning, and `35_fig3_biology.R:116` hard-codes "22,544 genes" in the Fig 3D facet label.
*Verification:* **CONFIRMED** from `rna/prot_variancepartition.csv` and the family map.
*Fix:* report the patient term as "line identity (nominally patient; not separable from residual at n=3 replicated patients)"; note the 3-level site caveat; log the 2 dropped genes.

### M9 — "does not translate into a difference in the number of genes detected" is an untested absence claim, and the medians do differ
**`scripts/01_rna_load_qc.R:101-103`** prints per-site medians only; no test is performed anywhere.
*Verification:* **CONFIRMED** — Mes-Masson (n=19) pseudoalign 92.2% / 19,817 genes; Huntsman (n=11) 88.1% / **20,765** genes; Huntsman/Vanderhyden (n=1) 87.9% / 20,119. The lower-alignment site detects ~950 **more** genes.
*Consequence:* an absence-of-effect conclusion drawn from a descriptive summary, in the direction opposite to the framing.
*Fix:* run a Wilcoxon test (19 vs 11) or reword to "the lower-alignment site does not detect fewer genes".

---

## 5. Minor issues and code-quality notes

**Counting / hygiene**
- `12:98`, `19:92` drop the single `NA`-symbol protein → 8,429 rows, so the manuscript's "Of **8,430** proteins, 6,856 (81.3%) … 1,573 (18.7%)" sums to 8,429 (8,430 − 6,856 = 1,574). Off by one.
- **71 of 8,430 proteins have `present_n_plex == 0`** — `NA` in all 31 analysis lines — and remain in the deposited `prot_abundance_matrix.csv`. `05:149` computes `n_any = 8,359` but never uses it.
- `05:113` `make.unique(abund$Symbol)` yields two rows literally named `NA.1` / `NA.2` (verified) plus 33 `SYMBOL.N` rows (5 × HLA-A, HLA-B, HLA-DRB1, CHTF8, POLR1D…). Symbol-keyed analyses (`10:96`, `11:63`, `12:98`, `13:87`, `19:92`) keep only base names, so those rows are orphaned; `06:137` `mat[match(mk$gene, feat$symbol), ]` silently takes the first row for a duplicated symbol.
- Peptide table has 8,617 accessions vs 8,430 abundance rows (`05:167-169` prints both, never reconciles).
- `07:112-121` rebuilds `HGVSp_Short` (empty in the MAFs) as `p.<ref><pos><alt>`, producing non-HGVS strings: `p.L639X`, `p.N210X`, `p.H214X`, `p.-156-157P`. The manuscript quotes `p.L639X` as a variant identifier; `X` for a frameshift will not match ClinVar/COSMIC. (`p.Q192*` is correct.)

**Weak or asserted claims**
- `14:230` and `build_report.py:41` print "same-patient families mostly co-cluster" as a **fixed string**. Verified: of the 3 multi-line families with 2 RNA lines each, 2295 and 3133 co-cluster; **1369 splits (clusters 1 vs 3) and still splits at k=2**. 2 of 3.
- `14`: PROGENy is called "orthogonal corroboration" (manuscript: "corroborated orthogonally by PROGENy") but is computed from the *same* log2-TPM matrix with different gene sets; the Hallmark-vs-PROGENy agreement is a printed contingency table with no statistic (no ARI/Rand index).
- `10:367-369`: `expression_consistent = "consistent"` for HGS/MMMT means only **absence** of a competing CC/MC program. Of the 26 "consistent": 15 HGS + 1 MMMT by absence, 2 SCCOHT from the same SWI/SNF evidence used to make the reclassification (circular), and only 5 CC + 3 MC on a positive lineage program. "26 expression-consistent" needs that breakdown.
- `10:129-140`, `:182-184`: `SMARCA2_loss` uses truncating status irrespective of tier, so Tier-3 calls (OV2295, TOV21G *SMARCA2*) feed `swisnf_deficient`.
- `11:163-166`: TOV2414's verdict hard-codes "SATB2-" from Sauriol 2020. In **these** data TOV2414's SATB2 z = **+0.17**, *higher* than VOA8762's (−0.59). The manuscript's "TOV2414 … (KRT7+/PAX8+/MUC5AC+/**SATB2−**…)" reads as four in-house measurements; three are (KRT7 z 2.58, PAX8 z 0.28, MUC5AC z 4.25), SATB2− is not.
- `11:149`: `ovarian_index = mean(KRT7_z, PAX8_z) − SATB2_z` is an unweighted ad-hoc composite mixing a protein z with RNA z's; no calibration, and the verdict thresholds (0.5, 1.0) are unjustified.
- `08:211` patient-level FGA = mean of a patient's lines; for CC the "median of 2 patients" is the mean of 0.671 and 0.071 (the manuscript does flag this).
- `08` reports FGA under "endometrioid 0.23 (n=1)" for TOV112D — the line the same paper reclassifies away from endometrioid.
- Abstract "TP53 mutation in all 11 high-grade serous patients": 11 is the WES subset of the resource's 16 HGS patients. Add the qualifier.
- `build_report.py:42` caption: passage "is 83% collinear with site and **not an independent structure driver**". `pass_site_r2 = 0.835` ✅, but `rna_passage_check.csv` gives `partial_r2_passage_after_site = 0.148` for PC1 (R²(site)=0.227, R²(passage)=0.078, R²(joint)=0.375 — classic suppression): passage adds ~15% of PC1 variance *beyond* site. The caption contradicts its own output. Passage is absent from the manuscript entirely.
- `17:354` writes 22 unadjusted `wilcox.test` p-values; the minimum attainable p for the n=2 groups is ~0.005.
- `17:326-331` (`get_vst`: single ENSG by max mean VST) and `04:54` (symbol-summed TPM) quantify the same markers differently, yet `rna_marker_effectsizes.csv` joins `lands_right` from the latter onto effect sizes from the former.
- `18:112-129`: our side collapses to Entrez by **summing linear TPM then logging**; DepMap by **averaging log2(TPM+1)** (`collapse_mean`). HVGs are also chosen on the *combined* 31+67 matrix, so platform-divergent genes are preferentially selected — this inflates the absolute Spearman (0.74–0.88) for self and non-self pairs alike, which is why the manuscript's "specificity is the signal" framing is the right one.
- `18:288`: the Cellosaurus parser falls back to the search top hit when there is no exact name match, which could attach a wrong accession. **Benign here** — all 30 hits are `high (exact name match)`, 12 `not found`, 0 low-confidence.
- `13:130-137`: RNA and protein subtype means are over different line sets (VOA6861 RNA-only, VOA14993 protein-only) with no note.

**Code quality**
- Hard-coded numbers in figure scripts (all currently correct, all verified, all liable to silent desynchronisation): `31:150-152` `25914`/`493`; `31:162` `ATM 82%->9%, ATR 77%->5%`; `31:109` `"1416"`; `32:36` `1416`, `6.9x`; `35:116-117` `22,544`/`6,856`.
- `31:87-90`: `ylab_text` assigned twice; the first assignment is dead.
- `07:180`: `n_hgs` assigned, never used.
- `07:430`: drift guard compares only `nrow`; a logic change preserving 6,036 rows would pass silently. Use a hash.
- `07:106-108`, `16:237-239`: `num_or0()` coerces multi-allelic AF strings (`"a,b"`) to `NA → 0`, i.e. "not common". The comment acknowledges this and relies on Mutect2's `common_variant` flag as a backstop.
- `07:330-331`, `16:112-121`: the hypermutator flag uses two different definitions (`>3× median` in `07`, `robust_z > 5 & fold > 3` in `16`). They agree on TOV21G, but there is one rule, two implementations.
- `10:378-379` `has_drv()` does regex matching over a comma-joined string when the structured table is in hand.
- `build_report.py:88-91` reads its CSS out of its **own previous output** (`SRC_HTML == OUT_HTML`); the HTML build cannot bootstrap from a clean tree.
- `reports/assets/f_auth_swisnf.png` is captioned in `build_report.py:38` but **no script generates it** (`10` writes only `figs/10_swisnf_panel.pdf`) — an orphaned figure in the HTML report.
- Nothing in the pipeline reads `data/`; it is empty apart from an empty `data/external/`. All source data live in `judy_archive/`, which is not part of the documented project layout.
- Manuscript figure numbering: **Fig. S8 is cited for two different things** — WES coverage metrics ("Fig. S8 is contingent on this") and within-HGSC strata. **Figs S1, S2, S3, S5 are generated (`37`, `33`) but never cited** in the manuscript text.
- Code Availability says `variancePartition` "is unavailable in this R build"; it **is** installed (1.38.1) and fails at run time — `17:162-166` describes it correctly.
- Manuscript attributes the marker effect sizes to "script `04`"; they are produced by `17` (`rna_marker_effectsizes.csv`). §4's "no BC Cancer VOA line exists publicly" and the ConsensusOV surfacing are likewise attributed to the wrong scripts (`10` vs `18`).

---

## 6. Reproducibility audit

**Seeds.** `set.seed(SEED)` in 20 of 30 R scripts. The 10 without (`00b`, `01`, `05`, `09`, `15`, `20`, `30`, `31`, `32`, `33`) contain no stochastic step — checked. Rtsne (`02:43-44`, `06:78-79`) is seeded; `progeny(perm = 1)` is deterministic; `lmer`/`bobyqa` under `mclapply` is deterministic. **One residual risk:** `fgsea()` (`03:197`) is Monte-Carlo and inherits the RNG state set at `03:40`, so a full top-to-bottom run of `03` is reproducible but a partial re-run gives slightly different GO p-values. Move `set.seed` inside `run_gsea()`.

**Paths.** All scripts resolve from `Sys.getenv("OVCAN_PROJ", unset = "/Users/dpcook/Analysis/ovcan_human")` — portable with the env var, personal default otherwise. `judy_archive/` is hard-wired as `DATA` in `00_setup.R:18` but appears nowhere in the documented structure.

**Versions.** **No `renv.lock`, no `sessionInfo()` call anywhere.** `check_pkgs()` (`00_setup.R:41-55`) only warns about missing packages. Both Methods and Code Availability claim `00_setup.R` "records paths, **package versions**, and the seed" — it does not. Observed environment: R 4.5.2 ✅, DESeq2 1.48.2, tximport 1.36.1, fgsea 1.34.2, singscore 1.28.1, maftools 2.24.0, MutationalPatterns 3.18.0, BSgenome.Hsapiens.UCSC.hg38 1.4.5, progeny 1.30.0, lme4 1.1.38, variancePartition 1.38.1 (installed, broken), scarHRD/copynumber absent ✅ (consistent with `09`).

**Execution order.** Numeric order is *not* the dependency order: `15` writes `metadata/line_family_map.csv`, consumed by `07`/`08`. `20` correctly runs last. `30`–`37` correctly consume only materialised outputs. No script consumes an output produced later except the `15` inversion.

**Orphaned / undocumented artifacts.** `output/wes_mutations_filtered.csv` (C1, written by nothing). `reports/assets/f_auth_swisnf.png` (generated by nothing). `output/external/` (528 MB DepMap + Cellosaurus JSON) is required by `18` but downloaded by no committed script — the header describes it as "downloaded to output/external/ by this workstream", with no code.

**Silent-failure risks.** tximport drops 1.9% of index targets without assertion (M1). `vp_lmer` drops failing genes with `keep <- !vapply(rows, is.null, …)` and no count (M8) — and `do.call(rbind, …)` on named vectors would misalign if `VarCorr` ever reordered components (it does not for a fixed formula, but the code relies on it). `09:83-89` wraps the scarHRD install in `tryCatch` and records the message — good practice. `10:150` `gv()` returns `NA_real_` for absent lines rather than erroring — intentional. `19:220-221` `apply(…, 1, paste)` on `present_by_plex` — verified correct.

**Verdict on a clean re-run.** **It would not complete.** Blockers, in order: (1) `wes_mutations_filtered.csv` cannot be produced (C1) → `07` aborts and `10`, `11`, `16`, `18` fail; (2) `output/external/` must be fetched by hand; (3) `metadata/line_family_map.csv` must pre-exist or `15` must be run out of numeric order; (4) `build_report.py` needs its own prior HTML output; (5) no lockfile, so Bioconductor drift is unbounded. With those five fixed, the numbers themselves should reproduce: every number I recomputed from raw inputs matched.

---

## 7. Code ↔ claim traceability

Legend: ✅ traced and verified · ⚠ traceable but qualified/mis-attributed · ❌ mismatch or untraceable.

| # | Manuscript claim | Source | Status |
|---|---|---|---|
| 1 | 42 models, 34 patients, 7 histotypes | `15`; `line_family_map.csv` | ✅ 42 rows / 34 patient_id |
| 2 | HGSC 16pt/24 lines, CC 8/8, MC 3/3, EC 2/2, MMMT 2/2, SCCOHT 2/2, LGS 1/1 | `15`, `30:150-153` | ✅ exact |
| 3 | RNA 31, TMT 31, WES 23, 13 tri-omic | `01:29`, `05:58`, `08:97`, `20:85-87` | ✅ 31/31/23/13 |
| 4 | Single RNA run `NS.1676.003`; single 5-set TMT design | `run_info.json` `call`; `tmt.layout.xlsx` (55 rows = 5×11) | ✅ |
| 5 | kallisto v0.46 | `run_info.json` `0.46.0` | ✅ |
| 6 | index built from **Ensembl release 105** | `run_info.json` 2021-10-20; `01:36-40` | ❌ **M1** (index predates 105) |
| 7 | 39,568 genes from tximport | `rna_txi.rds` | ✅ 39,568 × 31 |
| 8 | 22,544 genes after ≥10-in-≥2 filter | `01:70-74`; `rna_vst.rds` | ✅ 22,544 |
| 9 | 8,430 proteins / 146,830 peptides | `05`; `prot_qc.csv`; `peptide_ratio.xlsx` | ✅ 8,430 / 146,830 |
| 10 | 6,856 (81.3%) complete; 1,573 (18.7%) set-conditional | `19`; `prot_block_missingness.csv` | ⚠ correct but ÷8,429; 8,430−6,856 = 1,574 |
| 11 | median 13 of 31 lines lost per affected protein | `19:259-260` | ✅ 13 |
| 12 | ≥50%-presence filter retains 7,734 | `05:148`; `prot_qc.csv` | ✅ 7,734 |
| 13 | ch1 PIS, ch10 daisy-chain bridge (1↔2↔3↔4↔5), ch9 empty | `tmt.layout.xlsx`; `05:84-90`; `19:266` | ✅ 45 data cols; 4 consecutive links |
| 14 | Missingness per-set, all-or-nothing | `prot_sample_qc.csv` | ✅ exactly 1 distinct `n_detected` per plex |
| 15 | 23 CNV / 22 MAF; TOV3121D no MAF | `07:83`, `08:118`; `line_family_map.csv` | ✅ |
| 16 | OV2295 MAF 25,914 rows → 493 PASS | `07:30` comment; `31:150-152` hard-coded | ✅ recomputed from the MAF: 25,914 / 493 |
| 17 | ~206 coding variants per line | `07:174-176` | ✅ median 206.5 (133–1,416) |
| 18 | 50 driver calls; Tier1 27 / Tier2 11 / Tier3 12 | `wes_driver_tiers.csv` | ✅ 50; 27/11/12 |
| 19 | MAF header says GRCh37 but calls are GRCh38 | `16:22-31`; MAF `NCBI_Build` | ✅ field literally reads `GRCh37` |
| 20 | chrX median log2 −1.23 to +0.83, sign inverting | `08:163-174`; `wes_cnv_fga.csv` | ✅ −1.228 to +0.826 |
| 21 | Two genomes autosome FGA > 0.7 (TOV3121D, TOV2929D) | `08:201` | ✅ |
| 22 | 5 public normals PRJNA339046, kit unconfirmed | `08:22-31`, `09:48-50` (`commands.txt`) | ⚠ from source-file inspection, not a written output |
| 23 | R 4.5.2; fixed seed 1234; `00_setup.R` records **package versions** | `00_setup.R` | ⚠ R/seed ✅; **versions ❌** (no lockfile, no `sessionInfo`) |
| 24 | Pseudoalignment median 91.1%, range 85.8–93.1% | `rna_qc_metrics.csv` | ✅ |
| 25 | Median library size **56.5M reads** | `01:92` | ❌ **M2** (gene counts; true median 64.3M fragments) |
| 26 | ~20,119 genes detected per line | `rna_qc_metrics.csv` | ✅ 20,119 |
| 27 | Mes-Masson 92.2% vs Huntsman 88.1%, no detection difference | `01:101-103` | ⚠ **M9** (medians differ +948 genes, untested) |
| 28 | Bridge Pearson 0.991–0.994; median CV 5.3% | `prot_bridge_cor.csv`, `prot_qc.csv` | ✅ 0.9914–0.9943; 5.31% |
| 29 | RNA PC1 20.7%, PC2 10.4% | `rna_pc_confounder_joint.csv` | ✅ |
| 30 | Commonality PC1: unique-subtype 42%, unique-site 0.2%, shared 31% | `17:115-133` | ✅ 0.4238 / 0.0019 / 0.3112 |
| 31 | Same holds on PC2–PC3 | same | ⚠ PC2 0.05% ✅; PC3 0.24% (> the stated ≤0.2%) |
| 32 | Within clear cell (2+5), site explains 4–6% | `17:143-151` | ⚠ **M3** (wrong axes; within-CC PCA gives 13–27%, all n.s.) |
| 33 | RNA subtype 5.9% vs site 3.5%; protein 8.6 / 0.0 / 0.9 | `rna/prot_variancepartition.csv` | ✅ 5.949/3.525; 8.555/0.000/0.852 |
| 34 | Line/patient identity the largest structured term | same | ⚠ **M8** (28 levels / 31 obs) |
| 35 | Silhouette HGS 0.16, CC 0.12, MC 0.15, MMMT 0.74, SCCOHT 0.82; EC do not co-cluster | `rna_silhouette.csv` | ✅ 0.162/0.121/0.154/0.738/0.818; EC −0.014 |
| 36 | 16 of 22 markers land (top-2-of-6) | `rna_markers_summary.csv` | ✅ 16/22 |
| 37 | Median oriented AUC 0.69, median \|d\| 0.72, 8/22 AUC ≥ 0.80 | `rna_marker_effectsizes.csv` | ✅ 0.69 / 0.72 / 8 |
| 38 | KRT20 d 3.1, CDX2 2.8, SPP1 1.5, HNF1B 1.4 | same | ✅ 3.082/2.774/1.461/1.364 |
| 39 | PAX8 d 0.53, MUC2 −0.29, MECOM ≈0 | same | ✅ 0.529/−0.294/−0.064 |
| 40 | Effect sizes over 28 patient reps, "script 04" | `17:300-308` | ⚠ n=28 ✅; produced by `17`, not `04` |
| 41 | GO: HGS DNA-repair / CC oxidative-glutathione / SCCOHT cell-cycle recovered; MC suggestive; MMMT not | `rna_de_gsea_recovery.csv` | ✅ all five statuses exact |
| 42 | Strongest HGS GO axis is OXPHOS/ribosome biogenesis | `rna_de_gsea_go.csv` | ✅ top 8 are RNP/ribosome/mito-translation/OXPHOS |
| 43 | 30 lines, 8,212 shared genes; per-gene median 0.40 | `12`; `integ_rnaprot_cor*.csv` | ✅ 30; recomputed 8,212; 0.3971 |
| 44 | Per-line median 0.41, IQR 0.36–0.44 | same | ✅ 0.4076, 0.3636–0.4390 |
| 45 | External benchmarks 0.38/0.45/0.42/0.48/0.58, ceiling 0.72 | `12:66-74` literature table | ⚠ literature constants, not computed |
| 46 | Protein spread ≈0.30× RNA | `prot_dynamic_range.csv` | ✅ 0.299 (n-matched, 30 lines) |
| 47 | Protein-IQR terciles: 0.30 / 0.40 / 0.53 | `19:163-167` | ✅ 0.302/0.397/0.526 |
| 48 | ConsensusOV DIF 6, MES 4, IMR 3, PRO 2; 7/15 MES-or-IMR | `18:253-263` | ⚠ **M5** counts ✅ but parsed from notes; attributed to `10` |
| 49 | TP53 in 17/17 HGSC lines = 11/11 patients | `wes_driver_freq_patient.csv` | ✅ |
| 50 | CDK12 6/17 lines (35%) → 2–3/11 patients (18–27%) | same | ✅ 35.3% → 27.3% / Tier1-2 18.2% |
| 51 | Defensible somatic BRCA1/2 = 0; both BRCA2 Tier 3 (quiet genome; VAF 0.038) | `wes_driver_tiers.csv` | ✅ TOV81D p.R2845X; TOV3133D p.A22E VAF 0.038 — **but see C3** |
| 52 | KRAS G12/G13 hotspots in CC and MC; SMARCA4 truncation in TOV112D | same | ✅ TOV21G G13C, TOV3392D G12S, TOV2414 G12A; SMARCA4 p.L639X Tier1 |
| 53 | Arm freq per patient: 3q 82, 20q 91, 17p 82, 8q 73, 13q 64, 19q 55 | `wes_cnv_arm_freq_patient.csv` | ✅ all six exact |
| 54 | FGA HGSC 0.62 (11) > CC 0.37 (2) > MC 0.32 (1) > EC 0.23 (1) > LGS 0.02 (1) | `08:211-215` | ✅ 0.622/0.371/0.321/0.226/0.021 |
| 55 | chrX exclusion barely changes **high-FGA HGSC** genomes (0.61→0.62) | `08:216-218` | ❌ 0.61→0.62 is the **all-23-line** median; HGS-only is 0.627→0.635 |
| 56 | TOV81D 0.073 → 0.021 | `wes_cnv_fga.csv` | ✅ |
| 57 | 5 DepMap lines; Spearman 0.74–0.88; rank 1 of 67; reciprocal-best | `external_ccle_concordance.csv` | ✅ 0.739–0.879; all rank 1; reciprocal verified both ways; overlap exactly 5 |
| 58 | DepMap recovers TOV21G hypermutation (568 vs 7–18) | `output/external/depmap_damaging_overlap5.csv` | ✅ 568 / 7,17,17,18 |
| 59 | DepMap: TOV112D R175H + SMARCA4 damaging, OV90 SMAD4, all 5 TOV21G CC drivers | `18:226-237` | ⚠ printed to console only; no output file |
| 60 | TOV112D KRAS-A59T and OV90 BRAF in-frame indel not corroborated | `18:241-244`; `wes_driver_tiers.csv` | ⚠ A59T Tier2 ✅ in tiers; DepMap side console-only |
| 61 | COV434 & BIN67 SMARCA4/SMARCA2 loss; Cellosaurus flags COV434 | `auth_swisnf_panel.csv`, `cellosaurus_str_status.csv` | ✅ both deficient; flag text "Misclassified … SCCOHT" |
| 62 | BIN67 SMARCA4 mRNA retained, protein second-lowest | `auth_swisnf_panel.csv` | ✅ RNA rank 23/31; protein rank 2 |
| 63 | TOV112D p.L639X, low protein, silenced SMARCA2 mRNA, TP53 R175H | `auth_swisnf_panel.csv`, `wes_mutations_filtered.csv` | ✅ protein rank 4; SMARCA2 RNA rank 4; R175H |
| 64 | STR documented for 30 of 42; 12 without = 11 VOA + TOV3121D | `cellosaurus_str_status.csv` | ✅ exactly; all 30 exact-name matches |
| 65 | 26 expression-consistent, 20 genomics-consistent | `auth_perline_table.csv` | ⚠ counts ✅ as written, but **C2** → 20 becomes 19; and 26 rests largely on absence |
| 66 | OV90 PAX8/WT1/SOX17/CK7 ≈0, lowest serous profile among HGS | `10:296-298`; `rna_tpm.csv` | ✅ log2 TPM 1.14/0.48/0.00/1.87; serous z −1.25 = lowest |
| 67 | VOA5436 strong clear-cell program | `auth_perline_table.csv` | ✅ CC_prog_z 2.04 |
| 68 | VOA4841 lowest SMARCA4 of all lines | `auth_swisnf_panel.csv` | ✅ RNA rank 1, protein rank 1 |
| 69 | VOA4395 becomes sole EC | `10:450-451` | ✅ (hard-coded flag; supported by counts) |
| 70 | TOV2414 KRT7+/PAX8+/MUC5AC+/**SATB2−**, KRAS G12A | `auth_mucinous.csv` | ⚠ KRT7 z 2.58, PAX8 0.28, MUC5AC 4.25, KRAS G12A ✅; **SATB2 z = +0.17** (literature, not data) |
| 71 | VOA8762/8771 CK7-low, PAX8-low, CDX2-high, no provenance, no STR | `auth_mucinous.csv`, `cellosaurus_str_status.csv` | ✅ |
| 72 | Genomic HRD not computable (no allele-specific CN) | `09`; `wes_hrd_feasibility.md` | ✅ no BAF columns, no `--vcf`, no tumour BAMs; scarHRD/copynumber absent |
| 73 | Table S1: one row per line incl. **tiered** driver calls | `20`; `supplement_per_line.csv` | ❌ **C3** (42 × 41 ✅, but no tier column; Tier-3 calls unflagged) |
| 74 | FOLR1 bimodal within HGSC, RNA range 0.06–9.5; top TOV3133G | `13:181-186`; `adc_expression.csv` | ✅ 0.057–9.479; top 3 all HGS, TOV3133G first |
| 75 | MSLN protein range 1.8 vs RNA 9.3; protein narrowed 3–5× | `prot_dynamic_range.csv` | ✅ 1.756 / 9.299; ratios 0.19–0.31 |
| 76 | MSLN highest in HGSC; HER2 higher in CC and MC | `adc_subtype_summary.csv` | ✅ both assays |
| 77 | HER2→TOV3392D, NaPi2b→VOA12539, CDH6→OV3331 cross-assay consistent | `adc_expression.csv` | ✅ each is top line in **both** RNA and protein |
| 78 | TOV21G 1,416 coding, 6.9× median, 3.4× next; indel-enriched | `wes_mutation_load.csv` | ✅ 1,416; 6.86×; 3.43×; indel frac 0.279 |
| 79 | SBS6 cosine 0.88 + SBS44/SBS15/SBS20; no POLE signal | `wes_msi_mmr.csv` | ✅ 0.877, top3 SBS6;SBS44;SBS15; cos_pole 0.315 — but **M7** |
| 80 | No MMR-enzyme coding mutation present | `wes_msi_mmr.csv` | ⚠ **M6** (`EPCAM(truncating)` recorded, omitted from text) |
| 81 | TOV3392D not hypermutated → line-specific | `wes_mutation_load.csv` | ✅ 196 coding, `no hypermutation signal` |
| 82 | Biallelic ARID1A truncation, PIK3CA H1047Y, KRAS G13C in TOV21G | `wes_driver_tiers.csv` | ⚠ ARID1A truncating + PIK3CA + KRAS G13C ✅; "**biallelic**" is not computable from total-CN data |
| 83 | 15 HGSC lines → 3 strata; PROGENy corroborates; families mostly co-cluster | `hgs_heterogeneity.csv` | ⚠ 15/3 ✅; PROGENy not orthogonal + no statistic; 2 of 3 families co-cluster (1369 splits) |
| 84 | 2 SCCOHT (COV434, BIN67) and 2 carcinosarcoma (VOA5217, VOA5436) models | `line_family_map.csv` | ✅ |
| 85 | Passage / batch structure (Fig S2) | `rna_passage_check.csv` | ⚠ generated but never cited; asset caption contradicts `partial_r2 = 0.148` |

**Totals: 85 quantitative claims traced. 57 verified ✅ · 24 qualified ⚠ · 4 mismatched ❌ (#6, #25, #55, #73).** No claim was wholly untraceable, but #45 is literature, #22/#59/#60 exist only in console output, and #48 is inherited metadata.

---

## 8. What I could not check

1. **Upstream provenance of the archived intermediates.** Sarek/Mutect2 parameters, the WES capture kit, the CNVkit pooled-normal build, and the entire Morin-lab TMT search pipeline (PIS ratio normalisation, 1% FDR q-values, "CV replicates") are documented in script comments from source-file inspection but not re-executed. The claims that ch1 is a PIS and that abundances are PIS-ratio-normalised rest on the vendor `Readme.md`; I confirmed only that ch1 and ch9 are absent from the 45 data columns, consistent with the description. Needs the raw MS and the search parameters.
2. **Whether the archived `wes_mutations_filtered.csv` was produced by *this* filter logic.** I confirmed the current logic reproduces its row count exactly (6,036) from raw MAFs, which is strong but not conclusive — the file predates the current script. A per-row diff after uncommenting `07:429` would settle it.
3. **DE / GSEA p-values under patient collapse.** I quantified the pseudoreplication effect on PCA/silhouette (M4) but did not re-run DESeq2 on the 28 representatives. Re-running `03` with `dds[, rep_keep]` would settle how much of the HGS GO recovery survives.
4. **Whether the deposited matrices match `output/`.** Accessions are placeholders; nothing to compare.
5. **Genome-build verification beyond spot checks.** I confirmed the MAF `NCBI_Build` field reads `GRCh37` and trust the script's contig-length/hotspot argument for GRCh38, but did not independently re-liftover coordinates.
6. **Figure content.** I read every figure script and verified the data feeding each panel, but did not visually inspect the rendered PDFs/PNGs; a panel could still be mislabelled in a way the code does not reveal. Fig 1's coverage matrix, Fig 4's oncoprint geometry, and Fig 6's dual-ramp heatmap would benefit from a visual pass.
7. **Cellosaurus/DepMap currency.** The cached snapshots (DepMap 24Q4, Cellosaurus JSON of 23 Jul) were taken as given; I verified internal consistency (RRIDs match Model.csv, 5-line overlap, 67 expression models) but not that they are current.
8. **`reports/` narrative documents** — out of scope by instruction, so I cannot say whether the HTML/report text repeats or corrects the issues above. Note that `reports/assets/f_auth_swisnf.png` is referenced by `build_report.py` with no generating script, so at least one report figure is untraceable.
