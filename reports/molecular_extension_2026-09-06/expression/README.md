# Patient-aware expression and target-selection extension

6 September 2026. This is a separate exploratory analysis of the current matrices. It adds patient-aware concordance, explicit raw-unit/availability records, joint-ranking sensitivity and individual model leads to the existing script 13 ADC atlas. It does not repeat the atlas's subtype means or FOLR1 distribution modelling. No clinical receptor, MMR, SWI/SNF, amplification or treatment-response categories are assigned.

## Main findings

Seventeen of 19 targets have protein measurements in all 31 proteomics models. DPEP3 and PGR have no quantified protein feature. Each gene has RNA in 31 models; the modalities overlap in 30 models from 27 patients. All 42 resource models remain in the expression table and atlas, including models without either assay.

FOLR1, MSLN and TACSTD2 retain strong positive RNA/protein concordance after selecting one model per patient. ERBB2, SLC34A2, CDH6, CD276 and VTCN1 also show positive concordance, with differing strength across histotypes. ESR1 is discordant, and the smaller correlations for SMARCA4, ARID1A and CCNE1 do not support using RNA alone as a reliable substitute for their supplied protein measurements.

| Target | Patient representatives, n=27 | HGSC representatives, n=12 | Patient means of same-model pairs, n=27 |
| --- | ---: | ---: | ---: |
| ERBB2 | 0.524 | 0.692 | 0.564 |
| FOLR1 | 0.924 | 0.874 | 0.894 |
| MSLN | 0.839 | 0.816 | 0.830 |
| TACSTD2 | 0.882 | 0.979 | 0.889 |
| SLC34A2 | 0.737 | 0.797 | 0.747 |
| CDH6 | 0.686 | 0.601 | 0.722 |
| CD276 | 0.623 | 0.420 | 0.643 |
| VTCN1 | 0.778 | 0.881 | 0.786 |
| ESR1 | -0.257 | -0.294 | -0.254 |
| MLH1 | 0.555 | 0.427 | 0.548 |
| PMS2 | 0.546 | 0.545 | 0.531 |
| MSH2 | 0.382 | 0.049 | 0.420 |
| MSH6 | 0.274 | 0.140 | 0.156 |
| SMARCA4 | 0.141 | 0.056 | 0.158 |
| SMARCA2 | 0.487 | 0.517 | 0.513 |
| ARID1A | 0.226 | 0.336 | 0.203 |
| CCNE1 | 0.132 | 0.119 | 0.090 |

Values are Spearman correlations, not prediction accuracy. DPEP3 and PGR have zero complete pairs and no correlation. `concordance.csv` contains all values, denominators, model-level comparisons, six-patient clear-cell comparisons and patient-bootstrap intervals. For example, the FOLR1 interval is 0.824-0.961, whereas ESR1 spans -0.603 to 0.135. These resampling intervals describe sensitivity to the patients represented in this convenience panel; they do not establish population-wide validation.

## Model leads and important discordance

The table below identifies focused follow-up leads from the 27 preselected patient representatives. Joint rank combines within-cohort RNA and protein ranks, with methods detailed below. It is relative to this panel and does not imply an expression threshold. Raw protein units are supplied log2 normalised abundance and must not be compared across different proteins as absolute concentrations.

| Target | Leading model | RNA TPM | Protein | Follow-up interpretation |
| --- | --- | ---: | ---: | --- |
| ERBB2 | TOV3392D | 216.377 | 13.904 | Highest joint rank under all four assay-weighting methods and patient averaging; isoDoping flagged |
| FOLR1 | OV2295 | 277.290 | 13.235 | Stable with assay weights, but patient 2295 falls to rank 8 after averaging its measured sublines |
| MSLN | OV4485 | 385.762 | 13.081 | Joint rank 1-2 across assay-weighting methods |
| TACSTD2 | OV4485 | 930.943 | 14.974 | Joint rank 1-2; OV4453 is a close alternative with protein 15.184 and RNA 807.272 |
| SLC34A2 | VOA12539 | 425.436 | 13.575 | First under all four weighting methods and patient averaging |
| CDH6 | OV3331 | 580.406 | 14.717 | First under all four weighting methods and patient averaging |
| CD276 | VOA8762 | 153.795 | 11.947 | First under all four weighting methods and patient averaging |
| VTCN1 | OV3331 | 99.890 | 12.726 | First under all four weighting methods and patient averaging; isoDoping flagged |
| ESR1 | OV2295 | 18.696 | 12.662 | Joint lead, but the across-patient RNA/protein correlation is negative and protein ranks do not establish endocrine function |

The FOLR1 subline effect is substantive: OV2295-R2 has 30.955 TPM and protein 11.577, compared with OV2295's 277.290 TPM and 13.235. It is appropriate to choose OV2295 as a particular high-expression model, but not to describe patient 2295's models as uniformly high. VOA6861 has the highest measured FOLR1 RNA, 712.639 TPM, but lacks proteomics and therefore receives no joint score.

For the MMR and SWI/SNF panel, the ranking direction is deliberately **lower** abundance, to prioritise expression discrepancies for follow-up. TOV21G is lowest jointly for MLH1 (1.685 TPM, protein 12.164). TOV3291G is lowest jointly for MSH2 (1.652 TPM, 13.868) and MSH6 (5.803 TPM, 14.128). TOV2414 is lowest jointly for PMS2 (5.520 TPM, 10.287). These findings identify measurements worth checking; they do not diagnose mismatch-repair deficiency or MSI. TOV21G illustrates the distinction: its PMS2 protein is second-lowest, while PMS2 RNA is 15.988 TPM rather than at the bottom of the panel.

VOA4841 is lowest jointly for SMARCA4 (4.995 TPM, 13.306), and OV2085 for SMARCA2 (0.121 TPM, 10.961). BIN67 has low SMARCA4 protein (13.315; second-lowest among 31 models) but RNA 94.707 TPM. Its joint low-expression rank shifts from 5 to 14 as the weighting/weakest-assay rule changes. Thus a combined RNA/protein score should not conceal single-modality evidence, and relative protein abundance does not establish protein absence. The ARID1A joint low-expression lead is VOA5436; ARID1A is isoDoping flagged and its patient-level RNA/protein correlation is only 0.226. No functional loss or variant mechanism is inferred.

DPEP3 and PGR remain RNA-only leads: DPEP3 is highest in OV1369-R2 at 23.133 TPM, followed by the related TOV1369 at 2.176 TPM. PGR is highest in OV2295 at 4.816 TPM. DPEP3 has 25 measured zero-TPM models and PGR 21; the remaining missing assay rows are distinct from those zeros. Neither gene has a quantified protein row, so protein-negative or receptor-negative language is unsupported.

The CCNE1 comparison is intentionally retained as a reference. OV4453 and VOA12539 tie for the leading joint score, while TOV3291G's rank changes from 3 to 8 across weighting methods. This does not replace the earlier RNA/CNV review: different modalities support different model-selection questions. OV3291 and TOV3291G remain separate model rows despite their shared patient origin.

## Methods and denominators

**Inputs and mappings.** Script 42 reads `metadata/resource_models.csv`, the corrected release-93 `rna_tpm.csv` and `rna_gene_annotation.csv`, `prot_abundance_matrix.csv`, and `prot_qc.csv`. It sums TPM only where multiple retained primary-assembly Ensembl genes share a symbol, before calculating log2(TPM+1). Protein mapping uses the existing unique symbol-representative row. No live annotation, imputation, rescaling or new normalisation is introduced. The raw table includes identifiers, passages, TMT plex/channel, gene-level peptide counts and isoDoping flags. Protein passage is not inferred from RNA passage.

**Availability.** A missing assay, an unquantified protein target and a measured RNA value of zero have separate status fields. The 42-model table is complete for every target even though 11 models lack RNA and 11 lack protein. Protein feature annotations such as isoDoping and global peptide count can be present on rows with no protein assay; the abundance remains missing. Global peptide counts are not sample-specific peptide support.

**Concordance.** Complete pairs are assembled by exact model name before any patient operation. Primary analysis retains the globally preselected model per patient, without choosing the highest target expresser or substituting another subline if the selected model lacks a pair. The second patient analysis averages log2(TPM+1) and protein log2 abundance across each patient's available same-model pairs. It does not match an RNA-only model to a different protein-only model. The accompanying `rna_tpm` field in the mean cohort is the arithmetic mean TPM; the correlation/ranking uses the separately labelled mean log2(TPM+1), not log2(mean TPM+1).

Correlations are unadjusted. HGSC-only and clear-cell-only results address some histotype composition sensitivity, but centre, TMT plex, culture and passage remain uncontrolled. No p-values or discovery claims are produced. Bootstrap intervals resample selected patient rows 2,000 times with a fixed seed, using percentile limits; intervals are computed only with at least six paired patients and a nonconstant correlation. Valid resample counts are saved. These intervals do not represent assay precision.

**Joint ranking.** Within each paired cohort, each assay is converted to an average-tie rank percentile `(rank - 1)/(n - 1)`. The preferred direction is higher for surface/endocrine/CCNE1 targets and lower for MMR/SWI/SNF targets. Larger preference percentiles mean stronger support in that declared direction. Four scores are compared: equal-weight mean, minimum of the two percentiles (the weaker assay), 25% RNA/75% protein, and 75% RNA/25% protein. Descending joint ranks use minimum rank for ties. No score is calculated if an assay is missing or a scale is constant. The `both_assays_same_direction` field only indicates whether both preference percentiles lie in their respective cohort's upper half; it is not a biological or clinical cutoff.

The top-three-rank lead table retains ties, so a target can have more than three entries. It also reports the same patient's rank under within-patient averaging and within HGSC or clear-cell representatives. These are sensitivity comparisons with different denominators. Percentile ranks are convenient for combining scales but discard distances: a very small protein difference can span many ranks. Raw units and the assay-specific percentiles must therefore accompany any selected lead. The figure selects one leading row alphabetically when equal scores tie; the CSV retains all ties.

**Relationship to isoDoping.** ERBB2, VTCN1, ESR1, ARID1A and CCNE1 are flagged in the supplied QC. The flag does not automatically invalidate a measurement, and observed discordance cannot be attributed to it without the upstream processing evidence. See the [independent CCNE1 assay review](../../clinical_classification_2026-09-06/CCNE1_EXPRESSION_REVIEW.md) for source-workbook and published-method checks. No gene is discarded or given a quantitative correction merely because it is flagged.

## Files and reproduction

Run from the repository root, with Python containing numpy/pandas and the existing R figure dependencies:

```bash
python3 scripts/42_expression_target_extension.py
Rscript scripts/43_expression_extension_figures.R
```

`OVCAN_PROJ` can select the project root. The analysis uses processed matrices and does not need cluster alignments or raw proteomics. Script 43 sources only the shared figure theme, avoiding the canonical setup's session-log rewrites.

| File | Contents |
| --- | --- |
| `panel.csv` | The 19 genes, groups, rank direction and interpretation limits |
| `model_expression.csv` | 798 rows: 42 models x 19 genes; raw units, status, ranks and provenance fields |
| `target_summary.csv` | Measurement counts, zero/positive RNA counts, raw-unit ranges and peptide/isoDoping annotations |
| `unpaired_expression.csv` | Rows with exactly one measured modality; no joint score |
| `cohort_subjects.csv` | Exact rows/values underlying each correlation cohort |
| `concordance.csv` | Five cohort comparisons per target, correlation and bootstrap denominators |
| `joint_ranking.csv` | Complete pairs, all four ranking rules and cohort-specific denominators |
| `model_leads.csv` | Top-three-rank representatives and cross-method/cross-cohort sensitivity |
| `discordant_models.csv` | Largest RNA/protein preference-percentile gaps per gene among representatives |
| `analysis_qa.json` | Input hashes, versions, seeds, invariants and existing-atlas agreement |
| `visual_qa.json` | PDF export hashes, dimensions, fonts and render-review results |

## Figure legends

**Expression concordance and rank sensitivity** (`output/pdf/expression_concordance_rank_sensitivity.pdf`). **A:** Spearman correlation of same-model RNA log2(TPM+1) with supplied protein log2 normalised abundance for 27 selected patient models (rust circles) and the 12 HGSC representatives (slate diamonds). Pale horizontal intervals are the 2.5th and 97.5th percentiles from 2,000 patient-resampling replicates for the 27-patient cohort. They describe resampling sensitivity within this convenience panel. DPEP3 and PGR lack protein and have no correlation. **B:** For each target, the labelled representative is the leading equal-weight joint-rank model; ties are resolved alphabetically for display. Rust circles show its primary joint rank. Horizontal bars span ranks under equal weighting, the weaker-assay score, and 25% or 75% RNA weighting. Open diamonds show the same patient's rank after averaging available same-model pairs within each patient. All ranks use 27 patients. Higher expression is prioritised for surface/endocrine/CCNE1 targets; lower expression is prioritised for MMR/SWI/SNF targets. The CCNE1 equal-score lead VOA12539 is also retained in the data table. Asterisks denote a supplied isoDoping flag. Group separators divide surface, endocrine, MMR, SWI/SNF and CCNE1 targets. Correlations and ranks are exploratory and do not define clinical positivity, protein loss, localisation or response.

**Expression target model atlas** (`output/pdf/expression_target_model_atlas.pdf`). Rows show all 42 models, grouped by recorded histotype and alphabetically within histotype. Columns show the same 19 targets in **A**, RNA, and **B**, protein. Colour represents the within-target average-rank percentile among the 31 measured models for each assay, calculated as 100 x (rank - 1)/(n - 1), with higher abundance toward rust for every gene, including MMR/SWI/SNF. It is not the direction-dependent joint score from the other figure. Crosses mark unavailable values: an absent assay or an unquantified protein target, distinguished in the companion table. DPEP3 and PGR have no quantified protein. Measured zero TPM remains a ranked observation; tied zeros share an average rank rather than all receiving percentile zero. VOA6861 has RNA only and VOA14993 protein only. Raw TPM and supplied protein log2 normalised abundance are in the accompanying table. Asterisks mark genes whose protein feature has an isoDoping annotation. Colour does not compare absolute abundance between proteins and does not establish a receptor threshold or functional protein loss.

## Existing records worth retrieving next

These analyses run locally now; a new whole-alignment or cluster expression analysis is not required. If available in the original proteomics project, retrieve the exact search/quantification scripts and settings, PSM-level reporter-ion quantities with interference/isotope corrections, the isoDoping peptide list and channel design, and the protein rollup/scaling definition. Those records would permit targeted follow-up of ESR1/CCNE1/ARID1A discordance and distinguishing true low reporter signal from compression or rollup effects. Match assay aliquots and passages from laboratory records; do not invent a proteomics passage from the RNA label. Local supplied files contain only abundance/peptide workbooks, layout and an empty README, so these upstream records are not yet documented here.

Functional follow-up is a separate experiment or stock-record request: surface-specific measurement for putative surface-target models, appropriate protein localisation/loss assays for MMR or SWI/SNF candidates, and genotype interpretation with suitable allelic evidence. None is replaced by a percentile score. The parent molecular-extension report coordinates the WES/cluster-specific to-do list.

## Validation

All 42 models and 34 patients remain represented; all 19 target tables have 42 unique rows. Exact model pairing precedes patient averaging, and no missing values are imputed. All 527 finite comparisons against the existing ADC atlas agree within its three-decimal rounding (maximum absolute difference 0.00049919). The source matrices and metadata are recorded by SHA-256 in `analysis_qa.json`. Both PDFs were rendered with Poppler and inspected, including all model labels, missing-data marks, axis labels, legends and group boundaries. The figure script uses the shared rust/slate theme, Arial and Quartz PDF output. No canonical manuscript, release, prior CCNE1 output or scientific matrix was changed.
