# RNA, proteomics, integration, and confounder audit

Audit date: 5 September 2026. Scope: scripts 01–06, 12–14, 17, 19, and 21; their input data and current outputs; supplementary figure builder 37; the v5 manuscript and its Word comments. Raw observations and the v5 manuscript were preserved. The reference repair required coordinated rebuilding of downstream authentication, external-comparison, and figure outputs by the other agents.

## Main finding: the transcript-map problem is repaired

The previous analysis inferred an “Ensembl-104-era” reference from the October 2021 quantification date. That inference was incorrect: the date does not identify the reference release. The release-105 map omitted 3,529 of the 185,299 indexed transcripts, discarding 1.60–3.41% of library TPM (median 2.22%). This was not only a presentation caveat: individual genes could be incompletely quantified or absent.

A focused search of adjacent local analysis references recovered the Ensembl release-93 GRCh38 cDNA FASTA. Independent parsing established an exact match for **all 185,299 versioned transcript IDs and all corresponding sequence lengths in each of the 31 libraries**. The recovered FASTA also matches Ensembl's official release-93 checksum record: BSD sum **19566**, **65198** blocks. The original binary kallisto index is not available, so an original index sequence checksum cannot be compared; the evidence identifies the matched transcript reference through its complete ID/length fingerprint and the verified public FASTA.

The FASTA and a freshly derived map are now pinned in `data/reference/`. The deterministic builder is `scripts/build_matched_rna_reference.py`; its SHA-256 assertions reproduce the map byte for byte. Provenance and checksums are recorded in `data/reference/rna_reference_provenance.json`. The public source is the [Ensembl release-93 cDNA directory](https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/), with its [checksum record](https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/CHECKSUMS).

Script 01 now uses **exact versioned mapping**, verifies the complete target set and lengths before import, and asserts conservation of both estimated counts and TPM after gene summarization. All 31 libraries passed. Re-summarizing existing transcript estimates repairs this mapping error; raw-read re-quantification is not necessary for the assignment correction. No live annotation lookup remains in this import path.

| Measure | Before | Corrected |
|---|---:|---:|
| Transcripts without a gene map | 3,529 | 0 |
| Gene-level TPM rows | 39,568 | 39,733 |
| Genes retained by ≥10 estimated counts in ≥2 models | 22,544 | 22,678 |
| Median genes detected, ≥1 count among retained genes | 20,119 | 20,270 |
| Median pseudoalignment | 91.1% | 91.1% |
| Median processed fragments | 64.3 million | 64.3 million |

The before/after gene sets share 39,433 Ensembl IDs; the corrected map introduces 300 IDs and removes 135 IDs relative to the later mismatched annotation. Reference-matched symbols also differ because the two annotation releases use different symbol assignments. Stable Ensembl IDs, reference release, and locus annotation therefore accompany the exported data.

The principal expression structure is stable: across shared symbols, median sample-profile Spearman correlation before versus after repair is **0.99822** (minimum 0.99795), and absolute correlation between old and corrected PC1 scores is **0.999965**. Individual corrections can nevertheless be important. *FOXL2* was absent from the former symbol matrix and is now quantified in every model; its corrected TPM is 0.0219 in COV434 and 0.0429 in BIN67. High-abundance restored transcript contributions include *RPS27*, *GSTP1*, *RPL26*, and *RPS10*. The absence of a gene from the old table must not be interpreted as absent expression.

All original reference loci are retained, including 17,608 alternative-locus/patch transcripts. The map records `seq_region` and `primary_assembly`; Ensembl gene rows remain separate, while symbol queries sum TPM over identically named gene rows, as throughout the existing analysis. This avoids silently losing quantification assigned to original index targets. A primary-assembly-only PCA sensitivity, reselecting the 2,000 most variable primary-assembly features from the existing VST, gives absolute PC1–3 score correlations of **0.99970, 0.99903, and 0.99781** against the all-locus analysis. This is a feature-selection sensitivity, not a new raw-read quantification.

Evidence: `scripts/01_rna_load_qc.R:36`, `:62`, and `:88`; `output/rna_reference_reconciliation.csv`; `output/rna_gene_annotation.csv`; `reports/audit_2026-09-05/check_rna_reference.R`; `rna_reference_before_after.csv`, `rna_reference_profile_sensitivity.csv`, `rna_restored_gene_summary.csv`, `rna_marker_reference_delta.csv`, and `rna_primary_assembly_pca_sensitivity.csv` in this audit directory. Prior outputs and their SHA-256 manifest are under `rna_before/` and `rna_before_sha256.json`. The obsolete working release-105 map was moved out of `output/`; its pinned historical reference remains available for the comparison.

## Corrected methodological and interpretation issues

### High: dependent sublines were still receiving independent-observation significance language

Script 03 called HGS differential expression “formal” solely because its group size exceeded ten, although the 15 HGS models represent 12 patients and all share one contributing centre. Sample size alone does not license that interpretation. The model-level DE tables and figures now explicitly describe **panel associations**, and script 21's patient-representative tables are explicitly exploratory associations with the source-confounding limitation retained.

The previous comment that adding centre was unidentifiable was also too strong. The additive design is full-rank because some non-HGS histotypes span centres. Nevertheless, absent histotype–centre combinations prevent empirical checking of a common additive centre effect. A full-rank model cannot distinguish a source-aligned biological effect from an unmeasured culture/processing effect. The [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) describes these design-rank and confounding issues; here, the source limitation is scientific identification rather than a universal matrix-rank failure.

Script 12 likewise calculated gene-wise correlation q-values across 30 models and described significant inverse correlations as established associations without carrying donor dependence into those claims. A new `integ_rnaprot_patientrep_cor.csv` contains gene-wise estimates and asymptotic Spearman tests on **27 patient representatives**, with BH correction across the eligible patient-level genes. The combined all-model table preserves the descriptive correlations and labels its inferential fields exploratory. Exact ±1 correlations no longer receive missing p-values from the previous division-by-zero guard.

| Concordance summary | Corrected value |
|---|---:|
| Shared RNA/protein gene symbols | 8,358 |
| Genes with ≥10 paired models and nonzero variance | 8,033 |
| Median per-gene Spearman, model-level description | 0.3998 |
| Complete-case genes across all 30 paired models | 6,811 |
| Complete-case median per-gene Spearman | 0.4105 |
| Median sample-wise Spearman across genes | 0.4130 |
| Eligible genes on 27 patient representatives | 7,969 |
| Median patient-representative per-gene Spearman | 0.4182 |
| Inverse associations with BH q<0.05, patient representatives | 21/7,969 |

These are cross-model associations, not demonstrations of post-transcriptional mechanisms, predictive performance, or accuracy of individual protein measurements. Histotype, source, and plex structure have not been removed from the gene-wise correlations. See `scripts/12_rna_protein_concordance.R` and `output/integ_rnaprot_{cor_summary,patientrep_cor,patientrep_sensitivity,negative_genes}.csv`.

### Medium: the manuscript compared adjusted and unadjusted commonality statistics

V5 juxtaposed an adjusted unique-histotype component of approximately 0.393 in 31 models with an unadjusted component of approximately 0.523 after selecting patient representatives. It also described a negative adjusted component as though it were a variance share. Script 21 now writes both definitions explicitly. The main figure revision uses **raw components consistently**.

| Corrected PC1 statistic | 31 models | 28 patient representatives, HVGs refitted |
|---|---:|---:|
| Variance represented by PC1 | 20.632% | 19.858% |
| Raw unique histotype R² | 0.423852 | 0.520132 |
| Raw unique centre R² | 0.001913 | 0.002689 |
| Raw shared R² | 0.309954 | 0.269474 |
| Adjusted unique histotype increment | 0.392571 | 0.505662 |

Adjusted increments can be negative because of the complexity penalty; they are not physical negative variance components. Low unique-centre R² means little additional prediction conditional on the recorded histotype labels. It does not establish that centre has no effect.

The existing 1,000 all-model label permutations were also not patient-aware. Script 21 now adds **9,999 reduced-model residual permutations** on the 28 representatives: residuals are shuffled within levels of the reduced categorical factor and added back to fitted values, keeping the full design fixed. The patient-level PC1 histotype test gives F=10.0167, df=5,20, p=0.0011; BH correction across PC1–3 gives q=0.0033. PC1 centre p=0.5146. These tests condition on within-group exchangeability and remain limited by sparse and empty design combinations. They are evidence of histotype-associated expression structure, not removal of all source confounding.

Evidence: `scripts/21_rna_sensitivity.R:124` and `:154`; `output/sensitivity_patient_reps_pca.csv`; `output/sensitivity_patient_reps_pc_permutation.csv`. Script 17's all-model permutation file now explicitly redirects primary inference to the independent-patient result.

### Medium: marker validation needed source conditioning and consistent analysis units

Joint marker-label permutations already preserved correlation among marker genes. They did not preserve the source composition of histotypes, however. A new patient-representative permutation restricted **within contributing centre** conditions on that design structure. The corrected map retains **15/25** recovered markers. Against 20,000 source-restricted draws, the null mean is **8.4602**, its 95th percentile is **13**, and p=**0.01005**. Unrestricted patient permutations give null mean 6.8223 and p=0.00480. The restricted result is the stronger basis for the manuscript's source-aware association claim.

The old independent binomial p field is retired (NA); neither correlated markers nor unequal histotype sample sizes justify independent Bernoulli(1/3) trials. Script 17 also previously joined rankings computed on 31 models to effect sizes computed on 28 patient representatives. The rankings and recovery flags are now recomputed on those same 28 representatives. The small-group bootstrap caveat was corrected: an interval may be poorly calibrated at n=2 even if every resample produces a finite number.

Evidence: `scripts/04_rna_markers_genesets.R:190` and `:233`; `output/rna_marker_recovery_permutation.csv`; `scripts/17_variance_confounders.R` marker-ranking block; `output/rna_marker_effectsizes.csv`.

### Medium: one variance sensitivity invented correlations among unrelated patients

The former “family_multi_line_only” random-effect model pooled every unrelated singleton patient into a single random-intercept group. That is not a replacement for unavailable replication; it gives unrelated models a shared covariance. This model was removed. The patient-dropped and passage-adjusted sensitivities remain.

The full per-feature variance decomposition remains weakly identified: 31 models represent 28 patients, with only three replicated patient pairs. The estimate separating patient variance from residual variance is therefore supported by just those pairs. Its results remain supplementary and descriptive; the main figure was revised to show the better-supported PC commonality comparison. The corrected RNA model fitted 22,674 genes and explicitly lists four fit failures. Do not interpret the largest random component as a biological discovery or the residual as pure technical noise.

### Medium: cross-assay spread and external correlation values were given causal meanings

Script 12 formerly plotted a universal-looking “reproducibility ceiling” at 0.72 and stopped execution if the observed gene-wise median fell outside an arbitrary expected range. Both were removed. The verified CPTAC ovarian values, mean 0.38 and median 0.45, may provide context but cannot calibrate accuracy across different assays, feature sets, and populations. These values are reported by [Zhang et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4967013/). The interpretation in [Upadhya and Ryan, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9499981/) concerns measurement reproducibility and attenuation; it does not supply a dataset-specific ceiling for this experiment.

Script 19's output and exploratory figure language now describes **observed cross-assay spread**. The corrected set contains 8,035 genes with matched measurements in at least ten models; the median protein/RNA IQR ratio is **0.2990** among the 8,022 with nonzero RNA IQR (13 zero-IQR genes have no defined ratio). Two constant-RNA genes enter the spread table but cannot enter the correlation table, explaining 8,035 versus 8,033. This ratio does not isolate TMT ratio compression from regulation, turnover, assay processing, feature detection, or the log2(TPM+1) transformation. IQR is less sensitive than the range to isolated extremes; neither IQR nor SD is universally insensitive to a detection floor or transformation. Legacy console prose that attributed differences between distinct zero-RNA/nonzero-RNA gene sets causally to the floor was removed. Bridge-derived SD is labelled an approximate precision estimate rather than a theoretical noise floor.

The supplied protein values retain large protein-specific baselines. Therefore, the defensible scale label is **supplied log2 protein abundance**. The exact upstream normalization and scaling formula still requires confirmation; do not assert that these are simple zero-centred log2 ratios to the pooled internal standard. This uncertainty does not prevent describing the observed within-protein differences, but it precludes recovering the upstream pipeline from the processed matrix.

### Medium: within-HGS cluster labels were not stable to patient selection

Script 14 now writes a patient-representative clustering sensitivity for k=2,3,4. At k=3, the 15-model partition and a newly standardized/clustered 12-patient partition have an **adjusted Rand index of 0.189** on the shared representatives. Mean silhouettes are 0.142 and 0.180, respectively. The k=2 and k=4 adjusted Rand indices are 0.183 and 0.095. Thus the named strata should remain an illustration of querying pathway scores, not stable molecular subtypes or Technical Validation. The output carries patient IDs, representative flags, and an explicit exploratory interpretation.

The retained FOLR1 Gaussian-mixture exercise is also explicitly exploratory: two fitted Gaussian components do not necessarily imply a bimodal density, and 15 HGS models represent only 12 patients. Neither that exercise nor total RNA/protein abundance establishes surface-antigen density or ADC response.

Evidence: `output/hgs_cluster_patient_sensitivity.csv`; `output/hgs_heterogeneity.csv`; `output/adc_folr1_bimodality.csv`; scripts 13 and 14.

## Ancillary GO analyses and supplementary figures

The final patient-representative GO comparison changes the HGS repair-program interpretation. The strongest positive term matching the prespecified repair-program pattern has BH q=0.0457 across all 31 models, but **q=0.0646 among the 28 patient representatives**. The best matching term also changes. These are related GO terms, not an exact same-term sensitivity test. CC remains q=0.00470, MC q=0.150, MMMT q=0.569, and SCCOHT q=0.0000519 in the patient analysis. The three-seed model-level precision check was stable (HGS adjusted-p range 0.0454–0.0472); numerical precision does not resolve the difference induced by selecting independent patients. The HGS program therefore does not meet the 0.05 adjusted threshold on independent patients. The main manuscript appropriately omits these ancillary recovery grades. See `output/sensitivity_patient_reps_go_comparison.csv` and the full GO tables, which carry the tested term names and results.

Supplementary figures S1–S4, S7, and S8 were rebuilt in script 37. The final panels keep concise scientific captions by default (`OVCAN_FIG_PLAIN=0` only when the caller has not set a preference). Specific corrections include descriptive passage associations without raw p-value claims; source and patient-origin annotations without authentication claims; all 32 missingness patterns, including the 70 features absent from model channels in all plexes; marker recovery/effect sizes based on the same 28 patient representatives; and the HGS k=3 instability prominently reported as ARI=0.19. PROGENy theme agreement is described as agreement between scores computed from the same RNA, not independent validation. Dynamic captions report the actual group medians instead of a hardcoded Hypoxia ordering.

Visual QA inspected all six PNGs and separately rendered the two complex layout panels S3/S4 from PDF. A clipped dendrogram, point/tile row alignment, crowded passage annotation, and old internal revision prose were corrected. The final PDFs embed Arial/Arial Bold; S3 and S7 were checked with `pdffonts`. The final script-37 log contains no warnings. ComplexHeatmap's internal base-PDF measurement uses a standard metric-compatible Helvetica alias for Arial; final Cairo/ragg exports use the actual Arial font. PDF QA renders are retained as `figs3_pdf_qa.png` and `figs4_pdf_qa.png`.

## Existing methods checked and retained

* RNA inference uses estimated counts and tximport transcript-length normalization factors, not TPM supplied to DESeq2. The global expression filter is documented, and the VST is used for unsupervised display.
* The protein matrix preserves structural plex-level missingness. Complete-case features are used for PCA; the broader matrix retains feature-presence and zero-plex flags. Deterministic symbol representatives prevent arbitrary selection of duplicate accessions.
* The earlier bridge precision correction is already implemented and remains appropriate as an approximation: SD(single)=SD(pairwise difference)/√2 under equal independent errors, followed by the lognormal CV conversion. Current bridge-derived CVs are **10.35–13.21%**. These are four different adjacent-plex bridge samples, not four aliquots of a universal reference. High across-protein correlation does not alone establish agreement.
* The ADC atlas reports both modality-specific line sets and within-feature relative expression. It remains a model-querying example, with no evidence here for ADC sensitivity.
* RNA-versus-WES passage differences are documented for the available paired passage records; protein passage matching cannot be assumed from those records.

## Validation and remaining boundaries

The exact map builder reproduced its asserted checksum; every library passed complete versioned-ID/length matching and count/TPM conservation. Independent before/after and primary-assembly sensitivity calculations are retained in this directory. An independent `stats::cor.test` check for 120 eligible patient-representative genes agreed with exported Spearman coefficients and p-values to maximum absolute differences of 1.11×10⁻¹⁶ and 3.89×10⁻¹⁶, respectively (`rna_patient_correlation_verification.csv`). The RNA/protein scripts parse successfully; run logs document the affected analyses and the reused R 4.5.2 environment. The parent agent ran script 17 and coordinated the dependent figure and external-authentication rebuilds.

The completed audit does not turn this selected panel into a population-representative experiment. Source centre remains partly aligned with histotype; rare groups remain small; total protein abundance is not surface expression; upstream proteomics metadata and raw-spectrum provenance require collaborator confirmation. Those are bounded study-design or upstream-data limitations, not reasons to retain the now-repaired transcript-map mismatch.

Scripts 03 and 21 completed their full verified reruns; all five GO recovery grades were stable across the three model-level seeds. Scripts 19 and 37 completed their final clean rebuilds. A final wording-only GO plot refresh used the deposited statistics; no statistical computations were changed. The only warning in the full 03/21 runs notes that the installed S4Vectors binary was built under R 4.5.3 while this runtime is R 4.5.2; both runs completed and the independent numerical checks passed. Earlier interrupted streaming runs are retained in the logs; verified reruns use `source()` to parse each complete script before execution. Later comment/console-only corrections clarify that discordant TOV112D molecular evidence warrants annotation review, without treating it as a demonstrated histological reclassification.
