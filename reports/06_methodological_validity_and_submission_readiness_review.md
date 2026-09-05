# Methodological validity and Scientific Data submission-readiness review

**Project:** Multiomic characterization of ovarian cancer cell models  
**Review date:** 27 July 2026  
**Review scope:** methodological choices, statistical interpretation, data-descriptor fit, and manuscript readiness  
**Primary evidence reviewed:** `reports/01_multiomic_characterization_results.pdf`, its Markdown source, `reports/02_manuscript_outline.md`, the analysis scripts and tabular outputs, the workflow reproducibility reports, and manuscript drafts through `docs/manuscript/v4/OvCAN_data_descriptor_v4.md`.

## Executive assessment

The project supports a useful Scientific Data descriptor, but it is not yet submission-ready. The strongest parts are the documented sample/assay inventory, RNA quality control, internally consistent RNA–protein integration, transparent patient-family mapping, and canonical positive controls such as retention of *TP53* alterations in all 11 high-grade serous patients with exome data. The processed matrices and extensive per-model annotations would be reusable if deposited with sufficiently complete metadata.

The main concerns are not that the scripts failed to produce the reported outputs. They are that some statistical summaries answer a narrower question than the prose implies, and that several interpretations exceed what the study design can identify. The largest examples are:

1. histotype and source centre are partly confounded, so the data support histotype-associated structure after statistical adjustment but cannot establish that the structure is biological rather than a source, culture, or processing effect aligned with histotype;
2. several inferential RNA analyses use cell models as independent observations even though 13 models belong to five donor families; patient-representative analyses should be primary for inference;
3. proteomic bridge variability was converted to a “CV” without the replicate-difference adjustment required for a per-measurement estimate, and the comparison of limits-of-agreement width with a cross-model IQR mixes unlike dispersion measures;
4. tumour-only exomes remain unsuitable for interpreting burden or retained variants as definitively somatic, despite useful filtering and strong positive controls;
5. several validation statements—gene-set recovery, model authentication, ratio compression, and a mismatch-repair signature—require more qualified language.

These issues are repairable in the paper. They do not invalidate the resource. They change the appropriate claim from “the multiomic profiles faithfully establish subtype biology and genetic identity” to “the processed profiles pass multiple internal and external checks, while known design and assay limitations define how they should be reused.”

## Review boundaries and evidence integrity

This review accepts the user’s premise that generated outputs faithfully reflect the scripts. I therefore did not perform a second full workflow rerun, which would also risk touching current outputs. Instead, I:

- inspected all 36 pages of the primary results PDF for presentation defects;
- parsed the results report, outline, manuscript drafts, scripts, output tables, and existing clean-room reproducibility reports;
- confirmed that all 35 R scripts parse, `scripts/run_all.sh` passes `bash -n`, and `renv.lock` records R 4.5.2 and 288 packages;
- compared current outputs with the preserved post-fix snapshot: the differences were the expected updated authentication table, expanded WES bootstrap results, session/package timestamps, minor numerical/path changes, and one added protein-matrix README;
- ran one new, read-only statistical interrogation of the marker-recovery null model, described below.

No existing project file was modified or deleted. The auxiliary permutation script and its transient outputs are under `tmp/`.

## Overall claim-validity matrix

| Proposed claim | Assessment | Defensible interpretation |
|---|---|---|
| The resource contains complementary RNA, protein, mutation, and copy-number profiles across a diverse ovarian cancer model panel. | **Supported** | Coverage is uneven but clearly documented: 42 models/34 patients overall; RNA and protein each cover 31 models, exome mutations 22, copy number 23, and 13 models have all three principal layers. |
| RNA data are technically high quality. | **Supported with one reconciliation caveat** | Median pseudoalignment is 91.1%, median depth is 64.3 million paired fragments, and no sample is an obvious multimetric outlier. The pre-release-105 kallisto index and release-105 transcript map mismatch drops 1.90% of index targets and a median 2.22% of TPM; both index and mapping provenance must be deposited, and a matched-index sensitivity analysis would strengthen the claim. |
| Histotype, rather than centre, is the dominant source of RNA variation. | **Overstated** | Histotype retains unique explanatory signal after adjustment and in patient-representative sensitivity analysis. However, all high-grade serous models come from one centre and rare histotypes from another, so source/culture effects aligned with histotype cannot be excluded. |
| Known histotype markers and programmes are recovered. | **Partly supported** | Marker-direction recovery is stronger than a correlation-preserving permutation null. After patient collapse, 15/25 markers occupy the expected extreme categories (joint permutation *p* = 0.0052). Only the clear-cell and SCCOHT gene-set programmes remain convincingly multiple-testing significant; the other programmes should not be presented as equivalently recovered. |
| Proteomic bridge data show 15.7–20.4% repeatability CV and technical noise 2.4–3.1 times biological spread. | **Not supported as stated** | The bridge differences show useful but nontrivial repeatability limits. The published conversion treats the SD of a pairwise difference as a single-measurement SD. After division by √2, the approximate per-measurement lognormal CV is 10.4–13.2%. Comparing 95% limits-of-agreement width with a cross-model IQR is not a like-for-like noise-to-signal ratio. |
| A 3.34-fold TMT ratio compression sets the RNA–protein concordance ceiling. | **Over-attributed** | The proteomic matrix has a narrower cross-model dynamic range than RNA, consistent with TMT ratio compression. Protein turnover, homeostasis, normalization, missingness, and platform scale also contribute; these data cannot isolate a 3.34-fold causal compression factor or a concordance ceiling. |
| RNA–protein concordance is biologically and technically credible. | **Supported as descriptive validation** | Median per-model Spearman correlation is 0.408 across 30 paired models and median per-gene correlation is 0.397 across 7,894 genes. These are descriptive concordance estimates, not proof of accuracy for individual proteins. |
| Tumour-only exome filtering recovers canonical ovarian cancer genetics. | **Supported for selected positive controls only** | The filtering cascade removes many likely germline/background calls and retains *TP53* alterations in all 17 high-grade serous models/11 patients. Retained calls are still “coding candidates,” not confirmed somatic variants, and burdens/frequencies should not be interpreted as tumour mutation burden. |
| Copy-number profiles are quantitatively reliable. | **Provisionally supported for coarse patterns** | Canonical high-grade serous instability and low low-grade serous alteration are recovered, and threshold sensitivity is reported. Unknown capture-kit/BED compatibility with the five public pooled normals prevents strong claims about quantitative FGA or arm-frequency precision. |
| `TOV21G` is MSI-high/MMR deficient. | **Candidate only** | Burden, indel composition, signature fitting, and DepMap ordering converge on a strong candidate. Tumour-only calls, exome capture, lack of exome-to-genome normalization, correlated COSMIC signatures, and absent orthogonal testing prevent attribution. MMR IHC and MSI-PCR are required. |
| Cell-model identity has been authenticated. | **Not for the full panel** | Five models self-match their DepMap namesake at rank 1/67. Thirty of 42 have a public Cellosaurus STR profile, but that documents a reference profile rather than matching the present stocks. No in-house STR comparison or mycoplasma testing is available. |
| The data independently recover two published reclassifications. | **Partly incorrect** | The relevant paper re-assigned `COV434` as SCCOHT and `TOV112D` as dedifferentiated ovarian carcinoma; `BIN67` is an established SCCOHT model, not a reclassification in that study. The present profiles are molecularly consistent with those assignments but do not independently perform the pathology-based reclassification. See [Karnezis et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33328126/). |

## Major methodological findings

### 1. Centre–histotype confounding limits causal interpretation

The commonality analysis is a useful improvement over a simple PCA colour plot. For RNA PC1, unique histotype contribution remains substantial (raw 0.424; adjusted 0.393), unique centre contribution is near zero (raw 0.002; adjusted −0.025), and shared contribution is 0.311; the permutation result for unique histotype is significant. Patient-representative analysis similarly retains a large unique histotype component.

That does not make centre irrelevant. A shared variance component of 0.311 is large, the top RNA dendrogram split follows centre, all high-grade serous models are from CHUM, and several rare histotypes are from BC Cancer. The within-clear-cell comparison has only two CHUM and five BC Cancer models and is underpowered. Statistical adjustment cannot identify whether a predictor-aligned effect is histotype biology, tissue handling, growth conditions, passage history, extraction, library preparation, or another site-linked factor.

**Required manuscript change:** replace “histotype, not batch, drives the transcriptome” with language such as:

> Histotype labels are associated with the major transcriptomic axes and retain unique explanatory signal in adjusted and patient-representative analyses. Because histotype coverage is imbalanced by source centre, source- and culture-associated effects aligned with histotype cannot be excluded.

### 2. The patient, not the derived line, is the inferential unit

The overall panel contains 42 models from 34 patients, including 13 sublines from five donor families. In the RNA cohort, 31 models represent 28 patients. Treating all 31 line models as independent artificially increases the precision of every one-versus-rest comparison, including groups that do not themselves contain sublines, because duplicate donors also enter the reference set.

The existing sensitivity analysis demonstrates the practical consequence: collapsing to one RNA model per patient reduces high-grade serous from 15 to 12 observations and materially reduces differential-expression counts across several contrasts. Clear cell remains seven observations but loses 435 of 1,256 significant genes because the “rest” group changes.

**Required manuscript change:** use the 28 patient representatives as the primary unit for inferential marker, PCA-permutation, and differential-expression validation. Retain 31-model analyses as descriptive resource views. If all sublines are scientifically valuable, a donor-aware mixed model is another option, but with only three replicated RNA donor families it will not estimate donor variance precisely.

### 3. Marker validation survives a more appropriate null model

The current exact-binomial test treats the 25 literature markers as independent Bernoulli trials. They are not independent: markers within a histotype can be co-regulated, and a label permutation can move several markers together. I therefore reran the directional recovery test with 20,000 joint permutations of histotype labels while preserving the expression matrix and marker correlation.

| Analysis unit | Observed markers in expected top/bottom two groups | Permutation-null mean | Null 95th percentile | Empirical *p* |
|---|---:|---:|---:|---:|
| 31 line models | 16/25 | 6.68 | 11 | 0.0019 |
| 28 patient representatives | 15/25 | 6.81 | 12 | 0.0052 |

The qualitative conclusion is strengthened: marker recovery is unlikely under jointly permuted labels. The corrected patient-level *p* value is less extreme than the binomial *p* = 0.0017, as expected when marker dependence is preserved.

**Required workflow change before publication:** incorporate this joint-permutation test into a numbered analysis script, fix the seed, and deposit the null summary. Until that occurs, treat the values above as an independent audit result rather than a canonical pipeline output.

### 4. Gene-set “recovery” is uneven

The pathway search selects a favourable term for each regex-defined programme from a larger GO result set. Only the clear-cell and SCCOHT programmes remain Benjamini–Hochberg significant in the patient-representative analysis. High-grade serous and mucinous are nominal/suggestive, and MMMT is not recovered. Calling all five programmes recovered conflates confirmatory and exploratory evidence and does not fully account for selection of the best matching term.

**Required manuscript change:** report clear-cell and SCCOHT as robust pathway checks. Describe the remaining patterns as exploratory, or omit them from the main Technical Validation.

### 5. Proteomic repeatability should remain on the log2-difference scale

For an adjacent-plex bridge pair, the SD of the log2 difference contains error from both measurements. If the two measurements have equal, independent error variance, the per-measurement SD is:

`SD_single = SD_difference / sqrt(2)`

The current `100 × (2^SD_difference − 1)` calculation produces 15.7–20.4% and labels it “CV.” Applying the replicate-difference adjustment and the standard lognormal conversion,

`CV = 100 × sqrt(exp((ln(2) × SD_single)^2) − 1)`,

gives approximate per-measurement CVs of 10.4–13.2%. This remains an estimate based on four non-identical bridge samples rather than replicated aliquots of one common reference.

More importantly, the 95% limits-of-agreement span (0.826–1.052 log2) is a two-sided interval width, whereas the reported biological 0.344 log2 value is a median cross-model IQR. Their ratio does not show that “technical noise exceeds biological spread.” A like-for-like rough comparison instead gives a median observed cross-model protein SD of about 0.287 log2 versus a per-measurement bridge SD of 0.149–0.190 log2, an observed cross-model-to-estimated-technical SD ratio of approximately 1.5–1.9. The cross-model SD itself still contains technical variation and is not a pure biological variance estimate.

**Required manuscript change:** report bridge mean differences, SD of differences, and 95% limits of agreement directly. Do not label the existing conversion as a CV, and remove the 2.4–3.1-fold claim.

### 6. Narrower protein spread is multifactorial

The median protein-to-RNA IQR ratio of 0.299 is a clear property of these processed matrices. It does not isolate TMT ratio compression because RNA and protein are measured on different scales, with different normalization, detection, biological regulation, and missingness. The result is “consistent with known TMT ratio compression,” not proof that compression is exactly 3.34-fold or that it sets a hard concordance ceiling.

### 7. Proteomic histotype structure is weaker and retains plex/site structure

The protein analysis makes an appropriate choice for PCA by using the 6,855 proteins quantified in all five plexes, thereby avoiding imputation and preventing plex-structured missingness from directly determining the components. This is a sound primary complete-case analysis, though it preferentially retains widely quantified proteins and should not be treated as a view of all 8,427 matrix rows.

The resulting structure is not as uniformly biological as the summary language suggests:

- protein PC1 has adjusted R² of 0.464 for histotype and 0.315 for site;
- protein PC2 has adjusted R² of 0.404 for histotype and 0.215 for TMT plex;
- clear-cell mean silhouette is −0.003 and SCCOHT is 0.028, indicating essentially no separation from the other annotated groups in this protein representation;
- all silhouette values use line models, including same-patient sublines, and groups of two are not stable estimates.

The high bridge correlations show that protein rank order across thousands of proteins is reproducible within each re-run sample, but correlation is dominated by the large between-protein abundance range and does not by itself rule out smaller plex effects in between-model comparisons. The Bland–Altman quantities and leading-PC plex associations are the more relevant evidence.

**Required interpretation change:** the protein data show useful cross-plex comparability and some histotype-associated structure, but do not independently recover every histotype. Site and plex remain material on individual leading components.

### 8. RNA–protein concordance is descriptive, not an accuracy certificate

The two concordance views answer different questions:

- the per-model correlation (median 0.408) correlates thousands of genes within one model and is strongly influenced by stable gene-to-gene abundance differences;
- the per-gene correlation (median 0.397 across 7,894 genes) asks whether model-to-model differences agree between RNA and protein and is the more relevant reuse statistic.

The agreement between the two medians is reassuring descriptively, but it is not an “internal consistency” replication because the estimands differ. Per-gene estimates also use variable sample counts (minimum 10), structured protein missingness, and three repeated-patient families. A patient-representative sensitivity analysis would show whether same-donor sublines inflate the cross-model correlations.

The statement that 813 genes with negative point estimates are ones for which RNA “is not a proxy in any direction” is too strong without uncertainty intervals or significance/reliability filtering. With 10–30 paired observations per gene, some negative estimates will be unstable. Deposit the estimates and sample counts, but describe this as a set of candidates for discordant regulation rather than a validated negative-proxy list.

### 9. Tumour-only exomes are useful but intrinsically qualified

The filter cascade from 557,392 raw records to 6,036 retained coding nonsynonymous candidates is transparent and substantially improves interpretability. Population frequency filters and a normal panel reduce common germline and technical background, but they do not recreate matched normals. Rare germline variants, clonal culture-specific changes, and systematic artefacts can remain. The confidence tiers are heuristic and informed by prior gene- and variant-class knowledge, so they should not be presented as a validated somatic classifier.

The strongest genomic validation is the *TP53* positive control in high-grade serous models. Mutation burden, recurrent driver frequency, and inferred absence of somatic *BRCA1/2* require more cautious phrasing:

- use “retained coding candidate count,” not tumour mutation burden;
- use “candidate variant” or “high-confidence candidate,” not “somatic mutation,” unless orthogonally established;
- state that no defensible *somatic* *BRCA1/2* call can be made from this dataset, rather than that no *BRCA1/2* alteration exists.

The tier rules are sensible literature- and variant-class-based triage: hotspot or truncating variants in selected canonical genes receive higher confidence, whereas rare *BRCA1/2* calls are deliberately retained at Tier 3. They do not use each sample’s labelled histotype directly, which avoids a more serious circularity. Nevertheless, the rules are project-specific, are not a validated somatic classifier, and sometimes treat a gene class as higher confidence despite the lack of a matched normal. “Somatic-confidence tier” should therefore be defined as a prioritization heuristic.

### 10. Copy-number conclusions depend on undocumented capture compatibility

CNVkit with pooled public normals is a reasonable salvage strategy for cell-line exomes without matched normals. However, the capture kit and intervals are currently unknown. Target-set mismatch can create systematic log-ratio structure, segmentation artefacts, and arm calls. The canonical direction of FGA by histotype and threshold sensitivity support qualitative reuse, but the paper should not imply precision comparable to a matched, capture-compatible copy-number study.

Per-sample autosomal median centring also removes an unknown baseline in highly aneuploid models. Without purity/ploidy or allele-specific information, the segments are relative total-copy profiles, not absolute copy states. This is appropriate for coarse within-sample gains/losses and FGA under a declared threshold; it is not sufficient for whole-genome doubling, loss of heterozygosity, or HRD scoring. The 30-fold HGS-versus-LGS FGA contrast compares a well-sampled HGS group with one LGS model and should be presented as a panel range/positive control, not a histotype effect estimate.

### 11. Identity evidence is uneven

Expression self-matching is persuasive for the five models with DepMap counterparts: each ranks first among 67 ovarian models. It does not authenticate the remaining 37. Likewise, a Cellosaurus STR record supplies a reference profile; without STR profiling of the stocks used here, it does not demonstrate a match. The current v4 opening sentence says STR profiles and mycoplasma clearance “were obtained,” which conflicts with the project record and should be removed.

For publication, in-house STR matching and mycoplasma clearance are high-priority, preferably for all 42 models and at minimum for the 12 without Cellosaurus profiles and the provenance-flagged VOA models.

The project’s molecular consistency rules are useful as a transparent flagging system, not as authentication:

- HGS/MMMT “expression-consistent” generally means absence of a strong competing clear-cell or mucinous programme, not presence of a positive serous/carcinosarcoma programme;
- SCCOHT consistency reuses SWI/SNF loss, the same evidence used for the molecular label, and is therefore circular as an independent validation;
- clear-cell and mucinous use small, hand-selected marker panels and fixed z thresholds;
- only eight of 42 models have a positive lineage programme supporting the label.

The six-model “SWI/SNF-deficient” classification also combines project-specific cutoffs: *SMARCA4* protein z ≤ −1, *SMARCA2* RNA rank ≤4 or z ≤ −1.5, and selected Tier 1–2 truncations. This is biologically motivated and useful for prioritization, but the thresholds are not externally calibrated and the umbrella term combines loss of different subunits with different disease implications. Report the gene- and layer-specific evidence first; call newly flagged models “candidate [subunit]-loss” models pending immunoblot/IHC or orthogonal confirmation.

The mucinous ovarian-versus-GI analysis is commendably explicit that its mixed RNA/protein composite and thresholds are ad hoc. Its per-marker values are informative; the composite verdict is not a trained origin classifier. `VOA8762` and `VOA8771` should remain provenance flags requiring STR and histotype IHC, not molecular reclassifications.

### 12. Variance partitioning and passage effects are descriptive

Patient identity has 28 levels among 31 RNA models and only three replicated donor families, so a patient variance term is close to unidentifiable. The existing analysis appropriately labels this limitation in places, but the paper should keep the result supplementary and descriptive. Passage-number associations are similarly weak and are complicated by assay-specific passages that differ by −17 to +20 passages for some models.

### 13. Exploratory reuse analyses are useful when labelled as such

The ADC target atlas is a defensible usage illustration because it reports deposited expression values and already states the key translational limitation: bulk RNA and total TMT protein do not measure cell-surface abundance, antigen density, internalization, linker/payload sensitivity, or drug response. Per-line expression can shortlist models but cannot predict ADC efficacy. Subtype means for groups of two or three and RNA/protein means computed from slightly different line sets should remain descriptive. The small-sample `FOLR1` mixture/bimodality exercise is not needed for the resource claim; at 15 HGS lines, including related sublines, it is sensitive to one or two extreme models.

The within-HGS three-cluster analysis is substantially more exploratory. The value of `k = 3` is chosen for description rather than selected by a prespecified criterion; no bootstrap or consensus-clustering stability is shown; 15 lines represent 12 patients; cluster labels are assigned from the same Hallmark score space used to create the clusters; and PROGENy agreement is partial rather than a clean independent corroboration. This is acceptable as an explicitly illustrative supplement, but the named strata should not appear as established subtypes or as a main validation result.

`consensusOV` was developed for bulk tumours containing stromal and immune components. Ambiguous calls and sensitivity to the input set are expected in pure cultures. These calls do not validate model identity and are better omitted from the main descriptor or deposited as an exploratory annotation with a strong warning.

## Layer-by-layer readiness

| Layer | Methodological readiness | Appropriate reuse now | Main condition before publication |
|---|---|---|---|
| Sample roster and patient-family map | **High** | Coverage queries, donor-aware selection, provenance review | Add complete culture, derivation, consent/ethics, passage, and current-stock QC metadata. |
| RNA counts/TPM | **High with reference caveat** | Expression lookup, unsupervised exploration, patient-aware comparisons | Deposit the exact kallisto index/map provenance and add a matched-reference sensitivity analysis or re-quantification. Recover library-preparation details. |
| Protein relative-abundance matrix | **Moderate** | Relative abundance, target shortlisting, broad pathway/protein comparisons with plex awareness | Recover acquisition/search/channel-map metadata and raw-data status; correct repeatability language; preserve structural missingness flags. |
| Tumour-only mutation candidates | **Moderate for hypothesis generation; low for definitive somatic inference** | Canonical positive controls, candidate prioritization, cross-resource comparison | Deposit complete calling/filter details and capture metadata; keep all calls explicitly tumour-only/candidate. |
| Relative total-copy profiles | **Moderate for coarse patterns** | Broad gain/loss and within-model hypothesis generation | Confirm capture compatibility; describe centring/threshold sensitivity; avoid absolute-copy, LOH, and HRD claims. |
| Identity/lineage annotations | **Low to moderate** | Flagging stocks for follow-up; five-model DepMap self-match | Current-stock STR, mycoplasma results, and orthogonal validation of newly inferred lineage/SWI-SNF findings. |
| ADC and HGS usage examples | **Exploratory** | Demonstrating how a reuser might query the records | Keep clearly separate from Technical Validation; do not imply response prediction or stable molecular subtypes. |

## Statistical and reporting principles for the final resource paper

1. **State the unit next to every denominator.** “Model,” “patient,” “protein row,” “gene symbol,” “WES model,” and “paired model” are not interchangeable. The existing denominator inventory is excellent and should control all text and figure legends.
2. **Separate description from inference.** Groups with two or three models can be displayed and compared descriptively, but they cannot support generalizable subtype signatures. The `n ≥ 10` threshold used to label HGS “formal” is a project convention, not a statistical guarantee.
3. **Make patient-representative results primary for population claims.** Line-level values remain the correct resource unit for browsing and model selection.
4. **Do not interpret adjusted regression as deconfounding an unsupported design cell.** A model cannot estimate the HGS-at-BC-Cancer counterfactual when no such observations exist.
5. **Report effect sizes and uncertainty, not only feature counts.** Thousands of DE genes can arise in a high-dimensional, strongly separated dataset; effect sizes, donor sensitivity, and known-marker/pathway checks are more meaningful validation.
6. **Treat data-driven thresholds as conventions.** FGA cutoffs, hypermutator rules, z-score loss rules, lineage scores, and cluster number are useful operational definitions. Sensitivity analyses should accompany any conclusion that changes across reasonable values.
7. **Avoid causal assay language.** “Consistent with ratio compression,” “associated with histotype,” and “candidate MMR phenotype” correctly match the evidence.
8. **Do not use external references as if measured here.** Cellosaurus profiles, originator IHC, DepMap variants, and in-house omics should remain in distinct fields and sentences.

## Scientific Data genre and reporting fit

The current draft has the correct broad genre, but it needs several structural changes to match the current Scientific Data submission guidance:

- keep the title under 110 characters and the abstract at or below the recommended 170 words;
- include the required sections `Data Availability`, `Code Availability`, `References`, `Author Contributions`, and `Competing Interests`, plus `Funding` and a human-data ethics statement where relevant;
- make `Data Records` describe deposited files, identifiers, schemas, and access, not only planned deposition;
- reserve `Technical Validation` for checks of data quality and reliability rather than biological discovery;
- if retained, limit `Data Overview` to one paragraph and one or two display items;
- keep `Usage Notes` technical and practical rather than using it for a worked model-selection case study or discovery conclusions;
- provide anonymous reviewer access to the complete data at initial submission and formal public repository records no later than revision/publication.

The reviewed guidance was: [Scientific Data submission guidelines](https://www.nature.com/sdata/publish/submission-guidelines), [data policies](https://www.nature.com/sdata/policies/data-policies), [referee guidelines](https://www.nature.com/sdata/policies/for-referees), [editorial and publishing policies](https://www.nature.com/sdata/policies/editorial-and-publishing-policies), and the journal’s [artificial-intelligence policy](https://www.nature.com/sdata/policies/ai). The data policy also requires a Human Data Checklist for every human-derived dataset, even when the shared data are judged nonsensitive.

The current title is acceptable at 93 characters. The v4 abstract is 176 words and should be shortened. `consensusOV` calls on pure cell cultures, the HGS pathway-strata case study, and most “model-selection” material should be removed from the main text or placed in a clearly exploratory supplement. An ADC target-expression atlas can remain as the optional, concise Data Overview because it directly illustrates reuse of deposited values.

## Submission blockers

### Critical before submission

1. **Repository deposition and reviewer access.** Supply stable, anonymous reviewer-access URLs and accession placeholders for every claimed data record. At publication, make the records public in established repositories.
2. **Proteomics acquisition and upstream processing.** Recover the TMT channel map and details of sample preparation, digestion, labelling, fractionation, LC–MS/MS instrument/acquisition, database search, peptide/protein FDR, protein inference, and reporter-ion normalization. Deposit raw mass spectra if available. If raw data and core methods cannot be recovered, describe the protein matrix explicitly as a processed secondary dataset with limited upstream reproducibility.
3. **Cell authentication and contamination testing.** Add current-stock STR comparison and mycoplasma results or state the limitation prominently. A broad “authenticated” claim is not defensible without them.
4. **Ethics and provenance.** Provide the original human-sample consent/ethics approvals or a precise institutional determination explaining why additional review is not applicable; complete the journal’s human-data reporting requirements.
5. **Core culture metadata.** Supply media, supplements, atmosphere, temperature, passage procedures, and assay passage information. These variables matter because source centre is entangled with histotype.
6. **Exome capture metadata.** Identify the kit, target BED, read layout, reference/build details, and compatibility of pooled normals.

### Major before acceptance

7. Reframe patient representatives as the inferential unit and update marker/pathway validation accordingly.
8. Correct the protein bridge interpretation and remove the “technical noise exceeds biology” and hard “3.34-fold compression ceiling” claims.
9. Preserve `TOV21G` as a candidate only; ideally add MSI-PCR and MMR immunohistochemistry.
10. Re-quantify RNA against a transcript index matched to the deposited transcript-to-gene map, or provide a focused sensitivity analysis showing that key QC, PCA, and marker results are robust to the mismatch.
11. Add complete per-file data dictionaries, checksums, licenses, and links between raw/intermediate/processed objects.
12. Add code repository URL, release tag/commit, environment restoration instructions, and an explicit statement of any generative-AI assistance used in manuscript preparation.

## Recommended claim language

Prefer:

- “histotype-associated transcriptomic structure persisted after adjustment, although source-aligned effects cannot be excluded”;
- “15 of 25 prespecified markers occupied the expected extreme categories after patient collapse; joint label permutation *p* = 0.0052”;
- “the proteomic matrix had a narrower cross-model spread than the RNA matrix, consistent with TMT ratio compression and post-transcriptional regulation”;
- “retained tumour-only coding candidates”;
- “candidate mismatch-repair-deficient/MSI-high phenotype”;
- “five DepMap-overlapping models showed reciprocal expression self-matches; current-stock STR matching was not available”;
- “copy-number profiles support coarse, qualitative comparisons pending confirmation of capture compatibility.”

Avoid:

- “biology, not batch”;
- “all histotype programmes were recovered”;
- “technical noise exceeds biological spread by 2.4–3.1×”;
- “3.34-fold compression sets the concordance ceiling”;
- “somatic burden” or “somatic mutation frequency” for tumour-only candidates;
- “the panel was authenticated”;
- “TOV21G is MSI-high” without orthogonal validation;
- “two published SCCOHT reclassifications” when counting `BIN67`.

## Presentation and reproducibility review

The 36-page results PDF is visually clean. I found no clipped tables, overlapping elements, missing figures, or unreadable pages. Several multi-panel figures and legends are dense, but they are suitable as a technical report. For a journal submission, legends should be shortened and figure panels should be limited to validations needed to understand the deposited records.

The reproducibility framework is stronger than is typical for a mature exploratory project: pinned R environment, session records, deterministic scripts, a clean-room workflow check, post-fix output snapshot, and extensive parameter reporting. Remaining reproducibility gaps are mostly upstream experimental metadata and external resources rather than script execution.

## Bottom line

The project can become a strong Scientific Data descriptor if it is positioned as a carefully documented resource with explicitly bounded validity, not as a definitive biological classification study. RNA is the most submission-ready layer. Proteomics is reusable at the processed-matrix level once upstream acquisition/search details and raw-data status are resolved. Tumour-only mutation and pooled-normal copy-number layers are valuable for hypothesis generation and canonical positive controls but require persistent caveats. The revised manuscript should lead with coverage, deposited records, provenance, quality controls, and reuse constraints; biological findings should serve as validation examples rather than the paper’s endpoint.
