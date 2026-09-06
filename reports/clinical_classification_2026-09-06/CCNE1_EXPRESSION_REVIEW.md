# CCNE1 expression evidence review

Date: 6 September 2026. Scope: independent read-only review of current RNA/protein matrices, source proteomics workbooks and model/patient metadata. No canonical analysis, metadata or biological calls were changed.

## Interpretation for the exploratory classification

TOV3291G and OV2085 have the highest CCNE1 RNA abundance among the 15 measured HGSC models and remain first and second after reducing these models to 12 patients. Their supplied protein values are measurable but do not independently corroborate a protein-high classification. The `isoDoping=TRUE` annotation should remain visible; it is an assay annotation, not sufficient evidence that the protein measurement is invalid or that synthetic peptide abundance has been assigned to a biological model.

The strongest same-model combination among these expression leads is OV2085: the current separate CNV analysis reports a strong relative CCNE1 gain and its RNA ranks second within HGSC. WES and RNA were measured at different passages. The OV3291 CNV result and TOV3291G RNA result belong to distinct models from the same patient; combining them is patient-related evidence, not same-model multiomic confirmation. Neither the expression ranking nor relative CNV establishes a clinical amplification or protein-overexpression category.

## RNA ranking with explicit patient denominators

CCNE1 maps to one retained gene, ENSG00000105173, in the corrected release-93 RNA import. Rankings below use supplied TPM from `output/rna_tpm.csv`, descending within HGSC. The patient analysis uses the existing `selected_patient_model` flag in `metadata/resource_models.csv`; it does not select the highest-expressing subline per patient.

| Selected HGSC model | CCNE1 TPM | RNA rank among 15 HGSC models | RNA rank among 12 selected patient models |
| --- | ---: | ---: | ---: |
| TOV3291G | 158.022 | 1 | 1 |
| OV2085 | 124.433 | 2 | 2 |
| OV4453 | 95.830 | 3 | 3 |
| OV4485 | 75.292 | 4 | 4 |
| OV866-2 | 65.882 | 5 | 5 |
| OV2295 | 63.682 | 6 | 6 |
| TOV3041G | 39.520 | 7 | 7 |
| OV3133-R | 28.623 | 9 | 8 |
| OV3331 | 25.002 | 11 | 9 |
| OV90 | 22.298 | 12 | 10 |
| OV1369-R2 | 20.889 | 13 | 11 |
| OV1946 | 11.168 | 15 | 12 |

The median across these 12 representatives is 51.600895 TPM. TOV3291G and OV2085 are respectively 3.062 and 2.411 times this descriptive median. A sensitivity check averaging TPM across all available HGSC models within each patient preserves the first and second positions, as well as the full patient ordering. That alternative median is 42.1936035 TPM. These baselines are distinct from a median of patient means on the log2(TPM+1) scale used in the accompanying all-model classification table; their values must not be substituted for one another.

Across all 31 RNA models, TOV3291G ranks first, VOA12539 (clear-cell histotype) second at 146.954010 TPM, and OV2085 third. Thus the HGSC ranks are not all-histotype ranks. No RNA-high cutoff was fitted, and related models were not counted as independent patients for the patient summary.

## What the isoDoping annotation means

The published method from the Morin/Negri group adds synthetic proteotypic peptides to an isobaric channel to increase precursor selection and improve detection of specified proteins. Biological reporter channels remain separately quantified. That implementation also describes isotope-impurity and isodoped-peptide corrections. It therefore supports interpreting isoDoping as targeted detection enhancement, rather than a blanket measurement-failure flag. [Asleh et al., 2022, *Nature Communications*, Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8850446/).

The local `tmt.layout.xlsx` labels PIS as channel 1 and a separate `SM+iD` channel as channel 11 in every plex. The cited publication doped its PIS channel, so its exact correction formula and acquisition settings cannot be assumed to describe these supplied workbooks. Script 05 explicitly excludes the five `SM.iD` columns from the biological matrix. All 31 CCNE1 biological values were independently mapped by layout back to `protein_relative_abundance.xlsx` and match exactly, with maximum absolute difference zero. There is no observed channel-import error.

The local proteomics README contains no usable method description. The workbooks do not establish which individual CCNE1 peptides were synthetic targets, the applied interference corrections, raw reporter signal precision, or the exact protein scaling/rollup implementation. Consequently, the flag alone neither proves nor rules out a quantitative bias. It also cannot identify the cause of the observed RNA/protein disagreement. Keep the supplied protein values as descriptive measurements and retain the protocol-confirmation request.

## Protein evidence and peptide support

| Quantity | Verified result | Interpretation |
| --- | --- | --- |
| Protein identifier | CCNE1 / P24864 | One source protein row |
| Biological coverage | 31 of 31 models; all five plexes | Complete supplied feature |
| Source peptide counts | Five peptides; five unique; five quantified | Counts over the source dataset, not five values per model |
| Peptides quantified in every model | Three | Two additional peptides occur only in plex 2 |
| Source `qvalue` | 0 as stored | Does not establish literally zero identification error |
| Source `CV replicates` | 10.6048 | The precise replicate and calculation definition is unavailable; not a per-model standard error |
| Source `isoDoping` | TRUE | Retain as an assay flag |
| Protein minimum–maximum | 12.336876–13.171516 | Supplied log2 normalised abundance |
| Protein median | 12.803085 | Across all 31 measured models |
| Protein range | 0.834639 log2 units | About 1.783-fold on the supplied scale |
| TOV3291G protein | 12.871291; rank 13/31 overall, 7/15 HGSC | Does not match its first-place RNA ranking |
| OV2085 protein | 12.750837; rank 20/31 overall, 11/15 HGSC | Does not provide protein-high corroboration |
| Highest supplied CCNE1 protein | TOV3133G, 13.171516 | Its RNA is 38.151979 TPM, rank 8/15 HGSC |

Scanning all 146,831 rows of `peptide_ratio.xlsx` yielded exactly five P24864 peptide rows. DQHFLEQHPLLQPK, DSLDLLDKAR and GVADEDAHNIQTHR have values in all 31 models. DTMKEDGGAEFSAR and WMVPFAMVIR have values only in OV2085, OV4485, TOV1369, TOV3291G, VOA12539 and VOA295, the six biological models in plex 2. This is relevant when interpreting a global five-peptide annotation as sample-specific evidence.

For the three complete peptides, Spearman correlations with the supplied protein row across 31 models are 0.752, 0.694 and 0.576, respectively. Their correlations with RNA across the 30 same-model pairs are 0.087, −0.077 and −0.012. These checks do not reveal a consistent peptide-level RNA signal concealed by the protein rollup. Peptide-level values are retained on their supplied scale; a filename containing “ratio” is insufficient documentation to equate this scale with the offset protein abundance values.

The stored protein QC row is at `output/prot_qc.csv:1853`. Its other source fields (`Sum PEP=32.424`, `Average SN=12.84845080819294`) lack sufficient local definitions to translate into error probabilities or raw reporter-signal thresholds; no such conversion was made.

## Cross-assay spread and concordance

For the 30 same-model RNA/protein pairs, independent Spearman correlation is 0.1230256. The existing representative-per-patient output gives 0.1324786 across 27 patients. Restricting to HGSC gives 0.0857143 across 15 models and 0.1188811 across 12 selected patient models. These are descriptive associations without adjustment for centre, plex, passage or histotype; they do not establish regulation, an assay defect or predictive validity.

For these same 30 pairs, `output/prot_dynamic_range.csv:1023` reports RNA IQR 1.6871317 on log2(TPM+1), protein IQR 0.2483911 on its supplied log2 normalised scale, and protein/RNA IQR ratio 0.1472269. The corresponding ranges are 4.4235801 and 0.8346394. The protein values vary less across this panel, but the IQR ratio is a transformation-dependent comparison of observed biological/technical spread. It is not a calibrated TMT ratio-compression factor and cannot be attributed to isoDoping from these records alone. The all-31 protein IQR (0.2306995) differs because its denominator includes the protein-only model.

## Model identity and passage limits

| Model | Patient | WES | RNA | Protein |
| --- | --- | --- | --- | --- |
| OV3291 | 3291 | Passage P46 | Unavailable | Unavailable |
| TOV3291G | 3291 | Unavailable | TOV3291Gp62, passage p62 | Plex 2/channel 8; passage undocumented |
| OV2085 | OV2085 | Passage P64 | OV2085p68, passage p68 | Plex 2/channel 2; passage undocumented |

`metadata/samples.csv` explicitly describes OV3291 and TOV3291G as distinct OV/TOV isolates. Their shared patient identifier supports grouping for patient counts, not substitution between assays. In particular, the OV3291 CNV category must not be assigned to TOV3291G. Likewise, RNA passage must not be copied into an undocumented protein-passage field. The CCNE1-high-CNV leads TOV2929D, TOV2881EP and TOV3133D have no RNA/protein measurements in this resource; missing expression is not evidence of low expression.

## Recommended fields and wording

Retain exact model and patient identifiers, assay availability, RNA/WES passages, unknown protein passage, TMT plex/channel, CCNE1 identifiers, raw TPM and supplied protein values, descending ranks with explicit denominators, and the existing isoDoping/unique-peptide annotations. If peptide support is surfaced, distinguish the global count (five) from quantified peptides per model (three, or five in plex 2).

Suitable interpretations are: “OV2085: strong relative CCNE1 gain with high same-model RNA abundance at a different passage; supplied protein does not independently confirm overexpression”; and “TOV3291G: highest HGSC RNA abundance; a distinct model from the same patient, OV3291, has a relative CCNE1 gain.” Expression ranks can prioritise follow-up, but do not establish clinical protein positivity or drug response. This isoDoping interpretation also applies to other flagged targets such as ERBB2 and ESR1: retain the flag and assess the actual measurements rather than discarding the entire feature class.

## Reproducibility evidence

Scratch extraction tables and the five complete peptide source rows are in `tmp/clinical_classification_2026-09-06/`; these are review aids, not replacements for canonical outputs. Local checks used the bundled Python/pandas/openpyxl runtime, exact layout joins, descending average ranks, and Pearson correlation of ranks after pairwise missing-value removal. No values were imputed.

| Primary input or current matrix | SHA-256 |
| --- | --- |
| `judy_archive/data/proteomics/protein_relative_abundance.xlsx` | `a454b15b8cc6bca9ef9ac1921d8ab271bb30ef248c34f13883953a9afffabe9c` |
| `judy_archive/data/proteomics/peptide_ratio.xlsx` | `7c0034f94ab7f395c1730aef532060edacbc030107b4e17fb0a1b370918a0d41` |
| `judy_archive/data/proteomics/tmt.layout.xlsx` | `e260a7478dad5c07d89c38d37e83a0a06447c499ff7b835df7f38976fdd8e264` |
| `output/rna_tpm.csv` | `ccd06b0bfd6350618a23d78f0fdaa7859ca06aa6d92aea67943ae0a4e4e55975` |
| `output/prot_abundance_matrix.csv` | `1f8484522eca3fc80891183204823012d43284e1c6562e262e5eb68c2db8a3c9` |
| `metadata/samples.csv` | `e9403bd8e32b425bcc8ae1a57c2847331102d1e1aa562ee9c34c632c68650e8d` |
| `metadata/resource_models.csv` | `60962a78bbdebdda2f5744eaff679764218e124de45a73bbf218d972f41ef7c7` |
