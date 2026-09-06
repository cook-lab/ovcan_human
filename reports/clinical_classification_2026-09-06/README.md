# Clinical and molecular annotation opportunities

Assessment: 6 September 2026. Current manuscript: v9. This exploration adds model-selection evidence; it does not change the processed release or assign treatment eligibility to cell lines.

## What can be added now

The resource supports a useful annotation layer separating **variant significance**, **allele/copy-number state**, **expression**, **genomic scars**, and **functional phenotype**. These are different measurements. In particular, the existing research Tier 1–3 system prioritises candidate variants and uncertainty of somatic origin; it is not an AMP/ASCO/CAP clinical classification. A germline pathogenic allele can matter in a model even when tumour-only sequencing cannot establish its origin.

The strongest immediate additions are exact curated pathogenic/oncogenic allele annotations; CCNE1 relative copy number and RNA abundance; and model-specific leads for MMR, ERBB2, FOLR1 and pathway alterations. Clinical amplification, HRD, MSI, IHC positivity and drug response require additional evidence.

## CCNE1: concrete model leads

The values below come from the corrected target-only CNV segments and the release-93 RNA matrix. Gene coordinates are GRCh38 chr19:29,811,991–29,824,312, converted to 0-based half-open intervals for overlap. Every WES model has **seven retained, positive-depth target bins overlapping CCNE1**. Those bins provide evidence independently of the larger segment average. [NCBI Gene 898](https://www.ncbi.nlm.nih.gov/gene/898/).

| Model | Gene-overlapping segment log2 ratio | CCNE1 target-bin median log2 ratio | RNA evidence | Interpretation |
| --- | ---: | ---: | --- | --- |
| TOV2929D | 3.044 | 2.619 | Not measured | Strongest relative DNA gain; prioritise absolute-copy validation and expression measurement |
| OV2085 | 2.239 | 2.414 | 124.43 TPM; second of 15 HGSC models | Strong same-model DNA/RNA lead, measured at different passages |
| TOV2881EP | 1.352 | 1.336 | Not measured | Elevated DNA; only 0.278 log2 above its 19q median, so much is broad arm gain |
| TOV3133D | 1.242 | 1.457 | Not measured | Relative DNA gain; useful within-patient comparison with the other 3133 models |
| OV3291 | 1.241 | 1.560 | Not measured | Localised relative gain, 1.361 log2 above its 19q median |
| TOV2835EP | 0.449 | 1.198 | Not measured | Additional gene-bin signal diluted by the larger segment; inspect local depth and validate before classifying |
| TOV3291G | Not measured | Not measured | 158.02 TPM; first of 15 HGSC models | Highest HGSC RNA; distinct from the WES-profiled OV3291 isolate |

The first five pass an explicitly exploratory screen of gene-overlap segment log2 ≥1. This selects follow-up candidates; it is not a validated amplification threshold. TOV2929D and OV2085 correspond to approximately 8.25-fold and 4.72-fold ratios to their fitted autosomal baselines. Those ratios are **not absolute copy counts**. Their locus-minus-19q-median contrasts are also large (3.262 and 1.946 log2), so their ranking is supported by local contrast as well as autosomal centring. The segments span approximately 4.21 and 4.97 Mb; these data do not establish CCNE1 as the only gene driving the gained region.

OV90 has a moderate relative gain (0.694 segment log2; target-bin median 0.494), which should not be presented as the same strength of evidence as the top two. No WES measurement in a model means unknown copy number, not a negative result. Patient-related models retain separate assay entries.

Clinical/tumour studies use assay-specific CCNE1 amplification and Cyclin E1 protein criteria. For example, Chan et al. evaluated >8 copies by CISH together with an IHC threshold, neither of which is measured by this centred WES ratio or bulk TMT abundance. We therefore recommend “candidate CCNE1 amplification” or “strong relative CCNE1 gain” until absolute copy number or an orthogonal assay is available. [Chan et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7578325/).

### RNA and protein must remain separate

TOV3291G and OV2085 remain first and second after reducing HGSC RNA to 12 patients. VOA12539 is also RNA-high (146.95 TPM, first of seven RNA-profiled clear-cell models), but has no WES. The exact RNA values, model ranks, denominator counts and patient-weighted reference values are in the accompanying table. There is no matched normal expression baseline or calibrated clinical overexpression cutoff.

CCNE1 protein is quantified in all 31 protein-profiled models. Its `isoDoping` flag denotes targeted detection enhancement, not automatic invalidity: independent mapping to the original workbook confirms all 31 imported values exactly and excludes the SM+iD channels. However, the actual correction/rollup protocol remains undocumented. The five global unique peptides are not five per model: three are measured throughout, with two additional peptides confined to plex 2.

The supplied protein ranks do not corroborate the strongest RNA ranks: OV2085 is 20th of 31 and TOV3291G 13th of 31. RNA/protein Spearman correlation is 0.123 across 30 same-model pairs. Thus “CCNE1 RNA-high” is supported for these models, whereas protein overexpression has not been established. The mismatch cannot be attributed specifically to isoDoping, TMT interference or biological regulation from these records. [Independent expression and peptide review](CCNE1_EXPRESSION_REVIEW.md).

OV2085 has WES at P64 and RNA at p68; protein passage is undocumented. OV3291 WES P46 and TOV3291G RNA p62 come from distinct same-patient isolates. Do not transfer the OV3291 CNV result to TOV3291G.

## HRD and BRCA/HRR annotations

No validated HRD score can be calculated directly from the present total-copy segments. LOH and telomeric allelic imbalance require allele-specific states; fraction of genome altered and SBS3 refitting cannot replace them. A future WES workflow can produce exploratory genomic-scar measurements, but the current VCFs were not generated with systematic germline/PoN-site genotyping. Retained common SNPs alone do not prove an unbiased SNP input. The preferred next step is input recovery/preparation and a small, reviewed tumour-only allele-specific CN pilot, followed by scar scoring only for credible fits. [Detailed feasibility and input counts](HRD_FEASIBILITY.md).

**TOV81D is a concrete pathogenic-allele finding:** its VCF supports BRCA2 NM_000059.4:c.8537_8538del (p.Glu2846fs), ClinVar Variation 9328, classified Pathogenic by the ENIGMA expert panel. This refines the rough archived protein label p.R2845fs. The source VCF has reference/alternate read counts 75/45, and caller AF 0.392. AD-derived fraction (45/120 = 0.375) and model-estimated AF are different quantities. The legacy Tier 3 rationale concerns unconfirmed somatic origin; it must not be read as a benign or uncertain-pathogenicity verdict. [ClinVar VCV000009328.54](https://www.ncbi.nlm.nih.gov/clinvar/variation/9328/).

This is not evidence that TOV81D is BRCA2-null or HR-deficient. A primary historical study of TOV-81D reports retention and expression of the wild-type allele alongside the 8765delAG allele. Its current low relative FGA is also not a functional HR assay. Current-stock allele-specific status and functional RAD51 testing would resolve different questions. [Allelic transcript dosage study, 2012](https://aacrjournals.org/cancerpreventionresearch/article-abstract/5/5/765/32041).

The 3133 family has retained CDK12 frameshift candidates, and OV3291 a splice-site candidate; OV2295 has ATM candidates. These warrant exact allele/read review and second-hit assessment, not automatic pooling into a BRCA-like HRD group. Genomic scars can persist after HR repair restoration, so functional competence and treatment response should be recorded separately.

## Other useful classification axes

| Axis | Immediate model lead | What is supportable now | What remains needed |
| --- | --- | --- | --- |
| MMR/MSI | TOV21G | Prioritise its existing mutation/signature/literature evidence | Validated MSI analysis or locus assay and MMR-protein assessment; no MSI-high call from candidate count alone |
| ERBB2/HER2 | TOV3392D for expression; TOV2835EP for DNA | TOV3392D ranks first for ERBB2 RNA and protein. TOV2835EP has a local relative gain (segment 0.985; target-bin median 0.845 log2) | HER2 cell-surface/IHC confirmation and absolute-copy assessment; TOV3392D's broad segment is not gene-specific amplification evidence |
| FOLR1/FRalpha | OV2295, TOV3133G, OV3133-R; VOA14993 protein-only | Rank measured RNA/protein and retain their disagreement/missingness | Membrane expression and assay-specific positivity; bulk abundance is not the clinical IHC score |
| MAPK/PI3K | KRAS/PIK3CA variants in TOV21G; KRAS variants in TOV112D/TOV2414/TOV3392D; OV90 BRAF deletion | Exact allele-level model annotations; independently reported non-V600 BRAF deletion in OV90 | Variant-specific evidence and intended experimental context; do not relabel a non-V600 allele as BRAF V600E |
| SWI/SNF | Existing ARID1A/SMARCA4/SMARCA2 panels | Retain cross-assay evidence, candidate truncations and source-stock literature | Distinguish protein loss from sequence variation and from current-stock authentication |
| Endocrine / surface targets | ESR1-high RNA in OV2295/OV4485; PGR protein unavailable | Exploratory expression leads | Receptor localisation, protein assay and functional dependence |
| TMB / POLE-POLD1 / fusions | No new validated category assigned | Document present candidates and assay coverage | Callable-region and somatic filtering validation for TMB; pathogenic polymerase allele/signature support; read-level fusion analysis where justified |

The source-supported rationale, exact current clinical assay boundaries, primary OV90 BRAF evidence and additional model details are in [classification options](CLASSIFICATION_OPTIONS.md). This table does not assert treatment response or add an unmeasured clinical category.

## Recommended addition to the Data Descriptor

Add a supplementary **model-selection annotation table** with explicit evidence levels rather than one flat clinical label. Suggested fields are: exact model/patient, assay and passage, variant HGVS and source/version, pathogenicity/oncogenicity, unconfirmed origin where applicable, second-allele status, relative versus absolute CN, RNA/protein measurements, functional-test status, and verification needed. Keep “not measured”, “not established” and “no retained candidate” distinct.

The current exploration provides this starting evidence without rewriting manuscript v9 or the validated release. HRD remains a worthwhile extension, but its score and any binary research category should enter the paper only after an allele-specific fit and sensitivity review. No clinical HRD threshold is transferred directly to these cell lines. Orthogonal CCNE1 copy-number/protein testing and functional HR testing would be more informative than adding several weak proxy scores.

## Files and reproduction

- [All 42 models: CCNE1 DNA/RNA/protein](ccne1_model_summary.csv).
- [CCNE1 and ERBB2 loci](cnv_locus_summary.csv), [per-target evidence](locus_target_bin_evidence.csv), [coordinate sources](locus_definitions.csv).
- [36-gene RNA/protein panel](marker_expression.csv); 1,512 rows, retaining absent measurements.
- [57 retained variants in that panel](biomarker_variant_candidates.csv); inherited research tiers are not clinical tiers or fresh ClinVar classifications.
- [Curated variant significance](curated_variant_classifications.csv), [HRR VCF evidence](hrr_candidate_vcf_evidence.csv), [HRD input counts](hrd_vcf_input_counts.csv).
- [CCNE1 figure](../../output/pdf/ccne1_exploration.pdf) and [external caption](CCNE1_FIGURE_CAPTION.md).
- [Cluster follow-up plan](../../docs/cluster/CLINICAL_CLASSIFICATION_NEXT_STEPS.md).

Run `python3 scripts/40_clinical_biomarker_exploration.py` with numpy/pandas and the restored archived CNR files, then `Rscript scripts/41_ccne1_exploration_figure.R` with the current figure environment. `OVCAN_PROJ` and `OVCAN_DATA` are supported. Script 40 hashes source matrices/CNRs, checks all 23 current target-only CNS hashes and coordinate coverage, and writes separate exploratory outputs. It does not call variants, refit CNV, impute expression or mutate canonical data. The single curated ClinVar record and literature review are explicitly reviewed annotations, not outputs of a general clinical classifier. `analysis_validation.json` records the numerical source hashes; `review_validation.json` records independent and visual checks.
