# Analysis audit and manuscript revision

Audit date: 5 September 2026. Target: a Scientific Data Data Descriptor.

## Assessment

The resource supports a Data Descriptor centred on assay coverage, reusable matrices, documented provenance, and technical and biological consistency checks. The audit found two recoverable data omissions, several statistical and interpretive problems, and inconsistencies between analysis outputs, the manuscript, and the browser. The principal data omissions have been repaired: all quantified RNA targets now map to the matching transcript reference, and all 23 exome-profiled models now have both variant and copy-number records.

The revision preserves the distinction between models and patients, descriptive patterns and statistical evidence, measured features and inferred phenotypes, and historical labels and current biological interpretation. Experimental details that require collaborator confirmation remain explicit author actions. The package is not yet a repository submission.

## Scope and approach

Three parallel reviews covered RNA/proteomics/integration, genomics/provenance, and the v5 manuscript. The coordinating review checked figures, release metadata, browser consistency, and cross-file claims. The Word document was treated as authoritative: its tracked edit and 34 anchored comments were extracted before revision, rather than assuming the adjacent Markdown file contained the latest author input.

The work included code inspection, direct examination of archived and processed data, targeted numerical checks, reruns of affected downstream analyses, source verification, and visual review of regenerated figures and the revised Word document. It did not rerun kallisto or raw-read WES alignment/calling, nor did it reinterpret literature protocols as confirmed experimental records.

## High-priority corrections

| Finding | Consequence | Correction and evidence |
|---|---|---|
| The release-105 transcript map did not match the kallisto target set. | 3,529 quantified transcripts and 1.60–3.41% of sample TPM were omitted from gene summaries. | Recovered the Ensembl release-93 cDNA FASTA; its official Ensembl checksum matches. All 185,299 versioned target IDs and sequence lengths match every library. Derived a pinned map and rebuilt gene matrices and dependent analyses, with explicit estimated-count and TPM conservation checks. |
| Mutation loading scanned only existing MAF files. | TOV3121D was incorrectly treated as lacking a mutation record. | Recovered its archived annotated VCF with the same vcf2maf version. Conversion was validated by exactly reproducing 29,323 rows of another model's archived MAF for key annotation fields. Added 158 retained candidates; the original 6,036 retained candidates are unchanged. |
| Chromosome-arm summaries assigned entire segments by midpoint. | Segments crossing centromeres were attributed to one arm, distorting assessed lengths and some frequencies. | Intersected segments with pinned hg38 chromosome-arm intervals, excluded centromeric bases, and recalculated primary and threshold-sensitivity summaries. FGA uses its existing segment-length calculation and is unchanged by this arm fix. |
| Line-level significance tests did not consistently respect related sublines and centre structure. | Nominal inferential support could exceed what the sampling design warrants. | Added patient-based RNA–protein correlation inference, patient-based reduced-model permutation checks, and centre-restricted marker permutations. Model-level summaries remain available for descriptive model comparison. |

## Medium-priority corrections

| Finding | Correction |
|---|---|
| The draft compared an adjusted unique histotype component with an unadjusted patient-based component; a negative adjusted component read as a negative variance share. | Use comparable raw commonality decompositions in the main figure and clearly distinguish adjusted fit penalties from variance shares. |
| High VAF alone promoted CDKN2A A86D in OV90 and implied LOH. | Downgraded the call to unresolved Tier 3 and removed VAF-only LOH/pathogenicity language. |
| Main-figure patient bars included Tier 3 although the prose excluded Tier 3 from supported frequency statements. | Retained all tiers in the alteration matrix but restricted the adjacent patient counts to Tier 1–2 calls. |
| A published RNA–protein reproducibility estimate was shown as a dataset-specific correlation ceiling. | Removed the ceiling and acceptance-band graphics; correlations describe concordance, not measurement accuracy or a demonstrated TMT ceiling. |
| Mixed-model patient terms were weakly identified, and one sensitivity analysis pooled unrelated singleton patients as if related. | Removed that invalid sensitivity and moved emphasis in the main figure to interpretable histotype/centre commonality. |
| Exploratory within-HGSC clusters could read as stable biological strata. | Added a patient-selection sensitivity showing unstable assignments and kept these clusters out of the main Technical Validation claims. |
| The public-facing cohort required an inclusion filter on a mixed internal inventory. | Created an explicit 42-model release table, with patient relationships and assay availability checked against the packaged matrices. |
| The local viewer used the obsolete reference, selected one highest-mean gene per symbol, and labelled kallisto data as Salmon data. | Rebuilt it with the matched reference and the analysis convention of summing TPM across genes sharing a symbol; retained all contributing gene IDs and added source hashes. Public deployment remains a separate author action. |
| Report-style notes crowded the manuscript figures and obscured some data. | Regenerated the main figures with the existing Cook Lab colours, Arial typography, simpler labels, and methods in the legends. |

## Interpretation retained and bounded

The resource contains 42 models from 34 patients. RNA and protein each cover 31 models; their intersection contains 30 and their union 32. Thirteen models have RNA, protein, and exome data. The variant and copy-number layers now both cover 23 models.

The new variant cascade contains 582,474 input records, 16,081 PASS records, 15,995 population-filtered records, and 6,194 retained coding candidates. TP53 alterations are present in all 18 HGSC models with exome data, representing 11 patients. These are tumour-only candidates; the positive control does not prove somatic origin for every retained alteration.

A general claim that all non-HGSC models have fewer copy-number alterations is unsupported: clear-cell TOV3392D has autosomal FGA 0.671, above the HGSC model median of 0.635. The main HGSC chromosome-arm examples remain present after interval correction, but their thresholds, full-arm denominators, and patient aggregation rule remain part of the result.

The exact RNA reference produces 39,733 gene rows and 22,678 genes after the count filter. All transcript estimated counts and TPM are conserved to numerical tolerance. The reference includes alternative-locus transcripts; this composition is recorded rather than hidden by a new gene-exclusion rule. The original binary kallisto index checksum is not available, so the verified evidence is a complete match of target identifiers and lengths plus the official reference FASTA checksum.

Median gene-wise RNA–protein Spearman correlation is approximately 0.400 across 8,033 eligible genes. The patient-based estimate is approximately 0.418 across 7,969 eligible genes from 27 paired patient representatives. These analyses use different feature eligibility after patient selection; their agreement supports descriptive stability, not an independent validation of measurement accuracy.

Fifteen of 25 evaluated canonical markers retain the expected group ordering. The joint label permutation remains favourable, and the patient-based, centre-restricted permutation provides the more relevant design-aware check. Centre and histotype remain partly confounded, so these tests do not estimate an unobserved histotype-by-centre combination.

The four bridge samples support abundance-dependent estimates of cross-plex agreement. Their repeatability CV requires equal, independent measurement-error assumptions and does not establish a calibration curve for technical noise across every model. Differences in RNA and protein dispersion can arise from biology and assay properties; the dataset does not identify their separate contributions.

The exploratory HGSC DNA-repair GO result is sensitive to the analysis unit: its best targeted program has adjusted P = 0.0457 across all 15 models and 0.0646 after selecting the 12 patient representatives. It should remain suggestive at patient level. These ancillary gene-set grades are not main manuscript validation claims.

The supplied protein values are described as log2 abundances normalised using a pooled internal standard. Their feature-specific baselines do not support calling them simple log2 ratios to that standard without the generating laboratory's transformation formula.

The recovered CNVkit target footprint permits an exome-opportunity sensitivity analysis of mutational signatures. For TOV21G, the fitted MMR-associated fraction changes from 0.733 with the unadjusted full reference to 0.447 with the target-adjusted full reference (0.475 with the restricted target-adjusted set). These exploratory tumour-only fits support follow-up; they do not establish MSI status, clinical eligibility, or a precise biological exposure. The derived target footprint is not a sample-specific callable mask.

## Manuscript and author input

The v6 manuscript incorporates the supplied RNA library/sequencing methods, correctly separates newly generated RNA/proteomics from reprocessed WES, uses clearer patient terminology and Oxford commas, and removes draft-process commentary from the scientific narrative. Citation metadata were checked; author-and-DOI citations support the author's intended reference-manager workflow.

The author confirmed that the external normal exomes used the same capture kit; earlier prose stating unconfirmed capture compatibility has been corrected. The author also clarified that contributing laboratories may have performed STR and mycoplasma testing. Pending documentation is therefore treated as a confirmation task, never as evidence that testing was not performed.

Published culture methods and a relevant Negri/Morin proteomics protocol provide concrete drafts for collaborator review. They are recorded with their sources and unresolved choices, including differing oxygen conditions across cell-line derivation papers. None is silently asserted to be the protocol used for these samples.

## Reproducibility and release

The processed-data release supplies a cohort-only model table, principal matrices, feature/QC annotations, exact reference records, a data dictionary, a column inventory, source catalogue, checksums, and a rebuilt local viewer. The builder checks assay-model sets, patient counts, protein row alignment, zero-plex flags, and viewer source hashes. The canonical analysis order remains in `scripts/run_all.sh`; the new release builder runs after the scientific analyses.

The four primary DepMap input files were independently hashed and match the sizes and MD5 values in the official [DepMap 24Q4 Public version-1 record](https://doi.org/10.25452/figshare.plus.27993248.v1). The external-reference manifest pins these and the cached Cellosaurus responses. The latter preserve individual entry versions; the original response headers and database release are unavailable. The July 23, 2026 snapshot date comes from existing project documentation, not filesystem timestamps.

Earlier reports and the original v5 files remain preserved. They are historical records, not compatible numerical companions to the repaired outputs. Before-change snapshots and detailed logs are stored in this audit directory. The current entry point is `reports/README.md`.

## Verification and review records

- All 31 libraries passed transcript-reference reconciliation and abundance-conservation checks. Before/after and primary-locus sensitivity estimates are recorded in the [RNA/proteomics audit](rna_proteomics_audit.md).
- Independent checks reproduced the original 6,036 retained candidate identities, validated the recovered conversion, and recomputed all 897 model/arm calls. See the [genomics audit](genomics_audit.md) and `genomics_validation.json`.
- The processed package passed its cohort, assay, protein-row, missingness, external-reference, and viewer checks. All 49 files listed in `release/SHA256SUMS` were independently checked; `release_validation.json` records the result.
- Every viewer RNA value (1,112,342) and measured protein value (238,259) agrees with the processed CSVs within the declared display precision; 22,978 missing protein values remain missing. Browser checks covered restored FOXL2, multiple gene-ID links, and scale controls. The preview server was stopped; the hosted site was not changed.
- Main and supplementary figure exports were visually reviewed. The main figures use Arial with the Cook Lab palettes. The [manuscript review](manuscript_review.md) records page-by-page DOCX verification, and the [comment register](manuscript_comment_responses.md) addresses all 34 original Word comments.

These are checks of the repaired processed workflow and its reported claims. Raw-read alignment/calling and raw-spectrum reprocessing were outside this audit; their upstream experimental records still require completion.

## Remaining author actions

Update after the audit: the author has access to the original WES cluster analysis. The [cluster retrieval checklist](../../docs/manuscript/v6/wes_cluster_retrieval.md) prioritizes run scripts/configuration, sequencing and coverage QC, capture/reference records, and an inventory of final alignments. Historical cluster path hints have been extracted from the archived commands and VCF headers. Cluster access may resolve the local WES gaps; file availability and any additional analyses have not yet been verified.

1. Confirm source-stock STR and mycoplasma records with contributing laboratories, including test dates and correspondence to the analysed stocks.
2. Confirm culture conditions, passages, and the source-specific derivation references; review the drafted proteomics method with the generating laboratory.
3. Recover the exact Sarek release, capture-kit identity and original vendor design, and relevant read-level quality summaries from the cluster records. GATK 4.5.0.0 and VEP 113.0 were recovered from all 23 annotated VCF headers. The CNVkit-derived target footprint (290,475 bins; 63,465,385 union bases) has also been recovered. Same-kit compatibility is already author-confirmed.
4. Complete ethics/consent provenance, authorship, contributions, competing interests, funding, and acknowledgements.
5. Deposit raw and processed records with licences, accessions, and anonymous reviewer access; archive the code release and update the hosted viewer to the reviewed local build.

The detailed comment-response register and collaborator confirmation document accompany the revised manuscript. These actions require experimental or administrative information and are not replaced by computational inference.
