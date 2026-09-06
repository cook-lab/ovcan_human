# Current status and remaining work

The current manuscript is **v9**. The September 2026 audit and figure redesign are preserved. Both recovered WES handoffs have been independently checked and integrated into the methods, technical validation and dependent figures. Target-only CNV resegmentation corrects the antitarget artifact exposed by the second handoff. The public repository [cook-lab/ovcan_human](https://github.com/cook-lab/ovcan_human) remains public by author choice.

## Available now

- 42 models from 34 patients; RNA/protein each cover 31 models, 30 are paired, WES covers 23 models from 16 patients, and 13 have all three assays.
- All 23 variant and CNV profiles, with the same 6,194 coding candidates. CNV profiles now use target-only resegmentation after the second handoff exposed overlapping antitarget bins with extreme ratios; see the dated coverage update for the measured impact.
- Sarek 3.5.1 / Nextflow 24.10.2 successful-run evidence; acquisition records identify HiSeq 4000 PE100 and SeqCap EZ Exome v3. Twenty-two provider rows match names/passages; the TOV81D pooled alias needs confirmation.
- Independently parsed QC for all 23 WES models: read depth, alignment, duplicates, capture coverage and contamination, with denominators and source hashes.
- Independently verified coverage for all five CNV-reference exomes: mean target depths 75.85–95.72× (median 92.40×). All 204,706 reference-mask-passing targets have positive mean depth in at least three normals.
- Reconciled vendor/lifted/derived target definitions. The 63,709,951 bp sum used by mosdepth differs from the 63,514,049 bp interval union and the 63,465,385 bp standard-chromosome subset used for signature sensitivity.
- Corrected interpretation of MAF PASS (includes converter-added population flags) and of target-supported CNV segment spans (not direct continuous genomic coverage).
- Manuscript v9, 14 figure exports, the prior validated 49-file release, local browser build and source scripts. The browser has not been redeployed by this update; new QC companion records have not yet been added to the deposition package.

Read [the latest coverage and CNV correction](../reports/wes_cnv_coverage_2026-09-06/COVERAGE_UPDATE.md), [the initial WES completion report](../reports/wes_completion_2026-09-05/WES_COMPLETION.md), [manuscript readiness review](../reports/wes_completion_2026-09-05/manuscript_readiness.md) and [author confirmations](manuscript/v9/author_confirmation.md).

## Remaining paper requirements, before deposition

| Owner | Required records or decision |
| --- | --- |
| WES laboratory/provider | DNA extraction and detailed library preparation; confirm `TPV81D_23-pool` corresponds to TOV81D and document pooling |
| WES cluster | Manual execution/input-alignment evidence; targeted workflow/reference/liftover/converter provenance gaps where recoverable; additional existing normal alignment QC if available. All five normal coverage pairs are received and verified. |
| Contributing laboratories | Culture and harvest conditions and relation between assay aliquots (M01); RNA extraction/integrity (M02); actual proteomics acquisition/search/scaling/CV definition (M03) |
| Contributing laboratories | Stock-linked STR and mycoplasma records (M09) |
| Authors/institutions | Ethics/consent/data-sharing determination (M07); authors, affiliations, contributions, interests and funding (A01, A03–A06) |
| Authors | Final code licence/release identifier and choice of additional analysis figures |

The actionable cluster request is now [the targeted follow-up](cluster/recovery/2026-09-06/FOLLOWUP.md). The original broad search checklist is historical context, not a list of still-missing results. Missing from this tar is not proof of absence from the cluster.

The author's 13 current Word comments are addressed in v9; see [comment responses](manuscript/v9/user_comment_responses.md), including the PC1-decomposition explanation. A02 is resolved by confirmation that all three datasets are new. M06 is resolved by the analysis-team Cellosaurus provenance audit with explicit limits. M08 wording is approved. Culture and proteomics protocol drafts are directly in highlighted Methods for collaborator confirmation. The source v7 Word file remains the author's working copy.

## Exploratory molecular annotation

The author subsequently requested an exploration of clinical/molecular classifications. The [6 September assessment](../reports/clinical_classification_2026-09-06/README.md) adds separate CCNE1/ERBB2 locus and expression evidence, curated pathogenic BRCA2 significance in TOV81D, and a review of HRD/MSI/receptor annotation options. The [new focused cluster plan](cluster/CLINICAL_CLASSIFICATION_NEXT_STEPS.md) describes input preparation and a proposed allele-specific pilot. No HRD score or clinical positive/negative classification has been assigned, and manuscript v9 and the validated 49-file release are unchanged by this exploration.

The subsequent [expanded analysis](../reports/molecular_extension_2026-09-06/README.md) covers 52 CN loci, 19 expression targets, 65 exact variant candidates and SBS3 reference/filter/bootstrap sensitivity. AKT2/OV3331 is a concordant DNA/RNA/protein lead; NF1/CDKN2A depletion needs exon/read validation. No model has dictionary-robust positive SBS3 support. The current operational cluster handoff is [MEX01–MEX11](cluster/molecular_extension_2026-09-06/AGENT_TASK.md), with model aliases, candidate regions, input/status templates and a portable overlay builder. Discovery and proposed pilots are prepared; no cluster jobs have been run from this workstation. The author-approved [NCBI/ClinVar query](../reports/molecular_extension_2026-09-06/variants/clinvar/README.md) is complete: 23 of 53 alleles match exact records, covering 29 model observations. Updated ATR/PTEN/BRAF annotations and conflicting classifications are retained with their separate evidence scopes.

## Separate deposition and release phase

Sequence and proteomics accessions, reviewer access, controlled-access decisions, immutable data/code archives, licence selection and deployment of the corrected hosted browser remain. No new allele-specific CNV, HRD, MSI or signature-discovery analysis is required to support the current descriptive paper; those would need a separate justified analysis and suitable inputs. The raw handoff tar and extracted provider/run records remain outside Git.
