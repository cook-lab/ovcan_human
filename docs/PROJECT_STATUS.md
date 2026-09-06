# Current status and remaining work

The current manuscript is **v8**. The September 2026 audit and figure redesign are preserved. The recovered WES handoff has now been independently checked and integrated into the methods, technical validation, Figure 2 and Figure 4. The public repository [cook-lab/ovcan_human](https://github.com/cook-lab/ovcan_human) remains public by author choice.

## Available now

- 42 models from 34 patients; RNA/protein each cover 31 models, 30 are paired, WES covers 23 models from 16 patients, and 13 have all three assays.
- All 23 variant and CNV profiles, with the same 6,194 coding candidates and patient-level summaries as the audited release.
- Sarek 3.5.1 / Nextflow 24.10.2 successful-run evidence; acquisition records identify HiSeq 4000 PE100 and SeqCap EZ Exome v3. Twenty-two provider rows match names/passages; the TOV81D pooled alias needs confirmation.
- Independently parsed QC for all 23 WES models: read depth, alignment, duplicates, capture coverage and contamination, with denominators and source hashes.
- Reconciled vendor/lifted/derived target definitions. The 63,709,951 bp sum used by mosdepth differs from the 63,514,049 bp interval union and the 63,465,385 bp standard-chromosome subset used for signature sensitivity.
- Corrected interpretation of MAF PASS (includes converter-added population flags) and of target-supported CNV segment spans (not direct continuous genomic coverage).
- Manuscript v8, 14 figure exports, the prior validated 49-file release, local browser build and source scripts. The browser has not been redeployed by this update; new QC companion records have not yet been added to the deposition package.

Read [the WES completion report](../reports/wes_completion_2026-09-05/WES_COMPLETION.md), [manuscript readiness review](../reports/wes_completion_2026-09-05/manuscript_readiness.md) and [author confirmations](manuscript/v8/author_confirmation.md).

## Remaining paper requirements, before deposition

| Owner | Required records or decision |
| --- | --- |
| WES laboratory/provider | DNA extraction and detailed library preparation; confirm `TPV81D_23-pool` corresponds to TOV81D and document pooling |
| WES cluster | Five-normal coverage/QC files, manual execution evidence; targeted workflow/reference/liftover/converter provenance gaps where recoverable |
| Contributing laboratories | Culture and harvest conditions and relation between assay aliquots (M01); RNA extraction/integrity (M02); actual proteomics acquisition/search/scaling/CV definition (M03) |
| Contributing laboratories | Stock-linked STR and mycoplasma records (M09) |
| Authors/institutions | Ethics/consent/data-sharing determination (M07); AI-use verification (M08); authors, affiliations, prior data uses, contributions, interests and funding (A01–A06) |
| Authors | Final code licence/release identifier and choice of additional analysis figures; original Cellosaurus retrieval documentation if retained |

The actionable cluster request is now [the targeted follow-up](cluster/recovery/2026-09-05/FOLLOWUP.md). The original broad search checklist is historical context, not a list of still-missing results. Missing from this tar is not proof of absence from the cluster.

## Separate deposition and release phase

Sequence and proteomics accessions, reviewer access, controlled-access decisions, immutable data/code archives, licence selection and deployment of the corrected hosted browser remain. No new allele-specific CNV, HRD, MSI or signature-discovery analysis is required to support the current descriptive paper; those would need a separate justified analysis and suitable inputs. The raw handoff tar and extracted provider/run records remain outside Git.
