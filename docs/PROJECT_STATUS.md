# Current status and remaining work

The current manuscript is v7. The September 2026 audit corrected the RNA reference, recovered the 23rd variant profile, revised genomic interpretation and patient-level summaries, and rebuilt the data browser and processed release locally. All 14 figures were subsequently redesigned and checked. The private working repository [cook-lab/ovcan_human](https://github.com/cook-lab/ovcan_human) supports collaboration and cluster recovery; public deposition and journal submission remain separate steps.

## What is already available

- 42 models from 34 patients; 31 RNA profiles, 31 proteomic profiles, 30 paired models, 23 WES models from 16 patients and 13 models with all three modalities.
- All 23 processed variant and copy-number profiles, including the recovered TOV3121D MAF; 6,194 retained coding candidates.
- GATK 4.5.0.0, VEP 113.0 and CNVkit 0.9.10 evidence; 290,475 derived CNVkit target intervals; a 63,465,385 bp merged target footprint; five public CNV-reference accessions from PRJNA339046.
- The matched Ensembl release-93 map and RNA reconciliation records; verified DepMap 24Q4 Public version-1 metadata and current comparison results.
- The validated processed release, manuscript v7, clean figure exports, external legends and source scripts.

The [audit](../reports/audit_2026-09-05/AUDIT_REPORT.md), [WES input manifest](../output/wes_input_manifest.csv) and [data inventory](data/README.md) distinguish these recovered/local records from genuinely missing evidence. Missing from Git and missing from the study are different states.

## Priority retrievals from the WES cluster

| Need | Why it matters | Detailed task |
| --- | --- | --- |
| Sarek/Nextflow run version, revision, launch command, samplesheet, resolved parameters, configs and containers | Link the retained results to the actual workflow and references; do not infer the release from a tool version | [WES recovery](cluster/WES_RECOVERY.md) |
| Existing QC for 23 models and five reference exomes, including exported tables | Complete WES technical validation: depth, target coverage, duplicates, alignment and contamination with correct denominators | [WES recovery](cluster/WES_RECOVERY.md) |
| Vendor capture design, kit revision and preparation of calling/target/antitarget intervals | Distinguish the original design from derived CNVkit bins; same-kit compatibility is already author-confirmed | [WES recovery](cluster/WES_RECOVERY.md) |
| Custom filtering, VCF-to-MAF and CNV-reference scripts; logs and command records | Establish how outputs were produced and why they differ across runs | [WES recovery](cluster/WES_RECOVERY.md) |
| Original reference resources and checksums, reference-normal processing and variant panel of normals | Separate the two normal resources and make the methods reproducible | [WES recovery](cluster/WES_RECOVERY.md) |
| BAM/CRAM/FASTQ locations, indexes, sample/read-group mapping and existing checksums | Support any later QC calculation, targeted inspection or deposition without transferring large files now | [Cluster task prompt](cluster/CLAUDE_TASK.md) |

Retrieval should first recover existing evidence. New allele-specific CNV, HRD or MSI work is optional and needs a separate assessment of inputs and suitability. The present data do not establish those endpoints.

## Author and laboratory records still needed

These may not be on the WES cluster. Keep their ownership separate from a filesystem search.

| Identifier | Outstanding information |
| --- | --- |
| M01 | Actual culture and harvest conditions, model derivation and relationship between assay aliquots |
| M02 | RNA extraction, integrity measurements/criteria and pre-library checks; original quantification index/build provenance if available |
| M03 | Actual proteomics acquisition, search, filtering, channel/normalisation workflow and meaning of the supplied CV |
| M04–M05 | WES run/acquisition/capture/reference/QC records listed above |
| M06 | Original Cellosaurus release/retrieval documentation, if retained; DepMap provenance is resolved |
| M09 | Stock-specific STR and mycoplasma records tied to the profiled material; absent documentation does not mean tests were not performed |
| M07–M08 | Ethics/consent/data-sharing determination and author confirmation of the AI-assistance statement |
| A01–A06 | Authors, affiliations, prior uses of these exact datasets, contributions, interests, acknowledgements and funding |
| D01 | Data accessions, reviewer access, licences, access conditions and immutable deposited checksums |
| D02 | Verify/deploy the corrected hosted browser and archive that version; local rebuild alone is not deployment |
| D03 | Private working repository created at [cook-lab/ovcan_human](https://github.com/cook-lab/ovcan_human); public code availability and an immutable DOI archive remain pending |

See [author confirmations](manuscript/v7/author_confirmation.md) for the full questions and known evidence. Do not fill a manuscript TODO from a related publication unless the authors confirm that it describes this dataset.

## Recovery output expected

Provide a dated report linking each request to actual files, their generating run and verified model/sample aliases. Label items `found`, `partially found`, `not found in searched locations`, `inaccessible` or `not produced` only when the evidence supports that distinction. Include search roots and limitations. Keep checksums and original paths in a manifest, retain small scripts/configs with provenance, and list proposed follow-up calculations separately from recovered results.
