# Targeted WES follow-up after the first handoff

The 5 September 2026 handoff was received and independently compared with the retained study results. Read [the completion report](../../../../reports/wes_completion_2026-09-05/WES_COMPLETION.md) before searching. The original broad recovery checklist is now historical context. Do not repeat the completed 23-model QC extraction or merge the 22 additional legacy exomes into this cohort.

## Received and resolved

- Sarek 3.5.1 / Nextflow 24.10.2 run records, including the completed 23 October 2025 run and its executed variant-processing tasks.
- Model-level fastp, Picard, samtools, mosdepth and GATK contamination records for all 23 exomes. Picard files are nested under `reports/markduplicates/<sample>/`; the handoff generator missed them, but the local validator reads them.
- Provider acquisition workbooks, original SeqCap EZ Exome v3 target/probe definitions, converted GRCh38 intervals and CNV target bins.
- Retained variant identities/filtering reconciliation and 69 matching manual CNV segment files. The MAF PASS count already includes converter population-frequency flags.
- Evidence that current recalibrated alignments are interval restricted and later launches wrote to the same result tree. They must not be described as byte-identical retained-VCF inputs without a separate match.

## Files or records still worth retrieving

| Priority | Find | Why / acceptance evidence |
| --- | --- | --- |
| High | All five normal `*.targetcoverage.cnn` / `*.antitargetcoverage.cnn`, existing reference/QC plots or metrics, model coverage files if available | The pooled `reference.cnn` already exists in the local archive (35 identical copies; SHA-256 `848c052274264ac544897977648860455ac742d89bc1de5edc7ece1ca185eeda`). The bundle has no per-normal coverage `.cnn` files. Normal depths in the recovery summary are reported values, not independently verified primary evidence. Preserve exact paths, checksums and model/normal aliases. |
| High | Executed manual CNVkit job/script, `.command*` or shell/scheduler logs, tool versions, input alignment links | Establish which alignment generation was used and how the empirical reference, target bins, low-coverage dropping and segmentation were executed. The segment outputs match, but saved command text alone does not prove execution. Almost all model antitarget bins have zero depth; inspect the existing normal/model coverage tables for processing comparability. Do not treat missing off-target depth as a biological deletion. |
| Medium | Original 22-model VCF-to-MAF conversion script, command/version or job log | The output behaviour matches `common_variant` when any of seven gnomAD-exome subgroup frequencies exceeds 0.0004. This does not establish the original software version. TOV3121D's local recovery with vendored vcf2maf 1.6.22 is already documented. |
| Medium | The successful Sarek local workflow source directory/archive; exact reference FASTA/resource checksum manifests; container image identities/digests if retained | The 32-character Nextflow revision field is not a verified Git commit. The dictionary and index identify the reference, but do not substitute for the FASTA file hash. Use existing manifests first; do not checksum multi-terabyte trees. |
| Medium | Target-liftover command, tool version and chain file/checksum | Input and output BEDs reconcile, including 17 unmapped targets and split mappings. The exact transform invocation is still missing. |
| Author/provider | DNA extraction and library-preparation SOP; explanation of `TPV81D_23-pool` | The 2017 sheet has `TOV81D_P23` and `TOV81D_P23_-2` sample records and one differently named pooled sequencing record. Confirm its mapping and what was pooled. Similar provider/pipeline read counts alone are not proof. The other 22 sequencing rows match names/passages. |

Library record clues are `mcgill_r004741_provenance/WES 2017 Mes-Masson.xlsx` and `WES 2018 Mes-Masson.xlsx`, `HiSeqRead Tab`. The 2014 workbook describes a different acquisition and must not be substituted. All 23 provider sequencing rows specify HiSeq 4000 PE100; total-cycle fields are not read lengths.

## Keep deposition work separate

Original provider BAM/CRAM/FASTQ availability, exact final sequence inputs, CRAM references, permissions and repository packaging will matter for deposition. Inventory locations and existing hashes now if useful; do not transfer or publish large alignment/read files. This follow-up does not authorise new alignments, calling, HRD/MSI jobs, deletion, permission changes or new data deposition. Read-only targeted discovery and small evidence retrieval are sufficient.

Return a small dated inventory with `found`, `partially found`, `not found in searched locations`, or `inaccessible`, source paths and hashes. Do not convert absence from this tar into absence from the cluster. Keep copied raw logs private; curate nonsensitive evidence for this public repository.
