# Remaining WES follow-up after the coverage archive

The first WES handoff and the 6 September CNV coverage archive have both been inspected. Read [the coverage update](../../../../reports/wes_cnv_coverage_2026-09-06/COVERAGE_UPDATE.md) for the new checks and [the initial completion report](../../../../reports/wes_completion_2026-09-05/WES_COMPLETION.md) for the broader context. This file supersedes the 5 September follow-up. The later [cluster return and local review](../2026-09-06-molecular_extension/LOCAL_REVIEW.md) have now resolved the input-generation question and supplied full reference/resource hashes. New computation is specified separately in the [execution plan](../../molecular_extension_2026-09-06/EXECUTION_PLAN.md).

## Received: do not request these again

- Target and antitarget coverage CNNs for all 23 models and all five reference normals: SRR4039087, SRR4039088, SRR4039089, SRR4039096 and SRR4039097.
- Pooled `reference.cnn`, derived target BED and all three command files plus `run_cnvkit.sh`.
- Forty-six historical CNV plots. There is no README/handoff narrative in the new archive.
- Existing Sarek model sequencing QC, variants, CNV segment records, vendor capture design, lifted target intervals and acquisition workbooks from the first handoff.

All 108 payload files match the new archive's checksum manifest. Its self-entry for `SHA256SUMS` does not match the completed manifest; this is recorded separately and does not invalidate the verified payload hashes. Future manifests should exclude their own output file. A replacement archive is not needed solely for this packaging defect.

The normal target means are now verified. Reference-mask-passing CNV targets have coverage in at least three normals. The model/reference antitarget discrepancy is directly documented; the four positive model antitargets overlap captured targets and do not supply independent off-target evidence. Their extreme ratios inflated two archived chromosome 1 segments (OV2295 and TOV2835EP). The local correction removes antitarget rows and resegments the unchanged target ratios; its inputs, settings and comparison are in `output/wes_cnv_target_only/` and the coverage update. The later input-chain recovery establishes the upstream cause of near-empty model antitarget coverage: CNVkit received interval-restricted recalibrated alignments, while the public normals retained off-target reads. Complete duplicate-marked model CRAMs remain available. The separate cause of four antitarget intervals overlapping capture targets still needs construction provenance.

## Outstanding records to seek

| Priority | Existing records | Purpose |
| --- | --- | --- |
| If newly located | Successful manual CNVkit invocation logs or source image digests | Input paths and the Oct-24 recal-to-BAM chain are received for all 23 models; commands/output timestamps supply indirect execution evidence. Exact invocation logs were not found in the searched project/archive and old scheduler records have expired. Do not repeat those same searches; a separately supplied original location may help. |
| Medium, workstation records | Original 22-model VCF-to-MAF command/version | Not found in searched cluster roots. TOV3121D's local vcf2maf 1.6.22 reconstruction is already documented. |
| If recoverable | Remaining workflow modules/configuration and original container digests | The reported local-script revision agrees with released Sarek `main.nf`; this does not certify the inaccessible full source tree or its other files. Container tags were recovered; most historical image digests remain unavailable. |
| Targeted reference gap | `Homo_sapiens_assembly38.known_indels.vcf.gz` under the historical GATKBundle/beta path | Full FASTA, dbSNP146, Mills/1000G, variant PoN and AF-only gnomAD hashes are received. The beta known-indels resource is absent from the searched tar; do not request the five recovered hashes again. |
| Medium | Actual target-liftover command and chain identity/checksum | Input/output intervals reconcile. A liftOver v483 binary was found, but the executed transform may have used the web service. The binary alone does not identify the original command or chain. |
| Medium | CNV antitarget-construction command/log and exact target/accessibility BEDs used | Explain why four antitarget intervals overlap final capture targets. The near-zero model off-target coverage is already explained by interval-restricted input alignments. |
| If available | Existing normal alignment/duplicate QC | Bwa-mem2 without duplicate marking/BQSR is now documented; all five CNN coverage pairs are received. Additional mapping/insert-size/duplicate metrics remain useful only if already present or generated under the separate analysis plan. |
| Provider/laboratory | DNA extraction/library-preparation SOP and `TPV81D_23-pool` confirmation | Confirm the pooled preparation/stock/passage. Other acquisition rows, instrument, read configuration and capture identity are already reconciled locally. |

Return a small dated inventory with source paths, file roles, sample/run identities and available hashes, classifying records as found, partially found, not found in searched locations, or inaccessible. Keep raw run logs private and curate nonsensitive evidence for the public repository. Historical filenames and matching parameters alone do not prove byte-identical alignment inputs.

This is targeted retrieval of existing records. Do not submit alignment/calling/HRD/MSI jobs, change permissions, delete files, transfer large sequencing inputs or deposit additional data under this request. Sequence-input selection and data deposition remain a separate phase. If the requested historical evidence is unrecoverable, document that limit explicitly rather than repeating a broad search indefinitely.
