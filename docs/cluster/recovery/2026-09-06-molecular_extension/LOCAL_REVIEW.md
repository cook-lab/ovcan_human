# Local integration of the cluster return

Reviewed 6 September 2026. Incoming file: `data/ovcan_human_recovery_docs_ab00837.patch`, 2,280,112 bytes; SHA-256 `20566dcea61382ad76157aa24e22ef00e03715121a2dd168c962744815321cc1`. The patch identifies cluster commit `ab00837bb5d63cf7789cab835cc3054c43a95b3e` and contains **23 new files**, including inventories, historical recovery reports and five proposed job scripts. It contains no new MSI, allele-specific CN or HRD results. The cluster reports no submitted analysis jobs.

All 23 files were inspected and their extracted Git blob IDs verified before incorporation. Their bytes are preserved under their original repository paths; [received_files.json](received_files.json) records SHA-256 values. The original patch remains unchanged and is ignored as a transport artifact. No incoming script was executed, no cluster connection was made, and no raw data or manuscript results were changed. Curated records contain research model IDs, public reference accessions, cluster paths and lab usernames. Raw sequencing inputs, private logs and provider workbooks are not included.

## Accepted findings and qualifications

| Finding | Evidence and scope |
| --- | --- |
| All 23 complete md CRAMs, restricted recal generations and recorded CNV-input alignments located | Cluster-reported observations in [input_inventory.tsv](input_inventory.tsv). Readability/index and read-group checks must be repeated at execution. The workstation has not inspected these alignments. Matching names, sizes and mtimes alone do not establish byte-identical generations. |
| CNV source/target-only CNR consistency | All 23 source-CNR and reconstructed-target-CNR hash strings agree with the local canonical manifest. All 23 reported corrected CNS hashes agree with actual local file hashes. This independently checks the returned manifest against the checkout; it does not rehash remote alignments/CNR files. |
| Missing off-target model CNV depth explained | The recovered ApplyBQSR/CNV-input chain identifies interval-restricted model alignments, whereas public normal alignments retained off-target reads. This resolves the upstream cause of near-empty antitargets. The reason four historical antitarget intervals overlap capture targets still needs construction provenance. Canonical script 29 target-only profiles remain unchanged. |
| Five reference normals have different preprocessing | Bwa-mem2, no duplicate marking or BQSR in recovered scripts. Same kit is author-confirmed; process matching is not. The incoming normal-pilot phrase “process-matched approximation” is superseded. |
| No eligible germline-inclusive SNP, MSI or allele-specific results recovered | Within the searched WES roots. New pilots are justified; negative findings are not a claim that files never existed elsewhere. RNA/proteomics locations remain separate retrieval questions. |
| Reference hashes recovered | Five full reference/resource hashes, including the FASTA, are reported in [the follow-up manifest](../2026-09-06/followup_inventory_reference_hashes.tsv). These were computed on the cluster and not recomputed from multi-GB resources here. The beta known-indels resource is still unlocated in the searched archive. |
| Sarek source evidence strengthened | The reported Nextflow 24.10.2 local-script revision agrees with the released `main.nf`; the full recorded revision is `3954909713023f4328e976337e6e2cb9`. This supports top-level file consistency. It does not certify every workflow/module/config file as unmodified or establish a complete pipeline Git commit. Original source/image directories remain inaccessible. |
| Manual execution and liftover gaps narrowed | Input chain and indirect execution evidence were found. Exact manual-CNVkit invocation logs, original 22-model vcf2maf execution records, and original liftover chain/command remain absent from searched locations. The available liftOver v483 binary is not proof it executed the historical transform. Avoid repeating broad searches of the same roots. |

The inventory contains **204 rows: 202 labelled `verified`, one RNA `not_found` in WES roots and one proteomics `not_searched`**. The `verified` label needs context: the CNVkit container row explicitly says the pull was in progress. Check completion and image hash before use. The old claim of “204 verified rows” should not be repeated.

## What changed in the proposed analysis

The current [execution plan](../../molecular_extension_2026-09-06/EXECUTION_PLAN.md) and [corrected recipes](../../molecular_extension_2026-09-06/RECIPES.md) supersede the proposed scripts and technical questions in this source return. Main corrections:

- Complete md CRAMs are the default for SNP input as well as MSI/read review, with a fixed-SNP comparison to the recorded recal generation.
- PureCN receives CNR through `--tumor` and the full filtered VCF, without the incompatible Mutect2 stats argument or unrelated calling-interval file. High-purity homozygous-site retention and explicit base-quality/padding settings are included.
- Purity comparisons hold ploidy/max-copy settings fixed. The small differently processed normal panel is conditional sensitivity; beta-binomial calibration is not assumed. BAF-aware segmentation and competing-fit review precede LOH/scars.
- Read review expands from 25 to 30 records to include the five strong CDK12 loss candidates omitted from the original BED. Per-exon work preserves zeros, isoforms, interval definitions and processing differences.
- Input manifests, runtime/failure records, shared-file race prevention and checksum self-exclusion are required before submission.

`python3 scripts/51_prepare_cluster_execution_inputs.py` validates source hashes, model/passages, all CNR/CNS manifest relationships, cohort denominators and the focused request set, then writes the portable execution inputs. See [execution_input_validation.json](../../molecular_extension_2026-09-06/execution_input_validation.json).

## Historical documents and manuscript precedence

The added `2026-09-05/` files are the original cluster discovery record, not a replacement manuscript or analysis. Their statements about unrecorded instrument/capture details, archived CNS use, source versions and other open items predate local acquisition reconciliation and target-only resegmentation. The current `docs/PROJECT_STATUS.md`, manuscript v9 and later verified reports take precedence. Additional legacy exomes remain outside the 23-model analysis and are not silently added to the paper.

This return strengthens the future WES Methods/provenance update (interval-restricted historical CNV inputs, available complete alignments, reference hashes and remaining execution limits). It does not alter the 6,194 candidate calls, 23-model cohort, figures, v9 DOCX or validated 49-file release. Integrate new molecular results into a later manuscript revision only after accepted pilot outputs are returned.
