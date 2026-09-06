# Remaining WES follow-up after the coverage archive

The first WES handoff and the 6 September CNV coverage archive have both been inspected. Read [the coverage update](../../../../reports/wes_cnv_coverage_2026-09-06/COVERAGE_UPDATE.md) for the new checks and [the initial completion report](../../../../reports/wes_completion_2026-09-05/WES_COMPLETION.md) for the broader context. This file supersedes the 5 September follow-up.

## Received: do not request these again

- Target and antitarget coverage CNNs for all 23 models and all five reference normals: SRR4039087, SRR4039088, SRR4039089, SRR4039096 and SRR4039097.
- Pooled `reference.cnn`, derived target BED and all three command files plus `run_cnvkit.sh`.
- Forty-six historical CNV plots. There is no README/handoff narrative in the new archive.
- Existing Sarek model sequencing QC, variants, CNV segment records, vendor capture design, lifted target intervals and acquisition workbooks from the first handoff.

All 108 payload files match the new archive's checksum manifest. Its self-entry for `SHA256SUMS` does not match the completed manifest; this is recorded separately and does not invalidate the verified payload hashes. Future manifests should exclude their own output file. A replacement archive is not needed solely for this packaging defect.

The normal target means are now verified. Reference-mask-passing CNV targets have coverage in at least three normals. The model/reference antitarget discrepancy is directly documented; the four positive model antitargets overlap captured targets and do not supply independent off-target evidence. Their extreme ratios inflated two archived chromosome 1 segments (OV2295 and TOV2835EP). The local correction removes antitarget rows and resegments the unchanged target ratios; its inputs, settings and comparison are in `output/wes_cnv_target_only/` and the coverage update. This correction does not establish the upstream cause of the coverage discrepancy.

## Outstanding records to seek

| Priority | Existing records | Purpose |
| --- | --- | --- |
| High | Successful manual CNVkit Slurm/shell logs, `.command*` records if applicable, input alignment paths/links and headers, exact tool versions | The four supplied command/wrapper files are identical to the first handoff; they still do not prove which invocation completed or which alignment generation produced each output. Preserve sample/passage mappings and distinguish model and normal preprocessing, especially whether input model alignments were restricted to capture intervals while normal alignments retained off-target reads. This is a hypothesis to check in records, not an established cause. |
| Medium | Original 22-model VCF-to-MAF command/version or job log | The frequency-flag behaviour has been reconstructed, but the original converter version remains unknown. TOV3121D's local vcf2maf 1.6.22 reconstruction is already documented. |
| Medium | Successful Sarek local source archive/revision and local modifications; existing FASTA/resource checksum manifests; container identities/digests | The release and output-to-task links are verified. A 32-character Nextflow script hash is not an established Git commit, and sidecar hashes are not the full FASTA checksum. |
| Medium | Target-liftover executable/version, exact command and chain-file identity/checksum | Input/output intervals and unmapped targets reconcile, but the original transform invocation is still unavailable. |
| Medium | CNV antitarget-construction command/log and target/accessibility BED identities used at that step | Explain the four antitarget intervals that overlap the final capture targets, and whether the target/antitarget definitions were built from different intermediate interval sets. The interval files themselves are already present. |
| If available | Reference-normal alignment/duplicate QC logs and metadata | Mean bin-depth files do not establish mapping, duplication, insert-size or per-base threshold coverage. Record existing values with their actual denominators; absence of those reports is not a reason to invent or rerun QC automatically. |
| Provider/laboratory | DNA extraction and library-preparation SOP; confirmation of `TPV81D_23-pool` | Confirm the provider alias, constituent preparations/stock and passage. The other 22 acquisition rows match their model/passage names; instrument, read configuration and capture identity are resolved. |

Return a small dated inventory with source paths, file roles, sample/run identities and available hashes, classifying records as found, partially found, not found in searched locations, or inaccessible. Keep raw run logs private and curate nonsensitive evidence for the public repository. Historical filenames and matching parameters alone do not prove byte-identical alignment inputs.

This is targeted retrieval of existing records. Do not submit alignment/calling/HRD/MSI jobs, change permissions, delete files, transfer large sequencing inputs or deposit additional data under this request. Sequence-input selection and data deposition remain a separate phase. If the requested historical evidence is unrecoverable, document that limit explicitly rather than repeating a broad search indefinitely.
