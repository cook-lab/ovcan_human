# WES recovery from the original cluster project

This is the working handoff for the cluster instance of Claude Code. It consolidates the September 2026 audit and the manuscript-v7 retrieval list. The task is to **recover existing provenance, scripts and QC, and inventory existing sequencing inputs**. It is not a request to rerun the pipeline or initiate HRD/MSI analysis. The author has said the original cluster analysis is accessible; its paths and contents have not yet been verified from this checkout.

Use [CLAUDE_TASK.md](CLAUDE_TASK.md) as the copy-pastable task prompt. Manuscript items M04/M05 in [author_confirmation.md](../manuscript/v7/author_confirmation.md) describe the publication context. The detailed scientific corrections are in [genomics_audit.md](../../reports/audit_2026-09-05/genomics_audit.md). Earlier inventory prose is historical: statements about only 22 mutation profiles, incompatible capture kits, proven chrX-reference artefacts, or a guaranteed BAM-to-HRD solution have been superseded.

## What the Git checkout contains

The public repository includes code, documentation, current processed results and small source-evidence files. It is not a complete copy of the analysis workstation or cluster. An absent `judy_archive/data/`, BAM/CRAM, FASTQ or large VCF in Git does not mean it is absent on the cluster. The [archived-input inventory](../data/archived_input_inventory.tsv) lists the 906 excluded source-data files and sizes; its paths are relative to the original `judy_archive/data/` directory. Original paths recorded in manifests can refer to files intentionally excluded from Git. See [REPRODUCIBILITY.md](../REPRODUCIBILITY.md) for the checkout and external-input contract. Do not run the full R pipeline merely to test a fresh checkout: it depends on those external inputs and writes derived outputs.

Read these maintained inputs first:

| Repository path | Use |
|---|---|
| [metadata/samples.csv](../../metadata/samples.csv) | Canonical model names, assay availability and assay-specific passages. Filter the generated models with `wes_mut == Y` / `wes_cnv == Y`; do not use older summary spreadsheets as the WES inventory. |
| [metadata/line_family_map.csv](../../metadata/line_family_map.csv) | Patient groups and related models; model count is not patient count. |
| [wes_cluster_models.csv](../../reports/audit_2026-09-05/wes_cluster_models.csv) | Checklist of 23 WES models, passages and historical CNV alignment tokens. |
| [wes_cluster_path_hints.csv](../../reports/audit_2026-09-05/wes_cluster_path_hints.csv) | 322 literal path hints, their roles and original local evidence sources. These are not verified cluster locations. |
| [wes_input_manifest.csv](../../output/wes_input_manifest.csv) | The 23 archived/recovered MAF and source-VCF identities, including SHA-256 hashes for duplicate detection. |
| [wes_vcf_header_provenance.csv](../../output/wes_vcf_header_provenance.csv) and [wes_pipeline_parameters.csv](../../output/wes_pipeline_parameters.csv) | Recovered versions, resource labels and CNV command details. “Not recoverable” in the latter means not recovered from the original local archive, not absent from every cluster location. |
| [wes_filter_cascade.csv](../../output/wes_filter_cascade.csv) and [wes_recovery_validation.json](../../reports/audit_2026-09-05/wes_recovery_validation.json) | Baseline counts and checks on the recovered TOV3121D conversion. |
| [cnvkit_commands.txt](evidence/cnvkit_commands.txt), [cnvkit_command_diagram.txt](evidence/cnvkit_command_diagram.txt), [evidence manifest](evidence/manifest.csv) | Byte-preserved command snapshots with original archive paths and SHA-256 hashes. Preserve the distinction between the original archive path and this repository copy. |
| [wes_cnvkit_target_intervals.bed](../../output/wes_cnvkit_target_intervals.bed) | Recovered derived target bins, not the original vendor capture design. |

The [original WES notebook](<../../judy_archive/notebooks/All Inclusive Final Code WES.Rmd>) is retained as historical code evidence. It predates the corrected scripts and should be inspected for provenance/path clues, not executed as the current analysis.

If a listed supporting file is represented by a packaging placeholder, use its recorded source path/hash; do not silently manufacture an empty replacement. Report any required handoff input that is absent from the checkout.

## Already recovered: do not repeat these as missing results

- **All 23 WES models have both mutation and CNV inputs.** TOV3121D's annotated VCF was present; its MAF conversion was omitted. [wes_maf_inputs.R](../../scripts/lib/wes_maf_inputs.R) recovers it with the vendored vcf2maf 1.6.22 converter and `--inhibit-vep`, preserving the archived annotations. No variant-calling or VEP rerun was needed. Its historical source is `judy_archive/data/wes - old/mutect2/TOV3121D_P68/TOV3121D_P68.mutect2.filtered_VEP.ann.vcf.gz`; the workstation output is `output/recovered_maf/TOV3121D_P68.maf`.
- The audited cascade is **582,474 raw records → 16,081 PASS → 15,995 population-filtered → 6,194 coding candidates**. TOV3121D contributes 25,082 → 389 → 386 → 158, including TP53 p.H214fs. The original 6,036 candidates from the other 22 models were preserved. The WES cohort is 23 models from 16 patients; the HGSC subset is 18 models from 11 patients.
- All primary VCF headers establish **GATK Mutect2 and FilterMutectCalls 4.5.0.0, VEP v113.0, GRCh38.p14**, with VEP annotation labels gnomAD exome/genome v4.1, ClinVar 202404 and dbSNP 156. These labels do not identify the full workflow revision, exact container digest or all reference-file bytes.
- CNVkit **0.9.10** commands and the five public reference exomes **SRR4039087, SRR4039088, SRR4039089, SRR4039096 and SRR4039097** (PRJNA339046) are documented. These are unmatched public exomes used for CNV normalization, not matched normals for the model-specific variant calls. The author confirms use of the **same capture kit** for the models and reference exomes; its name/revision still needs documentary completion.
- The derived CNVkit target BED has **290,475 intervals**, recovered from 35 byte-identical archive copies. SHA-256: `2177970ffdbc933e488067f812e3ae4760310bd9a40fadf4e2fa5389480d955d`. Its standard-chromosome union is 63,465,385 bp and has already supported target-opportunity signature sensitivity. Recovering the original vendor/padded intervals still matters; these derived bins do not identify the kit by themselves.

The raw MAF `NCBI_Build` value sometimes says GRCh37 because of an incorrect conversion default. The audited coordinates are GRCh38, supported by VCF contig lengths, reference-resource paths and known loci. **Do not lift over the data to reconcile that label.** Record reference FASTA identifiers/checksums and treat header correction as a provenance-preserving downstream change, not a new coordinate conversion.

## Prioritized recovery targets

### High: required provenance and technical-validation gaps

| ID | Find and preserve existing evidence | What it resolves / acceptance evidence |
|---|---|---|
| H1 — Run identity and sample mapping | Original launch/submission scripts; scheduler logs; input samplesheets; parameter JSON/YAML; custom Nextflow configs; workflow source revision and modifications; `pipeline_info/`, `.nextflow.log*`, trace/report/timeline and software-version files. Include independently launched CNVkit/normal-reference runs. | The actual workflow identity/release/commit, Nextflow version, invocation, run date, completion state, model-to-library mapping and departures from defaults. Tool versions alone do not establish an nf-core/Sarek release. Tie each retained output to a successful task/run, rather than treating a draft command as execution evidence. |
| H2 — Capture, reference and acquisition provenance | Vendor capture BED/design/interval-list files; kit name/revision; target padding/merging/preparation scripts; exact genome FASTA `.fai`/`.dict`/existing checksum; library/flowcell/read-group manifest; sequencing/provider reports. Also locate Mutect2 germline/variant-PoN references and CNV normal-reference recipes. | Distinguishes original vendor targets, padded calling intervals, interval scatter shards and CNVkit-derived bins. Records reference identity, coordinate conventions, library kit, instrument/read configuration and assay-specific passages. Confirm the documented same-kit statement with available records; do not recast it as an incompatibility finding. |
| H3 — Existing WES QC | `multiqc/` and its `multiqc_data/` tables; sequencing/alignment/duplicate/insert-size/depth/capture reports; available contamination summaries; the underlying numerical files, not HTML alone. Search for FastQC, samtools, Picard, mosdepth and other modules only as filename hints, not assumptions that they ran. | A per-model QC table for all 23 models and a separate five-reference-exome table. Retain units, metric definitions and target denominators. Useful fields include reads/fragments, mapped/properly-paired fractions, duplicate fraction, insert-size summaries, on-target fraction, mean/median target depth, coverage-threshold fractions and contamination estimates. Missing metrics remain missing, never zero or an invented pass. Positive-control driver recovery does not establish exome-wide sensitivity. |
| H4 — Input preservation inventory | Final `*.recal.bam` or `*.recal.cram` and indexes; original FASTQs and existing checksum manifests; normal-exome alignments; symlink destinations. Inventory large inputs rather than copying them. | Whether the profiled libraries can be revisited or deposited. For each model: original path, resolved path, accessible/broken status, size, available checksum, model/passage/read-group evidence, alignment index and CRAM reference. Keep alternate versions separate until matched to the audited outputs. |

### Medium: recover if already present, without launching new analysis

| ID | Existing files to seek | Interpretation / limit |
|---|---|---|
| M1 — Filtering and annotation chain | Per-sample `*.mutect2.contamination.table`, `*.mutect2.segmentation.table`, `*.mutect2.artifactprior.tar.gz`; pileup/orientation-bias summaries; filtering `.stats`; interval merging/filtering scripts; VEP/vcf2maf commands; final and unfiltered VCFs/indexes. Relevant task `.command.sh`, `.command.run`, `.command.log`, `.command.err`, `.exitcode`, `versions.yml`, linked to the run trace. | These exact per-sample filtering filenames are named in the archived VCF headers. Establish actual thresholds and resources; distinguish VEP annotation gnomAD labels from the separate Mutect2 germline-resource VCF. A `mutect2.segmentation.table` is not automatically a major/minor copy-number solution. |
| M2 — CNV reference detail | `reference.cnn`, per-normal/per-model `*.targetcoverage.cnn` and `*.antitargetcoverage.cnn`, target/antitarget generation commands, normal preprocessing and reference sex configuration, segmentation/call options and fit/QC reports. | Documents the empirical pooled reference, normal processing, coverage behavior and relative-copy-number assumptions. Keep the **Mutect2 variant panel of normals** separate from the **five CNV reference exomes**. Do not infer absolute tumour ploidy or a causal explanation for chrX shifts from a relative ratio alone. |
| M3 — Previously completed specialized work | Existing allele-specific major/minor copy-number, purity/ploidy/BAF solutions, MSI analyses or orthogonal validation results, with inputs, commands and fit diagnostics. Search exact available workflow names before broad filenames. | Inventory and assess suitability only. BAM access, BAF fields, integer total-copy calls or a file named “segmentation” do not establish a usable HRD score. Current mutation signatures remain exploratory and MSI unconfirmed. No new HRD/MSI job is part of this handoff. |

Stock STR/mycoplasma records and ethics/capture-provider confirmations may sit outside the compute project. Note a relevant record if found, but do not let those laboratory requests block recovery of H1–H4. Caller contamination estimates are not stock authentication.

## Literal cluster clues and name traps

The following paths are copied from local evidence, **not invented locations or claims of current existence**:

- `/scratch/asmab/ovcan/` — CNV command reference/interval root.
- `/scratch/asmab/ovcan/intervals_sorted.bed` and `/scratch/asmab/ovcan/intervals.antitarget.bed`.
- `/scratch/asmab/ovcan/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta`.
- `/scratch/asmab/nf-core/cnvkit-0.9.10--pyhdfd78af_0.img` — a recorded image path, not a verified digest.
- `/project/6090753/active/ovcan_gq_analysis/work/` — staged Mutect2 reference paths; for example `7b/925b54083ad554059d4fc85b402ae6/1000g_pon.hg38.vcf.gz` and `af-only-gnomad.hg38.vcf.gz` for TOV3392D.
- `../work/<two-character-prefix>/<task-hash>/<sample>_P<passage>.recal.bam` and `normal_samples/all_bams/*.bam` in CNV commands. Resolve these relative to the original command's working directory established from run records, not relative to the Git checkout or workstation archive.
- VCF headers name staged `<sample>_P<passage>.recal.cram` inputs. BAM and CRAM may be different representations/stages; compare read groups, reference and generating tasks rather than assuming byte equality.

Important matching rules:

1. `OV` and `TOV` prefixes are meaningful. `OV2295`, `OV2295-R2` and `TOV2295-R` are three different models from one patient. `OV3133-R`, `OV3133-R2`, `TOV3133D` and `TOV3133G` are four different models. `OV3291` is not `TOV3291G`.
2. Separators and suffixes vary: canonical `OV1369-R2` maps to WES `OV1369_R2_P66`; CNV directories append `_new`; MAF filenames sometimes omit passage. Remove separators only for a candidate join, then verify uniqueness and sample/passages. Do not merge `R` with `R2`, `D` with `EP/G`, or samples merely sharing a patient number.
3. RNA and WES passages differ. For example, OV1369-R2 has RNA `p69` but WES `P66`; TOV112D has RNA `p63` but WES `P83`. The WES-specific fields in the model checklist take precedence for this search.
4. TOV3121D and TOV3121EP both have WES P68 but remain distinct models. TOV3121D's VCF is compressed; the archived alternatives are mostly plain `.vcf` despite surviving `.vcf.gz.tbi` filenames. A leftover index does not prove its compressed VCF is present or usable. Inventory compression/index relationships; do not regenerate files during discovery.
5. `wes - old` is the mutation module and `cnvkit wes - new` the CNV module. They are not two competing CNV call sets. Within the CNV archive, “All samples” is only a subset; top-level directories and `20251208 - new samples/` hold additional models, while CC/EC/HGS folders contain duplicates. Match hashes and model IDs before counting files.
6. Old summary spreadsheets overstate WES availability. The baseline is 23 models, not every model in the 42-model resource. If additional WES libraries are found, list them separately with evidence; do not silently expand the manuscript cohort.

### Corrected path-table extraction issue

Six lines in the original `commands.txt` contain an extra bare sample ID immediately after `cnvkit.py batch`, before the actual BAM token. An earlier handoff-table parser captured that ID as the alignment path. Both maintained CSVs now contain the actual `.recal.bam` token for all 23 models. The corrected rows are:

| Model | Original command line | Recorded BAM token |
|---|---:|---|
| OV1369-R2 | 18 | `../work/a5/a947cdda816a6ac440c35f5a99d7c7/OV1369_R2_P66.recal.bam` |
| TOV81D | 19 | `../work/a5/c66150f956d90d786643f4b3d08453/TOV81D_P23.recal.bam` |
| TOV2414 | 20 | `../work/ad/34d92ef2d6e1ba0272c5fd188d3a96/TOV2414_P65.recal.bam` |
| TOV2295-R | 21 | `../work/b6/ecdb564b1632e569b8ba3584395884/TOV2295_R_P57.recal.bam` |
| TOV2929D | 22 | `../work/e4/cf2272b3a898107c8fd5d1b2d00216/TOV2929D_P57.recal.bam` |
| TOV2881EP | 23 | `../work/ec/03f6bb30a3faac68cadaac88840ace/TOV2881EP_P64.recal.bam` |

Source: original `judy_archive/data/cnvkit wes - new/commands.txt`, retained as [cnvkit_commands.txt](evidence/cnvkit_commands.txt) with the same line numbering. Assertions verified 23 model rows, 23 `.recal.bam` tokens, 322 total path hints and exactly six corrections in each handoff CSV. Recover the executed task commands/logs to resolve the extra positional identifiers; **do not execute or “repair and rerun” these historical lines**. The table repair proves what token was recorded, not that the file exists now or that this text was the successful invocation.

## Read-only discovery sequence

1. Record the Git commit, cluster hostname, current directory and actual discovery time. Read root `CLAUDE.md`/`AGENTS.md` if present. Test the literal roots above and record readable, absent or permission-denied separately. Avoid a whole-cluster or `/` scan.
2. Inspect the known project roots and the shallow run/output directories first. Use `rg --files --hidden --maxdepth 4` with large-data extensions and deep `work/` excluded, then search the resulting names for `pipeline_info`, `multiqc`, `samplesheet`, `params`, `.config`, `.nextflow.log`, `trace`, `report`, `versions`, scheduler scripts and QC files. These names are discovery patterns; module existence must be demonstrated from the actual files.
3. Read matching small scripts/configs/logs as text. Do not source them. Connect the run record, input samplesheet, successful trace entries and final-output identifiers. A directory timestamp or date-looking folder name is not a verified run date.
4. Follow exact path hints and task hashes next. Expand to a bounded search within the identified project's `work/` only when needed. Preserve `.command.*` and `.exitcode` alongside each other and associate them with the trace. Do not assume a staged symlink target survives because the symlink itself exists.
5. Inventory large input files with path/size/accessibility and existing checksums. Available `samtools view -H` or `samtools quickcheck` can provide lightweight header/integrity checks without a full read scan. If tools are unavailable, record that; do not install software for discovery. Do not call `flagstat`, depth/coverage tools, variant callers or whole-file hashing across the cohort merely to fill the table. Calculate SHA-256 for the small evidence actually copied; use existing large-file checksum manifests where available and leave missing large-file hashes explicit.
6. Copy only the selected small evidence into a separate retrieval staging directory, preserving per-run structure and original paths in a manifest. Keep credentials out of copied configs. Leave FASTQs, alignments, large VCFs and reference genomes on the cluster; do not add them to Git. Do not overwrite the existing audited `output/` tables, raw calls or manuscript during this recovery pass.

No `sbatch`, Nextflow launch/resume, container execution of historical pipeline commands, alignment/variant/annotation rerun, new MSI/HRD job or bulk transfer is required. If existing QC is insufficient, return a specific proposed later computation with its required inputs and expected output; do not launch it as part of this task.

## Deliverables and acceptance checks

Write a small recovery report and machine-readable inventories under a new dated directory under `docs/cluster/recovery/` (a new handoff output location, not a claimed historical path). Large or provider-original evidence can stay in an ignored local staging directory such as `data/cluster_wes_retrieval/`, with repository-safe summaries and manifest pointers. Use the repository's final packaging policy for small originals.

Required outputs:

- `RECOVERY_REPORT.md`: found/not found/partial findings by H1–H4 and M1–M3; roots actually searched; run identities; resolved alias conflicts; completed versus still-open manuscript details; a precise next action for each remaining gap.
- `file_inventory.tsv`: `record_id`, `model_or_reference`, `sample_id`, `passage`, `run_id`, `role`, `original_path`, `resolved_path`, `status`, `bytes`, `checksum_algorithm`, `checksum`, `checksum_evidence`, `copied_relative_path`, `source_evidence`, `observed_at`, `notes`. Separate fresh observation time, historical run date and copy time. Do not infer any of them from mtime.
- `model_status.tsv`: exactly the 23 baseline models, with canonical/passaged aliases, matched run/library/read-group evidence, alignment/index availability, QC coverage and explicit unresolved items. Put extra discovered models in a separate section/table.
- `qc_metrics.tsv`: a long-form table with `model_or_reference`, `run_id`, `metric`, `value`, `unit`, `denominator`, `target_file_id`, `source_path`, `status` and `notes`. Keep model and public-reference exome groups distinguishable. Record absent/not-computed metrics as such.
- `run_reference_manifest.tsv` (or equally explicit JSON): workflow release/commit, Nextflow and tool versions, invocation/config source, reference and interval identities/hashes, padding/coordinate conventions, normal-reference membership and executed-command evidence. Unknown fields remain unknown.

The recovery pass is complete when every baseline model and priority item has an evidence-backed status. Success does not require pretending every file is recoverable. Acceptance checks are:

1. All 23 canonical WES models are accounted for exactly once; 16 patients and 11 HGSC patients are preserved. Each of the five public CNV reference exomes is accounted for separately. No alias or passage conflict is silently resolved.
2. Path-table clues are marked as hints until the cluster file is observed. A broken symlink, unreadable directory and a negative search in limited roots have distinct statuses; “not found in the searched roots” does not become “never existed.”
3. Retrieved run/command/config evidence can be linked to the retained calls. Final VCF candidates are compared with the audited manifest hashes where possible; a changed checksum is investigated, not overwritten or assumed to represent a better run.
4. Capture/vendor, calling/shard and CNV-bin intervals are distinguished. The recovered derived BED hash/count remains the baseline; sequence reference identity is documented without lifting over the audited coordinates.
5. QC has traceable numerical sources, units and denominators. Missing depth/coverage metrics stay open until existing reports or a separately authorized calculation resolves them.
6. Small copied files have hashes and provenance; raw/large data have inventory pointers. The retrieval work has not launched compute jobs or altered raw data, audited results, cohort membership or biological claims.

The maintained analysis entry points are [07_wes_mutations.R](../../scripts/07_wes_mutations.R), [08_wes_cnv.R](../../scripts/08_wes_cnv.R), [09_wes_hrd.R](../../scripts/09_wes_hrd.R), [16_wes_signatures_msi.R](../../scripts/16_wes_signatures_msi.R) and [22_wes_signature_refit.R](../../scripts/22_wes_signature_refit.R). Read them to understand what evidence downstream analyses consume; running them is outside the discovery task.
