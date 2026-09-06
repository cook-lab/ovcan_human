# Targeted WES follow-up — inventory, 2026-09-06

Response to `../2026-09-05/FOLLOWUP.md`. Read-only discovery on Nibi login node `l5`, plus
one trivial local Nextflow test run (a two-line workflow, to learn what the "revision" field
means). No job submitted, no source file modified. Item-level statuses, paths and hashes are
in `followup_inventory.tsv`; reference-file hashes in `followup_inventory_reference_hashes.tsv`;
the new bundle's per-file checksums in `cnv_coverage_bundle_SHA256SUMS.txt`.

## Bundle for local validation

```
/scratch/dcook/ovcan_human_wes_cnv_coverage_2026-09-06.tar.gz
303,052,628 bytes   sha256 c5d79958ae05115c975d63e3858fae7cd372b37c593295e0e2472bbf38c6706a
```

109 files (688 MB unpacked): `normals/` (the five `SRR*.sorted.{target,antitarget}coverage.cnn`,
`reference.cnn`, `intervals_sorted.target.bed`), `models/` (23 × target and antitarget
coverage `.cnn`), `plots/` (23 diagram PDFs and scatter PNGs), the four command/script files,
and `SHA256SUMS`. The normal coverage files and `reference.cnn` are byte-identical in all 23
model output directories (verified from the project manifest), so one copy each is included.
Permanent copies remain in `/project/6090753/active/ovcan_human/wes/cnvkit_0.9.10_manual/`.

## Status by requested item

| Item | Status | Summary |
|---|---|---|
| Five normal `*coverage.cnn`, model coverage files, plots | **found** | All present and bundled. Normal target depth 75.9–95.7×; antitarget 0.45–0.63× with ~2.5 % zero bins. Model target depth 78.6–98.7×; **model antitarget bins 42,705/42,709 zero in every model**. |
| Why model antitargets are empty | **found (explained)** | The manual CNVkit inputs were BAMs converted from Sarek's recalibrated CRAMs, and Sarek wrote those CRAMs by per-interval-shard ApplyBQSR, so off-target reads were physically absent from every tumour input. The normals were whole-genome bwa-mem2 alignments and keep off-target depth. Sarek's own CNVKIT_BATCH log shows antitarget `#reads=492, mean=0.0115` for OV1369_R2_P66. The zero-depth pattern is a processing artefact, not deletion; profiles are target-supported. |
| Manual CNVkit execution logs | **not found in searched locations** | No `slurm-*` from `run_cnvkit.sh` in the analysis directory or the scratch tar. Slurm accounting on Nibi holds no records before 2026-09-01 (PrivateData is `none`, so this is retention, not permission), so the Nov–Dec 2025 jobs cannot be recovered from `sacct`. Asma's home was not searched. |
| Manual CNVkit execution, indirect evidence | **partially found** | `commands.txt` (2025-12-03 12:58), `commands1.txt` (six corrected lines, 13:11), `run_cnvkit.sh` pointing at `commands1.txt` (13:12); 20 output directories written 2025-11-14 16:01–16:08, the six corrected models rewritten 2025-12-03 (last 17:52); PDF creation stamps agree; identical `reference.cnn` and normal coverage across all 23 directories. Saved commands plus timestamped outputs, no invocation log. |
| CNVkit input alignments | **found** | Each `work/<hash>/<model>.recal.bam` (6.8–8.6 GB, 2025-10-24 19:04–19:07) was created by `samtools view -T` from `<model>.recal.cram`, a symlink to the **Oct-24 recalibrated CRAM in `work/`**, and all 23 of those source CRAMs still exist (4.4–5.6 GB). For the 13 models whose published CRAM was rewritten on Oct 31, the Oct-24 work copy is the exact CNVkit input and remains available. |
| Manual CNVkit tool version | **partially found** | `cnvkit-0.9.10--pyhdfd78af_0.img` (biocontainers tag) named in the commands; the image file is in `/scratch/asmab/nf-core` (permission denied) and not in the tar. Plot metadata: Matplotlib 3.7.0, ReportLab. |
| Original 22-model vcf2maf command/version | **not found in searched locations** | Nothing MAF-related in the analysis directory, `work/` task names, or the tar index. Consistent with workstation-side conversion. |
| Sarek workflow source identity | **found (verified)** | Nextflow's "revision" for a path-launched script is a content hash, not MD5 or a Git commit. Nextflow 24.10.2 assigns the released nf-core/sarek 3.5.1 `main.nf` revision `3954909713`, identical to all eight launch logs, so the launched `main.nf` was the unmodified release file (3.5.0 and 3.5.1 share it; the release is fixed by the `3_5_1` directory and the `nf-core/sarek: v3.5.1` manifest line). The source directory itself is inaccessible and not in the tar. |
| Container image identities | **found (tags, not digests)** | Fourteen Singularity images resolved from `/scratch/asmab/nf-core/` are listed in `followup_inventory.tsv` (gatk4 4.5.0.0, ensembl-vep 113.0, fastp 0.23.4, fastqc 0.12.1, samtools 1.21, htslib 1.20, bcftools 1.20, mosdepth 0.3.8, vcftools 0.1.16, multiqc 1.25.1, gawk 5.1.0, two mulled images, one Seqera blob with sha256). Image files are inaccessible. |
| Reference FASTA and resource checksums | **found (computed)** | Streamed from the tar: `Homo_sapiens_assembly38.fasta` 3,249,912,778 B, MD5 `7ff134953dcca8c8997453bbb80b6b5e` (the widely published Broad hg38 bundle value; confirm against the bundle's own `.md5` if desired), SHA-256 `93157a16…5564`; plus `1000g_pon.hg38`, `af-only-gnomad.hg38`, `dbsnp_146.hg38`, `Mills_and_1000G…hg38` (full values in the companion TSV; byte counts match the tar index). `Homo_sapiens_assembly38.known_indels.vcf.gz` was staged from `GATKBundle/beta/`, which is **not in the tar**: not found in searched locations. |
| Liftover command, version, chain | **partially found** | Local binary `liftOver`, kent source version 483 (file dated 2025-06-27). Output names (`hglft_genome_2576bc_8f0.bed`, `#Deleted in new` list) match UCSC hgLiftOver web output; timeline 2025-10-03 11:07 → 12:59 → 13:20 → 14:58. No chain file and no command line in any searched root; the name/score columns were added after lifting. |
| `TPV81D_23-pool`, DNA/library SOP | not a cluster item | Only the samplesheet (lane 6) and provider BAM exist here. |

## Notes

- Nothing here changes cohort membership, variant candidates, CNV segments or the QC tables.
- Raw logs stay in `/project/6090753/active/ovcan_human/wes/sarek_3.5.1_run/logs/` (private);
  this directory carries only curated values.
- Hashing the five reference files read ~9.9 GB from the tar on the login node in one pass;
  no other heavy operation was run.
