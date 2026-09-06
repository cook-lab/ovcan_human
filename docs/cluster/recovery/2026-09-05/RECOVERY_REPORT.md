# WES recovery report — Nibi cluster, 2026-09-05

Discovery performed 2026-09-05 (17:59–20:30 UTC) on login node `l5.nibi.sharcnet` from the
public repository clone `~/ovcan_human` (commit `ed01287`). Read-only throughout: no jobs
submitted, no pipeline run, no file in the source directories modified or deleted. Small
evidence was copied to a new directory owned by `dcook`,
`/project/6090753/active/ovcan_human/` (6.0 GB, 2,139 files, `MANIFEST.sha256`), and the
originals were left where they are. Paths below written as `wes/...` are relative to that
new directory; `SRC` = `/project/6090753/active/ovcan_gq_analysis` (owner `asmab`);
`R004741` = `/project/6090753/active/ovcan_gq/ammmasson_OvCAN_WES_bioinformatics_analysis_R004741`
(owner `dcook`).

Companion files in this directory: `file_inventory.tsv` (2,409 records), `model_status.tsv`
(23 baseline models), `additional_wes_models.tsv` (22 further samples with WES-derived data,
kept out of the cohort), `qc_metrics.tsv` (1,494 long-form rows), `run_reference_manifest.tsv`.

## Roots searched and their state

| Root | State | What it holds |
|---|---|---|
| `SRC` = `/project/6090753/active/ovcan_gq_analysis` | readable (group rwx) | the definitive nf-core/sarek 3.5.1 run (`wes_sarek_results_20251022/`, `work/` 99,167 task dirs, 8 launch logs), the manual CNVkit analysis (`cnvkit/`), and `ovcan.tar.gz` |
| `SRC/ovcan.tar.gz` | readable; **truncated** plain tar misnamed `.tar.gz`; 315.5 GB, 5,663 members indexed before `Unexpected EOF` | copy of the old scratch root `/scratch/asmab/ovcan/` made 2026-05-06: samplesheet, design/interval files, `liftOver` binary, references (63 GB), superseded earlier Sarek attempts (`wes_new/` 211 GB), a 3-sample scratch copy of `wes_sarek_results_20251022` |
| `/scratch/asmab/ovcan/` | **permission denied** (mode 700) | historical root named in every command; partially preserved in the tar |
| `R004741` | readable (mine) | McGill Genome Centre delivery: 23 provider BAMs (Sarek input), SeqCap EZ Exome v3 design files, 2022 GenPipes/PCGR GRCh37 results for 32 samples, spreadsheets `WES 2014/2017/2018 Mes-Masson.xlsx` |
| `/project/6090753/dcook/projects/ovcan_gq` | empty | — |

Two reads were refused by the session's tool-permission classifier and were **not** done:
parsing the three `WES 20xx Mes-Masson.xlsx` spreadsheets, and `samtools view -H` on a
provider BAM. Both are trivial to run by hand (commands at the end).

## Headline findings

1. **The retained variant calls are fully traceable.** All 23 `*.mutect2.filtered_VEP.ann.vcf.gz`
   in `SRC/wes_sarek_results_20251022/annotation/mutect2/` hash-match (uncompressed SHA-256)
   the `vcf_sha256` values in `output/wes_input_manifest.csv`. They were produced by one
   successful launch: nf-core/sarek **v3.5.1**, Nextflow **24.10.2**, run name `angry_allen`,
   session `26cf85bf-0ac3-4237-8f48-c36b2ab20bd2`, Slurm job 3273686, started 2025-10-23 15:35
   EDT, completed 23:48:40 EDT, 38,541 tasks succeeded. Launch script, config, resolved
   `params_*.json`, traces, software versions and the `.nextflow/history` are all present.
2. **The capture kit is documented.** The McGill delivery ships `SeqCap_EZ_Exome_v3.targets.bed`
   (242,232 hg19 targets, 63.56 Mb) and `.probes.bed` (368,146 probes): Roche NimbleGen
   SeqCap EZ Exome v3. The scratch root holds the liftOver chain of evidence from that file to
   the hg38 `intervals_sorted.bed` used by Sarek (242,421 intervals; 17 not lifted). The
   repository's 290,475-bin CNVkit BED is the autobin split of that file (byte-identical copies
   found). The "same kit" author statement now has a documentary basis; what remains
   undocumented on the cluster is library-prep protocol, kit lot, instrument, flowcell and run
   dates (expected in the unread spreadsheets and McGill QC reports).
3. **Sequencing was delivered as aligned BAMs, not FASTQs.** The Sarek samplesheet lists a
   `bam` column pointing at `R004741/bam/<sample>.bam` (GRCh37 alignments from McGill). Sarek
   converted them to FASTQ, split them (50 M reads), and re-aligned to GRCh38 with bwa 0.7.18.
   No FASTQ exists anywhere searched. The 23 provider BAMs (6.4–8.1 GB each, ~154 GB) are
   present and readable; they are the primary input for any deposition.
4. **Complete technical QC exists** for all 23 models: MultiQC (2025-10-23) with fastp,
   FastQC, MarkDuplicates, samtools stats, mosdepth over the target BED, bcftools/vcftools
   stats and VEP summaries; per-sample source files were copied. `qc_metrics.tsv` carries
   1,469 measured values with units and denominators. The five public reference exomes have
   **no** alignment QC (aligned outside Sarek); only their CNVkit coverage files exist.
5. **The retained recalibrated CRAMs are not the byte-identical inputs of the retained VCFs.**
   Later launches (Oct 24 with `cnvkit` added; Oct 28–31, all killed by the 12h59 walltime
   after a lost Nextflow session) re-executed preprocessing and overwrote the 23
   `*.recal.cram`: 10 now date from Oct 24, 13 from Oct 31. Same inputs, same pipeline and
   parameters, but equivalence is by construction, not by checksum. All 23 pass
   `samtools quickcheck`.
6. **Recalibrated CRAMs are interval-restricted.** Sarek applied BQSR per interval shard
   (405 shards) and merged, so `recal.cram` holds only reads overlapping the shards (e.g.
   OV2085_P64: 77.4 M of 113.5 M reads). Genome-wide denominators (mapped %, duplicates,
   on-target fraction) must come from the `md.cram` statistics, which are also retained.

## Status by priority item

### H1 — Run identity and sample mapping: **found**

- Workflow: nf-core/sarek v3.5.1 launched by path from a local copy
  (`/scratch/asmab/nf-core/nf-core-sarek_3.5.1/3_5_1`), script revision
  `3954909713023f4328e976337e6e2cb9`; Nextflow 24.10.2 build 5932; plugins nf-schema 2.3.0,
  nf-amazon 2.9.2; Apptainer containers (digests not recorded).
- Invocation (`wes/sarek_3.5.1_run/launch/run_sarek.sh`): `--wes --tools mutect2,vep
  --genome GATK.GRCh38 --igenomes_base /scratch/asmab/ovcan/references --intervals
  /scratch/asmab/ovcan/intervals_sorted.bed --nucleotides_per_second 600000 -resume`, config
  `obcf-graham.cfg` (Slurm executor, queue 80, `-x c32`, two memory overrides). Despite its
  name the config describes "The OBCF Nibi DRAC profile"; the run was on Nibi.
- Eight launches 2025-10-22..31 across three sessions; only `angry_allen` completed. The
  Oct-24 `angry_kilby` launch (tools `mutect2,vep,cnvkit`) completed CNVKIT_BATCH for all 23
  with Sarek's flat reference, then failed at GETPILEUPSUMMARIES (TOV3133D_P66) after Slurm
  submission failures. Oct-28 `fervent_escher` started a new session (cache lost) and failed
  in BWAMEM1_MEM (TOV2414_P65); Oct-29/30/31 launches were killed at the walltime. Full
  table in `run_reference_manifest.tsv` (`launch_history`).
- Earlier scratch attempts (Oct 3–24; outdirs `wes_new`, `wes_sarek_results_20251007`,
  a 3-sample `wes_sarek_results_20251022`, `wes_sarek_results_20251024`; tools at times
  including strelka and controlfreec) are recorded from the tar and are superseded. Do not
  confuse the scratch `wes_sarek_results_20251022` (3 samples, 26.7 GB) with the project one.
- Sample mapping: `samplesheet.csv` (23 rows, `patient==sample`, `sex=XX`, `status=1`,
  tumour-only); lane = provider BAM `PU` tag via `get_lane_info.sh`; read groups
  `ID:<sample>_<lane> PU:<lane> SM:<sample>_<sample> LB:<sample> PL:ILLUMINA`. Lanes per
  model are in `model_status.tsv`. No alias conflicts: every Sarek sample id maps 1:1 to a
  repository model and passage (`OV1369_R2_P66`↔`OV1369-R2` P66, etc.); TOV3121D_P68 and
  TOV3121EP_P68 are distinct throughout.

### H2 — Capture, reference and acquisition provenance: **found (design, reference); not found (acquisition details)**

- Design: SeqCap EZ Exome v3 (`wes/mcgill_r004741_provenance/SeqCap_EZ_Exome_v3.{targets,probes}.bed`).
- Interval derivation: targets (hg19, 242,232) → UCSC `liftOver` → `hglft_genome_2576bc_8f0.bed`
  (242,421; `hglft_genome_ignored.txt` lists 17 intervals "Deleted in new"; 208 intervals
  landed on 14 alt/random/Un contigs) → `intervals_sorted.bed` (63,709,951 bp summed,
  63,514,049 bp merged) → 405 Sarek shards. No padding was applied. The Picard
  `interval_list` copy served only the first scratch attempts. CNVkit split the same file into
  290,475 target bins and 42,709 antitarget bins.
- Reference: GATK GRCh38 iGenomes bundle `Homo_sapiens_assembly38.fasta` (dict and fai
  extracted and hashed; FASTA itself only inside the truncated tar); known sites
  dbsnp_146, Mills_and_1000G, Homo_sapiens_assembly38.known_indels; Mutect2
  `--panel-of-normals 1000g_pon.hg38.vcf.gz` and `--germline-resource af-only-gnomad.hg38.vcf.gz`
  (the latter also for GetPileupSummaries); VEP cache 113_GRCh38.
- Not in cluster records: sequencing centre confirmation, instrument, flowcell/run IDs and
  dates, library kit lot, read configuration beyond 2×100 bp. The provider BAM headers and
  the three spreadsheets are the obvious next source (see "Next actions").

### H3 — Existing WES QC: **found (23 models); not produced (5 reference exomes)**

Per model (all in `qc_metrics.tsv`, MultiQC of 2025-10-23): fastp read counts before/after
filtering, Q30, GC, duplication, insert-size peak; MarkDuplicates pair counts, duplication
fraction, library size; samtools stats for md and recal CRAMs (mapped, properly paired, MQ0,
error rate, insert size); mosdepth mean depth genome-wide and over target, fraction of target
bases ≥1/10/20/30/50/100x (md and recal); GATK contamination estimate; CNVkit
reads-in-target; bcftools/vcftools counts. Summary over the 23 models:

| Metric (md CRAM unless noted) | Range |
|---|---|
| total reads after fastp (R1+R2) | 108.4–136.8 M |
| reads mapped | 99.83–99.98 % |
| MarkDuplicates duplication | 14.9–24.0 % |
| mosdepth mean target depth | 69.7–88.8× |
| target bases ≥30× | 66–84 % |
| GATK contamination estimate | see table (OV2085_P64 = 0.0039) |
| CNVkit target depth, manual run (derived) | 78.6–98.7× |

Caveats recorded in the table: recal-CRAM statistics are interval-restricted; the Sarek
CNVkit "percent reads in regions" (61–68 %) is computed on that restricted BAM, so it is not
a genome-wide on-target fraction; mosdepth counts duplicate reads. Normal exomes: only
CNVkit coverage (target depth 75.9–95.7×, derived here from `targetcoverage.cnn`); every
other metric is marked `not produced`.

Variant counts: vcftools `PASS` records per sample sum to 19,816 and TOV3121D_P68 has 547,
versus the audited 16,081 / 389. The units differ (VCF records including multi-allelic
sites vs MAF rows after vcf2maf), so reconcile before quoting either.

### H4 — Input preservation inventory: **found**

| Item | Count / size | Location | Status |
|---|---|---|---|
| provider BAM + index (Sarek input) | 23 × 6.4–8.1 GB (~154 GB) | `R004741/bam/` | accessible; no checksums |
| recalibrated CRAM + crai | 23 × 4.4–5.6 GB (~115 GB) | `SRC/wes_sarek_results_20251022/preprocessing/recalibrated/` | accessible; quickcheck OK; re-executions (Oct 24/31) |
| duplicate-marked CRAM + crai | 23 × 3.1–4.0 GB (~81 GB) | `.../preprocessing/markduplicates/` | accessible |
| CNVkit input BAMs (from Oct-24 CRAMs) | 23 × 6.8–8.6 GB (~175 GB) | `SRC/work/<hash>/` | accessible; all 23 historical `../work/...` tokens resolve |
| normal exomes: SRA, FASTQ×2, SAM, BAM, sorted BAM+bai | 5 × ~85 GB (~435 GB) | `SRC/cnvkit/normal_samples/SRR*/` | accessible; only `*.sorted.bam` (36 GB) is used |
| mosdepth per-base depth | 46 files, 17 GB | `.../reports/mosdepth/` | accessible; not copied |
| `ovcan.tar.gz` | 315.5 GB | `SRC/` | truncated; indexed (`wes/inventory/ovcan_tar_gz_index_2026-09-05.txt`) |
| 9 orphan `.bai` (no BAM) | — | `R004741/bam/` | BAMs not found in searched roots |

No checksum manifest exists for any large file in the searched roots. Checksums were
computed only for the 2,139 copied files.

### M1 — Filtering and annotation chain: **found**

Per-sample `mutect2.vcf.gz(.stats)`, `filtered.vcf.gz`, `filteringStats.tsv`,
`contamination.table`, `segmentation.table`, `artifactprior.tar.gz`, `pileups.table` and the
VEP-annotated VCF are retained and copied. Executed `.command.sh` files for FASTP, BWAMEM1_MEM,
MARKDUPLICATES, BASERECALIBRATOR, APPLYBQSR, MERGE_CRAM, MUTECT2, GETPILEUPSUMMARIES,
CALCULATECONTAMINATION, LEARNREADORIENTATIONMODEL, MERGEMUTECTSTATS, FILTERMUTECTCALLS,
ENSEMBLVEP_VEP, MOSDEPTH, SAMTOOLS_STATS, BCFTOOLS_STATS, FASTQC and the interval
preparation are copied under `wes/sarek_3.5.1_run/task_evidence/` with exit codes and
directory listings. FilterMutectCalls ran with GATK defaults (listed in the VCF header);
VEP ran `--everything --filter_common --per_gene --total_length --offline`. No custom
filtering happened on the cluster; the cascade to 6,194 candidates is repository-side.
The vcf2maf step for the 22 archived MAFs was done on the workstation; no MAF exists on the
cluster.

### M2 — CNV reference detail: **found**

`commands.txt` (byte-identical to the repository copy), `commands1.txt` (the six corrected
lines without the stray sample-ID token, i.e. the re-runs for OV1369_R2_P66, TOV81D_P23,
TOV2414_P65, TOV2295_R_P57, TOV2929D_P57, TOV2881EP_P64), `command_diagram.txt`,
`run_cnvkit.sh`, the normal-sample download/align/sort/index scripts and their Slurm logs
(2025-11-13/14), and all 23 output directories (`.cnr`, `.cns`, `.call.cns`, `.bintest.cns`,
coverage `.cnn`, `reference.cnn`, diagram, scatter) are copied. The pooled `reference.cnn` is
identical (SHA-256) across the 23 directories, as are the per-normal coverage files.
Command: `cnvkit.py batch <recal.bam> --drop-low-coverage -n normal_samples/all_bams/*.bam
-f <FASTA> -t intervals_sorted.bed -a intervals.antitarget.bed --diagram --scatter`; no
`--segment-method`, no `-y` (CNVkit defaults). The earlier hypothesis that the stray token
lines were never valid is resolved: those six models were re-run from `commands1.txt` and
their outputs exist (dated 2025-11-14..12-03). The Slurm `.out/.err` of the manual CNVkit
array jobs were not found. The Sarek-internal CNVkit run (Oct 24, flat reference, 663 MB)
also exists and is superseded; it is inventoried, not copied.

### M3 — Previously completed specialised work: **found (other pipeline, GRCh37); inventory only**

`R004741/USB/` holds the McGill GenPipes DNAseq results (HaplotypeCaller + VQSR joint VCF
split per sample, CNVkit germline CNA tables) and PCGR 0.9.2 reports (tumour-only, GRCh37,
with `--estimate_msi_status --estimate_signatures --estimate_tmb`) produced by a McGill
analyst in April 2022 for 25 of 32 samples; seven of the 23 baseline models
(OV2085_P64, OV3291_P46, OV90_P63, TOV112D_P83, TOV2414_P65, TOV3392D_P69, TOV81D_P23) are
listed in `missing_reports`. These are germline-caller-based, GRCh37, PCGR-tiered outputs:
useful as historical context, not as somatic or HRD evidence. No ASCAT, allele-specific CN,
purity/ploidy or HRD analysis exists on the cluster. Nothing here changes the manuscript's
position that MSI and HRD are not established.

## Additional WES-derived material outside the 23-model cohort (kept separate)

`additional_wes_models.tsv` lists 22 further samples with data in `R004741`:

- 9 Mes-Masson lines sequenced with the cohort (index files, GenPipes VCFs, CNA and PCGR
  reports present; **BAMs absent** in searched roots): OV1946_P49, OV2978_P67, OV4453_P63,
  OV4485_P60, TOV1946_P49, TOV2223G_P69, TOV2978G_P67, TOV3041G_P52, TOV3291G_P65.
- BIN67 and CSCOE (SCCOHT): VCF, CNA, PCGR; capture **SureSelect Human All Exon V7**, so not
  comparable to the SeqCap cohort.
- 11 VOA/OV866-2 lines (VOA10816p35 … VOA8771p26, OV866-2p70): VCF, CNA, PCGR only.
- `Public_Data/` CL-1…CL-15: identities not established.

**Author decision (David Cook, 2026-09-05): omit all of these from the manuscript and the
data record.** Reason: the McGill raw data (FASTQ/BAM) for these lines were lost in a storage
failure and were never retained locally; only GRCh37 GenPipes/PCGR outputs survive, produced
with a germline-design caller (HaplotypeCaller + VQSR), different build, annotation and CNV
reference, so they are not comparable with the Sarek/Mutect2 cohort. The 42-model,
34-patient, 23-WES counts in v7 are unchanged. The files stay where they are in `R004741`
for the lab's own reference.

## Acceptance checks

1. 23 canonical models accounted for exactly once; 16 patients and 11 HGSC patients
   preserved (asserted in the generator). Five reference exomes listed separately. No alias
   merged.
2. Every historical `../work/...` hint was observed on disk; `/scratch/asmab/ovcan` is
   recorded as permission-denied, the tar as truncated, and the nine extra models' BAMs as
   "not found in searched roots".
3. Run → outputs linkage established by trace, log, mtime and 23/23 VCF hash matches.
4. Vendor targets (242,232 hg19), lifted calling intervals (242,421 hg38), Sarek shards (405)
   and CNVkit bins (290,475 / 42,709) are distinguished; no coordinates were lifted or altered.
5. QC values carry source file, unit and denominator; missing metrics are `not produced`.
6. Copied files have SHA-256 in `/project/6090753/active/ovcan_human/MANIFEST.sha256`;
   large files have inventory pointers. No job launched, no raw data altered, no cohort or
   result changed.

## What the manuscript can now state (WES methods)

Exome capture with SeqCap EZ Exome v3 (Roche NimbleGen), 2×100 bp Illumina paired-end
sequencing (centre and instrument to confirm from provider records), reads re-extracted from
provider GRCh37 BAMs and processed with nf-core/sarek v3.5.1 (Nextflow 24.10.2): fastp 0.23.4,
bwa 0.7.18 to GATK GRCh38 (`Homo_sapiens_assembly38`), GATK 4.5.0.0 MarkDuplicates and BQSR
(dbSNP 146, Mills, known indels), tumour-only Mutect2 4.5.0.0 with the GATK 1000 Genomes
panel of normals and gnomAD germline resource, FilterMutectCalls with contamination,
segmentation and orientation-bias models, Ensembl VEP 113 (GRCh38.p14 cache; gnomAD v4.1,
ClinVar 202404, dbSNP 156). Copy number: CNVkit 0.9.10 hybrid method against a pooled
reference of five public exomes (PRJNA339046) aligned with bwa-mem2, target bins derived from
the same lifted SeqCap intervals. QC via MultiQC (values in `qc_metrics.tsv`).

## Gaps that remain and the precise next action for each

| Gap | Next action |
|---|---|
| Sequencing centre, instrument, flowcell/run IDs, dates, library kit lot | Parse `R004741/WES 2014|2017|2018 Mes-Masson.xlsx` (copied to `wes/mcgill_r004741_provenance/`) and read one provider BAM header: `module load samtools; samtools view -H "R004741/bam/OV2085_P64.bam" | grep -E '^@(RG|PG)'`. Both were blocked for the agent, not for you. |
| Same-kit confirmation for the 5 public normals | Check SRA/BioProject PRJNA339046 metadata for the capture kit; cluster records only show they cover the SeqCap targets at 76–96×. |
| Byte identity of alignments behind the VCFs | Not recoverable (Oct-23 CRAMs overwritten). State in methods that retained CRAMs are re-executions of the identical workflow. |
| Genome-wide on-target fraction | Not in existing reports (mosdepth ran with `--by`; Picard HsMetrics not run). Proposed later job, not launched: `mosdepth --by intervals_sorted.bed` summary already gives target vs genome mean; for read-level on-target use `samtools view -c -L intervals_sorted.bed <md.cram>` per model (23 × ~2 min, one small Slurm job). |
| PASS-count reconciliation (19,816 vcftools vs 16,081 audited) | Count PASS records in the copied `filtered.vcf.gz` with the repository's vcf2maf-equivalent normalisation before quoting either number. |
| Manual CNVkit job logs | Not found in searched roots; ask Asma whether `slurm-*_[1-27].out` from `run_cnvkit.sh` were kept elsewhere. |
| Container digests, liftOver chain/version, bwa-mem2 version | Not recorded anywhere found; report as "not recorded". |

## Housekeeping observations (no action taken)

- `R004741/sftp` (124 B) may contain provider credentials and sits in a group-readable
  project directory. It was neither read nor copied. Consider moving it to `~` with mode 600.
- Large regenerable intermediates in `SRC` (Asma's directory, not mine to touch): normal-sample
  SAM (~170 GB), unsorted BAM (~56 GB), uncompressed FASTQ (~162 GB) and 27 `tmp.*.bam`;
  the truncated 315 GB tar (its only unique content is the reference bundle and superseded
  runs); `work/` (99,167 task directories). Roughly 1 TB could be reclaimed, but every item
  is her call.
- Quota after copying: `/project` 8,117 GiB of 35 TiB, 3.03 M of 10 M files.

## Commands you can run yourself for the two blocked reads

```bash
# spreadsheets (script in the session scratchpad; pure standard library)
python3 /tmp/claude-3055330/-home-dcook/b0c11c8f-b271-4cf7-9676-37f8757b0bda/scratchpad/dump_xlsx.py \
  "/project/6090753/active/ovcan_human/wes/mcgill_r004741_provenance/WES 2014 Mes-Masson.xlsx" \
  "/project/6090753/active/ovcan_human/wes/mcgill_r004741_provenance/WES 2017 Mes-Masson.xlsx" \
  "/project/6090753/active/ovcan_human/wes/mcgill_r004741_provenance/WES 2018 Mes-Masson.xlsx"

# provider BAM header (read groups and aligner)
module load StdEnv/2023 samtools/1.22.1
samtools view -H "/project/6090753/active/ovcan_gq/ammmasson_OvCAN_WES_bioinformatics_analysis_R004741/bam/OV2085_P64.bam" | grep -E '^@(RG|PG|CO)'
```
