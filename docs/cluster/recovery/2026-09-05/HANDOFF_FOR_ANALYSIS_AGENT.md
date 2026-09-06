# WES provenance handoff for the manuscript analysis agent

Written 2026-09-05 on Nibi after read-only recovery of the original cluster records.
Everything below is evidence-backed; sources are the files in this directory and the
bundle `ovcan_human_wes_handoff_2026-09-05.tar.gz` (layout at the end). Where something is
**not recorded**, say so in the manuscript rather than guessing.

## Cohort (unchanged)

42 models / 34 patients; WES = 23 models from 16 patients (11 HGSC patients). 22 further
samples with legacy GRCh37 outputs exist but were **deliberately excluded by the author**
(raw data lost; incompatible pipeline). Do not mention them in the manuscript.

## What the WES methods can now state

| Step | Established value | Source |
|---|---|---|
| Capture | Roche NimbleGen SeqCap EZ Exome v3 (hg19 design, 242,232 targets, 63.56 Mb; 368,146 probes) | `mcgill_r004741_provenance/SeqCap_EZ_Exome_v3.{targets,probes}.bed` |
| Sequencing | Illumina paired-end 2×100 bp; one lane per library (lane 5–8); delivered by McGill Genome Centre as GRCh37 BAMs (`R004741`). Instrument, flowcell IDs, run dates, library kit lot: **not recorded on the cluster** (check `WES 2014/2017/2018 Mes-Masson.xlsx`) | `scratch_archive_extract/samplesheet.csv`, `get_lane_info.sh`, BWA read groups |
| Re-processing | nf-core/sarek **v3.5.1**, Nextflow **24.10.2**, Apptainer; reads re-extracted from provider BAMs (samtools collate/fastq), split at 50 M reads | `sarek_3.5.1_run/launch/*`, `pipeline_info/params_2025-10-23_15-35-28.json` |
| Read QC/trim | fastp 0.23.4, `--detect_adapter_for_pe --disable_adapter_trimming --length_required 15` (no trimming) | fastp JSON `command` field |
| Alignment | bwa 0.7.18 (`-K 100000000 -Y -B 3`) to GATK GRCh38 bundle `Homo_sapiens_assembly38.fasta`; samtools 1.20/1.21 | `task_evidence/oct23_BWAMEM1_MEM_*` |
| Duplicates/BQSR | GATK 4.5.0.0 MarkDuplicates (marked, not removed); BaseRecalibrator/ApplyBQSR per interval shard with dbSNP 146, Mills & 1000G indels, `Homo_sapiens_assembly38.known_indels` | `task_evidence/oct23_GATK4_*` |
| Intervals | hg19 targets lifted with UCSC liftOver → 242,421 hg38 intervals (17 unmapped) → `intervals_sorted.bed` (63.71 Mb summed, 63.51 Mb merged) → 405 Sarek shards, no padding | `scratch_archive_extract/hglft_genome_*`, `intervals_sorted.bed` |
| Somatic calling | Mutect2 4.5.0.0 tumour-only, `--panel-of-normals 1000g_pon.hg38.vcf.gz --germline-resource af-only-gnomad.hg38.vcf.gz`, F1R2 collection; GetPileupSummaries (gnomAD AF-only) → CalculateContamination with tumour segmentation; LearnReadOrientationModel; FilterMutectCalls with GATK defaults (`OPTIMAL_F_SCORE`, FDR 0.05, max-events-in-region 2, max-alt-allele-count 1, min median base quality 20) | `task_evidence/oct23_MUTECT2_*`, `oct23_FILTERMUTECTCALLS_*`, VCF headers |
| Annotation | Ensembl VEP 113.0, GRCh38.p14 cache 113, `--everything --filter_common --per_gene --total_length --offline`; databases in header: dbSNP 156, ClinVar 202404, COSMIC 99, gnomADe/gnomADg v4.1, GENCODE 47, PolyPhen 2.2.3, SIFT 6.2.1 | `task_evidence/oct23_ENSEMBLVEP_VEP_*`, annotated VCF header |
| Copy number | CNVkit 0.9.10 `batch --method hybrid` (default), `--drop-low-coverage`, pooled reference from 5 public exomes SRR4039087/88/89/96/97 (PRJNA339046) aligned with bwa-mem2 to the same GRCh38 bundle; targets = autobin split of `intervals_sorted.bed` (290,475 bins), antitargets 42,709; default CBS segmentation; reference sex inferred (no `-y`) | `cnvkit_0.9.10_manual/commands.txt`, `normal_samples/*.sh|*.txt` |
| QC | MultiQC (2025-10-23) aggregating fastp, FastQC, MarkDuplicates, samtools stats, mosdepth (`--by intervals_sorted.bed`), bcftools/vcftools stats, VEP summaries | `sarek_3.5.1_run/multiqc/multiqc_data/`, `qc_metrics.tsv` |

Run identity to cite internally: run `angry_allen`, session `26cf85bf-0ac3-4237-8f48-c36b2ab20bd2`,
Slurm 3273686, completed 2025-10-23 23:48 EDT, 38,541 tasks. The 23 VEP-annotated VCFs from
that run are byte-identical (SHA-256 of the uncompressed stream) to the VCFs in
`output/wes_input_manifest.csv`, so the audited variant set is exactly the cluster output.

## Technical-validation numbers (23 models, from `qc_metrics.tsv`)

Use the `.md` (duplicate-marked, genome-wide) values for read-level denominators; use the
mosdepth `*_target` values for depth. Ranges:

- fastp reads after filtering: 108.4–136.8 M (R1+R2); Q30 ≥ 0.88; GC 0.42–0.47
- reads mapped: 99.83–99.98 %; properly paired ≥ 98.6 %
- MarkDuplicates duplication: 0.149–0.240
- mean target depth (mosdepth, md CRAM): 69.7–88.8×; target bases ≥30×: 0.66–0.84; ≥10×: 0.91–0.97
- GATK contamination estimates: all small (e.g. OV2085_P64 = 0.0039); see per-model rows
- normal exomes (derived from CNVkit `targetcoverage.cnn`): mean target depth 75.9–95.7×; **no** alignment QC exists for them
- do **not** report `cnvkit_percent_mapped_reads_in_target_bins` (61–68 %) as on-target rate: it was computed on interval-restricted recal BAMs

Variant-count caveat: vcftools PASS per sample (sum 19,816; TOV3121D_P68 = 547) counts VCF
records; the audited 16,081 / 389 count MAF rows after vcf2maf. Keep using the audited
cascade; if you cite the VCF-level counts, label the unit.

## Caveats to bake into the text

1. Alignments were re-derived from provider GRCh37 BAMs; **no FASTQ exists**. Any deposition
   of reads must start from those 23 BAMs (154 GB, `R004741/bam/`, no provider checksums).
2. The retained recalibrated CRAMs were overwritten by later, walltime-killed Sarek launches
   (10 dated 2025-10-24, 13 dated 2025-10-31). They are re-executions of the identical
   workflow, not the byte-identical inputs of the VCFs. Describe them as such if deposited.
3. Recalibrated CRAMs contain only reads overlapping the 405 interval shards (Sarek per-shard
   ApplyBQSR). Genome-wide statistics come from the `.md` CRAMs.
4. Not recorded anywhere: container digests, liftOver chain/version, bwa-mem2 version, the
   manual CNVkit Slurm logs, sequencing instrument/flowcell/dates. Write "not recorded".
5. The Mutect2 variant panel of normals (GATK 1000 Genomes PoN) and the five CNVkit
   reference exomes are different resources; keep them distinct in the text.
6. `GRCh37` labels in some archived MAFs are a vcf2maf header default; coordinates are GRCh38.

## Bundle layout (`ovcan_human_wes_handoff_2026-09-05.tar.gz`)

```
recovery/                       this directory (report, TSVs, generator)
sarek_3.5.1_run/launch/         run_sarek.sh, obcf-graham.cfg, nextflow_history.tsv, interval_list
sarek_3.5.1_run/pipeline_info/  params_*.json, execution_trace_*.txt, software versions (HTML reports omitted)
sarek_3.5.1_run/csv/            Sarek sample CSVs
sarek_3.5.1_run/multiqc/multiqc_data/   all MultiQC tables
sarek_3.5.1_run/task_evidence/  executed .command.sh/.err, exit codes, versions.yml per task
sarek_3.5.1_run/variant_calling/mutect2/<s>/   contamination, segmentation, filteringStats, stats (VCFs omitted: you have them)
sarek_3.5.1_run/reports/{markduplicates,samtools,mosdepth summary+dist,bcftools,vcftools}
cnvkit_0.9.10_manual/           commands*.txt, run_cnvkit.sh, normal_samples scripts+logs, <s>_new/{.cns,.call.cns,.bintest.cns}
scratch_archive_extract/        samplesheet*.csv, SeqCap targets, hglft_*.bed, ignored list, intervals_sorted*.bed, antitarget bed, get_lane_info.sh, run_test.sh, params_2025-10-07, obcf-graham.cfg (scratch version), reference .dict/.fai
mcgill_r004741_provenance/      README, PCGR/GenPipes scripts, SeqCap probes+targets, SureSelect V7 bed, WES 20xx spreadsheets, USB lists, GRCh37 CNA tables
```

Full 6.0 GB evidence set with SHA-256 manifest: `/project/6090753/active/ovcan_human/` on Nibi.
