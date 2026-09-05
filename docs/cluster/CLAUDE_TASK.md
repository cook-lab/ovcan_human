# Task prompt for Claude Code on the cluster

Copy the following prompt into Claude Code from the public repository checkout. No cluster root needs to be guessed in advance: the evidence table contains the historical roots, and the task begins by checking them.

---

We are preparing the OvCAN human ovarian cell-model resource for a Scientific Data Data Descriptor. The local audit recovered the final WES mutation and CNV inputs for all 23 models, but the original cluster project may contain the missing pipeline provenance, acquisition records, scripts and technical QC needed for the manuscript.

Your task is to recover **existing** WES evidence and inventory existing sequencing inputs. Read root `CLAUDE.md` and `AGENTS.md` if present, then `docs/cluster/WES_RECOVERY.md`. That document is the current scope and acceptance specification. Read these small inputs before searching:

- `metadata/samples.csv` and `metadata/line_family_map.csv`;
- `reports/audit_2026-09-05/wes_cluster_models.csv` and `wes_cluster_path_hints.csv`;
- `output/wes_input_manifest.csv`, `wes_vcf_header_provenance.csv`, `wes_pipeline_parameters.csv` and `wes_filter_cascade.csv`;
- `docs/cluster/evidence/cnvkit_commands.txt`, `cnvkit_command_diagram.txt` and their `manifest.csv`;
- `docs/data/archived_input_inventory.tsv` for the excluded source-data filenames/sizes (relative to `judy_archive/data/`), and `docs/REPRODUCIBILITY.md` for checkout expectations;
- `reports/audit_2026-09-05/genomics_audit.md` for resolved scientific issues.

The Git checkout intentionally omits raw archives and large data. A missing file in Git is not evidence that it is absent on the cluster. Historical roots include `/scratch/asmab/ovcan/` and `/project/6090753/active/ovcan_gq_analysis/`; verify accessibility and actual paths rather than assuming they still exist. Resolve relative `../work/` and normal-alignment paths against the original run's working directory. First inspect shallow run directories and exact recorded task hashes; avoid a whole-cluster scan.

Prioritize the H1–H4 items in WES_RECOVERY.md: actual workflow identity/version/commit and successful run configuration; model/library/passaged sample mapping; original capture design and reference provenance; existing numerical WES QC; and alignment/FASTQ/index inventory. Recover the M1–M3 items when already present: filtering/annotation support files and scripts, pooled CNV-reference details, and any previously completed allele-specific or MSI analyses with their diagnostics. Keep the Mutect2 variant panel of normals separate from the five public CNV-reference exomes.

Do not reopen resolved gaps. TOV3121D's VCF-to-MAF conversion is complete, all 23 models have mutation and CNV results, and the audited cohort has 6,194 coding candidates across 16 patients. GATK 4.5.0.0, VEP 113.0, CNVkit 0.9.10 and the five public normal accessions are already documented. The original capture-kit name/revision and vendor design remain to retrieve, but same-kit use is author-confirmed. The 290,475 derived CNVkit target bins are already recovered. GRCh37 labels in some raw MAFs are a conversion-header error; do not lift over the GRCh38 coordinates.

Pay particular attention to sample names and passages. OV/TOV, R/R2, D/EP/G and shared patient numbers do not make two samples interchangeable. RNA passage can differ from WES passage. The six corrected BAM-path rows in the handoff CSVs came from command lines containing a bare sample ID before the actual `.recal.bam` token; read their executed task logs before assuming those historical command lines were valid successful invocations. TOV3121D has a compressed annotated VCF, while most other archived VCFs are uncompressed with stale compressed-file index names. Preserve these distinctions.

Proceed autonomously with read-only discovery and copying selected small evidence into a separate staging directory. Read discovered scripts/configuration as text; do not source or execute them. Do not submit scheduler jobs, launch/resume Nextflow, rerun alignment/calling/annotation/QC, install software, initiate HRD/MSI analysis, scan entire alignments, hash the full large-data collection, or transfer large reads/alignments/VCFs. Existing checksums and lightweight file/header checks are sufficient for this inventory. Keep raw files and the audited `output/` tables unchanged. Do not commit raw genomic data, large files or credentials to Git.

Write the report and inventories specified in WES_RECOVERY.md under a new dated directory under `docs/cluster/recovery/`: `RECOVERY_REPORT.md`, `file_inventory.tsv`, `model_status.tsv`, `qc_metrics.tsv`, and `run_reference_manifest.tsv` or JSON. Account for all 23 baseline models and all five public-reference exomes, and list any additional models separately. Every claim needs its source file/path and observed status. Preserve units, target denominators, original paths, run identity and available checksums; separate historical run date from actual discovery/retrieval time. Do not infer provenance from mtime or fill missing metrics with zero.

Finish with a concise account of what was recovered, what manuscript details are now resolved, which precise gaps remain and why, and whether existing evidence supports a specific later computation. “Not found in the searched roots” is a valid documented outcome. A proposed later job should identify inputs and outputs, but this task does not authorize launching it. Continue useful independent recovery work if a directory is missing or a particular item needs laboratory confirmation; do not stop at the first unavailable path.

---
