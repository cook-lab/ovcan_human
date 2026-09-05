# WES cluster retrieval checklist

Added 5 September 2026 after confirmation that the original cluster analysis is accessible. The audit's WES gaps refer to the local archive; cluster availability has not yet been checked. This checklist supports manuscript items M04/M05 and stronger genomic technical validation.

## First retrieval: run records, methods, and QC

| Priority | Files or records to retrieve | What they resolve |
|---|---|---|
| 1 | Sarek launch script, scheduler submission script and job log, exact pipeline release/commit, Nextflow version, input samplesheets, parameter JSON/YAML, custom `.config` files, and any local workflow modifications | Reconstruct the actual run, sample assignments, tumour-only configuration, and departures from pipeline defaults. |
| 1 | `pipeline_info/`, `.nextflow.log*`, execution trace/report/timeline, software-version files, and container image identifiers/digests | Link the final outputs to successful tasks, versions, and exact commands. Include the separate CNVkit and normal-reference runs, if they were launched independently. |
| 1 | `multiqc/` including its underlying data directory; individual sequencing/alignment/depth/duplicate QC reports; capture metrics and contamination estimates, if produced | Build per-model WES technical validation for all 23 models and the five public reference exomes. The exported numerical tables are needed alongside HTML reports. |
| 1 | Original vendor capture BED/design files and kit name/revision; the exact BED/interval list used for alignment QC and variant calling; CNVkit target/antitarget preparation scripts | Resolve capture provenance and distinguish vendor targets, padded calling intervals, and derived CNVkit bins. Record coordinate conventions and padding. Same-kit compatibility is already author-confirmed. |
| 1 | Sample/library/flowcell manifest, sequencing instrument and read length, read-group mapping, and any service-provider run reports | Complete acquisition methods and link model names/passages to the sequenced libraries and lanes. |
| 1 | Reference file manifest: genome FASTA with `.fai`/sequence dictionary and checksum; germline-resource and Mutect2 variant-panel-of-normals VCF identities/checksums; normal-sample inputs and pooled CNV-reference recipe | Establish the exact reference sequences and resources. The five normal exomes used for CNV normalization and the Mutect2 variant panel of normals are different resources and need separate provenance. |

Useful QC fields include read/fragment counts, mapped and properly paired fractions, duplicate fraction, insert-size summaries, on-target fraction, target-depth distribution, target bases above available coverage thresholds, and contamination estimates. Preserve metric definitions, target denominators, and the intervals used. A caller's contamination estimate is not an independent stock-authentication result. If a metric was not computed, record that; it can be evaluated later from the alignments.

The filenames above are search hints, not evidence that each module ran. The [Sarek output documentation](https://nf-co.re/sarek/docs/output) describes the usual `pipeline_info`, `multiqc`, `reports`, and preprocessing outputs. Use the documentation for the recovered pipeline version when interpreting the actual run; do not apply current defaults retrospectively.

## Variant-filtering and custom analysis files

Retrieve the original scripts/notebooks and configuration used to merge interval-level calls, run FilterMutectCalls, annotate variants, convert VCF to MAF, construct the normal reference, and run CNVkit. Include versioned source or a commit identifier, not just rendered notebook output. If the commands survive only in Nextflow `work/`, retain the relevant `.command.sh`, `.command.run`, `.command.log`, `.command.err`, `.exitcode`, and software-version records for those tasks, linked to their execution-trace entries.

The archived FilterMutectCalls headers explicitly name per-sample `*.mutect2.contamination.table`, `*.mutect2.segmentation.table`, and `*.mutect2.artifactprior.tar.gz`. These are particularly useful retrieval targets. Also retain available filtering statistics, pileup summaries, orientation-bias summaries, and final/unfiltered VCFs with indexes. Compare final VCF checksums with the local input manifest before transferring duplicates.

The `mutect2.segmentation.table` used during filtering is not, by its filename alone, a major/minor copy-number solution. Preserve its schema and generating command.

## Alignments and larger data: inventory first

For every model, locate the final recalibrated BAM or CRAM and its index. Record full path, file size, sample/read-group identifiers, checksum if already available, and the reference needed to decode CRAM. Include the normal-exome alignments used by CNVkit. If reads or alignments are reached through symlinks into `work/` or scratch storage, record the resolved target and confirm that it still exists.

These files would support missing coverage/QC calculations and targeted review of variants or copy-number calls. They can remain on the cluster for those computations; the first local transfer can be the small reports and manifests. Also inventory the original FASTQs, lane/sample mapping, and existing checksums for eventual data deposition.

If allele-specific copy-number, B-allele-frequency, purity/ploidy, MSI, or related analyses were already run, collect their outputs, input identities, commands, and fit diagnostics. Their existence and suitability need checking. Access to alignments makes additional analyses possible to assess; it does not by itself establish that a reliable HRD or MSI score can be produced. New HRD/MSI work is optional follow-up, separate from completing the current descriptor's WES provenance and QC.

## Path clues already recovered locally

- Archived CNVkit commands reference `/scratch/asmab/ovcan/`, including `intervals_sorted.bed`, `intervals.antitarget.bed`, and the `references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/` directory.
- Mutect2 headers also reference `/project/6090753/active/ovcan_gq_analysis/work/` for staged reference resources.
- All 23 CNVkit batch commands name a `*.recal.bam`; many paths begin `../work/`. Resolve those relative to the original command's working directory, not the local archive.
- The [path-hint table](../../../reports/audit_2026-09-05/wes_cluster_path_hints.csv) contains 322 recorded paths with their exact local provenance. These are historical hints, not verified cluster locations. The [model checklist](../../../reports/audit_2026-09-05/wes_cluster_models.csv) links the 23 models to their recorded CNV alignment inputs.

## Already resolved in the local audit

The TOV3121D conversion, all 23 final variant/CNV model records, GATK 4.5.0.0 and VEP 113.0 header versions, CNVkit 0.9.10 commands, five public-normal accession IDs, and the 290,475 derived target bins have been recovered. Cluster retrieval should confirm run provenance and supply missing upstream/QC records rather than treating these as absent results. The original vendor capture design and full pipeline release remain distinct from those recovered details.

## Recording the handoff

Keep each run's relative directory structure and supply a manifest with original cluster path, local destination, sample/run identifier, size, and checksum. Record actual retrieval time separately from the historical run date. A suitable local landing directory is `data/cluster_wes_retrieval/`, with separate `run_provenance`, `qc`, `references`, and `analysis_scripts` subdirectories; large alignment files can be represented by their cluster inventory until needed.
