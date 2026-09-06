# Current WES results and evidence

The current CNV profiles are the **target-only reconstruction** produced by script 29 and summarised by script 08. The earlier recovered-provenance tables describe historical inputs; they are not a pointer back to canonical archived CNS calls. Read the [6 September update](../reports/wes_cnv_coverage_2026-09-06/COVERAGE_UPDATE.md) and [current cluster follow-up](../docs/cluster/recovery/2026-09-06/FOLLOWUP.md).

| Files | Role |
| --- | --- |
| `wes_cnv_target_only/*.cns`, `manifest.csv`, `runtime.json` | Native CNVkit 0.9.10 target-only CBS outputs, input/output hashes, baseline replay and runtime evidence |
| `wes_cnv_segments.csv`, `wes_cnv_fga.csv`, `wes_cnv_arm_*` | Current centred relative segments, FGA, arm boundaries/calls/frequencies and threshold sensitivity |
| `wes_cnv_coverage_sample_summary.csv`, `wes_cnv_coverage_metric_summary.csv` | Verified primary CNN coverage for 23 models and five normals; bin means are not per-base threshold coverage |
| `wes_cnv_coverage_reference_support.csv`, `wes_cnv_coverage_chromosome_summary.csv` | Cross-normal support for reference-mask-passing targets and chromosome-level coverage patterns |
| `wes_cnv_coverage_sources.csv`, `wes_cnv_coverage_validation.json`, `wes_cnv_coverage_normal_comparison.csv` | New coverage source hashes, checks and reconciliation with the first handoff's rounded normal summaries |
| `wes_qc_model_summary.csv`, `wes_qc_*` | First-handoff model sequencing/duplicate/alignment/contamination/mosdepth QC; unchanged by the CNV correction |
| `wes_recovered_provenance*`, `wes_acquisition*` | Historical execution/input/variant-filter reconciliation and provider acquisition evidence |
| `wes_mutations_filtered.csv`, `wes_driver_tiers.csv`, `wes_mutation_load.csv` | Current variant candidates and summaries; the 6,194 coding-candidate set is unchanged |

The first handoff's normal-depth rows were initially marked unverified because their primary CNNs were absent. Their verification is now supplied by `wes_cnv_coverage_*`; the dated historical parser outputs remain intact. The original manual CNV execution log and full model/normal input-alignment provenance are still unavailable. Neither a compatible capture design nor the new coverage table alone proves identical preprocessing.

All copy-number measurements are relative. Target-only segment spans interpolate across gaps; their length does not equal directly observed exonic bases. FGA depends on the centring baseline and threshold, particularly for OV1369-R2. These records do not establish absolute ploidy, allele-specific copy number, LOH or HRD.

In the native CNR/CNS files, the column named `gene` preserves the source interval labels (often strings such as `chr1:14426-14627`). These labels are not a newly generated gene-level copy-number annotation.
