# Recovered WES QC validation

Validated all **23 baseline models from 16 patients** against the recovered primary reports and contemporaneous MultiQC records. No alignments were read, analyses rerun, or existing scientific result tables modified.

## Quantitative summary

All rows below have n=23. Fractions are on 0–1 scales. Full machine precision and per-model extrema are in `output/wes_qc_metric_summary.csv`; integer numerators remain in `output/wes_qc_model_summary.csv`.

| Metric | Minimum | Median | Maximum | Unit |
| --- | ---: | ---: | ---: | --- |
| `fastp_reads_after` | 108,382,950 | 123,188,938 | 136,781,554 | read ends (R1+R2) |
| `fastp_pairs_after` | 54,191,475 | 61,594,469 | 68,390,777 | paired fragments |
| `fastp_q30_fraction_after` | 0.881333074 | 0.919675632 | 0.930489356 | fraction |
| `fastp_read_retention_fraction` | 0.970112527 | 0.985629943 | 0.988750451 | fraction |
| `md_mapped_fraction` | 0.998307289 | 0.9994912 | 0.999840454 | fraction |
| `md_properly_paired_fraction` | 0.985724023 | 0.990897876 | 0.992481714 | fraction |
| `picard_duplication_fraction` | 0.149101 | 0.197735 | 0.239987 | fraction |
| `mosdepth_md_mean_target_depth_x` | 69.7151035 | 78.7567263 | 88.7499923 | x |
| `mosdepth_md_fraction_target_ge_10x` | 0.91 | 0.95 | 0.97 | fraction |
| `mosdepth_md_fraction_target_ge_20x` | 0.79 | 0.88 | 0.92 | fraction |
| `mosdepth_md_fraction_target_ge_30x` | 0.66 | 0.78 | 0.84 | fraction |
| `gatk_contamination_fraction` | 0 | 0.00266769106 | 0.00488455358 | fraction |
| `recal_primary_read_retention_fraction` | 0.681890997 | 0.738430476 | 0.769000874 | fraction |
| `md_mq0_fraction_of_mapped` | 0.122580977 | 0.131946742 | 0.159538016 | fraction |

## Denominators and source semantics

- fastp counts are R1+R2 read ends. Division by two is supported by equal first/last read-end counts in all 23 genome-wide samtools reports. Raw fastp summary/filtering/duplication fields agree between the MultiQC TSV and JSON. Original standalone fastp JSONs were not included in the handoff.
- Mapping and proper-pair fractions use genome-wide duplicate-marked (`md`) primary read ends, excluding secondary/supplementary alignments from that denominator. Picard duplication is a fraction despite the field name `PERCENT_DUPLICATION`; its numerator and denominator were independently rebuilt from paired/unpaired counts and reconciled to samtools.
- mosdepth 0.3.8 ran `--by intervals_sorted.bed` with no flag, mapping-quality or fast-mode override. Its [version-pinned implementation](https://github.com/brentp/mosdepth/blob/v0.3.8/mosdepth.nim#L219-L254) excludes FLAG 1796 (unmapped, secondary, QC-failed and duplicate alignments), uses MAPQ>=0 and corrects overlapping mates. Supplementary alignments are not excluded by 1796. This corrects the handoff prose claiming duplicate-inclusive depth.
- The [BED-mode distribution code](https://github.com/brentp/mosdepth/blob/v0.3.8/mosdepth.nim#L612-L627) accumulates per-base depth, rather than the means used for numeric windows. The md denominator is **63,709,951 summed BED bases across 242,421 intervals**, not the **63,514,049 bp union**; 195,902 bp are counted again because of overlaps. All 23 md summaries include the full summed length. Recal summaries omit 507–3,843 bp on no-read contigs; their sample-specific lengths are retained. Every stage denominator was reconciled to its included per-contig lengths.
- Means were reconstructed from integer depth sums divided by integer target lengths. The distribution files retain only two decimals: coverage fractions are approximate to 0.01, and exact threshold-covered base counts cannot be reconstructed. No spurious exact counts were calculated from those rounded fractions.
- Archived MultiQC per-contig values were reproduced as sums of the printed cumulative fractions excluding the zero-depth level. Because the depth distribution is rounded and omits unchanged high-depth levels, these plotted values differ from the authoritative integer-sum/length means. The latter are used in the new model summary and manuscript.
- Recalibrated CRAM QC describes alignments retained after interval-sharded ApplyBQSR. Its lower primary-read count is not a quality-filter failure. A target-read percentage calculated from those restricted alignments is not a genome-wide on-target fraction. Later report/CRAM rewrites are distinguished from the Oct 23 calling run; agreement with contemporaneous MultiQC establishes numerical QC consistency, not alignment byte identity.
- GATK contamination tables map uniquely through their repeated patient/sample IDs for all 23 models. The estimate and reported error field are retained without calling the error a confidence interval or inferring stock authenticity.

## Source checks and corrections

- Checked 197 present source-file hashes against the recovery inventory. All matched. `output/wes_qc_sources.csv` contains local SHA-256 values, cluster paths, inventory hashes and explicit missing-file status.
- Cross-checks: {'fastp_tsv_json_fields': 69, 'samtools_summary_fields': 897, 'mosdepth_coverage_points': 90, 'mosdepth_cdf_derived_contig_means': 874, 'picard_fields': 207}. Every compared value passed. The metric definitions, checks and source IDs are also recorded in `qc_validation.json`.
- The recovery generator searched the wrong Picard path, omitting all 23 model duplication records. This builder reads the actual nested sample directories and supplies the missing values.

## Five CNV reference exomes

Primary coverage CNNs independently validated: **0/5**. The recovery table reports length-weighted means of 75.85–95.72x (median 92.40x), using 63,514,049 bp in 290,475 CNVkit bins. These are a distinct method/denominator from model mosdepth. Missing primary CNNs prevent independently endorsing their reported arithmetic; missing paths and inventory hashes are preserved. No normal alignment, duplication or contamination QC was supplied.

## Actionable interpretation

TOV2929D is consistently the lowest-depth/lowest-coverage model (69.7151x mean; rounded fractions 0.91/0.79/0.66 at 10/20/30x). Locus-level depth should be checked before interpreting an absent variant in this model. These aggregate summaries alone do not justify exclusion or demonstrate a biological cause. All contamination estimates are below 0.005; the largest is OV3291. A reference-bundle mapping fraction includes MAPQ0 alignments (12.26–15.95% of mapped primary read ends), so high alignment rate is not synonymous with uniquely usable target coverage.

## Figure suggestion and reproduction

Use two compact panels: all 23 per-model mean target depths, and rounded cumulative target-coverage curves (or the 10/20/30x distributions), with a consistent model order. Mean-depth and coverage filters/denominators belong in the external legend. Source-ready curves for both stages are in `output/wes_qc_coverage_profile.csv`; the primary manuscript should use `stage=md`.

```sh
python3 scripts/24_wes_recovered_qc.py
```

The script uses only the Python standard library and validates source hashes, aliases, denominators, filtering accounting and cross-report agreement before writing its new outputs. It does not modify metadata, existing analysis outputs, figures, manuscript or release records.
