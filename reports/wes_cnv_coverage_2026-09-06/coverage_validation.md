# CNV coverage validation — 6 September 2026

Independently parsed all **5 normal and 23 model target/antitarget pairs** (56 primary CNNs). All 58 primary coverage/design/reference source hashes match the supplied archive manifest. The pooled reference is byte-identical to the previously audited reference. Historical September 5 QC, CNV calls and manuscript/figure files were not rewritten.

## Target depth and zero-bin summaries

| Cohort | n | Length-weighted mean depth, min / median / max (×) | Zero-depth target bins, min / median / max |
| --- | ---: | --- | --- |
| cnv_reference_exomes | 5 | 75.849018 / 92.397735 / 95.722376 | 31,560 / 41,795 / 48,010 |
| baseline_models | 23 | 78.644238 / 89.867763 / 98.667319 | 1,420 / 1,517 / 1,695 |

Every target file has **290,475 bins spanning 63,514,049 bp**, with no within-design overlaps. Every antitarget file has **42,709 bins spanning 2,580,384,051 bp**, also without within-design overlaps. All 28 samples have identical ordered coordinates and gene labels within each bin group; the target CNNs also match the supplied target BED. The target and antitarget designs overlap by **23,356 bp** and must not be described as disjoint. The pooled reference contains precisely their combined bin coordinate set.

The denominator differs from the 63,709,951 summed bases of the overlapping mosdepth calling intervals. The CNVkit target length is the 63,514,049 bp union after merging/splitting. Directly pooling CNN and mosdepth metrics would conflate interval definitions, alignment stages and coverage implementations.

## Recovered normal means

| Reference exome | Verified length-weighted mean (×) | Prior reported mean (×) |
| --- | ---: | ---: |
| SRR4039087 | 95.467332221 | 95.47 |
| SRR4039088 | 83.750941443 | 83.75 |
| SRR4039089 | 75.849018118 | 75.85 |
| SRR4039096 | 92.397735015 | 92.40 |
| SRR4039097 | 95.722376204 | 95.72 |

All five values agree at the original two-decimal precision. They are now independently validated from primary CNNs; the earlier dated report correctly retains its historical reported-only status.

## Normal target coverage gaps

| Normal | Zero-depth bins (%) | bp in zero-depth bins (%) | Positive bins after reference mask / 204,706 |
| --- | ---: | ---: | ---: |
| SRR4039087 | 41,795 (14.389%) | 9,475,321 (14.918%) | 202,822 |
| SRR4039088 | 48,010 (16.528%) | 10,867,183 (17.110%) | 200,349 |
| SRR4039089 | 31,560 (10.865%) | 7,161,334 (11.275%) | 204,418 |
| SRR4039096 | 40,929 (14.090%) | 9,283,199 (14.616%) | 203,149 |
| SRR4039097 | 46,259 (15.925%) | 10,477,068 (16.496%) | 200,604 |

All five normals have zero-depth target bins on every autosome. Within each normal, chr19 has the largest zero-bin fraction (24.43–37.13% across normals), and chr13 the smallest (5.57–7.55%). Thus the gaps are distributed unevenly rather than confined to a wholly missing chromosome. CNN summaries do not identify their cause. Per-chromosome means, counts and base footprints are recorded in `wes_cnv_coverage_chromosome_summary.csv`.

## Reference support and processing compatibility

- Target: 23,406 bins have zero depth in all five normals; 62,536 have zero depth in at least one. Of 204,706 reference-mask-passing bins, 8,833 lack positive depth in at least one normal. Detailed counts and base lengths by number of supporting normals are in `wes_cnv_coverage_reference_support.csv`.
- Antitarget: 1,016 bins have zero depth in all five normals; 1,235 have zero depth in at least one. Of 40,923 reference-mask-passing bins, 131 lack positive depth in at least one normal. Detailed counts and base lengths by number of supporting normals are in `wes_cnv_coverage_reference_support.csv`.

Normal antitarget mean depths range from 0.451933 to 0.632583× (median 0.587532×). Each model instead has exactly 42,705 zero-depth antitarget bins and only four positive bins. The primary files confirm that the usable CNV signal remains predominantly on-target. Coordinate agreement confirms a common bin design; it does not by itself establish equivalent upstream processing or justify normalising away the asymmetric off-target support.

All retained target bins have positive depth in at least three normals: 3,355 have support from three, 5,478 from four, and 195,873 from all five. The mask excludes every bin with no positive depth in any normal. Its 204,706 retained targets span 44,729,346 bp. This is a support summary, not an independent calibration of CNV accuracy.

The pooled-reference mask retains 204,706 target and 40,923 antitarget bins. Reconstruction uses the previously validated union of |log2|>5, spread>1, depth=0, and GC outside 0.3–0.7. All 184 model bin-count/support comparisons with the September 5 audit pass. Mask criteria overlap and are not additive. This validation does not refilter or regenerate CNV profiles.

The pooled reference contains 25 tiny negative depth values (minimum −2.77556×10⁻¹⁷), consistent with numerical cancellation in robust aggregation. Two target bins pass the literal historical depth==0 mask: chr17:63963626–63963868 and chr20:2396442–2396684 (BED coordinates). Their coordinates and values are recorded in the JSON audit. These are not physically negative read depths. The historical mask is preserved exactly; no rounded-zero rule or additional bin exclusion has been applied.

## Units, precision and limits

CNVkit [coverage implementation v0.9.10](https://raw.githubusercontent.com/etal/cnvkit/v0.9.10/cnvlib/coverage.py) defines depth as mean pileup per bin and computes log2 from that mean. For all positive-depth input rows, the maximum discrepancy between log2(depth) and the stored log2 is 0.00005683, consistent with text rounding; every zero-depth row uses the -20 sentinel. The [reference implementation](https://raw.githubusercontent.com/etal/cnvkit/v0.9.10/cnvlib/reference.py) separately robustly aggregates absolute depths and bias-corrected log2 values, so the reference log2/depth columns are not interchangeable.

- CNN depth columns contain rounded bin means. Weighted totals computed from them are not exact original integer read-base counts.
- A positive mean bin does not establish coverage at every base. Fractions of bases at >=10x, >=20x or >=30x cannot be reconstructed from these CNNs; none are calculated.
- Base length in zero-mean bins is exact from the interval coordinates, but is only a lower bound on all uncovered bases because positive-mean bins may contain uncovered positions.
- The reference depth column is a robust aggregate across normals, not the arithmetic mean or a sixth independent sample. Corrected reference log2 values are not log2 of its absolute depth column.
- Model CNNs describe interval-restricted recalibrated inputs from the recorded manual CNVkit workflow. Normal CNNs describe separately processed reference exomes. Matching interval coordinates does not establish matching read-processing filters or genome-wide coverage.
- All 23 model antitarget files have only four positive-depth bins. The target and antitarget designs overlap at a small set of bases; positive antitarget signals are not necessarily independent off-target support. The broad normal antitarget coverage does not supply missing model off-target measurements; CNV profiles remain predominantly target-supported.
- No normal alignment/mapping, duplication, contamination, read-level integrity or per-base coverage reports were supplied in this archive. Those checks remain unavailable, rather than failed.
- The reference mask is reproduced for validation only; no bin thresholds, normal selection, CNV segmentation, copy-number interpretation or historical scientific outputs are changed.

## Reproduction and outputs

```sh
python3 scripts/26_wes_cnv_coverage_qc.py
```

`output/wes_cnv_coverage_sample_summary.csv` contains 56 sample/bin-group rows. The metric summary gives n, min, median, max, units and source IDs; the normal comparison records prior rounding agreement; the reference support table records cross-normal support; the source table records primary checksums. Machine-readable validation is saved both beside this report and under output/. Source inputs remain read-only.
