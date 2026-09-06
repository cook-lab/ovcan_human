# CNV coverage provenance and compatibility assessment — 6 September 2026

The new return closes the missing **per-normal target/antitarget coverage** request and independently reproduces the existing model coverage inputs. It also reveals a consequential limitation beyond the earlier “target-supported” caveat: four residual antitarget bins have technical extreme values that inflate archived segment means. Two small chromosome 1 segments are materially affected. Target-only resegmentation is warranted; this assessment does not replace any canonical CNV outputs.

## Evidence and reproduction

Run `python3 scripts/27_cnv_coverage_provenance.py` with the two returned bundles and original workstation CNV inputs restored. The standard-library script reads historical files but does not execute their commands. It writes only these proof tables and [the JSON assessment](provenance_assessment.json):

- [File identities](provenance_file_identity.csv): 62 returned files compared with prior files or recorded historical hashes.
- [Antitarget/segment contexts](provenance_antitarget_context.csv): 92 positive antitarget records, with 94 coordinate-overlap contexts.
- [Segment influence](provenance_segment_influence.csv): 71 overlapping segments, 69 of which actually contain antitarget contributors.
- [Bin-weight comparison](provenance_model_bin_weights.csv): 23 model summaries.
- [Archived centering baseline](provenance_archived_centering_baseline.csv): 23 original autosomal centers, pinned before correction with original canonical-table and CNS hashes. Subsequent canonical changes do not redefine this historical diagnostic.

The independent [coverage validator](coverage_validation.md) and [archive validator](archive_validation.json) provide complementary numerical and archive-integrity checks. The supplied `SHA256SUMS` has an incorrect self-entry; all 108 actual payload files pass its recorded hashes. That self-entry problem does not explain any coverage pattern.

## Input and command identity

Here `N` denotes `data/cluster_wes_retrieval/2026-09-06/ovcan_human_wes_cnv_coverage_2026-09-06/` and `P` denotes `data/cluster_wes_retrieval/2026-09-05/ovcan_human_wes_handoff_2026-09-05/`.

All **46 model CNNs** in `N` are byte-identical to every corresponding copy under `judy_archive/data/cnvkit wes - new/`. The **10 normal CNNs** match 230 filename/checksum records in `P/recovery/file_inventory.tsv`: each normal/bin-group file was inventoried in 23 model directories. Those earlier checksums were explicitly computed on recovered copies, not freshly on each original cluster path. The identity table preserves that distinction.

`N/normals/reference.cnn` is 26,875,054 bytes, SHA-256 `848c052274264ac544897977648860455ac742d89bc1de5edc7ece1ca185eeda`, matching all 35 archived copies. The target BED also matches all 35 archived copies, SHA-256 `2177970ffdbc933e488067f812e3ae4760310bd9a40fadf4e2fa5389480d955d`. These are established inputs, not new replacements.

`commands.txt`, `commands1.txt`, `command_diagram.txt` and `run_cnvkit.sh` are byte-identical to the prior return. They add **no successful-execution log or newly verified runtime version**. The named container remains CNVkit 0.9.10. The wrapper specifies array 1–27 but currently selects the six-line `commands1.txt`; this is a current recipe, not proof that 27 samples or 23 successful jobs ran. The earlier six malformed command lines are already resolved by the corrected text and need not be requested again. No new mapping, alignment header, duplication or contamination report is present.

## Normal coverage and model/reference compatibility

The five primary normal CNN pairs reproduce their earlier rounded target-depth means: 75.849–95.722×, median 92.398×. All 28 normal/model samples share the same target and antitarget coordinates. Each normal nevertheless has 31,560–48,010 zero-mean target bins (10.865–16.528% of bins). Zeros occur across all 22 autosomes; no wholly omitted chromosome is demonstrated. Across normals, 23,406 target bins are zero in all five and 62,536 in at least one. Crucially, **every one of the 204,706 reference-retained target bins has positive source coverage in at least three normals**. See [cross-normal support](../../output/wes_cnv_coverage_reference_support.csv).

These observations support incomplete or uneven capture coverage as a plausible technical explanation rather than a missing-path assertion. They do not establish GC as the cause or certify every alignment/preparation step. The reference mask already removes unsupported or unstable bins, and there is no evidence here to change that mask merely because original unfiltered normal coverage contains zeros. Tiny near-zero signed values in the robust pooled depth column are documented by the numerical validator; they are not physical negative read depth.

Normal antitarget depth means are 0.452–0.633×, whereas every model has zero depth in 42,705/42,709 antitarget bins. This is consistent with the independently documented interval-restricted model recalibration workflow. Same-kit compatibility is author-confirmed, and coordinate identity is verified; neither establishes equivalent upstream duplicate handling or filtering. The archived normal BWA recipe and sort/index logs do not demonstrate a duplicate-marking invocation. Absence of that record is not proof duplicate marking was omitted. CNVkit recommends compatible reference/sample preparation and duplicate flags for coverage counting. [CNVkit pipeline documentation](https://cnvkit.readthedocs.io/en/stable/pipeline.html).

## Technical antitarget artifacts and their segment influence

All four positive model antitarget coordinates overlap target bins. Coordinates below are zero-based, half-open:

| Antitarget interval | Overlapping target bins | Target-overlap bases |
| --- | ---: | ---: |
| chr1:145242279–145290970 | 2 | 265 |
| chr9:67727116–67844434 | 1 | 88 |
| chrX:135764982–135776495 | 44 | 11,513 |
| chrX:135782267–135793757 | 44 | 11,490 |

Thus none provides clean independent off-target validation. Across 23 models the 92 normalized CNR values range from **+12.4774 to +19.7489 log2**. Their weights are **0.690–0.979**, comparable to model target-bin median weights of 0.807–0.851. They cannot be dismissed as negligible-weight observations or interpreted directly as biological amplifications.

CNVkit normalizes target and antitarget coverage separately. Centering an antitarget group dominated by zero-depth sentinel values plausibly explains why its few nonzero observations become extreme; this mechanism is an inference from the observed input/output pattern and the implementation, not a newly recovered execution log. [CNVkit `fix` implementation](https://raw.githubusercontent.com/etal/cnvkit/v0.9.10/cnvlib/fix.py).

To test influence without re-running segmentation, script 27 partitions the ordered positive-depth CNR rows using each archived CNS probe count. The reconstructed bin-weighted means match **all 5,670 archived CNS means**, maximum absolute difference **4.998×10⁻⁶**, consistent with printed precision. This validates actual contributing-bin membership; simple interval overlap would incorrectly assign the chromosome 1 antitarget to two adjacent segments in OV2295 and TOV2835EP. Recorded CNS aggregate weights also include field-transfer effects and must not substitute for the actual contributing-bin weight sums.

The four antitargets contribute to 69 segments. Most effects are small, but the following two 1,889,317-bp spans require correction:

| Model / archived chr1 span | Probes | Archived raw log2 | Target-only mean, same contributing targets | Archived centered log2 | Target-only log2 at unchanged center |
| --- | ---: | ---: | ---: | ---: | ---: |
| OV2295, 143401653–145290970 | 74 | 0.514619 | 0.206097 | 0.460517 | 0.151995 |
| TOV2835EP, 143401653–145290970 | 74 | 0.569084 | 0.263667 | 0.419764 | 0.114347 |

Each antitarget supplies approximately 1.8% of its segment's contributing weight, enough to raise its mean by about 0.31 log2. Both spans switch from gain to neutral at the existing 0.20 threshold under this fixed-membership, fixed-center diagnostic. The next-largest affected autosomal amplitude difference is 0.05536 in TOV2929D chr9. These are diagnostic estimates, not replacement calls: omitting technical bins can also change breakpoint selection and subsequent centering. It is therefore inaccurate to claim all CNV summaries are unaffected simply because most segment amplitudes are moderate.

## Corrective scope and remaining requests

The appropriate bounded correction is to resegment **target CNR bins only across all 23 models**, preserving the original reference-masked target measurements, their normalized ratios/weights, model identities and the documented segmentation settings. Do not selectively repair only the two conspicuous segment means, normalize the four artifacts into apparent biological signal, or rerun WES alignment/coverage to close this issue. Keep original outputs and a comparison of FGA, arm summaries, hotspots and related manuscript numbers. Canonical integration is owned by the main analysis task; this report records the pre-correction evidence.

Per-normal CNN retrieval is **closed**. Remaining focused provenance requests are the successful manual CNV execution log/container identity and normal alignment headers or existing duplicate/mapping QC sufficient to establish preprocessing compatibility. The supplied antitarget BED construction history remains useful for explaining its target overlap. These are read-only recovery requests, not demands for new sequencing or expensive cluster jobs. They should not be conflated with the now-demonstrated local antitarget artifact, which can be corrected from available CNR data.

The 23 models, 16 independent patients and 6,194 coding candidates are unchanged by this audit. No mutation, cohort or normal-selection changes were made.

## Corrective execution completed

[Script 29](../../scripts/29_wes_cnv_target_only.py) subsequently replayed the original OV2295 and TOV2835EP cases with native CNVkit 0.9.10. All 202 boundaries/probe counts/log2 means in each case matched exactly. It then resegmented target CNR rows for all 23 models, producing 6,073 segments across all retained contigs, supported by 4,692,292 positive target bins. The source reference mask and target ratios/weights were preserved; all original inputs remain unchanged. See [the correction README](../../output/wes_cnv_target_only/README.md), [manifest](../../output/wes_cnv_target_only/manifest.csv) and [runtime record](../../output/wes_cnv_target_only/runtime.json). The prior fixed-membership assessment remains historical evidence, while the corrected segments provide the inputs for downstream FGA/arm/figure rebuilding.

## OV1369-R2: centering sensitivity after correction

Reproduce this bounded check with `python3 scripts/29c_cnv_centering_sensitivity.py`. [Script 29c](../../scripts/29c_cnv_centering_sensitivity.py) verifies original/current input hashes, recomputes both segment centers and the two bin-level medians using base R/matrixStats, checks the archived center against its pinned baseline, and writes the four diagnostic tables/JSON without changing canonical results. `OVCAN_DATA` supports an externally restored archive.

The independent span arithmetic reproduces corrected autosomal FGA **0.8591733609961716**. Its change from 0.6414685 is largely a relative-baseline effect. The segment/probe-weighted median changes from −0.1894744675 to −0.1178736413, a shift of +0.0716008262. Removing antitargets changes the native arm/gap partition and segmentation; it does not add new measured target ratios. Assessed segment span also decreases from 2,741.891 to 2,652.975 Mb because segment spans no longer bridge the same gaps. These span denominators are not the 63.514-Mb capture union.

| Segments used | Archived center | Corrected center |
| --- | ---: | ---: |
| Archived hybrid | 0.641469 | 0.827249 |
| Corrected target-only | 0.671549 | 0.859173 |

With corrected boundaries fixed, changing the center moves 515.184 Mb in 24 spans across the loss threshold. For example, chr4:49197114–190060527 has raw log2 −0.359791. Its centered value changes from −0.170317 to −0.241917, crossing −0.20 without any new biological measurement. See [the 2×2 check](ov1369_centering_2x2.csv), [threshold crossings](ov1369_centering_threshold_crossings.csv), and [source hashes/definitions](ov1369_centering_diagnostic.json).

A further check uses the same 195,544 positive-depth autosomal target CNR bins, excluding antitargets. Base R `median(log2)` is **−0.0137164500**; `matrixStats::weightedMedian(log2, weight)` is **−0.0201551582**. Applied to the corrected segments, these centers give FGA 0.911403 and 0.938344, respectively; using the existing CNR baseline with no additional centering gives 0.906069. The [alternative-center table](ov1369_alternative_centers.csv) records all definitions.

The bin-level alternatives remove dependence of the baseline on segment grouping, but remain affected by probe density, measurement noise and a highly altered genome. None establishes a biologically neutral state or absolute ploidy. The observed sensitivity warrants explicit qualification of precise FGA/state labels, not selecting a center to recover an earlier preferred value. Retaining the declared segment/probe-weighted convention with these sensitivity results is defensible for a descriptive relative-CNV resource. A change of convention would require a consistently defined comparison across all models, rather than an ad hoc adjustment of OV1369-R2 alone.

The [centering-sensitivity plot caption](ov1369_centering_sensitivity_caption.md) explains the archived/corrected curves and highlighted centers. [Script 39](../../scripts/39_cnv_centering_sensitivity.R) regenerates the local report PDF/PNG from existing profiles; [its curve data](ov1369_centering_curve.csv) are tracked independently of the rendered artifacts.
