# SBS3 is strongly dependent on the input screen and reference model

The extension completed **276 point fits and 36,800 bootstrap refits** across **23 models from 16 patient groups**. No model has consistently nonzero SBS3 support across the target-adjusted full and restricted dictionaries under both candidate screens. This does not identify any model as HR-proficient: it shows that these tumour-only exome spectra do not support a robust categorical SBS3/HRD assignment with the tested approach.

The new tables extend the existing analysis with target-adjusted bootstrap uncertainty, matched-resample dictionary/opportunity contrasts, stricter candidate screens, per-patient descriptive summaries and a reconstruction test excluding SBS3. All 46 existing genome-reference SBS3 point estimates were independently reproduced to floating-point precision; original scientific outputs remain unchanged.

## What changes when the assumptions change

Each row below represents the same 23 models. “Nonzero estimate” is an algorithmic contribution, not an HRD-positive designation. “Lower bound >0” refers to the 2.5th percentile of conditional bootstrap fractions, not a clinically validated decision threshold.

| Candidate screen | Opportunity model | Dictionary | Models with nonzero SBS3 point estimate | Models with bootstrap lower bound >0 |
|---|---|---|---:|---:|
| Baseline | Genome | Full, 60 signatures | 7 | 0 |
| Baseline | Genome | Restricted, 22 signatures | 21 | 12 |
| Baseline | Target-adjusted | Full, 60 signatures | 2 | 0 |
| Baseline | Target-adjusted | Restricted, 22 signatures | 6 | 0 |
| Stricter | Genome | Full, 60 signatures | 10 | 0 |
| Stricter | Genome | Restricted, 22 signatures | 21 | 11 |
| Stricter | Target-adjusted | Full, 60 signatures | 6 | 0 |
| Stricter | Target-adjusted | Restricted, 22 signatures | 12 | 1 |

Machine-readable values are in [dictionary_summary.csv](dictionary_summary.csv). The intermediate population-frequency-only screen and every model-specific result are included in [sbs3_sensitivity_summary.csv](sbs3_sensitivity_summary.csv).

The effect is substantial, rather than a small rescaling of the previous estimates. In **OV1369-R2**, the baseline full-dictionary point estimate changes from **37.7% using the genome reference to 0% after capture adjustment**; the target-adjusted bootstrap range is still **0–52.2%**. Its restricted-dictionary target-adjusted estimate is **18.4%**, with range **0–38.2%**. The data do not distinguish a unique stable allocation to SBS3.

**TOV1369** is the only model with a nonzero target-adjusted bootstrap lower bound under the stricter screen and restricted dictionary: **61.1%, range 23.0–69.9%**, selected in 200/200 resamples. With the full dictionary on the same candidates, its point estimate is **10.2%, range 0–47.1%**, selected in 122/200 resamples. This is a model-dependent lead for further evaluation, not a robust classifier result. TOV1369 and OV1369-R2 belong to the same patient family.

**TOV81D**, which carries the separately curated pathogenic BRCA2 deletion, has a target-adjusted baseline SBS3 point estimate of zero in both dictionaries. The stricter screen produces estimates of **14.8%** (full dictionary, range **0–38.6%**) and **14.5%** (restricted, range **0–46.1%**), using only **171 SNVs**, compared with 608 at baseline. Neither the mutation nor these unstable fits establish biallelic BRCA2 loss or present HR function. The exact allele and prior wild-type transcript evidence remain in the [HRD feasibility review](../../clinical_classification_2026-09-06/HRD_FEASIBILITY.md).

## Why apparently precise fits can be misleading

With the target-adjusted full dictionary, removing SBS3 changes the best unpenalized reconstruction cosine by at most **0.000427** at baseline or **0.001032** under the stricter screen. Other components can approximate the spectra closely without SBS3. Conversely, restricting the dictionary can force more contribution onto the remaining components. Neither dictionary is established as the biological truth; their disagreement demonstrates dependence on the assumed alternatives.

Bootstrap contribution trade-offs make this ambiguity visible. In the baseline target-adjusted full fit for OV1369-R2, SBS3 count contributions correlate negatively with SBS39 (**r = −0.659**) and SBS8 (**r = −0.624**) across resamples. Similar trade-offs occur in TOV1369. These are fitting relationships, not evidence that one biological repair mechanism suppresses another. Full component diagnostics are in [sbs3_bootstrap_competitors.csv](sbs3_bootstrap_competitors.csv).

The stricter screen retains **6,738/12,984 SNVs** and leaves eight models below 200 SNVs. It uses a lower maximum available population frequency together with minimum read support and VAF; it is not a matched-normal somatic filter. **7,969 baseline candidates lack all four population-frequency annotations** and remain eligible under the documented rule. New bootstrap precision cannot remove uncertainty caused by residual germline variants, ascertainment or the lack of model-specific callable masks.

## Deliverables and reproduction

- [Methods and limitations](METHODS.md), [external figure legend](FIGURE_LEGEND.md).
- [Per-model estimates and uncertainty](sbs3_sensitivity_summary.csv), [candidate burdens](substrate_burden.csv), [dictionary summary](dictionary_summary.csv).
- [Matched-resample reference contrasts](paired_reference_sensitivity.csv), [patient descriptive means](patient_descriptive_sensitivity.csv), [all component point fits](all_signature_point_exposures.csv).
- [Validation](validation.json), [source hashes](source_manifest.csv), [runtime manifest](run_manifest.json).

From the repository root, with private inputs restored under their archived layout or selected by `OVCAN_DATA`:

```bash
Rscript --vanilla scripts/44_wes_sbs3_sensitivity.R
python3 reports/molecular_extension_2026-09-06/signatures/validate_sbs3.py
Rscript --vanilla scripts/44b_fig_sbs3_sensitivity.R
```

The analysis defaults to 200 resamples and three local workers; `OVCAN_SBS3_BOOTS` and `OVCAN_SBS3_WORKERS` can override these explicitly. Matching per-job caches are reusable and intentionally untracked. The reviewed figure PDF, source, numeric tables and external legend are included in the repository; the PNG remains a reconstructable local preview. No package installation or cluster job submission is performed by these scripts.

Validation checks all model/reference combinations, bootstrap percentiles and denominators, source/script hashes, component fraction sums, the nested least-squares residual inequality within numerical tolerance, all previous genome-reference point estimates, and patient-level averaging. The maximum absolute difference from previous point estimates is **4.44 × 10⁻¹⁶**. Repeated input-stage warnings were checked in a first-model probe and concern deprecated `S4Vectors:::anyMissing()` calls in the installed genomic-range dependencies; the rebuilt input matrix and transformed reference both match their canonical counterparts exactly. The final PDF was rendered with Poppler and inspected: equal-width comparison panels, matching axes, readable model labels/counts, visible zero estimates and no clipped interval endpoints. The plot uses Arial and the shared rust/slate theme.

The next useful data addition remains germline-inclusive SNP/allelic counts and a defensible allele-specific copy-number fit, followed by scar metrics only if fit quality and ploidy ambiguity permit. Orthogonal current-stock HR assays would address function. The [cluster next-step request](../../../docs/cluster/CLINICAL_CLASSIFICATION_NEXT_STEPS.md) remains the place to authorize and document that work; these diagnostic fractions should not be entered as clinical HRD labels.
