# SBS3 sensitivity extension

This standalone analysis is implemented in `scripts/44_wes_sbs3_sensitivity.R`. It writes only this report directory and does not change the canonical mutation, copy-number or signature outputs. `scripts/44b_fig_sbs3_sensitivity.R` makes the companion diagnostic plot. No HRD class, clinical probability, tumour mutational burden or functional repair status is assigned.

## Inputs and screens

All 23 MAF/VCF pairs are matched to `output/wes_input_manifest.csv` and SHA256-verified. The starting candidates are MAF PASS single-base substitutions on chromosomes 1–22, X and Y, with maximum available `gnomADe_AF`, `AF`, `AA_AF` and `EA_AF` ≤0.001. This is **MAF PASS**, including converter-added filtering, rather than an assertion that all Mutect2 PASS variants enter. Sites are restricted to the same derived target union used by script22: 63,465,385 bp on standard chromosomes. Coordinates, alleles and hg38 reference bases are checked when creating SBS96 contexts. The rebuilt baseline count matrix exactly matches `output/wes_sbs_context_target_restricted.csv`.

Three nested screens are compared:

| Screen | Definition | Models | SNVs, range (median) | Total |
|---|---|---:|---:|---:|
| Baseline | Above canonical candidate screen | 23 | 303–2,417 (430) | 12,984 |
| Lower population-frequency threshold | Maximum available population AF ≤0.00001 | 23 | 221–2,129 (340) | 10,548 |
| Lower population-frequency threshold plus read support | Previous screen, AD reference + alternate ≥20, alternate AD ≥5, caller AF ≥0.05 | 23 | 127–1,668 (216) | 6,738 |

Population frequencies absent from all four fields remain eligible, consistently with the canonical filtering rule. This applies to 7,969/12,984 baseline candidates; absence from an annotation is not evidence of a measured zero population frequency. The stricter screen removes potential artifacts and more population-observed variants but does **not** establish somatic origin. There is no upper-VAF or heterozygous-VAF exclusion: either would remove plausible real clonal variants in pure/aneuploid cultures. The read-support settings are transparent exploratory thresholds, not newly validated caller cutoffs. They reduce counts and can change the biological mixture. Eight models have fewer than 200 variants after the combined screen; 200 is a descriptive count band, not a validated eligibility threshold.

## References and fitting

COSMIC v3.2 GRCh38 is the package's bundled 60-signature dictionary; the restricted dictionary retains the same 22 signatures specified in script22. These are two competing modelling assumptions, not an assertion that other processes cannot occur. Exact membership and pairwise similarity to SBS3 are in `reference_dictionaries.csv`.

For each dictionary, both its genome-reference spectrum and a target-opportunity-adjusted spectrum are fitted to the same target-restricted variants. Each reference entry is multiplied by the target/genome trinucleotide opportunity ratio and each signature column is then normalized to sum one. This transform exactly reproduces `output/wes_cosmic_target_normalized.csv`. The previously verified target opportunity table and BED are hashed in `source_manifest.csv`; no live capture or genome annotation is substituted. The target union is a capture-footprint proxy, not a sample-specific callable mask. It cannot correct coverage or variant-detection differences among individual models.

The point fit uses MutationalPatterns strict backward selection with `max_delta=0.004`. The procedure repeatedly removes small fitted components while controlling the deterioration in reconstruction. Dictionary restriction and bootstrap resampling are supported diagnostics of attribution stability; similarly shaped signatures can exchange contributions. [MutationalPatterns methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8845394/).

There are 276 point fits: 23 models × three screens × two dictionaries × two opportunity models. The baseline and combined stricter screens additionally receive 200 bootstrap resamples each, giving 184 bootstrapped fits and 36,800 resample fits. The intermediate AF-only screen supplies a point-fit sensitivity check. The model/substrate seed is held identical across all four reference conditions, so dictionary/opportunity contrasts use matched resampled mutation spectra. Never-selected signatures are explicitly zero-filled before summarization. Full seed, parameter, version and source details are saved in `run_manifest.json` and `session_info.txt`.

## Quantities and interpretation

`sbs3_sensitivity_summary.csv` records the fitted SBS3 count and its fraction of all fitted contributions, total observed and fitted counts, reconstruction cosine, and the number of selected components. For bootstrap fits, it includes the fraction selected, median fraction, and 2.5th/97.5th percentiles of fitted counts and fractions. The fraction uses each replicate's sum of fitted contributions; it is not approximated by dividing a count interval by observed N. Individual SBS3 replicates and paired reference differences are retained as compact CSVs.

These percentile ranges measure **conditional sampling variability in the retained candidates**. They do not incorporate uncertain germline origin, caller selection, reference-dictionary correctness, sample callability or biological stock history. Selection frequency is neither a posterior probability of HRD nor a diagnostic sensitivity/specificity estimate. A zero estimate or interval containing zero cannot establish HR proficiency. Bootstrap anticorrelations with other components are recorded to expose attribution trade-offs; they are not biological anticorrelations between pathways.

As an independent identifiability diagnostic, unpenalized non-negative least-squares fits are made with and without SBS3. The loss of squared reconstruction accuracy on exclusion is divided by observed N² for scale comparability. The nested fits are checked for nondecreasing squared residual error after SBS3 removal. A tiny loss means that other dictionary components approximate the spectrum well without SBS3. This is not a likelihood-ratio test, a p-value or a clinical decision rule. A high overall reconstruction cosine does not establish that any individual signature is correctly identified.

`patient_descriptive_sensitivity.csv` averages available model fractions within each of the 16 patient groups. Sublines remain separately visible, and no inferential patient-frequency comparison is performed. The 23 cultures are not 23 independent patients or independent validations of a process shared by one family.

## Why no new classifier is substituted

Methods such as SigMA can improve signature detection with sparse counts using cancer/assay-specific background models and calibrated classifiers. Their published performance does not transfer automatically to these tumour-only, filtered cell-line spectra. A separately justified application would need an appropriate trained model, assay/context compatibility and a defensible variant input, rather than treating strict-fit SBS3 exposure as a calibrated probability. [SigMA primary study](https://www.nature.com/articles/s41588-019-0390-2), [author software and parameter guidance](https://github.com/parklab/SigMA/wiki).

The next evidence gap remains suitable germline-inclusive allele counts for an allele-specific copy-number fit, followed by genomic-scar assessment only for accepted fits; see `reports/clinical_classification_2026-09-06/HRD_FEASIBILITY.md` and `docs/cluster/CLINICAL_CLASSIFICATION_NEXT_STEPS.md`. This extension provides a quantified description of signature uncertainty, not a replacement for those data or a RAD51 functional assay.
