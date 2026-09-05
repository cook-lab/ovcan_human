# Figure revision report — half B (Figs 4, 5, 6; Supp S5, S6)

**Scripts owned and edited:** `scripts/31_fig4_genomics.R`, `scripts/32_fig5_rare.R`,
`scripts/33_supp_genomics.R`, `scripts/36_fig6_adc.R`.
**Specifications for the manuscript writer:** `docs/manuscript/v2/FIGURE_SPECS_B.md`.
**Rendered and inspected in the isolated sandbox** (`OVCAN_PROJ=<sandbox>`); the real
`output/` was never read or written. Every one of the five figures was opened and
looked at after each edit round, not just re-rendered.

No script edited outside the four above. `00b_figure_theme.R` untouched.

---

## What changed, per figure

### Figure 4 — `31_fig4_genomics.R`

| Defect | Fix |
|---|---|
| `25914` / `493` hard-coded at 31:150-152; ATM/ATR percentages hard-coded at 31:162 | **All derived at render time** from the read-only MAFs in `judy_archive/data/wes - old/mutect2/` using 07's own rule (FILTER==PASS, popAF ≤ 0.001, coding non-synonymous). Stage 3 is `stopifnot`-checked against `wes_mutation_load.csv` per line and in total. Cost ≈ 1 s (7 columns × 557k rows). |
| "artefact collapse ATM 82%→9%" undefined | Defined on the panel: **% of the 22 WES lines with ≥1 coding non-syn. call in that gene**, with line counts (ATM 18→2, ATR 17→1). |
| Panel A was a single line (OV2295) presented as general | **Panel-wide series added beside it** (557,392 → 15,692 → 6,036), so the reader sees both. |
| **Panel B (5×5 ρ submatrix) could not support "rank 1 of 67, reciprocal-best"** | **Replaced** with a self-match-vs-best-non-self margin plot across **all 67 DepMap ovarian lines**: 66 non-self ticks per row, named best non-self, self-match, the Δρ segment, and z vs the non-self distribution. The BIN67↔COV434 SCCOHT-twin nuance is stated on the panel. The 5×5 matrix survives as `reports/assets/f_external_concordance.png` with the colour scale extended past 0.8 (it previously made 0.82 and 0.88 identical). |
| Panel C "coding cand." track had no axis | In-panel axis added (dotted references + labels 0/200/400), truncation at 460 made explicit with a double-slash break glyph on TOV21G and its value 1,416. |
| Panel C right-hand "patients" bars had no axis | Axis added (0/4/8/12/16) with the title `patients with ≥1 call (of 16)`. |
| Panel C "Variant class" legend = two identical greys keyed to nothing | Variant class is now a **real encoding on a second channel**: tier = full-cell neutral fill, variant class = centre bar in `variant_class_colours` (truncating/missense/in-frame/multi-hit, all four realised). Legend titles name the channel. |
| Tier colours colliding with the subtype strip | Verified the new neutral ramp took effect. Also found and fixed a **local family palette** (`31:69`) that reused four subtype hues; now the locked plum ramp, guarded by `check_palette_keys()`. |
| CDK12 frequency reads as somatic | `†` on every gene with **no Tier-1 call in any line** (ARID1B, BRAF, BRCA2, CDK12, CDKN2A, NF1, SMARCA2) plus an explicit line: *CDK12: 6 lines / 3 patients, every call Tier 2–3*. |
| TP53 row reads 12 vs an HGSC n of 11 | Denominator printed (`of 16` WES patients) and the composition spelled out in the spec: 11/11 HGSC patients + TOV112D (EC). |
| Panel D call rule invisible | Rule printed on the panel, read from the CSV: |log2| > 0.20 over > 50% of arm length, per-sample autosome probe-weighted median centring; acrocentric p-arms named as omitted; the sensitivity sweep cited. Six headline arms carry their values. |
| n / unit missing | Every panel states them (`22 WES lines / 16 patients`, `% HGSC patients (n=11)`, `n=66 per line`). |

### Figure 5 — `32_fig5_rare.R`

* **New Panel F: SWI/SNF.** The paper cited fig5 for the SWI/SNF reclassifications and
  no figure anywhere contained one. Panel F shows 7 lines × (3 genes × 2 assays), the
  rank of 31 printed in each cell, z-fill against the same 31 lines, and a WES-call
  column carrying the tier. It supports SMARCA4/SMARCA2 loss in both SCCOHT lines,
  TOV112D's Tier-1 `p.L639fs` (written correctly, not `p.L639X`), and — printed on the
  panel — **BIN67's SMARCA4 mRNA retained at rank 23/31 while protein is 2/31, margin
  0.105 log2**, i.e. why RNA-only authentication misses it. `‡` marks the Tier-3-only
  call that was withdrawn (deficient count 7 → 6).
* **New Panel D: refit with uncertainty.** MMR-d group relative exposure for all 22
  lines (TOV21G 0.733 vs 0.000–0.292), plus SBS6's bootstrap 95% interval converted to
  the same relative axis, selection frequency 99.5%/200 boots, and reconstruction
  cosine 0.977 (highest of 22).
* **Panel B labelled**: the line name and variant count are now in the axis title
  (`TOV21G, 2,417 exome-wide PASS SNVs`), with the total asserted against
  `wes_msi_mmr.csv`.
* **Panel A indel key**: limits rounded outward to 0.04–0.28 so the **endpoints render**
  (breaks placed exactly on limits are silently dropped by `guide_colourbar`); the
  hypermutator's 0.279 and the flag's rule (robust z > 5 AND > 3× median) are printed.
* **Panel C**: labelled as the **group MAXIMUM** with the winning signature named, the
  corrected top three (SBS6 0.877, SBS44 0.835, SBS15 0.809), and an explicit
  "screen, not attribution — SBS1/5/6/15/44 are mutually similar" note.
* **Panel E**: column headers are now neutral line names; the **z reference set is on
  the panel** (per gene across all 31 RNA / 31 protein lines, *not* across the 3
  columns); both assays shown with grey for not-measured; the CDX2-alone and
  literature-vs-measured SATB2 facts moved into the spec as measured values only. The
  GI/ovarian partition is derived from the recorded `CDX2 z > 1.00` convention, not
  from hard-coded line names.
* `1416` and `6.9×` (32:36) now read from `wes_mutation_load.csv`; the hypermutator is
  identified from `is_hypermutator`, never named by hand.

### Supplementary S5, S6 — `33_supp_genomics.R`

* **S5 made self-sufficient** so it can be cited: 23 rows asserted in code; bin size
  (10 Mb), segment resolution (median 185, range 88–481) and centring stated; a 50 Mb
  ruler row added inside each chromosome; the driver-label colour coding given its own
  **"Top labels"** legend (rust = gain locus, navy = loss locus); `‡` marks the CNV-only
  line (TOV3121D, whose sibling TOV3121EP *is* in Fig 4C); a per-line **FGA bar** added
  on the right, which makes the near-blank row quantitative — **TOV81D, FGA 0.021**,
  named on the panel. The diverging scale is symmetric about 0 and the legend now
  labels the clamp (`≤ -1.5` / `≥ 1.5`), which previously implied a ±2 range in which
  1.5 and 2.0 rendered identically. Local family palette replaced with the locked plum ramp.
* **S6 rebuilt around the probabilities.** Three sub-panels: stacked per-class
  posteriors (sum-to-1 asserted) with the call and its top probability; the
  top-vs-second margin with the 0.10 band and the **median 0.160** marked; the
  intrinsic-stratum strip. The five thin margins (0.022, 0.044, 0.072, 0.078, 0.088) and
  the two input-set-unstable calls (OV2085 IMR→PRO, TOV3041G MES→DIF) are both visible.
  "Intrinsic stratum" is defined on the panel. Local `call_pal` and `strat_pal`
  replaced with locked `consensusov_colours` / `stratum_colours`.

### Figure 6 — `36_fig6_adc.R`

* **Rebuilt as native ggplot** (was a bare `ComplexHeatmap` draw) so the two new panels
  compose.
* **Both blocks row-scaled** (z per feature across the 30 drawn lines) — the protein
  block is now legible instead of nearly uniform. Absolute information kept as a
  per-row **cross-line IQR in log2** column, so the reader can see what one z unit is
  worth. Clamped at ±3 with `squish` and the clamp stated (observed max |z| 4.94, from
  DPEP3 RNA, whose IQR is 0.00 — flagged on the panel as a single-line signal).
* **"not detected" has a legend entry** — a pale tile **plus a `×` glyph**. A pale tile
  alone was indistinguishable from z ≈ 0, which would have reproduced the original
  defect in a new form.
* Subtype strip separated by a gap and a hairline border; block gap widened so
  `MMMT (n=2)` and `SCCOHT (n=2)` no longer collide.
* **n per subtype block printed** (HGS 15, CC 6, EC 2, MC 3, MMMT 2, SCCOHT 2 = 30
  drawn). The different per-modality line sets (RNA 31 / protein 31; VOA6861 RNA-only,
  VOA14993 protein-only, both CC) are stated in the footnote.
* Protein unit defined: log2 TMT reporter intensity **relative to the Pooled Internal
  Standard (PIS, channel 1 of every plex)**.
* **New Panel B (FOLR1)**: per-line dots, RNA and protein, with BC 0.440 vs the 0.555
  threshold and mclust ΔBIC −1.06 → one component. Described as a wide range with a
  near-zero cluster of 4/15 — **not bimodal**.
* **New Panel C (compression)**: range / IQR / SD ratios for all 8 targets, sorted, with
  MSLN's 0.189 / 0.104 / 0.141 called out and the point made that the
  floor-insensitive statistics show *more* compression.

---

## Numbers that disagreed with the brief (reported, not forced)

1. **"median cross-line protein IQR is 0.34 log2".** Correct — but that is the
   **whole measured proteome** (8,357 rows × 31 lines: median IQR **0.344**). For the
   **8 ADC target proteins specifically** the median IQR is **0.728**. Both are stated
   in the figure footnote and the spec, because 0.34 alone would understate the spread
   of the proteins actually drawn while 0.73 alone would understate the general
   compression. Sources: `prot_dynamic_range.csv$prot_iqr` (per gene) and
   `prot_compression_floor_check.csv` (whole-proteome, floor-stratified).
2. **MSLN IQR ratio.** The brief gives 0.104; measured **0.10372**, which rounds to
   0.104 ✓ — but only when computed on the **30 paired lines** as
   `prot_dynamic_range.csv` does. Computing it from `adc_expression.csv` with 31 lines
   per modality gives 0.106. The figure uses `prot_dynamic_range.csv` (n_paired = 30).
   The brief cited `prot_compression_floor_check.csv` for these per-gene values; that
   file has no per-gene rows (2 rows, whole-proteome, stratified by whether the gene has
   an RNA zero), so the per-gene source is `prot_dynamic_range.csv`.
3. **figs5 FGA value depends on the column.** The near-blank row TOV81D reads
   **0.021** on `fga_auto_0.2` (autosome-restricted, which is what the figure plots) and
   0.073 on `fga_0.2` (all chromosomes). The figure and spec use the autosome column and
   name it.
4. **"the reclassified line" in the TP53 count** is **TOV112D**, an **endometrioid (EC)**
   line — the same line that carries the Tier-1 `SMARCA4 p.L639fs`. The brief's phrasing
   is consistent with the data; I mention it only so the writer does not read
   "reclassified" as "reclassified out of HGSC".
5. Every other number in the brief reproduced **exactly**: 25,914/493/170; ATM 82%→9%,
   ATR 77%→5%; DepMap margins 0.434 / 0.231 / 0.204 / 0.139 / 0.036 with z 11.5 and 4.0;
   CDK12 6 lines / 3 patients all Tier 2–3; TP53 12 patients; arm frequencies 3q 82%,
   20q 91%, 17p 82%, 8q 73%, 13q 64%, 19q 55% at n = 11; indel fraction 0.279; SBS6
   0.877 / SBS44 0.835 / SBS15 0.809 / SBS20 0.565; rel_mmr_d 0.733 vs 0.000–0.292;
   SBS6 selected 99.5%/200; reconstruction cosine 0.9765; SATB2 z +0.17 rank 19/31;
   SMARCA4 mRNA 23/31 vs protein 2/31, margin 0.105 log2; 6 deficient lines;
   DIF 5 / MES 5 / IMR 3 / PRO 2, IMR+MES 8/15, five margins < 0.10 with median 0.160,
   2 unstable; FOLR1 BC 0.440 vs 0.555 with ΔBIC −1.06; MSLN range 0.189 / SD 0.141.

---

## Palette additions needed in `00b_figure_theme.R` (not edited by me)

1. **`stratum_colours` is missing the key the data actually carries.**
   `hgs_heterogeneity.csv$cluster_label` is **`"Inflammatory/NF-kB-EMT"`**, but the
   locked palette key is `"Inflammatory"`. `check_palette_keys()` would fail on the raw
   label. `33_supp_genomics.R` therefore maps label → palette key through an explicit
   `STRAT_KEY` lookup (never redefining a colour), and asserts the mapping is total.
   **Recommendation:** add `"Inflammatory/NF-kB-EMT"` as an alias key pointing at
   `#C2410C`, exactly as `site_colours` already carries both display and metadata
   spellings. `37_supp_rnaprot.R` (other agent) has the same label and currently defines
   it locally at `37:220`, so a single alias fixes both call sites.
2. **No new colours are otherwise required.** `variant_class_colours` covers all four
   realised classes; `family_colours`, `subtype_colours`, `tier_colours`,
   `consensusov_colours` all had complete key coverage under `check_palette_keys()`.

---

## Not fixed, and why

1. **Residual FDR in the driver landscape.** Tumour-only WES without a matched normal
   cannot be cleaned by filtering; the figures now label tiers and flag non-Tier-1-only
   frequencies, but no figure change makes a Tier-2/3 call reliable. Orthogonal
   validation is the only fix.
2. **BIN67's SMARCA4 protein margin (0.105 log2) is genuinely thin.** The panel plots it
   and prints it rather than smoothing it. The SCCOHT protein calls are
   consistent-with, not proof-of, complete loss; IHC would settle it.
3. **SBS-96 exome renormalisation not applied** — the capture kit is not recoverable
   from the archive. Recorded in `wes_pipeline_parameters.csv`, stated in the spec, and
   flagged as a PI question, not something a figure can repair.
4. **Hartigan's dip test for FOLR1 is absent** because `diptest` is not installed in this
   R build. The bimodality coefficient and mclust ΔBIC both point the same way (not
   bimodal), so the conclusion does not hang on it, but the spec says the test was not
   computed rather than implying it was.
5. **Fig 4 is 9.1 in tall.** Four panels each now carrying its own axis and legend will
   not compress further without dropping an axis a referee asked for. If the journal
   caps figure height, panel D (CNV arms) is the natural candidate to move to
   supplementary, since S5 already carries the per-line profiles.
6. **`adc_subtype_summary.csv` top-mean subtypes for FOLR1 do not match the expected
   HGSC association** (`top_mean_matches = FALSE`, top subtype MC on n = 3). The figure
   shows the data and prints the per-subtype n; the recovery table is printed to the
   console for the writer. That is a result, not a figure defect.

---

## Incidental fix

`33_supp_genomics.R` and `36_fig6_adc.R` were emitting **512 and ~28**
`font family 'Arial' not found in PostScript font database` warnings per run (plus a
stray `Rplots.pdf`), because `ComplexHeatmap` measures text through `pdf(NULL)`, whose
font database has no Arial. `33` now aliases Arial to the Helvetica metrics in the
`pdfFonts()`/`postscriptFonts()` databases (the two share metrics, so measurement is
exact; the real exports go through `cairo_pdf`/`ragg`, which resolve Arial properly).
Warnings for that script: **0**. The remaining warnings in `33` are pre-existing
`S4Vectors:::anyMissing` notices from `GenomicRanges` built under R 4.5.3 while running
on 4.5.2 — not introduced here and harmless.

---

## Verification performed

* Each of `fig4.png`, `fig5.png`, `figs5.png`, `figs6.png`, `fig6.png` was read as an
  image after every edit round; all reported fixes were confirmed visually, including
  the ones that only became visible after re-rendering (clipped log-scale bars in
  Fig 4A, a pale "not detected" tile in Fig 6A that was indistinguishable from z ≈ 0,
  an outlier row stretching Fig 6A's z scale, inverted gene/assay headers in Fig 5F).
* `stopifnot` guards added where a number could silently desynchronise: waterfall stage 3
  vs `wes_mutation_load.csv` (per line and pooled), SBS-96 total vs `n_snv_used`,
  ConsensusOV probabilities summing to 1, 23 rows in figs5, 22 MAF lines, the mucinous
  CDX2 partition being 2 + 1, all variant classes mapped, and
  `check_palette_keys()` on every categorical fill.
* All four scripts render cleanly with `OVCAN_PROJ=<sandbox> Rscript scripts/<file>`
  and no errors.
