# fig_report_A.md — figure revision report (Agent A)

**Scope:** `scripts/30_fig1_overview.R` (Fig. 1), `scripts/34_fig2_qc.R` (Fig. 2),
`scripts/35_fig3_biology.R` (Fig. 3), `scripts/37_supp_rnaprot.R` (Supp. S1, S2, S3, S4, S7,
S8). Nine figures. Date 2026-07-24. R 4.5.2.

**How this was rendered.** Everything was rendered through an isolated sandbox
(`OVCAN_PROJ=<sandbox>`) holding a frozen copy of `output/`; the real `output/` was never read
or written. Script edits went to the real `scripts/`. See §6 for the isolation audit.

**Two directories in the real project are NOT mine and must not be read as stale or as my
output:**

- **`figs/`** currently holds diagnostic figures from the clean-room pipeline run that was
  executing while these figures were revised. Nothing in it came from my four scripts (they
  write only to `docs/manuscript/figures/` and `reports/assets/`).
- **`reports/assets/`** is shared. Four asset filenames are written by **both** a pipeline
  analysis script and one of my figure scripts, so whichever runs last wins and the HTML
  report will show that version:

  | asset | also written by | my script |
  |---|---|---|
  | `f_prot_compression.png` | `19_proteomics_dynamic_range.R` | `34_fig2_qc.R` |
  | `f_variance_partition.png` | `17_variance_confounders.R` | `35_fig3_biology.R` |
  | `f_passage_check.png` | `17_variance_confounders.R` | `37_supp_rnaprot.R` |
  | `f_passage_check_crossassay.png` | `17_variance_confounders.R` | `37_supp_rnaprot.R` |

  This collision predates my changes, but my versions of all four now differ substantially
  from the analysis-script versions (the passage panel gained a third model; the compression
  panel gained the criterion definition; the variance panel gained means and the n_features
  fix). **The final re-render must run 34/35/37 after 17/19**, or the report will silently
  display the un-revised panels. Renaming the figure-script outputs (e.g. a `f_ms_` prefix)
  would remove the ordering dependency permanently. All four scripts run clean end to end
with no errors. Every PNG was viewed at render size before being called done.

The only remaining warning is pre-existing and benign: `37_supp_rnaprot.R` emits
`grid.Call(C_textBounds, ...): font family 'Arial' not found in PostScript font database`
while `ComplexHeatmap` measures text for the Supp. S3 layout. It was present in the baseline
run, the rendered `figs3.pdf`/`.png` embed Arial correctly (cairo_pdf and ragg both resolve
it), and the annotation-strip proportions are correct in the output. I left it alone rather
than change the S3 device path and risk shifting a layout that now renders correctly.

**Panel-letter changes** (three figures gained panels; the specs file documents each):

| figure | before | after |
|---|---|---|
| Fig. 2 | A–D | **A–F** (added a library-depth panel and a CV-by-abundance panel; the bridge panel moved to F) |
| Fig. 3 | A–E | **A–F** (added a per-line marker panel F) |
| Supp. S1 | A–B | **A–C** (added a silhouette panel C) |
| Supp. S8 | A–B | **A–C** (proliferation moved out of A's size channel into its own panel B; PROGENy is now C) |

---

## 1. Numbers that disagreed with what the brief told me to expect

These are the important ones. Nothing in the brief was badly wrong, but four items differ and
one of them changes a sentence the paper might otherwise write.

1. **Supp. S3 correlation range.** The brief said "state the true range (~0.75–0.90)". The
   actual off-diagonal Spearman range is **0.757–0.969** (median 0.834). The high end is
   0.969, not ~0.90 — it comes from the tightest same-family/adjacent pairs. The colour scale
   now spans the observed range exactly and the footnote quotes 0.757–0.969. **If the
   manuscript text says "~0.75–0.90" it is wrong.**

2. **Fig. 2A depth difference of "~14 M fragments".** Both numbers exist and they are
   different quantities: **assigned gene counts** differ by **+13.9 M** (54.4 → 68.3 M,
   Wilcoxon p = 0.77) and **processed fragments** differ by **+17.7 M** (60.0 → 77.8 M,
   p = 0.37). The "~14 M" in the brief is the assigned-gene-counts figure. Both are printed
   in the panel with both p values, because **neither per-site depth difference is
   significant** — depth is a plausible mechanism, not a demonstrated one, and the panel must
   not be written as if it were.

3. **Supp. S7 widest bootstrap CI.** The brief quoted KRT20's *d* CI as [1.72, 18.97], which
   is correct. But **MUC5AC's is wider still: [−0.44, 20.64]** — and unlike KRT20 it crosses
   zero. If the paper wants to make the "CIs are the point" argument, MUC5AC is the stronger
   example: a large point estimate (d = 1.92) with an interval spanning no effect to an
   effect of 20. Three CIs exceed the ±6 axis clip (CDX2 → 9.5, KRT20 → 19.0, MUC5AC → 20.6)
   and all three are printed as labelled arrows.

4. **Supp. S1 protein CC silhouette.** The brief said **−0.003**; `silhouette_by_modality.csv`
   stores **−0.0035**, which the figure prints as **−0.004**. Trivial rounding, flagged so the
   text and the figure do not differ in the third decimal.

5. **Supp. S8 PROGENy agreement count depends on the threshold, and I had to choose one.**
   The brief said "only 2 of 5 hypoxic-stratum lines show elevated PROGENy Hypoxia". That is
   true at **z > +1** (OV3331, OV4453). At z > +0.5 it is 3 of 5 (TOV3291G joins). I fixed the
   threshold at **z > +1**, stated it in the panel, and added three threshold-free statistics
   so the claim does not rest on the cut: **Spearman ρ = 0.63 (p = 0.014, n = 15)**, and the
   inflammatory stratum's **median** PROGENy Hypoxia (+0.64) is **higher** than the hypoxic
   stratum's (+0.56). The honest summary is *partial, axis-dependent* agreement — the NF-κB
   axis agrees well (ρ = 0.82, p = 0.0003), the hypoxia axis only moderately.

6. **Fig. 3 panel E "MKI67 boxed under SCCOHT".** MKI67 was **never boxed** — the existing
   code already filtered the proliferation control out of `expected_cell`. What the reviewers
   most likely saw was the group-separator rule meeting the panel edge next to MKI67's rust
   SCCOHT cell. No bug to fix, but the ambiguity was real, so MKI67's axis label is now
   `MKI67 (control)` and the panel states that it has no expected cell and is not graded.

7. **Fig. 1 WES-SNV marks were already correct.** The brief asked me to check whether SNV was
   marked for both TOV3121EP and TOV3121D. It was not: the script already read `has_wes_maf`,
   and `metadata/line_family_map.csv` correctly gives TOV3121D `has_wes_cnv = TRUE,
   has_wes_maf = FALSE`. I added a `stopifnot` asserting that TOV3121D is the **only** line
   where the two columns differ, so a future metadata edit cannot silently break it. The 23 /
   22 split is now printed under the exome box in panel A and beneath the columns in panel B.

8. **Supp. S7 EC markers were not silently dropping out.** Adding `"EC"` to `grp_lvl` was
   necessary and is done, but the failure mode was milder than the brief expected: `grp` is
   used only in `arrange()`, so with `grp = NA` the three EC markers sorted to the *end* of
   the axis (rendering above SCCOHT, outside their block) rather than disappearing. They were
   present but misplaced and unlabelled. The script now asserts `!any(is.na(me$grp))`.

---

## 2. What changed, per figure

### Fig. 1 — `30_fig1_overview.R`
- **Every count is now read from `metadata/line_family_map.csv`**, and the script asserts
  42 lines / 34 patients / 13 all-three models, that the assay column totals equal the panel A
  box numbers, and that TOV3121D is the sole CNV-without-MAF line. Panel A's box labels and
  the per-centre sub-label (`CHUM 29 · BC Cancer 12 · OHRI 1`) are derived, not typed.
- **Panel B now has all four legends**, including the previously missing **subtype** key. The
  subtype block labels on the right are additionally drawn *in their subtype colour*, so the
  strip is keyed twice.
- **Family strip moved to the locked plum ramp** (`family_colours`) and gained a **second
  channel**: a bracket spanning each family's rows with the patient id. The script asserts
  same-family rows are contiguous before drawing the brackets.
- **The 13 all-three models are marked** in a dedicated `All 3` column (ink dot, white halo),
  with the count in the legend title and under the column.
- **Three-variable colour separation verified:** subtype = hue wheel, family = plum ramp,
  presence = ink/pale grey. The closest remaining pair is assay-present ink `#1E2A44` and the
  SCCOHT subtype navy `#0F172A`; they are separated by a gap and a hairline rule and never
  abut. See the palette note in §4.
- **Panel C prints both quantities** as `lines / patients` per row, and the diamond has a real
  legend key instead of a glyph buried in the axis title. Sums asserted (42, 34).
- Height 6.2 → 6.6 in to fit the four-legend block.

### Fig. 2 — `34_fig2_qc.R`
- **Panel A no longer contradicts the text.** It now carries the per-site medians
  (19,817 / 20,765 / 20,119), the relative magnitude (**+948 genes, +4.8%, at −4.1 pp**), the
  **Wilcoxon p = 0.0082 (n = 30)**, `n = 31`, and — new — a **horizontal reference at the
  median genes detected** alongside the pseudoalignment reference. Both truncated axes are
  declared.
- **New panel B: library depth.** Genes detected vs processed fragments, with
  r = 0.50 (p = 0.0038) against r(pseudoalignment, detected) = −0.63, per-site depth medians,
  and both depth p values (0.37, 0.77).
- **Subtype encoding dropped from A/B** in favour of site (fill + shape). Six subtype hues
  plus three site shapes on a 7-pt mark was the main legibility failure, and the panels' claim
  is about centre. Per-line subtype is available in `output/rna_qc_metrics.csv` and Fig. 3A.
- **Panel C: "98% of genes compressed" is now defined** in the panel (protein IQR < RNA IQR
  for the same gene, ratio < 1: 98.0% of 7,896 genes), with the floor-insensitive statistics
  preferred (median IQR ratio 0.30, median SD ratio 0.32, 99.2% < 1) and the
  different-normalisation-scales caveat stated on the panel and in the spec.
- **Panel D sums to 8,427 and says why** (`8,430 search rows − 3 with no gene symbol`), and
  **the 0-plex bar is labelled `70 identified but quantified in no plex (retained, no
  measurement)`** with a leader line. Asserted against `prot_feature_accounting.csv`.
- **New panel E: CV by abundance decile**, both vendor replicates (11.0% → 3.0%, with IQR
  band) and the bridge (34.7% → 7.6%), the latter converted to `100 × (2^SD − 1)` so the two
  series share an axis.
- **Panel F replaces the correlation scatter with Bland–Altman.** Each facet names its cell
  line, `[external]` replaces the unexplained asterisk and is defined in-panel (VOA3993, a
  Carey LGS proteomics-only sample, `provenance = external`, `analysis_include = N`, not one
  of the 42 models), and each facet reports n, bias, SD, repeatability CV and LoA span. **The
  rust band is the median cross-line protein IQR (0.34 log2) and the panel states that the
  LoA span is 2.4–3.1× it.** Pearson r survives only as a caption aside, explicitly labelled
  as uninformative. The existing assertion that r matches `prot_bridge_cor.csv` is retained,
  and a new one asserts every LoA span exceeds the protein IQR.
- Height 5.1 → 6.8 in (three rows).

### Fig. 3 — `35_fig3_biology.R`
- **A/B state n and the unit** (`31 cell lines, one point per LINE, not per patient
  representative`) and the PCA input (top 2,000 VST genes, centred, unscaled).
- **A names the two EC lines** (TOV112D, VOA4395) with leader lines and quotes the
  silhouettes (overall 0.220, **EC −0.014**, SCCOHT 0.819).
- **B shows all three site levels** — the previous 2-level collapse folded BIN67 into
  Huntsman without saying so — and **reports all three commonality components raw AND
  adjusted** (unique subtype 42.4 / 39.3%, unique site 0.2 / −2.5%, **shared 31.1 / 28.9%**,
  joint 73.7 / 65.7%) plus the permutation p values. The annotation moved out of the data
  area into a caption block, which also removes the point overlap.
- **C's n corrected from 8,212 to 7,894** with the ≥10-paired-lines rule stated, and the
  **10.3% negatively-correlated fraction is now shaded and labelled**.
- **D's facet labels are read from the CSVs**: `RNA (22,542 genes)` / `Protein (6,855
  proteins)` — the hard-coded 22,544 / 6,856 at the old lines 116–117 are gone.
- **D plots means beside medians** (open diamond + label below), and annotates `Patient*` as a
  **design artefact** (28 levels on 31 observations, 3 replicated patients, family-restricted
  median **0.76%**). It also states that per-feature percentages need not sum to 100 and that
  RNA has no batch term while protein has TMT plex.
- **E: the EC block is outlined and labelled `EC (not recovered)`** with the failing absolute
  values printed (ESR1 EC mean 0.03, PGR 0.00), and the script asserts all three EC markers
  have `lands_right = FALSE`.
- **E: boxes now encode direction** — solid = expected HIGH, dashed = expected LOW — and loss
  markers carry `↓` in the axis label, so the dark boxed SCCOHT cells no longer read as
  failures. `MKI67 (control)` is labelled as ungraded.
- **New panel F gives the per-line values E cannot show:** all 31 lines as individual points
  on an absolute `log2(TPM + 1)` axis, intended-subtype lines coloured, others grey, with the
  expression-floor reference at 1. The per-line matrix is rebuilt exactly as script 04 builds
  it and **asserted equal to `rna_markers_summary.csv`** — a re-render, not a recomputation.
- Height 5.6 → 7.4 in.

### Supp. S1 — protein PCA
- **Plex moved to `plex_colours`** (ordinal slate); verified that A and B no longer share a
  hue with each other or with any other locked palette.
- n = 31 and the PCA input stated; the panel says plainly that there is **essentially no
  subtype structure** and labels the CC line at the PC2 extreme (**VOA10816**).
- **New panel C: mean silhouette by subtype for both modalities**, with the definition, the
  zero reference, and the per-subtype n, so the negative is quantified rather than asserted
  (protein CC −0.004, SCCOHT 0.028, all 0.135).
- W15 → W2, height 2.6 → 3.1 in.

### Supp. S2 — passage
- **Three nested models instead of two.** Adding `passage | subtype + site` (from
  `unique_passage_beyond_subtype_site`) is what makes the confounding visible:
  **PC1 7.8% → 14.8% → 0.4%**.
- **Uncertainty added** as the `passage alone` p value under each PC, with an explicit note
  that `output/` carries no p or CI for the two conditional quantities.
- **PC5's bar-less conditional value is explained** (raw R² = 0.00004; conditioning can
  increase an R²).
- **The per-gene consequence is stated**: adding passage as a fixed covariate moves the
  **site median from 3.53% to 0.00%**.
- Panel B gains n = 13 **and why only those lines**, both axes labelled `passage number`, the
  diagonal labelled `y = x (identical passage)`, and **Pearson r = −0.33 [−0.75, 0.27],
  p = 0.26** plus the −17 to +20 discordance range. The n is asserted against
  `rna_passage_discordance.csv`.
- W15 → W2, height 2.7 → 3.5 in.

### Supp. S3 — line × line correlation
- **Site annotation strip added** (the single most informative annotation, previously absent);
  the footnote states that the top-level split is essentially by centre.
- **Colour scale now spans the observed off-diagonal range 0.757–0.969** instead of a rounded
  0.70 floor; the footnote gives the range and the median (0.834) and notes the diagonal is 1
  by construction.
- **Family harmonised to `family_colours`.** Feature set and n stated (all 22,544 VST genes,
  31 lines, Spearman, Ward.D2 on Euclidean distance of r).
- **Family co-clustering stated:** 3133 co-clusters as an adjacent pair; **1369 and 2295 do
  not**. Height 6.2 → 6.6 in to hold the footnote.

### Supp. S4 — presence patterns
- **The 70 zero-plex proteins are now a row** (32 patterns, not 31), labelled in rust, and the
  script asserts the pattern counts sum to 8,427 and that the zero row is the zero-plex set.
- **The count arithmetic is surfaced:** 1,572 absent from ≥1 plex **includes** the 70, so
  1,502 are partially observed.
- **Panel A gained row labels** (the protein count per pattern) — it previously had neither
  labels nor counts.
- **Panel B switched from bars on a log10 axis to points**, since bar length is meaningless
  once the axis is logged; every point is value-labelled and the zero-plex point is rust.
- W15 kept, height 4.2 → 5.0 in.

### Supp. S7 — marker effect sizes
- **`"EC"` added to `grp_lvl`** (was line 186) plus an assertion that no marker has
  `grp = NA`. 25 markers in 6 blocks; the EC block is in place and labelled as failing.
- **Bootstrap 95% CIs added as error bars in both panels** (`cohens_d_lo/hi`,
  `auc_oriented_lo/hi`, 2,000 resamples) — the direct referee request.
- **BH significance marked with `*`: only 4 of 25 survive at 0.05** (HNF1B, SPP1, KRT20,
  SMARCA2), with the count printed. **The floor-limited marker is marked with `†`**
  (SMARCA2) and the reason given.
- **The signed-d vs oriented-AUC mismatch is annotated rather than hidden:** axis titles say
  `SIGNED` and `ORIENTED`, `↓` marks loss markers, and the caption states that a loss marker
  succeeds at negative d in A but above 0.5 in B.
- **n per intended subtype printed** (HGS 12, CC 7, EC 2, MC 3, MMMT 2, SCCOHT 2 of 28) with
  the unit made explicit (28 patient representatives, not 31 lines) — this is the only figure
  of mine on that unit.
- Group brackets with n on the right of panel B. W15 → W2, height 4.0 → 5.4 in.

### Supp. S8 — within-HGSC strata
- **Panel A is presented as a partition, not a result.** The panel states that the strata come
  from ward.D2 on the **full 50-set Hallmark z matrix** and that the two plotted axes are the
  theme means used to **name** the clusters by a greedy rule — so the partition is not
  circular but the labels are, and a point's position cannot test its own label. (This is more
  precise than the brief's "points plotted against the two scores that define the strata":
  the scores define the *labels*, not the *partition*.)
- **Stratum colours locked** via `stratum_colours`, mapping the deposited
  `Inflammatory/NF-kB-EMT` string onto the `Inflammatory` key. Consistent with S6.
- **The size encoding is gone:** proliferation z moved to its own panel B with a visual zero
  and a sign. The panel notes proliferation was not used to name the strata.
- **Leader lines added** to all 15 labels in A (`min.segment.length = 0`, increased padding).
- **Panel C quantifies the PROGENy agreement** instead of implying corroboration: ρ = 0.63
  (p = 0.014), 2/5 vs 1/4 above z = +1, and the inflammatory stratum's higher median; plus
  the better-agreeing NF-κB axis (ρ = 0.82, p = 0.0003). n = 15 stated, and the panel states
  that **columns are ordered by assigned stratum and are NOT clustered** while rows are.
- Height 4.2 → 6.6 in (three panels).

---

## 3. What I could not fix, and why

1. **No CI or p value for the conditional passage R² in Supp. S2A.** `output/` carries only
   point estimates for `partial_r2_passage_after_site` and
   `unique_passage_beyond_subtype_site`. Bootstrapping them inside a figure script would be
   new analysis, not rendering, so the panel prints the `passage alone` p values and states
   plainly that no uncertainty is available for the conditional quantities. **Recommended
   upstream fix:** have `scripts/17_variance_confounders.R` bootstrap or permute the
   commonality components and add `*_lo` / `*_hi` columns to
   `output/rna_pc_confounder_joint.csv`, as it already does for `unique_*_perm_p`.

2. **Two Spearman ρ values in Supp. S8C are computed in the figure script.** No upstream CSV
   carries the PROGENy-vs-Hallmark agreement statistics, so `37_supp_rnaprot.R` computes them
   directly from the two deposited columns of `output/hgs_heterogeneity.csv`. This is a
   read-off of a deposited table rather than new modelling, but it is the one place where a
   figure of mine originates a statistic. **Recommended upstream fix:** move the ρ values, the
   per-stratum medians and the above-threshold counts into `scripts/14_hgs_heterogeneity.R`
   and write them to a small CSV, and have the figure read them. The same applies to the
   Pearson r/CI in Supp. S2B (computed from the deposited `rna_passage_check.csv` columns).

3. **Fig. 3F does not use a deposited per-line marker matrix**, because none exists. It
   rebuilds the symbol-level `log2(TPM + 1)` matrix exactly as `scripts/04` does and asserts
   the resulting subtype means equal `rna_markers_summary.csv` to the stored 2 dp (max
   absolute difference 0). That guard makes divergence impossible without a render failure,
   but a deposited per-line marker table would be cleaner. **Recommended upstream fix:** have
   script 04 write `output/rna_markers_perline.csv` (26 markers × 31 lines).

4. **In-panel text is heavier than the lab convention prefers.** The house rule keeps in-panel
   text to axis titles, keys and essential values. Several of these panels now carry 4–10 line
   caption blocks at 5.2–5.4 pt because the brief requires n, unit, thresholds, tests and
   caveats to be legible *in the figure*, and because four of the supplements must be
   self-explanatory enough to cite. Everything in those blocks is also in
   `docs/manuscript/v2/FIGURE_SPECS_A.md`. **If the journal's caption style makes the panels
   redundant, the caption blocks are the first thing to trim** — they are single `labs(caption
   = ...)` calls per panel and can be deleted without touching any geometry. I did not trim
   them myself because the reviewers' complaint was under-specification, not clutter.

5. **Fig. 2A/2B no longer show subtype.** This is a deliberate trade, not a limitation I
   worked around: a reader can no longer check from Fig. 2 whether a QC outlier is one
   subtype. Reinstating subtype means either a third visual channel on a 7-pt mark or a fourth
   legend in row 1. If the reviewers want it back, the cheapest route is a small
   subtype-coloured strip panel rather than re-encoding the points.

6. **Panel-letter drift.** Three of my figures gained panels, so letters shifted relative to
   the previous draft (Fig. 2's bridge panel is now **F**, not D; Fig. 3's marker heatmap is
   **E** with a new **F**; S8's PROGENy heatmap is **C**, not B). Any surviving cross-reference
   in the manuscript body or in `docs/manuscript/v2/` written against the old letters will be
   wrong. The specs file lists every panel by its rendered letter.

---

## 4. Palette additions / changes requested in `00b_figure_theme.R`

I did not edit the theme file. Three requests, in priority order:

1. **Add the deposited stratum label string as a key.** `stratum_colours` has
   `"Inflammatory"`, but `output/hgs_heterogeneity.csv` stores `"Inflammatory/NF-kB-EMT"`.
   I mapped it in `37_supp_rnaprot.R`, and `33_supp_genomics.R` (not mine) reads the same
   column, so both figures are carrying the same local mapping. Adding
   `"Inflammatory/NF-kB-EMT" = "#C2410C"` as a fourth entry would let
   `check_palette_keys()` guard the deposited string directly, which is the behaviour the
   guard is for.

2. **Consider a `singleton` key for `family_colours`.** Both Fig. 1B (blank tiles) and
   Supp. S3 (`singleton = cook_grey`) need a "no multi-line family" level and each invents it
   locally. A locked `"singleton" = "#E2E8F0"` entry would make the two figures agree by
   construction.

3. **Nothing else is needed.** `plex_colours`, `site_colours` + `site_shapes`,
   `family_colours`, `present_colours` and `subtype_colours` covered every variable I plot,
   and `check_palette_keys()` caught nothing at render time after the mappings above — which
   is the intended outcome.

**One near-collision worth recording, not a change request.** In Fig. 1B the assay-present
ink (`#1E2A44`) and the SCCOHT subtype navy (`#0F172A`) are ~1 ΔE apart in practice. They
never abut (the subtype strip is at x = 1, the assay block starts at x = 4.0, with a gap and a
hairline rule between), and both are keyed in separate legends, so I did not treat it as a
defect. If a reviewer raises it, lightening `present_colours["Yes"]` toward `#2B3A55` would
separate them without touching the subtype palette.

---

## 5. Verification performed

- **Every panel viewed as a rendered PNG at final size** before being called done; several
  rounds of caption re-wrapping were needed to eliminate text clipped at the print width.
  Two real clipping bugs were found and fixed this way: the outer patchwork `&
  theme(plot.margin = ...)` was silently overriding per-panel right margins (which clipped
  Supp. S7's group labels), and Supp. S4B's top row needed `clip = "off"` because
  `expand = c(0, 0)` (required to keep its rows aligned with panel A's tile grid) put the
  6,855 point half outside the panel.
- **Assertions added, all of which pass:** Fig. 1 — 42 / 34 / 13, assay totals vs panel A,
  TOV3121D as the only CNV-without-MAF line, family-row contiguity, panel C sums.
  Fig. 2 — n = 31, bridge r vs `prot_bridge_cor.csv`, agreement table completeness, LoA span
  > median protein IQR, panel D bars sum to the analysis matrix, zero-plex and complete-case
  bars vs `prot_feature_accounting.csv`, 10 CV deciles. Fig. 3 — n = 31, one PC1 row, panel C
  rows == `n_pergene_reported`, single `n_features` per assay, one family-restricted row, EC
  markers all `lands_right = FALSE`, **recomputed subtype means == `rna_markers_summary.csv`**.
  Supp. — n = 31 protein lines, palette keys for subtype/plex/site/family/stratum, all three
  passage models present for every PC, cross-assay n vs the discordance summary, presence
  patterns sum to 8,427, no `grp = NA` among markers, bootstrap CI columns present.
- **Console summaries** were extended in all four scripts to print every figure-facing number,
  so a future re-render against the real `output/` will show any drift on stdout without
  opening a PDF.
- **Renders confirmed in the sandbox:** `fig1`, `fig2`, `fig3`, `figs1`, `figs2`, `figs3`,
  `figs4`, `figs7`, `figs8` — `.pdf` and `.png` for each. A final re-render against the real
  `output/` is someone else's step; nothing in these scripts depends on the sandbox.

---

## 6. Isolation audit

Recorded because a clean-room pipeline run was certifying the real `output/` concurrently.

**Every invocation of the four figure scripts had `OVCAN_PROJ` set to the sandbox in the same
shell command as the `Rscript` call.** Confirmed three ways: every run banner printed
`project: /private/tmp/.../figsandbox`; the sandbox `output/session_info.txt` records
`caller: 37_supp_rnaprot.R` with `project /private/tmp/.../figsandbox`; and no file in the real
`docs/manuscript/figures/` has an mtime after 09:48 (my session began 15:33).

**Two Rscript calls in this session did run without `OVCAN_PROJ`** — throwaway verification
scripts `chk2.R` (15:40:24) and `chk3.R` (15:41:30) in the scratchpad, used to check numbers
before editing. They could not have written to the real project and did not:

- neither contains `source()`, `00_setup.R`, or any reference to `PROJ` — so the `PROJ`
  default in `00_setup.R` was never evaluated and `record_versions()` was never called;
- neither contains any `write*`, `save*`, `ggsave`, `dir.create`, `pdf()`, `png()`, `agg_*` or
  `sink` call — they only `cat()` to stdout;
- both addressed the sandbox by hard-coded absolute path;
- **zero files in the real `output/` have an mtime in the 15:35–15:45 window.**

**The real `output/` writes that do exist are the pipeline's own.**
`output/package_versions.csv` and `output/session_info.txt` are both stamped 16:11:42, and
`session_info.txt` self-identifies as written by `21_rna_sensitivity.R` for
`project /Users/dpcook/Analysis/ovcan_human`. The four recent `reports/assets/` files match the
same pipeline scripts (17 at 16:11:42, 19 at 16:07:31) — see the collision table above. The
`Rplots.pdf` in the real project root dates to 09:26, before this session.

**Sandbox `reports/` and `reports/assets/` are real directories, not symlinks** (only
`judy_archive` and `scripts` are symlinked), so the report assets my scripts write went to
`<sandbox>/reports/assets/` — confirmed by their 16:20:33–16:20:46 mtimes there.
