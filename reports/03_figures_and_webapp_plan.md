# Figure Re-engineering, Manuscript Draft & Explorer-App — Orchestration Plan

**Prepared:** 2026-07-24 · **Status:** In progress (sub-agents launched)
**Shared spec:** this document is the single source of truth for the figure/branding contract; every figure agent reads §A + its rows in §B.

---

## §A. Visual branding contract (ALL manuscript figures)

**Foundation module:** `scripts/00b_figure_theme.R` — source it after `00_setup.R` in every figure script:
```r
source("scripts/00_setup.R"); source("scripts/00b_figure_theme.R")
```
It provides: `theme_ovcan()`, locked palettes, scale helpers, ComplexHeatmap col-funs, and `save_fig()`. **Do not define palettes locally.**

**Palette decision table** (which palette for which data type):

| Data type | Palette | Helper |
|---|---|---|
| Centered / signed (z-score, log2FC, log2 copy-ratio, PROGENy) | diverging **navy↔rust** | `scale_fill_cook_div(midpoint=0)` · `cook_div_colfun(max_abs)` |
| Non-centered continuous (absolute expression, %, load) | sequential **rust ramp** | `scale_fill_cook_seq()` · `cook_seq_colfun(limits)` |
| Dense continuous needing perceptual uniformity | viridis/magma (state why) | `scale_*_viridis_c()` |
| Subtype (fixed) | `subtype_colours` | `scale_fill_subtype()` / `scale_colour_subtype()` |
| Source site (fixed) | `site_colours` | `scale_*_manual(values=site_colours)` |
| Somatic tier (fixed, ordered) | `tier_colours` | `scale_*_manual(values=tier_colours)` |
| Present/absent | `present_colours` | `scale_*_manual(values=present_colours)` |
| Generic categorical ≤8 | `cook_categorical` (in order, never skip) | `scale_*_manual(values=cook_categorical)` |
| Single series | `cook_rust` | — |

**NEVER:** rainbow/jet, red-green, diverging map on non-centered data.

**Fonts:** standard sans — `FIG_FONT` resolves to **Arial** (deliberate deviation from the Inter brand font; journal norm + portable embedding). Do not use Manrope/Inter in figures.

**Text minimalism (per PI):** default to **no plot title/subtitle** — the figure caption carries narrative. In-panel text = axis titles, legend, essential stat/value annotations, and panel tags (A, B, …) only. No sentence-length titles.

**Density mandate (per PI):** minimize whitespace, maximize legible content density. Use `theme_ovcan(base_size = 7–8)` for composite panels. Tight `plot.margin`, compact legends (`legend.key.size` already small), trimmed axis expansion (`expand = expansion(mult = 0.02)` where sensible), no oversized point/line sizes. **Legibility is the constraint:** every axis tick, annotation, and legend label must be readable at final print width — test by viewing the exported PNG at 100%.

**Export:** widths `W1=3.46in` (single col), `W15=5.51in` (1.5 col), `W2=7.20in` (double col). Use `save_fig(p, file, w, h)` → cairo_pdf for `.pdf`, ragg for `.png` (dpi 400). **Export at true final size** — do not export large and scale down.

**Output locations:**
- Composite manuscript figures → `docs/manuscript/figures/fig{1..6}.{pdf,png}` and `figs{1..7}.{pdf,png}`.
- Keep/refresh each single panel `reports/assets/f_*.png` (restyled) so the branded HTML report (`build_report.py`) stays current.

**Panel assembly:** patchwork (`p1 + p2 + plot_layout(...)`, `plot_annotation(tag_levels='A')`). Wrap grid-graphics panels (ComplexHeatmap/oncoPrint) with `ggplotify::as.ggplot(grid.grabExpr(draw(ht)))` so they compose in patchwork.

---

## §B. Figure architecture (6 main + supplementary)

> Refinements vs the outline: **Fig 1 replaces Table 1** (PI) — a branded schematic conveys cohort/samples/data types. The **WES filtering waterfall moves into Fig 4** (all WES viz in one figure, clean ownership). **Within-HGSC heterogeneity is demoted to supplementary** (model-selection example only); ADC atlas is the sole main reuse figure (Fig 6).

| Fig | Content (panels) | Source script(s) | Key palettes | Owner |
|---|---|---|---|---|
| **1** | Resource overview — (A) workflow schematic tissue→line→{RNA,TMT,WES}→processing→deposition; (B) **sample×assay coverage matrix** (42 lines, grouped by subtype, patient-family brackets); (C) subtype×patient counts. Programmatic + branded, condensed to figure space. **Replaces Table 1.** | NEW `21_figure1_overview.R` (from `samples.csv`, `line_family_map.csv`) | subtype, present/absent, categorical | **A** |
| **2** | Sequencing & proteome quality — RNA QC (genes/pseudoalign); TMT bridge reproducibility; TMT compression + block-missingness | 01, 05, 19 | subtype, sequential, single-series | **B** |
| **3** | Data recapitulate subtype biology — RNA PCA (subtype + site); variance decomposition; marker recovery; RNA–protein concordance | 02, 17, 04, 12 | subtype, diverging (z), single-series | **B** |
| **4** | WES: raw→validated genomics & identity — filtering waterfall; canonical-driver oncoprint (tier fill, family/subtype tracks); autosome CNV landscape; DepMap external concordance | 07, 08, 18 | tier, cnv-diverging, subtype/family | **A** |
| **5** | Rare-subtype & flagged models — TOV21G hypermutation/MSI (load, indel frac, SBS-96, COSMIC cosine); mucinous ovarian-vs-GI | 16, 11 | sequential, diverging, subtype | **A** |
| **6** | Reuse example — subtype-resolved ADC-target atlas (RNA + protein) | 13 | sequential/diverging (rework heatmap) | **B** |

**Supplementary**

| Fig | Content | Source | Owner |
|---|---|---|---|
| S1 | Proteomic PCA by subtype | 06 | B |
| S2 | Passage sensitivity (within- + cross-assay) | 17 | B |
| S3 | **NEW** RNA sample–sample correlation heatmap (family/subtype annotated; identity/no-swap QC) | `rna_vst.rds` | B |
| S4 | **NEW** per-plex protein presence matrix (structural block-missingness; ggplot presence-grid, *not* an UpSet pkg — none installed) | `prot_block_missingness.csv` | B |
| S5 | **NEW** genome-wide per-line CNV heatmap (lines × bins) | `wes_cnv_segments.csv` | A |
| S6 | **NEW** ConsensusOV calls vs intrinsic HGSC strata (visualize TME-label caveat) | `consensusov_calls.csv` (+ strata) | A |
| S7 | **NEW/expand** marker effect-size panel (all 22; d + AUC) | `rna_marker_effectsizes.csv` | B |
| S8 | Within-HGSC heterogeneity (demoted from main) | 14 | B |

**Ownership (disjoint — do NOT edit another agent's scripts):**
- **Agent A** scripts: `07, 08, 10, 11, 16, 18` + NEW `21`. Figures 1, 4, 5; supp S5, S6.
- **Agent B** scripts: `01, 02, 04, 05, 06, 12, 13, 14, 17, 19`. Figures 2, 3, 6; supp S1, S2, S3, S4, S7, S8.
- Shared read-only: `00_setup.R`, `00b_figure_theme.R`, everything in `output/`, `metadata/`.

---

## §C. Package-default visualizations to re-brand (per PI)

Make these match branding instead of package defaults:
- **Oncoprint (Fig 4):** prefer `ComplexHeatmap::oncoPrint` with `tier_colours` alt fills, thin hairline grid, Arial, compact tracks (patient-family + subtype) using `subtype_colours`; TMB top-annotation as a branded bar. Avoid maftools default rainbow.
- **CNV heatmaps (Fig 4 landscape, S5 genome-wide):** `cook_div_colfun()` (navy↔rust centered at neutral copy ratio); chromosome dividers as light hairlines; subtype/family row/col annotations from locked palettes.
- **ADC atlas heatmap (Fig 6):** rebuild with `cook_seq_colfun()` (absolute expression) or `cook_div_colfun()` (if z-scored); square cells, hairline borders, Arial; RNA and protein sub-panels visually consistent.
- Apply `ht_opt_cook()` once per session before drawing ComplexHeatmaps.

---

## §D. Manuscript draft (Agent: descriptor-draft)

- **Output:** `reports/05_scientific_data_descriptor_draft.md` — full first draft in *Scientific Data* Data Descriptor structure.
- **Structure (fixed):** Title · Abstract (≤170 w) · Background & Summary · Methods (subsectioned per assay + analysis) · Data Records · Technical Validation · Usage Notes · Code Availability. **No Results/Discussion.**
- **Reframing:** follow `reports/02_manuscript_outline.md` §0 mapping — QC/biology-recapitulation/identity → Technical Validation; MSI line / ADC atlas / heterogeneity → Usage Notes. Recast every "we found" as "the data recover / reusers can select."
- **Inputs:** `reports/02_manuscript_outline.md`, `reports/01_multiomic_characterization_results.md`, `reports/00_synthesis_and_recommendations.md`, `reports/lit_review/` (for Background citations), `PROJECT_SPEC.md`, `output/supplement_per_line.csv`.
- **Style:** invoke the **`academic-prose`** skill and write to that standard. Weave real statistics/quantifications inline (effect sizes, n, ρ, %). Calibrated language. Figures referenced by the §B numbering.
- **Deposition:** leave GEO/PRIDE/figshare accessions and any deposition specifics as clearly-marked `[PLACEHOLDER: …]` (PI handles behind the scenes).

---

## §E. Explorer web app (Agent: webapp-explorer)

- **Goal:** a small, hostable app to look up a gene and see its per-line RNA + protein expression across the panel (precedent in spirit: `/Users/dpcook/Analysis/cldn6_ovcan_models/reports/CLDN6_report.{html,pdf}` — branded, self-contained).
- **Data (small enough for client-side):** `output/rna_tpm.csv` (+ `rna_vst.rds`), `output/prot_abundance_matrix.csv`, `metadata/samples.csv`, `metadata/line_family_map.csv`.
- **Deliverables:** (1) `reports/04_webapp_feasibility.md` (feasibility + hosting options: GitHub Pages / Netlify / Cloudflare Pages / shinyapps.io / Claude Artifact — recommend one); (2) a **working self-contained static prototype** under `app/` (data embedded/compressed, in-browser lookup, no server), gene search → per-line RNA + protein, coloured by `subtype_colours`, grouped/sortable; (3) build notes.
- **Branding:** the web app MAY use the full web brand (Inter/Manrope via `~/Lab/Branding/tokens/colors_and_type.css`) — the Arial rule is figures-only. Rust/navy palette; match the CLDN6 report look.
- Keep it **feasibility-first**: prove the smallest useful version works before adding features.

---

## §F. Agent roster & parallelism

Foundation (`00b_figure_theme.R`) is built + tested. All four agents run in parallel; the two figure agents have disjoint script ownership (§B).

1. **fig-overview-genomics** (A) — Figs 1, 4, 5 + S5, S6.
2. **fig-rna-prot-reuse** (B) — Figs 2, 3, 6 + S1–S4, S7, S8.
3. **descriptor-draft** — the SD draft (§D).
4. **webapp-explorer** — feasibility + prototype (§E).
