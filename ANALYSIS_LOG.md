# Analysis Log: OvCAN Human Ovarian Cancer Cell-Line Multi-Omic Resource

## Status
**Current stage:** Phase 8 — descriptor authoring + figure re-engineering COMPLETE (pending PI review). PI approved the outline and directed: draft the manuscript (sub-agent + `academic-prose`), re-engineer all figures programmatically (patchwork/cowplot) with the updated branding, and explore a gene-explorer web app. All three delivered + verified. **Foundation built + tested:** `scripts/00b_figure_theme.R` (shared theme/palettes: new **navy↔rust diverging** map per branding commit 7b1f3dd, locked subtype/site/tier palettes, Arial figure font per PI, text-minimalism + density mandate, `save_fig()` + column widths). Plan documented in `reports/03_figures_and_webapp_plan.md`.
**Current task:** ALL 4 Phase-8 agents COMPLETE + orchestrator-verified. Deliverables: **6 main figures + 8 supplements** (`docs/manuscript/figures/fig1–6`, `figs1–8`; branded via `00b_figure_theme.R`, pure Arial in all 14 PDFs — theme-wide cairo_pdf font bug found + fixed); **descriptor draft** (`reports/05_...draft.md`); **explorer app** (`app/` + `reports/04_webapp_feasibility.md`, renamed "Gene Explorer: OvCAN Collection", white bg). Everything ready for PI review. **Optional later:** small cosmetic figure-polish pass (legend spacing, 1 annotation, 1 caption note); make white the default bg in the HTML report too (open Q). **Still last:** `renv::snapshot`. **Standing PI decisions:** HRD; mucinous VOA8762/8771; deposition specifics (placeholders in draft); STR/IHC + TOV21G MSI confirmation; author list + code host. **Still deferred to end:** `renv::snapshot`. **Standing PI decisions unchanged:** HRD; mucinous VOA8762/8771; deposition (accessions as placeholders per PI — team handles behind the scenes; MAF GRCh37→GRCh38 header; strip third-party normals; PRIDE raw MS); STR/IHC + TOV21G MSI-confirmation; author list + code host.
**Last updated:** 2026-07-27 (methodological audit response — five statistical corrections + two factual corrections applied across scripts, outputs, the results report and the v4 manuscript; see the top log entry. `reports/01_multiomic_characterization_results.pdf` now 38 pp. **Outstanding before submission:** a fresh clean-room re-execution, and gaps 16–22 in the v4 register — RNA re-quantification against a matched index, deposition with data dictionaries/checksums/reviewer access, ethics and Human Data Checklist, author contributions/competing interests/funding, locked references.)

## Quick Summary
Re-scoping an MSc thesis (Judy Sobh) into a *Scientific Data* Data Descriptor for a multi-omic (bulk RNA-seq + TMT proteomics + WES) resource of human ovarian cancer cell-line models (OvCAN / Ovarian Cancer Canada collection; CHUM Mes-Masson TOV/OV + BC Cancer Huntsman VOA + OHRI lines). Phase 1 (full assessment of Judy's archive) and Phase 2 (6-cluster literature review) are essentially complete. Emerging reframe: **descriptor, not discovery**; report coverage per-assay (only **13 lines are truly tri-omic**); WES **mutations are tumor-only → artifact-prone** (CNV is sound); several claims (ConsensusOV subtyping; "novel" EC–SCCOHT convergence; expression-"HRD") must be dropped/relabeled. Source-of-truth sample sheet built (`metadata/samples.csv`).

---

## Log

### [2026-07-27] — Methodological audit response: five statistical corrections, two factual corrections, five new deposited tables
**Type:** Analysis / Report / Manuscript · **Phase:** Pre-submission validity audit · **Status:** Complete
**Trigger:** `reports/06_methodological_validity_and_submission_readiness_review.md` (independent methodological validity and *Scientific Data* submission-readiness review)
**Scripts edited:** `01_rna_load_qc.R`, `04_rna_markers_genesets.R`, `12_rna_protein_concordance.R`, `19_proteomics_dynamic_range.R`, `34_fig2_qc.R`, `docs/manuscript/v4/build_docx.sh`
**New outputs:** `rna_reference_sensitivity.csv`, `rna_marker_recovery_permutation.csv`, `prot_noise_vs_biology.csv`, `integ_rnaprot_negative_genes.csv`, `integ_rnaprot_patientrep_sensitivity.csv` (output/ now 107 CSV + 6 rds + 13 txt + 2 md)
**Documents updated:** `reports/01_multiomic_characterization_results.{md,html,pdf}` (38 pp), `docs/manuscript/v4/OvCAN_data_descriptor_v4.{md,docx}`, `reports/02_manuscript_outline.md` and `reports/05_scientific_data_descriptor_draft.md` (both marked SUPERSEDED with inline corrections)

**What was done — the review's substantive findings were accepted and implemented:**

1. **Proteomic bridge repeatability was on the wrong scale.** `sd(primary − bridge)` is the SD of a *difference of two measurements*, so a per-measurement SD needs the replicate-difference adjustment `/√2`. The old `100 × (2^sd_diff − 1)` = 15.7–20.4% applied neither that nor the standard lognormal conversion. Corrected: per-measurement SD **0.149–0.190 log2**, CV **10.4–13.2%**. Legacy column retained as `sd_diff_cv_pct_legacy`, flagged deprecated.
2. **The "technical noise exceeds biological spread by 2.4–3.1×" claim was an artefact and reverses on correction.** It divided a 95% LoA span (3.92 SD, paired-difference scale) by a cross-model IQR (1.35 SD, single-measurement scale) — a ~2.9× inflation from the choice of dispersion measure alone. Like-for-like (SD vs SD): observed cross-model protein SD 0.287 log2 is **1.5–1.9×** the technical SD, so technical error is **27–44%** of observed cross-model variance and single-measurement reliability is **0.56–0.73**. Better for the 8 ADC targets (technical share 6–10%). New table `prot_noise_vs_biology.csv`. Fig. 2F's reference band is now the 95% span of a cross-model difference (1.59 log2), so both quantities are difference intervals — the technical span is **0.52–0.66×** it, and the LoA lines now sit *inside* the band. Fig. 2E's bridge series is now per-measurement (21.3% → 5.2% across deciles, was 34.7% → 7.6%).
3. **Marker recovery now uses a correlation-preserving joint permutation null.** The exact binomial treated 25 markers as independent Bernoulli trials; same-histotype markers are co-regulated, so one relabelling moves several together. 20,000 joint permutations of the histotype labels at seed 1234: **16/25 on 31 models, p = 0.0013**; **15/25 on 28 patient representatives, p = 0.0046** (null mean 6.7–6.8, i.e. ~27% not the binomial's 33%, because the up rule's expression floor makes landing harder). Conclusion strengthened. Deposited to `rna_marker_recovery_permutation.csv`.
4. **The 28 patient representatives are now the explicit inferential unit** in both the report and the manuscript, with the 31 line models the descriptive unit; the split is stated in §1.2 of the report and in Methods, and the reuse-defaults table now names them separately.
5. **The "3.34× compression sets the concordance ceiling" claim was retired.** The spread ratio is a real property of the two matrices but cannot be attributed to reporter-ion compression alone — post-transcriptional buffering, turnover, pooled-standard ratioing and structured missingness all narrow the protein axis, and the layers sit on different scales. Reported as "consistent with".
6. **The 813 negative transcript–protein correlations were reported as established anti-correlation.** With n = 10–30 the SE of ρ is ~1/√(n−3), so most are indistinguishable from zero: **51 nominally significant, 30 (0.38%) BH-significant**, against 3,626 BH-significant positive. Per-gene p and q now deposited in `integ_rnaprot_cor.csv`; tiers in `integ_rnaprot_negative_genes.csv`. Reframed as candidates for discordant regulation.
7. **Histotype-versus-centre language is now associational throughout.** The commonality analysis establishes only that centre adds nothing beyond histotype; the PC1 shared component of 0.311 is not apportionable because no design cell has the two varying independently. New limitation 5b.
8. **Two factual corrections.** `BIN67` was **not** reclassified by Karnezis 2021 (that paper reassigned `COV434` and `TOV112D`) and is no longer counted as a recovered reclassification — the count is one corroborated, not two independently recovered. And the v4 draft contained a sentence stating that STR profiles and mycoplasma clearance "were obtained for the panel"; **neither was performed**, and the sentence was removed and replaced with the limitation.
9. **New sensitivity analysis, which returned a positive rather than a null result.** Testing whether the Ensembl-104-index / release-105-map mismatch is design-neutral: the per-model dropped-TPM fraction has **R² = 0.330 on contributing centre (p = 0.0033) and 0.388 on histotype (p = 0.010)**. The loss is small and bounded (1.60–3.41% of TPM, a 1.81-pp range) and the reference pair is identical for every model, but it is *not* design-neutral — so matched-index re-quantification was promoted to a "would change the science" open item (§12 item 5) rather than argued away. Deposited to `rna_reference_sensitivity.csv`.
10. **Other wording corrections accepted:** somatic-confidence tiers defined as a prioritisation heuristic rather than a validated classifier; "no defensible *somatic* *BRCA1/2* call can be made" rather than "defensible somatic BRCA1/2 is zero"; the 30-fold FGA contrast presented as a panel range (18 models vs 1) not a histotype effect estimate; copy number described as relative total-copy with unconfirmed pooled-normal capture compatibility; GO programme grades disclosed as regex-selected best terms (programme prespecified, term not); SWI/SNF calls described as candidate subunit-loss flags with uncalibrated thresholds; consensusOV explicitly an exploratory annotation; the 98.0%/97.8% compression discrepancy between report and manuscript resolved (the figure script had `na.rm` silently dropping 14 non-finite ratios).
11. ***Scientific Data* genre compliance.** Added the required **Data Availability, Ethics, Author Contributions, Competing Interests, Funding, Use of Generative AI** and **References** sections; abstract cut **176 → 164 words** (title unchanged at 93 chars); Data Overview reduced and Fig. 6B marked explicitly exploratory with a gap item to relocate it to the Supplementary Information; gap register grown from 15 to 22 items; four new supplementary tables (S14–S17) specified. `build_docx.sh`'s header strip re-anchored on the Authors line so it no longer breaks when a provenance note is added.

**Not changed, with reasons:** the review's points 4, 7, 9, 10 and 12 were largely already satisfied by the existing text (two-of-five pathway grading, protein plex/site structure, tumour-only caveats, capture-kit disclosure, the patient variance term flagged as a design artefact); those sections received wording tightening rather than substantive revision. The reviewer's own permutation p-values (0.0019 / 0.0052) differ from the deposited ones (0.0013 / 0.0046) only through the seed; both give the same conclusion, and the project seed 1234 is now canonical.

**Verification:** all four edited analysis scripts rerun to exit 0; Fig. 2 regenerated in both plain and manuscript modes and the report rebuilt (38 pp, 14 SVG figures embedded). A full clean-room re-execution has **not** been rerun since these edits and should be, before submission — the previously recorded 122-output figure is now stale at 128.

**ROUND 2, same date — six refinements from a second review pass, all accepted.** Recorded in `reports/06_changelog_methodological_validity_and_submission_readiness_review.md` §12.

1. **`reliability_icc` was not an ICC.** It mixed a bridge error distribution pooled *across proteins* (from four different samples) with an observed SD computed per gene and then medianed *across genes*, under an unstated homoscedasticity assumption. Renamed **`approx_reliability_ratio`**, with the four separating approximations enumerated in `19_proteomics_dynamic_range.R` §4b-ii and in both documents. Phrasing everywhere is now "under a homoscedastic independent-error approximation, bridge variability corresponds to roughly 27–44% of the variance of a typical protein's observed cross-model measurements" — never "reliability 0.56–0.73".
2. **Dropping homoscedasticity replaced the global figure with a gradient, and this is the substantive gain.** Matching technical and observed spread *within* mean-abundance deciles (both decile-local) gives a technical variance share of **67.5% in decile 1 falling to 7.5% in decile 10** (approx. variance ratio 0.33 → 0.92). So whether a single protein in a single model is interpretable is **abundance-conditional**, and the global 27–44% understates the problem at low abundance and overstates it at high. New columns in `prot_cv_by_abundance.csv`.
3. **The ADC comparison is now genuinely per-class.** Round 1 swapped in the ADC-target observed SD but kept the *global* technical SD. Bridge differences restricted to the ADC target rows give technical SD **0.134 log2 from 32 differences** against observed 0.600 → technical share **≈5%** (was reported as 6–10%). Quoted with its n; the abundance-matched deciles agree.
4. **The reference analysis was mislabelled.** It is a **design-alignment diagnostic**, not a matched-index sensitivity analysis: the association it finds (dropped-TPM fraction on centre/histotype) is equally consistent with genuine histotype differences in expression of the omitted transcripts, so it flags a risk and cannot demonstrate a harm. The round-1 changelog sentence "I implemented it and it did the opposite" is retracted. Matched-index re-quantification remains the actual, **undone** sensitivity analysis (report §12 item 5, manuscript gap 16).
5. **Spearman inference and "proxy" terminology narrowed.** Switched from a Fisher z transform scaled by 1/√(n−3) — a Pearson-derived approximation applied to a rank correlation — to the **conventional asymptotic Spearman t approximation (df = n − 2)**, which `cor.test` uses under ties. Counts move **51 → 50** nominal and **30 → 29** FDR (FDR-positive 3,626 → 3,599); nothing substantive turns on it, but the conventional test is reproducible with `cor.test()`. Separately, the FDR-supported subset is now described as an **inverse association**, not as genes for which transcript "fails as a proxy" — an inverse relationship can still be predictive, and proxy value additionally depends on predictive error and external validation, neither assessed here.
6. **A protein-key join defect, repaired.** `19_proteomics_dynamic_range.R` §4c keyed the abundance table on the decorated matrix row id (`SYMBOL|UNIPROT` for the 31 non-representative rows) but keyed the bridge differences on the raw `abund$Symbol` — a miss for every decorated row and a many-to-many fan-out for the **44 workbook rows sharing 13 symbols**, with `relationship = "many-to-many"` silencing the warning that would have caught it. `(Symbol, Uniprot)` is unique on both sides and is now the join key, resolved to the canonical row id in a new §4b-i and asserted one-to-one. Effect is small as predicted (decile-1 SD 0.430 → 0.429, decile-10 0.106 → 0.105; Fig. 2E endpoints unchanged at 21.3% → 5.2%) but it was wrong.

**Also corrected: one mischaracterisation of my own.** Round 1 filed the commonality-analysis discussion under "Not accepted", describing the review as saying the analysis "adds nothing". The review had explicitly called it "a useful improvement over a simple PCA colour plot". Recategorised as *accepted, evidential boundary clarified*. Related over-read fixed: **donor collapse rules out duplicate donor families, not a centre-aligned batch effect** — the 28 representatives come from the same centres in the same proportions — and the commonality partition excludes centre as a source of *additional* variance, not a centre-aligned batch effect, which lives in the shared component by construction. Report §4.2 and the summary callout now say so, so the two narrow negative results cannot accumulate into an implied positive one.

**Round-2 verification:** `19` and `12` rerun to exit 0, Fig. 2 regenerated in both modes, report rebuilt (38 pp), `.docx` rebuilt (37 gap markers). The second reviewer independently reproduced two results — the conventional Spearman test (50/29) and the donor-collapse comparison on the identical 7,832-gene set (0.3976 → 0.4168, inverse fraction 0.1023 → 0.0970) — both confirming the deposited values. **The clean-room rerun is now the first outstanding item**; the analytical content is settled.

### [2026-07-26] — Consolidated results report rebuilt (37 pp) covering the full analysis, with clean-figure mode
**Type:** Report / Script · **Phase:** Post-v4 consolidation · **Status:** Complete
**Outputs:** `reports/01_multiomic_characterization_results.{md,html,pdf}` · `reports/assets2/` (14 svg embedded + 14 pdf + 14 png + `panels/` 18 png) · `scripts/build_report_full.py` · previous report + edited scripts archived to `reports/_archive/2026-07-26/`

**What was done:**
- **Rewrote the internal results report from scratch** to cover everything examined since the 2026-07-24 version (which predated v2–v4 of the manuscript, the evidence dossier, and the 2026-07-25 fixes). 12 sections + 2 appendices, ~11,500 words: composition and the ten denominators; per-layer methods and what is not recoverable; per-layer quality; what structures the expression data (commonality, permutation nulls, within-CC control, passage, subline collapse); biology recovered (DE grading, GSEA grades, marker effect sizes); the transcript–protein relationship; genomic fidelity; signatures and the TOV21G MSI candidate; identity and authentication; the reuse layer; the 23-item limitations register with a number on every row; and open items framed by what each one would change. Bullets carry the numbers, prose carries the interpretation.
- **Numbers re-verified against `output/` directly**, not transcribed: 30+ headline quantities recomputed in R (RNA QC medians, concordance, FGA by histotype, driver tiers + Tier-3-only genes, silhouettes by modality, PC1 commonality raw/adjusted, marker counts incl. the 4 BH survivors, arm frequencies + unit, DepMap margins, compression, bridge LoA, burden, GSEA padj, FOLR1 bimodality, HGS strata, annotation-support decomposition). All matched. The v4 correction (8 positive-lineage + 16 absence-only + 2 circular = 26) was confirmed from `auth_perline_table.csv`.
- **Figure set regenerated in a new "plain" mode** so the report does not inherit the manuscript figures' in-panel methodological footnotes (PI: "weird annotation text"). Added an **env-gated** section 8 to `scripts/00b_figure_theme.R`: with `OVCAN_FIG_PLAIN=1`, `save_fig()` strips every plot caption from the object tree (patchwork-aware) and a shadowed `annotate()` drops long/multi-line in-panel note blocks; three ComplexHeatmap/leader-line cases needed explicit guards (`33:draw_s5`, `37:draw_s3`, `34` panel-D leader). Geometry, palettes, locked colour mappings, panel tags and value labels are untouched — 11 note blocks and all plot captions suppressed, nothing needed to read a panel removed.
- **Proved the manuscript path is unaffected:** re-ran `34`, `33`, `37` with the flag unset in an isolated sandbox and compared to the pre-edit backup — `fig2`, `figs1–figs5`, `figs7`, `figs8` all **pixel-identical** (max delta 0). `docs/manuscript/figures/` and `reports/assets/` were never written (mtimes unchanged); all figure runs used `OVCAN_PROJ` pointed at a scratch sandbox with symlinked inputs and real copies of the two files `00_setup.R` rewrites.
- **New builder `scripts/build_report_full.py`** (the old `build_report.py` is untouched and still drives the legacy report). Markdown stays the single source of truth and remains readable standalone; three constructs are expanded at build time — `{{stats}}`, `{{callout: …}}`, `{{figure: id | caption}}` — with figure labels derived from the id so `fig3 → Figure 3` and `figs7 → Figure S7`, matching the `Fig. N` references in the prose. Branded Cook Lab sheet (rustNavy tokens, Manrope/Inter/JetBrains Mono) matched to `00b_figure_theme.R`.

**Decisions:**
- Update the report **in place** at `reports/01_...` rather than creating a parallel file, so there is one canonical path; the 2026-07-24 md/html/pdf, the pre-edit figure scripts and both figure directories are archived under `reports/_archive/2026-07-26/`.
- Report figures live in `reports/assets2/` (not `assets/`), so the legacy report keeps rendering unchanged.
- Use the **14 composite figures** rather than the 18 single-panel exports: the composites are the canonical branded set and splitting them would lose panels. Panels are still copied to `reports/assets2/panels/` if a section ever wants one.
- **`academic-prose` editing pass over the report text** (PI request). Audit → rewrite → verify, in eight batches of exact-match replacements so nothing could drift. Removed: **em-dashes 109 → 14** in body prose (the 14 remaining are structural: title subtitle, stat-card separators the builder renders as line breaks, empty table cells, and one verbatim quoted file string); 10 self-referential emphasis scaffolds ("worth stating", "deserves a note", "should not be buried"); 3 dramatic colon reveals ("The consequence is blunt:", "The honest statement is that", "So:"); all `Notably`/`Critically`/`Note also` boilerplate; decorative `decisively`/`precisely`. Split 14 clause-stacked sentences past 40 words. Word count 14,099 → 13,788. Kept deliberately: the directive register appropriate to an internal reference ("do not use the word bimodal", "Report probabilities and margins, never bare labels"), evidence-tied hedging (`candidate`, `consistent with`, `may be n = 1`), and every technical term. **Verified no content moved:** numeric-token multiset identical at 2,054 tokens (zero lost, zero gained), and gene symbols, backticked identifiers, § cross-references, figure references, citations and all 24 table shapes unchanged.
- **Heading revision** (PI: headings were "ambiguous and flowery", and the same pattern recurred wherever the text made a declarative claim). Every heading is now either a plain description of content or a declarative statement with the finding in it, on the manuscript convention that a reader should be able to scan the headings alone and know the results. Examples: "Screening and attribution are different analyses, and the words matter" → "Cosine screening is inconclusive; the bootstrap refit supports a mismatch-repair exposure in TOV21G alone"; "The mechanism is consistent with what an exome cannot see" → "No causal coding lesion is present, consistent with epigenetic *MLH1* silencing"; "Commonality analysis carries the weight" → "Commonality analysis: histotype explains the leading components, centre adds nothing"; "The collinearity should be seen, not only modelled" → "Unsupervised clustering of all 22,544 genes splits by contributing centre". 49 headings revised, plus 14 in-body claim openers carrying the same obliqueness ("Two things this is **not**", "This divergence should not be averaged away", "One trap is easy to propagate", "The partition is not circular; the labels are"). Convention set: gene symbols italic in headings, cell-line names plain (so the §8 title no longer italicises TOV21G); `LGSC` replaced by the abbreviation the report actually defines. **Verified:** every number newly surfaced in a heading already occurred verbatim in the body, and no data value anywhere in the body changed — the only body-text token additions are one `BIN67` model mention and one §4.2 cross-reference.
- Print layout: letter, 10 mm margins, figure image capped at **200 mm** — above ~205 mm the figure block exceeds one page, `break-inside: avoid` cannot place it and Chrome emits a near-blank page. Fig 3 was moved to the end of §4.2 and Fig S4 kept under its own heading for the same reason. Final: **36 pp, mean vertical fill 76%, no page below 28%.**
- **Figures embedded as vector SVG, not PNG** (PI question). Tested rather than assumed: `pdftocairo -svg` converts the `cairo_pdf` originals, Chrome preserves SVG as vector through print-to-PDF, and the result is **3.4 MB vs 6.4 MB for the PNG build at the same pagination** — vector is both sharper *and* smaller here, because 400 dpi rasters of 14 full-width figures cost more than the paths do. Raster content in the output PDF drops from **14 images / 95.7 Mpx to 22 images / 4.6 Mpx**, the remainder being the layers deliberately rasterised in R with `ggrastr` (the bridge and concordance scatters), which is the correct behaviour. Print render time is unchanged (3.0 s vs 3.8 s). `build_report_full.py` regenerates a stale SVG automatically when its PDF is newer, and `--figs=png` falls back. **Trade-off recorded:** pdftocairo converts glyphs to paths, so figure text is font-independent and crisp but not selectable in the PDF (body text still is — 15,282 words extractable). Since vector figures are zoomable, the printed-size/legibility trade-off that governed the raster build no longer binds; pagination proved insensitive to the image cap between 160 and 200 mm, so the larger value was kept.

**Issues / notes:**
- `mdls kMDItemNumberOfPages` returns a stale Spotlight value after overwriting a PDF in place — it reported 23 pp for a 39-pp file. Use `pdfinfo`.
- The report records that `cooklab.ca/ovcan_viewer` already exposes the full RNA + protein matrices publicly, and flags reconciling that with the deposition/embargo plan as an open item.

### [2026-07-24 17:05] — Integration follow-up: tx2gene relocated to data/reference/ with a checksum assertion
**Type:** Decision / Script · **Phase:** Revision (post-repair integration) · **Status:** Complete
**Script:** `scripts/01_rna_load_qc.R:49-115`

**What was done:**
- **Established what actually happened to the tx2gene map during the clean run.** It was NOT preserved by me: in BOTH passes `01` fell back to a **live biomaRt query** against `dec2021.archive.ensembl.org`, because the map lived only in `output/` — the directory documented as "safe to regenerate". Pass-2 log line 6 is explicit. So the certified verdict was "reproducible **with network access** to the Ensembl 105 archive", not offline. Two claims, and only the weaker one had been earned.
- **The result itself is a strong reproducibility datum:** two independent live queries ~30 min apart returned a file byte-identical to the hand-cached copy and the baseline snapshot — md5 `119dfe0ab4f856e6b81efc6ce78f4ba7`, 13,884,511 bytes, 266,615 transcripts, all four copies.
- **Relocated** the map to `data/reference/tx2gene_ensembl_rel105.csv`. `01` now loads it from there, **asserts md5 + 266,615 rows + the 4 column names**, and refreshes the `output/` working copy by **byte copy** so the other 7 consumers (`03, 04, 12, 13, 17, 19, 21`) need no change. biomaRt fallback retained but now emits a loud `warning()`.
- **Verified offline:** with `output/tx2gene_ensembl_rel105.csv` deleted, `01` loads from `data/reference/`, logs `tx2gene verified: md5 119dfe0a...`, reports 39,568 genes → 22,544 filtered, exit 0, **no network query**. All 8 files `01` writes are byte-identical to the certified pass-2 copies, so nothing moved.

**Decisions:**
- Assert on the **file on disk** and propagate with `file.copy()`, never `write_csv()`: my first attempt checksummed a re-serialised tibble and failed against the very file it verified, because a `read_csv`→`write_csv` round-trip is not byte-preserving here (empty `external_gene_name` → `NA`).
- Keep a working copy in `output/` rather than repointing 7 scripts: smaller diff, and `output/` stays the single place downstream scripts look.

**Issues:**
- **RETRACTED a wrong claim of mine.** I had reported the concurrent figure-script processes as contamination and said the figure work was void. They were running with `OVCAN_PROJ` pointed at an isolated sandbox, so they never touched the real `output/`/`metadata/`/`docs/` (all 28 files in `docs/manuscript/figures/` stayed byte-identical to the pre-revision snapshot). My exclusivity guard keyed on the **script name**, which is over-broad; the correct test keys on the **target root** (`OVCAN_PROJ` via `ps -Eww`, or `lsof +D output/`). Cost was ~10 s of waiting. The pass-1 contamination by the RNA workstream was real and is confirmed by that agent.
- Also verified none of `00`–`22` sources `00b_figure_theme.R` or uses its palettes (`17`/`18` define local `save_fig()`), so concurrent theme edits cannot affect the chain.
- Write coverage now fully closed: the 23 apparently-unreferenced outputs are constructed filenames (`sprintf` in `03`'s subtype loop; `write_session_info()`), and with tx2gene resolved **no file in `output/` is written by no script**.

### [2026-07-24 16:40] — Integration test: clean-room reproducibility run of the full pipeline
**Type:** Verification / Script · **Phase:** Revision (post-repair integration) · **Status:** Complete
**Scripts:** all of `scripts/` except `3*` · **Reports:** `reports/integration_report.md`, `reports/integration_diff.md`, `reports/integration_logs/`

**What was done:**
- Emptied `output/` (122 files moved aside, `output/external/` never touched) and deleted `metadata/line_family_map.csv`, then ran the whole chain sequentially in dependency order with per-script logs and timings. **VERDICT: a clean run from raw inputs reproduces the full output set** — 23 invocations, 0 non-zero exits, all 122 files regenerated, 1,775 s (~30 min) total.
- Comparison vs `reports/_snapshot_postfix/output/`: **109 identical, 12 timestamped `session_info*.txt`, 1 explained (bootstrap cache), 0 missing, 0 value/column differences.** All **25** headline numbers PASS (557,392→15,692→15,609→6,036, median 206.5, TOV21G 1,416, 22,542 genes, 8,427 = 6,855 + 1,572, 70 zero-plex, 42 supplement rows, 19 genomics-consistent, arm freqs 82/91/82/73/64/55). `07`'s sha256 content guard reproduces exactly, so the C1 fix genuinely works from an empty `output/`.

**Decisions:**
- **Zero-plex approach (the item the RNA agent could not finish):** keep `prot_abundance_matrix.csv` numeric-only and add ONE shared loader `read_prot_matrix()` (`00_setup.R:136-181`) used by `10`, `11`, `12`, `13`, `19`. Rationale: a flag column would coerce the matrix to character for every consumer doing `as.matrix(x[, -1])`, including external reusers. The loader *derives* the 70 zero-plex rows and asserts them against `prot_qc.csv$zero_plex`, so the two files cannot drift. No denominator moved.
- **Abandoned pass 1 as evidence** because the RNA repair agent was concurrently executing scripts into the same `output/` (session_info mtimes 15:20-15:29 during my `03` run). Re-ran a second, exclusive pass. Cost: one extra ~30 min run; benefit: the verdict is certifiable.

**Key outputs:**
- File: `reports/integration_report.md` (verdict, failure log, cross-half reconciliation, regeneration table, non-determinism, headline table, runtime, residual risks, canonical **Execution order** for Methods/Code Availability)
- File: `reports/integration_diff.md` (per-file 122-row comparison)
- File: `reports/integration_logs/` (per-script logs, `timings_pass1.tsv`, `timings_pass2.tsv`, `compare_outputs.R`, `verify_headlines.R`, `pass1_08_wes_cnv_familymap_guard.log`)

**Issues:**
- **Fixed — `ensure_family_map()` guard was incomplete:** `17`, `20`, `21` read the generated `metadata/line_family_map.csv` bare (same class of bug as the original `wes_mutations_filtered.csv` defect). Added at `17:85`, `21:56`, `20:46-49`.
- **Fixed — bootstrap cache reusable while incomplete:** `22`'s validity test checked only cache keys, so a short cache passed and the join's `replace_na(boot_selected_frac, 0)` would silently DROP signatures from `wes_signature_refit_exposures.csv`. Added a complete-grid requirement (`22:202-213`) + a coverage `stopifnot` (`22:262-265`); both reject-and-recompute and reuse paths tested.
- **Non-determinism: none found.** No seed fix needed. Two independent from-scratch 200-replicate bootstraps produced byte-identical caches.
- **Runtime estimates were far too pessimistic:** `22` = 570 s from scratch (not ~40 min), `17` = 243 s (not ~40 min). Whole pipeline ~30 min, so a full clean run is a practical pre-submission check.
- **UNRESOLVED (not mine):** stale hard-coded figure labels remain — `35_fig3_biology.R` (22,544→22,542; 6,856→6,855; 26 marker rows; panel C n=7,894) and `37_supp_rnaprot.R:186` (add `"EC"` to `grp_lvl`). Figure scripts were run by another agent against a half-populated `output/` during this test and should be rebuilt.
- **UNRESOLVED:** `output/` has no lock; concurrent agents interleave writes (this cost a full re-run). A lockfile in `00_setup.R` would fix it.

### [2026-07-24] — Phase 8: descriptor draft exported to Word for PI review
**Type:** Script · **Status:** Complete
- Converted `reports/05_scientific_data_descriptor_draft.md` → **`docs/manuscript/OvCAN_data_descriptor_draft.docx`** (pandoc; no branding per PI — plain Word with native Heading styles + a real table for Data Records, for comments/track-changes). Preprocess: dropped HTML build comments, promoted the real title to H1; body verbatim incl. inline `[PLACEHOLDER]`/`[CHECK]`/`[NOTE]` flags. Verified valid docx (30 KB, "Microsoft Word 2007+") and all 9 sections + title round-trip correctly. Markdown at `reports/05_...` remains the source of truth.

### [2026-07-24] — Phase 8: explorer app — single-file build, cooklab hosting, default gene, DEPLOYED
**Type:** Script · **Status:** Complete — committed + pushed; **LIVE at cooklab.ca/ovcan_viewer**
- **Hosting path resolved:** `cooklab.ca` is Astro→Vercel; Astro serves `public/` at the site root, so the app deploys at **cooklab.ca/ovcan_viewer** via `public/ovcan_viewer/index.html`. New `app/build_single.py` inlines `data.js` into `index.html` → self-contained `app/ovcan_viewer_standalone.html` (3.47 MB; no external refs → robust to trailing-slash / relative-path, no server config). Copied to `/Users/dpcook/Projects/cooklab/public/ovcan_viewer/index.html` (untracked/local; deploy = user's `git push`, Vercel auto-deploys). Flagged: a public URL releases the full unpublished RNA+protein matrices → align with deposition/embargo. Deploy steps in `app/README.md`.
- **PI request:** default highlighted gene **CLDN6 → WT1** (ENSG00000184937) — "no reason to draw attention to CLDN6"; also demoted CLDN6 from first "Try" chip + placeholder to last. Changed in source `app/index.html`, rebuilt single-file, recopied. Verified default renders WT1.
- **Deployed (user-authorized):** the app was already live from the user's earlier push (commit `cc94e83`, CLDN6 default); the WT1 fix was committed (`4410329`) + pushed to `origin/main` → Vercel auto-deployed. Verified the live URL `https://cooklab.ca/ovcan_viewer/` serves `show(urlGene || "WT1")`. **The full RNA+protein matrices are now public** (user's decision) — relevant to descriptor deposition timing.

### [2026-07-24] — Phase 8: Agent A figures COMPLETE + theme-wide Arial fix — ALL figures done
**Type:** Script / Verification · **Phase:** 8 · **Status:** Complete (all 4 Phase-8 agents done)
- **fig-overview-genomics** delivered scripts `30–33` → `fig1` (resource overview — REPLACES Table 1: workflow schematic + 42×4 coverage matrix + subtype line/patient counts), `fig4` (WES: filtering waterfall 25,914→493→170 + oncoprint + DepMap 5×5 self-match + per-arm CNV gain/loss), `fig5` (TOV21G hypermutation/SBS-96/COSMIC + mucinous ovarian-vs-GI), `figs5` (genome-wide per-line CNV heatmap), `figs6` (ConsensusOV vs intrinsic strata). Refreshed 5 `reports/assets/` panels.
- **Oncoprint rebuilt as native ggplot** (ComplexHeatmap::oncoPrint clips under patchwork `grid.grabExpr`): tier fill + truncating/missense box-height + subtype/family/coding-load tracks + patient-count bar; composes cleanly. CNV split non-redundant: Fig4 = per-arm gain/loss freq (foregrounds collapse-robustness); figs5 = full per-line genome heatmap. Orchestrator viewed fig1/fig4 at true size: high quality, branded, legible.
- **THEME-WIDE ARIAL FIX (cross-cutting, flagged by Agent A):** `theme_ovcan()` `%+replace%` had dropped `base_family` from the root `text` element → cairo_pdf fell back to **Helvetica** in PDFs (PNGs were fine). Fixed at source in `scripts/00b_figure_theme.R`: (1) `family=base_family` on root text; (2) `update_geom_defaults("text"/"label"/repel, family=FIG_FONT)` for in-plot text geoms; (3) explicit `family=FIG_FONT` on the two `ggrepel::geom_text_repel` calls in script 37 (ggrepel `::`-called/not attached, so the global default didn't reach it). Re-ran all 8 figure scripts. **Verified via pdffonts: 14/14 manuscript PDFs pure Arial (0 Helvetica).**
- Minor cosmetic notes for an optional later polish pass (none affect readability): Fig6 subtype-legend crowds at MMMT/SCCOHT; Fig3B annotation sits near a point; Fig4 waterfall caption should specify "representative line OV2295 + pop-AF removal".
- **Phase-8 figure set COMPLETE:** 6 main (`fig1–6`) + 8 supp (`figs1–8`) in `docs/manuscript/figures/` (pdf+png), all branded via `00b`, pure Arial; `reports/assets/` refreshed (18 panels) so the HTML report restyles. Descriptor draft + explorer app complete (prior entries). All 4 Phase-8 agents idle.

### [2026-07-24] — Phase 8: Agent B figures COMPLETE + verified (Figs 2,3,6 + S1–4,7,8)
**Type:** Script / Verification · **Phase:** 8 · **Status:** Complete
- **fig-rna-prot-reuse** delivered scripts `34–37` (read `output/` only; analysis `00–20` untouched). Main: `docs/manuscript/figures/fig2` (RNA/proteome QC), `fig3` (subtype biology), `fig6` (ADC atlas). Supp: `figs1–figs4`, `figs7`, `figs8` (pdf+png). Refreshed 13 `reports/assets/` panels so the HTML report restyles. All Arial, locked palettes, tight density, panel tags, no narrative titles; PDFs embed Arial via cairo.
- **Orchestrator spot-checked fig2/fig3/fig6 at true size:** high quality, on-brand, numbers reconcile (pseudoalign median 91.1%→"91", concordance 0.40 w/ 0.72 ceiling + 0.38–0.48 benchmark band, PC1 20.7/PC2 10.4, site 0.2% of PC1, variance subtype≥site, 6856 complete, bridge r 0.991–0.994, 16/22 markers boxed). **ADC heatmap successfully rebranded** (cook_seq rust ramps, subtype top bar, square cells, hairline grid, DPEP3 NA grey) — replaces old magma+Dark2.
- Accepted Agent-B style calls: Fig3-D variance = median+IQR point-range (not violins; `output/` has only summary stats, no heavy recompute); Fig6 = absolute expression (not z; model-selection intent). Minor cosmetic notes for a later polish pass: Fig6 subtype-legend spacing near SCCOHT; Fig3B annotation near a point.
**Still running:** fig-overview-genomics — Fig 5 done; Fig 4 + S5/S6 in progress.

### [2026-07-24] — Phase 8: explorer app — PI-requested UI fixes
**Type:** Script · **Status:** Complete (verified in headless Chrome)
- `app/index.html`: (1) page background cream `#F0EEE9` → **white `#FFFFFF`** (PI prefers clean white; aligns with design-system note that warm tone is editorial-only, not default bg); (2) RNA/protein chart headers relabeled to separate the printed value from the bar encoding — **"value = TPM · bar = log₁₀ (per-gene/shared axis)"** (was the ambiguous "TPM, log₁₀ rel/abs"), since the number is absolute TPM while the bar length is log-scaled; (3) **renamed the app to "Gene Explorer: OvCAN Collection"** (h1 + `<title>`; eyebrow "Cook Lab · Ovarian Cancer Canada"; lede now defines the acronym) to reference the collection officially. Open question for PI: whether to also make white the default background in the branded HTML report (`build_report.py`) and other lab HTML outputs. Hosting: guidance provided (static → any free host; Netlify Drop / GitHub Pages / Cloudflare Pages / private Artifact); NOT deployed — flagged that a public URL releases the full unpublished RNA+protein matrices, so hosting should track the deposition/embargo plan.

### [2026-07-24] — Phase 8: descriptor draft + explorer app COMPLETE (verified)
**Type:** Report / Verification
**Phase:** 8 (descriptor authoring)
**Status:** Complete (2 of 4 Phase-8 agents done + independently verified)
**What was done:**
- **descriptor-draft** → `reports/05_scientific_data_descriptor_draft.md` (~5,700 w). Verified: all 9 *Scientific Data* sections in policy order, **no Results/Discussion**, reframing correct (biology-recovery framed as data-quality evidence). Cohort composition cross-checked against `line_family_map.csv` — **exact match**: HGS 24 lines/16 patients, CC 8/8, MC 3/3, EC 2/2, MMMT 2/2 (VOA5217, VOA5436), SCCOHT 2/2 (BIN67, COV434), LGS 1/1 → 42/34; all 5 families are HGS. All quantitative claims match report 01. Placeholders (deposition) + CHECK (RNA read config, WES capture kit) + NOTE (HRD, table renumber, deposition hygiene) flagged for PI. Real lit-review citations w/ DOIs; academic-prose applied (no em-dashes in body).
- **webapp-explorer** → `reports/04_webapp_feasibility.md` + `app/` (index.html 23 KB + data.js 3.4 MB gzip+base64 + build_payload.py + README). Fully client-side static, no server; 2.6 MB gzipped payload; in-browser `DecompressionStream`; 28,901 RNA symbols + 8,427 proteins × 32-line union. Verified: CLDN6 (ENSG00000184697) per-line TPM matches raw `rna_tpm.csv` EXACTLY (VOA8771 1072.2 / OV3331 892.1 / OV90 84.7); headless-Chrome render confirms Cook branding + subtype palette matching the figures + honest n.d./no-RNA states + ‡ label-conflict flags. Recommended host: GitHub Pages. NOT deployed (PI's call).
**Key outputs:** `reports/05_scientific_data_descriptor_draft.md`; `reports/04_webapp_feasibility.md`; `app/{index.html,data.js,build_payload.py,README.md}`.
**Decisions:** web app = client-side static (no Shiny/server); explorer reuses `subtype_colours` from `00b_figure_theme.R` for cross-deliverable cohesion.
**Issues:** none blocking. Web-app gene-ID caveats (1,835 multi-mapping symbols → max-mean-TPM; ~17% Ensembl IDs unsymboled dropped; protein block-missingness) surfaced in UI + doc.

### [2026-07-24] — Phase 8: figure re-engineering + draft + web-app orchestration
**Type:** Decision / Orchestration / Script
**Phase:** 8 (descriptor authoring)
**Status:** In progress (4 sub-agents running)
**What was done:**
- **Recon of branding:** read `~/Lab/Branding` tokens (`tokens.jsx` DATAVIZ, `colors_and_type.css`) + the `visualization` skill. Confirmed the updated default **diverging colormap = navy/slate-blue↔rust** `['#1E3A5F','#5A7AA3','#BCC9DC','#F8F6F1','#FABEA0','#EB6235','#C2410C']` (git 7b1f3dd) — this OVERRIDES the skill's `scico::vik` for this project. Probed R env: patchwork/cowplot/ggplotify/scico/ComplexHeatmap/circlize present; **only Arial/Helvetica/Inter fonts** available (Manrope/JetBrains missing); no UpSet pkg.
- **Built the shared foundation** `scripts/00b_figure_theme.R` (tested — renders, font=Arial, diverging maps −2→navy/0→pale/+2→rust): `theme_ovcan()` (dense, base 7–8), locked `subtype_colours`/`site_colours`/`tier_colours`/`present_colours`, `scale_*_subtype()`, `scale_fill_cook_div()` (midpoint-centered), `scale_fill_cook_seq()`, `cook_div_colfun()`/`cook_seq_colfun()` for ComplexHeatmap, `ht_opt_cook()`, `save_fig()`, column widths `W1/W15/W2`. Encodes a **palette decision table** (centered→diverging, non-centered→sequential, fixed-categorical→locked) + PI's rules: **Arial** (deliberate deviation from Inter brand font), **text minimalism** (no narrative titles; captions narrate), **max legible density**.
- **PI directives captured:** (a) figures laid out programmatically (patchwork/cowplot), re-optimized for density/legibility; (b) new navy↔rust diverging map everywhere; (c) rebrand package-default viz (oncoprint, CNV/ADC heatmaps); (d) cohesive palettes mindful of data type; (e) deposition = placeholders for now; (f) **Fig 1 schematic replaces Table 1**; (g) explore a gene-explorer web app (CLDN6-report precedent).
- **Wrote the orchestration plan** `reports/03_figures_and_webapp_plan.md` (agent-consumable shared spec: branding contract §A, figure architecture §B with disjoint script ownership, package-default reworks §C, draft §D, web-app §E).
- **Launched 4 parallel sub-agents** (dedicated NEW figure scripts read `output/`; analysis scripts 00–20 untouched → zero collision): fig-overview-genomics (scripts 30–33 → Figs 1/4/5, S5/S6), fig-rna-prot-reuse (scripts 34–37 → Figs 2/3/6, S1–S4/S7/S8), descriptor-draft (→ `reports/05_...draft.md`, uses `academic-prose`), webapp-explorer (→ `reports/04_webapp_feasibility.md` + `app/`).
**Key outputs (so far):** `scripts/00b_figure_theme.R`; `reports/03_figures_and_webapp_plan.md`. Agent outputs pending.
**Decisions:**
- Diverging navy↔rust overrides `scico::vik` for this project (PI branding update). Recorded in the theme module.
- Manuscript figures use **Arial**, not the Inter brand font (PI; journal norm + embedding). HTML report keeps Inter via CSS.
- Fig 1 (branded schematic + coverage matrix) replaces Table 1.
- Dedicated figure scripts (read `output/`) instead of editing analysis scripts — faster, no heavy re-runs, no cross-agent collisions.
**Issues:**
- No UpSet package installed → S4 uses a ggplot presence-matrix (avoid new dependency). Non-blocking.

### [2026-07-24] — Phase 8: Scientific Data manuscript outline + figure/table plan
**Type:** Decision / Report
**Phase:** 8 (descriptor authoring)
**Script:** —
**Status:** Complete (awaiting PI review)
**What was done:**
- Created `reports/02_manuscript_outline.md` — a full Data Descriptor outline built from `reports/01_...results.md`, `PROJECT_SPEC.md`, and the 19 existing `f_*` assets.
- **Core reframing** (documented as a mapping table): a Descriptor has NO Results/Discussion. F2 QC → **Technical Validation §1**; F3 "recapitulates biology" → **Technical Validation §2** (data-quality evidence, not discovery); F4 genetics/identity → **Technical Validation §3–4**; F4 TOV21G-MSI + F5 ADC atlas + within-HGSC strata → **Usage Notes** (reuse features). Every "we found" recast as "the data recover / reusers can select."
- **Figure plan:** 6 main (recommend adding a NEW **Fig 1** resource-overview: workflow schematic 🎨 + sample×assay coverage matrix 🔧 + subtype/patient bars; `PROJECT_SPEC` targeted 5 — noted the merge fallback) + 8 supplementary. Flagged which are ready ✅ / generatable-now 🔧 / need-raw-data 🔒 (S8 WES coverage), plus alternative viz options (Bland-Altman concordance, RNA UMAP, per-arm CNV lollipop).
- **Table plan:** Table 1 cohort inventory + Table 2 Data Records (main); S1–S6 supp (all already materialized in `output/`).
- **Gating items surfaced:** raw-data deposition availability (WES FASTQs? raw MS from Morin?) — determines what Data Records can claim; accessions; the standing HRD/mucinous/STR/authorship decisions; landscape-citation lock.
**Key outputs:**
- File: `reports/02_manuscript_outline.md`
**Decisions:**
- Add a resource-overview Fig 1 (near-mandatory for a Descriptor; currently missing): rationale recorded in the outline.
- Hold figure generation + prose drafting until PI signs off on the outline structure.
**Issues:**
- WES raw-read/BAM availability unconfirmed — if only processed MAF/CNV are archivable, Data Records is processed-only with a stated limitation. UNRESOLVED (PI/deposition).

### [2026-07-23 ~23:45] — Phase 7: five workstreams integrated; report + branded PDF regenerated
**Type:** Script / Report
**Status:** Complete (renv pending, deferred to end)
**What was done:** All five parallel agents completed, verified their outputs, and returned report-ready deltas (`reports/reviews/deltas/*.md`). No file collisions (strict ownership held; shared MAF untouched). Integrated everything into `reports/01_multiomic_characterization_results.md` and regenerated the branded HTML/PDF (22 pp).
**Key results integrated:**
- **Patient-level genomics** (`wes-drivers`, scripts 07/08): TP53 **11/11 HGSC patients (100%)** survives collapse; CDK12 35% lines → **18% patients defensible** (3133 frameshift was ×4 pseudoreplication); somatic-confidence tiers (27 T1/11 T2/12 T3) replace `germline_like_vaf`; **defensible somatic BRCA2 = 0** (TOV81D/TOV3133D both Tier3); **FGA now autosome-restricted** (chrX was silently included — a sex artifact; TOV81D 0.073→0.021); CNV arm events robust to collapse (asymmetry stated). New: `wes_driver_freq_patient.csv`, `wes_driver_tiers.csv`, `wes_cnv_arm_freq_patient.csv`; regen `f_wes_oncoplot`/`f_wes_cnv` w/ family tracks.
- **TOV21G = candidate MSI-high/MMR-d OCCC** (`wes-signatures`, script 16): 1,416 coding = 6.9× median; SBS6/44/15 (not POLE); indel-rich; MLH1-methylation-type (WES-invisible). Reuse feature. **Deposition catch: archived MAF `NCBI_Build`=GRCh37 is WRONG — data are GRCh38.** `wes_mutation_load.csv`, `wes_msi_mmr.csv`, `wes_sbs_context.csv`, `f_wes_hypermutation`.
- **Confounder resolved** (`rna-variance`, script 17): joint model → **site adds ≤0.2% of PC1 beyond subtype**; CC (sole cross-site) 4–6%. variancePartition (lme4 fallback; pkg broken). Passage 83% collinear w/ site, not independent. Marker effect sizes (median AUC 0.69; SMARCA2 d=−2.8 vs SMARCA4 d=−0.7 backs post-transcriptional loss). `rna_pc_confounder_joint.csv`, `rna_variancepartition.csv`, `prot_variancepartition.csv`, `rna_passage_check.csv`, `rna_marker_effectsizes.csv`, `f_variance_partition`.
- **External identity** (`external-benchmark`, script 18): 5/42 in DepMap, each self-matches rank 1/67 (ρ 0.74–0.88); **COV434/BIN67 SMARCA4-damaging in DepMap** + Cellosaurus COV434 "misclassified→SCCOHT" flag → SCCOHT calls externally corroborated; ConsensusOV surfaced w/ TME caveat; **STR: 30/42 citable, 12 none** (11 VOA + TOV3121D). `external_ccle_concordance.csv`, `consensusov_calls.csv`, `cellosaurus_str_status.csv`, `f_external_concordance`.
- **Proteomics depth** (`prot-depth`, script 19): 0.40 concordance is a **compression ceiling** (protein spread 0.30× RNA; ~1.7× noise floor); **bridge is a PIS common-reference hub** (corrects reviewer's chain worry); 18.7% proteins block-missing per-plex. `prot_dynamic_range.csv`, `prot_block_missingness.csv`, `f_prot_compression`.
- **Consolidated supplement** (script 20): `output/supplement_per_line.csv` (42×41, Table S1 backbone).
**Infra:** `scripts/15_patient_family_map.R` (family map), `16–20`, `scripts/build_report.py` (md→branded-HTML→PDF, reproducible). Script 01 fixes: lazy biomaRt, header 104→105.
**Decisions:** "authenticate"→"corroborate/flag"; freqs per patient; somatic tiers; FGA autosome-only; TOV21G MSI reframed as a feature.
**Issues:** `variancePartition` broken in this R (findbars/reformulas) → lme4 REML fallback (documented). MAF build-header mislabel → deposition fix. Both non-blocking.

### [2026-07-23 ~23:00] — Phase 7: two independent peer reviews + comprehensive revision orchestration launched
**Type:** Decision / Orchestration
**Status:** In progress (5 agents running)
**What was done:**
- Ran two cold (context-free) Opus reviewers on the report: `sci-reviewer` (scientific quality) and `data-sci-reviewer` (code↔report fidelity). Full reviews: `reports/reviews/{scientific_review.md, data_science_review.md}`.
- **Fidelity verdict: HIGH** (~60 claims checked, **0 Critical / 0 Major** discrepancies; the Mutect2 FILTER re-filtering was independently re-derived from raw MAFs; HRD confirmed not-fabricated). **Scientific verdict: valuable + honest but MAJOR REVISION** — all numbers reconcile; issues are interpretive/statistical.
- **#1 scientific issue = pseudoreplication:** genomic frequencies were per-LINE not per-PATIENT (same-patient sublines carry identical variants: 3133→Q192*+CDK12, 2295→I195T, 1369→G244C). TP53 100% and CNV freqs survive collapse; point-mutation freqs (e.g. "CDK12 35%") do not.
- Built the shared dependency `metadata/line_family_map.csv` (`scripts/15_patient_family_map.R`): 42 generated lines → **34 independent patients** (5 multi-line families: 1369/2295/3133/3291/3121); **HGS-with-WES-MAF: 17 lines → 11 patients** (matches both reviewers).
- PI directive: execute **Tier A + Tier B + relevant extras**, go comprehensive (beyond Scientific Data length), `renv` LAST. Launched 5 file-partitioned sub-agents (see Status). Report/log integration held single-threaded to avoid collisions.
**Key outputs so far:** `reports/reviews/{scientific_review,data_science_review}.md`; `metadata/line_family_map.csv`; `scripts/15_patient_family_map.R`. Agent deltas will land in `reports/reviews/deltas/`.
**Decisions:** "authenticate" → "molecularly corroborate/flag" (STR not done — journal requires the statement); somatic-confidence tiering to replace the uninformative `germline_like_vaf` flag; drop chrX from the CNV figure; report freqs at patient level with family tracks.

### [2026-07-23 ~22:00] — Phase 6: consolidated results report compiled
**Type:** Report
**Status:** Complete
**What was done:** Wrote `reports/01_multiomic_characterization_results.md` — consolidated Phase 1–5 results organized under the proposed 5 descriptor figures (F1 cohort/provenance; F2 QC; F3 recapitulates-known-biology; F4 genomics+authentication; F5 usage/ADC+heterogeneity). Every claim backed by a stat + figure reference; honest caveats; outstanding-items list. This is the descriptor backbone.
**Decisions applied (defaults, flagged for PI):** genomic HRD documented as not-feasible/pending (needs BAMs); mucinous VOA8762/8771 flag-and-keep with GI-origin caveat.
**Outputs:** `reports/01_multiomic_characterization_results.{md,html,pdf}` — branded HTML→PDF via headless Chrome (Cook Lab report template + tokens), 13 pp Letter; figures rasterized (200 dpi) to `reports/assets/`. Branding + full-width figure legibility verified (masthead/stats/callout/tables + oncoplot/CNV).

### [2026-07-23 ~21:45] — Phase 4 FINALIZED + cross-verified (mucinous detail in)
**Type:** Script / Decision
**Status:** Complete (both authentication agents independently agree)
**Finalized findings (scripts 10–11; output/auth_mucinous.csv; figs/10–11):**
- **Mucinous authenticity — a real catch:** TOV2414 genuinely ovarian (KRAS G12A, CK7+/PAX8+/MUC5AC+/SATB2−; Sauriol 2020). **VOA8762 & VOA8771: ovarian origin NOT supported — intestinal-leaning (CK7-low, PAX8-low, CDX2+); SATB2 low so not clearly colorectal → intestinal/upper-GI-NOS. VOA8771 is the most GI-leaning** (CDX2≈41 TPM; the protein-z softens it only because TMT under-quantifies TFs). → **genuine ovarian mucinous may be n=1 (TOV2414)** pending BC Cancer STR + IHC (CK7/SATB2/PAX8/CDX2) + KRAS/APC/SMAD4.
- **OV90 refinement:** soften from clean 'HGS' → **'HGS-family carcinoma, serous identity NOT confirmed by expression'** (PAX8/WT1/SOX17/CK7 all ≈0 — lowest serous of the HGS lines; SMAD4-nonsense + BRAF; lowest HGS FGA). OV3331 = clean HGS. (Explains the thesis-era "OV90 resembles no subtype" flag.)
- SWI/SNF: COV434 & BIN67 lack WES → their SMARCA4 loss rests on protein alone → **SMARCA4 IHC is the confirmatory next test**.
- **Consolidated STR/IHC targets for BC Cancer/OVCARE:** VOA8762, VOA8771 (ovarian-vs-GI), VOA5436 (MMMT-vs-clear-cell), VOA4841 (CC atypical SMARCA4-loss), VOA4395 (EC markers), OV90 (serous identity).
**Infra:** first authentication agent died in the outage; authentication-2 completed scripts 10–11; the recovered original independently re-ran + corroborated (no overwrite). Cross-verified.

### [2026-07-23 ~21:40] — Phases 4 & 5 (authentication + integration) — COMPLETE / near-complete
**Type:** Script
**Status:** Phase 5 Complete; Phase 4 core complete (mucinous detail + figures finalizing)
**Infra note:** a transient internet outage (PI machine) killed the first agent pair mid-run (no partial writes); re-launched as `authentication-2`/`integration-2`; the originals recovered when the connection returned and coordinated via the shared task list (no duplication). Scripts 10–14 written.

**Phase 5 (integration) — COMPLETE:**
- **RNA–protein concordance** (n=30 lines with both): per-line Spearman median **0.41** (IQR 0.36–0.44), per-gene median 0.40. MATCHES CPTAC ovarian (0.38/0.45) + cell-line resources (0.42–0.58); ceiling ~0.72. Replaces the thesis's 2-line figure. (`output/integ_rnaprot_cor*.csv`; `figs/12_*`)
- **ADC-target atlas** (RNA+protein by subtype; `figs/13_*` verified): MSLN→HGS ✓, ERBB2/HER2→CC/MC ✓; NOTABLE exception — FOLR1 highest in MC/CC not HGS in our lines (cell-line vs tumour difference; report honestly). Per-line table = reusable product (`output/adc_expression.csv`).
- **Within-HGSC heterogeneity** (n=15; Hallmark + PROGENy): 3 subgroups — Inflammatory/NF-kB-EMT, Low-signaling, Hypoxic-glycolytic. Descriptive model-selection example. (`output/hgs_heterogeneity.csv`; `figs/14_*`)

**Phase 4 (authentication) — core complete:**
- **Multi-omic SWI/SNF panel** (`output/auth_swisnf_panel.csv`): TOV112D (SMARCA4 WES-trunc L639X + protein-low; SMARCA2 RNA-low), COV434 (SMARCA4 protein-low; SMARCA2 RNA-low), BIN67 (SMARCA4 protein z=−2.0 despite retained mRNA; SMARCA2 RNA rank 1) → all SWI/SNF-deficient, confirming Karnezis-2021. Also surfaced: TOV21G ARID1A-truncating (canonical clear-cell ✓), OV2085/OV2295 SMARCA2 variants, VOA4841 atypical SMARCA4 protein-loss (verify).
- **Per-line histotype-consistency** (`output/auth_perline_table.csv`): most concordant. Resolutions/flags: **OV90 & OV3331 → HGS supported over 'Adenocarcinoma'** (TP53-mut + high CNV); **TOV112D → dedifferentiated carcinoma (SWI/SNF-null)**, discordant with 'EC'; **NEW: VOA5436 (labeled MMMT) expression resembles clear-cell (HNF1B/NAPSA-high) — flag for review**; VOA4841 (CC) atypical SWI/SNF — verify.
- **Mucinous** (script 11, finalizing): TOV2414 externally authenticated (Sauriol 2020, KRAS G12A); VOA8762/VOA8771 ovarian-vs-GI unresolved (no external provenance) → STR/IHC from BC Cancer. STR column present (empty) for later.
**Outputs:** `output/{integ_,adc_,hgs_,auth_}*`; `scripts/10–14`; `figs/10–14` (13 ADC atlas verified).
**Decisions:** RNA–protein concordance now benchmarked to CPTAC (normal); ADC atlas honest re FOLR1; authentication recovers the two published reclassifications + resolves OV90/OV3331 + surfaces VOA5436 as a new identity flag.

### [2026-07-23 ~19:15] — Phase 1 RNA known-biology validation — COMPLETE (rna-biology agent)
**Type:** Script
**Status:** Complete
**Scripts:** `03_rna_de_signatures.R`, `04_rna_markers_genesets.R`
**Key results:**
- **Marker recovery (primary validation): 16/22 canonical markers land in the expected subtype** (MC CDX2/TFF1/TFF3/KRT20/MUC5AC; CC HNF1B/NAPSA/SPP1/GCLC; HGS SOX17/WT1/MUC16). The 6 "misses" are coherent (PAX8 pan-Müllerian; ZEB1/CDH2 peak in EC via TOV112D-dediff; MUC2 off in 2D lines; GPX3 nonspecific).
- GO recovery (fgsea): HGS DNA-repair ✓, CC oxidative/glutathione ✓, SCCOHT cell-cycle ✓, MC glycan suggestive, MMMT EMT not recovered (n=2). Hallmark singscore corroborates (CC PI3K/AKT/MTOR; SCCOHT/EC MYC/E2F/G2M; MC fatty-acid).
- **CAVEAT on record:** one-vs-rest among uniformly-proliferating lines is stringent; HGS's top GO axis (OXPHOS/ribosome-biogenesis) is confounded with site (all 15 HGS = Mes-Masson) and was NOT "corrected" (unidentifiable). Marker recovery (lineage-specific, not a batch artifact) is the load-bearing evidence.
- SWI/SNF (feeds Phase 4): TOV112D dual-low SMARCA4(rank2)/SMARCA2(rank4) → dedifferentiated; COV434 dual-low → SCCOHT; BIN67 SMARCA2 rank1 but **SMARCA4 mRNA RETAINED** (SCCOHT SMARCA4 loss is post-transcriptional → protein/WES definitive; don't over-read BIN67 SMARCA4 RNA).
- CC/MC gene-symbol bug FIXED (stopifnot-enforced). Rare-subtype DE/GO descriptive only.
- **Cross-omic:** SWI/SNF deficiency in TOV112D & COV434 now supported by RNA + protein + WES (SMARCA4 mutation in TOV112D) — a triple-omic authentication result.
**Outputs:** `output/rna_de_*`, `rna_signatures_*`, `rna_de_gsea_*`, `rna_markers_summary.csv`, `rna_swisnf.csv`, `rna_geneset_scores.csv`; `figs/03_*`, `04_*`.

### [2026-07-23 ~19:00] — Phases 2 & 3 (proteomics + WES) — COMPLETE (parallel agents)
**Type:** Script
**Status:** Complete
**Scripts:** `05_proteomics_load_qc.R`, `06_proteomics_separation.R`; `07_wes_mutations.R`, `08_wes_cnv.R`, `09_wes_hrd.R`

**Proteomics (Phase 2) — QC gate PASSES:**
- 8,430 proteins / 146,830 peptides; 31 generated lines mapped via tmt.layout; 5 Carey LGS excluded. Missingness is PER-PLEX (structural). Presence≥50% filter → 7,734 proteins (+12.8% vs na.omit's 6,856).
- Bridge-replicate reproducibility **Pearson 0.99 / CV median 5.3%** → strong; largely satisfies Technical Validation without the raw pipeline (raw MS still needed for PRIDE).
- Separation subtype-driven but weaker than RNA (subtype R² > plex & site on all top PCs; residual plex on PC2 adj-R²=0.21). Silhouette: HGS 0.18 carries it; CC ~0 (not proteomically cohesive), SCCOHT 0.03; n=2 unreliable.
- Markers recover: WFDC2/MUC16 HGS-high; SMARCA4/SMARCA2 lowest in SCCOHT; EPCAM/KRT7 low in TOV112D. TF markers (PAX8/WT1/HNF1B) flat — MS under-quantifies low-abundance TFs.
- Flags: subtype↔site confounded by design; 82 isoDoping proteins flagged; don't over-claim CC protein separation.

**WES (Phase 3):**
- **ROOT CAUSE of the ATM/ATR/BRCA2 artifacts:** the archived analysis IGNORED the Mutect2 FILTER column (MAFs hold ALL calls, pre-flagged germline/PoN/common; OV2295 25,914 rows → 493 PASS). Fix = FILTER==PASS + gnomAD/1000G/ESP pop-AF>0.001.
- Positive control HOLDS: **TP53 100% (17/17 HGS)**; credible hotspots consistent within patient families. Artifacts GONE: **ATM 82→9%, ATR 77→5%, BRCA2 68→9%, POLE 27→0%**. Coherent: KRAS in CC/MC/EC, CTNNB1/ARID1A in CC/EC, SMARCA4 in TOV112D.
- CNV (CNVkit "new"): canonical HGSC events recovered — 3q26 89%, 8q24/MYC 78%, 19q12/CCNE1 72%, 20q 72%, 17p/TP53 78%, 13q/RB1 56%, 10q/PTEN 28%. FGA HGS 0.63 > CC 0.40 > MC 0.36 > EC 0.27 > LGS 0.07 (TOV81D quietest).
- **HRD: NOT FEASIBLE from archived data** (CNVkit total-CN only, no allele-specific/BAF; 0/23 used --vcf; no tumor recal BAMs archived; scarHRD dep `copynumber` also unavailable ≥Bioc 3.18). No score fabricated. Path: recover WES recal BAMs → Sequenza → scarHRD.
- Flags: residual germline persists (~206 coding/line vs ~50–80 real) → drivers DESCRIPTIVE + tumor-only caveat, NOT a mutation-burden metric; TOV21G TMB outlier 1416 (check MMR/MSI vs artifact); frequencies inflated by RELATED ISOLATES (2295 trio / 3133 quartet / 1369 pair) → note patient-family non-independence (TP53 control unaffected); systematic chrX gain likely pooled-normal sex artifact; capture-kit concordance for the 5 public normals unconfirmed (deposition).
**Outputs:** `output/prot_*`, `output/wes_*`; `figs/05_*,06_*,07_*,08_*` (07 oncoplot + 08 CNV landscape visually verified).
**Decisions:** genomic HRD deferred (infeasible from archived data); WES drivers descriptive with tumor-only + patient-family caveats; CNV is the sound genomic layer.

### [2026-07-23 ~15:15] — Phase 1 (RNA): load+QC + subtype separation — COMPLETE
**Type:** Script
**Status:** Complete
**Scripts:** `scripts/00_setup.R`, `scripts/01_rna_load_qc.R`, `scripts/02_rna_separation.R`
**What was done:**
- tx2gene **pinned to Ensembl 105** (dec2021 archive; nearest to the release-104 index era — 104 archive retired), cached to `output/tx2gene_ensembl_rel105.csv` — replaces the notebook's non-reproducible live biomaRt.
- `tximport` (kallisto, ignoreTxVersion) → 39,568 genes × 31 lines; `DESeqDataSetFromTximport` (counts + avgTxLength offset — the correct path, resolving the TPM-vs-counts ambiguity); filtered to **22,544** expressed genes; `vst` for separation.
- QC: pseudoalign median **91.1%** (85.8–93.1); ~**20,119** genes detected/line; site pattern (Mes-Masson 92.2% vs Huntsman 88.1% pseudoalign) minor, detected-gene counts similar.
- Separation: PCA (PC1 20.7%, PC2 10.4%), t-SNE (perplexity 5), silhouette, Spearman heatmap.
**Key results:**
- **Confounder check (R² per PC): subtype ≫ site on all top PCs** (PC1 0.73 vs 0.31; PC2 0.69 vs 0.21; PC3 0.58 vs 0.05) → separation is biology-driven, not batch. Visually, Mes-Masson CC lines (TOV21G/TOV3392D) co-cluster with Huntsman CC lines, not with Mes-Masson HGS — confirms biology > site.
- Silhouette: HGS 0.16 · CC 0.12 · MC 0.15 (modest); MMMT 0.74 · SCCOHT 0.82 (inflated, n=2); **EC −0.01** (the two "EC" lines do not co-cluster).
- **TOV112D clusters with SCCOHT (BIN67/COV434) on PCA** → independent support for the Karnezis-2021 reclassification (dedifferentiated / SWI-SNF-null); feeds Phase-4 authentication + Fig 4.
**Outputs:** `output/rna_{txi,dds,vst}.rds`, `rna_{tpm,counts,qc_metrics,silhouette,pc_confounder}.csv`, `rna_pca.rds`; `figs/01_rna_qc_detected.pdf`, `figs/02_rna_{pca_subtype,pca_site,tsne_subtype,spearman_heatmap}.pdf`.
**Decisions:** Ensembl 105 pin (104 retired); PCA on vst top-2000 HVG; formal DE deferred to next step (HGS-restricted / descriptive).

### [2026-07-23 ~14:45] — Phase 2→3: lit-review report, synthesis, anchor docs, PI decisions
**Type:** Decision / Script
**Status:** Complete

**What was done:**
- Wrote the literature-review report (`reports/lit_review/ovarian-cell-line-multiomic-resource-literature-review.md`; 8 themes, ~80 refs, gaps, methods) and the consolidated synthesis (`reports/00_synthesis_and_recommendations.md`; reframe + proposed 5-figure structure + 6-phase plan + decision list).
- Citation verification: 36 load-bearing citations resolved via OpenAlex (titles/authors/venues match; none retracted). Publisher-paywall limitation disclosed.
- Wrote `PROJECT_SPEC.md` and `ANALYSIS_PLAN.md`.

**Decisions (PI):**
- Framing: **lean Data Descriptor** + ADC-target expression atlas retained as a featured **usage** example.
- WES: **re-filter (gnomAD/PoN) + canonical drivers only** + genomic HRD via scarHRD; CNV as-is.
- LGSC/Carey: **EXTERNAL** (published Carey/OVCARE data, not generated here) → **removed from resource**; cite Shrestha 2021. Sample-sheet `provenance` update delegated to data-inventory.

### [2026-07-23 ~14:10] — Literature review: refinement round (Phase 2)
**Type:** Exploration
**Status:** Complete
**What was done:**
- Ran OpenAlex co-citation chaining (`references/cochain.py`) over 94 collected identifiers (86 DOIs + 8 PMIDs); 88/94 resolved (existence-verification signal). 3,358 backward candidates; 93 co-cited by ≥4 seeds. Triaged top-40 → ~9 foundational primary papers to add (Zhang 2016 CPTAC ovarian proteogenomics; Ahmed 2010 TP53-ubiquity; Kuo 2009 CC PIK3CA; Hess 2004 mucinous; Konecny 2014 + Helland 2011 subtype classifiers; Langdon 1988 + van den Berg-Bakker 1993 historical cell-line landmarks incl. COV434 "granulosa" origin; Ganzfried 2013 curatedOvarianData).
- Spawned 3 refinement agents: `lit-refine-novelty` (confirm no competing ovarian multi-omic cell-line resource; carcinosarcoma-model gap), `lit-refine-mucinous` (mucinous ovarian line authenticity — GI-contaminant risk for TOV2414/VOA8762/VOA8771), `lit-chain-explore` (web-ground + verify the 9 chained papers).
**Key outputs:** `reports/lit_review/_all_findings.md`, `_cochain_candidates.{md,json}`; raw agent findings `reports/lit_review/raw/0[1-6].md` (07–09 pending).

### [2026-07-23 ~14:05] — Literature review: 6-cluster parallel search (Phase 2)
**Type:** Exploration
**Status:** Complete (core clusters)
**What was done:** 6 web-grounded search agents (≥5 searches each; ~100 queries total) across: (1) ovarian cell-line resources + STR authentication; (2) CCLE/ProCan/proteomics + *Scientific Data* descriptor norms; (3) molecular subtyping + TME-confounding critique; (4) rare-subtype biology + SWI/SNF; (5) ADC targets; (6) genomic HRD + tumor-only WES methods.
**Key findings (feed synthesis):**
- **Subtyping:** TCGA/ConsensusOV subtypes are substantially TME-driven (Zhang 2019 IHC; Schwede 2020; Olbrecht 2021; Tanis 2026). Mesenchymal + immunoreactive defined by stroma/immune cells absent from pure lines → **drop/heavily caveat ConsensusOV on our lines** (inference, not directly cited).
- **RNA–protein correlation norm** ~0.4–0.6 (CCLE 0.48; ProCan 0.42; Jarnuczak [not Frejno] 0.58; ceiling ~0.72) → our 0.34–0.46 slightly low; **recompute properly**.
- **Descriptor templates:** Kalocsay 2023 (60 breast lines TMT, *Sci Data*) = QC/normalization/deposition template; LL-100 (Quentmeier 2019) for WES/RNA arms. **No ovarian multi-omic cell-line descriptor exists** → novelty.
- **Rare subtypes:** Karnezis 2021 reclassifies **COV434→SCCOHT** (matches our label) and **TOV112D→dedifferentiated carcinoma (SMARCA4/SMARCA2-null)** (conflicts with our "endometrioid"). EC–SCCOHT "convergence" is the established SWI/SNF-deficient concept, likely anchored on a mislabeled line → reframe as validation, not discovery.
- **ADC:** FOLR1→HGSC (mirvetuximab approved), HER2→clear cell (amplified in mucinous), MSLN→serous. NaPi2b/MSLN/DPEP3 had high expression but NEGATIVE pivotal trials → frame atlas as model-selection/hypothesis-generation (expression necessary, not sufficient). Add CDH6.
- **HRD/WES:** genomic HRD from WES → `scarHRD` (HRDetect/CHORD are WGS-only). Peng signature is expression, **not genomic HRD** (relabel). Tumor-only calling: ~224 private germline/sample, ~50–70% FDR (Halperin 2017); worst in ~100%-pure lines (UNMASC) → mechanism for inflated ATM/ATR/BRCA2. CNVkit reference must be **same capture kit** (Talevich 2016). Takamatsu 2024 (*Sci Data*): HRD measurable in ovarian lines but tracks resistance, not sensitivity.

### [2026-07-23 ~13:55] — DECISION: exclude Carey LGSC from analysis (PI)
**Type:** Decision
**Status:** Complete
**Decision:** Per PI, remove all Carey low-grade serous samples from analysis (8 RNA `LGSOC_P#` + 5 proteomics `VOA*`, incl. contested VOA10841). Kept in `samples.csv` documented (`analysis_include=N`), not deleted.
**Rationale:** PI judges the batch effect unrecoverable. Note (on record): the LGSOC RNA-seq is also corrupted at processing (FASTQ mis-merge → ~20% pseudoalignment), so that RNA "batch effect" is at least partly artifact. **Consequence:** LGS effectively drops from the resource (sole remainder = TOV81D, Mes-Masson, WES-only — PI to confirm keep/drop). Background framing must no longer claim LGS coverage.

### [2026-07-23 ~13:52] — metadata/samples.csv created (source of truth)
**Type:** Script/Data
**Status:** Complete
**What was done:** Built `metadata/samples.csv` from the verified inventory — 55 rows × 21 cols (incl. per-assay presence, passage-per-assay, TMT plex/channel, subtype_status conflict flags, analysis_include + exclusion_reason). 42 included / 13 excluded. Flags: OV90/OV3331 defaulted HGS (Adenocarcinoma conflict recorded — **PI to ratify**); TOV112D reclass conflict (Karnezis 2021); VOA10841 LGS/CC conflict.
**Key outputs:** `metadata/samples.csv`; generator `scratchpad/make_samplesheet.py`.

### [2026-07-23 ~13:50] — Phase 1: archive assessment (4 parallel agents)
**Type:** Exploration
**Status:** Complete
**What was done:** 4 agents produced `reports/assessment/`:
- `01_thesis_manuscript_synthesis.md` — narrative reads as discovery, not descriptor; "31×3 omics" over-claim; ATM/ATR/BRCA2 rates are no-matched-normal artifacts; manuscript dropped thesis caveats.
- `02_data_inventory.md` — verified inventory: RNA 39 samples (31 clean cell lines + 8 corrupted LGSOC), proteomics 36 lines (5× TMT-11plex, 8,430 proteins), WES 23 lines (CNVkit CN + Mutect2 SNV). **13 lines tri-omic** (all Mes-Masson). CNVkit reference = 5 public normal exomes (PRJNA339046); mutations tumor-only.
- `03_methods_notebooks.md` — 3 .Rmd only; proteomics processing + several matrices generated by unarchived code; hardcoded paths, live biomaRt, CC/MC signature Gene bug, one-vs-rest DESeq2 across unreplicated n=2 groups.
- `04_figures_narrative.md` — thesis (16 figs) vs manuscript (4 dense figs); manuscript's 4-figure spine is descriptor-appropriate; missing proteomics/WES QC, STR authentication, batch/source figures.
**Decisions:** WES CNV canonical = the single normal-referenced CNVkit set ("old/new" are just Mutect2 vs CNVkit modules, verified from command files, confirmed by 2 agents).

### [2026-07-23 ~13:40] — Session start / anchoring
**Type:** Decision
**Status:** Complete
**What was done:** No prior anchor docs existed (fresh project). Surveyed directory; loaded Cook Lab branding + analysis standards. Established 3-phase plan: assess archive → literature review → synthesis/plan. Created project scaffold (`scripts/`, `output/`, `figs/`, `metadata/`, `reports/`, `data/external/`, `docs/manuscript/figures/`, `shellscripts/`).

---

## Decisions for PI

**RESOLVED (2026-07-23):**
1. ✅ **Framing** — lean *Scientific Data* Data Descriptor + ADC-target atlas as a featured usage example.
2. ✅ **WES mutations** — re-filter (gnomAD/PoN) + canonical drivers only + genomic HRD via scarHRD; CNV as-is.
3. ✅ **Cohort / deposition** — deposit full *generated* dataset (flagged); LGSC/Carey removed as EXTERNAL (published Carey/OVCARE data; cite Shrestha 2021).

**Open / secondary (proceeding on defaults):**
4. **OV90/OV3331 canonical subtype** — HGS (default applied) vs. Adenocarcinoma.
5. **TOV112D** — verify Karnezis-2021 dedifferentiated reassignment from our SMARCA4/SMARCA2 RNA+protein; drops EC to n=1 if confirmed.
6. **TOV81D** — kept in generated deposit (CHUM WES); LGS not featured as a subtype.
7. **Proteomics provenance** — obtain TMT pipeline + per-sample passages + QC from Morin lab (GSC).
8. **VOA8762/VOA8771** — request STR + histotype IHC + mutation from BC Cancer/OVCARE (no external provenance).
9. **Provenance check** — confirm no other samples besides LGSC/Carey are external; confirm whether proteomics LGS channels were physically run in our TMT.
