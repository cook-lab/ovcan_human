# Figure layout and PDF typography refinement

This pass addresses the author's follow-up on Figure 1C, PDF axis-title lettering, Figure 3A legend spacing and Figure 3D marker outlines. It follows the initial September redesign; the analytical results and manuscript text are unchanged.

## Changes

- **Figure 1C:** replaced the wide assay tiles with square cells and rotated assay headings. Model names, histotype strips and patient brackets remain visible. Panels A and B now stack beside C, using the space released by the compact table. The whole figure is 7.2 × 4.0 inches, down from the previous nominal 7.2 × 5.55 inches (about 28% less area). All 42 models, 126 availability cells, 34 patients and 31/31/23 assay counts are preserved.
- **PDF lettering:** the original `grDevices::cairo_pdf` export on this workstation rounded small Arial glyph advances to approximately whole PDF points. In a controlled 8 pt `Per-gene Spearman r` label, the advance after `e` was approximately 4.00 pt; Quartz retains approximately 4.45 pt. This explains the visibly crowded/uneven axis titles in the previous exports. The problem was observable in glyph positioning, not missing font embedding: both versions embed Arial. The shared `figure_pdf()` helper now uses Quartz on macOS, with fractional character positioning and the requested font size. All 14 PDFs were re-exported, including the two direct ComplexHeatmap device calls. Text remains extractable and fonts remain embedded. PNGs still use ragg at 400 dpi.
- **Figure 3A:** the gap between Histotype and Centre legend blocks is 10 pt, with 2 pt key-row spacing inside Histotype. The block hierarchy is now visually clear.
- **Figure 3D:** removed the grey horizontal class separators from the heatmap. Solid and dashed expected-cell outlines remain intact. The corresponding panel E row separators and marker order are unchanged.

The revised main figures are embedded in manuscript v7. The combined PDF, individual-file ZIP and export manifest were regenerated from the canonical exports. The previous figures remain available in Git commit `ed012870b2caf0acf1d68c2dd7d20e550975a84a`.

## Validation

All 14 final PDFs were rendered with Poppler and visually inspected. Main figures 1–6 and all supplementary figures S1–S8 have readable labels and no newly clipped or overlapping elements. Figure 1C's fixed coordinate aspect guarantees equal cell width and height; all 42 model labels are extractable from its PDF. The actual Figure 3 PDF shows separate legend blocks and unobstructed marker-cell outlines.

[pdf_validation.csv](pdf_validation.csv) records page dimensions, embedded fonts, text-glyph counts and hashes. Every individual PDF has exactly one page, embedded fonts, extractable text and no PDF annotation objects. Glyph origins lie within page bounds. [font_advance_probe.csv](font_advance_probe.csv) records the controlled small-text comparison. The independent [supplementary review](supplementary_visual_qa.md) also confirms unchanged extracted word/number/punctuation multisets across all eight supplementary PDFs.

The Word document was rebuilt from unchanged Markdown, rendered with the bundled LibreOffice/Poppler workflow, and all 17 pages reviewed. The manuscript text and tables remain unchanged; only the embedded Figures 1 and 3 and the build manifest changed. See [the manuscript review](manuscript_visual_qa.md).

A SHA256 comparison of the pre-refinement `output/`, `metadata/` and `release/` trees confirmed that only the two normal runtime records (`output/package_versions.csv` and `output/session_info.txt`) were refreshed. No quantitative tables, metadata, release files or original inputs changed. The dependency-light checkout check and all 49 release checksums pass. Only figure builders and artifact packaging were run; no statistical analysis was rerun. A diagnostic isolated the S5/S6 builder warnings to the installed S4Vectors version notice and deprecated `anyMissing()` calls inside the existing genomic-interval stack; there were no font warnings or missing-glyph errors. [render_warnings.csv](render_warnings.csv) records their types/counts. No dependency versions were changed for this figure correction.

## Reproduction and platform limits

The canonical PDFs were checked on macOS with R 4.5.2, Quartz and Arial. On platforms without Quartz, `figure_pdf()` uses the Cairo R package, which also retained fractional advances in the local diagnostic. Device font metrics differ: the local Cairo-package probe used an effective 8.5 pt size for a requested 8 pt label. A Linux/cluster rendering has not been visually validated, so inspect rebuilt layouts and fonts before replacing publication exports. The R package is already represented in `renv.lock`; no packages or fonts were installed for this pass.

The underlying distinction between fractional metrics and device-grid rounding is described in the [Cairo font-options documentation](https://cairographics.org/manual/cairo-cairo-font-options-t.html). The specific numbers above are measurements from this project's diagnostic PDFs, not a general claim that every Cairo build behaves identically.

Run the figure builders with `bash scripts/run_all.sh --figures` when the documented inputs and environment are available, then `python3 docs/manuscript/figures/build_figure_bundle.py` and `python3 docs/manuscript/v7/build_docx_v7.py`. Preserve the separate legends. Render and inspect outputs again after any future layout/device change.
