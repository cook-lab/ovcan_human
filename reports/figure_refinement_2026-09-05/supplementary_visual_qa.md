# Supplementary figure PDF QA, 5 September 2026

**Result: all eight supplementary figures pass visual and text-preservation checks. No corrective edits are requested.**

Reviewed the complete, freshly exported single-page PDFs `docs/manuscript/figures/figs1.pdf` through `figs8.pdf` using their Poppler-rendered PNGs in `tmp/figure_refinement_2026-09-05/final_pdf/`. This was a read-only review of the figures and their data; only this QA report was written.

| Figure | Visual checks and outcome |
| --- | --- |
| S1 | Both proteomic PCA panels, histotype and TMT-plex legends, silhouette labels/counts and modality key are readable. Axis titles are separated from tick labels; the VOA10816 label is intact. |
| S2 | Passage-variance bars and the complete passage-pair display remain legible. Superscript R² and ΔR² legend text is intact. Panel B model labels, axis title and modality key do not collide. |
| S3 | The complete 31-by-31 correlation heatmap, both dendrograms, rotated column labels, row labels, annotation-strip names, three categorical legends and numerical colour scale are visible. The long model names and right-side legend text fit without clipping or overlap. |
| S4 | Missingness-pattern rows align with their count display. All counts, including 6,855 and the rust 70, fit inside the export. Quantification legend, plex labels and log-scale axis title are clear. |
| S5 | The chromosome heatmap retains all model/histotype labels, chromosome numbers 1–22, aligned FGA bars and scale, patient-family key and missing-segment key. The log2-copy-ratio legend and ≤/≥ endpoint symbols are intact. Dense chromosome labels are distinguishable and the legends do not overlap the heatmap. |
| S6 | All classifier probabilities and model names are readable, including the OV2085 asterisk. Maximum-probability cell outlines, probability colour bar, margin display and RNA-cluster column are intact. The upper axis title clears its tick labels. |
| S7 | Marker labels, arrows, asterisks and dagger are visible. Effect-size intervals and annotated off-scale upper limits remain readable; those annotations are intentional, not clipped labels. AUC scale, patient counts and recovery-rule key fit. |
| S8 | All labelled model points, proliferation labels, illustrative-group key and complete pathway heatmap are readable. The heatmap's rotated model names and PROGENy colour scale fit; panel and axis titles are clearly separated from other labels. |

## Export and content checks

- `pdfinfo` identifies all eight files as one-page **Quartz R Device** exports. `pdffonts` confirms embedded Arial fonts; bold Arial is embedded where used. No missing-glyph boxes or replacement symbols were observed.
- Extracted each current PDF's text and the corresponding committed PDF at `ed01287` with `pdftotext -layout`. After Unicode NFKC normalisation, the complete multisets of words, numbers and punctuation are **identical for all eight figures**. This includes numerical summary labels, probability tables, counts, model names and plot keys. Text order was not treated as scientific content because it can change with PDF layout extraction.
- Compared all 255 entries in `tmp/figure_refinement_2026-09-05/before_data_sha256.json` against their current SHA-256 hashes. **253 are unchanged**. Only `output/package_versions.csv` and `output/session_info.txt` differ, as expected when figure builders source the shared R setup. No analytical tables, matrices, metadata or release records in that snapshot changed.

This check establishes preservation and legibility of the re-exported supplementary figures. It does not re-run or independently revalidate their underlying statistical analyses.
