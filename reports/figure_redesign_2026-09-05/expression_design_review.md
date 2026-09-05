# Expression figure redesign and QA

Scope: scripts `34_fig2_qc.R` and `35_fig3_biology.R`, canonical Figure 2/3 PDF and PNG outputs, and their directly generated report assets. Statistical output tables were not edited. The parent agent preserved previous figures and scripts in `before/` and provided the softened common theme.

## Figure 2

The original layout allocated six panels over three uneven bands, with four bridge facets competing for width in the bottom strip. The new composition uses a full-height bridge column at right and three balanced bands at left: matched RNA QC plots, spread/coverage summaries, and a wide precision-by-abundance plot. The four bridge facets now have larger labels, consistent limits, and repeated vertical tick labels. A shared centre key sits outside the RNA data. Both RNA plots use exactly the same vertical scale.

Protein coverage uses counts encoded by horizontal position on a logarithmic axis, with values printed alongside all six points. This preserves visibility of the 70-feature class and the 6,855-feature class in the same panel. The zero-plex class has a rust accent. The paired-spread panel uses a quieter slate/rust pair and only median labels; methods and interpretation are external. The CV panel uses a wide horizontal axis with all ten deciles and clearly separated keys.

A new consistency assertion caught a pre-existing bridge display error: the old points were bridge-minus-primary, but the deposited bias and limits were primary-minus-bridge. The new display uses **primary-minus-bridge**, explicitly labelled on the axis and in the legend. Plotted sample counts, mean differences and SDs are checked against the deposited table; the maximum permitted numeric discrepancy is 1e-10. The statistical table itself was preserved. All 29,349 paired protein observations enter those statistics; the display preserves the previous +/-1.6 limit and excludes 0.191% of points visually.

## Figure 3

The two duplicate PCA panels are combined into one larger view: histotype is encoded by fill and centre by shape. Separate keys sit below the data; the histotype key uses two rows so every label remains visible. This releases room for an explicit comparison of all-model and patient-representative PC1 commonality. Matched dot-plot facets place model-level and patient-level percentages in separate columns, making the small unique-centre increments legible alongside the larger components. Direct labels show percentages without narrative notes. Each component has one labelled point per facet, preventing adjacent model/patient values from colliding or touching the panel bounds.

Patient-mean and individual-model marker displays now occupy the full lower band and share row centres and group boundaries. All marker names appear once. Histotype patient counts are printed under the heatmap. The model-level points use light grey for other histotypes and coloured points for the intended histotype; fixed-seed vertical jitter reveals overlap without changing expression values. The heatmap uses the shared slate-white-rust scale, a compact colour key, and solid/dashed expected-histotype outlines. The centre-confounding, small-group and model/patient qualifications are in the external legend.

## Quality checks

- Rebuilt canonical PDF and 400-dpi PNG files from materialized analysis results.
- Inspected rasterizations of the actual PDFs made with `pdftoppm`, rather than relying on the separately exported PNGs.
- Revised the first exports to fix an overlong PCA key, colliding commonality labels, and tight coverage-count padding.
- Preserved validation of the reconstructed individual-model RNA marker matrix against saved histotype means and the patient means against the marker effect-size table.
- Confirmed bridge sample counts and the sign/mean/SD agreement numerically.
- Final log and font/geometry checks are stored alongside the PDF review renders.

Final panel mappings and complete publication legends are in `legends_expression.md`. Figures contain short panel titles, axes, keys and numerical labels only; no explanatory paragraphs depend on theme-level suppression.

## Final verification status

Both final rebuilds completed without warnings or errors. Figure 2 is 7.2 by 6.15 inches (Cairo rounds the PDF media box to 518 by 442 points); Figure 3 is 7.2 by 7.0 inches (518 by 504 points). Body labels are approximately 7-8 points, panel titles 8.5 points and panel letters 11 points. Both PDFs embed Arial and Arial Bold; no custom branding font is used. The final PDF renders were inspected after the final changes: keys and endpoint labels are complete, patient/model commonality values occupy separate facets, and all 26 marker rows align across the two lower panels.
