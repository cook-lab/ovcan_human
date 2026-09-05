# OvCAN figure exports

The September 2026 redesign provides 14 clean scientific figure pages with standard sans serif type, rust and slate data colours, and panel-specific layouts. Explanatory prose is in [figure_legends.md](figure_legends.md), separate from the figures.

- [Combined figure PDF](OvCAN_figures_redesigned.pdf): one figure per page, with PDF bookmarks. Pages retain the dimensions of the individual exports.
- [Individual figure package](OvCAN_figure_exports.zip): each figure as vector PDF and 400 dpi PNG, plus legends and checksums.
- [Manuscript v7](../v7/OvCAN_Scientific_Data_draft_v7.docx): the four main figures, revised legends and updated panel citations.

Figures 1–4 are the main figures in the current Data Descriptor draft. Figures 5–6 and S1–S8 remain additional analysis figures; this export collection does not change their manuscript status or imply that exploratory analyses are independent validation.

| File | Content |
| --- | --- |
| fig1 | Resource overview and assay coverage |
| fig2 | RNA and protein quality assessment |
| fig3 | Expression patterns and assay agreement |
| fig4 | Genomic filtering, candidate alterations and external matching |
| fig5 | Exploratory rare-histotype and mutational features |
| fig6 | ADC-target abundance atlas |
| figs1 | Proteomic principal components and separation |
| figs2 | Passage-related exploratory comparisons |
| figs3 | RNA sample correlations |
| figs4 | Protein detection patterns across plexes |
| figs5 | Genome-wide copy-number profiles |
| figs6 | ConsensusOV outputs and exploratory RNA groups |
| figs7 | Patient-level marker checks |
| figs8 | Exploratory within-HGSC expression patterns |

Run the R figure scripts 30–37 from the project root, then run `build_figure_bundle.py` with a Python environment containing pypdf. The figure builders read the validated outputs and metadata. The package script checks that every PDF has one page and no PDF annotation objects; `figure_export_manifest.json` records dimensions and SHA256 checksums. The before-redesign snapshot and design review are in `reports/figure_redesign_2026-09-05` at the project root. The subsequent compact Figure 1 layout, Figure 3 legend/separator fixes and PDF typography checks are documented in [the figure refinement report](../../../reports/figure_refinement_2026-09-05/FIGURE_REFINEMENT.md).

PDFs use embedded standard sans serif fonts. The shared `figure_pdf()` helper uses Quartz on macOS to preserve fractional glyph positioning; elsewhere it uses the Cairo R package. The original `grDevices::cairo_pdf` exporter rounded small Arial glyph advances on the analysis workstation. The current canonical PDFs were rendered with Quartz; verify text and panel spacing after rendering with a different platform/device. PNGs continue to use ragg at 400 dpi.
