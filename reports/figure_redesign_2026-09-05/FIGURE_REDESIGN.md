# Figure redesign review

All 14 canonical figures were evaluated against the analysis outputs and current manuscript, redesigned, rebuilt and visually checked from their PDF exports. The main changes are clearer panel proportions, larger labels, quieter rust/slate contrasts and separation of figure artwork from explanatory prose. The four main figures are incorporated in manuscript v7 with revised panel references and legends.

## Deliverables

- `docs/manuscript/figures/OvCAN_figures_redesigned.pdf`: 14 bookmarked figure pages at their native dimensions.
- `docs/manuscript/figures/OvCAN_figure_exports.zip`: individual PDFs and 400 dpi PNGs, complete legends, index and checksums.
- `docs/manuscript/figures/figure_legends.md`: all 14 external legends, including methods and interpretive qualifications removed from the artwork.
- `docs/manuscript/v7/OvCAN_Scientific_Data_draft_v7.docx` and Markdown source: the current four-figure manuscript presentation.

The preceding figures, builders and v6 document records are preserved in `before/`; `before_manifest.json` records their checksums. The original v5 and v6 manuscript artifacts remain available. Figures 5–6 and S1–S8 retain their status as additional analysis figures; they were not added as new main figures or promoted to independent biological validation.

## Changes by figure

| Figure | Redesign |
| --- | --- |
| 1 | Replaced the 42-row narrow matrix with two aligned blocks and larger model names. Combined the identical variant/CNV model sets into one exome column after an explicit equality assertion. Moved compact model/patient counts beside the workflow. Labelled shared patients with brackets. |
| 2 | Gave the four bridge comparisons a full-height right column. Aligned the RNA QC scales, widened the abundance/precision plot and used log-position counts for protein coverage. Corrected the bridge display sign as described below. |
| 3 | Combined redundant PCA views using histotype fill and centre shape. Compared commonality components in paired small panels. Aligned patient-mean marker heatmaps and individual-model distributions across the lower band. Increased figure height where required for readable marker labels. |
| 4 | Used filtering trajectories and a compact reference-matching panel. Replaced competing variant-class colours with shapes, retaining slate priority-tier fills. Labelled patient groups directly and enlarged the oncoprint. |
| 5 | Put the model burden comparison beside spectrum and refit panels; aligned the expression and molecular-evidence panels underneath. Retained essential WES labels, ranks and missing-value symbols. |
| 6 | Tightened the paired target atlas and its lower panels. Used gene symbols consistently, retained absolute IQR values beside row-standardised heatmaps, and distinguished unquantified protein cells with crosses. |
| S1–S2 | Removed footnotes from the panels, enlarged PCA labels and replaced crowded passage scatter labels with a paired passage display. |
| S3–S4 | Moved correlation annotation keys to the side, improved dendrogram and scale margins, and aligned all 32 protein-presence patterns with a compact count display. |
| S5–S6 | Removed oversized locus callouts and repeated annotation keys, added a clear FGA label and separated chromosome labels. Replaced dense classifier bars with a probability matrix and aligned margin/cluster columns. |
| S7–S8 | Enlarged marker labels, showed out-of-view interval bounds explicitly, improved sample-label placement and tightened pathway heatmaps. |

## Scientific correction found during design review

**Medium priority, corrected:** the previous Figure 2 bridge plot showed bridge-minus-primary values while its stored mean differences and limits of agreement were primary-minus-bridge. The plotting direction and axis label now match the deposited summaries. New assertions compare each plotted pair's count, mean difference and SD with the stored agreement table at a tolerance of 1e-10. This fixes the display; the validated statistical table did not change.

Visual simplifications retain their scientific qualifications in the legends. In particular, priority tiers do not establish somatic origin; tumour-only signatures remain exploratory; relative copy-number profiles do not establish absolute ploidy or LOH; related models do not become independent patient replicates; and RNA/protein abundance differences do not establish assay mechanisms or drug response. The manuscript's pending provenance and cluster-retrieval items are carried forward in `v7/author_confirmation.md` and `v7/wes_cluster_retrieval.md`.

## Visual system and export checks

The common theme uses Arial and embedded Arial Bold/Italic where needed, with brand rust `#C2410C`, slate `#64748B`, light neutral fills and white pages. Signed heatmaps use a slate–white–rust scale. Text uses dark slate for contrast. Panel-specific categories retain stable mappings; variant classes now use shapes rather than a second fill palette. Explanatory paragraphs are no longer drawn by the plotting code. The old text-length suppression heuristic was removed so that it cannot silently discard a scientific label; the export helper suppresses captions only.

All 14 actual PDF exports were rasterized and inspected, with final source/render hashes checked after the last edits. Each canonical PDF has one page, embedded Arial-family fonts and no PDF annotation objects. Review iterations corrected clipping, overlapping numerical labels, crowded chromosome labels and a cut-off colour key. Sizes were chosen by content: 12 figures occupy less page area than before, while Figures 3 and S6 gain modest space for readability. `final_dimensions.json` records exact exported sizes and changes. PNGs are 400 dpi; PDFs retain vector artwork and text.

The 49-file processed release still passes its saved SHA256 checksums, and all 26 live output/metadata files matched to release copies are byte-identical. The browser data and source artifacts also match the preceding audit. Figure builders refresh runtime version/session logs as part of setup; no quantitative analysis tables were edited in this figure revision.

Full reviewer notes and validation records are in `expression_design_review.md`, `genomics_design_review.md`, `supp_expression_design_review.md`, `final_pdf_render/validation.json` and the preservation-check JSON files. Manuscript v7 was rendered separately to check the figure sizes and the updated legend pagination.
