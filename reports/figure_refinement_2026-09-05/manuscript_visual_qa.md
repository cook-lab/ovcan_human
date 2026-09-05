# Manuscript v7 visual QA

Reviewed 5 September 2026. Source: `docs/manuscript/v7/OvCAN_Scientific_Data_draft_v7.docx` (SHA-256 `da411b8ea08a967e1f14c363137691d39c8e950e9d68746a8c283a655f35e1ee`). Reviewed the existing bundled-runtime render in `tmp/figure_refinement_2026-09-05/docx_render/`; no additional rendering, Word edits or script changes were needed.

**Result: pass.** The PDF has **17 pages**. Every text page, **1–13 inclusive**, was opened and visually inspected at readable resolution (source PNGs: 1547 × 2002 px). Figure pages 14–17 are covered separately by the parent review.

| Pages | Inspection result |
| --- | --- |
| 1–3 | Title, abstract, background and opening methods fit within the page margins. Headings stay with their following text; paragraph continuations are normal. |
| 4–6 | Methods, provenance, ethics and data-record text are readable. No clipping, overlapping text, broken glyphs or isolated headings. |
| 7–9 | Technical validation, usage notes and availability statements retain clear spacing and normal pagination. |
| 10 | References fit on the page, with legible links and aligned hanging indents. |
| 11 | Both tables fit intact; all column headers, cells and table notes remain readable. No split rows or clipped text. |
| 12 | Figure legends 1–3 fit intact. The existing Figure 3d wording specifying means of marker log2(TPM + 1) fits without overflow or a new page break. |
| 13 | Figure 4's legend is intact. Remaining space is the normal end of the legend section before separately paginated figure pages. |

As an additional comparison, page-image SHA-256 hashes for pages 1–11 and 13 are identical to an earlier pre-import v7 render. Its page 12 predates the Figure 3d legend clarification, which was already present in Git commit `ed01287` before this figure refinement. A direct comparison of the current Word paragraphs and table cells against `ed01287` confirms identical manuscript text. All thirteen pages were visually checked regardless of hash identity.

No material layout issues or new pagination issues were found. Existing highlighted author-confirmation placeholders remain as intended in the working draft.
