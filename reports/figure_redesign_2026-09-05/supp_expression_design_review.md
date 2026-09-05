# Expression supplement design review

Revised 5 September 2026. Scope: `scripts/37_supp_rnaprot.R`, Supplementary Figures S1–S4, S7 and S8, and their external legends. The existing figures are preserved in the centrally created `before/` backup. This revision changes presentation only; it does not change quantitative analysis tables or the manuscript.

## Design decisions

All explanatory paragraphs have been removed from plot captions, annotations and heatmap footers. The former `OVCAN_FIG_PLAIN=0` override is gone. Complete legends, including methods, statistical units, missingness and inferential qualifications, are generated reproducibly in `legends_supp_expression.md`. The figure PDFs contain axes, sample labels, keys, compact group/sample-size labels, interval limits and panel letters only.

The shared lab palette is used without editing the theme: slate for neutral marks, rust for emphasis, and the common slate–white–rust scale for signed pathway values. The positive RNA correlation heatmap uses a light-to-rust sequential scale ending at brand rust; it does not use a diverging scale for positive data or a near-black high endpoint. Dendrograms are slate. Established histotype, centre, plex and patient colours remain consistent with the other resource figures. Arial is embedded in the PDFs; normal plot labels are 7–8 pt, with 11 pt panel letters.

| Figure | Previous size (in) | Revised size (in) | Main changes |
| --- | --- | --- | --- |
| S1 | 7.2 × 3.1 | 7.2 × 2.85 | Three aligned panels with larger labels; PCA keys moved below the data; one informative outlier label; RNA/protein silhouette key outside the plot. |
| S2 | 7.2 × 3.5 | 7.2 × 3.2 | Replaced the crowded cross-assay scatter with a paired passage plot. Every eligible model has a labelled row, and each pair is directly connected. Nested-model bar keys sit outside the data. |
| S3 | 6.9 × 7.3 | 7.2 × 5.65 | Larger 7 pt model labels, three compact annotation keys at right, short scale below, lighter continuous ramp and slate dendrograms. Removed the large explanatory footer. |
| S4 | 5.51 × 5.0 | 5.51 × 4.65 | Kept a natural 1.5-column width for this simple display. Presence tiles and count points are aligned across all 32 rows; exact counts appear once. Slate replaces the dark presence fill. |
| S7 | 7.2 × 5.4 | 7.2 × 5.1 | All 25 marker labels enlarged and grouped in the usual histotype order. The linear d range is −4 to 6; three upper interval bounds outside the range are printed with arrows. Oriented AUC remains 0–1, with patient-group sizes outside the data. One shared recovery-shape key replaces redundant legends. |
| S8 | 7.2 × 6.6 | 7.2 × 5.8 | Rebalanced scatter, proliferation and pathway panels. All 15 model names remain visible. The illustrative-group key sits between panel rows; the PROGENy scale stays beside its heatmap. Repelled labels have more space around the central models. |

## Scientific information retained

- S1 retains both views of the same 31-model protein PCA and both modality silhouette summaries. Its legend states that related sublines remain and that visual overlap does not remove batch concerns.
- S2 retains all three variance summaries for PC1–PC5 and all 13 models with paired RNA/WES passage records. Its legend distinguishes incremental total R² from residual-scaled partial R² and makes the descriptive, noncausal interpretation explicit.
- S3 retains every 31 × 31 correlation, the diagonal of 1, all 22,678 input genes, and centre, histotype and patient annotations. Its legend gives the observed range and warns that high global correlations alone do not authenticate stocks.
- S4 retains all 8,427 feature rows represented by 32 patterns, including 70 rows unquantified across all five plexes. Counts sum to the complete analysis matrix. Exact counts are shown only in the right panel, with positions on a logarithmic axis. The legend distinguishes missing quantification from absent expression.
- S7 retains all 25 checks, the 28-patient unit, bootstrap intervals, three large upper interval bounds, low-marker arrows, BH flags and the attainable-p-value flag. The external legend includes the abundance threshold in the recovery rule and states the limitations from very small groups and unadjusted centre effects.
- S8 retains all 15 HGS models, the 14 pathways and the fixed descriptive grouping. Its external legend states the 12-patient sensitivity result (ARI 0.189), dependence on the same RNA data, and the limitation that these are illustrative groups rather than established molecular subtypes.

## Verification

The final script was run against the frozen shared theme. The R log contains no missing-row, dropped-label, clipping or rendering warnings. Source-level checks guard the model, marker, patient and feature totals. The script writes only canonical figure files, four existing report image assets and the external legend; the standard setup script also refreshes its runtime package/session provenance files.

Each final PDF is one vector page and was rasterized with `pdftoppm -png -r 180` for inspection at readable resolution. The rendered pages are in `qa_supp_expression/`. Arial embedding and PDF dimensions were checked independently. The final PDF checksums and PNG dimensions are recorded in `supp_expression_artifacts.csv`.

| Page | Visual QA |
| --- | --- |
| S1 | All points, the labelled PC2 extreme, axis values, silhouettes and external keys are legible. No key covers data; no clipping or overlaps. |
| S2 | All 13 model rows and paired endpoints fit. Equal passage records may have coincident symbols. No model labels overlap; bar keys are clear and do not cover bars. |
| S3 | All row/column labels and annotation strips fit. Both dendrograms, all scale ticks and categorical keys have clear margins after increasing the compact canvas height. No footer prose remains. |
| S4 | All 32 pattern rows and count positions align, including the final zero-plex row. Count labels are printed once and fit within the plot. No truncation or omitted patterns. |
| S7 | All 25 markers and both intervals are visible. Arrow labels 9.5, 19.0 and 20.6 preserve the three upper bounds beyond the d axis. Group sizes, symbols, footnote marks and the recovery key are readable and unclipped. |
| S8 | All 15 names and proliferation rows fit. Central labels were separated with deterministic repulsion and a targeted label nudge. All pathway and model labels, three group strips and the symmetric colour key are readable. |

The first compact S3 revision clipped a dendrogram tip and the lowest scale labels; that render was rejected and the final height increased to 5.65 in. S8 label padding and the TOV3041G nudge were adjusted after inspection. No final QA issues remain.

## Rebuild

From the project root:

```sh
/usr/local/bin/Rscript scripts/37_supp_rnaprot.R
```

Canonical PDFs and 400 dpi PNGs are written to `docs/manuscript/figures/figs1`, `figs2`, `figs3`, `figs4`, `figs7` and `figs8`. Complete narrative legends are in `reports/figure_redesign_2026-09-05/legends_supp_expression.md`; they should travel with these figures in manuscript and supplementary-material exports.
