# Genomic figure redesign and final PDF review

## Scope and outcome

Rebuilt scripts `31_fig4_genomics.R`, `32_fig5_rare.R` and `33_supp_genomics.R` for Figure 4, Figure 5, Supplementary Figure S5 and Supplementary Figure S6. Canonical PDFs and 400-dpi PNGs are in `docs/manuscript/figures/`; companion report assets were refreshed. Full methods, definitions and caveats are in `legends_genomics.md`. The earlier designs are preserved in the centrally created `before/` backup.

The builders no longer add explanatory paragraphs, captions or legacy `OVCAN_FIG_PLAIN=0` defaults. They use the shared slate/rust theme without modifying it. All scientific values are read from the validated analysis tables; no underlying analysis, variant filtering, patient assignment, CNV calling or classifier result was changed for this redesign.

## Design decisions

| Figure | Final size | Principal changes |
|---|---|---|
| Figure 4 | 7.2 × 8.2 in | Replaced log-scale bars with labelled point trajectories. External matching uses its full comparison distribution and two simple endpoint markers. The oncoprint uses marker shape for variant class and fill for priority tier, eliminating competing colour channels. Related patient groups are directly bracketed. A logarithmic burden track removes the old capped TOV21G bar. The right bars retain the Tier 1–2 patient denominator; all 39 assessed CNV arms remain visible. |
| Figure 5 | 7.2 × 6.9 in | An asymmetric layout puts a vertical, log-axis model burden plot beside the spectrum and two three-row comparisons. This removes 23 rotated labels and the large empty region created by one extreme bar. RNA/protein marker and SWI/SNF panels share one z-score legend. Ranks, WES tiers and missing-measurement crosses remain legible. The figure is 21% shorter than its previous 8.7-in form. |
| Supplementary S5 | 7.2 × 4.2 in | Removed oversized generic gene-locus callouts and the repeated coordinate ruler. All 23 model profiles and 300 genomic bins remain. Direct histotype row labels replace a redundant colour strip and its legend. A compact family strip, horizontal legends, readable chromosome numbers and a clearly labelled FGA track use the available width. The figure is 26% shorter than its previous 5.7-in form. |
| Supplementary S6 | 7.2 × 3.8 in | Replaced stacked bars and long class/cluster keys with a labelled 15 × 4 probability matrix, an aligned margin plot and a numbered cluster strip. Winning probabilities have an outline; the input-set-sensitive model has an asterisk. Probabilities are explicitly labelled uncalibrated. Cluster-name interpretation and instability are described externally. |

Text uses Arial, with regular, italic and bold variants embedded as needed. Data labels and axes are approximately 7–8.5 pt at final size; compact track ticks are 6.7 pt and panel letters 11 pt. Gene labels are italic where practical. Data marks use slate (`#64748B`), rust and the shared white-centered expression/CNV ramp; charcoal-like fills were removed. The indel-fraction key now ends at rust instead of the darkest brown of the longer sequential palette.

The standalone `f_external_concordance.png` report asset now shows the complete self/non-self comparison plot corresponding to Figure 4B, replacing its former five-by-five submatrix. Its appropriate explanation is Figure 4B's external legend.

## Final panel mappings

- Figure 4: A, filtering trajectories; B, DepMap self/non-self correlations; C, driver-candidate oncoprint with coding burden and patient counts; D, HGSC arm-level CNV frequencies.
- Figure 5: A, coding-candidate burden; B, TOV21G SBS-96 spectrum; C, signature-group cosine maxima; D, target/reference refit sensitivity; E, mucinous-model marker expression; F, SWI/SNF RNA/protein ranks and WES tier labels.
- Supplementary S5: one integrated CNV heatmap with a patient-family strip and FGA track.
- Supplementary S6: A, class-probability matrix; B, top-minus-second probability margin; C, illustrative RNA cluster number. Cluster key: 1, hypoxic/glycolytic; 2, low-signalling; 3, inflammatory/NF-κB/EMT.

## Validation and visual quality assurance

All three builders executed successfully. Each final PDF was rendered with `pdftoppm` and inspected directly rather than relying on the separately exported PNG. The final review renders are `fig4_pdf_review.png`, `fig5_pdf_review.png`, `figs5_pdf_review.png` and `figs6_pdf_review.png` in this directory.

The review corrected the following issues before final export: labels touching Figure 4A trajectories and insufficient headroom for its largest count; the narrow LGS strip label; Figure 5 assay-header and colour-key clearance; merged chromosome 19–22 labels in S5; annotation labels and the missing-segment key colliding with S5's bottom legend/FGA ticks; and half-clipped first/last margin points in S6. S6 now asserts identical y-axis ranges across its three panels so margin points align with the corresponding heatmap rows. The final rendered pages have no identified clipping, overlap, blank panel regions or embedded narrative paragraphs.

`genomics_figure_qa.json` records dimensions, source-table hashes, PDF hashes and embedded fonts. `pdffonts` confirms that all four final PDFs use embedded Arial subsets. `pdfinfo` confirms one page per figure. PDF text checks recovered the expected filtering counts and signature-refit values, and confirmed removal of the former explanatory paragraphs. Data-level checks retain 23 WES models, 6,194 coding candidates, 51 driver-variant records represented in 46 gene/model cells, 14 driver genes and 16 mutation patients. Builder assertions preserve the 11-patient HGSC CNV denominator, 2,417 TOV21G signature SNVs, 23 CNV models and 15 HGSC RNA models.

Script 33 reports installed-package warnings about the S4Vectors build version and deprecated `anyMissing()` calls in interval operations. These are unchanged library warnings, not missing observations, failed coordinate mapping or font substitutions. No additional analysis rerun was required for this presentation-only revision.
