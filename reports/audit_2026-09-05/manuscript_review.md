# Manuscript review on 5 September 2026

The revision is based on `docs/manuscript/v5/OvCAN_Scientific_Data_draft_v5.docx`, including its 34 comments and one tracked replacement. The Word document differs from its Markdown source only in the tracked replacement of “reprocessed” with “processed” in the abstract. This was accepted as the starting text, then clarified to distinguish new RNA/protein data from reprocessed exomes. The complete v5 directory was preserved.

## Deliverables

- `docs/manuscript/v6/OvCAN_Scientific_Data_draft_v6.md` — revised manuscript source.
- `docs/manuscript/v6/OvCAN_Scientific_Data_draft_v6.docx` — readable manuscript with four current main figures embedded for review.
- `docs/manuscript/v6/build_docx_v6.py` — reproducible DOCX builder using the Markdown and current main figure PNGs.
- `docs/manuscript/v6/author_confirmation.md` — source-specific culture/proteomics drafts and a consolidated author TODO list.
- `reports/audit_2026-09-05/manuscript_comment_responses.md` — response to every Word comment.
- `reports/audit_2026-09-05/v5_comments_anchored.json` — full original comments and exact text anchors.
- `reports/audit_2026-09-05/v5_word_vs_markdown.diff` — proof of the substantive difference between the two v5 sources.
- `reports/audit_2026-09-05/reference_verification.json` — original DOI metadata checked through Crossref and Europe PMC.

## Scientific corrections integrated

| Priority | Issue in v5 | Revision and evidence |
|---|---|---|
| High | RNA and proteomics described as reprocessed rather than newly generated | Corrected in abstract, background, and workflow description using the user's Word comment. |
| High | 22 mutation profiles versus 23 copy-number profiles | Recovered TOV3121D adds 158 retained coding candidates. Both exome analyses now describe 23 models. Cascade: 582,474 total records, 16,081 Mutect2 PASS, 15,995 population-filtered, and 6,194 retained coding/splice records. |
| High | Avoidable transcript-to-gene annotation mismatch | Exact Ensembl93 reference recovered and mapped, with all 185,299 targets retained. Recomputed gene counts are 39,733 total and 22,678 count-filtered. The old release105 caveat is removed because the processing problem was corrected. |
| High | Unsupported statement that stock testing was not performed | Removed; explicit TODO requests source laboratory STR/mycoplasma records. External-reference availability is described separately. |
| High | Ambiguous unit of patient-level inference | Expression uses representatives; mutation and CNV arm frequencies count each patient once if any related model has the event; patient FGA uses the mean of that patient's model values. These are now distinguished explicitly. |
| Medium | Mixed adjusted/unadjusted variance partition interpretation | Main display and prose use raw R-squared components. Unique histotype/centre/shared fractions are 0.424/0.002/0.310 across 31 models and 0.520/0.003/0.269 after refitting in 28 patient representatives. |
| Medium | Histotype marker testing did not preserve centre composition | Main result uses the new within-centre joint permutation: 15 of 25 checks, null mean 8.46, P = 0.0101. The unrestricted test remains a sensitivity result. |
| Medium | Obsolete RNA–protein and QC figures | Gene detection median 20,270; correlations 8,033 genes/30 models, median 0.400; patient sensitivity 7,969 genes/27 models, median 0.418. Spread panel contains 8,035 genes because two constant-RNA genes are ineligible for correlation. |
| Medium | Protein units too specifically described as log2 ratios | Revised to supplied log2 protein abundances normalised using a pooled internal standard. Per-feature baselines remain in the supplied matrix; exact upstream formula/scaling is an author TODO. |
| Medium | 14 genes described as the predefined candidate panel | Corrected to 19 selected genes; 14 is the number with retained candidates, not the prespecified list size. Tier descriptions explain functional prioritisation and retain the tumour-only caveat. |
| Medium | Broad copy-number interpretation and genomic overclaims | Retained established positive controls and removed broad non-HGSC comparisons. Clear cell TOV3392D has FGA 0.671, exceeding the HGSC model median. Corrected arm intersections and patient aggregation are described. |
| Medium | Missing capture provenance presented as established incompatibility | Accepted the user's same-kit statement. Included recovered normal accessions and 290,475 CNVkit target bins. Original kit name/version/vendor BED remain clearly labelled TODOs. |

The analysis reports describe the underlying code changes and sensitivity results in detail. This manuscript review integrates the corrected outputs without treating exploratory signatures, candidate-driver tiers, or expression matching as proof of clinical validity or complete stock authentication.

## Writing and presentation

The manuscript now leads with the data and their intended use, preserves the Data Descriptor section structure, and uses past tense for performed procedures and present tense for resource contents. Patient relationships are described without “donor-family” terminology. Serial commas are applied throughout. Working filenames are removed from the narrative; the release manifest supplies file-specific documentation.

Technical Validation now describes what each check measures and reports representative evidence. It avoids blanket assertions that the data “passed” unspecified thresholds, novelty claims, repeated defensive caveats, and causal explanations unsupported by the design. Small histotype groups and confounding are explained at their relevant points. TOV112D's historical group label remains visible with its published reassignment, rather than being silently treated as confirmed endometrioid carcinoma.

References use author/year and DOI in the text as requested. All 17 original DOI records were verified; omitted punctuation was restored in the Beaufort, Soneson, Garcia, and Talevich titles. Two references became unused after removing the external correlation benchmark comparison and were omitted from v6's bibliography. The verified McLaren et al. VEP methods reference was added when annotation provenance was recovered. The general Ghandi et al. CCLE reference was replaced with the precise DepMap 24Q4 Public version-1 dataset citation after all four primary local files matched official checksums and byte counts. The later reference-manager conversion remains an author task.

The [Scientific Data submission guidelines](https://www.nature.com/sdata/publish/submission-guidelines) were checked on 5 September 2026. The revised title is 55 characters, the abstract is 144 words, and the required Data Descriptor sections are present. Four main figures are retained. The manuscript foregrounds data generation, records, and technical validation. Exact prior use of these datasets and repository access remain author confirmations under the journal's [editorial policies](https://www.nature.com/sdata/policies/editorial-and-publishing-policies) and [data policies](https://www.nature.com/sdata/policies/data-policies).

## Author work still required

The remaining TODOs concern facts that cannot be established by editing: actual culture/harvest conditions; RNA extraction/integrity; acquisition-level proteomics methods and scaling; WES run provenance and vendor capture design; the remaining Cellosaurus release/retrieval documentation; STR/mycoplasma records; ethics and consent; author declarations; public data/code deposits; and the hosted browser's release version. They are recorded at the point of use and consolidated in `author_confirmation.md`.

Published culture and proteomics methods were researched rather than left as empty requests. The confirmation document supplies draft prose, model-specific derivation references, and the conflicts that require collaborator input. These protocols are explicitly provisional: evidence that a laboratory used a procedure in another paper does not prove it was used for the present samples.

## Verification

The final bounded provenance update resolved the DepMap DOI, input filenames, and checksums locally. Cellosaurus is described as the documented 23 July 2026 snapshot; only the unsaved database release and original retrieval documentation remain in M06.

The v6 DOCX was built from the revised Markdown and four final main figure PNGs with the bundled Python runtime. Rendering used the documents skill renderer and bundled LibreOffice override at `/Users/dpcook/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/soffice`. The original v5 source files were preserved.

Structural verification found two tables, four inline figures, no unresolved build tokens, and no decorative paragraph borders. The title has 55 characters; the abstract has 144 words. The build manifest records the Markdown, builder, and figure SHA-256 values, and all inputs matched it. Main text, including author TODOs but excluding the abstract, references, tables, and figure legends, contains 3,753 whitespace-delimited words.

The first render exposed an inherited coloured rule below the title and wrapping inside three Table 1 header words. Both were corrected. Figure legends were expanded only where needed to define abbreviations, colour/line encodings, plotted intervals, and denominators. Body text uses Arial 11 pt; references use 9.5 pt; legends use 10 pt; table cells use 9 pt. All figures retain their supplied aspect ratios.

Final page-by-page visual inspection was performed at the renderer's original 1547 × 2000 pixel resolution. The final text pages are free of clipping, missing glyphs, overlaps, orphaned headings, and unintended table splits:

| Page | Content inspected | Result |
|---|---|---|
| 1 | Title, abstract, background opening | Pass; title rule removed |
| 2 | Background conclusion, study design, culture TODO, RNA preparation | Pass |
| 3 | RNA quantification, proteomics methods and TODO | Pass |
| 4 | Exome methods, candidate tiers, copy-number methods | Pass |
| 5 | Copy-number definitions, validation analysis methods | Pass |
| 6 | External references, ethics/AI TODOs, data records | Pass |
| 7 | Data deposit TODO, sequencing/protein and expression validation | Pass |
| 8 | RNA–protein interpretation, genomic checks, external matching | Pass |
| 9 | Usage notes, availability, contributions | Pass |
| 10 | Author declarations and 16 references | Pass |
| 11 | Both tables and notes | Pass; all column headings fit |
| 12 | Four complete figure legends | Pass |
| 13 | Figure 1 resource overview | Pass; final patient/variant labels verified |
| 14 | Figure 2 quality assessment | Pass; axes, legends, and bridge panels readable |
| 15 | Figure 3 expression validation | Pass; patient-mean marker display and PC1 components verified |
| 16 | Figure 4 genomics | Pass; full figure and rotated model/arm labels fit |

After the final figure-label updates, the document was rebuilt and rendered to `reports/audit_2026-09-05/docx_render_v6_verified`. Pages 1–12 were byte-identical to their already inspected page PNGs; all four final figure pages were inspected again at original resolution. The page comparison is recorded in `v6_render_comparison.json`. The subsequent M06 provenance cleanup changed only pages 6 and 10. Both were inspected again at original resolution and passed; the other 14 page PNGs were byte-identical to the already inspected version. The final 16-page render is in `docx_render_v6_final_provenance`, with comparison evidence in `v6_provenance_render_comparison.json`. The final DOCX SHA-256 is `dee5d807b88eb3d38210e88972f97acc896aa419395f7b4335a7a16c4b1cd221`.

Render PDFs and page PNGs are retained in the audit directory as QA intermediates. They are not additional manuscript deliverables.
