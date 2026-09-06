# Scientific Data manuscript v9

This draft incorporates the independently checked CNV coverage handoff received on 6 September 2026. It preserves v8 and the established figure designs, while updating Figure 4D and Supplementary Figure S5 for target-only CNV resegmentation.

- [Manuscript source](OvCAN_Scientific_Data_draft_v9.md) · [Word draft](OvCAN_Scientific_Data_draft_v9.docx) · [PDF for review](OvCAN_Scientific_Data_draft_v9.pdf).
- [Coverage update and remaining work](../../../reports/wes_cnv_coverage_2026-09-06/COVERAGE_UPDATE.md).
- [Author confirmations](author_confirmation.md) and [targeted cluster follow-up](../../cluster/recovery/2026-09-06/FOLLOWUP.md).
- [Responses to all 13 current Word comments](user_comment_responses.md), including a plain explanation of the PC1 variance partition.

The five reference-normal depths are now verified from the primary coverage tables. Methods and technical validation distinguish mean bin depth from per-base sequencing coverage, describe reference-target support, and identify the absence of independent off-target validation. The final CNV profiles exclude artifactual antitarget bins; the relative baseline sensitivity in OV1369-R2 is documented. Missing normal coverage files are no longer an outstanding request. Manual execution/alignment provenance and the separate laboratory/author questions remain open.

Build with `python3 docs/manuscript/v9/build_docx_v9.py` in an environment containing python-docx. Render and inspect every output page after changes. The build manifest pins the manuscript and four figure inputs.

The author's new-data and AI-wording confirmations are incorporated. Cellosaurus provenance is documented by the analysis team. Proposed culture and proteomics protocols appear directly in highlighted Methods for collaborator review; they remain provisional where the current samples' records are unavailable.
