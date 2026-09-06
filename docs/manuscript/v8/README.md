# Scientific Data manuscript v8

This draft integrates independently checked evidence from the WES cluster handoff. It preserves v7 and its figure refinements.

- [Manuscript source](OvCAN_Scientific_Data_draft_v8.md) and [Word draft](OvCAN_Scientific_Data_draft_v8.docx) · [PDF for review](OvCAN_Scientific_Data_draft_v8.pdf).
- [Author confirmations](author_confirmation.md): remaining laboratory, provenance and submission questions.
- [WES completion report](../../../reports/wes_completion_2026-09-05/WES_COMPLETION.md): verified changes and remaining work other than deposition.
- [Targeted cluster follow-up](../../cluster/recovery/2026-09-05/FOLLOWUP.md).

Changes include acquisition/platform details, the complete tumour-only variant filtering chain, 23-model sequencing QC, target-coordinate definitions and explicit limits of the target-supported copy-number profiles. Figure 2 adds WES depth/coverage panels; Figure 4 now identifies its intermediate count as MAF PASS. The 6,194 coding candidates, patient-level summaries and existing 49-file processed release remain unchanged.

The Word build uses `python-docx`; run `python3 docs/manuscript/v8/build_docx_v8.py` from the project environment. Render the resulting document and inspect every page after changes. Its build manifest records source/figure hashes. TODOs remain where laboratory or author records are still required; this is not yet a submission-ready manuscript.
