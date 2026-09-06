# Working in this repository

## Read first

1. `README.md`
2. `docs/PROJECT_STATUS.md`
3. For cluster work: `docs/cluster/recovery/2026-09-06/FOLLOWUP.md` first, then `docs/cluster/CLAUDE_TASK.md` and `docs/cluster/WES_RECOVERY.md` as historical search guidance
4. For analysis changes: `reports/audit_2026-09-05/AUDIT_REPORT.md` and the relevant domain audit

`docs/manuscript/v9/` is the current manuscript. Scripts under `scripts/` are canonical; numeric script order is not dependency order. Use `scripts/run_all.sh` to understand dependencies. Original notebooks and earlier reports are historical evidence.

## Cluster recovery

Start by locating existing files in the run/project directories supplied by the operator. Do not launch a whole-filesystem scan, run the full analysis, submit jobs, change permissions, delete files or copy large sequencing files as part of initial discovery. Inspect existing logs, commands, reports and manifests first. A later explicit operator request can authorise additional work.

Keep original run directories read-only. Record exact paths, symlink targets, sample/run IDs, sizes and checksums when available. A historical path is not proof that a file still exists. Report permission-denied and broken-symlink cases separately from files not found in searched locations. Do not calculate checksums of multi-terabyte trees merely to build the first inventory.

Write curated recovery findings and small nonsensitive evidence in a new dated directory under `docs/cluster/recovery/`, with original paths and hashes. Put bulk retrievals under `data/cluster_wes_retrieval/` or managed storage; this directory is ignored by Git. Do not commit credentials, signed download URLs, raw sequence/alignment files, clinical identifiers or unreviewed logs. The author explicitly chose public GitHub visibility. Preserve public visibility unless the author requests a change; deposition of additional data still requires author direction.

Use the acceptance checks in `WES_RECOVERY.md`. Return an inventory, an item-by-item status table and a short account of what remains missing. Never replace a missing result with a pipeline default or an assumed protocol.

## Scientific invariants

- The resource has 42 models from 34 patients; related sublines are not independent patients. RNA/protein each cover 31 models; 30 are paired; WES covers 23 models from 16 patients; 13 have all three modalities.
- TOV3121D variants were recovered from its archived annotated VCF. There are 6,194 retained coding candidates across 23 models. The earlier 22-model statement is obsolete.
- Variant tiers are prioritisation categories for tumour-only calls, not confirmation of somatic origin. Relative CNV profiles do not establish absolute ploidy, LOH or HRD.
- Current CNV summaries use script 29's target-only resegmentation, followed by script 08. Do not restore archived CNS profiles as canonical: four overlapping antitarget bins carry extreme ratios and inflated two chromosome 1 segments. The 6 September coverage report documents the correction and verifies all five reference-normal coverage pairs.
- The five CNV-reference exomes and the Mutect2 variant panel of normals are different resources. Same-kit compatibility of model/reference exomes is author-confirmed. The recovered design is SeqCap EZ Exome v3; see the September WES completion report for verified run/QC evidence and remaining targeted requests. MAF PASS includes a converter-added common_variant flag. CNV profiles are supported almost entirely by target bins; segment spans interpolate between targets.
- RNA quantification uses the matched Ensembl release-93 transcript map. Do not substitute live annotations or the old release-105 map.
- Protein values are supplied log2 normalised abundance, not simply zero-centred log ratios to the internal standard. Bridge plots use primary-minus-bridge differences.
- Preserve historical histotype labels in metadata while retaining their interpretive qualifications. Do not infer stock authenticity, MSI, protein loss or drug response from an exploratory plot.

## Editing and validation

Run `python3 scripts/check_checkout.py` for a dependency-light repository/release check. Read `docs/REPRODUCIBILITY.md` before executing analysis. `OVCAN_PROJ` selects the project root and `OVCAN_DATA` can select an external source-data directory with the archived internal layout. R setup refreshes version/session logs, so even sourcing it is not a strictly read-only check.

Keep original inputs unchanged. Make focused changes and record methods, denominators, transformations, limitations and affected artifacts. Rebuild only dependent results; do not alter data to make a figure cleaner. Figures use standard sans serif fonts and the shared rust/slate theme in `scripts/00b_figure_theme.R`. Narrative explanations belong in external legends. Render and inspect PDF/DOCX exports after changing their layout.

Do not add an open-source licence, assert completed ethics/stock tests, fill repository accessions or mark submission requirements resolved without supporting records and author authority.
