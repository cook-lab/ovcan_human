# OvCAN human ovarian cancer multiomic resource

Analysis, validated data records and a working Scientific Data Data Descriptor for 42 ovarian cancer cell models from 34 patients. The repository contains current processed results, reproducible analysis/figure scripts, the manuscript and an evidence-based checklist for recovering missing WES records from the cluster.

**Current state:** manuscript **v9** incorporates both September 2026 WES handoffs, verified five-normal coverage and a target-only CNV correction. The analysis audit and figure redesign are preserved. Several specific provenance, laboratory-method, stock-testing and submission details still need records; see the [coverage update](reports/wes_cnv_coverage_2026-09-06/COVERAGE_UPDATE.md) and [initial WES completion report](reports/wes_completion_2026-09-05/WES_COMPLETION.md). This working repository is public by author choice; journal data deposition and an immutable code archive remain separate steps.

## Start here on the cluster

For the new molecular-analysis extension, give the agent [the reviewed execution plan](docs/cluster/molecular_extension_2026-09-06/EXECUTION_PLAN.md). The new cluster inventory resolves alignment availability; the plan corrects the proposed scripts and supplies explicit model paths, 30 variant-review requests, MSI/allele-specific CN pilot gates and conditional HRD-scar work. A portable overlay can be prepared with `python3 scripts/49_prepare_molecular_handoff.py` when the local extension files have not yet been published.

For remaining provenance retrieval only:

1. Read the [updated targeted follow-up](docs/cluster/recovery/2026-09-06/FOLLOWUP.md), then [the agent guide](AGENTS.md) and [current project status](docs/PROJECT_STATUS.md).
2. Give Claude Code [the cluster task prompt](docs/cluster/CLAUDE_TASK.md). Its detailed checklist is [WES recovery](docs/cluster/WES_RECOVERY.md).
3. Start with the [23-model checklist](reports/audit_2026-09-05/wes_cluster_models.csv), [322 recorded path hints](reports/audit_2026-09-05/wes_cluster_path_hints.csv) and [archived CNV commands](docs/cluster/evidence/cnvkit_commands.txt).
4. Run the lightweight checkout check from the clone root:

```bash
python3 scripts/check_checkout.py
```

The initial task is to locate existing records, preserve provenance and report what is available. It does not require R, re-running the analysis, downloading alignments or submitting new cluster jobs. The historical paths in this repository are clues, not verified current locations.

## Current results and manuscript

- [Expanded molecular characterization](reports/molecular_extension_2026-09-06/README.md): 52-locus relative CN screen, AKT2 and deletion leads, 19-target patient-aware expression, exact variant read review and SBS3 stability; separate from manuscript v9 and the release.
- [Clinical and molecular annotation exploration](reports/clinical_classification_2026-09-06/README.md): CCNE1 DNA/RNA/protein evidence, curated BRCA2 significance, HRD feasibility and targeted follow-up. These exploratory records are separate from manuscript v9 and the processed release.
- [Latest WES coverage and CNV correction](reports/wes_cnv_coverage_2026-09-06/COVERAGE_UPDATE.md)
- [Initial WES completion and remaining paper requirements](reports/wes_completion_2026-09-05/WES_COMPLETION.md)
- [Analysis audit and corrections](reports/audit_2026-09-05/AUDIT_REPORT.md)
- [Manuscript v9](docs/manuscript/v9/OvCAN_Scientific_Data_draft_v9.md) · [Word file](docs/manuscript/v9/OvCAN_Scientific_Data_draft_v9.docx)
- [Redesigned figures and complete legends](docs/manuscript/figures/README.md)
- [Figure design review](reports/figure_redesign_2026-09-05/FIGURE_REDESIGN.md)
- [Processed-data release and dictionaries](release/README.md)
- [Author confirmations](docs/manuscript/v9/author_confirmation.md)
- [Reproduction and input restoration](docs/REPRODUCIBILITY.md)

## Repository layout

| Path | Purpose |
| --- | --- |
| `scripts/` | Canonical analysis, figure and packaging builders; `run_all.sh` defines dependency order |
| `metadata/` | Model inclusion, assay availability, aliases, passages and patient families |
| `output/` | Current processed matrices, statistical summaries and provenance; large reconstructable/downloadable objects excluded |
| `release/` | The validated 49-file processed-data package, including its saved checksums |
| `docs/manuscript/v9/` | Current draft, author questions and Word builder |
| `docs/cluster/` | Prioritised recovery instructions, Claude task and copied source-command evidence |
| `docs/data/` | Inventory and storage boundary for inputs intentionally omitted from Git |
| `reports/audit_2026-09-05/` | Audit findings, validation records and historical path clues |
| `reports/figure_redesign_2026-09-05/` | Figure review decisions, source checks and export manifests |
| `app/` | Local data browser and builders; hosted deployment has not been updated by creating this repo |
| `judy_archive/notebooks/` | Three original notebooks, retained as historical provenance |

Original sequencing/analysis inputs under `judy_archive/data/`, large download caches, old result snapshots and render scratch files are intentionally outside Git. See [data inventory and restoration](docs/data/README.md). No input files were deleted when preparing the repository. A clone supports review and cluster discovery immediately; it does not contain everything needed for a raw-input rerun.

The earlier `ANALYSIS_LOG.md`, plans, notebooks, v5/v6 manuscripts and older reports are historical records. The September audit supersedes their quantitative results and interpretations where they disagree. Do not combine old RNA results with the corrected Ensembl 93 import or reinstate the former 22-model variant count.
