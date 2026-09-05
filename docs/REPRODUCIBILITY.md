# Checkout and reproduction

A clean clone supports manuscript review, inspection of processed results and cluster recovery immediately. A full raw-input rerun needs the archived inputs and large reference files listed in [data restoration](data/README.md), plus the analysis environment. Repository creation does not certify a fresh-cluster end-to-end rerun.

## Review and recovery without R

```bash
python3 scripts/check_checkout.py
```

This read-only standard-library check validates the handoff files, 42-model/34-patient metadata, all 23 historical WES BAM aliases, 322 path hints, copied command hashes and the processed release's 49 checksums. It neither scans the cluster nor executes analysis. Follow [the cluster task](cluster/CLAUDE_TASK.md) for discovery.

The current figures, manuscript Word file and browser payload are already committed. Open them directly when reviewing. The dated audit/design manifests document the local builds; some refer to preserved workstation snapshots intentionally omitted from Git. Use the processed release's own checksum file to check the data package in this clone.

## Analysis environment

The recorded environment is R 4.5.2 with Bioconductor 3.21. `renv.lock` records the July environment; per-script September session records under `output/session_info_*.txt` and `output/package_versions.csv` provide later execution evidence. Treat these as reproducibility records, and check package differences rather than assuming that a fresh restore has already been validated.

In an appropriate cluster software environment, after choosing a library location and provisioning the pinned inputs:

```bash
# Install renv separately if it is not available in the chosen R environment.
Rscript -e 'renv::restore(lockfile = "renv.lock", prompt = FALSE)'
```

Do not perform package installation or computational analysis on a login node contrary to the cluster's policies. Exact scheduling/module commands depend on the cluster and are not guessed here. WES recovery additionally has a vendored `vcf2maf` utility; see `scripts/vendor/README.md`. Python builders use pandas/numpy for the browser and release packaging, python-docx for Word, and pypdf for the figure bundle. Scientific figure rendering uses the R packages in the scripts; the shared theme chooses Arial or an available standard sans serif fallback. The shared PDF device uses Quartz on macOS and the Cairo R package elsewhere to retain fractional text positioning; current canonical PDFs were checked with Quartz. Check font embedding, glyphs and layout again if rebuilding with another platform/device. PNG export uses ragg.

## Paths and input restoration

Run from the clone root, or set `OVCAN_PROJ` explicitly. `OVCAN_DATA` selects a source-data directory outside the clone while preserving the archived internal structure. The default remains `judy_archive/data/`. No analysis should modify that source tree.

```bash
export OVCAN_PROJ="$PWD"
export OVCAN_DATA="/actual/managed/path/to/archived/data"
python3 scripts/check_checkout.py --inputs
```

Reference resources under `data/reference/` and external comparison inputs under `output/external/` have separate locations. Restore the exact omitted source files using their provenance records. The corrected RNA import uses `tx2gene_ensembl_rel93.csv`; never substitute the previous release-105 map or a live annotation query.

## Analysis and exports

After inputs and packages are available, run the dependency order, not numeric filename order:

```bash
bash scripts/run_all.sh --no-fetch
```

Without `--no-fetch`, the runner also invokes the external-data fetcher. Review its pinned resources and saved provenance before any download. The scripts overwrite their derived outputs and runtime logs; run in a branch or separate clone when evaluating a change. Do not run the pipeline simply to locate missing cluster files.

Figure-only rebuilding is available via `bash scripts/run_all.sh --figures`, but it still requires some archived inputs: Figure 2 reconstructs bridge differences from the original proteomics workbooks. The clone therefore includes ready-made figures without claiming every figure can be rebuilt from committed summaries alone.

After a reviewed analytical change, relevant packaging commands are:

```bash
python3 app/build_payload.py
python3 app/build_single.py
python3 scripts/build_release.py
python3 docs/manuscript/figures/build_figure_bundle.py
python3 docs/manuscript/v7/build_docx_v7.py
```

Recheck data counts and units, release checksums, manuscript values and all affected panel references. Render and inspect changed PDF/DOCX pages. Publication, repository visibility, hosted browser deployment and archival DOI creation are separate actions.

## Repository preparation changes

The Git import preserves quantitative results and original inputs. Three standalone R scripts formerly defaulted to the workstation's absolute home path; their defaults now use the working directory, consistent with the shared setup. `OVCAN_DATA` was added to the shared setup, and the WES feasibility inventory follows that selected source tree. Defaults retain the original local layout. The six malformed historical BAM hints were corrected from archived command tokens; this changes handoff metadata, not WES calls. No scientific computation was re-run for the repository import.
