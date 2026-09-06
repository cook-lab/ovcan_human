# CNVkit 0.9.10 local CBS runtime

Prepared only for the authorised target-only resegmentation from existing CNR inputs. No alignments, normal reference construction, coverage calculation, HMM analysis or cluster jobs are involved.

Executable in the prepared working copy: `tmp/cnvkit_runtime_0.9.10/bin/cnvkit.py`. Create any project-local venv with the recipe below; its absolute path is not an algorithm requirement.

Python: same directory, `python` (3.12.14). Rscript: `/usr/local/bin/Rscript` (R 4.5.2, DNAcopy 1.82.0). The historical cluster DNAcopy version was not recovered, so original-CNR baseline replay is required before interpreting the target-only comparison.

The venv is isolated from global Python packages. [`scripts/requirements-cnvkit-0.9.10.txt`](../../scripts/requirements-cnvkit-0.9.10.txt) lists every installed package. Numerical packages are pinned to NumPy 1.26.4, pandas 2.2.3 and SciPy 1.13.1; NumPy 2 and pandas 3 were avoided because CNVkit 0.9.10 uses removed/changed APIs. The HMM dependency pomegranate 1.1.2 satisfies the package requirement and imports, but its HMM API was not validated and is not used by CBS.

## Checks completed

- `cnvkit.py version` returns 0.9.10.
- Native `segment --help` succeeds (available from the installed CLI).
- `pip check` reports no broken requirements.
- Rscript with CNVkit's `--no-restore --no-environ` flags loads DNAcopy 1.82.0.
- Eight installed segmentation/filtering/bin-grouping/core source files are byte-identical to official v0.9.10 files. SHA-256 values and source URLs are in `runtime_source_verification.json`.
- No patches or substitutions to CNVkit algorithms or package source were made.

## Native invocation

```sh
MPLCONFIGDIR="$PWD/tmp/cnvkit_runtime_probe_2026-09-06/mpl" \
  ./tmp/cnvkit_runtime_0.9.10/bin/cnvkit.py segment INPUT.cnr \
  --method cbs --threshold 0.0001 --drop-low-coverage \
  --drop-outliers 10 --rscript-path /usr/local/bin/Rscript \
  --processes 1 --output OUTPUT.cns
```

Do not supply `--smooth-cbs`, a VCF or recenter the CNR input for this replay. The coordinating analysis script owns input hashes, target filtering, baseline replay, output comparisons and downstream integration. No original files should be replaced.

## Version-pinned behavior

CNVkit first splits the input using `by_arm()` before low-depth filtering. It infers the largest internal chromosome gap of at least 100,000 bp, searching outside margins of `max(50, round(0.1*n_bins))`; this is not the separate cytoband arm definition used by the downstream scientific analysis.

Within each inferred arm: remove bins with log2<-15 or depth==0; remove rolling outliers using width 50, quantile **0.95**, factor 10; remove zero/missing weights. The code's 0.95 is authoritative despite the outlier docstring mentioning a 90th quantile. Native smoothing handles mirrored boundaries and Savitzky–Golay residuals, so a simple rolling-median substitute is not equivalent.

CBS receives starts converted to one-based coordinates and text-rounded data (`%.6g`). It uses `set.seed(0xA5EED)` separately for each arm and weighted `DNAcopy::segment(alpha=0.0001)`, without optional `smooth.CNA`. Remaining local DNAcopy defaults are nperm=10000, hybrid permutation method, min.width=2, kmax=25, nmin=200, eta=0.05, trim=0.025, undo.splits='none'. It recomputes segment means from the original filtered bin log2 values and weights, restores original first/last bin endpoints, and converts starts back to zero-based. Python `transfer_fields` adds aggregate annotations and stretches outer arm endpoints to the full original arm limits; it does not recenter means. A native CBS run preserves these details.

Original-CNR replay must precede the target-only comparison because removing all antitarget rows can change inferred arm boundaries as well as removing their influence. Full replay validation and all23 corrected exports are delegated to the genomics analysis agent; this runtime preparation itself did not perform segmentation.

## Recreate the tested environment

The lock describes a tested macOS arm64 / Python 3.12 environment, not a claim of testing every operating system. The original cluster R/DNAcopy dependency record remains unavailable. Use a separate environment and validate the original-CNR baseline if the Python/R platform or versions differ.

```sh
python3.12 -m venv .venv-cnvkit-0.9.10
.venv-cnvkit-0.9.10/bin/python -m pip install -r scripts/requirements-cnvkit-0.9.10.txt
.venv-cnvkit-0.9.10/bin/python -m pip check
MPLCONFIGDIR="$PWD/tmp/cnvkit-mpl" .venv-cnvkit-0.9.10/bin/cnvkit.py version
Rscript --no-restore --no-environ -e 'library(DNAcopy); packageVersion("DNAcopy")'
```

DNAcopy must already be present in the chosen R library; the Python lock does not install it. Pass the environment-specific CNVkit/Rscript paths to script 29 through its documented arguments. First-import matplotlib font caching can take around a minute on this workstation. Package installation and source probes did not alter global Python/R libraries, original CNR files or historical segment calls.
