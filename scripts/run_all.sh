#!/usr/bin/env bash
# =============================================================================
# run_all.sh — canonical execution order for the OvCAN human multi-omic resource
#
# Numeric script order is NOT dependency order. This file is the dependency
# order, derived from the inputs each script reads and the outputs it writes.
# The two places where the numbering misleads:
#
#   * 18_external_benchmarking.R must run BEFORE 10_authentication.R, because
#     10 joins the Cellosaurus reference-profile status that 18 parses. Neither
#     reads the other's outputs, so the ordering is free.
#   * 15_patient_family_map.R must run BEFORE the figure scripts, because it
#     writes metadata/line_family_map.csv, which 30_fig1_overview.R reads.
#
# Usage, from the project root:
#     bash scripts/run_all.sh              # everything
#     bash scripts/run_all.sh --figures    # figure scripts only
#     bash scripts/run_all.sh --no-fetch   # skip the external-data download
#     python3 app/build_payload.py        # refresh viewer data (pandas, numpy)
#     python3 app/build_single.py          # refresh standalone viewer
#     python3 scripts/build_release.py    # validate/package principal records
#
# Environment: R 4.5.2 / Bioconductor 3.21. Restore the recorded package
# versions first with  Rscript -e 'renv::restore(lockfile = "renv.lock")'.
# The matched Ensembl 93 transcript reference is pinned under data/reference,
# and the original data inputs are stored outside Git. See docs/REPRODUCIBILITY.md
# for restoring inputs before running; a clone alone is for review/recovery.
# =============================================================================

set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p output logs

FETCH=1
FIGURES_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --no-fetch)   FETCH=0 ;;
    --figures)    FIGURES_ONLY=1 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

run () {
  local script="scripts/$1"
  echo "── $1"
  Rscript "$script" 2>&1 | tee "logs/$(basename "$1" .R).log"
}

# -----------------------------------------------------------------------------
# Stage 0 — environment record and external data
# 00_setup.R and 00b_figure_theme.R are sourced by every downstream script;
# running 00 directly writes package_versions.csv and session_info.txt.
# -----------------------------------------------------------------------------
if [ "$FIGURES_ONLY" -eq 0 ]; then
  run 00_setup.R
  [ "$FETCH" -eq 1 ] && run fetch_external_data.R

  # ---------------------------------------------------------------------------
  # Stage 1 — per-layer load and QC (independent of one another)
  # ---------------------------------------------------------------------------
  run 01_rna_load_qc.R          # -> rna_*, tx2gene_matched.csv
  run 05_proteomics_load_qc.R   # -> prot_*
  run 07_wes_mutations.R        # -> wes_mutations_filtered.csv, wes_driver_tiers.csv
  run 08_wes_cnv.R              # -> wes_cnv_*
  run 09_wes_hrd.R              # -> wes_hrd_feasibility.md, wes_pipeline_parameters.csv
  run 15_patient_family_map.R   # -> metadata/line_family_map.csv

  # ---------------------------------------------------------------------------
  # Stage 2 — per-layer structure and derived features
  # ---------------------------------------------------------------------------
  run 02_rna_separation.R       # <- 01
  run 03_rna_de_signatures.R    # <- 01
  run 04_rna_markers_genesets.R # <- 01
  run 06_proteomics_separation.R  # <- 05, 02 (rna_silhouette.csv)
  run 16_wes_signatures_msi.R   # <- 07
  run 13_adc_atlas.R            # <- 01
  run 12_rna_protein_concordance.R  # <- 01, 05

  # ---------------------------------------------------------------------------
  # Stage 3 — cross-layer analysis
  # ---------------------------------------------------------------------------
  run 14_hgs_heterogeneity.R    # <- 01, 04 (rna_geneset_scores.csv)
  run 17_variance_confounders.R # <- 01, 02, 04, 05
  run 19_proteomics_dynamic_range.R  # <- 01, 05, 12, 13
  run 21_rna_sensitivity.R      # <- 01, 03
  run 22_wes_signature_refit.R  # <- 16

  # ---------------------------------------------------------------------------
  # Stage 4 — external comparison, then authentication
  # 18 BEFORE 10: 10 reads output/cellosaurus_str_status.csv.
  # ---------------------------------------------------------------------------
  run 18_external_benchmarking.R    # <- 01, 07, output/external/
  run 10_authentication.R           # <- 01, 05, 07, 08, 18
  run 11_mucinous_authenticity.R    # <- 01, 07

  # ---------------------------------------------------------------------------
  # Stage 5 — master per-model table (Table S1)
  # ---------------------------------------------------------------------------
  run 20_supplement_table.R     # <- 01, 05, 07, 08, 10, 16, 18
fi

# -----------------------------------------------------------------------------
# Stage 6 — figures. Part of the certified run, not a separate pass.
# -----------------------------------------------------------------------------
run 30_fig1_overview.R    # <- metadata/line_family_map.csv (15)
run 34_fig2_qc.R          # <- 01, 05, 19
run 35_fig3_biology.R     # <- 01, 02, 04, 12, 17
run 31_fig4_genomics.R    # <- 07, 08, 16, 18
run 32_fig5_rare.R        # <- 10, 11, 16, 22
run 36_fig6_adc.R         # <- 13, 19
run 33_supp_genomics.R    # <- 08, 14, 18
run 37_supp_rnaprot.R     # <- 01, 05, 06, 14, 17

echo
echo "Complete. Logs in logs/, results in output/, figures in figs/."
