# =============================================================================
# Script: 00_setup.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Environment, paths, packages, and shared helpers for the re-analysis.
#          Sourced by every numbered script. Uses the system R installation.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# =============================================================================
# Anchors: PROJECT_SPEC.md, ANALYSIS_PLAN.md, metadata/samples.csv
# Reproducibility notes:
#   - Do NOT use live biomaRt for the transcript->gene map (the archived RNA
#     notebook did; it is non-reproducible). Use a pinned org.Hs.eg.db / EnsDb.
#   - Restrict formal DE to HGS or present descriptively (n=2 rare subtypes).
#   - Analysis set = provenance == "generated" & analysis_include == "Y".

# Paths -----------------------------------------------------------------------
PROJ  <- normalizePath(Sys.getenv("OVCAN_PROJ", unset = getwd()), mustWork = TRUE)
DATA  <- normalizePath(Sys.getenv("OVCAN_DATA", unset = file.path(PROJ, "judy_archive", "data")),
                       mustWork = FALSE)  # source data (read-only; may live outside the clone)
OUT   <- file.path(PROJ, "output")                  # reproducible intermediates
FIGS  <- file.path(PROJ, "figs")                    # exploratory figures
MSFIG <- file.path(PROJ, "docs", "manuscript", "figures")  # final paper figures
META  <- file.path(PROJ, "metadata")
for (d in c(OUT, FIGS, MSFIG)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

SAMPLE_SHEET <- file.path(META, "samples.csv")
stopifnot("samples.csv not found; set OVCAN_PROJ or run from project root" =
            file.exists(SAMPLE_SHEET))

# Packages --------------------------------------------------------------------
.required <- c(
  "tidyverse", "readxl", "matrixStats",                       # data wrangling
  "tximport", "DESeq2", "edgeR", "limma",                     # RNA quant / DE
  "AnnotationDbi", "org.Hs.eg.db",                            # pinned annotation
  "Rtsne", "cluster", "factoextra",                          # dim-reduction / clustering
  "ComplexHeatmap", "pheatmap", "ggrepel", "RColorBrewer", "viridisLite", # viz
  "singscore", "fgsea",                                       # gene-set scoring
  "maftools"                                                  # WES mutations
)
.optional <- c("progeny", "GenVisR", "scarHRD")               # install when needed (Phases 3/5)

check_pkgs <- function(load = FALSE) {
  ip   <- rownames(installed.packages())
  miss <- setdiff(.required, ip)
  if (length(miss))
    warning("Missing required packages: ", paste(miss, collapse = ", "),
            "\n  Install: BiocManager::install(c(",
            paste0('\"', miss, '\"', collapse = ", "), "))", call. = FALSE)
  optmiss <- setdiff(.optional, ip)
  if (length(optmiss))
    message("Optional not yet installed (Phase 3/5): ", paste(optmiss, collapse = ", "))
  if (load) invisible(lapply(intersect(.required, ip),
                             function(p) suppressPackageStartupMessages(
                               library(p, character.only = TRUE))))
  invisible(miss)
}

# Version recording ------------------------------------------------------------
# [review revision] Methods and Code Availability both claim this script "records
# paths, package versions, and the seed". check_pkgs() only WARNED about missing
# packages — nothing recorded a version, and there is no renv.lock or sessionInfo()
# anywhere in the pipeline. record_versions() closes that gap: it is called at the
# bottom of this file, so every script that sources 00_setup.R refreshes
#   output/package_versions.csv  (one row per analysis package + R/Bioc/platform)
#   output/session_info.txt      (sessionInfo() + seed + the calling script)
# session_info.txt reflects the LAST script sourced (it is rewritten each time);
# package_versions.csv is script-independent. Neither is a substitute for a
# lockfile — record renv::snapshot() before submission — but they make the
# environment recoverable, which nothing previously did.
.version_pkgs <- unique(c(
  .required, .optional,
  # WES / signature / external-benchmark stack (used by 07-11, 16, 18, 22)
  "data.table", "digest", "circlize", "GenomicRanges", "GenomeInfoDb", "Biostrings",
  "MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg38", "consensusOV", "jsonlite",
  # RNA / proteomics / variance stack (01-06, 12-14, 17, 19)
  "patchwork", "ggrepel", "lme4", "variancePartition", "Rtsne", "cluster",
  "ragg", "remotes", "BiocManager"))

record_versions <- function(quiet = TRUE) {
  ver <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                              error = function(e) NA_character_)
  pv <- data.frame(package = .version_pkgs,
                   version = vapply(.version_pkgs, ver, character(1)),
                   installed = !is.na(vapply(.version_pkgs, ver, character(1))),
                   row.names = NULL, stringsAsFactors = FALSE)
  pv <- pv[order(pv$package), ]
  # R / Bioconductor / platform go in the same table so one file answers "what ran"
  meta <- data.frame(
    package = c("R", "Bioconductor", "platform", "OS", "seed", "recorded_utc"),
    version = c(as.character(getRversion()),
                tryCatch(as.character(BiocManager::version()), error = function(e) NA_character_),
                R.version$platform, utils::sessionInfo()$running,
                as.character(SEED), format(Sys.time(), tz = "UTC", usetz = TRUE)),
    installed = TRUE, stringsAsFactors = FALSE)
  utils::write.csv(rbind(meta, pv), file.path(OUT, "package_versions.csv"),
                   row.names = FALSE)

  # sessionInfo() captures whatever the CALLING script has attached at source time
  # (base + anything already loaded); the package table above is the complete record.
  caller <- {
    a <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(a)) basename(sub("^--file=", "", a[1])) else "interactive"
  }
  si <- c(sprintf("# OvCAN session record  |  written by 00_setup.R  |  caller: %s", caller),
          sprintf("# %s  |  R %s  |  seed %d  |  project %s",
                  format(Sys.time(), tz = "UTC", usetz = TRUE), getRversion(), SEED, PROJ),
          "# Rewritten on every source() of 00_setup.R; see output/package_versions.csv",
          "# for the full analysis-package version table (script-independent).", "",
          utils::capture.output(utils::sessionInfo()))
  writeLines(si, file.path(OUT, "session_info.txt"))
  if (!quiet) message("Recorded output/package_versions.csv + output/session_info.txt")
  invisible(pv)
}

# Patient/subline family map ----------------------------------------------------
# [review revision] DEPENDENCY-ORDER GUARD. metadata/line_family_map.csv is written
# by 15_patient_family_map.R but consumed by 07/08/16/20 — i.e. numeric script order
# is NOT dependency order, and a clean checkout only worked because the file was
# committed. Rather than renumber 15 (the manuscript, reports/ and docs/ all cite it
# as "script 15"), any consumer calls ensure_family_map() first: it regenerates the
# map from metadata/samples.csv if absent. 15 depends on nothing but samples.csv, so
# this is always safe and never overwrites an existing map.
ensure_family_map <- function() {
  f <- file.path(META, "line_family_map.csv")
  if (!file.exists(f)) {
    message("metadata/line_family_map.csv absent -> running 15_patient_family_map.R")
    src <- file.path(PROJ, "scripts", "15_patient_family_map.R")
    stopifnot("15_patient_family_map.R not found; cannot build the family map" =
                file.exists(src))
    local(source(src, local = TRUE))
    stopifnot("15_patient_family_map.R did not produce line_family_map.csv" =
                file.exists(f))
  }
  invisible(f)
}

# Protein abundance matrix + zero-plex contract ---------------------------------
# [integration revision] ONE LOADER for output/prot_abundance_matrix.csv, used by
# 10, 11, 12, 13 and 19. Two things made this necessary:
#   (a) 05_proteomics_load_qc.R now drops the 3 no-symbol search rows (8,430 ->
#       8,427) and names non-representative duplicate-symbol rows SYMBOL|UNIPROT.
#       Each consumer had its own ad-hoc guard against the OLD shape; 10 and 11
#       still carried a `filter(!is.na(protein), !duplicated(protein))` that is now
#       dead code, so the five scripts no longer agreed on what they were reading.
#   (b) 70 rows are ZERO-PLEX: identified in the search output but NA in all 31
#       lines. They are deliberately RETAINED in the deposited matrix (deleting a
#       row loses the identification), which means every consumer holds 70 rows
#       that carry no measurement. Nothing must compute a statistic from them.
# The flag deliberately lives in prot_qc.csv, NOT as an extra column in the matrix:
# the matrix is a deposited resource that reusers read as `as.matrix(x[, -1])`
# (which is exactly what all five scripts do), and a non-numeric column would break
# that contract for them and for us. Instead the zero-plex set is DERIVED from the
# matrix here and cross-checked against prot_qc.csv$zero_plex, so the two files
# cannot silently drift. It is returned as attr(m, "zero_plex").
read_prot_matrix <- function(drop_zero_plex = FALSE, quiet = FALSE) {
  f <- file.path(OUT, "prot_abundance_matrix.csv")
  stopifnot("output/prot_abundance_matrix.csv missing — run 05_proteomics_load_qc.R first" =
              file.exists(f))
  pa <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
  stopifnot("no-symbol protein rows should have been dropped in script 05" =
              !any(is.na(pa[[1]]) | pa[[1]] == ""),
            "protein row names must be unique (05 assigns SYMBOL|UNIPROT to non-representatives)" =
              !any(duplicated(pa[[1]])))
  m <- as.matrix(pa[, -1]); rownames(m) <- pa[[1]]
  zp <- rownames(m)[rowSums(!is.na(m)) == 0L]        # derived from the data itself
  # Cross-check against the flag 05 wrote. If these ever disagree, one of the two
  # files was regenerated without the other and the run must not continue.
  fq <- file.path(OUT, "prot_qc.csv")
  if (file.exists(fq)) {
    q <- readr::read_csv(fq, show_col_types = FALSE, progress = FALSE)
    if (all(c("row", "zero_plex") %in% names(q)))
      stopifnot("prot_qc.csv$zero_plex disagrees with the all-NA rows of prot_abundance_matrix.csv — regenerate 05" =
                  setequal(zp, q$row[q$zero_plex]))
  }
  attr(m, "zero_plex") <- zp
  if (drop_zero_plex) m <- m[!rownames(m) %in% zp, , drop = FALSE]
  if (!quiet)
    message(sprintf("Protein matrix: %d proteins x %d lines (%.1f%% NA) | %d zero-plex rows %s",
                    nrow(m), ncol(m), 100 * mean(is.na(m)), length(zp),
                    if (drop_zero_plex) "DROPPED" else "retained (carry no measurement)"))
  m
}

# Hypermutator rule -------------------------------------------------------------
# [review revision] ONE RULE, ONE IMPLEMENTATION. 07 previously used ">3x median"
# and 16 used "robust_z > 5 & fold > 3"; they agree on TOV21G but were two separate
# definitions of one manuscript-facing flag. Both now call this helper. Requiring
# BOTH a median-absolute-deviation outlier AND a >3x fold keeps a 22-line panel from
# flagging a merely-high line; TOV21G's robust_z (~21) is far beyond the next line
# (~3.6), so the call is insensitive to either cutoff.
HYPERMUT_ROBUST_Z <- 5      # (n - median) / MAD
HYPERMUT_FOLD     <- 3      # n / median
hypermutator_stats <- function(n_coding, cell_line = NULL) {
  med <- stats::median(n_coding); mad_n <- stats::mad(n_coding)
  fold <- n_coding / med
  rz   <- if (mad_n > 0) (n_coding - med) / mad_n else rep(NA_real_, length(n_coding))
  data.frame(cell_line = if (is.null(cell_line)) seq_along(n_coding) else cell_line,
             n_coding = n_coding, median_load = med, mad_load = mad_n,
             fold_over_median = fold, robust_z = rz,
             is_hypermutator = !is.na(rz) & rz > HYPERMUT_ROBUST_Z & fold > HYPERMUT_FOLD,
             stringsAsFactors = FALSE)
}

# Cook Lab branding (see ~/Lab/Branding; visualization skill) ------------------
COOK_RUST <- "#C2410C"; COOK_NAVY <- "#0F172A"
# subtype_colours <- c(...)  # define once the final subtype set is locked

# Reproducibility -------------------------------------------------------------
options(stringsAsFactors = FALSE)
SEED <- 1234                      # deterministic tSNE / k-means / sampling
set.seed(SEED)

# Record the environment (see record_versions() above) — this is what makes the
# Methods/Code-Availability claim that setup "records package versions" true.
try(record_versions(), silent = TRUE)

message("00_setup.R loaded | R ", getRversion(), " | project: ", PROJ,
        " | versions -> output/package_versions.csv")
