#!/usr/bin/env Rscript
# =============================================================================
# Script: fetch_external_data.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Make output/external/ REPRODUCIBLE. 18_external_benchmarking.R depends on
#          ~650 MB of DepMap and Cellosaurus data that was downloaded by hand and by
#          no committed script — one of the five blockers to a clean re-run.
# Author:  Cook Lab (analyst: Claude)  |  Date: 2026-07-24
# =============================================================================
#
# *** DEFAULT MODE IS VERIFY, NOT DOWNLOAD ***
#   Run with no arguments and this script only CHECKS the files already on disk
#   against the manifest below (size, md5, row/column counts) and reports. It writes
#   nothing into output/external/ and cannot clobber the hand-downloaded snapshot.
#
#     Rscript scripts/fetch_external_data.R                 # verify (default)
#     Rscript scripts/fetch_external_data.R --download      # fetch anything MISSING
#     Rscript scripts/fetch_external_data.R --download --force   # re-fetch everything
#
#   `--download` skips any file that already exists unless `--force` is also given.
#   The manifest is the contract: a downloaded file that does not match its expected
#   md5 is reported as a MISMATCH (DepMap re-releases silently under the same names,
#   so a mismatch usually means the release moved, not that the download failed).
#
# WHAT IS FETCHED
#   A. DepMap Public 24Q4 (four primary files). DepMap distributes releases through
#      figshare; the article for 24Q4 is 27993248. Per-file download URLs are figshare
#      file IDs that change between releases, so this script resolves them through the
#      figshare REST API by FILE NAME rather than hard-coding IDs that would rot:
#        GET https://api.figshare.com/v2/articles/27993248/files
#      Portal landing page (human route): https://depmap.org/portal/data_page/?tab=allData
#   B. Cellosaurus records, via the public REST API (no key required):
#        GET https://api.cellosaurus.org/search/cell-line?q=<name>&format=json
#        GET https://api.cellosaurus.org/cell-line/<CVCL_xxxx>?format=json
#      One search JSON per line name in output/external/cellosaurus/names.txt.
#   C. Four DERIVED subsets used by 18. These are NOT downloads — they are computed
#      from A, so this script can regenerate them deterministically.
#
# STILL NOT REPRODUCIBLE FROM HERE (flagged for the PI, not solvable in code):
#   - The WES exome CAPTURE KIT behind intervals_sorted.bed (see
#     output/wes_pipeline_parameters.csv).
#   - The proteomics search parameters: judy_archive/data/proteomics/Readme.md is a
#     1-byte empty file, so the search engine, database, FDR and quantification
#     scheme are genuinely unavailable and must come from the Morin laboratory.
# =============================================================================

PROJ <- Sys.getenv("OVCAN_PROJ", unset = getwd())
EXT   <- file.path(PROJ, "output", "external")
CELLO <- file.path(EXT, "cellosaurus")

args     <- commandArgs(trailingOnly = TRUE)
DOWNLOAD <- "--download" %in% args
FORCE    <- "--force" %in% args
if (FORCE && !DOWNLOAD)
  stop("--force only makes sense together with --download")
message(sprintf("mode: %s%s | target: %s",
                if (DOWNLOAD) "DOWNLOAD" else "VERIFY ONLY (nothing will be written)",
                if (FORCE) " + FORCE (re-fetch existing files)" else "", EXT))

DEPMAP_RELEASE  <- "DepMap Public 24Q4"
DEPMAP_FIGSHARE <- "27993248"                       # figshare article id for 24Q4
FIGSHARE_API    <- sprintf("https://api.figshare.com/v2/articles/%s/files", DEPMAP_FIGSHARE)
CELLO_SEARCH    <- "https://api.cellosaurus.org/search/cell-line?q=%s&format=json"
CELLO_RECORD    <- "https://api.cellosaurus.org/cell-line/%s?format=json"

# ---------------------------------------------------------------------------
# MANIFEST — the contract. Values were read off the snapshot currently on disk
# (2026-07-23 download) and are asserted on every run.
#   kind: "download" = fetched from source; "derived" = computed from downloads
#   rows / cols: as counted by data.table::fread (rows EXCLUDE the header)
# ---------------------------------------------------------------------------
MANIFEST <- list(
  list(path = "Model.csv", kind = "download", source = "depmap",
       bytes = 645696L, md5 = "675210d17675f3517b0ce39a3c274f16",
       rows = 2105L, cols = 47L,
       desc = "DepMap cell-model metadata incl. RRID (Cellosaurus accession) and Oncotree annotation"),
  list(path = "OmicsExpressionProteinCodingGenesTPMLogp1.csv", kind = "download", source = "depmap",
       bytes = 506628654L, md5 = "71794802b750ce77c422dad0720a40af",
       rows = 1673L, cols = 19194L,
       desc = "log2(TPM+1) protein-coding gene expression; rows = models, cols = 'SYMBOL (entrez)'"),
  list(path = "OmicsSomaticMutationsMatrixHotspot.csv", kind = "download", source = "depmap",
       bytes = 4210723L, md5 = "a5aeb1deef897ead3c955c372148d840",
       rows = 1929L, cols = 543L,
       desc = "per-model hotspot-mutation count per gene"),
  list(path = "OmicsSomaticMutationsMatrixDamaging.csv", kind = "download", source = "depmap",
       bytes = 147655356L, md5 = "cb20fdbe1cf3b9b0d8ed4f53e1f399b6",
       rows = 1929L, cols = 19098L,
       desc = "per-model damaging-mutation count per gene"),
  list(path = "ovarian_ach.txt", kind = "derived", source = "Model.csv",
       bytes = 825L, md5 = "18afa3a0d81d5999c703a4bbd1d4daee",
       rows = 75L, cols = 1L,
       desc = "ModelIDs with OncotreeLineage == 'Ovary/Fallopian Tube' (one per line)"),
  list(path = "depmap_ovarian_models.csv", kind = "derived", source = "Model.csv",
       bytes = 6015L, md5 = "48425c6800aeb57043229309622a91dc",
       rows = 75L, cols = 5L,
       desc = "metadata for the ovarian models: ModelID, StrippedCellLineName, RRID, OncotreePrimaryDisease, OncotreeSubtype"),
  list(path = "depmap_expr_ovarian.csv", kind = "derived", source = "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
       bytes = 20926714L, md5 = "483f246a9f5f137ad4886143ba4cd1d3",
       rows = 67L, cols = 19194L,
       desc = "expression rows for the ovarian models present in the expression matrix (67 of 75)"),
  list(path = "depmap_expr_overlap5.csv", kind = "derived", source = "depmap_expr_ovarian.csv",
       bytes = 1818015L, md5 = "98bbfd5d1f0a1aebb1c5ba142675c5ea",
       rows = 5L, cols = 19194L,
       desc = "expression for the 5 models that overlap this resource (OV90, TOV21G, TOV112D, BIN67, COV434)"),
  list(path = "depmap_hotspot_overlap5.csv", kind = "derived", source = "OmicsSomaticMutationsMatrixHotspot.csv",
       bytes = 18327L, md5 = "fca40291552a55fa72bc847e381eeda7",
       rows = 5L, cols = 543L,
       desc = "hotspot counts for the 5 overlap models"),
  list(path = "depmap_damaging_overlap5.csv", kind = "derived", source = "OmicsSomaticMutationsMatrixDamaging.csv",
       bytes = 663680L, md5 = "f605adfefb01b61d91d22a84e46b75c0",
       rows = 5L, cols = 19098L,
       desc = "damaging counts for the 5 overlap models")
)
# Cellosaurus: one search JSON per line, plus the two accession records that were
# pulled directly. names.txt is the input list (42 generated lines).
CELLO_EXPECT_SEARCH_JSON <- 42L   # search_<name>.json, one per line in names.txt
CELLO_EXPECT_RECORD_JSON <- 2L    # cvcl_2010.json (COV434), cvcl_a1sr.json (TOV2414)
OVERLAP5 <- c("OV90", "TOV21G", "TOV112D", "BIN67", "COV434")

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

md5_of <- function(f) {
  if (has_pkg("digest")) unname(digest::digest(file = f, algo = "md5"))
  else unname(tools::md5sum(f))
}
dims_of <- function(f) {
  if (grepl("\\.txt$", f)) {
    n <- length(readLines(f, warn = FALSE)); return(c(rows = n, cols = 1L))
  }
  # nrows = -1 reads the whole file; these are the sizes we are asserting, so no
  # shortcut. The 507 MB expression matrix takes ~10-20 s.
  d <- data.table::fread(f, nrows = 0L)                 # header only -> ncol
  nl <- length(readLines(f, warn = FALSE)) - 1L         # rows excluding header
  c(rows = nl, cols = ncol(d))
}

verify_one <- function(m) {
  f <- file.path(EXT, m$path)
  if (!file.exists(f))
    return(data.frame(file = m$path, kind = m$kind, status = "MISSING",
                      detail = "not on disk", stringsAsFactors = FALSE))
  b <- as.integer(file.info(f)$size)
  h <- md5_of(f)
  d <- dims_of(f)
  probs <- c(
    if (!identical(b, m$bytes)) sprintf("bytes %d != %d", b, m$bytes),
    if (!identical(h, m$md5))   sprintf("md5 %s != %s", h, m$md5),
    if (!is.na(m$rows) && d[["rows"]] != m$rows) sprintf("rows %d != %d", d[["rows"]], m$rows),
    if (!is.na(m$cols) && d[["cols"]] != m$cols) sprintf("cols %d != %d", d[["cols"]], m$cols))
  data.frame(file = m$path, kind = m$kind,
             status = if (length(probs)) "MISMATCH" else "OK",
             detail = if (length(probs)) paste(probs, collapse = "; ")
                      else sprintf("%d rows x %s cols, %.1f MB, md5 %s",
                                   d[["rows"]], ifelse(is.na(m$cols), d[["cols"]], m$cols),
                                   b / 1e6, substr(h, 1, 8)),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# DOWNLOAD (only with --download; never overwrites unless --force)
# ---------------------------------------------------------------------------
figshare_urls <- function() {
  stopifnot("jsonlite is required to resolve figshare file URLs" = has_pkg("jsonlite"))
  fl <- jsonlite::fromJSON(FIGSHARE_API)
  setNames(fl$download_url, fl$name)
}

fetch_depmap <- function() {
  urls <- figshare_urls()
  for (m in MANIFEST) {
    if (m$kind != "download") next
    f <- file.path(EXT, m$path)
    if (file.exists(f) && !FORCE) { message("  skip (exists): ", m$path); next }
    u <- urls[[m$path]]
    if (is.null(u) || is.na(u))
      stop("figshare article ", DEPMAP_FIGSHARE, " has no file named '", m$path,
           "'. The release has probably moved; check ",
           "https://depmap.org/portal/data_page/?tab=allData and update ",
           "DEPMAP_FIGSHARE / MANIFEST together.")
    message("  downloading ", m$path, " <- ", u)
    utils::download.file(u, f, mode = "wb", quiet = FALSE)
  }
}

fetch_cellosaurus <- function() {
  nm_file <- file.path(CELLO, "names.txt")
  stopifnot("output/external/cellosaurus/names.txt is required (the 42 line names)" =
              file.exists(nm_file))
  nms <- trimws(readLines(nm_file, warn = FALSE)); nms <- nms[nzchar(nms)]
  for (n in nms) {
    f <- file.path(CELLO, paste0("search_", gsub("[/ ]", "_", n), ".json"))
    if (file.exists(f) && !FORCE) next
    u <- sprintf(CELLO_SEARCH, utils::URLencode(n, reserved = TRUE))
    message("  cellosaurus search: ", n)
    try(utils::download.file(u, f, mode = "wb", quiet = TRUE), silent = TRUE)
    Sys.sleep(0.5)                                     # be polite to the public API
  }
  for (acc in c("CVCL_2010", "CVCL_A1SR")) {           # the two records pulled directly
    f <- file.path(CELLO, paste0(tolower(acc), ".json"))
    if (file.exists(f) && !FORCE) next
    message("  cellosaurus record: ", acc)
    try(utils::download.file(sprintf(CELLO_RECORD, acc), f, mode = "wb", quiet = TRUE),
        silent = TRUE)
  }
}

# ---------------------------------------------------------------------------
# DERIVE (the four subsets 18 actually reads)
# ---------------------------------------------------------------------------
derive_subsets <- function() {
  model <- data.table::fread(file.path(EXT, "Model.csv"))
  ov <- model[OncotreeLineage == "Ovary/Fallopian Tube"]
  wr <- function(x, p) {
    f <- file.path(EXT, p)
    if (file.exists(f) && !FORCE) { message("  skip (exists): ", p); return(invisible()) }
    message("  writing derived ", p)
    data.table::fwrite(x, f)
  }
  wt <- function(v, p) {
    f <- file.path(EXT, p)
    if (file.exists(f) && !FORCE) { message("  skip (exists): ", p); return(invisible()) }
    message("  writing derived ", p); writeLines(v, f)
  }
  wt(ov$ModelID, "ovarian_ach.txt")
  wr(ov[, .(ModelID, StrippedCellLineName, RRID, OncotreePrimaryDisease, OncotreeSubtype)],
     "depmap_ovarian_models.csv")

  ach5 <- ov$ModelID[toupper(gsub("[^A-Za-z0-9]", "", ov$StrippedCellLineName)) %in%
                       toupper(gsub("[^A-Za-z0-9]", "", OVERLAP5))]
  sub_rows <- function(src, ids, out) {
    f <- file.path(EXT, out)
    if (file.exists(f) && !FORCE) { message("  skip (exists): ", out); return(invisible()) }
    d <- data.table::fread(file.path(EXT, src))
    setnames(d, 1, "ModelID")
    message("  writing derived ", out)
    data.table::fwrite(d[ModelID %in% ids], f)
  }
  sub_rows("OmicsExpressionProteinCodingGenesTPMLogp1.csv", ov$ModelID, "depmap_expr_ovarian.csv")
  sub_rows("depmap_expr_ovarian.csv",                       ach5,       "depmap_expr_overlap5.csv")
  sub_rows("OmicsSomaticMutationsMatrixHotspot.csv",        ach5,       "depmap_hotspot_overlap5.csv")
  sub_rows("OmicsSomaticMutationsMatrixDamaging.csv",       ach5,       "depmap_damaging_overlap5.csv")
}

# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------
if (DOWNLOAD) {
  if (!dir.exists(CELLO)) dir.create(CELLO, recursive = TRUE)
  message("\n-- DepMap (", DEPMAP_RELEASE, ", figshare ", DEPMAP_FIGSHARE, ") --")
  fetch_depmap()
  message("\n-- Cellosaurus REST API --")
  fetch_cellosaurus()
  message("\n-- derived subsets --")
  derive_subsets()
}

message("\n-- verifying against the manifest --")
res <- do.call(rbind, lapply(MANIFEST, verify_one))
print(res, row.names = FALSE, right = FALSE)

nsearch <- length(list.files(CELLO, pattern = "^search_.*\\.json$"))
nrecord <- length(list.files(CELLO, pattern = "^cvcl_.*\\.json$"))
cat(sprintf("\ncellosaurus/: %d search JSON (expected %d) | %d record JSON (expected %d) | names.txt %s\n",
            nsearch, CELLO_EXPECT_SEARCH_JSON, nrecord, CELLO_EXPECT_RECORD_JSON,
            if (file.exists(file.path(CELLO, "names.txt"))) "present" else "MISSING"))

tot <- sum(file.info(list.files(EXT, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE)
cat(sprintf("output/external/ total: %.0f MB across %d files\n", tot / 1e6,
            length(list.files(EXT, recursive = TRUE))))

bad <- res$status != "OK"
if (any(bad)) {
  cat("\n", sum(bad), " file(s) not matching the manifest:\n", sep = "")
  print(res[bad, ], row.names = FALSE, right = FALSE)
  if (!DOWNLOAD)
    cat("Re-run with --download to fetch anything MISSING. A MISMATCH usually means\n",
        "the DepMap release moved: verify the release, then update DEPMAP_FIGSHARE and\n",
        "the MANIFEST entry together (never one without the other).\n")
} else {
  cat("\nAll manifest files verified: sizes, md5 checksums and row/column counts match.\n")
}
stopifnot("output/external/ does not match the manifest (see the table above)" = !any(bad),
          "unexpected number of Cellosaurus search JSON files" =
            nsearch == CELLO_EXPECT_SEARCH_JSON,
          "unexpected number of Cellosaurus record JSON files" =
            nrecord == CELLO_EXPECT_RECORD_JSON)
message("\nfetch_external_data.R complete (", if (DOWNLOAD) "download + verify" else "verify only", ").")
