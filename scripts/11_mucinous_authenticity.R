# =============================================================================
# Script: 11_mucinous_authenticity.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Mucinous-line authentication — the highest-risk histotype, because
#          many "mucinous ovarian" lines are actually metastatic GI carcinomas
#          (Korch 2012; Meagher 2025). For the three OvCAN mucinous lines
#          (TOV2414, VOA8762, VOA8771) we apply the ovarian-vs-GI immunomarker
#          discriminators from RNA + protein, plus WES drivers where available,
#          and give a per-line verdict: genuinely ovarian mucinous vs possible GI.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   4 (authentication) — script 2 of 2
# =============================================================================
# Discriminators (Meagher 2025 10.1002/path.6407; Cheasley 2019 10.1038/
# s41467-019-11862-x; lit review Theme 2):
#   - Ovarian mucinous (MOC):  KRT7+ (diffuse) / SATB2- / PAX8 focal-positive;
#     KRAS-mutant; NO APC or SMAD4 alteration.
#   - Colorectal metastasis:   KRT7- / SATB2+ / CDX2+ / PAX8- ; APC-mutant.
#   - Pancreatic metastasis:   SMAD4 loss.
#   Intestinal-type differentiation (CDX2, MUC2, KRT20, gastric MUC5AC) occurs in
#   BOTH ovarian MOC and GI, so CDX2 ALONE is NOT discriminating — the useful
#   contrast is KRT7/PAX8 (Mullerian) HIGH + SATB2 LOW (ovarian) vs the reverse.
#
# Assay-aware caveats:
#   - These markers are validated as IHC on tissue. We approximate with mRNA
#     (log2 TPM, z across the 31 RNA lines) and TMT protein (log2 abundance, z
#     across 31 protein lines). Cultured lines can down-regulate cytokeratins, so
#     a LOW RNA KRT7 is weaker evidence than a LOW protein KRT7; we show both.
#   - SATB2, CDX2, MUC2, TFF3 are NOT in the proteomics panel -> RNA only.
#   - WES exists for TOV2414 ONLY (tumor-only; canonical hotspots trustworthy).
#     VOA8762/VOA8771 have NO WES and NO external provenance (absent from
#     Cellosaurus; no primary paper) -> our -omics is the only in-house check,
#     which can raise a flag but cannot definitively call origin. STR + IHC from
#     BC Cancer/OVCARE are required and are flagged as pending.
#   - TOV2414 is externally authenticated ovarian mucinous (Sauriol 2020: KRAS
#     G12A, SATB2-, focal PAX8+) -> it serves as our positive anchor.
# =============================================================================

source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({ library(tidyverse); library(matrixStats) })
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
set.seed(SEED)

mc_lines  <- c("TOV2414","VOA8762","VOA8771")
# comparators: serous (PAX8/KRT7 high, ovarian) + a clear-cell line, for context
comparators <- c("OV90","TOV21G","VOA4395")

# --- VERDICT THRESHOLDS, NAMED AND SWEPT [review revision] -------------------
# These three cut-points were bare literals inside the case_when below, with no
# stated basis and no sensitivity check, yet they decide whether two lines are
# described as possibly-GI — a call that determines whether anyone uses them again.
# They are marker-z cut-points on a 31-line reference distribution, i.e. conventions,
# not calibrated decision boundaries. They are named here, written into the output,
# and swept in auth_mucinous_sensitivity.csv so the referee can see which verdicts
# depend on the convention and which do not.
THR_SATB2_POS   <- 0.5    # SATB2 z above this = "SATB2-positive" (colorectal)
THR_CDX2_POS    <- 1.0    # CDX2  z above this = "intestinal CDX2-positive"
THR_OV_INDEX    <- 0.0    # ovarian_index below this = ambiguous/GI-leaning
THR_SATB2_SWEEP <- c(0, 0.25, 0.5, 0.75, 1.0)
THR_CDX2_SWEEP  <- c(0.5, 0.75, 1.0, 1.5)
THR_OVIDX_SWEEP <- c(-0.5, 0, 0.5)

# -----------------------------------------------------------------------------
# 1. Marker matrices (RNA log2 TPM z; protein log2 abundance z)
# -----------------------------------------------------------------------------
tpm  <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g  <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"), show_col_types = FALSE)
g2e  <- t2g %>% distinct(ensembl_gene_id, external_gene_name) %>%
  filter(!is.na(external_gene_name), external_gene_name != "")
m_id <- as.matrix(tpm[, -1]); rownames(m_id) <- tpm$gene_id
id2sym <- setNames(g2e$external_gene_name, g2e$ensembl_gene_id)
common <- intersect(rownames(m_id), names(id2sym))
logtpm <- log2(rowsum(m_id[common, , drop = FALSE], group = id2sym[common]) + 1)
rna_z  <- t(scale(t(logtpm)))

# [integration revision] Shared loader (see 00_setup.R). The local
# `filter(!is.na(protein), !duplicated(protein))` guard was written against the old
# 8,430-row matrix and is handled at source by script 05 since the revision.
# Zero-plex rows are NA in all 31 lines, so t(scale(t())) would return NaN for them;
# they are excluded from the z matrix so no marker lookup can return a z computed
# from no data. None of the 10 mucinous markers is zero-plex (asserted below).
prot <- read_prot_matrix()
prot_zero_plex <- attr(prot, "zero_plex")
prot_z <- t(scale(t(prot[!rownames(prot) %in% prot_zero_plex, , drop = FALSE])))

# --- STATE THE z-REFERENCE SET EXPLICITLY [review revision] ------------------
# The referee could not tell from the figure whether the mucinous z-scores were
# computed across the 3 mucinous columns (which would make one-high/two-low true by
# construction and the panel circular) or across the whole panel. They are computed
# across EVERY line in each assay, independently per assay, and that is now recorded
# in the output and printed here.
Z_REF_RNA  <- sprintf("z per gene across all %d RNA lines (log2 TPM+1)", ncol(rna_z))
Z_REF_PROT <- sprintf("z per gene across all %d protein lines (log2 TMT relative abundance)",
                      ncol(prot_z))
message("z reference sets: ", Z_REF_RNA, " | ", Z_REF_PROT,
        "  (NOT the 3 mucinous columns — the panel is not circular)")

# marker panel with expected direction for OVARIAN identity
markers <- tribble(
  ~symbol,  ~group,               ~ovarian_dir, ~note,
  "KRT7",   "Mullerian (ovarian)", "high", "diffuse in ovarian MOC; often lost in GI",
  "PAX8",   "Mullerian (ovarian)", "high", "Mullerian TF; focal in MOC, negative in GI",
  "WT1",    "Mullerian (ovarian)", "high", "serous>MOC; low overall in MOC",
  "SATB2",  "intestinal (GI)",     "low",  "colorectal marker; NEGATIVE in ovarian MOC",
  "CDX2",   "intestinal (shared)", "n/a",  "intestinal TF; in BOTH MOC and GI (not discriminating alone)",
  "KRT20",  "intestinal (shared)", "n/a",  "intestinal keratin; shared",
  "MUC2",   "intestinal (shared)", "n/a",  "intestinal mucin; shared",
  "MUC5AC", "gastric (mucinous)",  "high", "gastric mucin; mucinous differentiation",
  "TFF1",   "gastric (mucinous)",  "n/a",  "trefoil; secretory/gastric",
  "TFF3",   "gastric (mucinous)",  "n/a",  "trefoil; secretory/intestinal"
)

# A zero-plex protein carries no measurement in any line, so a marker z built from
# one would reflect missingness. None of the 10 markers is zero-plex. [integration revision]
stopifnot("a mucinous marker protein is zero-plex (NA in all 31 lines) — its protein z would be meaningless" =
            !any(markers$symbol %in% prot_zero_plex))

show_lines <- c(mc_lines, comparators)
get_z <- function(mat, g, cl) if (g %in% rownames(mat) && cl %in% colnames(mat)) mat[g, cl] else NA_real_
get_raw <- function(mat, g, cl) if (g %in% rownames(mat) && cl %in% colnames(mat)) mat[g, cl] else NA_real_

marker_long <- tidyr::crossing(symbol = markers$symbol, cell_line = show_lines) %>%
  left_join(markers, by = "symbol") %>%
  rowwise() %>%
  mutate(
    rna_logtpm = round(get_raw(logtpm, symbol, cell_line), 2),
    rna_z      = round(get_z(rna_z, symbol, cell_line), 2),
    prot_log2  = round(get_raw(prot, symbol, cell_line), 2),
    prot_z     = ifelse(symbol %in% rownames(prot), round(get_z(prot_z, symbol, cell_line), 2), NA_real_),
    in_protein = symbol %in% rownames(prot)
  ) %>% ungroup()

# -----------------------------------------------------------------------------
# 2. WES drivers relevant to mucinous origin (TOV2414 only) --------------------
# -----------------------------------------------------------------------------
stopifnot("output/wes_mutations_filtered.csv missing — run 07_wes_mutations.R first" =
            file.exists(file.path(OUT, "wes_mutations_filtered.csv")))
wm <- readr::read_csv(file.path(OUT, "wes_mutations_filtered.csv"), show_col_types = FALSE)
origin_genes <- c("KRAS","APC","SMAD4","TP53","CDKN2A","ERBB2","BRAF","GNAS","RNF43","TGFBR2")
mc_wes <- wm %>% filter(cell_line %in% mc_lines, Hugo_Symbol %in% origin_genes) %>%
  transmute(cell_line, Hugo_Symbol,
            variant = paste0(gsub("^p\\.", "", HGVSp_Short), " (", Variant_Classification, ")"),
            vaf, is_driver) %>%
  arrange(cell_line, Hugo_Symbol)
cat("\n=== Mucinous-origin WES variants (WES = TOV2414 only) ===\n")
print(as.data.frame(mc_wes), row.names = FALSE)
wes_avail <- mc_lines[mc_lines %in% wm$cell_line]

# helper: does a line carry a mutation in gene g?
has_mut <- function(cl, g) cl %in% wm$cell_line &&
  any(wm$cell_line == cl & wm$Hugo_Symbol == g)
kras_change <- function(cl) {
  h <- wm %>% filter(cell_line == cl, Hugo_Symbol == "KRAS")
  if (nrow(h) == 0) NA_character_ else paste(gsub("^p\\.", "", h$HGVSp_Short), collapse = ";")
}

# -----------------------------------------------------------------------------
# 3. Per-line ovarian-vs-GI scoring + verdict ---------------------------------
# -----------------------------------------------------------------------------
# Combine the discriminating markers (KRT7, PAX8 -> ovarian; SATB2 -> GI) into a
# transparent index. We DO NOT use CDX2/KRT20/MUC2 in the index (shared), but
# report them.
# ASSAY-AWARE marker readout (important): the nuclear transcription factors
# PAX8/WT1/CDX2/SATB2/HNF1B are LOW-ABUNDANCE and are UNDER-QUANTIFIED by TMT
# (the Phase-2 QC flagged PAX8/WT1/HNF1B as near-flat across channels), so their
# PROTEIN z is an artifact and RNA is authoritative for them. Abundant
# structural/secretory markers (KRT7, MUC5AC, KRT20) are well measured by TMT, so
# for those we prefer protein z (culture-robust) and fall back to RNA. Using
# protein z for PAX8 would spuriously soften a GI-leaning line (e.g. VOA8771,
# whose flat TMT PAX8 reads ~0 while its mRNA is clearly low).
tf_markers <- c("PAX8","WT1","CDX2","SATB2","HNF1B","NAPSA")     # RNA authoritative
z_best <- function(cl, g) {
  if (g %in% tf_markers) return(get_z(rna_z, g, cl))            # TF/nuclear -> RNA
  p <- get_z(prot_z, g, cl)                                     # structural/secretory -> protein if present
  if (!is.na(p) && g %in% rownames(prot)) p else get_z(rna_z, g, cl)
}
verdict_tbl <- tibble(cell_line = mc_lines) %>%
  rowwise() %>%
  mutate(
    KRT7_z   = round(z_best(cell_line, "KRT7"), 2),      # abundant cytokeratin -> protein-preferred
    PAX8_z   = round(z_best(cell_line, "PAX8"), 2),      # TF -> RNA (TMT under-quantifies)
    WT1_z    = round(z_best(cell_line, "WT1"), 2),       # TF -> RNA
    SATB2_rna_z = round(get_z(rna_z, "SATB2", cell_line), 2),
    CDX2_rna_z  = round(get_z(rna_z, "CDX2",  cell_line), 2),
    MUC5AC_z = round(z_best(cell_line, "MUC5AC"), 2),    # secretory mucin -> protein-preferred
    # --- ovarian_index: AN AD-HOC, UNWEIGHTED COMPOSITE. Document it as such. ---
    # ovarian_index = mean(KRT7_z, PAX8_z) - SATB2_rna_z
    #   i.e. (Mullerian evidence) minus (colorectal evidence), on a z scale.
    # [review revision] Three properties a reader must know before quoting it:
    #   (1) It MIXES ASSAYS. KRT7 and MUC5AC take PROTEIN z (abundant, well measured
    #       by TMT); PAX8/WT1/CDX2/SATB2 take RNA z (low-abundance TFs that TMT
    #       under-quantifies). So the index adds a protein z to an RNA z. That choice
    #       is assay-justified (see the z_best note above) but it is not a single
    #       measurement scale, and the two assays have different spreads.
    #   (2) It is UNWEIGHTED and UNCALIBRATED — no training set, no external anchor.
    #       KRT7 and PAX8 contribute half each; SATB2 contributes a full unit. There
    #       is no evidence that 1 z of SATB2 offsets 1 z of Mullerian marker.
    #   (3) Its verdict thresholds are conventions (see THR_* above), swept below.
    # It is a transparent summary of the three discriminating markers, NOT a
    # classifier. The per-marker values and their panel-wide ranks (written to
    # auth_mucinous_marker_ranks.csv) are the primary evidence.
    ovarian_index = round(mean(c(KRT7_z, PAX8_z), na.rm = TRUE) - SATB2_rna_z, 2),
    ovarian_index_definition =
      "mean(KRT7_z, PAX8_z) - SATB2_rna_z ; ad-hoc unweighted composite mixing protein z (KRT7) with RNA z (PAX8, SATB2); uncalibrated",
    KRAS   = kras_change(cell_line),
    APC    = if (cell_line %in% wm$cell_line) (if (has_mut(cell_line, "APC")) "candidate present" else "no retained candidate") else "no WES",
    SMAD4  = if (cell_line %in% wm$cell_line) (if (has_mut(cell_line, "SMAD4")) "candidate present" else "no retained candidate") else "no WES",
    # EXTERNAL (literature/registry) status, kept STRICTLY SEPARATE from anything
    # measured here and explicitly attributed. [review revision] This string
    # previously read "...KRAS G12A, SATB2-", and the manuscript reproduced
    # "KRT7+/PAX8+/MUC5AC+/SATB2-" as though all four were in-house measurements.
    # Three of them are; SATB2-negative is Sauriol's IHC, and in THESE data
    # TOV2414's SATB2 is not low (see literature_vs_measured below).
    external_provenance = ifelse(cell_line == "TOV2414",
      "Cellosaurus CVCL_A1SR; primary paper Sauriol 2020 REPORTS (their IHC/sequencing, not measured here): KRAS G12A, SATB2-negative, focal PAX8+",
      "NONE (absent from Cellosaurus; no primary paper)")
  ) %>% ungroup()

# --- LITERATURE vs MEASURED, per line [review revision] ----------------------
# The one place where the two disagree, stated in the output so the text cannot
# silently substitute the literature value for the measurement.
satb2_rank <- rank(rna_z["SATB2", ], ties.method = "min")   # 1 = lowest of the panel
verdict_tbl <- verdict_tbl %>%
  rowwise() %>%
  mutate(SATB2_rna_rank = as.integer(satb2_rank[cell_line]),
         SATB2_rna_rank_of = length(satb2_rank),
         literature_vs_measured = if (cell_line == "TOV2414")
           sprintf(paste0("DISCORDANT: Sauriol 2020 reports SATB2-NEGATIVE by IHC, but measured SATB2 mRNA z here is %+.2f ",
                          "(rank %d of %d, i.e. mid-panel) — HIGHER than VOA8762 (%+.2f) and VOA8771 (%+.2f). ",
                          "The ovarian call for TOV2414 rests on measured KRT7/PAX8/MUC5AC + KRAS G12A + published derivation/reference STR record, ",
                          "NOT on a measured SATB2-low. Do not quote SATB2- as an in-house measurement. ",
                          "Possible explanations: IHC protein vs mRNA discordance, or 2D-culture drift."),
                   SATB2_rna_z, as.integer(satb2_rank[cell_line]), length(satb2_rank),
                   round(rna_z["SATB2", "VOA8762"], 2), round(rna_z["SATB2", "VOA8771"], 2))
         else "no external literature value to compare (no primary paper)") %>%
  ungroup()

# --- verdict logic: MEASURED EXPRESSION FIRST, provenance stated separately ---
# [review revision] TOV2414 previously had a HARD-CODED verdict branch that bypassed
# its own measured markers and asserted the Sauriol phenotype. Now every line goes
# through the SAME data-derived rule (expression_verdict), and the external status is
# composed on afterwards. TOV2414 still comes out ovarian-compatible — but on its
# measured KRT7/PAX8/MUC5AC, which is the point.
expr_verdict <- function(SATB2_rna_z, CDX2_rna_z, PAX8_z, KRT7_z, ovarian_index,
                         t_satb2 = THR_SATB2_POS, t_cdx2 = THR_CDX2_POS,
                         t_ovidx = THR_OV_INDEX) {
  if (!is.na(SATB2_rna_z) && !is.na(CDX2_rna_z) &&
      SATB2_rna_z > t_satb2 && CDX2_rna_z > t_cdx2)
    "GI/COLORECTAL-LEANING (SATB2+/CDX2+)"
  else if (!is.na(PAX8_z) && !is.na(KRT7_z) && !is.na(CDX2_rna_z) &&
           PAX8_z < 0 && KRT7_z < 0 && CDX2_rna_z > t_cdx2)
    "GI-LEANING (Mullerian-low PAX8-/KRT7-, intestinal CDX2+)"
  else if (!is.na(ovarian_index) && ovarian_index < t_ovidx)
    "AMBIGUOUS/GI-leaning (ovarian index below threshold)"
  else "OVARIAN-COMPATIBLE expression (Mullerian-positive, not SATB2/CDX2-driven)"
}
verdict_tbl <- verdict_tbl %>%
  rowwise() %>%
  mutate(
    expression_verdict = expr_verdict(SATB2_rna_z, CDX2_rna_z, PAX8_z, KRT7_z, ovarian_index),
    verdict = paste0(
      expression_verdict,
      sprintf(" [measured: KRT7 z=%+.2f, PAX8 z=%+.2f, SATB2 z=%+.2f, CDX2 z=%+.2f, MUC5AC z=%+.2f, ovarian_index=%+.2f]",
              KRT7_z, PAX8_z, SATB2_rna_z, CDX2_rna_z, MUC5AC_z, ovarian_index),
      if (cell_line == "TOV2414")
        paste0(". Published ovarian derivation (Cellosaurus CVCL_A1SR; Sauriol 2020) and a matching KRAS G12A candidate observed here, ",
               "so the ovarian call does not rest on expression alone. Carries a tumour-only candidate SMAD4 frameshift + low SMAD4 mRNA ",
               "— a co-occurring event (SMAD4 loss occurs in a subset of MOC), NOT evidence of pancreatic origin given the ",
               "published derivation/reference STR record. NB its measured SATB2 is NOT low (see literature_vs_measured)")
      else
        paste0(". NO external provenance (absent from Cellosaurus, no primary paper), so expression is the only in-house ",
               "check and it can raise a flag but cannot call origin. STR + IHC required"),
      "."),
    needs = ifelse(cell_line == "TOV2414", "Reference STR profile available; obtain current-stock STR/mycoplasma documentation",
                   "STR profile + histotype IHC (CK7/SATB2/PAX8/WT1) + WES (KRAS/APC/SMAD4) from BC Cancer/OVCARE")
  ) %>% ungroup()

# --- THRESHOLD SENSITIVITY [review revision] ---------------------------------
# Which verdicts are conclusions and which are consequences of a convention?
sens_grid <- tidyr::expand_grid(t_satb2 = THR_SATB2_SWEEP, t_cdx2 = THR_CDX2_SWEEP,
                                t_ovidx = THR_OVIDX_SWEEP,
                                cell_line = mc_lines) %>%
  left_join(verdict_tbl %>% select(cell_line, SATB2_rna_z, CDX2_rna_z, PAX8_z,
                                   KRT7_z, ovarian_index), by = "cell_line") %>%
  rowwise() %>%
  mutate(verdict_class = expr_verdict(SATB2_rna_z, CDX2_rna_z, PAX8_z, KRT7_z,
                                      ovarian_index, t_satb2, t_cdx2, t_ovidx),
         is_headline = t_satb2 == THR_SATB2_POS & t_cdx2 == THR_CDX2_POS &
                       t_ovidx == THR_OV_INDEX) %>%
  ungroup()
readr::write_csv(sens_grid, file.path(OUT, "auth_mucinous_sensitivity.csv"))
cat("\n=== Verdict-threshold SENSITIVITY: how many of the ", nrow(sens_grid) / length(mc_lines),
    " threshold combinations give each verdict? ===\n", sep = "")
print(as.data.frame(sens_grid %>% dplyr::count(cell_line, verdict_class) %>%
  group_by(cell_line) %>%
  mutate(pct_of_grid = round(100 * n / sum(n))) %>% ungroup()), row.names = FALSE)
stab <- sens_grid %>% group_by(cell_line) %>%
  summarise(n_distinct_verdicts = dplyr::n_distinct(verdict_class),
            headline = verdict_class[is_headline][1], .groups = "drop")
print(as.data.frame(stab), row.names = FALSE)
cat("  => a line with n_distinct_verdicts == 1 is a conclusion; > 1 means the verdict\n",
    "     depends on the (uncalibrated) threshold convention and must be stated as such.\n")

auth_mucinous <- verdict_tbl %>%
  mutate(z_reference_rna = Z_REF_RNA, z_reference_protein = Z_REF_PROT,
         thresholds = sprintf("SATB2>%.2f, CDX2>%.2f, ovarian_index<%.2f (conventions; see auth_mucinous_sensitivity.csv)",
                              THR_SATB2_POS, THR_CDX2_POS, THR_OV_INDEX)) %>%
  select(cell_line, KRT7_z, PAX8_z, WT1_z, SATB2_rna_z, SATB2_rna_rank,
         SATB2_rna_rank_of, CDX2_rna_z, MUC5AC_z,
         ovarian_index, ovarian_index_definition, KRAS, APC, SMAD4,
         expression_verdict, external_provenance, literature_vs_measured,
         z_reference_rna, z_reference_protein, thresholds, verdict, needs)
readr::write_csv(auth_mucinous, file.path(OUT, "auth_mucinous.csv"))

# --- PER-MARKER VALUES AGAINST THE FULL PANEL DISTRIBUTION [review revision] ---
# The referee's concern was that a 3-column z-score would be circular by
# construction. It is not (see Z_REF_* above), but a 3-column FIGURE cannot show
# that. This table gives, for every marker, each mucinous line's value, its z, its
# RANK among all lines in that assay, and the panel min/median/max — so the reader
# can see where the three lines sit in the whole distribution.
panel_stats <- function(mat, matz, assay, unit) {
  zref <- if (identical(assay, "RNA")) Z_REF_RNA else Z_REF_PROT
  as_tibble(mat[intersect(markers$symbol, rownames(mat)), , drop = FALSE], rownames = "symbol") %>%
    pivot_longer(-symbol, names_to = "cell_line", values_to = "value") %>%
    left_join(as_tibble(matz[intersect(markers$symbol, rownames(matz)), , drop = FALSE],
                        rownames = "symbol") %>%
                pivot_longer(-symbol, names_to = "cell_line", values_to = "z"),
              by = c("symbol","cell_line")) %>%
    group_by(symbol) %>%
    mutate(rank_low_is_lowest = rank(value, ties.method = "min"),
           n_lines_in_assay = dplyr::n(),
           panel_min = min(value), panel_median = median(value), panel_max = max(value)) %>%
    ungroup() %>%
    mutate(assay = assay, unit = unit, z_reference_set = zref)
}
marker_ranks <- bind_rows(panel_stats(logtpm, rna_z,  "RNA",     "log2(TPM+1)"),
                          panel_stats(prot,   prot_z, "protein", "log2 TMT relative abundance")) %>%
  left_join(markers, by = "symbol") %>%
  mutate(is_mucinous_line = cell_line %in% mc_lines,
         percentile = round(100 * (rank_low_is_lowest - 0.5) / n_lines_in_assay, 1)) %>%
  select(symbol, group, ovarian_dir, cell_line, is_mucinous_line, assay, unit, value, z,
         rank_low_is_lowest, n_lines_in_assay, percentile,
         panel_min, panel_median, panel_max, z_reference_set, note) %>%
  arrange(symbol, assay, desc(z))
readr::write_csv(marker_ranks, file.path(OUT, "auth_mucinous_marker_ranks.csv"))

cat("\n=== The 3 mucinous lines against the FULL panel distribution (not 3 columns) ===\n")
print(as.data.frame(marker_ranks %>% filter(is_mucinous_line) %>%
  transmute(symbol, assay, cell_line, z = round(z, 2),
            rank = sprintf("%d/%d", rank_low_is_lowest, n_lines_in_assay),
            pctile = percentile) %>%
  pivot_wider(names_from = cell_line, values_from = c(z, rank, pctile),
              names_glue = "{cell_line}_{.value}") %>%
  arrange(assay, symbol)), row.names = FALSE)
message("Wrote output/auth_mucinous_marker_ranks.csv (", nrow(marker_ranks),
        " marker x line x assay rows) and output/auth_mucinous_sensitivity.csv")

cat("\n=== Mucinous authentication verdicts ===\n")
print(as.data.frame(auth_mucinous %>%
  select(cell_line, KRT7_z, PAX8_z, SATB2_rna_z, CDX2_rna_z, MUC5AC_z,
         ovarian_index, KRAS, APC, SMAD4)), row.names = FALSE)
for (i in seq_len(nrow(auth_mucinous)))
  cat(sprintf("\n%s: %s\n   provenance: %s\n", auth_mucinous$cell_line[i],
              auth_mucinous$verdict[i], auth_mucinous$external_provenance[i]))

# -----------------------------------------------------------------------------
# 4. Figure — marker panel (RNA + protein z) for the 3 MC lines + comparators
# -----------------------------------------------------------------------------
plot_long <- marker_long %>%
  select(symbol, group, cell_line, ovarian_dir, rna_z, prot_z) %>%
  pivot_longer(c(rna_z, prot_z), names_to = "assay", values_to = "z") %>%
  mutate(assay = recode(assay, rna_z = "RNA (z log2 TPM)", prot_z = "Protein (z log2 abund.)"),
         assay = factor(assay, levels = c("RNA (z log2 TPM)","Protein (z log2 abund.)")),
         symbol = factor(symbol, levels = rev(markers$symbol)),
         cell_line = factor(cell_line, levels = show_lines),
         is_mc = cell_line %in% mc_lines,
         z = pmax(pmin(z, 2), -2))

# annotate the ovarian-relevant direction on the y-axis via group facet
ann_rows <- markers %>% mutate(symbol = factor(symbol, levels = rev(markers$symbol)))

xcol  <- ifelse(show_lines %in% mc_lines, COOK_RUST, "grey30")
xface <- ifelse(show_lines %in% mc_lines, "bold", "plain")

p_mk <- ggplot(plot_long, aes(cell_line, symbol, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = subset(plot_long, is.na(z)),
            aes(label = "NA"), size = 2.3, colour = "grey60") +
  scale_fill_gradient2(low = COOK_NAVY, mid = "white", high = COOK_RUST,
                       midpoint = 0, limits = c(-2, 2), name = "z", na.value = "grey92") +
  facet_grid(group ~ assay, scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Mucinous-line authentication: ovarian-vs-GI markers",
       subtitle = paste0("Ovarian MOC = KRT7+/PAX8+/SATB2- ; colorectal = SATB2+/CDX2+/KRT7-/PAX8-.  ",
                        "Rust = the 3 MC lines.\nCDX2/KRT20/MUC2 shared (not discriminating alone).  ",
                        "z across 31 RNA / 31 protein lines; SATB2/CDX2/MUC2/TFF3 not in TMT panel.\n",
                        "NB: PAX8/WT1 are TFs under-quantified by TMT (flat) -> RNA is authoritative for the verdict index."),
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = xcol, face = xface),
        strip.text.y.left = element_text(angle = 0, size = 7.5),
        strip.placement = "outside", panel.grid = element_blank(),
        legend.position = "right")
ggsave(file.path(FIGS, "11_mucinous_markers.pdf"), p_mk, width = 8.5, height = 5.5)

# --- companion: KRT7/PAX8 (ovarian) vs CDX2/SATB2 (GI) scatter, MC lines labeled
scat_all <- tibble(
    cell_line = colnames(rna_z),
    mullerian     = sapply(colnames(rna_z), function(cl) mean(c(z_best(cl, "KRT7"), z_best(cl, "PAX8")))),
    intestinal_GI = sapply(colnames(rna_z), function(cl) mean(c(get_z(rna_z, "CDX2", cl),
                                                                get_z(rna_z, "SATB2", cl))))) %>%
  mutate(grp = dplyr::case_when(cell_line == "TOV2414" ~ "TOV2414 (auth. ovarian)",
                                cell_line %in% c("VOA8762","VOA8771") ~ "VOA MC (unverified)",
                                TRUE ~ "other lines"))
p_sc <- ggplot(scat_all, aes(mullerian, intestinal_GI)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  geom_point(aes(colour = grp, size = grp)) +
  ggrepel::geom_text_repel(data = subset(scat_all, cell_line %in% mc_lines),
                           aes(label = cell_line), size = 3, min.segment.length = 0) +
  scale_colour_manual(values = c("TOV2414 (auth. ovarian)" = "#0072B2",
                                 "VOA MC (unverified)" = COOK_RUST,
                                 "other lines" = "grey75"), name = NULL) +
  scale_size_manual(values = c("TOV2414 (auth. ovarian)" = 3,
                               "VOA MC (unverified)" = 3, "other lines" = 1.6), guide = "none") +
  labs(title = "Mullerian (KRT7/PAX8) vs intestinal-GI (CDX2/SATB2) balance",
       subtitle = "ovarian mucinous expected upper-left-ish (Mullerian+); GI-leaning lower-right (intestinal+, Mullerian-)",
       x = "Mullerian score (mean z: KRT7, PAX8)",
       y = "Intestinal/GI score (mean z: CDX2, SATB2)") +
  theme_bw(base_size = 9) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(FIGS, "11_mucinous_scatter.pdf"), p_sc, width = 6.5, height = 5.5)

cat("\nOutputs: output/auth_mucinous.csv\n")
cat("Figures: figs/11_mucinous_markers.pdf, figs/11_mucinous_scatter.pdf\n")
message("\n11_mucinous_authenticity.R complete.")
