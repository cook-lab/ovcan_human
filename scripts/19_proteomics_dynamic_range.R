# =============================================================================
# Script: 19_proteomics_dynamic_range.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Peer-review response (§3.7) — quantify the assay-specific limits of
#          the TMT proteomics that shape how the resource should be reused:
#            (1) CROSS-ASSAY SPREAD — cross-line dynamic range at protein vs
#                RNA (range / IQR / SD per gene, and the protein/RNA ratio).
#                The spread ratio is descriptive. It does not identify a causal
#                compression factor, establish a correlation ceiling, or measure
#                model-selection discrimination.
#            (2) BRIDGE DESIGN — the inter-plex bridge is a daisy-CHAIN
#                (plex1<->2<->3<->4<->5) not a common-reference HUB; describe the
#                inter-plex normalization ACTUALLY implemented (PIS-based).
#            (3) STRUCTURAL BLOCK-MISSINGNESS — TMT missingness is per-plex, so
#                some proteins are absent from ENTIRE plex-blocks (all lines in
#                that plex), not missing at random.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   Revision (peer-review response) — proteomics depth/limits
# -----------------------------------------------------------------------------
# READ-ONLY inputs & provenance:
#   output/prot_abundance_matrix.csv  — 8,427 proteins x 31 lines, log2 relative
#         abundance (NAs preserved). Values are supplied log2 abundances with feature-specific baselines.
#         PIS normalization is reported; the exact upstream scaling formula needs confirmation.
#   output/rna_tpm.csv + output/tx2gene_matched.csv — RNA, collapsed to
#         gene symbol (sum TPM) then log2(TPM+1), IDENTICAL to script 12 so the
#         dynamic-range result speaks directly to the concordance number.
#   output/integ_rnaprot_cor.csv — per-gene RNA-protein Spearman (from script 12)
#         used to describe association of spread with concordance.
#   output/prot_bridge_cor.csv — per-link inter-plex bridge correlations (05).
#   judy_archive/data/proteomics/{protein_relative_abundance,tmt.layout}.xlsx —
#         re-read ONLY to (a) confirm the PIS/bridge layout and (b) estimate the
#         bridge-derived precision estimates from repeat-sample differences.
#   metadata/samples.csv — plex/channel/subtype/site per line.
# Scripts 05/06 are READ-ONLY references for the normalization description.
# Writes: output/prot_dynamic_range.csv, output/prot_block_missingness.csv,
#         figs/f_prot_compression.pdf, reports/assets/f_prot_compression.png.
# =============================================================================

source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(matrixStats); library(patchwork)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
set.seed(SEED)

SUBTYPE_LEVELS <- c("HGS", "CC", "EC", "MC", "MMMT", "SCCOHT")
MIN_LINES <- 10L   # min paired lines for a stable per-gene spread estimate (matches script 12)

dir.create(file.path(PROJ, "reports", "assets"), showWarnings = FALSE, recursive = TRUE)

# Lab theme — identical to script 12 so this figure matches the concordance figs
theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text = element_text(colour = "black"),
      plot.title = element_text(size = rel(1.1), hjust = 0, face = "bold",
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = rel(0.85), hjust = 0, colour = "grey30",
                                   margin = margin(b = 6)),
      axis.title = element_text(size = rel(0.95)),
      axis.text = element_text(size = rel(0.8), colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      legend.title = element_text(size = rel(0.85)),
      legend.text = element_text(size = rel(0.8)),
      legend.background = element_blank(), legend.key = element_blank(),
      panel.background = element_blank(), panel.border = element_blank(),
      panel.grid = element_blank(), strip.background = element_blank(),
      strip.text = element_text(size = rel(0.9), face = "bold"),
      plot.margin = margin(8, 8, 8, 8))
}
COL_RNA  <- COOK_NAVY   # RNA in navy, protein in rust (no red-green, on-brand)
COL_PROT <- COOK_RUST

# =============================================================================
# 1. Build the SHARED gene x line RNA / protein matrices (mirror script 12)
# =============================================================================
tpm  <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g  <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"), show_col_types = FALSE)
gene_meta <- t2g %>% distinct(ensembl_gene_id, external_gene_name) %>%
  filter(!is.na(external_gene_name), external_gene_name != "")

m_id    <- as.matrix(tpm[, -1]); rownames(m_id) <- tpm$gene_id
id2sym  <- setNames(gene_meta$external_gene_name, gene_meta$ensembl_gene_id)
common  <- intersect(rownames(m_id), names(id2sym))
tpm_sym <- rowsum(m_id[common, , drop = FALSE], group = id2sym[common])
rna_log <- log2(tpm_sym + 1)                                   # symbols x RNA lines

# [integration revision] Shared loader (00_setup.R); feature hygiene handled at source
# by script 05. The 70 zero-plex rows are retained in the matrix (the missingness
# accounting in section 4 below has to see them) but are excluded explicitly from the
# abundance-decile stratification, where rowMeans() over an all-NA row returns NaN.
prot_mat <- read_prot_matrix()
prot_zero_plex <- attr(prot_mat, "zero_plex")

shared_lines <- intersect(colnames(rna_log), colnames(prot_mat))
shared_genes <- intersect(rownames(rna_log), rownames(prot_mat))
R <- rna_log[shared_genes, shared_lines]                        # RNA log2(TPM+1)
P <- prot_mat[shared_genes, shared_lines]                       # protein log2 rel. abundance
stopifnot(identical(dim(R), dim(P)), identical(dimnames(R), dimnames(P)))
message(sprintf("Paired set: %d shared genes x %d shared lines (RNA-only: %s; protein-only: %s)",
                length(shared_genes), length(shared_lines),
                paste(setdiff(colnames(rna_log), shared_lines), collapse = ","),
                paste(setdiff(colnames(prot_mat), shared_lines), collapse = ",")))

# =============================================================================
# 2. TASK 1 — cross-line dynamic range at protein vs RNA (per gene)
#    n-MATCHED per gene: for each gene, compute spread over the SAME lines (the
#    lines where protein is present; RNA is always finite). Range is sensitive
#    to n, so matching the line set removes that confound. Require >= MIN_LINES.
# =============================================================================
spread_stats <- function(x) {                                    # x = numeric vector
  x <- x[is.finite(x)]
  c(min = min(x), max = max(x), range = max(x) - min(x),
    iqr = as.numeric(diff(quantile(x, c(.25, .75)))), sd = sd(x))
}

dr <- lapply(shared_genes, function(g) {
  xp <- P[g, ]; xr <- R[g, ]
  ok <- is.finite(xp) & is.finite(xr)                            # == protein-present lines
  n  <- sum(ok)
  if (n < MIN_LINES) return(NULL)
  rs <- spread_stats(xr[ok]); ps <- spread_stats(xp[ok])
  tibble(gene = g, n_paired = n,
         rna_min = rs["min"],  rna_max = rs["max"],  rna_range = rs["range"],
         rna_iqr = rs["iqr"],  rna_sd  = rs["sd"],
         prot_min = ps["min"], prot_max = ps["max"], prot_range = ps["range"],
         prot_iqr = ps["iqr"], prot_sd  = ps["sd"])
}) %>% bind_rows()

# ratios (protein / RNA); guard zero-variance RNA denominators
# SPREAD RATIOS compare observed scales; they do not estimate causal compression.
# The two scales have different FLOORS: RNA is log2(TPM+1) with a hard floor at 0
# for an undetected transcript. The supplied protein scale retains feature
# baselines; its exact normalization formula needs confirmation, and an observed
# minimum is not a calibrated limit of detection. So a RANGE ratio partly measures the RNA zero-floor,
# not solely assay compression. IQR reduces sensitivity to isolated extremes,
# whereas SD and IQR can both depend on the transformation and abundance, so:
#     PRIMARY   = iqr_ratio and sd_ratio, descriptive and transformation-dependent
#     SECONDARY = range_ratio, reported but LABELLED floor-sensitive
# rna_has_zero flags the genes where the floor is actually engaged, so the
# floor-sensitivity of the range ratio can be shown rather than argued.
dr <- dr %>%
  mutate(range_ratio_floor_sensitive = ifelse(rna_range > 0, prot_range / rna_range, NA_real_),
         range_ratio = range_ratio_floor_sensitive,   # kept: downstream/figure name
         iqr_ratio   = ifelse(rna_iqr   > 0, prot_iqr   / rna_iqr,   NA_real_),
         sd_ratio    = ifelse(rna_sd    > 0, prot_sd    / rna_sd,    NA_real_),
         complete30  = n_paired == length(shared_lines))
# does the RNA zero-floor touch this gene among its paired lines?
rna_zero_n <- vapply(dr$gene, function(g) {
  ok <- is.finite(P[g, ]) & is.finite(R[g, ]); sum(R[g, ok] == 0)
}, numeric(1))
dr <- dr %>% mutate(rna_n_zero_lines = as.integer(rna_zero_n),
                    rna_has_zero     = rna_n_zero_lines > 0)

# attach per-gene concordance (script 12) as a descriptive association
pg_cor <- readr::read_csv(file.path(OUT, "integ_rnaprot_cor.csv"), show_col_types = FALSE) %>%
  filter(level == "per_gene") %>% transmute(gene = id, concord_spearman = spearman)
dr <- dr %>% left_join(pg_cor, by = "gene")

# ADC targets (the reuse shortlist layer) — flag for consequence (b)
adc_genes <- readr::read_csv(file.path(OUT, "adc_expression.csv"),
                             show_col_types = FALSE) %>% pull(symbol) %>% unique()
dr <- dr %>% mutate(is_adc_target = gene %in% adc_genes)

readr::write_csv(dr, file.path(OUT, "prot_dynamic_range.csv"))

# --- headline metrics --------------------------------------------------------
cc30 <- dr %>% filter(complete30)                                # fully n-matched (30 lines each)
med  <- function(x) median(x, na.rm = TRUE)
frac_compressed <- mean(dr$iqr_ratio < 1, na.rm = TRUE)

# FOLR1 exemplar
folr1 <- dr %>% filter(gene == "FOLR1")

# Descriptive association between observed protein spread and concordance.
# These gene-level quantities cannot identify an assay-compression mechanism.
ceil_test <- dr %>% filter(is.finite(concord_spearman))
rho_prot_iqr <- cor(ceil_test$prot_iqr, ceil_test$concord_spearman, method = "spearman")
rho_rna_iqr  <- cor(ceil_test$rna_iqr,  ceil_test$concord_spearman, method = "spearman")
# concordance in low vs high protein-spread terciles
ceil_binned <- ceil_test %>%
  mutate(prot_iqr_tercile = cut(prot_iqr, quantile(prot_iqr, c(0, 1/3, 2/3, 1)),
                                labels = c("low", "mid", "high"), include.lowest = TRUE)) %>%
  group_by(prot_iqr_tercile) %>%
  summarise(n = n(), median_concord = median(concord_spearman), .groups = "drop")

cat("\n=== TASK 1: cross-line dynamic range (protein vs RNA) ===\n")
cat(sprintf("Genes with >= %d paired lines: %d (of %d shared); fully n-matched at 30 lines: %d\n",
            MIN_LINES, nrow(dr), length(shared_genes), nrow(cc30)))
cat("\nMedian per-gene cross-line spread (log2 units) [n-matched, complete-30 subset]:\n")
cat(sprintf("  IQR   : RNA %.2f  vs  protein %.2f   (protein/RNA median ratio = %.2f)  <- PRIMARY (less sensitive to isolated extremes)\n",
            med(cc30$rna_iqr), med(cc30$prot_iqr), med(cc30$iqr_ratio)))
cat(sprintf("  SD    : RNA %.2f  vs  protein %.2f   (protein/RNA median ratio = %.2f)  <- COMPLEMENTARY SPREAD SUMMARY\n",
            med(cc30$rna_sd), med(cc30$prot_sd), med(cc30$sd_ratio)))
cat(sprintf("  RANGE : RNA %.2f  vs  protein %.2f   (protein/RNA median ratio = %.2f)  <- SECONDARY, FLOOR-SENSITIVE\n",
            med(cc30$rna_range), med(cc30$prot_range), med(cc30$range_ratio)))
cat(sprintf("\nFull >=%d-line set: median IQR ratio %.2f; %.1f%% of genes have protein IQR < RNA IQR\n",
            MIN_LINES, med(dr$iqr_ratio), 100 * frac_compressed))

# Descriptive stratification by observed zero RNA values. These groups contain
# different genes, so their differences cannot isolate a transformation effect.
# All spread ratios remain abundance-, feature-, and transformation-dependent.
floor_split <- dr %>% group_by(rna_has_zero) %>%
  summarise(n_genes = dplyr::n(),
            median_range_ratio = med(range_ratio),
            median_iqr_ratio   = med(iqr_ratio),
            median_sd_ratio    = med(sd_ratio),
            median_rna_range   = med(rna_range),
            median_rna_iqr     = med(rna_iqr), .groups = "drop")
cat("\nFLOOR CHECK - genes split by whether any paired line has RNA == 0 (log2(TPM+1) floor):\n")
print(as.data.frame(floor_split %>% mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)
cat(sprintf("  %.1f%% of genes have at least one zero-RNA model.\n",
            100 * mean(dr$rna_has_zero)))
cat("  The two strata contain different genes; their spread differences do not isolate\n")
cat("  a causal effect of the RNA floor. IQR, SD and range can depend on abundance,\n")
cat("  feature selection and log(TPM+1). These are descriptive cross-assay spread ratios.\n")
readr::write_csv(floor_split, file.path(OUT, "prot_compression_floor_check.csv"))

# --- per-gene exemplars: both range and robust-spread summaries
exemplars <- dr %>% filter(gene %in% c("FOLR1", "MSLN", "ERBB2", "TACSTD2")) %>%
  select(gene, n_paired, rna_n_zero_lines,
         rna_range, prot_range, range_ratio,
         rna_iqr, prot_iqr, iqr_ratio, rna_sd, prot_sd, sd_ratio)
cat("\nExemplars - the MSLN 'protein 1.8 vs RNA 9.3' comparison is the floor-sensitive one:\n")
print(as.data.frame(exemplars %>% mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)
cat(sprintf("\nDescriptive association — per-gene spread vs RNA-protein concordance (n=%d genes):\n", nrow(ceil_test)))
cat(sprintf("  Spearman(protein IQR, concordance) = %.3f ; Spearman(RNA IQR, concordance) = %.3f\n",
            rho_prot_iqr, rho_rna_iqr))
cat("  Concordance by protein-spread tercile:\n"); print(as.data.frame(ceil_binned), row.names = FALSE)

# ADC-target dynamic range (consequence b)
adc_dr <- dr %>% filter(is_adc_target) %>%
  select(gene, n_paired, rna_range, prot_range, range_ratio, concord_spearman) %>%
  arrange(range_ratio)
cat("\nConsequence (b) — ADC-target cross-line dynamic range (protein vs RNA):\n")
print(as.data.frame(adc_dr %>% mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)

# =============================================================================
# 3. TASK 3 — structural block-missingness (per-plex)
#    Recompute plex coverage from the NA-preserving matrix + plex map, so the
#    block-missingness table is self-contained and internally consistent.
# =============================================================================
ss <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  filter(proteomics == "Y", provenance == "generated", analysis_include == "Y")
plex_of <- setNames(as.integer(ss$tmt_plex), ss$cell_line)
plexes  <- sort(unique(plex_of))
stopifnot("all 31 analysis lines present in matrix" =
            all(names(plex_of) %in% colnames(prot_mat)))

full_mat <- prot_mat[, names(plex_of), drop = FALSE]             # 31 analysis lines
present_by_plex <- sapply(plexes, function(p) {
  cols <- names(plex_of)[plex_of == p]
  rowSums(!is.na(full_mat[, cols, drop = FALSE])) > 0            # protein seen in >=1 line of plex p
})
colnames(present_by_plex) <- paste0("plex", plexes)
present_n_plex  <- rowSums(present_by_plex)
present_n_lines <- rowSums(!is.na(full_mat))
n_lines_total   <- ncol(full_mat)

plexes_present <- apply(present_by_plex, 1, function(v) paste(plexes[v], collapse = ";"))
plexes_absent  <- apply(present_by_plex, 1, function(v) paste(plexes[!v], collapse = ";"))

block <- tibble(
  protein         = rownames(full_mat),
  present_n_lines = present_n_lines,
  present_n_plex  = present_n_plex,
  plexes_present  = plexes_present,
  plexes_absent   = plexes_absent,
  absent_from_any_plex = present_n_plex < length(plexes),
  pct_lines_missing    = round(100 * (1 - present_n_lines / n_lines_total), 1),
  pass_presence50      = present_n_lines >= ceiling(0.5 * n_lines_total))
readr::write_csv(block, file.path(OUT, "prot_block_missingness.csv"))

# per-plex summary: proteins quantified + lines + subtype composition
plex_summary <- tibble(
  plex        = plexes,
  n_lines     = as.integer(table(plex_of)[as.character(plexes)]),
  n_proteins  = colSums(present_by_plex),
  subtypes    = sapply(plexes, function(p)
    paste(sort(table(ss$subtype[ss$tmt_plex == p])[
      table(ss$subtype[ss$tmt_plex == p]) > 0]) %>%
        { paste0(names(.), "=", .) }, collapse = ",")))

# NB: build the coverage counts BEFORE the tibble so `present_n_plex` resolves to
# the per-protein vector (tibble data-masking would otherwise shadow it with 0:5).
cover_counts <- as.integer(table(factor(present_n_plex, levels = 0:length(plexes))))
tier_tab <- tibble(present_n_plex = 0:length(plexes), n_proteins = cover_counts)
n_absent_any <- sum(block$absent_from_any_plex)
n_complete   <- sum(block$present_n_plex == length(plexes))

cat("\n=== TASK 3: structural block-missingness ===\n")
cat(sprintf("Proteins total: %d | present in all %d plexes (complete): %d (%.1f%%) | absent from >=1 whole plex: %d (%.1f%%)\n",
            nrow(block), length(plexes), n_complete, 100 * n_complete / nrow(block),
            n_absent_any, 100 * n_absent_any / nrow(block)))
cat("\nPer-protein plex-coverage distribution (# plexes a protein is quantified in):\n")
print(as.data.frame(tier_tab), row.names = FALSE)
cat("\nPer-plex summary (a protein absent from a plex is missing for ALL its lines):\n")
print(as.data.frame(plex_summary), row.names = FALSE)
cat(sprintf("\nMedian lines lost per non-complete protein: %d of %d\n",
            median(n_lines_total - block$present_n_lines[block$absent_from_any_plex]), n_lines_total))

# =============================================================================
# 4. TASK 2 — bridge chain topology + inter-plex normalization + precision
# =============================================================================
bc <- readr::read_csv(file.path(OUT, "prot_bridge_cor.csv"), show_col_types = FALSE)
is_chain <- all(bc$bridge_plex == bc$prim_plex + 1)              # consecutive links only
cat("\n=== TASK 2: inter-plex bridge design ===\n")
cat(sprintf("Bridge links (n=%d): %s\n", nrow(bc),
            paste(sprintf("plex%d<->%d[%s%s r=%.3f]", bc$prim_plex, bc$bridge_plex,
                          bc$cell_line, ifelse(bc$external, "*", ""),
                          bc$pearson), collapse = "  ")))
cat(sprintf("Topology: %s (each link connects CONSECUTIVE plexes -> daisy-CHAIN, not a common-reference hub)\n",
            ifelse(is_chain, "CHAIN confirmed", "NOT a simple chain")))
cat(sprintf("Per-link Pearson range: %.4f-%.4f (Spearman %.4f-%.4f). Max separation plex1<->plex%d = %d links.\n",
            min(bc$pearson), max(bc$pearson), min(bc$spearman), max(bc$spearman),
            length(plexes), length(plexes) - 1))
cat("* = external Carey LGS line VOA3993 (the plex3<->4 link). COV434 occupies plex1 ch10 as a\n",
    "  unique biological sample (no prior plex to bridge from), not a replicate.\n")

# Estimate bridge-derived precision from repeat-sample differences.
# sd(primary-bridge)/sqrt(2) approximates single-measurement SD under equal,
# independent measurement errors; it is pooled across proteins within each link.
PROT_DIR <- file.path(DATA, "proteomics")
abund <- read_excel(file.path(PROT_DIR, "protein_relative_abundance.xlsx"))
lay   <- read_excel(file.path(PROT_DIR, "tmt.layout.xlsx")) %>%
  transmute(id, plex, tmt_label = `TMT.label`, name)
data_cols <- names(abund)[12:ncol(abund)]
cmap <- tibble(data_col = data_cols) %>% left_join(lay, by = c("data_col" = "id")) %>%
  mutate(role = case_when(grepl("^SM\\.iD", data_col) ~ "SM.iD",
                          tmt_label == 10 & plex >= 2 ~ "bridge",
                          TRUE ~ "sample"))
brs  <- cmap %>% filter(role == "bridge")  %>% transmute(name, bridge_col = data_col)
prm  <- cmap %>% filter(role == "sample")  %>% transmute(name, prim_col = data_col)
bp   <- brs %>% left_join(prm, by = "name")
# The /sqrt(2) is the replicate-difference adjustment: sd(primary - bridge) is the
# SD of a DIFFERENCE of two measurements, so a per-measurement SD needs it. Section
# 4b deposits the same quantity per link and 4b-ii the like-for-like comparison
# against observed cross-line spread; this block is the console summary only.
noise_sd <- map_dbl(seq_len(nrow(bp)), function(i) {
  d <- as.numeric(abund[[bp$prim_col[i]]]) - as.numeric(abund[[bp$bridge_col[i]]])
  sd(d[is.finite(d)]) / sqrt(2)
})
tech_sd <- median(noise_sd)
cat(sprintf("\nBridge-derived precision: approximate per-measurement SD ~ %.3f log2 (per-link %s)\n",
            tech_sd, paste(sprintf("%.3f", noise_sd), collapse = ", ")))
cat(sprintf("Median observed cross-line protein SD = %.3f log2 -> observed:technical ~ %.1fx\n",
            med(dr$prot_sd), med(dr$prot_sd) / tech_sd))

# --- 4b. BRIDGE AGREEMENT, not bridge correlation ----------------------------
# A Pearson r on log2 abundances spanning ~12 log2 units is dominated by dynamic
# range, so high correlation does not establish agreement of the paired values.
# Compute the Bland-Altman/MA
# quantities that DO measure agreement, per link, with the n each rests on:
#   bias           = mean(primary - bridge)                (systematic offset)
#   sd_diff        = sd(primary - bridge)                  (repeatability, PAIRED scale)
#   loa_lower/upper= bias +/- 1.96*sd_diff                 (95% limits of agreement)
#   ma_slope       = slope of (primary-bridge) on mean     (abundance-dependent bias)
#
# TWO UNITS, AND THE CONVERSION BETWEEN THEM  [audit revision 2026-07-27]
# sd_diff is the SD of a DIFFERENCE of two measurements. If the two carry equal,
# independent error variance then Var(d) = 2*Var(single), so a PER-MEASUREMENT SD
# requires the replicate-difference adjustment
#       sd_single = sd_diff / sqrt(2).
# An earlier version of this script reported 100*(2^sd_diff - 1) and labelled it a
# "repeatability CV". That quantity is neither a per-measurement SD (no sqrt(2)
# adjustment) nor the standard lognormal CV, and it overstated per-measurement
# precision loss by ~1.5x. Both the adjusted SD and the correct lognormal CV,
#       cv_pct = 100 * sqrt(exp((ln2 * sd_single)^2) - 1),
# are now deposited; `sd_diff_cv_pct_legacy` is retained ONLY so the superseded
# number in earlier drafts is traceable, and is flagged deprecated in the file.
#
# WHAT THIS CAN AND CANNOT BE COMPARED WITH. The 95% LoA span is 2*1.96*sd_diff
# = 3.92 * sd_diff: an interval WIDTH on the paired-difference scale. A cross-line
# IQR is 1.35 * SD on the single-measurement scale. Dividing one by the other is
# not a noise-to-signal ratio - it inflates by ~2.9x purely from the choice of
# dispersion measure, which is exactly how the retired "technical noise exceeds
# biological spread by 2.4-3.1x" claim arose. The like-for-like comparison is
# SD against SD, and it is computed below.
ba <- purrr::map_dfr(seq_len(nrow(bp)), function(i) {
  x <- as.numeric(abund[[bp$prim_col[i]]]); y <- as.numeric(abund[[bp$bridge_col[i]]])
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  d <- x - y; m <- (x + y) / 2
  fit <- lm(d ~ m)
  sd_d <- sd(d); sd_1 <- sd_d / sqrt(2)
  tibble(cell_line = gsub("_", "-", bp$name[i]),
         n_proteins = length(d),
         bias = mean(d), sd_diff = sd_d,
         loa_lower = mean(d) - 1.96 * sd_d, loa_upper = mean(d) + 1.96 * sd_d,
         median_abs_diff = median(abs(d)),
         sd_single = sd_1,
         cv_pct_per_measurement = 100 * sqrt(exp((log(2) * sd_1)^2) - 1),
         sd_diff_cv_pct_legacy = 100 * (2^sd_d - 1),
         ma_slope = unname(coef(fit)[2]),
         ma_slope_p = summary(fit)$coefficients[2, 4],
         pearson = cor(x, y), spearman = cor(x, y, method = "spearman"))
})
ba <- bc %>% select(cell_line, prim_plex, bridge_plex, external) %>%
  left_join(ba, by = "cell_line") %>% arrange(bridge_plex) %>%
  mutate(loa_span = loa_upper - loa_lower,
         units_note = paste("sd_diff/loa_* are on the paired-difference scale;",
                            "sd_single and cv_pct_per_measurement are per measurement",
                            "(sd_single = sd_diff/sqrt(2));",
                            "sd_diff_cv_pct_legacy is DEPRECATED - superseded by cv_pct_per_measurement"),
         replicate_caveat = paste("four different bridge samples, one per link, each re-run in the",
                                  "adjacent plex; NOT replicated aliquots of one common reference"))
readr::write_csv(ba, file.path(OUT, "prot_bridge_agreement.csv"))
cat("\nBridge AGREEMENT (report these instead of, or alongside, Pearson r):\n")
print(as.data.frame(ba %>% select(cell_line, prim_plex, bridge_plex, external, n_proteins,
                                  bias, sd_diff, loa_span, sd_single, cv_pct_per_measurement,
                                  sd_diff_cv_pct_legacy, ma_slope_p, pearson) %>%
                      mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat(sprintf("Per-link n proteins: %d-%d. Bias %.3f to %+.3f log2 (negligible, no systematic offset).\n",
            min(ba$n_proteins), max(ba$n_proteins), min(ba$bias), max(ba$bias)))
cat(sprintf("  SD of paired differences %.3f-%.3f log2; 95%% LoA span %.3f-%.3f log2.\n",
            min(ba$sd_diff), max(ba$sd_diff), min(ba$loa_span), max(ba$loa_span)))
cat(sprintf("  Per-measurement SD %.3f-%.3f log2 (= %.1f-%.1f%% lognormal CV).\n",
            min(ba$sd_single), max(ba$sd_single),
            min(ba$cv_pct_per_measurement), max(ba$cv_pct_per_measurement)))
cat(sprintf("  [DEPRECATED, for traceability only: 100*(2^sd_diff-1) = %.1f-%.1f%%]\n",
            min(ba$sd_diff_cv_pct_legacy), max(ba$sd_diff_cv_pct_legacy)))

# --- 4b-i. BRIDGE DIFFERENCES, KEYED TO THE MATRIX ROW ID --------------------
# Built once here because both 4b-ii (ADC stratum) and 4c (abundance deciles) need it.
#
# PROTEIN KEY, FIXED  [audit revision 2026-07-27, second round]
# The matrix row identifier is `prot_qc.csv$row`: the bare symbol for the
# representative row of a symbol, and `SYMBOL|UNIPROT` for the 31 non-representative
# duplicate-symbol rows (see 05_proteomics_load_qc.R). The previous version keyed the
# abundance table on that decorated identifier but keyed the bridge differences on the
# RAW `abund$Symbol`, so the join was a MISS for every decorated row and a
# many-to-many FAN-OUT for the 44 workbook rows sharing one of 13 symbols, each
# difference being replicated across every row carrying that symbol. A
# `relationship = "many-to-many"` argument silenced the warning that would have caught
# it. `(Symbol, Uniprot)` is unique on both sides (asserted), so differences now
# resolve one-to-one to the canonical row id.
pq <- readr::read_csv(file.path(OUT, "prot_qc.csv"), show_col_types = FALSE)
sym2row <- pq %>% transmute(Symbol = symbol, Uniprot = uniprot, protein = row)
stopifnot("(symbol, uniprot) must be unique in prot_qc.csv" =
            !any(duplicated(paste(sym2row$Symbol, sym2row$Uniprot))),
          "(Symbol, Uniprot) must be unique in the source workbook" =
            !any(duplicated(paste(abund$Symbol, abund$Uniprot))))
bridge_diff <- purrr::map_dfr(seq_len(nrow(bp)), function(i) {
  x <- as.numeric(abund[[bp$prim_col[i]]]); y <- as.numeric(abund[[bp$bridge_col[i]]])
  tibble(Symbol = abund$Symbol, Uniprot = abund$Uniprot,
         link = gsub("_", "-", bp$name[i]), diff = x - y)
}) %>%
  filter(is.finite(diff)) %>%
  inner_join(sym2row, by = c("Symbol", "Uniprot"), relationship = "many-to-one") %>%
  select(protein, link, diff)
stopifnot("every bridge difference must resolve to exactly one matrix row" =
            !any(is.na(bridge_diff$protein)))

# --- 4b-ii. LIKE-FOR-LIKE technical vs observed spread, on ONE scale ---------
# Both quantities are per-measurement SDs in log2 units, so their ratio is
# interpretable - which the retired LoA-span-over-IQR ratio was not.
#
# WHAT THIS IS AND IS NOT  [audit revision 2026-07-27, second round]
# `approx_reliability_ratio` is a variance-ratio DIAGNOSTIC, not an intraclass
# correlation. It was previously named `reliability_icc`, which overclaimed: an ICC
# is estimated from a variance-components model fitted to replicated measurements of
# the SAME units, and nothing here is. Four approximations stand between this number
# and an ICC:
#   1. HOMOSCEDASTICITY. One pooled technical SD is applied to the median gene, but
#      technical error is a strong function of abundance (bridge SD of differences
#      falls ~4x from decile 1 to decile 10). Section 4c now reports the
#      abundance-MATCHED version, which drops this assumption and is the better
#      estimate; this global figure is an approximation to that gradient.
#   2. MISMATCHED GRANULARITY. The technical SD is pooled ACROSS proteins within a
#      link; the observed SD is per-gene and then medianed ACROSS genes. They are not
#      the same protein-level quantity.
#   3. NON-IDENTICAL REPLICATES. The four links are four different samples, one
#      re-run each, not replicated aliquots of one common reference.
#   4. INDEPENDENT, ADDITIVE ERROR. Var(observed) = Var(biological) + Var(technical)
#      is assumed, so the biological SD is reached by subtraction rather than
#      estimated. The observed SD is not a pure biological quantity.
# Defensible phrasing: "under a homoscedastic independent-error approximation, bridge
# variability corresponds to roughly X% of the variance of a typical protein's
# observed cross-line measurements." NOT: "a single protein in a single model has
# reliability X."
obs_sd   <- med(dr$prot_sd)
obs_iqr  <- med(dr$prot_iqr)
APPROX_NOTE <- paste(
  "variance-ratio diagnostic, NOT an intraclass correlation: pooled-across-proteins",
  "technical SD vs median-across-genes observed SD, under a homoscedastic",
  "independent-additive-error approximation; four non-identical bridge samples.",
  "See prot_cv_by_abundance.csv for the abundance-matched version.")
noise_tbl <- tibble(
  technical_sd_log2 = c(min(ba$sd_single), median(ba$sd_single), max(ba$sd_single)),
  bound = c("min link", "median link", "max link")) %>%
  mutate(stratum = "all paired genes", .before = 1) %>%
  mutate(technical_sd_basis = "pooled across proteins within a link",
         observed_cross_line_sd_log2 = obs_sd,
         observed_cross_line_iqr_log2 = obs_iqr,
         observed_to_technical_sd_ratio = obs_sd / technical_sd_log2,
         technical_variance_share_pct = 100 * technical_sd_log2^2 / obs_sd^2,
         implied_biological_sd_log2 = sqrt(pmax(obs_sd^2 - technical_sd_log2^2, 0)),
         approx_reliability_ratio = pmax(obs_sd^2 - technical_sd_log2^2, 0) / obs_sd^2,
         n_genes = nrow(dr),
         approximation_note = APPROX_NOTE)

# ADC TARGETS, with their OWN technical estimate  [second round]
# The earlier ADC comparison swapped in the ADC-target median observed SD but kept
# the GLOBAL technical SD, so it was not a per-class statement. The bridge
# differences restricted to the ADC target rows give a direct estimate; n is small
# (a few tens of differences), so it is reported WITH its n rather than as a
# precise value, and the abundance-matched decile figure in 4c is the cross-check.
adc_rows      <- dr$gene[dr$is_adc_target]
adc_diff_vals <- bridge_diff$diff[bridge_diff$protein %in% adc_rows]
adc_obs_sd    <- med(dr$prot_sd[dr$is_adc_target])
adc_sd_single <- sd(adc_diff_vals) / sqrt(2)
noise_tbl <- bind_rows(
  noise_tbl,
  tibble(stratum = sprintf("ADC targets (direct bridge estimate, n = %d differences)",
                           length(adc_diff_vals)),
         bound = "pooled over links",
         technical_sd_log2 = adc_sd_single,
         technical_sd_basis = "bridge differences restricted to the ADC target rows",
         observed_cross_line_sd_log2 = adc_obs_sd,
         observed_cross_line_iqr_log2 = med(dr$prot_iqr[dr$is_adc_target]),
         observed_to_technical_sd_ratio = adc_obs_sd / adc_sd_single,
         technical_variance_share_pct = 100 * adc_sd_single^2 / adc_obs_sd^2,
         implied_biological_sd_log2 = sqrt(pmax(adc_obs_sd^2 - adc_sd_single^2, 0)),
         approx_reliability_ratio =
           pmax(adc_obs_sd^2 - adc_sd_single^2, 0) / adc_obs_sd^2,
         n_genes = sum(dr$is_adc_target),
         approximation_note = paste("technical SD estimated from the ADC target rows",
                                    "themselves; few differences and non-identical bridge samples;",
                                    "variance-ratio approximation, not an ICC or a variance-component estimate")))
readr::write_csv(noise_tbl, file.path(OUT, "prot_noise_vs_biology.csv"))
cat("\nLIKE-FOR-LIKE technical vs observed cross-line spread (SD against SD, log2 units):\n")
print(as.data.frame(noise_tbl %>% select(stratum, bound, technical_sd_log2,
                                          observed_cross_line_sd_log2,
                                          observed_to_technical_sd_ratio,
                                          technical_variance_share_pct,
                                          implied_biological_sd_log2,
                                          approx_reliability_ratio) %>%
                      mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)
gl <- noise_tbl %>% filter(stratum == "all paired genes")
cat(sprintf("  All paired genes: observed cross-line SD is %.1f-%.1fx the per-measurement technical SD.\n",
            min(gl$observed_to_technical_sd_ratio), max(gl$observed_to_technical_sd_ratio)))
cat(sprintf("  Under a homoscedastic independent-error approximation, bridge variability corresponds to\n"))
cat(sprintf("  %.0f-%.0f%% of the variance of a typical protein's observed cross-line measurements.\n",
            min(gl$technical_variance_share_pct), max(gl$technical_variance_share_pct)))
cat(sprintf("  ADC targets, own technical estimate (n = %d differences): %.0f%%.\n",
            length(adc_diff_vals), noise_tbl$technical_variance_share_pct[nrow(noise_tbl)]))
cat("  This is a variance-ratio diagnostic, not an ICC. Section 4c drops the homoscedasticity\n")
cat("  assumption by matching the technical estimate to the observed spread within abundance deciles.\n")

# --- 4c. CV / noise stratified by ABUNDANCE ----------------------------------
# Observed spread and measurement variability can depend on abundance; a median CV hides
# the part a reuser needs. Two stratifications by mean protein abundance decile:
#   vendor_cv_pct      = the Morin-provided "CV replicates" column (prot_qc.csv)
#   bridge_sd_diff     = SD of the bridge primary-minus-bridge difference
# Both are reported per decile so a shortlist built at low abundance can be
# discounted appropriately.
# ZERO-PLEX EXCLUSION, MADE EXPLICIT [integration revision]. rowMeans(na.rm = TRUE)
# over a row that is NA in all 31 lines returns NaN, which the old
# `filter(is.finite(mean_abund))` then dropped — the right answer, but reached
# incidentally, so nothing would have noticed if a zero-plex row had ever acquired a
# value. The 70 zero-plex rows are now removed by name and the identity of the two
# sets is asserted, so the decile denominator (8,357 = 8,427 - 70) is deliberate.
prot_mean_ab <- rowMeans(full_mat[!rownames(full_mat) %in% prot_zero_plex, , drop = FALSE],
                         na.rm = TRUE)
stopifnot("dropping the zero-plex rows must leave every mean abundance finite" =
            all(is.finite(prot_mean_ab)),
          "the non-finite mean abundances must be exactly the zero-plex rows" =
            setequal(prot_zero_plex,
                     rownames(full_mat)[!is.finite(rowMeans(full_mat, na.rm = TRUE))]))
dec_of <- function(v) cut(v, quantile(v, seq(0, 1, 0.1), na.rm = TRUE),
                          labels = 1:10, include.lowest = TRUE)
ab_tbl <- tibble(protein = names(prot_mean_ab), mean_abund = prot_mean_ab) %>%
  mutate(abundance_decile = dec_of(mean_abund)) %>%
  left_join(pq %>% transmute(protein = row, vendor_cv_pct = cv_replicates), by = "protein")

# `bridge_diff` is keyed to the matrix row id in 4b-i, so the decile join below is
# one-to-one rather than the many-to-many fan-out it used to be.
bridge_by_decile <- bridge_diff %>%
  inner_join(ab_tbl %>% select(protein, abundance_decile), by = "protein",
             relationship = "many-to-one") %>%
  filter(!is.na(abundance_decile)) %>%
  group_by(abundance_decile) %>%
  summarise(n_bridge_diffs = dplyr::n(),
            bridge_sd_diff = sd(diff),
            bridge_median_abs_diff = median(abs(diff)), .groups = "drop") %>%
  # same replicate-difference adjustment as 4b: the decile SD is on the
  # paired-difference scale, so it needs /sqrt(2) before it is comparable
  # with a vendor per-measurement CV plotted on the same axis.
  mutate(bridge_sd_single = bridge_sd_diff / sqrt(2),
         bridge_cv_pct_per_measurement =
           100 * sqrt(exp((log(2) * bridge_sd_single)^2) - 1))

# ABUNDANCE-MATCHED noise-to-signal  [audit revision 2026-07-27, second round]
# The global comparison in 4b-ii applies ONE pooled technical SD to the median gene,
# which assumes homoscedastic error across the abundance range. It is not: the bridge
# SD of differences falls ~4-fold from decile 1 to decile 10. Matching the technical
# estimate to the observed spread WITHIN each decile removes that assumption, and it
# is the version a reuser should consult for a specific protein, because they know
# their protein's abundance. `prot_sd_median_decile` is the median cross-line SD of
# the paired genes falling in that decile, so both sides are decile-local.
prot_sd_by_row <- dr %>%
  transmute(protein = gene, prot_sd) %>%          # dr is keyed on the matrix row id
  inner_join(ab_tbl %>% select(protein, abundance_decile), by = "protein")
cv_by_decile <- ab_tbl %>% group_by(abundance_decile) %>%
  summarise(n_proteins = dplyr::n(),
            mean_abund_lo = min(mean_abund), mean_abund_hi = max(mean_abund),
            vendor_cv_median = median(vendor_cv_pct, na.rm = TRUE),
            vendor_cv_q25 = quantile(vendor_cv_pct, .25, na.rm = TRUE),
            vendor_cv_q75 = quantile(vendor_cv_pct, .75, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(bridge_by_decile, by = "abundance_decile") %>%
  left_join(prot_sd_by_row %>% group_by(abundance_decile) %>%
              summarise(n_paired_genes_decile = dplyr::n(),
                        prot_sd_median_decile = median(prot_sd), .groups = "drop"),
            by = "abundance_decile") %>%
  mutate(observed_to_technical_sd_ratio = prot_sd_median_decile / bridge_sd_single,
         technical_variance_share_pct = 100 * bridge_sd_single^2 / prot_sd_median_decile^2,
         approx_reliability_ratio =
           pmax(prot_sd_median_decile^2 - bridge_sd_single^2, 0) / prot_sd_median_decile^2)
readr::write_csv(cv_by_decile, file.path(OUT, "prot_cv_by_abundance.csv"))
cat("\nABUNDANCE-MATCHED technical vs observed spread (both decile-local, log2 SD):\n")
print(as.data.frame(cv_by_decile %>%
                      select(abundance_decile, n_paired_genes_decile, bridge_sd_single,
                             prot_sd_median_decile, observed_to_technical_sd_ratio,
                             technical_variance_share_pct, approx_reliability_ratio) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat(sprintf("  Technical share of observed cross-line variance: %.0f%% (decile 1) -> %.0f%% (decile 10);\n",
            cv_by_decile$technical_variance_share_pct[1],
            cv_by_decile$technical_variance_share_pct[10]))
cat("  the global figure in 4b-ii is a homoscedastic approximation to this gradient.\n")
cat("\nCV / bridge noise by mean-abundance decile (1 = lowest abundance):\n")
print(as.data.frame(cv_by_decile %>% mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)

# =============================================================================
# 5. FIGURE f_prot_compression (4 panels)
# =============================================================================
# Panel A: paired per-gene cross-line spread (IQR), RNA vs protein -------------
long_iqr <- dr %>% select(gene, RNA = rna_iqr, protein = prot_iqr) %>%
  pivot_longer(c(RNA, protein), names_to = "assay", values_to = "iqr") %>%
  mutate(assay = factor(assay, levels = c("RNA", "protein")))
medlab <- long_iqr %>% group_by(assay) %>% summarise(m = median(iqr), .groups = "drop")
pA <- ggplot(long_iqr, aes(assay, iqr, fill = assay)) +
  geom_violin(colour = NA, alpha = 0.85, width = 0.85, scale = "width") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", colour = "black", linewidth = 0.35) +
  geom_text(data = medlab, aes(assay, m, label = sprintf("%.2f", m)),
            hjust = -0.9, size = 3, colour = "black") +
  scale_fill_manual(values = c(RNA = COL_RNA, protein = COL_PROT), guide = "none") +
  coord_cartesian(ylim = c(0, quantile(long_iqr$iqr, 0.99, na.rm = TRUE))) +
  labs(title = "A  Cross-line spread is narrower at the protein level",
       subtitle = sprintf("per-gene IQR across matched lines; n = %d genes", nrow(dr)),
       x = NULL, y = expression("cross-line IQR ("*log[2]*" units)")) +
  theme_lab()

# Panel B: distribution of the protein/RNA IQR ratio ---------------------------
mr <- med(dr$iqr_ratio)
XMAX <- 1.5                                             # display cap (median 0.30; =1 reference visible)
frac_beyond <- mean(dr$iqr_ratio > XMAX, na.rm = TRUE)
pB <- ggplot(dr %>% filter(is.finite(iqr_ratio), iqr_ratio <= XMAX), aes(iqr_ratio)) +
  geom_histogram(binwidth = 0.025, fill = "grey78", colour = "white", boundary = 0) +
  geom_vline(xintercept = 1, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = mr, colour = COL_PROT, linewidth = 0.9) +
  annotate("text", x = mr, y = Inf, label = sprintf("median = %.2f", mr),
           vjust = 1.6, hjust = -0.08, colour = COL_PROT, size = 3) +
  annotate("text", x = 1, y = Inf, label = "equal spread (=1)",
           vjust = 1.6, hjust = 1.05, colour = "grey40", size = 2.7) +
  coord_cartesian(xlim = c(0, XMAX)) +
  labs(title = "B  Protein / RNA dynamic-range ratio",
       subtitle = sprintf("%.0f%% of genes: protein spread < RNA spread (%.0f%% beyond x-axis cap)",
                          100 * frac_compressed, 100 * frac_beyond),
       x = "protein IQR / RNA IQR (per gene)", y = "genes") +
  theme_lab()

# Panel C: ADC-target exemplars — RNA vs protein per-line spread (FOLR1 marked) -
adc_show <- intersect(c("FOLR1", "MSLN", "ERBB2", "TACSTD2", "CDH6", "SLC34A2", "CD276", "VTCN1", "DPEP3"),
                      rownames(P))
cdat <- lapply(adc_show, function(g) {
  ok <- is.finite(P[g, ]) & is.finite(R[g, ])
  bind_rows(
    tibble(gene = g, assay = "RNA",     line = shared_lines[ok], val = R[g, ok]),
    tibble(gene = g, assay = "protein", line = shared_lines[ok], val = P[g, ok]))
}) %>% bind_rows() %>%
  mutate(assay = factor(assay, levels = c("RNA", "protein")),
         gene  = factor(gene, levels = adc_show))
seg <- cdat %>% group_by(gene, assay) %>%
  summarise(lo = min(val), hi = max(val), .groups = "drop")
pC <- ggplot(cdat, aes(val, gene, colour = assay)) +
  geom_segment(data = seg, aes(x = lo, xend = hi, y = gene, yend = gene),
               linewidth = 0.5, alpha = 0.5,
               position = position_dodge(width = 0.55)) +
  geom_point(size = 1.5, alpha = 0.7, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(RNA = COL_RNA, protein = COL_PROT), name = NULL) +
  labs(title = "C  ADC-target range: wide in RNA, narrow in protein",
       subtitle = "each point = 1 model; assay scales differ; abundance is not surface expression",
       x = expression("expression ("*log[2]*" units, RNA & protein overlaid)"), y = NULL) +
  theme_lab() + theme(legend.position = "top", legend.justification = "left")

# Panel D: block-missingness — per-protein plex coverage -----------------------
pD <- ggplot(tier_tab %>% mutate(present_n_plex = factor(present_n_plex)),
             aes(present_n_plex, n_proteins)) +
  geom_col(width = 0.72, fill = COL_PROT) +
  geom_text(aes(label = n_proteins), vjust = -0.4, size = 2.8) +
  annotate("text", x = 0.6, y = Inf,
           label = sprintf("%d proteins (%.0f%%) absent from >=1 whole plex",
                           n_absent_any, 100 * n_absent_any / nrow(block)),
           vjust = 1.5, hjust = 0, size = 2.8, colour = COOK_NAVY) +
  labs(title = "D  Missingness is structural (per-plex blocks)",
       subtitle = "a protein absent from a plex is missing for ALL lines in it",
       x = "# plexes a protein is quantified in (0-5)", y = "proteins") +
  coord_cartesian(clip = "off") + theme_lab()

fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1.05))
ggsave(file.path(FIGS, "f_prot_compression.pdf"), fig, width = 11, height = 8.4)
ggsave(file.path(PROJ, "reports", "assets", "f_prot_compression.png"),
       fig, width = 11, height = 8.4, dpi = 200)

cat("\nOutputs written:\n")
cat("  output/prot_dynamic_range.csv     (", nrow(dr), "genes )\n")
cat("  output/prot_block_missingness.csv (", nrow(block), "proteins )\n")
cat("  output/prot_compression_floor_check.csv, prot_bridge_agreement.csv, prot_cv_by_abundance.csv\n")
cat("  figs/f_prot_compression.pdf + reports/assets/f_prot_compression.png\n")

# =============================================================================
# 6. Environment record
# =============================================================================
# 00_setup.R checks that required packages are PRESENT but records no versions
# and there is no renv lockfile yet, so each script writes its own session record
# beside its outputs. Reviewers asked for this explicitly.
write_session_info <- function(script) {
  si <- utils::sessionInfo(); pk <- c(si$otherPkgs, si$loadedOnly)
  f  <- file.path(OUT, paste0("session_info_", script, ".txt"))
  writeLines(c(sprintf("script: %s", script),
               sprintf("run_at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               sprintf("seed: %s", SEED), "", "--- sessionInfo() ---",
               utils::capture.output(print(si)),
               "", "--- package versions (attached + loaded) ---",
               sort(sprintf("%-28s %s", names(pk),
                            vapply(pk, function(p) as.character(p$Version), character(1))))), f)
  message("  session record -> ", f)
}
write_session_info("19_proteomics_dynamic_range")

message("\n19_proteomics_dynamic_range.R complete.")
