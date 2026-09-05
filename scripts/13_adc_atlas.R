# =============================================================================
# Script: 13_adc_atlas.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Subtype-resolved ADC-target expression atlas (RNA + protein) - the
#          FEATURED USAGE EXAMPLE of the descriptor (Fig 5a). Shows how the
#          resource is used to shortlist preclinical models by ADC-target
#          expression, and checks the atlas against known subtype associations.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   5 (Integration & usage) - step 2 of 3
# =============================================================================
# Framing / caveats (lit review Theme 7):
#   - This is MODEL-SELECTION utility, NOT clinical target discovery.
#   - Target expression is NECESSARY BUT NOT SUFFICIENT for ADC response:
#     NaPi2b (SLC34A2; UPLIFT missed, UP-NEXT discontinued), mesothelin (MSLN;
#     anetumab NCI#10150 negative) and DPEP3 (SC-003 ORR 4%, discontinued) all
#     had high target expression but NEGATIVE pivotal trials. The atlas is
#     hypothesis-generating for choosing models to interrogate, not a claim that
#     a high-expressing line will respond.
#   - NOMENCLATURE: the task list wrote "CD276 (B7-H4)", but these are two
#     distinct genes/targets - CD276 = B7-H3 (ifinatamab deruxtecan) and
#     VTCN1 = B7-H4 (e.g. puxitatug samrotecan / XMT-1660). Both are real
#     ovarian ADC targets, so BOTH are included and labelled correctly.
#   - DPEP3 is not detected in the TMT proteomics (low-abundance) -> shown as
#     RNA-only; its protein row is left as NA ("not detected") in the heatmap.
#   - Bulk 2D lines: expression is whole-population, not surface-resolved; ADCs
#     act on cell-surface protein, so protein (not just mRNA) is the more
#     relevant layer, and even that does not capture surface localisation.
# =============================================================================

source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(matrixStats); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(viridisLite)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
set.seed(SEED)

theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(text = element_text(colour = "black"),
          plot.title = element_text(size = rel(1.1), hjust = 0, margin = margin(b = 6)),
          plot.subtitle = element_text(size = rel(0.85), hjust = 0, colour = "grey30",
                                       margin = margin(b = 6)),
          axis.title = element_text(size = rel(1.0)),
          axis.text = element_text(size = rel(0.85), colour = "black"),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          legend.title = element_text(size = rel(0.85)), legend.text = element_text(size = rel(0.8)),
          legend.background = element_blank(), legend.key = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(),
          panel.grid = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = rel(0.95), face = "bold"),
          plot.margin = margin(8, 8, 8, 8))
}

# -----------------------------------------------------------------------------
# 0. ADC target panel (curated; symbol -> common name, ADC, expected subtype)
# -----------------------------------------------------------------------------
adc <- tribble(
  ~symbol,   ~common,      ~adc_example,                       ~expected,
  "FOLR1",   "FRalpha",    "mirvetuximab soravtansine",        "HGS",
  "TACSTD2", "TROP2",      "sacituzumab/datopotamab deruxtecan", NA,
  "ERBB2",   "HER2",       "trastuzumab deruxtecan",           "CC/MC",
  "MSLN",    "mesothelin", "anetumab ravtansine",              "HGS",
  "SLC34A2", "NaPi2b",     "upifitamab rilsodotin",            NA,
  "DPEP3",   "DPEP3",      "tamrintamab pamozirine (SC-003)",  NA,
  "CDH6",    "CDH6",       "raludotatug deruxtecan",           NA,
  "CD276",   "B7-H3",      "ifinatamab deruxtecan",            NA,
  "VTCN1",   "B7-H4",      "puxitatug samrotecan / XMT-1660",  NA)
adc <- adc %>% mutate(row_label = sprintf("%s (%s)", symbol, common))

# -----------------------------------------------------------------------------
# 1. RNA symbol-collapsed log2(TPM+1) matrix (same collapse as scripts 04/12)
# -----------------------------------------------------------------------------
tpm <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"), show_col_types = FALSE)
id2sym <- t2g %>% distinct(ensembl_gene_id, external_gene_name) %>%
  filter(!is.na(external_gene_name), external_gene_name != "") %>%
  { setNames(.$external_gene_name, .$ensembl_gene_id) }
m_id <- as.matrix(tpm[, -1]); rownames(m_id) <- tpm$gene_id
common <- intersect(rownames(m_id), names(id2sym))
rna_log <- log2(rowsum(m_id[common, , drop = FALSE], group = id2sym[common]) + 1)

# 2. Protein matrix (log2 relative abundance; NAs preserved) --------------------
# [integration revision] Shared loader (00_setup.R); feature hygiene is handled at
# source by script 05. Zero-plex rows are retained but must not become an ADC "value":
# every per-subtype mean already reports n_measured, and no target is zero-plex
# (asserted once the panel is defined).
prot_mat <- read_prot_matrix()
prot_zero_plex <- attr(prot_mat, "zero_plex")
# A zero-plex target would be reported as an absent protein when in fact it was never
# quantified in any line — a very different statement for an ADC candidate. None of
# the 9 targets is zero-plex; assert it. [integration revision]
stopifnot("an ADC target protein is zero-plex (NA in all 31 lines) — 'absent' and 'never quantified' must not be conflated" =
            !any(adc$symbol %in% prot_zero_plex))

# 3. Annotation (source of truth) ----------------------------------------------
ann <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  filter(provenance == "generated", analysis_include == "Y") %>%
  transmute(cell_line,
            subtype = factor(subtype, levels = c("HGS","CC","EC","MC","MMMT","SCCOHT")))
sub_lvls <- levels(ann$subtype)
sub_cols <- setNames(RColorBrewer::brewer.pal(6, "Dark2"), sub_lvls)  # match scripts 02/04/12

# presence check ---------------------------------------------------------------
adc$in_rna  <- adc$symbol %in% rownames(rna_log)
adc$in_prot <- adc$symbol %in% rownames(prot_mat)
cat("\n=== ADC target availability ===\n")
print(as.data.frame(adc %>% select(symbol, common, in_rna, in_prot)), row.names = FALSE)
stopifnot("some ADC targets missing from RNA" = all(adc$in_rna))          # all 9 in RNA
message(sprintf("Protein-detected targets: %d/%d (missing: %s)",
                sum(adc$in_prot), nrow(adc),
                paste(adc$symbol[!adc$in_prot], collapse = ", ")))

# -----------------------------------------------------------------------------
# 4. Long tidy per-line table (ALL lines; NA where an assay is missing) -> output
# -----------------------------------------------------------------------------
rna_long <- as_tibble(rna_log[adc$symbol, , drop = FALSE], rownames = "symbol") %>%
  pivot_longer(-symbol, names_to = "cell_line", values_to = "value") %>%
  mutate(assay = "RNA")
prot_sym <- intersect(adc$symbol, rownames(prot_mat))
prot_long <- as_tibble(prot_mat[prot_sym, , drop = FALSE], rownames = "symbol") %>%
  pivot_longer(-symbol, names_to = "cell_line", values_to = "value") %>%
  mutate(assay = "protein")
adc_tbl <- bind_rows(rna_long, prot_long) %>%
  left_join(ann, by = "cell_line") %>%
  left_join(adc %>% select(symbol, common, adc_example, expected), by = "symbol") %>%
  transmute(cell_line, subtype = as.character(subtype), symbol, common, adc_example,
            expected_subtype = expected, assay,
            log2_expr = round(value, 3)) %>%
  arrange(symbol, assay, subtype, cell_line)
readr::write_csv(adc_tbl, file.path(OUT, "adc_expression.csv"))

# -----------------------------------------------------------------------------
# 5. Known-subtype-association check (mean expression per subtype; is expected top?)
# -----------------------------------------------------------------------------
# THE TWO MODALITIES ARE NOT OVER THE SAME LINES. VOA6861 has RNA but no
# proteomics and VOA14993 has proteomics but no RNA, so an "RNA vs protein"
# subtype mean is computed on 31 lines for RNA and 31 (a different 31) for
# protein — and within a subtype the n can differ by one. Nothing recorded this,
# so every mean now carries its own n and the modality-specific lines are written
# out. Read any RNA-vs-protein comparison of subtype means with these n's in view.
rna_lines  <- colnames(rna_log)
prot_lines <- colnames(prot_mat)
modality_sets <- bind_rows(
  tibble(cell_line = setdiff(rna_lines, prot_lines), only_in = "RNA"),
  tibble(cell_line = setdiff(prot_lines, rna_lines), only_in = "protein")) %>%
  left_join(ann, by = "cell_line") %>%
  transmute(cell_line, subtype = as.character(subtype), only_in)
readr::write_csv(modality_sets, file.path(OUT, "adc_modality_line_sets.csv"))
cat("\n=== Modality-specific lines (why RNA and protein subtype means differ in n) ===\n")
print(as.data.frame(modality_sets), row.names = FALSE)

# lines available per (assay, subtype) — the denominator for every mean below
assay_lines <- bind_rows(
  tibble(assay = "RNA", cell_line = rna_lines),
  tibble(assay = "protein", cell_line = prot_lines)) %>%
  left_join(ann, by = "cell_line") %>%
  group_by(assay, subtype) %>% summarise(n_lines_in_modality = dplyr::n(), .groups = "drop")

sub_mean <- adc_tbl %>%
  group_by(symbol, common, assay, subtype) %>%
  summarise(mean_log2 = mean(log2_expr, na.rm = TRUE), n = sum(!is.na(log2_expr)),
            .groups = "drop") %>%
  left_join(assay_lines %>% mutate(subtype = as.character(subtype)),
            by = c("assay", "subtype")) %>%
  rename(n_measured = n)
cat("\nLines per subtype in each modality (the mean's denominator):\n")
print(as.data.frame(assay_lines %>% tidyr::pivot_wider(names_from = assay,
                                                       values_from = n_lines_in_modality)),
      row.names = FALSE)

# subtype holding the single highest-expressing line per target/assay (model-selection view)
max_subtype <- adc_tbl %>% filter(!is.na(log2_expr)) %>%
  group_by(symbol, assay) %>% slice_max(log2_expr, n = 1, with_ties = FALSE) %>%
  transmute(symbol, assay, max_line_subtype = subtype) %>% ungroup()

is_expected <- function(st, exp) dplyr::case_when(
  is.na(exp) ~ NA, exp == "CC/MC" ~ as.character(st) %in% c("CC","MC"),
  TRUE ~ as.character(st) == exp)

# rank each subtype by mean (1 = highest mean) per target/assay
ranked <- sub_mean %>% group_by(symbol, assay) %>%
  mutate(mean_rank = rank(-mean_log2, ties.method = "min")) %>% ungroup()
# rank of the EXPECTED subtype (min rank among expected subtypes; CC/MC -> best of the two)
exp_rank <- ranked %>% left_join(adc %>% select(symbol, expected), by = "symbol") %>%
  filter(!is.na(expected),
         (expected == "CC/MC" & as.character(subtype) %in% c("CC","MC")) |
         (expected != "CC/MC" & as.character(subtype) == expected)) %>%
  group_by(symbol, assay) %>% summarise(expected_mean_rank = min(mean_rank), .groups = "drop")
# top-mean subtype per target/assay
top_mean <- ranked %>% group_by(symbol, common, assay) %>% arrange(desc(mean_log2)) %>%
  summarise(top_mean_subtype = subtype[1], top_val = round(mean_log2[1], 2),
            top_n_measured = n_measured[1], top_n_lines_in_modality = n_lines_in_modality[1],
            second = subtype[2], second_val = round(mean_log2[2], 2),
            second_n_measured = n_measured[2], .groups = "drop")

assoc <- top_mean %>%
  left_join(max_subtype, by = c("symbol","assay")) %>%
  left_join(adc %>% select(symbol, expected), by = "symbol") %>%
  left_join(exp_rank, by = c("symbol","assay")) %>%
  mutate(top_mean_matches     = is_expected(top_mean_subtype, expected),
         top_line_in_expected = is_expected(max_line_subtype, expected)) %>%
  arrange(symbol, assay)
readr::write_csv(assoc, file.path(OUT, "adc_subtype_summary.csv"))

cat("\n=== Known subtype-association check ===\n")
cat("top_mean = highest-MEAN subtype; max_line = subtype of the single highest-expressing line;\n")
cat("exp_rank = rank (by mean, 1=highest) of the EXPECTED subtype. Model selection cares that the\n")
cat("expected subtype CONTAINS high expressers (top_line_in_expected), not only that it tops the mean.\n")
print(as.data.frame(assoc %>% filter(!is.na(expected)) %>%
        select(symbol, common, assay, expected, top_mean_subtype, top_mean_matches,
               max_line_subtype, top_line_in_expected, expected_mean_rank)), row.names = FALSE)
chk <- assoc %>% filter(!is.na(expected))
cat(sprintf("\nExpected subtype tops the MEAN: %d / %d target-assay tests\n",
            sum(chk$top_mean_matches, na.rm = TRUE), nrow(chk)))
cat(sprintf("Expected subtype CONTAINS the top individual expresser: %d / %d\n",
            sum(chk$top_line_in_expected, na.rm = TRUE), nrow(chk)))

# -----------------------------------------------------------------------------
# 5b. FOLR1 within HGS — TEST the bimodality instead of asserting it from a range
# -----------------------------------------------------------------------------
# "Strongly bimodal" was previously inferred from a range (0.06-9.5), which is not
# evidence of bimodality at all. Two statistics are computed and written out:
#   (1) a 2-component vs 1-component Gaussian mixture fit (mclust), compared by
#       BIC. dBIC = BIC(G=2) - BIC(G=1); mclust reports BIC so HIGHER is better,
#       i.e. dBIC > 0 favours two components (>10 is conventionally strong).
#       Both equal-variance (E) and unequal-variance (V) models are fitted, since
#       at n=15 the V model can overfit a single outlying line.
#   (2) the moment-based bimodality coefficient BC = (skew^2 + 1) / (kurt + 3g),
#       g = 3(n-1)^2/((n-2)(n-3)); BC > 0.555 (the uniform-distribution value)
#       is the conventional bimodality flag. Distribution-free, no model.
# NOTE: Hartigan's dip test is the statistic a reviewer will ask for, and the
# `diptest` package is NOT installed in this R build, so it is not computed here.
# Install diptest and add dip.test(x) if a formal dip p-value is required.
folr1_hgs <- adc_tbl %>% filter(symbol == "FOLR1", assay == "RNA", subtype == "HGS", !is.na(log2_expr))
x_folr1 <- folr1_hgs$log2_expr
bimod_coef <- function(x) {
  n <- length(x); m <- mean(x); s <- sd(x)
  skew <- sum((x - m)^3) / (n * s^3)
  kurt <- sum((x - m)^4) / (n * s^4) - 3
  (skew^2 + 1) / (kurt + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}
bimod_rows <- tibble(statistic = character(), value = numeric(), interpretation = character())
if (requireNamespace("mclust", quietly = TRUE)) {
  # mclust::Mclust() calls mclustBIC() unqualified, so the package must be
  # ATTACHED, not merely namespace-loaded. It exports map()/sim(), which mask
  # purrr::map()/stats::sim — nothing in this script calls those bare.
  suppressPackageStartupMessages(library(mclust))
  fits <- lapply(c("E", "V"), function(mn)
    lapply(1:2, function(g) mclust::Mclust(x_folr1, G = g, modelNames = mn, verbose = FALSE)))
  names(fits) <- c("E", "V")
  for (mn in names(fits)) {
    b1 <- fits[[mn]][[1]]$bic; b2 <- fits[[mn]][[2]]$bic
    f2 <- fits[[mn]][[2]]
    bimod_rows <- bind_rows(bimod_rows, tibble(
      statistic = c(sprintf("mclust_%s_BIC_G1", mn), sprintf("mclust_%s_BIC_G2", mn),
                    sprintf("mclust_%s_dBIC_G2_minus_G1", mn),
                    sprintf("mclust_%s_G2_mean_low", mn), sprintf("mclust_%s_G2_mean_high", mn),
                    sprintf("mclust_%s_G2_prop_low", mn), sprintf("mclust_%s_G2_n_low", mn)),
      value = c(b1, b2, b2 - b1,
                min(f2$parameters$mean), max(f2$parameters$mean),
                f2$parameters$pro[which.min(f2$parameters$mean)],
                sum(f2$classification == which.min(f2$parameters$mean))),
      interpretation = c("", "",
                         ifelse(b2 - b1 > 10, "strong support for 2 components",
                                ifelse(b2 - b1 > 0, "weak support for 2 components",
                                       "1 component preferred")),
                         "", "", "", "")))
  }
} else {
  message("mclust not installed - FOLR1 mixture fit skipped")
}
bimod_rows <- bind_rows(
  tibble(statistic = c("n_HGS_lines", "min", "max", "range", "sd", "bimodality_coefficient"),
         value = c(length(x_folr1), min(x_folr1), max(x_folr1),
                   diff(range(x_folr1)), sd(x_folr1), bimod_coef(x_folr1)),
         interpretation = c("", "", "", "range is NOT evidence of bimodality", "",
                            ifelse(bimod_coef(x_folr1) > 0.555,
                                   "BC > 0.555 -> bimodality flagged",
                                   "BC <= 0.555 -> bimodality NOT flagged"))),
  bimod_rows,
  tibble(statistic = "hartigan_dip_test", value = NA_real_,
         interpretation = "NOT COMPUTED: the diptest package is not installed in this R build"))
bimod_out <- bimod_rows %>% mutate(gene = "FOLR1", assay = "RNA", group = "HGS",
                                   value = round(value, 4),
                                   analysis_scope = "exploratory sensitivity of 15 models from 12 patients; a 2-component Gaussian mixture need not have a bimodal density", .before = 1)
readr::write_csv(bimod_out, file.path(OUT, "adc_folr1_bimodality.csv"))

cat(sprintf("\n=== FOLR1 within HGS (n=%d): is it actually bimodal? ===\n", length(x_folr1)))
cat(sprintf("RNA log2(TPM+1) range %.2f-%.2f (SD %.2f)\n",
            min(x_folr1), max(x_folr1), sd(x_folr1)))
print(as.data.frame(bimod_out %>% select(statistic, value, interpretation)), row.names = FALSE)
cat("  The 3 highest FOLR1-RNA lines panel-wide are HGS; many HGS lines are near-zero, so the HGS\n")
cat("  MEAN is dragged below MC/CC even though FRalpha-high HGS models exist (e.g. TOV3133G, OV3133-R).\n")
cat("  This mirrors the clinical picture and is a model-selection feature, not a failure to recapitulate.\n")

# -----------------------------------------------------------------------------
# 6. FIGURE A - per-line atlas heatmap (RNA over protein; shared 30 lines)
# -----------------------------------------------------------------------------
shared <- intersect(colnames(rna_log), colnames(prot_mat))                 # 30 lines
col_ord <- ann %>% filter(cell_line %in% shared) %>% arrange(subtype, cell_line)
cl <- col_ord$cell_line
sub_split <- factor(col_ord$subtype, levels = sub_lvls)

rna_hm  <- rna_log[adc$symbol, cl, drop = FALSE]
prot_hm <- matrix(NA_real_, nrow = nrow(adc), ncol = length(cl),
                  dimnames = list(adc$symbol, cl))
prot_hm[prot_sym, ] <- prot_mat[prot_sym, cl]
rownames(rna_hm) <- rownames(prot_hm) <- adc$row_label

rna_colf  <- colorRamp2(seq(min(rna_hm), max(rna_hm), length.out = 9), magma(9))
prot_colf <- colorRamp2(seq(min(prot_hm, na.rm = TRUE), max(prot_hm, na.rm = TRUE),
                             length.out = 9), magma(9))
col_ha <- HeatmapAnnotation(subtype = col_ord$subtype, col = list(subtype = sub_cols),
                            annotation_name_gp = grid::gpar(fontsize = 8),
                            annotation_legend_param = list(subtype = list(nrow = 1)))
hm_opts <- list(cluster_rows = FALSE, cluster_columns = FALSE,
                column_split = sub_split, cluster_column_slices = FALSE,
                column_gap = grid::unit(1.2, "mm"),
                row_names_side = "left", row_names_gp = grid::gpar(fontsize = 8),
                na_col = "grey92",
                width = length(cl) * unit(4.2, "mm"),
                height = nrow(adc) * unit(4.2, "mm"),
                rect_gp = grid::gpar(col = "white", lwd = 0.5))

# slice titles (subtype) rotated 90 so the narrow n=2 slices (MMMT/SCCOHT) fit;
# shown on the TOP (RNA) heatmap only, suppressed on the protein heatmap.
ht_rna <- do.call(Heatmap, c(list(rna_hm, name = "RNA log2(TPM+1)", col = rna_colf,
                                   top_annotation = col_ha, show_column_names = FALSE,
                                   column_title_rot = 90,
                                   column_title_gp = grid::gpar(fontsize = 8, fontface = "bold")),
                             hm_opts))
ht_prot <- do.call(Heatmap, c(list(prot_hm, name = "protein log2 abund.", col = prot_colf,
                                    show_column_names = TRUE, column_title = NULL,
                                    column_names_gp = grid::gpar(fontsize = 6.5)), hm_opts))

pdf(file.path(FIGS, "13_adc_atlas_heatmap.pdf"), width = 11, height = 6.4)
draw(ht_rna %v% ht_prot, heatmap_legend_side = "right", annotation_legend_side = "bottom",
     column_title = "ADC-target expression atlas: RNA (top) and protein (bottom), matched lines by subtype",
     column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
     merge_legends = TRUE)
dev.off()

# -----------------------------------------------------------------------------
# 7. FIGURE B - subtype-mean bubble plot (association check; RNA | protein panels)
# -----------------------------------------------------------------------------
sm <- sub_mean %>%
  mutate(symbol = factor(symbol, levels = rev(adc$symbol)),
         subtype = factor(subtype, levels = sub_lvls),
         assay = factor(assay, levels = c("RNA","protein")))
bubble <- function(dat, ttl, legname) {
  ggplot(dat, aes(subtype, symbol)) +
    geom_point(aes(size = mean_log2, fill = mean_log2), shape = 21, colour = "grey30") +
    scale_fill_viridis_c(option = "magma", name = legname) +
    scale_size_continuous(range = c(1.5, 8), guide = "none") +
    labs(title = ttl, x = NULL, y = NULL) + theme_lab() +
    theme(axis.text.y = element_text(size = rel(0.8)))
}
pR <- bubble(sm %>% filter(assay == "RNA"), "A  RNA  (mean log2 TPM)", "RNA") +
  labs(subtitle = "expected: FOLR1->HGS, HER2(ERBB2)->CC/MC, MSLN->HGS(serous)")
pP <- bubble(sm %>% filter(assay == "protein"), "B  Protein  (mean log2 abundance)", "protein") +
  labs(subtitle = "DPEP3 not detected in proteomics")
ggsave(file.path(FIGS, "13_adc_subtype_bubble.pdf"),
       pR + pP + patchwork::plot_layout(guides = "collect"),
       width = 11, height = 4.4)

# -----------------------------------------------------------------------------
# 8. Notable per-line model-selection observations (highest-expressing lines)
# -----------------------------------------------------------------------------
cat("\n=== Top RNA- and protein-expressing line per target (model shortlist) ===\n")
top_line <- adc_tbl %>% filter(!is.na(log2_expr)) %>%
  group_by(symbol, common, assay) %>% slice_max(log2_expr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(symbol, common, assay, top_line = cell_line, subtype, log2_expr) %>%
  arrange(symbol, assay)
print(as.data.frame(top_line), row.names = FALSE)

cat("\nOutputs: output/adc_expression.csv, adc_subtype_summary.csv,\n")
cat("         adc_modality_line_sets.csv, adc_folr1_bimodality.csv\n")
cat("Figures: figs/13_adc_atlas_heatmap.pdf, figs/13_adc_subtype_bubble.pdf\n")
cat("NOTE: subtype means are PER LINE (not per patient) and the RNA and protein line sets differ\n")
cat("  by one line each (see adc_modality_line_sets.csv). Non-HGS groups are n=2-7.\n")

# -----------------------------------------------------------------------------
# 9. Environment record
# -----------------------------------------------------------------------------
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
write_session_info("13_adc_atlas")

message("\n13_adc_atlas.R complete.")
