# =============================================================================
# Script: 12_rna_protein_concordance.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: RNA-protein concordance across the matched lines, benchmarked to
#          published proteogenomic references. This is a Technical-Validation
#          element of the descriptor (Fig 3): does the resource show the
#          EXPECTED, moderate mRNA-protein correlation seen in tumours and other
#          cell-line panels? It REPLACES the thesis's 2-line correlation with a
#          proper distribution across all shared lines/genes.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   5 (Integration & usage) - step 1 of 3
# =============================================================================
# Assay-aware notes / design decisions:
#   - TWO distinct metrics, often conflated. We report both and label them:
#       (a) PER-GENE (gene-wise) correlation across lines - for each gene,
#           correlate mRNA vs protein across the shared lines. Summarised over
#           genes (median/mean), THIS is the metric the published benchmarks
#           report (CPTAC, CCLE, ProCan, Jarnuczak). It is the headline number.
#       (b) PER-LINE (sample-wise) correlation across genes - for each line,
#           correlate mRNA vs protein across shared genes. Typically HIGHER
#           because it is dominated by the large dynamic range of expressed vs
#           silent genes; it is NOT directly comparable to the gene-wise
#           benchmarks and is reported as a complementary sample-QC read.
#   - Spearman is the primary statistic (rank-based; robust to the TPM<->TMT
#     scale/normalisation mismatch and invariant to the log transform). Pearson
#     on the log scale is reported alongside.
#   - RNA gene_id (Ensembl) -> symbol collapse by SUMMING TPM (same as script
#     04); protein is already gene-symbol-level log2 relative abundance with NAs
#     preserved -> pairwise-complete handling throughout.
#   - Cross-line correlation is partly driven by subtype structure (as in CPTAC,
#     where it is driven by tumour heterogeneity); this is expected, not a
#     confounder to remove.
#   - QC GATE: if symbol mapping loses most proteins, or the gene-wise median is
#     wildly off the 0.4-0.5 benchmark, STOP (likely a mapping bug) rather than
#     report a spurious number.
# =============================================================================

source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(matrixStats); library(patchwork)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
set.seed(SEED)

# Local lab theme (visualization skill standard; setup does not define one) -----
theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text = element_text(colour = "black"),
      plot.title = element_text(size = rel(1.15), hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0, colour = "grey30",
                                   margin = margin(b = 6)),
      axis.title = element_text(size = rel(1.0)),
      axis.text = element_text(size = rel(0.85), colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      legend.title = element_text(size = rel(0.85)),
      legend.text = element_text(size = rel(0.8)),
      legend.background = element_blank(), legend.key = element_blank(),
      panel.background = element_blank(), panel.border = element_blank(),
      panel.grid = element_blank(), strip.background = element_blank(),
      strip.text = element_text(size = rel(0.95), face = "bold"),
      plot.margin = margin(8, 8, 8, 8))
}
# Published benchmarks (lit review Theme 3; gene-wise Spearman unless noted) -----
BENCH <- tribble(
  ~label, ~value, ~note,
  "CPTAC ovarian (median)", 0.45, "Zhang et al. 2016; different cohort and assay; contextual only",
  "CPTAC ovarian (mean)", 0.38, "Zhang et al. 2016; mean across genes; contextual only")
# No universal reproducibility ceiling is estimated for this dataset. External
# correlations cannot act as a pass/fail accuracy standard across assays/cohorts.


# -----------------------------------------------------------------------------
# 1. RNA: TPM -> symbol-collapsed log2 matrix (same collapse as script 04)
# -----------------------------------------------------------------------------
tpm <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"),
                       show_col_types = FALSE)
gene_meta <- t2g %>% distinct(ensembl_gene_id, external_gene_name) %>%
  filter(!is.na(external_gene_name), external_gene_name != "")

m_id <- as.matrix(tpm[, -1]); rownames(m_id) <- tpm$gene_id
id2sym  <- setNames(gene_meta$external_gene_name, gene_meta$ensembl_gene_id)
common  <- intersect(rownames(m_id), names(id2sym))
tpm_sym <- rowsum(m_id[common, , drop = FALSE], group = id2sym[common])  # symbols x RNA lines
rna_log <- log2(tpm_sym + 1)
message(sprintf("RNA symbol-level matrix: %d symbols x %d lines",
                nrow(rna_log), ncol(rna_log)))

# -----------------------------------------------------------------------------
# 2. Protein: load (gene symbols; log2 relative abundance; NAs preserved)
# -----------------------------------------------------------------------------
# [integration revision] Shared loader (00_setup.R). Script 05 drops the no-symbol
# rows at source and names non-representative duplicate-symbol rows SYMBOL|UNIPROT,
# so the row count is the single denominator in prot_feature_accounting.csv. The 70
# zero-plex rows are RETAINED here on purpose: they are NA in all 31 lines, so they
# contribute 0 to n_prot_per_gene and are removed by the n>=10 requirement below,
# which keeps the shared-symbol denominator comparable to the deposited matrix.
prot_mat <- read_prot_matrix()
prot_zero_plex <- attr(prot_mat, "zero_plex")

# subtype annotation (source of truth) for colouring / featured lines
ann <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  filter(provenance == "generated", analysis_include == "Y") %>%
  transmute(cell_line,
            subtype = factor(subtype, levels = c("HGS","CC","EC","MC","MMMT","SCCOHT")))
sub_cols <- setNames(RColorBrewer::brewer.pal(6, "Dark2"), levels(ann$subtype)) # match scripts 02/04

# -----------------------------------------------------------------------------
# 3. Restrict to shared genes + shared lines  (+ QC GATE)
# -----------------------------------------------------------------------------
shared_lines <- intersect(colnames(rna_log), colnames(prot_mat))
shared_genes <- intersect(rownames(rna_log), rownames(prot_mat))
frac_mapped  <- length(shared_genes) / nrow(prot_mat)
cat(sprintf("\nShared lines: %d | shared genes: %d (%.1f%% of %d proteins map to an RNA symbol)\n",
            length(shared_lines), length(shared_genes), 100 * frac_mapped, nrow(prot_mat)))
# GATE: a correct HGNC<->Ensembl-symbol mapping should retain the large majority
# of proteins. Losing most of them would indicate a mapping bug -> STOP.
stopifnot("QC GATE FAILED: <70% of proteins mapped to an RNA symbol - suspect a mapping bug" =
            frac_mapped > 0.70)

R <- rna_log[shared_genes, shared_lines]      # RNA log2(TPM+1)
P <- prot_mat[shared_genes, shared_lines]     # protein log2 relative abundance
stopifnot(identical(dim(R), dim(P)), identical(dimnames(R), dimnames(P)))

# -----------------------------------------------------------------------------
# 4. PER-LINE (sample-wise) correlation across shared genes
# -----------------------------------------------------------------------------
per_line <- tibble(cell_line = shared_lines) %>%
  left_join(ann, by = "cell_line") %>%
  rowwise() %>%
  mutate(
    n_genes  = sum(is.finite(R[, cell_line]) & is.finite(P[, cell_line])),
    spearman = cor(R[, cell_line], P[, cell_line], method = "spearman",
                   use = "pairwise.complete.obs"),
    pearson  = cor(R[, cell_line], P[, cell_line], method = "pearson",
                   use = "pairwise.complete.obs")) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 5. PER-GENE (gene-wise) correlation across lines  [benchmark-comparable]
#    Require protein measured in >= MIN_LINES lines for a stable estimate.
# -----------------------------------------------------------------------------
MIN_LINES <- 10L
n_prot_per_gene <- rowSums(!is.na(P))
keep <- names(which(n_prot_per_gene >= MIN_LINES))
cat(sprintf("Per-gene set: %d / %d shared genes have protein in >= %d lines\n",
            length(keep), length(shared_genes), MIN_LINES))

gene_cor <- function(method) {
  vapply(keep, function(g) {
    x <- R[g, ]; y <- P[g, ]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < MIN_LINES || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = method))
  }, numeric(1))
}
per_gene <- tibble(gene = keep,
                   n_lines  = n_prot_per_gene[keep],
                   spearman = gene_cor("spearman"),
                   pearson  = gene_cor("pearson")) %>%
  filter(!is.na(spearman))
cat(sprintf("Per-gene correlations computed for %d genes (dropped %d zero-variance/NA)\n",
            nrow(per_gene), length(keep) - nrow(per_gene)))

# Sensitivity of the headline median to the MIN_LINES threshold ---------------
sens <- map_dfr(c(5L, 10L, 15L, 20L, 30L), function(m) {
  gk <- names(which(n_prot_per_gene >= m))
  vals <- vapply(gk, function(g) {
    x <- R[g, ]; y <- P[g, ]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < m || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  }, numeric(1))
  tibble(min_lines = m, n_genes = sum(!is.na(vals)),
         median_spearman = median(vals, na.rm = TRUE))
})

# -----------------------------------------------------------------------------
# 6. Distribution summaries (median, IQR, range) + benchmark comparison
# -----------------------------------------------------------------------------
summ <- function(x) {
  x <- x[is.finite(x)]
  tibble(n = length(x), median = median(x), mean = mean(x),
         q25 = quantile(x, .25), q75 = quantile(x, .75),
         min = min(x), max = max(x))
}
# -----------------------------------------------------------------------------
# 5b. DISCLOSE THE n STRUCTURE. The per-gene n is NOT 30 for every gene: TMT
#     block-missingness leaves it between MIN_LINES and 30, and a Spearman on
#     n=10 is nearly uninformative. Pooling those into one median without saying
#     so is not defensible, so the whole n distribution is written out and the
#     COMPLETE-CASE median is reported as the headline-comparable statistic
#     (same genes, same 30 lines, no n heterogeneity). Also record the two counts
#     that were previously conflated: shared symbols vs genes that actually
#     receive a correlation.
# -----------------------------------------------------------------------------
n_complete_lines <- length(shared_lines)
cc_gene   <- per_gene %>% filter(n_lines == n_complete_lines)
frac_neg  <- mean(per_gene$spearman < 0)
n_dist <- per_gene %>% dplyr::count(n_lines, name = "n_genes") %>% arrange(n_lines) %>%
  mutate(cum_genes = cumsum(n_genes),
         pct_genes = round(100 * n_genes / nrow(per_gene), 2))
n_thresh <- tibble(threshold = c(10L, 15L, 20L, 25L, n_complete_lines)) %>%
  rowwise() %>%
  mutate(n_genes_below = sum(per_gene$n_lines <  threshold),
         n_genes_atleast = sum(per_gene$n_lines >= threshold),
         median_spearman_atleast = median(per_gene$spearman[per_gene$n_lines >= threshold])) %>%
  ungroup()
readr::write_csv(n_dist,   file.path(OUT, "integ_rnaprot_n_distribution.csv"))
readr::write_csv(n_thresh, file.path(OUT, "integ_rnaprot_n_thresholds.csv"))

cat(sprintf("\n=== Per-gene n structure (n = lines with BOTH assays for that gene) ===\n"))
cat(sprintf("range %d-%d | genes at n<15: %d | n<20: %d | complete case (n=%d): %d\n",
            min(per_gene$n_lines), max(per_gene$n_lines),
            sum(per_gene$n_lines < 15), sum(per_gene$n_lines < 20),
            n_complete_lines, nrow(cc_gene)))
cat(sprintf("Genes with NEGATIVE RNA-protein correlation: %d (%.1f%%)\n",
            sum(per_gene$spearman < 0), 100 * frac_neg))
print(as.data.frame(n_thresh), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5c. UNCERTAINTY ON THE NEGATIVE GENES  [audit revision 2026-07-27]
#     A negative point estimate is not evidence of anti-correlation at these n.
#     With 10-30 paired models the sampling SE of a Spearman rho is roughly
#     1/sqrt(n-3), i.e. ~0.38 at n=10 and ~0.19 at n=30, so a substantial share of
#     the negative estimates are indistinguishable from zero. Earlier drafts wrote
#     that for these genes "a transcript measurement is not a proxy for protein in
#     any direction", which asserts a real negative relationship the estimates do
#     not establish. So: attach a two-sided asymptotic p to every gene, BH-adjust
#     across the whole per-gene set, and report the negative genes at three tiers
#     (any negative point estimate / nominally significant / BH-significant). The
#     deposited per-gene file carries p and q so a reuser can set their own filter.
# TEST CHOICE  [audit revision 2026-07-27, second round]
# Use the CONVENTIONAL asymptotic Spearman test - the t approximation on rho with
# n - 2 degrees of freedom, which is what stats::cor.test(method = "spearman") falls
# back to once n is past the exact range or ties are present (and ties are ubiquitous
# here). An earlier version used a Fisher z transform scaled by 1/sqrt(n-3); that is
# a Pearson-derived approximation applied to a rank correlation, and although it is
# in wide use it is the less standard choice. The two agree closely on this data
# (Fisher z gave 51 nominal / 30 BH-significant negatives; the t approximation gives
# the counts printed below), so nothing substantive turns on it - but the
# conventional test is the one a reader can reproduce with cor.test().
per_gene <- per_gene %>%
  mutate(
    t_stat = spearman * sqrt((n_lines - 2) / pmax(1 - spearman^2, 0)),
    p_value = ifelse(is.na(t_stat), NA_real_,
                     2 * stats::pt(-abs(t_stat), df = n_lines - 2)),
    q_value = stats::p.adjust(p_value, method = "BH"))
# TERMINOLOGY  [second round]. An inverse association is NOT an absence of proxy
# value: a gene whose protein reliably falls as its transcript rises is predictable,
# simply in the opposite direction. Proxy value additionally depends on predictive
# error and external validation, neither of which is assessed here. So these tiers
# are named for the ASSOCIATION they establish, and the interpretation field says
# what the counts do and do not license.
neg_tiers <- tibble(
  tier = c("inverse point estimate (rho < 0)",
           "inverse and nominally significant (p < 0.05)",
           "FDR-supported inverse association (q < 0.05)",
           "FDR-supported positive association (q < 0.05)"),
  n_genes = c(sum(per_gene$spearman < 0),
              sum(per_gene$spearman < 0 & per_gene$p_value < 0.05, na.rm = TRUE),
              sum(per_gene$spearman < 0 & per_gene$q_value < 0.05, na.rm = TRUE),
              sum(per_gene$spearman > 0 & per_gene$q_value < 0.05, na.rm = TRUE))) %>%
  mutate(pct_of_pergene_set = round(100 * n_genes / nrow(per_gene), 2),
         n_pergene_set = nrow(per_gene),
         most_negative_rho = min(per_gene$spearman),
         median_n_lines_negative = median(per_gene$n_lines[per_gene$spearman < 0]),
         test = "two-sided asymptotic Spearman t approximation, df = n - 2; BH across the per-gene set",
         analysis_unit = "model; exploratory because related sublines are not independent",
         interpretation = paste("Model-level q-values are exploratory; patient-level tests are deposited separately.",
                                "An inverse POINT ESTIMATE is a candidate for discordant",
                                "regulation; an association does not establish regulation or predictive validity.",
                                "ASSOCIATION, which is not the same as transcript having no",
                                "proxy value - an inverse relationship can still be",
                                "predictive. Predictive error and external validation are",
                                "not assessed here."))
readr::write_csv(neg_tiers, file.path(OUT, "integ_rnaprot_negative_genes.csv"))
cat("\n=== Inverse per-gene correlations, with uncertainty ===\n")
print(as.data.frame(neg_tiers %>% select(tier, n_genes, pct_of_pergene_set)), row.names = FALSE)
cat(sprintf("Most negative rho %.3f; median n for the inverse-estimate genes %.0f models.\n",
            min(per_gene$spearman), median(per_gene$n_lines[per_gene$spearman < 0])))
cat("These model-level tests are exploratory; use the patient-level tests written below.\n")
cat("an FDR-supported INVERSE ASSOCIATION, not 'transcript is not a proxy' - an inverse\n")
cat("relationship can still be predictive, in the other direction.\n")

# -----------------------------------------------------------------------------
# 5d. PATIENT-REPRESENTATIVE SENSITIVITY  [audit revision 2026-07-27]
#     Three of the 30 dual-layer models are same-donor sublines, so the per-model
#     and per-gene medians are computed on partially non-independent columns. This
#     recomputes both on one model per patient and reports the shift, so no reader
#     has to assume it is negligible. It is a sensitivity analysis, not a
#     replacement: the 30-model values remain the resource-level description.
FAM_PATH <- file.path(PROJ, "metadata", "line_family_map.csv")
stopifnot("metadata/line_family_map.csv must exist — run 15_patient_family_map.R first (see run_all.sh)" =
            file.exists(FAM_PATH))
fam_c <- readr::read_csv(FAM_PATH, show_col_types = FALSE)
stopifnot("line_family_map.csv must carry cell_line and patient_representative" =
            all(c("cell_line", "patient_representative") %in% names(fam_c)))
rep_lines_c <- intersect(fam_c$cell_line[fam_c$patient_representative], shared_lines)
stopifnot("collapsing donors must not increase the model count" =
            length(rep_lines_c) <= length(shared_lines))
cat(sprintf("\nPatient-representative dual-layer subset: %d of %d models (%d donors collapsed)\n",
            length(rep_lines_c), length(shared_lines),
            length(shared_lines) - length(rep_lines_c)))
Rr <- R[, rep_lines_c, drop = FALSE]; Pr <- P[, rep_lines_c, drop = FALSE]
per_line_rep <- vapply(rep_lines_c, function(l) {
  ok <- is.finite(Rr[, l]) & is.finite(Pr[, l])
  suppressWarnings(cor(Rr[ok, l], Pr[ok, l], method = "spearman"))
}, numeric(1))
n_prot_rep <- rowSums(!is.na(Pr))
keep_rep <- names(which(n_prot_rep >= MIN_LINES))
rep_cor <- function(g, method) {
  x <- Rr[g, ]; y <- Pr[g, ]; ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < MIN_LINES || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}
per_gene_rep_tbl <- tibble(
  gene = keep_rep, n_patients = n_prot_rep[keep_rep],
  spearman = vapply(keep_rep, rep_cor, numeric(1), method = "spearman"),
  pearson = vapply(keep_rep, rep_cor, numeric(1), method = "pearson")) %>%
  filter(is.finite(spearman)) %>%
  mutate(t_stat = spearman * sqrt((n_patients - 2) / pmax(1 - spearman^2, 0)),
         p_value = 2 * stats::pt(-abs(t_stat), df = n_patients - 2),
         q_value = p.adjust(p_value, "BH"),
         analysis_unit = "one model per patient",
         interpretation = "exploratory cross-model association; source/histotype/plex structure not removed; does not demonstrate regulation or predictive validity")
readr::write_csv(per_gene_rep_tbl, file.path(OUT, "integ_rnaprot_patientrep_cor.csv"))
per_gene_rep <- setNames(per_gene_rep_tbl$spearman, per_gene_rep_tbl$gene)
# Quantify inverse associations using INDEPENDENT PATIENT representatives, not
# q-values computed from all related sublines. The all-model correlations remain
# descriptive resource summaries in integ_rnaprot_cor.csv.
negative_summary <- function(tab, unit) {
  tibble(tier = c("inverse point estimate (rho < 0)", "inverse and nominally significant (p < 0.05)",
                  "inverse association with BH q < 0.05", "positive association with BH q < 0.05"),
         n_genes = c(sum(tab$spearman < 0), sum(tab$spearman < 0 & tab$p_value < 0.05),
                     sum(tab$spearman < 0 & tab$q_value < 0.05),
                     sum(tab$spearman > 0 & tab$q_value < 0.05)),
         analysis_unit = unit, n_pergene_set = nrow(tab),
         pct_of_pergene_set = 100 * n_genes / nrow(tab),
         test = "two-sided asymptotic Spearman t; BH across all eligible genes within analysis unit",
         interpretation = "patient-representative q-values are primary for association; all-model q-values exploratory; association does not imply regulation or predictive validity")
}
readr::write_csv(bind_rows(
  negative_summary(per_gene, "30 models (descriptive; related sublines)"),
  negative_summary(per_gene_rep_tbl, "27 patient representatives (primary association unit)")),
  file.path(OUT, "integ_rnaprot_negative_genes.csv"))

conc_sens <- tibble(
  statistic = c("per-model Spearman (median)", "per-gene Spearman (median)",
                "fraction of per-gene estimates negative"),
  all_models = c(median(per_line$spearman), median(per_gene$spearman), frac_neg),
  patient_representatives = c(median(per_line_rep), median(per_gene_rep),
                              mean(per_gene_rep < 0)),
  n_all = c(length(shared_lines), nrow(per_gene), nrow(per_gene)),
  n_reps = c(length(rep_lines_c), length(per_gene_rep), length(per_gene_rep))) %>%
  mutate(shift = patient_representatives - all_models,
         note = paste("MIN_LINES rule held at", MIN_LINES,
                      "in both columns, so the representative gene set is smaller"))
readr::write_csv(conc_sens, file.path(OUT, "integ_rnaprot_patientrep_sensitivity.csv"))
cat("\n=== Concordance: same-donor subline sensitivity ===\n")
print(as.data.frame(conc_sens %>% select(statistic, all_models, patient_representatives,
                                         shift, n_all, n_reps) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)

dist_summary <- bind_rows(
  summ(per_line$spearman) %>% mutate(metric = "per_line_spearman", .before = 1),
  summ(per_line$pearson)  %>% mutate(metric = "per_line_pearson",  .before = 1),
  summ(per_gene$spearman) %>% mutate(metric = "per_gene_spearman", .before = 1),
  summ(per_gene$pearson)  %>% mutate(metric = "per_gene_pearson",  .before = 1),
  summ(cc_gene$spearman)  %>% mutate(metric = "per_gene_spearman_complete_case", .before = 1),
  summ(cc_gene$pearson)   %>% mutate(metric = "per_gene_pearson_complete_case",  .before = 1)) %>%
  # the counts the manuscript conflated, carried in the same file as the medians
  mutate(n_shared_symbols   = length(shared_genes),
         n_protein_rows     = nrow(prot_mat),
         n_pergene_reported = nrow(per_gene),
         min_lines_rule     = MIN_LINES,
         frac_pergene_negative = round(frac_neg, 4))

cat("\n=== Concordance distribution summary ===\n")
print(as.data.frame(dist_summary %>% mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)
cat(sprintf("\nCOUNTS, stated separately because they are not the same number:\n"))
cat(sprintf("  shared gene symbols (RNA symbol present in the protein matrix): %d\n",
            length(shared_genes)))
cat(sprintf("  genes that RECEIVE a correlation (protein in >= %d lines, non-zero variance): %d\n",
            MIN_LINES, nrow(per_gene)))
cat(sprintf("  complete-case genes (protein in all %d shared lines): %d\n",
            n_complete_lines, nrow(cc_gene)))
cat(sprintf("The n attached to the reported median MUST be %d, not %d.\n",
            nrow(per_gene), length(shared_genes)))

gw_med <- dist_summary$median[dist_summary$metric == "per_gene_spearman"]
gw_mn  <- dist_summary$mean[dist_summary$metric   == "per_gene_spearman"]
cc_med <- dist_summary$median[dist_summary$metric == "per_gene_spearman_complete_case"]
cat(sprintf("\nHEADLINE (gene-wise Spearman): median %.3f / mean %.3f across %d genes (n per gene %d-%d)\n",
            gw_med, gw_mn, nrow(per_gene), min(per_gene$n_lines), max(per_gene$n_lines)))
cat(sprintf("HEADLINE-COMPARABLE (complete case, n=%d lines for every gene): median %.3f across %d genes\n",
            n_complete_lines, cc_med, nrow(cc_gene)))
cat("Benchmarks (gene-wise Spearman):\n")
print(as.data.frame(BENCH), row.names = FALSE)
cat(sprintf("\nObserved gene-wise median %.2f; CPTAC ovarian median 0.45 is contextual, not an accuracy threshold.\n", gw_med))
# GATE: sanity-check the headline is in a plausible moderate range.
stopifnot("Non-finite descriptive correlation summary" = is.finite(gw_med))

# -----------------------------------------------------------------------------
# 7. Output table (tidy: both levels)  -> output/integ_rnaprot_cor.csv
# -----------------------------------------------------------------------------
out_tbl <- bind_rows(
  per_line %>% transmute(level = "per_line", id = cell_line,
                         subtype = as.character(subtype),
                         n = n_genes, spearman, pearson),
  per_gene %>% transmute(level = "per_gene", id = gene,
                         subtype = NA_character_,
                         n = n_lines, spearman, pearson,
                         # deposited so a reuser can filter on reliability rather
                         # than on the sign of a point estimate (section 5c)
                         p_value, q_value,
                         inference_note = "model-level asymptotic tests are exploratory; use integ_rnaprot_patientrep_cor.csv for independent-patient associations")) %>%
  mutate(across(c(spearman, pearson), ~round(.x, 4)),
         across(c(p_value, q_value), ~signif(.x, 4)))
readr::write_csv(out_tbl, file.path(OUT, "integ_rnaprot_cor.csv"))
readr::write_csv(dist_summary %>% mutate(across(where(is.numeric), ~round(.x, 4))),
                 file.path(OUT, "integ_rnaprot_cor_summary.csv"))

# -----------------------------------------------------------------------------
# 8. Figures
# -----------------------------------------------------------------------------
# --- 8a. Distributions: per-gene density (with benchmarks) + per-line strip ---
bench_gw <- BENCH %>% filter(label != "CPTAC ovarian (mean)")  # plot contextual median only
pA <- ggplot(per_gene, aes(x = spearman)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.05,
                 fill = "grey80", colour = "white", boundary = 0) +
  geom_density(colour = COOK_NAVY, linewidth = 0.7) +
  geom_vline(xintercept = gw_med, colour = COOK_RUST, linewidth = 0.9) +
  geom_vline(data = bench_gw, aes(xintercept = value),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = bench_gw, aes(x = value, y = 0, label = sub(" \\(median\\)| \\(integrated\\)| integrated", "", label)),
    angle = 90, size = 2.5, colour = "grey30", direction = "x",
    segment.colour = NA, hjust = 0, nudge_y = 0.15) +
  annotate("text", x = gw_med, y = Inf, vjust = 1.6, hjust = -0.08,
           label = sprintf("ours = %.2f", gw_med), colour = COOK_RUST, size = 3) +
  labs(title = "A  Gene-wise concordance (external reference for context)",
       subtitle = sprintf("Spearman r per gene across %d matched lines; n = %d genes",
                          length(shared_lines), nrow(per_gene)),
       x = "Spearman r (mRNA vs protein, across lines)", y = "Density") +
  coord_cartesian(xlim = c(-0.6, 1)) + theme_lab()

pB <- per_line %>%
  select(cell_line, subtype, spearman, pearson) %>%
  pivot_longer(c(spearman, pearson), names_to = "stat", values_to = "r") %>%
  mutate(stat = factor(stat, c("spearman","pearson"))) %>%
  ggplot(aes(x = stat, y = r)) +
  geom_jitter(aes(colour = subtype), width = 0.15, height = 0, size = 2, alpha = 0.9) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5,
               colour = "black", linewidth = 0.4) +
  scale_colour_manual(values = sub_cols, name = "subtype", drop = FALSE) +
  labs(title = "B  Sample-wise concordance (per line)",
       subtitle = sprintf("r across %d shared genes; each point = 1 line (n = %d)",
                          length(shared_genes), length(shared_lines)),
       x = NULL, y = "correlation (mRNA vs protein, across genes)") +
  coord_cartesian(ylim = c(0, 1)) + theme_lab()

ggsave(file.path(FIGS, "12_concordance_distributions.pdf"),
       pA + pB + patchwork::plot_layout(widths = c(1.35, 1)),
       width = 11, height = 4.6)

# --- 8b. Representative scatters -------------------------------------------
# Gene-vs-gene (across lines): pick a well-measured gene nearest the gene-wise median.
full_genes <- per_gene %>% filter(n_lines == length(shared_lines))
rep_gene <- full_genes$gene[which.min(abs(full_genes$spearman - gw_med))]
rg_r <- per_gene$spearman[per_gene$gene == rep_gene]
gdf <- tibble(cell_line = shared_lines, rna = R[rep_gene, ], prot = P[rep_gene, ]) %>%
  left_join(ann, by = "cell_line")
pG <- ggplot(gdf, aes(rna, prot)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey60", linewidth = 0.5) +
  geom_point(aes(colour = subtype), size = 2.4, alpha = 0.9) +
  scale_colour_manual(values = sub_cols, name = "subtype", drop = FALSE) +
  labs(title = sprintf("A  Representative gene: %s", rep_gene),
       subtitle = sprintf("across %d lines; Spearman r = %.2f (~ gene-wise median)",
                          length(shared_lines), rg_r),
       x = expression(mRNA~log[2]*"(TPM+1)"),
       y = expression(protein~log[2]~relative~abundance)) + theme_lab()

# Line-vs-line (across genes): line nearest the per-line median Spearman.
pl_med <- median(per_line$spearman)
rep_line <- per_line$cell_line[which.min(abs(per_line$spearman - pl_med))]
rl_r <- per_line$spearman[per_line$cell_line == rep_line]
ldf <- tibble(rna = R[, rep_line], prot = P[, rep_line]) %>% filter(is.finite(rna), is.finite(prot))
pL <- ggplot(ldf, aes(rna, prot)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "genes") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey85", linewidth = 0.5) +
  labs(title = sprintf("B  Representative line: %s", rep_line),
       subtitle = sprintf("across %d shared genes; Spearman r = %.2f (~ per-line median)",
                          nrow(ldf), rl_r),
       x = expression(mRNA~log[2]*"(TPM+1)"),
       y = expression(protein~log[2]~relative~abundance)) + theme_lab()

ggsave(file.path(FIGS, "12_concordance_scatter_examples.pdf"),
       pG + pL + patchwork::plot_layout(widths = c(1, 1)),
       width = 10.5, height = 4.6)

cat("\n--- Sensitivity of gene-wise median Spearman to MIN_LINES threshold ---\n")
print(as.data.frame(sens %>% mutate(median_spearman = round(median_spearman, 3))),
      row.names = FALSE)
cat(sprintf("\nRepresentative gene: %s (r=%.2f); representative line: %s (r=%.2f)\n",
            rep_gene, rg_r, rep_line, rl_r))
cat("\nOutputs: output/integ_rnaprot_cor.csv, integ_rnaprot_cor_summary.csv,\n")
cat("         integ_rnaprot_n_distribution.csv, integ_rnaprot_n_thresholds.csv\n")
cat("Figures: figs/12_concordance_distributions.pdf, figs/12_concordance_scatter_examples.pdf\n")
cat("NOTE: this analysis is PER LINE (30 shared lines), not per patient.\n")

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
write_session_info("12_rna_protein_concordance")

message("\n12_rna_protein_concordance.R complete.")
