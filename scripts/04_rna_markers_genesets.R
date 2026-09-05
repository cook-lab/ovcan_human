# =============================================================================
# Script: 04_rna_markers_genesets.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Canonical subtype-marker recovery (primary validation), the SWI/SNF
#          check (feeds authentication), and Hallmark gene-set scoring
#          (singscore). Together with the GO recovery in script 03, this is the
#          "recapitulates known subtype biology" evidence for the descriptor.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   1 (RNA re-processing) — step 4 of 4
# =============================================================================
# Assay-aware notes:
#   - Marker "recovery" is judged on log2-TPM: z-scored across lines shows
#     ENRICHMENT (relative), but a marker only truly "lands" if it is also
#     ABSOLUTELY expressed in its subtype. We report both (heatmap = z; summary
#     table = mean log2-TPM in-subtype vs rest + rank), and annotate n.
#   - Cultured 2D lines lose some in-vivo markers and gain culture-adaptation
#     genes; VIM (mesenchymal) in particular is broadly expressed in culture, so
#     EMT markers are weaker discriminators here than in tissue. Stated, not hidden.
#   - SCCOHT is defined by LOSS (SMARCA4/SMARCA2), not a positive marker; RNA
#     loss is only partly informative because SCCOHT SMARCA4 inactivation is
#     often post-transcriptional (mutation/protein loss with retained mRNA),
#     whereas SMARCA2 is epigenetically silenced (mRNA low). Protein/WES arms are
#     the definitive test; RNA here is supportive-if-present.
# =============================================================================

source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(matrixStats); library(singscore)
  library(ComplexHeatmap); library(circlize); library(fgsea); library(GSEABase)
})
# ComplexHeatmap/GSEABase/AnnotationDbi mask several dplyr verbs — reassert them
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count  <- dplyr::count;  first  <- dplyr::first;  slice  <- dplyr::slice
set.seed(SEED)

# -----------------------------------------------------------------------------
# 1. TPM -> symbol-collapsed log2 matrix + ordered sample annotation
# -----------------------------------------------------------------------------
tpm <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"),
                       show_col_types = FALSE)
gene_meta <- t2g %>% distinct(ensembl_gene_id, external_gene_name) %>%
  filter(!is.na(external_gene_name), external_gene_name != "")

m_id <- as.matrix(tpm[, -1]); rownames(m_id) <- tpm$gene_id

# collapse gene_id -> symbol by SUMMING TPM (patch/scaffold duplicates collapse;
# TPM is additive so summing recovers the symbol-level abundance)
id2sym <- setNames(gene_meta$external_gene_name, gene_meta$ensembl_gene_id)
common <- intersect(rownames(m_id), names(id2sym))
sym_vec <- id2sym[common]
tpm_sym <- rowsum(m_id[common, , drop = FALSE], group = sym_vec)   # symbols x lines
logtpm  <- log2(tpm_sym + 1)
message(sprintf("Symbol-level matrix: %d symbols x %d lines", nrow(logtpm), ncol(logtpm)))

ann <- readr::read_csv(file.path(OUT, "rna_sample_annotation.csv"),
                       show_col_types = FALSE) %>%
  transmute(cell_line, subtype = factor(subtype, levels = c("HGS","CC","EC","MC","MMMT","SCCOHT")),
            site) %>%
  arrange(subtype, cell_line)
stopifnot(setequal(ann$cell_line, colnames(logtpm)))
logtpm <- logtpm[, ann$cell_line]                     # order cols by subtype
subtypes <- levels(ann$subtype)
n_by <- table(ann$subtype)

sub_cols  <- setNames(RColorBrewer::brewer.pal(6, "Dark2"), subtypes)   # match scripts 01/02
z_fun     <- colorRamp2(c(-2, 0, 2), c(COOK_NAVY, "white", COOK_RUST))  # CB-safe diverging
zscore    <- function(mat) t(scale(t(mat)))            # z across lines, per gene
sub_lab   <- setNames(sprintf("%s\n(n=%d)", subtypes, as.integer(n_by)), subtypes)

# -----------------------------------------------------------------------------
# 2. Canonical marker panel (curated from the literature)
# -----------------------------------------------------------------------------
markers <- tribble(
  ~marks,            ~symbol,   ~direction, ~note,
  "HGS",             "PAX8",    "up",   "Mullerian/serous TF",
  "HGS",             "WT1",     "up",   "serous marker",
  "HGS",             "MUC16",   "up",   "CA125",
  "HGS",             "MECOM",   "up",   "3q26 HGSC amplicon TF",
  "HGS",             "SOX17",   "up",   "serous TF",
  "CC (clear cell)", "HNF1B",   "up",   "clear-cell master TF",
  "CC (clear cell)", "NAPSA",   "up",   "napsin A",
  "CC (clear cell)", "SPP1",    "up",   "osteopontin",
  "CC (clear cell)", "GPX3",    "up",   "glutathione peroxidase",
  "CC (clear cell)", "GCLC",    "up",   "glutathione synthesis",
  "EC (endometrioid)", "ESR1",  "up",   "oestrogen receptor (Hollis 2020 EC panel)",
  "EC (endometrioid)", "PGR",   "up",   "progesterone receptor (Hollis 2020 EC panel)",
  "EC (endometrioid)", "ARID1A","down", "SWI/SNF; WEAK RNA proxy - EC loss is mutational, not transcriptional",
  "MC (mucinous)",   "CDX2",    "up",   "intestinal TF",
  "MC (mucinous)",   "MUC2",    "up",   "intestinal mucin",
  "MC (mucinous)",   "TFF1",    "up",   "trefoil factor",
  "MC (mucinous)",   "TFF3",    "up",   "trefoil factor",
  "MC (mucinous)",   "KRT20",   "up",   "intestinal keratin",
  "MC (mucinous)",   "MUC5AC",  "up",   "gastric mucin",
  "MMMT (EMT)",      "VIM",     "up",   "vimentin (broad in culture)",
  "MMMT (EMT)",      "ZEB1",    "up",   "EMT TF",
  "MMMT (EMT)",      "CDH2",    "up",   "N-cadherin",
  "MMMT (EMT)",      "SNAI2",   "up",   "SLUG",
  "SCCOHT (loss)",   "SMARCA4", "down", "SWI/SNF (loss; often protein-level)",
  "SCCOHT (loss)",   "SMARCA2", "down", "SWI/SNF (epigenetically silenced)",
  "proliferation",   "MKI67",   "up",   "proliferation control"
) %>% mutate(marks = factor(marks, levels = c("HGS","CC (clear cell)","EC (endometrioid)",
                                              "MC (mucinous)","MMMT (EMT)","SCCOHT (loss)",
                                              "proliferation")))
stopifnot("marker(s) absent from matrix" = all(markers$symbol %in% rownames(logtpm)),
          "marker symbols must be unique (sm_mat is keyed on symbol)" =
            !any(duplicated(markers$symbol)))
# ENDOMETRIOID SET (added at revision). The panel previously covered HGS/CC/MC/
# MMMT/SCCOHT and lacked EC markers. Include the EC group while retaining
# the supplied annotation and showing individual-model results.
# Hollis et al. 2020 names ARID1A, ESR1, PGR and vimentin. Three are added here;
# VIM is deliberately NOT duplicated because it is already scored in the MMMT/EMT
# group (a symbol may appear once — sm_mat is keyed on symbol) and it is broadly
# expressed in 2D culture, so it is a poor discriminator either way. Its EC mean
# is reported in section 2b below for completeness. Read the EC result with the n
# in mind: the annotated EC group is n=2, and discordant SWI/SNF evidence for
# TOV112D warrants annotation review. Section 2b therefore shows both models
# separately; these marker data do not establish a histological reclassification.

mk_mat <- logtpm[markers$symbol, , drop = FALSE]
mk_z   <- zscore(mk_mat)

# --- marker heatmap: markers (rows, split by intended subtype) x lines (cols by subtype)
col_ha <- HeatmapAnnotation(subtype = ann$subtype,
                            col = list(subtype = sub_cols),
                            annotation_name_gp = grid::gpar(fontsize = 8))
pdf(file.path(FIGS, "04_rna_markers_heatmap.pdf"), width = 11, height = 7)
draw(Heatmap(mk_z, name = "z(log2 TPM)", col = z_fun,
             top_annotation = col_ha,
             row_split = markers$marks, cluster_rows = FALSE, cluster_row_slices = FALSE,
             row_title_rot = 0, row_gap = grid::unit(1.5, "mm"),
             column_split = ann$subtype, cluster_columns = FALSE, cluster_column_slices = FALSE,
             row_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
             column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
             row_names_gp = grid::gpar(fontsize = 8),
             column_names_gp = grid::gpar(fontsize = 7),
             column_title = "Canonical subtype markers x cell lines (z-scored log2 TPM)"),
     heatmap_legend_side = "right")
dev.off()

# --- quantitative "lands right" summary (direction-aware; absolute + relative)
# NB: the row-panel labels ("CC (clear cell)") differ from the subtype CODES in
# the annotation ("CC"); map explicitly. Many canonical markers are SHARED
# Mullerian/serous lineage markers (e.g. PAX8 is high in HGS AND CC), so a
# rank-1-of-6 test understates recovery. We grade an UP marker as landing when
# its intended subtype is in the TOP-2 of the 6 subtype means AND it is actually
# expressed there (mean log2 TPM > 1); a DOWN marker (SCCOHT SWI/SNF) lands when
# the intended subtype is in the BOTTOM-2. MKI67 is a proliferation control
# (no single subtype) — reported descriptively.
code_map <- c("HGS" = "HGS", "CC (clear cell)" = "CC", "EC (endometrioid)" = "EC",
              "MC (mucinous)" = "MC",
              "MMMT (EMT)" = "MMMT", "SCCOHT (loss)" = "SCCOHT",
              "proliferation" = NA_character_)
# subtype-mean matrix (markers x 6 subtypes) on log2 TPM
sm_mat <- t(vapply(markers$symbol, function(g)
  tapply(mk_mat[g, ], factor(ann$subtype, levels = subtypes), mean)[subtypes],
  numeric(length(subtypes))))
colnames(sm_mat) <- subtypes

marker_summary <- markers %>%
  mutate(marks_code = unname(code_map[as.character(marks)])) %>%
  bind_cols(as_tibble(round(sm_mat, 2)) %>% rename_with(~ paste0(.x, "_mean"))) %>%
  rowwise() %>%
  mutate(
    mean_in     = if (is.na(marks_code)) NA_real_ else sm_mat[symbol, marks_code],
    rank_in     = if (is.na(marks_code)) NA_integer_
                  else as.integer(rank(-sm_mat[symbol, ], ties.method = "min")[marks_code]),
    top_subtype = colnames(sm_mat)[which.max(sm_mat[symbol, ])],
    lands_right = dplyr::case_when(
      is.na(marks_code)  ~ NA,
      direction == "up"  ~ rank_in <= 2L & mean_in > 1,
      direction == "down"~ rank_in >= 5L)
  ) %>% ungroup() %>%
  mutate(mean_in = round(mean_in, 2)) %>%
  select(marks, marks_code, symbol, direction, note,
         ends_with("_mean"), mean_in, rank_in, top_subtype, lands_right)
readr::write_csv(marker_summary, file.path(OUT, "rna_markers_summary.csv"))

cat("\n=== Canonical marker recovery (per-subtype mean log2 TPM) ===\n")
print(as.data.frame(marker_summary %>%
        select(marks_code, symbol, direction, ends_with("_mean"),
               rank_in, top_subtype, lands_right)), row.names = FALSE)
scored <- marker_summary %>% filter(!is.na(lands_right))
cat(sprintf("\nMarkers landing in the expected subtype: %d / %d graded (top-2 & expressed for 'up'; bottom-2 for 'down'; MKI67 excluded as control)\n",
            sum(scored$lands_right), nrow(scored)))
cat("Not landing (worth noting):",
    paste(scored$symbol[!scored$lands_right], collapse = ", "), "\n")
# Marker genes are correlated and histotype sizes are unequal. An independent
# Bernoulli(1/3) model is not a calibrated null. Use the joint permutation below.
bt <- function(k, n) NA_real_ # retired binomial field retained for schema compatibility

# -----------------------------------------------------------------------------
# 2a. JOINT-PERMUTATION NULL for marker recovery  [audit revision 2026-07-27]
# -----------------------------------------------------------------------------
# The retired binomial test treated the 25 graded markers as independent Bernoulli
# trials at p = 1/3. They are NOT independent: markers intended for the same
# histotype are co-regulated (PAX8/WT1/MUC16 all track a serous programme), and
# a single relabelling of the models can move several markers together. The
# binomial p is therefore anti-conservative by an unknown amount.
#
# The correct null permutes the HISTOTYPE LABELS and recomputes the whole landing
# count, which preserves (a) the marker-marker correlation structure, (b) the
# group sizes, and (c) the absolute-expression floor in the 'up' rule. The test
# statistic is the same landing count, so only the reference distribution changes.
#
# Two units are reported because they answer different questions (see the
# patient-representative note in 21_rna_sensitivity.R):
#   - 31 RNA models      : three patient pairs contribute related sublines, so
#                          the columns are not independent patient observations
#   - 28 patient reps    : the INFERENTIAL unit; one RNA model per patient, so
#                          source-restricted exchangeability is assessed below
# The permutation p is (1 + #{null >= observed}) / (n_perm + 1), so the floor at
# 20,000 draws is 5.0e-5.
N_PERM_MARKER <- 20000L
scored_idx <- which(!is.na(marker_summary$marks_code))
land_count <- function(mat, labels) {
  labels <- factor(labels, levels = subtypes)
  means <- vapply(subtypes,
                  function(st) rowMeans(mat[, labels == st, drop = FALSE]),
                  numeric(nrow(mat)))
  rownames(means) <- rownames(mat)
  sum(vapply(scored_idx, function(i) {
    sym <- marker_summary$symbol[i]
    vals <- means[sym, ]
    r <- rank(-vals, ties.method = "min")[marker_summary$marks_code[i]]
    if (marker_summary$direction[i] == "up") r <= 2L && vals[marker_summary$marks_code[i]] > 1
    else r >= 5L
  }, logical(1)))
}
perm_null <- function(mat, labels, n_perm = N_PERM_MARKER, strata = NULL) {
  obs <- land_count(mat, labels)
  set.seed(SEED)                                   # deterministic, seed = 1234
  null <- replicate(n_perm, {
    permuted <- labels
    if (is.null(strata)) permuted <- labels[sample.int(length(labels))]
    else for (g in unique(strata)) {
      ii <- which(strata == g); permuted[ii] <- labels[ii[sample.int(length(ii))]]
    }
    land_count(mat, permuted)
  })
  tibble(observed = obs, n_markers = length(scored_idx),
         null_mean = mean(null), null_median = median(null),
         null_p95 = unname(quantile(null, 0.95)), null_max = max(null),
         n_perm = n_perm,
         permutation_p = (1 + sum(null >= obs)) / (n_perm + 1))
}

fam_mk <- readr::read_csv(file.path(PROJ, "metadata", "line_family_map.csv"),
                          show_col_types = FALSE)
stopifnot("line_family_map.csv must carry patient_representative and has_rna" =
            all(c("cell_line", "patient_representative", "has_rna") %in% names(fam_mk)))
rep_lines <- fam_mk$cell_line[fam_mk$patient_representative & fam_mk$has_rna]
rep_keep  <- ann$cell_line %in% rep_lines
stopifnot("patient-representative RNA set must be the expected 28 models" =
            sum(rep_keep) == 28L)

marker_perm <- bind_rows(
  perm_null(mk_mat, ann$subtype) %>%
    mutate(unit = "31 line models", n_samples = ncol(mk_mat), .before = 1),
  perm_null(mk_mat[, rep_keep, drop = FALSE], ann$subtype[rep_keep]) %>%
    mutate(unit = "28 patient representatives", n_samples = sum(rep_keep), .before = 1),
  perm_null(mk_mat[, rep_keep, drop = FALSE], ann$subtype[rep_keep],
            strata = ann$site[rep_keep]) %>%
    mutate(unit = "28 patient representatives; within-centre permutation",
           n_samples = sum(rep_keep), .before = 1)
) %>%
  mutate(binomial_p = vapply(seq_len(dplyr::n()),
                             function(i) bt(observed[i], n_markers[i]), numeric(1)),
         null_model = ifelse(grepl("within-centre", unit),
           "joint histotype-label permutation within contributing centre; patient representatives; source composition preserved",
           "joint histotype-label permutation; marker correlation and group sizes preserved; source not conditioned"),
         seed = SEED)
readr::write_csv(marker_perm, file.path(OUT, "rna_marker_recovery_permutation.csv"))

cat("\n=== Marker recovery against a JOINT label-permutation null (correlation-preserving) ===\n")
print(as.data.frame(marker_perm %>% select(unit, n_samples, observed, n_markers,
                                           null_mean, null_p95, permutation_p, binomial_p) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
cat("Use the patient-representative, within-centre permutation for the source-aware association.\n")
cat("The retired binomial p assumed marker independence, which co-regulation violates.\n")

# -----------------------------------------------------------------------------
# 2b. ENDOMETRIOID markers, reported at the n that actually exists
# -----------------------------------------------------------------------------
# The annotated EC group is 2 RNA models (TOV112D, VOA4395). Discordant SWI/SNF
# evidence for TOV112D motivates reporting each model separately while retaining
# the supplied histotype. These are per-model values and z-scores across all 31
# models, not an inference about a histotype based on a single model.
ec_syms <- c(markers$symbol[markers$marks == "EC (endometrioid)"], "VIM")
ec_z    <- zscore(logtpm[ec_syms, , drop = FALSE])
ec_tbl  <- as_tibble(logtpm[ec_syms, , drop = FALSE], rownames = "symbol") %>%
  pivot_longer(-symbol, names_to = "cell_line", values_to = "log2_tpm") %>%
  left_join(as_tibble(ec_z, rownames = "symbol") %>%
              pivot_longer(-symbol, names_to = "cell_line", values_to = "z_all31"),
            by = c("symbol", "cell_line")) %>%
  left_join(ann %>% select(cell_line, subtype), by = "cell_line") %>%
  mutate(across(c(log2_tpm, z_all31), ~round(.x, 3)),
         in_hollis_ec_panel = TRUE,
         scored_group = ifelse(symbol == "VIM", "MMMT (EMT) - not double-counted for EC",
                               "EC (endometrioid)")) %>%
  arrange(symbol, subtype, cell_line)
readr::write_csv(ec_tbl, file.path(OUT, "rna_ec_markers.csv"))

cat("\n=== EC marker check (Hollis 2020 panel), reported at the n that exists ===\n")
cat("Annotated EC group n=2 (TOV112D + VOA4395); discordant SWI/SNF evidence for\n")
cat("TOV112D warrants histotype review. Individual-model z-scores use all 31 models.\n")
ec_wide <- ec_tbl %>% filter(cell_line %in% c("TOV112D", "VOA4395")) %>%
  select(symbol, cell_line, log2_tpm, z_all31) %>%
  pivot_wider(names_from = cell_line, values_from = c(log2_tpm, z_all31))
print(as.data.frame(ec_wide), row.names = FALSE)
cat("\nEC subtype means vs the best-scoring subtype for each EC marker:\n")
print(as.data.frame(marker_summary %>% filter(marks_code == "EC" | symbol == "VIM") %>%
        select(symbol, direction, EC_mean, mean_in, rank_in, top_subtype, lands_right)),
      row.names = FALSE)
cat("VIM is the 4th Hollis EC marker; it is scored under MMMT/EMT above and is NOT\n")
cat("  double-counted for EC (broadly expressed in 2D culture -> poor discriminator).\n")

# -----------------------------------------------------------------------------
# 3. SWI/SNF check (SCCOHT + TOV112D) — feeds authentication (Phase 4)
# -----------------------------------------------------------------------------
swi_genes <- c("SMARCA4","SMARCA2","ARID1A","SMARCB1")
stopifnot(all(swi_genes %in% rownames(logtpm)))
swi_log <- logtpm[swi_genes, , drop = FALSE]
swi_z   <- zscore(swi_log)

swisnf <- tibble(cell_line = colnames(logtpm),
                 subtype   = ann$subtype) %>%
  bind_cols(as_tibble(t(round(2^swi_log - 1, 2))) %>%
              rename_with(~ paste0(.x, "_tpm"))) %>%
  bind_cols(as_tibble(t(round(swi_z, 2))) %>%
              rename_with(~ paste0(.x, "_z")))
# ranks across all 31 lines (1 = lowest expression)
for (g in swi_genes)
  swisnf[[paste0(g, "_rank")]] <- rank(swisnf[[paste0(g, "_tpm")]], ties.method = "min")
swisnf <- swisnf %>% arrange(SMARCA4_tpm)
readr::write_csv(swisnf, file.path(OUT, "rna_swisnf.csv"))

key <- c("BIN67","COV434","TOV112D")
cat("\n=== SWI/SNF: key lines (rank 1 = lowest of 31) ===\n")
print(swisnf %>% filter(cell_line %in% key) %>%
        select(cell_line, subtype,
               SMARCA4_tpm, SMARCA4_rank, SMARCA2_tpm, SMARCA2_rank,
               ARID1A_tpm, SMARCB1_tpm) %>% as.data.frame(), row.names = FALSE)
cat("\nSMARCA4 5 lowest:", paste(head(swisnf$cell_line, 5), collapse = ", "), "\n")
cat("SMARCA2 5 lowest:",
    paste((swisnf %>% arrange(SMARCA2_tpm) %>% pull(cell_line))[1:5], collapse = ", "), "\n")

# --- SWI/SNF heatmap (genes x lines, ordered by subtype; key lines highlighted)
hl <- ifelse(colnames(swi_z) %in% key, "bold", "plain")
col_ha2 <- HeatmapAnnotation(subtype = ann$subtype, col = list(subtype = sub_cols),
                             annotation_name_gp = grid::gpar(fontsize = 8))
pdf(file.path(FIGS, "04_rna_swisnf.pdf"), width = 10, height = 3.6)
draw(Heatmap(swi_z, name = "z(log2 TPM)", col = z_fun,
             top_annotation = col_ha2,
             cluster_rows = FALSE, column_split = ann$subtype,
             cluster_columns = FALSE, cluster_column_slices = FALSE,
             row_names_gp = grid::gpar(fontsize = 9),
             column_names_gp = grid::gpar(fontsize = 7, fontface = hl),
             column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
             column_title = "SWI/SNF subunit expression (bold = SCCOHT-candidate + TOV112D)"),
     heatmap_legend_side = "right")
dev.off()

# -----------------------------------------------------------------------------
# 4. Hallmark gene-set scoring (singscore) — MSigDB Hallmark v7.4
# -----------------------------------------------------------------------------
gmt_file <- file.path(DATA, "h.all.v7.4.symbols.gmt.txt")     # read-only source
stopifnot("Hallmark GMT not found" = file.exists(gmt_file))
hallmark <- fgsea::gmtPathways(gmt_file)
message(sprintf("Hallmark sets: %d (v7.4)", length(hallmark)))

rankData <- rankGenes(logtpm)                       # rank-based; log vs raw TPM identical ranks
score_one <- function(gs) {
  gs <- intersect(gs, rownames(logtpm))
  simpleScore(rankData, upSet = gs)$TotalScore
}
sc <- vapply(hallmark, score_one, numeric(ncol(logtpm)))   # lines x sets
rownames(sc) <- colnames(logtpm)
scores_out <- as_tibble(sc, rownames = "cell_line") %>%
  left_join(ann %>% select(cell_line, subtype), by = "cell_line") %>%
  relocate(subtype, .after = cell_line)
readr::write_csv(scores_out, file.path(OUT, "rna_geneset_scores.csv"))

# z-score across lines per Hallmark set for the heatmap
sc_z <- zscore(t(sc))                               # sets x lines
rownames(sc_z) <- sub("^HALLMARK_", "", rownames(sc_z))

pdf(file.path(FIGS, "04_rna_hallmark_singscore_heatmap.pdf"), width = 11, height = 11)
draw(Heatmap(sc_z, name = "z(singscore)", col = z_fun,
             top_annotation = col_ha,
             column_split = ann$subtype, cluster_columns = TRUE, cluster_column_slices = FALSE,
             cluster_rows = TRUE,
             row_names_gp = grid::gpar(fontsize = 7),
             column_names_gp = grid::gpar(fontsize = 7),
             column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
             column_title = "Hallmark singscore (z across lines) x cell lines"),
     heatmap_legend_side = "right")
dev.off()

# --- which Hallmark sets track a subtype? (mean z per subtype; report extremes)
hz <- as_tibble(sc_z, rownames = "hallmark") %>%
  pivot_longer(-hallmark, names_to = "cell_line", values_to = "z") %>%
  left_join(ann %>% select(cell_line, subtype), by = "cell_line") %>%
  group_by(hallmark, subtype) %>% summarise(mz = mean(z), .groups = "drop")
top_sets <- hz %>% group_by(subtype) %>% slice_max(mz, n = 3, with_ties = FALSE) %>%
  arrange(subtype, desc(mz)) %>% ungroup()
cat("\n=== Top Hallmark sets by mean z per subtype (validation read; n annotated) ===\n")
for (st in subtypes) {
  hh <- top_sets %>% filter(subtype == st)
  cat(sprintf("[%s n=%d] %s\n", st, as.integer(n_by[st]),
              paste(sprintf("%s(%.2f)", hh$hallmark, hh$mz), collapse = ", ")))
}

cat("\nOutputs: output/rna_markers_summary.csv, rna_ec_markers.csv, rna_swisnf.csv, rna_geneset_scores.csv\n")
cat("Figures: figs/04_rna_markers_heatmap.pdf, 04_rna_swisnf.pdf, 04_rna_hallmark_singscore_heatmap.pdf\n")
cat("NOTE: Hallmark scoring here is PER LINE (n=31), not per patient.\n")

# -----------------------------------------------------------------------------
# 5. Environment record
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
write_session_info("04_rna_markers_genesets")

message("\n04_rna_markers_genesets.R complete.")
