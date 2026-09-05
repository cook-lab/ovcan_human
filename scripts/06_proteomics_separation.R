# =============================================================================
# Script: 06_proteomics_separation.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Does the TMT proteomics recapitulate histological subtypes, and is
#          any separation confounded by batch (TMT plex) or source site?
#          PCA / t-SNE / silhouette / correlation clustering with an explicit
#          subtype-vs-PLEX-vs-SITE confounder check (mirrors the RNA check in
#          output/rna_pc_confounder.csv, adding the plex axis). Plus canonical
#          subtype protein markers.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   2 (Proteomics) — step 2 of 2
# -----------------------------------------------------------------------------
# BENCHMARK EXPECTATION: proteomic subtype separation is typically WEAKER than
# RNA; the n=2 rare-subtype groups (EC, MMMT, SCCOHT) cannot support reliable
# silhouettes and may be negative. HGS is 100% Mes-Masson while Huntsman is
# 100% non-HGS, so subtype and site are partly confounded by design.
# MISSINGNESS: TMT missingness is per-plex. To avoid plex-structured missingness
# (or imputation) driving the components, the NA-intolerant multivariate methods
# here use the COMPLETE-CASE proteins (quantified in all 5 plexes). The
# documented reusable matrix (output/prot_abundance_matrix.csv) keeps all 8,427
# proteins; the recommended reuse filter is presence >=50% of lines (7,734).
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(tidyverse); library(matrixStats)
  library(Rtsne); library(cluster); library(ComplexHeatmap); library(circlize)
})

SUBTYPE_LEVELS <- c("HGS", "CC", "EC", "MC", "MMMT", "SCCOHT")

# 1. Load matrix + annotation --------------------------------------------------
obj  <- readRDS(file.path(OUT, "prot_matrix.rds"))
mat  <- obj$mat            # 8,427 proteins x 31 lines (log2 relative abundance, NAs)
feat <- obj$feat
ann  <- obj$sample_ann %>%
  transmute(cell_line, subtype = factor(subtype_ss, levels = SUBTYPE_LEVELS),
            site = factor(site_ss), plex = factor(plex)) %>%
  as.data.frame()
stopifnot(identical(colnames(mat), ann$cell_line))

sub_cols  <- setNames(RColorBrewer::brewer.pal(6, "Dark2"), SUBTYPE_LEVELS)
site_cols <- setNames(c("#C2410C", "#0F172A", "#6B7280"),
                      c("Mes-Masson", "Huntsman", "Huntsman/Vanderhyden"))
plex_cols <- setNames(viridisLite::viridis(5), as.character(1:5))

# 2. Complete-case matrix for multivariate analysis ----------------------------
cc  <- rowSums(is.na(mat)) == 0
mcc <- mat[cc, , drop = FALSE]
message(sprintf("Multivariate input: %d complete-case proteins (of %d) x %d lines",
                nrow(mcc), nrow(mat), ncol(mcc)))

# 3. Top variable proteins + PCA (center, no scale — mirrors RNA) --------------
rv  <- rowVars(mcc); top <- head(order(rv, decreasing = TRUE), 2000)
pca <- prcomp(t(mcc[top, ]), center = TRUE, scale. = FALSE)
pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pcd <- as_tibble(pca$x[, 1:10]) %>%
  mutate(cell_line = rownames(pca$x), subtype = ann$subtype,
         site = ann$site, plex = ann$plex)
saveRDS(pca, file.path(OUT, "prot_pca.rds"))

pc_plot <- function(colvar, cols, ttl)
  ggplot(pcd, aes(PC1, PC2, colour = .data[[colvar]])) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = cell_line), size = 2.3, max.overlaps = 20) +
    scale_colour_manual(values = cols) +
    labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2]),
         colour = colvar, title = ttl) +
    theme_minimal(base_size = 11)
ggsave(file.path(FIGS, "06_prot_pca_subtype.pdf"),
       pc_plot("subtype", sub_cols, "Protein PCA by subtype"), width = 7.5, height = 6)
ggsave(file.path(FIGS, "06_prot_pca_plex.pdf"),
       pc_plot("plex", plex_cols, "Protein PCA by TMT plex (batch)"), width = 7.5, height = 6)
ggsave(file.path(FIGS, "06_prot_pca_site.pdf"),
       pc_plot("site", site_cols, "Protein PCA by source site"), width = 7.5, height = 6)

# 4. t-SNE (perplexity 5; n=31) ------------------------------------------------
set.seed(SEED)
ts  <- Rtsne(pca$x[, 1:10], perplexity = 5, theta = 0.0, pca = FALSE)$Y
tsd <- tibble(t1 = ts[, 1], t2 = ts[, 2], subtype = ann$subtype,
              site = ann$site, cell_line = ann$cell_line)
ggsave(file.path(FIGS, "06_prot_tsne_subtype.pdf"),
       ggplot(tsd, aes(t1, t2, colour = subtype)) + geom_point(size = 3) +
         ggrepel::geom_text_repel(aes(label = cell_line), size = 2.3, max.overlaps = 20) +
         scale_colour_manual(values = sub_cols) +
         labs(title = "Protein t-SNE (perplexity 5)") + theme_minimal(base_size = 11),
       width = 7.5, height = 6)

# 5. Silhouette by subtype (Euclidean on PC1-10; annotate n) -------------------
d   <- dist(pca$x[, 1:10])
sil <- silhouette(as.integer(ann$subtype), d)
sil_tab <- tibble(subtype = levels(ann$subtype)[sil[, "cluster"]],
                  width = sil[, "sil_width"]) %>%
  group_by(subtype) %>% summarise(n = n(), mean_sil = mean(width), .groups = "drop") %>%
  mutate(subtype = factor(subtype, levels = SUBTYPE_LEVELS)) %>% arrange(subtype)
readr::write_csv(sil_tab, file.path(OUT, "prot_silhouette.csv"))

# 5b. BOTH modalities in one table, labelled ------------------------------------
#     The RNA silhouettes were the ones that reached the manuscript, unlabelled by
#     modality, and the protein values differ sharply in the direction that
#     matters (clear cell is NEGATIVE in protein, SCCOHT ~0). Writing them side by
#     side with an explicit `modality` column makes them impossible to confuse.
#     This script runs after 02, so rna_silhouette.csv is on disk.
rna_sil_f <- file.path(OUT, "rna_silhouette.csv")
stopifnot("rna_silhouette.csv missing - run 02_rna_separation.R first" = file.exists(rna_sil_f))
sil_both <- bind_rows(
  readr::read_csv(rna_sil_f, show_col_types = FALSE) %>% mutate(modality = "RNA"),
  sil_tab %>% mutate(subtype = as.character(subtype), modality = "protein"),
  tibble(subtype = "ALL", n = ncol(mcc),
         mean_sil = mean(sil[, "sil_width"]), modality = "protein")) %>%
  transmute(modality = factor(modality, levels = c("RNA", "protein")),
            subtype = factor(subtype, levels = c(SUBTYPE_LEVELS, "ALL")),
            n, mean_sil = round(mean_sil, 4),
            unit = "cell line (per-line, NOT per-patient)") %>%
  arrange(modality, subtype)
readr::write_csv(sil_both, file.path(OUT, "silhouette_by_modality.csv"))
cat("\nSilhouettes, BOTH modalities (subtype = ALL is the overall mean):\n")
print(as.data.frame(sil_both), row.names = FALSE)

# 6. Confounder check: R^2 of each PC ~ subtype vs plex vs site ----------------
#    Raw R^2 favours factors with more levels (subtype 5 df, plex 4 df, site 2
#    df), so adjusted R^2 is reported alongside. This mirrors + extends the RNA
#    check (which had only subtype + site) with the TMT plex axis.
r2  <- function(pc, f) summary(lm(pcd[[pc]] ~ f))$r.squared
ar2 <- function(pc, f) summary(lm(pcd[[pc]] ~ f))$adj.r.squared
conf <- tibble(
  PC = paste0("PC", 1:5), var_pct = pv[1:5],
  r2_subtype = sapply(paste0("PC", 1:5), r2,  f = ann$subtype),
  r2_plex    = sapply(paste0("PC", 1:5), r2,  f = ann$plex),
  r2_site    = sapply(paste0("PC", 1:5), r2,  f = ann$site),
  adjr2_subtype = sapply(paste0("PC", 1:5), ar2, f = ann$subtype),
  adjr2_plex    = sapply(paste0("PC", 1:5), ar2, f = ann$plex),
  adjr2_site    = sapply(paste0("PC", 1:5), ar2, f = ann$site))
readr::write_csv(conf, file.path(OUT, "prot_pc_confounder.csv"))

# 7. Sample-correlation heatmap (Spearman) with subtype+plex+site annotation ---
sc <- cor(mcc, method = "spearman")
ha <- HeatmapAnnotation(subtype = ann$subtype, plex = ann$plex, site = ann$site,
                        col = list(subtype = sub_cols, plex = plex_cols, site = site_cols))
pdf(file.path(FIGS, "06_prot_spearman_heatmap.pdf"), width = 8.5, height = 7.5)
draw(Heatmap(sc, name = "Spearman", top_annotation = ha,
             clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
             col = colorRamp2(quantile(sc, c(.5, .75, 1)), viridisLite::mako(3)),
             show_row_names = FALSE, column_names_gp = grid::gpar(fontsize = 7)))
dev.off()

# 8. Canonical subtype protein markers ----------------------------------------
#    Only markers present in the proteome are used (CDX2/MUC2/SATB2 not detected).
markers <- tribble(
  ~gene,     ~category,
  "PAX8",    "HGS/Mullerian", "WT1", "HGS/Mullerian", "MUC16", "HGS (CA125)", "WFDC2", "HGS (HE4)",
  "MECOM",   "HGS",
  "HNF1B",   "Clear cell",    "NAPSA", "Clear cell",
  "MUC5AC",  "Mucinous",      "MUC5B", "Mucinous",     "TFF1",  "Mucinous",   "KRT7", "Epithelial",
  "EPCAM",   "Epithelial",
  "SMARCA4", "SWI/SNF (SCCOHT/dediff)", "SMARCA2", "SWI/SNF (SCCOHT/dediff)",
  "ARID1A",  "SWI/SNF",       "ARID1B", "SWI/SNF")
mk <- markers %>% filter(gene %in% feat$symbol)
# 13 symbols carry more than one Uniprot accession, so a bare match() on symbol
# silently took whichever accession came first in the search table. Restrict the
# lookup to the deterministic symbol representative defined in script 05
# (most non-missing lines, then most peptides, then lowest q, then accession).
rep_feat <- feat[feat$symbol_representative, ]
stopifnot("symbol representative is not unique" = !any(duplicated(rep_feat$symbol)))
mk_mat <- mat[rep_feat$row[match(mk$gene, rep_feat$symbol)], , drop = FALSE]
rownames(mk_mat) <- mk$gene
# row z-score (per protein across lines), NA-aware (scale() is not)
zsc <- t(apply(mk_mat, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
ord <- order(ann$subtype, ann$cell_line)
ha2 <- HeatmapAnnotation(subtype = ann$subtype[ord], col = list(subtype = sub_cols))
pdf(file.path(FIGS, "06_prot_markers.pdf"), width = 9.5, height = 6)
draw(Heatmap(zsc[, ord], name = "row z\n(log2 abund)",
             top_annotation = ha2, cluster_columns = FALSE,
             column_split = factor(ann$subtype[ord], levels = SUBTYPE_LEVELS),
             column_title_gp = grid::gpar(fontsize = 9),
             row_split = factor(mk$category, levels = unique(mk$category)),
             row_title_gp = grid::gpar(fontsize = 8), row_title_rot = 0,
             cluster_rows = FALSE,
             col = colorRamp2(c(-2, 0, 2), c("#2166AC", "grey95", "#B2182B")),
             column_names_gp = grid::gpar(fontsize = 7),
             row_names_gp = grid::gpar(fontsize = 8), na_col = "grey80"))
dev.off()

# per-subtype mean marker abundance (validation table) -------------------------
mk_means <- as_tibble(mk_mat, rownames = "gene") %>%
  pivot_longer(-gene, names_to = "cell_line", values_to = "abund") %>%
  left_join(tibble(cell_line = ann$cell_line, subtype = ann$subtype), by = "cell_line") %>%
  group_by(gene, subtype) %>% summarise(mean_abund = mean(abund, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = subtype, values_from = mean_abund)

# 9. Console report + QC-gate check --------------------------------------------
cat("\n=== PROTEOMICS SEPARATION (script 06) ===\n")
cat(sprintf("Multivariate input: %d complete-case proteins x %d lines; top 2000 variable used for PCA\n",
            nrow(mcc), ncol(mcc)))
cat("\nPCA variance explained (PC1-5):\n"); print(pv[1:5])
cat("\nSilhouette by subtype (Euclidean, PC1-10; n=2 groups unreliable):\n")
print(as.data.frame(sil_tab))
cat(sprintf("Overall mean silhouette: %.3f\n", mean(sil[, "sil_width"])))
cat("\nConfounder check — variance in each PC explained by subtype vs plex vs site:\n")
print(as.data.frame(conf), digits = 3)
cat("\nCanonical marker mean log2 relative abundance by subtype:\n")
print(as.data.frame(mk_means), digits = 4)

# QC-GATE: does batch dominate subtype on the top PCs?
gate_flag <- with(conf, any(r2_plex[1:2] > r2_subtype[1:2]))
site_flag <- with(conf, any(r2_site[1:2] > r2_subtype[1:2]))
cat("\n--- QC GATE ---\n")
cat(sprintf("PLEX exceeds subtype on PC1/PC2: %s\n", gate_flag))
cat(sprintf("SITE exceeds subtype on PC1/PC2: %s (note: subtype/site confounded by design)\n", site_flag))
if (gate_flag)
  cat("!! FLAG: TMT plex explains more top-PC variance than subtype -> batch may dominate; interpret separation cautiously.\n")
if (!gate_flag)
  cat("OK: subtype structure is not dominated by TMT plex on the top PCs (consistent with effective PIS normalization).\n")
cat("\nNB for reporting: prot_pc_confounder.csv is the PC-level view and is NOT reassuring in the\n")
cat("  way the per-protein medians are (site and plex load substantially on the leading protein PCs).\n")
cat("  Report the per-protein medians and this table together, not one instead of the other.\n")

# 10. Environment record -------------------------------------------------------
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
write_session_info("06_proteomics_separation")

message("\n06_proteomics_separation.R complete. Figures in figs/06_* ; tables in output/prot_*.")
