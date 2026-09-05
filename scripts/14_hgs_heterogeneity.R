# =============================================================================
# Script: 14_hgs_heterogeneity.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: A USAGE EXAMPLE for the descriptor (Fig 5): within the best-sampled
#          subtype (HGS, n=15 RNA lines) show how a user can stratify models by
#          pathway/gene-set activity to pick a line matching a biology of
#          interest. Cluster on Hallmark singscore (primary) and corroborate with
#          PROGENy pathway activity. DESCRIPTIVE ONLY (n=15) - no survival, no
#          discovery, no "novel subtype" claims.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   5 (integration & usage) - step 3 of 3
# =============================================================================
# Assay-aware notes / decisions:
#   - Input = per-line Hallmark singscore (output/rna_geneset_scores.csv, from
#     script 04; MSigDB Hallmark v7.4). We cluster the 15 HGS lines on those
#     scores z-scored ACROSS the 15 lines (so each pathway contributes equally
#     and clustering reflects RELATIVE within-HGS differences).
#   - PROGENy (14 pathways, top=100 responsive genes, scaled across the 15 lines)
#     is computed on the symbol-level log2 TPM as an ORTHOGONAL corroboration of
#     the Hallmark-derived groups (different gene sets + a footprint model, not a
#     simple overlap of the same genes).
#   - Confounders: all 15 HGS lines are Mes-Masson/CHUM -> NO source-site
#     confound here (unlike cross-subtype analyses). We still report a
#     proliferation read (Hallmark G2M/E2F/MYC) so groups are not just "fast vs
#     slow cycling", and we note where same-patient isolate families land.
#   - k=3 is chosen for description (also report k=2/k=4 membership); with n=15
#     these are illustrative strata, not defined molecular subtypes.
# =============================================================================

source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(matrixStats); library(progeny)
  library(ComplexHeatmap); library(circlize); library(cluster)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; slice <- dplyr::slice
set.seed(SEED)

theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace% theme(
    text = element_text(colour = "black"),
    plot.title = element_text(size = rel(1.15), hjust = 0, margin = margin(b = 6)),
    plot.subtitle = element_text(size = rel(0.85), hjust = 0, colour = "grey30", margin = margin(b = 6)),
    axis.title = element_text(size = rel(1.0)), axis.text = element_text(size = rel(0.85), colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    legend.title = element_text(size = rel(0.85)), legend.text = element_text(size = rel(0.8)),
    legend.key = element_blank(), legend.background = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid = element_blank(), strip.background = element_blank(),
    strip.text = element_text(size = rel(0.95), face = "bold"), plot.margin = margin(8, 8, 8, 8))
}
okabe_ito <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
z_fun  <- colorRamp2(c(-2, 0, 2), c(COOK_NAVY, "white", COOK_RUST))     # match script 04
zrows  <- function(mat) t(scale(t(mat)))                               # z across cols, per row

# -----------------------------------------------------------------------------
# 1. HGS lines + Hallmark singscore matrix (z across the 15 lines)
# -----------------------------------------------------------------------------
gs <- readr::read_csv(file.path(OUT, "rna_geneset_scores.csv"), show_col_types = FALSE)
hgs <- gs %>% filter(subtype == "HGS")
hgs_lines <- hgs$cell_line
n_hgs <- length(hgs_lines)
stopifnot(n_hgs == 15)
message(sprintf("HGS lines: %d", n_hgs))

H <- as.matrix(hgs[, grep("^HALLMARK_", colnames(hgs))]); rownames(H) <- hgs_lines
colnames(H) <- sub("^HALLMARK_", "", colnames(H))
Hz <- scale(H)                                     # lines x sets, z per set across lines

# -----------------------------------------------------------------------------
# 2. PROGENy pathway activity on symbol-level log2 TPM (orthogonal view)
# -----------------------------------------------------------------------------
tpm <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"), show_col_types = FALSE)
gene_meta <- t2g %>% distinct(ensembl_gene_id, external_gene_name) %>%
  filter(!is.na(external_gene_name), external_gene_name != "")
m_id <- as.matrix(tpm[, -1]); rownames(m_id) <- tpm$gene_id
id2sym <- setNames(gene_meta$external_gene_name, gene_meta$ensembl_gene_id)
common <- intersect(rownames(m_id), names(id2sym))
logtpm <- log2(rowsum(m_id[common, , drop = FALSE], group = id2sym[common]) + 1)

pg <- progeny::progeny(logtpm[, hgs_lines], scale = TRUE, organism = "Human",
                       top = 100, perm = 1)          # lines x 14 pathways, z-scaled across lines
message(sprintf("PROGENy: %d lines x %d pathways", nrow(pg), ncol(pg)))

# -----------------------------------------------------------------------------
# 3. Cluster the 15 HGS lines on Hallmark z (ward.D2 / Euclidean); k=3 primary
# -----------------------------------------------------------------------------
d  <- dist(Hz)
hc <- hclust(d, method = "ward.D2")
cl3 <- cutree(hc, 3); cl2 <- cutree(hc, 2); cl4 <- cutree(hc, 4)

# Independent-patient sensitivity of the illustrative partition. These metrics
# assess dependence on related sublines; they do not validate molecular subtypes.
fm <- readr::read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)
rep_hgs <- intersect(hgs_lines, fm$cell_line[fm$patient_representative])
stopifnot(length(rep_hgs) == 12L)
Hz_rep <- scale(H[rep_hgs, , drop = FALSE])
d_rep <- dist(Hz_rep); hc_rep <- hclust(d_rep, method = "ward.D2")
ari <- function(a, b) {
  tab <- table(a, b); choose2 <- function(x) x * (x - 1) / 2
  n2 <- choose2(sum(tab)); expected <- sum(choose2(rowSums(tab))) * sum(choose2(colSums(tab))) / n2
  (sum(choose2(tab)) - expected) /
    ((sum(choose2(rowSums(tab))) + sum(choose2(colSums(tab)))) / 2 - expected)
}
cluster_sensitivity <- purrr::map_dfr(2:4, function(k) {
  a <- cutree(hc, k); b <- cutree(hc_rep, k)
  tibble(k = k, n_models = length(hgs_lines), n_patients = length(rep_hgs),
         silhouette_models = mean(cluster::silhouette(a, d)[, "sil_width"]),
         silhouette_patient_representatives = mean(cluster::silhouette(b, d_rep)[, "sil_width"]),
         adjusted_rand_on_shared_representatives = ari(a[rep_hgs], b[rep_hgs]),
         interpretation = "illustrative partitions; k was not prespecified or externally validated; does not establish stable molecular subtypes")
})
readr::write_csv(cluster_sensitivity, file.path(OUT, "hgs_cluster_patient_sensitivity.csv"))

# robustness: cluster independently on PROGENy and cross-tabulate with Hallmark k=3
hc_pg  <- hclust(dist(scale(pg)), method = "ward.D2")
cl3_pg <- cutree(hc_pg, 3)
xtab   <- table(hallmark = cl3[hgs_lines], progeny = cl3_pg[hgs_lines])

# -----------------------------------------------------------------------------
# 4. Characterise clusters; assign descriptive labels from theme scores
# -----------------------------------------------------------------------------
theme_sets <- list(
  inflammatory = c("TNFA_SIGNALING_VIA_NFKB","INFLAMMATORY_RESPONSE","IL6_JAK_STAT3_SIGNALING",
                   "INTERFERON_GAMMA_RESPONSE","EPITHELIAL_MESENCHYMAL_TRANSITION","IL2_STAT5_SIGNALING"),
  hypoxic_glycolytic = c("HYPOXIA","GLYCOLYSIS"),
  proliferation = c("G2M_CHECKPOINT","E2F_TARGETS","MYC_TARGETS_V1","MITOTIC_SPINDLE"))
theme_score <- function(sets) rowMeans(Hz[, intersect(sets, colnames(Hz)), drop = FALSE])
themes <- tibble(cell_line = hgs_lines,
                 inflammatory = theme_score(theme_sets$inflammatory),
                 hypoxic_glycolytic = theme_score(theme_sets$hypoxic_glycolytic),
                 proliferation = theme_score(theme_sets$proliferation),
                 cluster = cl3[hgs_lines])
clabel_by <- themes %>% group_by(cluster) %>%
  summarise(inflammatory = mean(inflammatory), hypoxic_glycolytic = mean(hypoxic_glycolytic),
            proliferation = mean(proliferation), .groups = "drop")
# greedy label assignment: inflammatory -> hypoxic-glycolytic -> remaining = low-signaling
lab <- setNames(rep(NA_character_, 3), clabel_by$cluster)
c_inf <- clabel_by$cluster[which.max(clabel_by$inflammatory)]; lab[as.character(c_inf)] <- "Inflammatory/NF-kB-EMT"
rem <- clabel_by %>% filter(cluster != c_inf)
c_hyp <- rem$cluster[which.max(rem$hypoxic_glycolytic)];       lab[as.character(c_hyp)] <- "Hypoxic-glycolytic"
c_low <- setdiff(clabel_by$cluster, c(c_inf, c_hyp));          lab[as.character(c_low)] <- "Low-signaling"
cluster_label <- lab[as.character(cl3[hgs_lines])]

# mean Hallmark z per cluster (characterisation table) + top sets
cl_means <- as_tibble(Hz, rownames = "cell_line") %>%
  mutate(cluster = cl3[cell_line]) %>%
  pivot_longer(-c(cell_line, cluster), names_to = "hallmark", values_to = "z") %>%
  group_by(cluster, hallmark) %>% summarise(mean_z = mean(z), .groups = "drop")
cl_means_wide <- cl_means %>% pivot_wider(names_from = cluster, values_from = mean_z,
                                          names_prefix = "cluster")
readr::write_csv(cl_means_wide %>% mutate(across(where(is.numeric), ~round(.x, 3))),
                 file.path(OUT, "hgs_hallmark_cluster_means.csv"))

# -----------------------------------------------------------------------------
# 5. Output: per-line assignment + PROGENy + proliferation  -> hgs_heterogeneity.csv
# -----------------------------------------------------------------------------
family <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  filter(cell_line %in% hgs_lines) %>% select(cell_line, notes)
out <- tibble(cell_line = hgs_lines,
              cluster = as.integer(cl3[hgs_lines]),
              cluster_label = cluster_label,
              cluster_k2 = as.integer(cl2[hgs_lines]),
              cluster_k4 = as.integer(cl4[hgs_lines]),
              proliferation_z = round(themes$proliferation, 3),
              inflammatory_z = round(themes$inflammatory, 3),
              hypoxic_glycolytic_z = round(themes$hypoxic_glycolytic, 3),
              patient_id = fm$patient_id[match(hgs_lines, fm$cell_line)],
              patient_representative = hgs_lines %in% rep_hgs,
              interpretation = "exploratory k=3 pathway-score partition; not a validated molecular subtype") %>%
  bind_cols(as_tibble(round(pg[hgs_lines, ], 3)) %>% rename_with(~ paste0("progeny_", .x))) %>%
  arrange(cluster, cell_line)
readr::write_csv(out, file.path(OUT, "hgs_heterogeneity.csv"))

# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------
cl_cols <- setNames(okabe_ito[1:3], as.character(sort(unique(cl3))))
lab_for_leg <- tapply(cluster_label, cl3[hgs_lines], function(x) x[1])
# order columns by cluster then name
col_ord <- out$cell_line
prolif_vec <- setNames(themes$proliferation, themes$cell_line)[col_ord]

# 6a. Hallmark heatmap (sets x lines), split by cluster
top_ann <- HeatmapAnnotation(
  cluster = factor(cl3[col_ord]),
  proliferation = prolif_vec,
  col = list(cluster = cl_cols,
             proliferation = colorRamp2(c(-1.5, 0, 1.5), c(COOK_NAVY, "white", COOK_RUST))),
  annotation_legend_param = list(cluster = list(title = "cluster")),
  annotation_name_gp = grid::gpar(fontsize = 8))
Hz_plot <- t(Hz)[, col_ord]                        # sets x lines
# order rows to surface structure
rv <- rowVars(Hz_plot); Hz_plot <- Hz_plot[order(rv, decreasing = TRUE), ]
pdf(file.path(FIGS, "14_hgs_hallmark_heatmap.pdf"), width = 9, height = 11)
draw(Heatmap(Hz_plot, name = "z(singscore)", col = z_fun, top_annotation = top_ann,
             column_split = factor(cl3[col_ord]), cluster_column_slices = FALSE,
             cluster_columns = TRUE, cluster_rows = TRUE,
             row_names_gp = grid::gpar(fontsize = 6.5), column_names_gp = grid::gpar(fontsize = 8),
             column_title = "Within-HGS heterogeneity: Hallmark singscore (z across 15 HGS lines)",
             column_title_gp = grid::gpar(fontsize = 10, fontface = "bold")),
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# 6b. PROGENy heatmap (14 pathways x lines), columns in the same cluster order
pgz <- t(pg)[, col_ord]                            # pathways x lines (already scaled)
pdf(file.path(FIGS, "14_hgs_progeny_heatmap.pdf"), width = 9, height = 5)
draw(Heatmap(pgz, name = "PROGENy z", col = z_fun, top_annotation = top_ann,
             column_split = factor(cl3[col_ord]), cluster_column_slices = FALSE,
             cluster_columns = TRUE, cluster_rows = TRUE,
             row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 8),
             column_title = "PROGENy pathway activity (orthogonal corroboration)",
             column_title_gp = grid::gpar(fontsize = 10, fontface = "bold")),
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# 6c. Theme-score scatter (inflammatory vs hypoxic-glycolytic), sized by proliferation
pC <- themes %>% mutate(cluster = factor(cl3[cell_line]), label = cluster_label) %>%
  ggplot(aes(inflammatory, hypoxic_glycolytic, colour = cluster)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(size = proliferation), alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = cell_line), size = 2.4, max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = cl_cols, name = "cluster",
                      labels = paste0(names(cl_cols), ": ", lab_for_leg[names(cl_cols)])) +
  scale_size_continuous(name = "proliferation z", range = c(1.5, 6)) +
  labs(title = "HGS strata by pathway theme (n = 15)",
       subtitle = "descriptive; point size = Hallmark proliferation (G2M/E2F/MYC/spindle) z",
       x = "Inflammatory / NF-kB-EMT score (mean Hallmark z)",
       y = "Hypoxic-glycolytic score (mean Hallmark z)") +
  theme_lab() + theme(legend.position = "right")
ggsave(file.path(FIGS, "14_hgs_theme_scatter.pdf"), pC, width = 8, height = 5.6)

# -----------------------------------------------------------------------------
# 7. Console report
# -----------------------------------------------------------------------------
cat("\n=== Within-HGS strata (Hallmark ward.D2, k=3; DESCRIPTIVE, n=15) ===\n")
for (k in sort(unique(cl3))) {
  ln <- names(cl3)[cl3 == k]
  cat(sprintf("[Cluster %d = %s | n=%d] %s\n", k, lab[as.character(k)], length(ln),
              paste(ln, collapse = ", ")))
}
cat("\n--- Top Hallmark sets per cluster (mean z) ---\n")
for (k in sort(unique(cl3))) {
  hh <- cl_means %>% filter(cluster == k) %>% slice_max(mean_z, n = 5)
  cat(sprintf("[C%d] %s\n", k, paste(sprintf("%s(%.2f)", hh$hallmark, hh$mean_z), collapse = ", ")))
}
cat("\n--- Cluster theme means (Hallmark z) ---\n"); print(as.data.frame(clabel_by %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
cat("\n--- PROGENy mean pathway activity per cluster (z) ---\n")
pg_means <- as_tibble(pg[hgs_lines, ], rownames = "cell_line") %>% mutate(cluster = cl3[cell_line]) %>%
  group_by(cluster) %>% summarise(across(where(is.numeric), mean), .groups = "drop")
print(as.data.frame(pg_means %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
cat("\n--- Robustness: Hallmark-k3 vs PROGENy-k3 clustering (contingency) ---\n"); print(xtab)
cat("\nAll 15 HGS lines are Mes-Masson/CHUM -> no source-site confound in this analysis.\n")
cat("Related sublines are not independent validation; patient-representative sensitivity is deposited.\n")
cat("\nOutputs: output/hgs_heterogeneity.csv, output/hgs_hallmark_cluster_means.csv\n")
cat("Figures: figs/14_hgs_hallmark_heatmap.pdf, 14_hgs_progeny_heatmap.pdf, 14_hgs_theme_scatter.pdf\n")
message("\n14_hgs_heterogeneity.R complete.")
