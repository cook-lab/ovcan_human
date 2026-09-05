# =============================================================================
# Script: 02_rna_separation.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Does the RNA data recapitulate histological subtypes? PCA / t-SNE /
#          silhouette / Spearman clustering — with an explicit SITE-vs-SUBTYPE
#          confounder check (subtype is partly entangled with source lab).
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   1 (RNA re-processing) — step 2 of 2
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(tidyverse); library(DESeq2); library(matrixStats)
  library(Rtsne); library(cluster); library(ComplexHeatmap); library(circlize)
})

vsd <- readRDS(file.path(OUT, "rna_vst.rds"))
ann <- readRDS(file.path(OUT, "rna_dds.rds")) |> colData() |> as.data.frame()
ann$subtype <- factor(ann$subtype, levels = c("HGS","CC","EC","MC","MMMT","SCCOHT"))
v   <- assay(vsd)
sub_cols  <- setNames(RColorBrewer::brewer.pal(6, "Dark2"), levels(ann$subtype))
site_cols <- setNames(c("#C2410C","#0F172A","#6B7280"),
                      c("Mes-Masson","Huntsman","Huntsman/Vanderhyden"))

# 1. Top variable genes + PCA --------------------------------------------------
rv  <- rowVars(v); top <- head(order(rv, decreasing = TRUE), 2000)
pca <- prcomp(t(v[top, ]), scale. = FALSE, center = TRUE)
pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pcd <- as_tibble(pca$x[, 1:10]) |> mutate(cell_line = rownames(pca$x),
                                          subtype = ann$subtype, site = ann$site)
saveRDS(pca, file.path(OUT, "rna_pca.rds"))

pc_plot <- function(colvar, cols, ttl)
  ggplot(pcd, aes(PC1, PC2, colour = .data[[colvar]])) +
    geom_point(size = 3) + ggrepel::geom_text_repel(aes(label = cell_line), size = 2.3, max.overlaps = 20) +
    scale_colour_manual(values = cols) +
    labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2]),
         colour = colvar, title = ttl) + theme_minimal(base_size = 11)
ggsave(file.path(FIGS,"02_rna_pca_subtype.pdf"), pc_plot("subtype", sub_cols, "RNA PCA by subtype"), width=7.5, height=6)
ggsave(file.path(FIGS,"02_rna_pca_site.pdf"),    pc_plot("site", site_cols, "RNA PCA by source site"), width=7.5, height=6)

# 2. t-SNE (perplexity 5; n=31) ------------------------------------------------
set.seed(SEED)
ts <- Rtsne(pca$x[, 1:10], perplexity = 5, theta = 0.0, pca = FALSE)$Y
tsd <- tibble(t1 = ts[,1], t2 = ts[,2], subtype = ann$subtype, site = ann$site, cell_line = ann$cell_line)
ggsave(file.path(FIGS,"02_rna_tsne_subtype.pdf"),
       ggplot(tsd, aes(t1,t2,colour=subtype)) + geom_point(size=3) +
         ggrepel::geom_text_repel(aes(label=cell_line), size=2.3, max.overlaps=20) +
         scale_colour_manual(values=sub_cols) + labs(title="RNA t-SNE (perplexity 5)") +
         theme_minimal(base_size=11), width=7.5, height=6)

# 3. Silhouette by subtype (Euclidean on PC1-10; annotate n) -------------------
d   <- dist(pca$x[, 1:10])
sil <- silhouette(as.integer(ann$subtype), d)
sil_tab <- tibble(subtype = levels(ann$subtype)[sil[, "cluster"]],
                  width = sil[, "sil_width"]) |>
  group_by(subtype) |> summarise(n = n(), mean_sil = mean(width)) |>
  mutate(subtype = factor(subtype, levels = levels(ann$subtype))) |> arrange(subtype)
# The overall mean was previously only printed; write it too, so script 06 can put
# the RNA and protein silhouettes into one modality-labelled table and so the
# manuscript's overall value has a file to be quoted from.
sil_tab <- bind_rows(sil_tab |> mutate(subtype = as.character(subtype)),
                     tibble(subtype = "ALL", n = nrow(pca$x),
                            mean_sil = mean(sil[, "sil_width"])))
readr::write_csv(sil_tab, file.path(OUT, "rna_silhouette.csv"))

# 4. SITE vs SUBTYPE confounder check (R^2 of each PC ~ factor) ----------------
r2 <- function(pc, f) summary(lm(pcd[[pc]] ~ f))$r.squared
conf <- tibble(PC = paste0("PC", 1:5),
               var_pct   = pv[1:5],
               r2_subtype = sapply(paste0("PC",1:5), r2, f = ann$subtype),
               r2_site    = sapply(paste0("PC",1:5), r2, f = ann$site))
readr::write_csv(conf, file.path(OUT, "rna_pc_confounder.csv"))
# HGS is entirely Mes-Masson, so within-HGS data cannot estimate a centre effect.
# Cross-centre comparisons within CC/EC/MC remain sparse and descriptive.

# 5. Spearman sample-correlation heatmap --------------------------------------
sc <- cor(v[top, ], method = "spearman")
ha <- HeatmapAnnotation(subtype = ann$subtype, site = ann$site,
                        col = list(subtype = sub_cols, site = site_cols))
pdf(file.path(FIGS, "02_rna_spearman_heatmap.pdf"), width = 8, height = 7)
draw(Heatmap(sc, name = "Spearman", top_annotation = ha,
             clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
             col = colorRamp2(quantile(sc, c(.5,.75,1)), viridisLite::mako(3)),
             show_row_names = FALSE, column_names_gp = grid::gpar(fontsize = 7)))
dev.off()

# 6. Report --------------------------------------------------------------------
cat("\n=== PCA variance explained (PC1-5) ===\n"); print(pv[1:5])
cat("\n=== Silhouette by subtype (Euclidean, PC1-10; watch n=2 groups) ===\n"); print(as.data.frame(sil_tab))
cat(sprintf("\nOverall mean silhouette: %.3f\n", mean(sil[, "sil_width"])))
cat("\n=== Confounder check: variance in each PC explained by subtype vs site ===\n")
print(as.data.frame(conf))
cat("\nNote: HGS(15)=all Mes-Masson, MMMT(2)=all Huntsman -> subtype/site partly confounded.\n")
cat("Note: every number here is PER LINE (n=31), not per patient. The patient-representative\n")
cat("  sensitivity analysis (n=28) lives in scripts/21_rna_sensitivity.R; quote both.\n")

# 7. Environment record --------------------------------------------------------
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
write_session_info("02_rna_separation")

message("\n02_rna_separation.R complete. Figures in figs/ ; tables in output/.")
