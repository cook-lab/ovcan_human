# =============================================================================
# Script: 37_supp_rnaprot.R
# Purpose: Expression supplementary figures S1-S4, S7-S8 from audited outputs.
# Design revision: 2026-09-05. All explanatory prose is in external legends.
# Canonical outputs: docs/manuscript/figures/figs{1,2,3,4,7,8}.{pdf,png}
# =============================================================================
source("scripts/00_setup.R")
source("scripts/00b_figure_theme.R")
suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(matrixStats)
  library(ComplexHeatmap); library(circlize); library(grid)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count <- dplyr::count; slice <- dplyr::slice
set.seed(SEED)
# ComplexHeatmap measures strings on an internal PDF device; the final devices
# embed actual Arial. Helvetica is the compatible measurement alias only.
if (FIG_FONT == "Arial" && !"Arial" %in% names(grDevices::pdfFonts()))
  grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])
ht_opt_cook()
SUB_LVL <- c("HGS", "LGS", "CC", "EC", "MC", "MMMT", "SCCOHT")
ASSETS <- file.path(PROJ, "reports", "assets")
DESIGN <- file.path(PROJ, "reports", "figure_redesign_2026-09-05")
dir.create(ASSETS, showWarnings = FALSE, recursive = TRUE)
dir.create(DESIGN, showWarnings = FALSE, recursive = TRUE)
pt_mm <- function(pt) pt / (72.27 / 25.4)
base_theme <- function() theme_ovcan(base_size = 8) + theme(
  axis.text = element_text(size = 7.3, colour = cook_ink),
  axis.title = element_text(size = 8),
  legend.text = element_text(size = 7), legend.title = element_text(size = 7.5),
  legend.position = "bottom", legend.key.size = unit(8, "pt"),
  legend.margin = margin(1, 0, 0, 0), legend.spacing.x = unit(3, "pt"),
  plot.margin = margin(4, 4, 3, 4))
tag_theme <- theme(plot.tag = element_text(size = 11, face = "bold", family = FIG_FONT,
                                          colour = cook_ink))
save_pair <- function(p, stem, height, width = W2) {
  save_fig(p, file.path(MSFIG, paste0(stem, ".pdf")), width, height)
  save_fig(p, file.path(MSFIG, paste0(stem, ".png")), width, height)
}

# S1: matched protein PCA views and modality-specific silhouette summaries.
ppca <- readRDS(file.path(OUT, "prot_pca.rds"))
pann <- read_csv(file.path(OUT, "prot_sample_annotation.csv"), show_col_types = FALSE) %>%
  transmute(cell_line, subtype = factor(subtype, levels = SUB_LVL), plex = factor(plex))
pvar <- 100 * ppca$sdev^2 / sum(ppca$sdev^2)
pscore <- as_tibble(ppca$x[, 1:2], rownames = "cell_line") %>% left_join(pann, by = "cell_line")
stopifnot(nrow(pscore) == 31L)
check_palette_keys(pscore$subtype, subtype_colours, "histotype")
check_palette_keys(pscore$plex, plex_colours, "TMT plex")
cc_drv <- pscore %>% filter(subtype == "CC") %>% slice_max(PC2, n = 1)
pca_base <- function() ggplot(pscore, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = cook_hair, linewidth = .25) +
  geom_vline(xintercept = 0, colour = cook_hair, linewidth = .25) +
  labs(x = sprintf("PC1 (%.1f%%)", pvar[1]), y = sprintf("PC2 (%.1f%%)", pvar[2])) +
  base_theme() + theme(legend.title = element_text(size = 7.2),
    legend.key.width = unit(7, "pt"), legend.spacing.x = unit(2, "pt"))
s1a <- pca_base() +
  geom_point(aes(fill = subtype), shape = 21, colour = cook_ink, size = 2.3, stroke = .25) +
  ggrepel::geom_text_repel(data = cc_drv, aes(label = cell_line), size = pt_mm(7),
    colour = cook_ink, family = FIG_FONT, min.segment.length = 0,
    segment.size = .25, box.padding = .35, seed = SEED) +
  scale_fill_subtype(name = "Histotype", drop = TRUE) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE, title.position = "top"))
s1b <- pca_base() +
  geom_point(aes(fill = plex), shape = 21, colour = cook_ink, size = 2.3, stroke = .25) +
  scale_fill_plex(name = "TMT plex") +
  guides(fill = guide_legend(nrow = 1, title.position = "top"))
silm <- read_csv(file.path(OUT, "silhouette_by_modality.csv"), show_col_types = FALSE) %>%
  mutate(modality = factor(ifelse(modality == "RNA", "RNA", "Protein"), c("RNA", "Protein")),
    lab = ifelse(subtype == "ALL", sprintf("All (%d)", n), sprintf("%s (%d)", subtype, n)))
sil_ord <- silm %>% filter(modality == "Protein") %>% arrange(mean_sil) %>% pull(lab)
silm$lab <- factor(silm$lab, sil_ord)
s1c <- ggplot(silm, aes(mean_sil, lab, colour = modality, shape = modality)) +
  geom_vline(xintercept = 0, colour = cook_hair, linewidth = .35) +
  geom_line(aes(group = lab), colour = cook_hair, linewidth = .5, orientation = "y") +
  geom_point(size = 2) +
  scale_colour_manual(values = c(RNA = cook_slate, Protein = cook_rust), name = NULL) +
  scale_shape_manual(values = c(RNA = 16, Protein = 17), name = NULL) +
  scale_x_continuous(breaks = c(0, .4, .8), expand = expansion(mult = .10)) +
  labs(x = "Mean silhouette width", y = "Histotype (models)") + base_theme() +
  theme(axis.title.y = element_blank())
figs1 <- s1a + s1b + s1c + plot_layout(widths = c(1, 1, 1.05)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_pair(figs1, "figs1", 2.85)
save_fig(s1a, file.path(ASSETS, "f_prot_pca.png"), 3.6, 3.1)

# S2: descriptive nested-model associations and paired assay passage records.
pc <- read_csv(file.path(OUT, "rna_passage_check.csv"), show_col_types = FALSE)
jnt <- read_csv(file.path(OUT, "rna_pc_confounder_joint.csv"), show_col_types = FALSE)
M_LVL <- c("Passage alone (R²)", "After centre (ΔR²)", "After histotype + centre (ΔR²)")
within <- pc %>% filter(check == "PC~passage regression") %>%
  transmute(PC = id, `Passage alone (R²)` = r2_passage,
    `After centre (ΔR²)` = partial_r2_passage_after_site) %>%
  left_join(jnt %>% transmute(PC, `After histotype + centre (ΔR²)` =
    unique_passage_beyond_subtype_site), by = "PC") %>%
  pivot_longer(all_of(M_LVL), names_to = "model", values_to = "r2") %>%
  mutate(model = factor(model, M_LVL))
stopifnot(!anyNA(within$r2))
s2a <- ggplot(within, aes(PC, 100 * r2, fill = model)) +
  geom_col(position = position_dodge(width = .8), width = .73) +
  scale_fill_manual(values = setNames(c(cook_grey, cook_slate, cook_rust), M_LVL), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .06))) +
  labs(x = NULL, y = "PC variance explained (%)") + base_theme() +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.justification = "left", legend.key.height = unit(7, "pt"))
disc <- read_csv(file.path(OUT, "rna_passage_discordance.csv"), show_col_types = FALSE)
d_sum <- function(k) {
  v <- disc$passage_diff[disc$block == "summary" & disc$cell_line == k]
  stopifnot(length(v) == 1); v
}
cross <- pc %>% filter(check == "cross-assay passage mismatch") %>%
  transmute(cell_line = id, rna_passage = as.numeric(rna_passage),
    wes_passage = as.numeric(wes_passage), passage_diff = as.numeric(passage_diff)) %>%
  arrange(rna_passage, cell_line) %>% mutate(cell_line = factor(cell_line, cell_line))
stopifnot(nrow(cross) == d_sum("n_lines_with_both"))
ct_p <- cor.test(cross$rna_passage, cross$wes_passage)
cross_long <- cross %>% pivot_longer(c(rna_passage, wes_passage), names_to = "assay", values_to = "passage") %>%
  mutate(assay = factor(assay, c("rna_passage", "wes_passage"), c("RNA-seq", "WES")))
s2b <- ggplot(cross, aes(y = cell_line)) +
  geom_segment(aes(x = rna_passage, xend = wes_passage, yend = cell_line),
    colour = cook_hair, linewidth = .7) +
  geom_point(data = cross_long, aes(x = passage, colour = assay, shape = assay), size = 2.2) +
  scale_colour_manual(values = c(`RNA-seq` = cook_rust, WES = cook_slate), name = NULL) +
  scale_shape_manual(values = c(`RNA-seq` = 16, WES = 17), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(.06, .07))) +
  labs(x = "Passage number", y = NULL) + base_theme() +
  theme(axis.ticks.y = element_blank(), panel.grid.major.y = element_line(colour = "#F1F5F9", linewidth = .25))
figs2 <- s2a + s2b + plot_layout(widths = c(1, 1.18)) + plot_annotation(tag_levels = "A") & tag_theme
save_pair(figs2, "figs2", 3.2)
save_fig(s2a, file.path(ASSETS, "f_passage_check.png"), 3.6, 3.2)
save_fig(s2b, file.path(ASSETS, "f_passage_check_crossassay.png"), 3.6, 3.2)

# S3: all-model global RNA correlation with source and patient annotations.
vm <- SummarizedExperiment::assay(readRDS(file.path(OUT, "rna_vst.rds")))
cmat <- cor(vm, method = "spearman")
off_diag <- cmat[upper.tri(cmat)]
R_MIN <- min(off_diag); R_MAX <- max(off_diag); R_MED <- median(off_diag)
fam <- read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE) %>%
  transmute(cell_line, subtype = factor(subtype, SUB_LVL), source_site,
    family = ifelse(is_multiline_family, family, NA_character_))
fam <- tibble(cell_line = colnames(cmat)) %>% left_join(fam, by = "cell_line") %>%
  mutate(subtype = droplevels(subtype))
stopifnot(!anyNA(fam$source_site))
multi_fams <- sort(unique(na.omit(fam$family)))
check_palette_keys(multi_fams, family_colours, "patient")
check_palette_keys(fam$source_site, site_colours, "centre")
fam_lab <- replace_na(fam$family, "Other patients")
fam_cols <- c(family_colours[multi_fams], `Other patients` = cook_grey)
SITE_S3 <- c("Mes-Masson" = "CHUM", "Huntsman" = "BC Cancer", "Huntsman/Vanderhyden" = "OHRI")
site_lab <- unname(SITE_S3[fam$source_site])
site_cols <- setNames(site_colours[names(SITE_S3)], unname(SITE_S3))
# Positive correlations use a sequential light-to-rust scale; diagonal remains 1.
col_cor <- circlize::colorRamp2(seq(R_MIN, 1, length.out = 6), cook_sequential[1:6])
lgnd <- function(...) list(title_gp = gpar(fontsize = 7.5, fontfamily = FIG_FONT, col = cook_ink),
  labels_gp = gpar(fontsize = 7, fontfamily = FIG_FONT, col = cook_ink),
  grid_width = unit(3, "mm"), grid_height = unit(3, "mm"), ...)
ra <- rowAnnotation(Centre = site_lab, Histotype = fam$subtype, Patient = fam_lab,
  col = list(Centre = site_cols, Histotype = subtype_colours[levels(fam$subtype)], Patient = fam_cols),
  annotation_name_gp = gpar(fontsize = 7, fontfamily = FIG_FONT, col = cook_ink),
  annotation_legend_param = list(Centre = lgnd(), Histotype = lgnd(), Patient = lgnd()),
  simple_anno_size = unit(3, "mm"), gap = unit(1, "mm"))
ht_s3 <- Heatmap(cmat, name = "Spearman r", col = col_cor, right_annotation = ra,
  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
  row_dend_width = unit(7, "mm"), column_dend_height = unit(7, "mm"),
  row_dend_gp = gpar(col = cook_slate), column_dend_gp = gpar(col = cook_slate),
  row_names_gp = gpar(fontsize = 7, fontfamily = FIG_FONT, col = cook_ink),
  column_names_gp = gpar(fontsize = 7, fontfamily = FIG_FONT, col = cook_ink),
  rect_gp = gpar(col = NA), border_gp = gpar(col = cook_hair, lwd = .6),
  heatmap_legend_param = list(direction = "horizontal", legend_width = unit(43, "mm"),
    at = c(round(R_MIN, 2), .8, .9, 1),
    title_gp = gpar(fontsize = 7.5, fontfamily = FIG_FONT, col = cook_ink),
    labels_gp = gpar(fontsize = 7, fontfamily = FIG_FONT, col = cook_ink)),
  width = ncol(cmat) * unit(3.3, "mm"), height = nrow(cmat) * unit(3.3, "mm"))
draw_s3 <- function() draw(ht_s3, heatmap_legend_side = "bottom", annotation_legend_side = "right",
  merge_legends = FALSE, legend_grouping = "original", padding = unit(c(2, 2, 2, 2), "mm"))
grDevices::cairo_pdf(file.path(MSFIG, "figs3.pdf"), width = W2, height = 5.65)
draw_s3(); dev.off()
ragg::agg_png(file.path(MSFIG, "figs3.png"), width = W2, height = 5.65, units = "in", res = 400)
draw_s3(); dev.off()

# S4: whole-plex presence patterns and counts; every feature appears once.
bm <- read_csv(file.path(OUT, "prot_block_missingness.csv"), show_col_types = FALSE)
fa <- read_csv(file.path(OUT, "prot_feature_accounting.csv"), show_col_types = FALSE)
fa_n <- function(lvl) { v <- fa$n[fa$level == lvl]; stopifnot(length(v) == 1); as.integer(v) }
N_MATRIX <- fa_n("analysis rows (matrix)")
N_ZEROPLEX <- fa_n("quantified in 0 lines (zero-plex)")
N_ABSENT1 <- fa_n("absent from >=1 whole plex")
pat <- bm %>% mutate(plexes_present = replace_na(plexes_present, "")) %>%
  count(plexes_present, present_n_plex, name = "n_proteins") %>%
  arrange(desc(present_n_plex), desc(n_proteins)) %>%
  mutate(pattern_id = row_number(), is_zero = present_n_plex == 0)
stopifnot(sum(pat$n_proteins) == N_MATRIX, any(pat$is_zero), pat$n_proteins[pat$is_zero] == N_ZEROPLEX)
pid_lv <- rev(pat$pattern_id)
grid_long <- pat %>% rowwise() %>%
  mutate(present = list(if (plexes_present == "") integer(0) else as.integer(strsplit(plexes_present, ";")[[1]]))) %>%
  ungroup() %>% tidyr::crossing(plex = 1:5) %>%
  mutate(state = ifelse(map2_lgl(present, plex, ~ .y %in% .x), "Quantified", "Not quantified"),
    pattern_id = factor(pattern_id, pid_lv))
pat_ord <- pat %>% mutate(pattern_id = factor(pattern_id, pid_lv))
s4a <- ggplot(grid_long, aes(factor(plex), pattern_id, fill = state)) +
  geom_tile(colour = "white", linewidth = .5) +
  scale_fill_manual(values = c(Quantified = cook_slate, `Not quantified` = cook_grey),
    name = NULL, breaks = c("Quantified", "Not quantified")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = expansion(add = .5), breaks = NULL) +
  labs(x = "TMT plex", y = NULL) + base_theme() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
s4b <- ggplot(pat_ord, aes(n_proteins, pattern_id)) +
  geom_hline(yintercept = seq_len(nrow(pat)), colour = "#F1F5F9", linewidth = .25) +
  geom_point(aes(colour = is_zero), size = 1.8) +
  geom_text(aes(label = format(n_proteins, big.mark = ",", trim = TRUE)),
    hjust = -.35, size = pt_mm(7.2), colour = cook_ink, family = FIG_FONT) +
  scale_colour_manual(values = c(`TRUE` = cook_rust, `FALSE` = cook_slate), guide = "none") +
  scale_x_log10(limits = c(10, 18000), breaks = c(10, 100, 1000, 10000),
    labels = scales::label_comma(), expand = expansion(mult = c(.02, .02))) +
  scale_y_discrete(expand = expansion(add = .5), breaks = NULL) +
  labs(x = "Protein features (log scale)", y = NULL) + base_theme() +
  theme(axis.ticks.y = element_blank(), panel.grid.major.x = element_line(colour = cook_hair, linewidth = .25))
figs4 <- s4a + s4b + plot_layout(widths = c(.78, 1.35), guides = "collect") +
  plot_annotation(tag_levels = "A") & tag_theme & theme(legend.position = "top")
save_pair(figs4, "figs4", 4.65, width = W15)

# S7: uncertainty for all 25 marker checks, using patient representatives.
me <- read_csv(file.path(OUT, "rna_marker_effectsizes.csv"), show_col_types = FALSE) %>%
  filter(!is.na(cohens_d)) %>%
  mutate(grp = factor(marks_code, c("HGS", "CC", "EC", "MC", "MMMT", "SCCOHT")),
    subtype = factor(marks_code, SUB_LVL), bh_sig = wilcox_p_bh < .05,
    floor_lim = wilcox_p_floor_limited) %>%
  arrange(grp, desc(cohens_d)) %>% mutate(sym_f = factor(symbol, rev(symbol)))
stopifnot(nrow(me) == 25, !anyNA(me$grp), length(unique(me$n_boot_valid)) == 1)
N_BOOT <- unique(me$n_boot_valid); N_UNITS <- unique(me$n_intended + me$n_rest)
stopifnot(length(N_UNITS) == 1)
y_lab <- setNames(sprintf("%s%s%s%s", me$symbol, ifelse(me$direction == "down", " ↓", ""),
  ifelse(me$bh_sig, " *", ""), ifelse(me$floor_lim, " †", "")), me$symbol)
grp_pos <- me %>% mutate(idx = as.integer(sym_f)) %>% group_by(grp) %>%
  summarise(ylo = min(idx) - .5, yhi = max(idx) + .5, ymid = mean(idx),
    n_int = first(n_intended), .groups = "drop")
sep7 <- head(sort(grp_pos$yhi), -1)
D_MIN <- -4; D_MAX <- 6
me <- me %>% mutate(d_hi_clip = pmin(cohens_d_hi, D_MAX), d_lo_clip = pmax(cohens_d_lo, D_MIN),
  d_over = cohens_d_hi > D_MAX)
mk_panel <- function(xvar, lo, hi, ref, xlab) ggplot(me, aes(.data[[xvar]], sym_f)) +
  geom_vline(xintercept = ref, colour = cook_hair, linewidth = .45) +
  geom_hline(yintercept = sep7, colour = cook_hair, linewidth = .3) +
  geom_linerange(aes(xmin = .data[[lo]], xmax = .data[[hi]], colour = subtype), linewidth = .55) +
  geom_point(aes(fill = subtype, shape = lands_right), colour = cook_ink, size = 2.1, stroke = .25) +
  scale_fill_subtype(guide = "none", drop = TRUE) + scale_colour_subtype(guide = "none", drop = TRUE) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), name = "Recovery rule",
    breaks = c("TRUE", "FALSE"), labels = c("Met", "Not met")) +
  guides(shape = guide_legend(nrow = 1, override.aes = list(fill = cook_slate, size = 2.1))) +
  labs(x = xlab, y = NULL) + base_theme() + theme(axis.ticks.y = element_blank())
s7a <- mk_panel("cohens_d", "d_lo_clip", "d_hi_clip", 0, "Signed Cohen's d") +
  geom_text(data = ~ filter(.x, d_over), aes(x = D_MAX, label = sprintf("→ %.1f", cohens_d_hi)),
    hjust = -.05, size = pt_mm(7), colour = cook_ink, family = FIG_FONT) +
  scale_x_continuous(limits = c(D_MIN, D_MAX), breaks = seq(-4, 6, 2), expand = c(0, 0)) +
  scale_y_discrete(labels = y_lab, expand = expansion(add = .6)) +
  coord_cartesian(clip = "off") + theme(plot.margin = margin(4, 34, 3, 4))
s7b <- mk_panel("auc_oriented", "auc_oriented_lo", "auc_oriented_hi", .5, "Oriented AUC") +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1), expand = c(0, 0)) +
  scale_y_discrete(labels = NULL, breaks = NULL, expand = expansion(add = .6)) +
  geom_segment(data = grp_pos, aes(x = 1.035, xend = 1.035, y = ylo, yend = yhi),
    inherit.aes = FALSE, colour = cook_hair, linewidth = .35) +
  geom_text(data = grp_pos, aes(x = 1.06, y = ymid, label = sprintf("%s\n(n = %d)", grp, n_int)),
    inherit.aes = FALSE, hjust = 0, size = pt_mm(7.2), lineheight = 1,
    colour = cook_ink, family = FIG_FONT) +
  coord_cartesian(xlim = c(0, 1), clip = "off") + theme(plot.margin = margin(4, 51, 3, 4))
figs7 <- s7a + s7b + plot_layout(widths = c(1.12, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") & tag_theme & theme(legend.position = "bottom")
save_pair(figs7, "figs7", 5.1)

# S8: illustrative HGS partitions; separate signed proliferation panel.
het <- read_csv(file.path(OUT, "hgs_heterogeneity.csv"), show_col_types = FALSE)
k3_sens <- read_csv(file.path(OUT, "hgs_cluster_patient_sensitivity.csv"), show_col_types = FALSE) %>% filter(k == 3L)
stopifnot(nrow(k3_sens) == 1)
STRAT_KEY <- c("Inflammatory/NF-kB-EMT" = "Inflammatory", "Hypoxic-glycolytic" = "Hypoxic-glycolytic",
  "Low-signaling" = "Low-signaling")
clab_cols <- setNames(stratum_colours[unname(STRAT_KEY)], names(STRAT_KEY))
check_palette_keys(het$cluster_label, clab_cols, "illustrative HGS group")
het <- het %>% mutate(cluster_label = factor(cluster_label, names(STRAT_KEY)))
strat_n <- het %>% count(cluster_label)
strat_labels <- setNames(sprintf("%s (n = %d)", c("Inflammatory", "Hypoxic-glycolytic", "Low-signalling"), strat_n$n),
  as.character(strat_n$cluster_label))
s8a <- ggplot(het, aes(inflammatory_z, hypoxic_glycolytic_z)) +
  geom_hline(yintercept = 0, colour = cook_hair, linewidth = .3) +
  geom_vline(xintercept = 0, colour = cook_hair, linewidth = .3) +
  geom_point(aes(fill = cluster_label), shape = 21, size = 2.5, colour = cook_ink, stroke = .3) +
  ggrepel::geom_text_repel(aes(label = cell_line), size = pt_mm(7), colour = cook_ink, family = FIG_FONT,
    max.overlaps = Inf, segment.size = .25, segment.colour = cook_slate, min.segment.length = 0,
    box.padding = .55, point.padding = .25, force = 2, max.iter = 10000,
    nudge_x = ifelse(het$cell_line == "TOV3041G", .45, 0),
    nudge_y = ifelse(het$cell_line == "TOV3041G", .25, 0), seed = SEED) +
  scale_fill_manual(values = clab_cols, name = "Illustrative group", labels = strat_labels) +
  scale_x_continuous(expand = expansion(mult = .18)) + scale_y_continuous(expand = expansion(mult = .18)) +
  labs(x = "Inflammatory theme (mean z)", y = "Hypoxic-glycolytic theme (mean z)") + base_theme() +
  guides(fill = guide_legend(nrow = 1, title.position = "top"))
prolif_ord <- het %>% arrange(proliferation_z) %>% pull(cell_line)
s8b <- ggplot(het %>% mutate(cell_line = factor(cell_line, prolif_ord)), aes(proliferation_z, cell_line)) +
  geom_vline(xintercept = 0, colour = cook_hair, linewidth = .4) +
  geom_segment(aes(x = 0, xend = proliferation_z, yend = cell_line, colour = cluster_label), linewidth = .6) +
  geom_point(aes(fill = cluster_label), shape = 21, size = 2.1, colour = cook_ink, stroke = .25) +
  scale_fill_manual(values = clab_cols, guide = "none") + scale_colour_manual(values = clab_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = .12)) +
  labs(x = "Proliferation theme (mean z)", y = NULL) + base_theme() + theme(axis.ticks.y = element_blank())
prog <- het %>% select(cell_line, cluster, cluster_label, starts_with("progeny_"))
pm <- prog %>% select(cell_line, starts_with("progeny_")) %>% column_to_rownames("cell_line") %>% as.matrix() %>% t()
rownames(pm) <- sub("^progeny_", "", rownames(pm))
ord <- prog %>% arrange(cluster, cell_line)
pm <- pm[, ord$cell_line]
path_order <- rownames(pm)[hclust(dist(pm))$order]
short_lab <- c("Inflammatory/NF-kB-EMT" = "Inflam.", "Hypoxic-glycolytic" = "Hypoxic", "Low-signaling" = "Low-sig.")
prog_long <- as_tibble(pm, rownames = "pathway") %>% pivot_longer(-pathway, names_to = "cell_line", values_to = "z") %>%
  left_join(prog %>% select(cell_line, cluster_label), by = "cell_line") %>%
  mutate(pathway = factor(pathway, path_order), cell_line = factor(cell_line, ord$cell_line),
    stratum = factor(short_lab[as.character(cluster_label)], short_lab[names(clab_cols)]))
rho_hyp <- suppressWarnings(cor.test(het$hypoxic_glycolytic_z, het$progeny_Hypoxia, method = "spearman"))
rho_nfkb <- suppressWarnings(cor.test(het$inflammatory_z, het$progeny_NFkB, method = "spearman"))
zc <- max(abs(prog_long$z), na.rm = TRUE)
s8c <- ggplot(prog_long, aes(cell_line, pathway, fill = z)) +
  geom_tile(colour = "white", linewidth = .4) + facet_grid(~stratum, scales = "free_x", space = "free_x") +
  scale_fill_cook_div(midpoint = 0, limits = c(-zc, zc), name = "PROGENy z") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) +
  base_theme() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 7.2), axis.text.y = element_text(size = 7.3),
    panel.spacing = unit(2, "pt"), strip.text = element_text(size = 7.3), legend.position = "right",
    legend.key.width = unit(8, "pt"), legend.key.height = unit(12, "pt"))
# Collect the categorical key below A/B only; C retains its continuous key.
row8 <- s8a + s8b + plot_layout(widths = c(1.4, 1), guides = "collect") & theme(legend.position = "bottom")
figs8 <- row8 / s8c + plot_layout(heights = c(1.38, 1)) + plot_annotation(tag_levels = "A") & tag_theme
save_pair(figs8, "figs8", 5.8)
save_fig(s8a, file.path(ASSETS, "f_hgs_het.png"), 4.6, 3.6)

# Reproducible external legends: narrative belongs beside, not inside, the PDFs.
legends <- c(
  "# Expression supplementary figure legends", "", "Generated by scripts/37_supp_rnaprot.R from the audited analysis outputs. Revised 5 September 2026.", "",
  "## Supplementary Figure S1. Proteomic ordination and histotype separation", "",
  sprintf("**A, B.** Principal components of the 2,000 most variable complete-case protein features across %d models, using the supplied normalized log2 abundances, centred without feature scaling. The same coordinates are coloured by annotated histotype (A) and TMT plex (B). PC1 and PC2 explain %.1f%% and %.1f%% of variance, respectively. The clear-cell model %s is labelled at the PC2 extreme. Visual overlap does not establish the absence of batch effects; plex and contributing centre should be considered in comparisons.", nrow(pscore), pvar[1], pvar[2], cc_drv$cell_line), "",
  "**C.** Mean silhouette widths for the annotated histotypes in RNA and protein data, using Euclidean distances on the first ten principal components of each modality. Lines connect summaries for the same histotype. Values near zero indicate overlapping groups; negative values indicate greater similarity to another group. Parentheses give model counts, and All denotes the overall model-weighted summary. These descriptive model-level estimates include related sublines; they are not independent-patient estimates. Histotype abbreviations: HGS, high-grade serous; CC, clear cell; EC, historical endometrioid; MC, mucinous; MMMT, malignant mixed Müllerian tumour; SCCOHT, small cell carcinoma of the ovary, hypercalcaemic type.", "",
  "## Supplementary Figure S2. Passage associations and cross-assay passage differences", "",
  "**A.** Associations between recorded passage number and the first five RNA principal components in 31 models, including related sublines. Light grey bars show R² for passage alone; slate and rust bars show the increment in total R² after adding passage to models containing centre, or histotype and centre, respectively. The conditional measures are ΔR², not partial R² rescaled by residual variance. Passage, centre and histotype are partly aligned; these descriptive associations do not estimate a causal effect of passage.", "",
  sprintf("**B.** Paired RNA-seq and WES passage records for the %d models with both values, ordered by RNA-seq passage. Connecting segments identify the two records from the same model; other models lack one or both values. WES minus RNA-seq passage ranges from %+d to %+d, with a median absolute difference of %d passages. The descriptive Pearson correlation is %.2f. Protein passage records are not represented.", nrow(cross), d_sum("min_diff"), d_sum("max_diff"), d_sum("median_abs_diff"), ct_p$estimate), "",
  "## Supplementary Figure S3. Global RNA similarity with sample provenance", "",
  sprintf("Spearman correlations between %d RNA profiles across %s variance-stabilized genes. Off-diagonal correlations range from %.3f to %.3f (median %.3f); the diagonal is 1 by construction. The continuous colour scale spans the observed minimum to 1. Rows and columns use Ward.D2 clustering on Euclidean distances between correlation profiles. Annotation strips identify contributing centre, annotated histotype and patient identifiers for models belonging to multi-model patient groups in the resource. Grey denotes the other distinct patients; the numbered patient colours are shared with the resource figures. Related sublines remain in this descriptive display. High global transcriptomic correlations alone do not authenticate stocks or separate biological from centre effects. CHUM, Centre hospitalier de l'Université de Montréal; OHRI, Ottawa Hospital Research Institute; histotype abbreviations are defined in Supplementary Figure S1.", ncol(vm), format(nrow(vm), big.mark = ","), R_MIN, R_MAX, R_MED), "",
  "## Supplementary Figure S4. Whole-plex protein quantification patterns", "",
  sprintf("**A.** Quantification patterns across five TMT plexes for all %s protein-feature rows. Each row is one of %d observed patterns, ordered first by decreasing number of plexes with quantification and then by decreasing feature count. Slate indicates quantification in model channels within that plex; light grey indicates no quantification. **B.** Corresponding feature counts in the same row order. Points encode count by position on a logarithmic axis, with exact counts printed once. The rust point in the final row marks %d features identified in search output but unquantified in model channels in all five plexes. These are included among %s features lacking quantification in at least one complete plex. The pattern quantified in all five plexes contains %s features. Missingness is structured by plex and does not establish absent or low protein expression; protein-feature rows need not correspond one-to-one with gene symbols.", format(N_MATRIX, big.mark = ","), nrow(pat), N_ZEROPLEX, format(N_ABSENT1, big.mark = ","), format(pat$n_proteins[pat$present_n_plex == 5], big.mark = ",")), "",
  "## Supplementary Figure S7. Patient-level marker effect sizes and uncertainty", "",
  sprintf("Twenty-five prespecified marker checks across %d patient representatives, using symbol-summed log2(TPM + 1) RNA abundance. MKI67 is an additional display control in the main figure and is excluded here. **A.** Signed Cohen's d for the intended histotype versus the remaining representatives. Negative effects are expected for markers annotated for low abundance (↓). **B.** AUC oriented to the expected direction: values above 0.5 indicate higher abundance for high markers or lower abundance for low markers in the intended histotype. Intervals are 95%% percentile intervals from %s stratified bootstrap resamples. Large intervals reflect small groups; intervals for groups of two or three patients are poorly calibrated. Upper Cohen's d interval bounds beyond the displayed maximum of 6 are shown by arrow labels with their actual values; no point estimates are omitted. Brackets give the number of intended-histotype patient representatives.", N_UNITS, format(N_BOOT, big.mark = ",")), "",
  sprintf("Circles mark recovery under the prespecified rank rule applied to histotype means in these same patient representatives: the intended histotype must rank among the top two for high markers, with mean log2(TPM + 1) above 1, or among the bottom two for low markers. Triangles mark checks outside that rule. An asterisk identifies BH-adjusted Wilcoxon q < 0.05 across the 25 markers (%d markers: %s). A dagger identifies a rank-sum comparison constrained by the minimum attainable p value (%s). Individual marker comparisons are not adjusted for centre; these intervals and tests should not be interpreted as isolating histotype from centre effects. The panel-wide centre-restricted recovery permutation analysis is a separate analysis described in the main text. Histotype abbreviations are defined in Supplementary Figure S1.", sum(me$bh_sig), paste(me$symbol[me$bh_sig], collapse = ", "), paste(me$symbol[me$floor_lim], collapse = ", ")), "",
  "## Supplementary Figure S8. Illustrative within-HGS expression patterns", "",
  sprintf("**A.** Inflammatory and hypoxic-glycolytic theme scores for %d HGS models from %d patients. Colours identify an illustrative k = 3 Ward.D2 partition of all 50 Hallmark gene-set scores standardized within HGS. The two displayed theme means were used to name the groups; they do not independently validate the partition. All models are labelled. Repeating clustering with one representative per patient changes the assignments (adjusted Rand index %.3f on shared representatives). The groups therefore should not be treated as stable molecular subtypes. **B.** Signed proliferation theme scores for the same models, ordered by value; zero denotes the mean within HGS. Colour identifies the group from A. Proliferation was not used to name the groups.", nrow(het), k3_sens$n_patients, k3_sens$adjusted_rand_on_shared_representatives), "",
  sprintf("**C.** PROGENy scores for %d pathways standardized within HGS. Columns are grouped by the illustrative partition and ordered by model name within groups; rows use hierarchical clustering with Euclidean distance and complete linkage. The abbreviated group labels correspond to the full key below A and B. Descriptive Spearman correlations between the Hallmark theme means and corresponding PROGENy scores are %.2f for Hypoxia and %.2f for NF-κB. Both scoring systems use the same RNA profiles, so agreement is not independent assay validation. The colour scale is symmetric about zero and covers the observed range. Related models and centre composition remain part of this exploratory display.", nrow(pm), rho_hyp$estimate, rho_nfkb$estimate), "")
writeLines(legends, file.path(DESIGN, "legends_supp_expression.md"), useBytes = TRUE)
message("Expression supplements and external legends complete: S1-S4, S7-S8.")
