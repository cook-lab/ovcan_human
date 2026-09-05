# Figure3: expression structure, descriptive concordance and marker checks.
# Redesign 2026-09-05. Reads existing results; statistical outputs are unchanged.
# External legend: reports/figure_redesign_2026-09-05/legends_expression.md
source("scripts/00_setup.R"); source("scripts/00b_figure_theme.R")
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count <- dplyr::count; slice <- dplyr::slice
set.seed(SEED)
ASSETS <- file.path(PROJ, "reports", "assets")
SUB_LVL <- c("HGS", "LGS", "CC", "EC", "MC", "MMMT", "SCCOHT")
fig_theme <- function() theme_ovcan(base_size = 8.3) + theme(
  plot.title = element_text(size = 8.5, face = "bold", margin = margin(b = 5)),
  plot.tag = element_text(size = 11, face = "bold", family = FIG_FONT),
  plot.tag.position = "topleft", axis.text = element_text(size = 7.2),
  axis.title = element_text(size = 8), legend.text = element_text(size = 7),
  legend.title = element_text(size = 7), legend.key.size = unit(8, "pt"),
  legend.margin = margin(0, 0, 0, 0), plot.margin = margin(5, 5, 5, 5))
SITE_LVL <- c("Mes-Masson", "Huntsman", "Huntsman/Vanderhyden")
SITE_DISP <- c("Mes-Masson" = "CHUM", "Huntsman" = "BC Cancer", "Huntsman/Vanderhyden" = "OHRI")
site_shapes_disp <- setNames(site_shapes[SITE_LVL], SITE_DISP[SITE_LVL])
pca <- readRDS(file.path(OUT, "rna_pca.rds"))
ann <- read_csv(file.path(OUT, "rna_sample_annotation.csv"), show_col_types = FALSE) %>%
  transmute(cell_line, subtype = factor(subtype, levels = SUB_LVL),
            site_f = factor(SITE_DISP[site], levels = SITE_DISP[SITE_LVL]))
var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
scores <- as_tibble(pca$x[, 1:2], rownames = "cell_line") %>% left_join(ann, by = "cell_line")
stopifnot(nrow(scores) == 31L, !anyNA(scores$site_f))
ec_lines <- scores %>% filter(subtype == "EC")
pA <- ggplot(scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = cook_hair, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = cook_hair, linewidth = 0.25) +
  geom_point(aes(fill = subtype, shape = site_f), colour = cook_ink, size = 2.6, stroke = 0.35) +
  ggrepel::geom_text_repel(data = ec_lines, aes(label = cell_line), size = 2.55,
    family = FIG_FONT, colour = cook_ink, min.segment.length = 0,
    segment.size = 0.25, segment.colour = cook_ink_muted, box.padding = 0.5, seed = SEED) +
  scale_fill_subtype(name = "Histotype", drop = TRUE) +
  scale_shape_manual(values = site_shapes_disp, name = "Centre") +
  scale_x_continuous(expand = expansion(mult = 0.08), n.breaks = 4) +
  scale_y_continuous(expand = expansion(mult = 0.08), n.breaks = 4) +
  labs(title = "RNA expression (31 models)", tag = "A",
       x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
  fig_theme() + theme(legend.position = "bottom", legend.box = "vertical", legend.justification = "left",
    legend.spacing.y = unit(1, "pt")) +
  guides(fill = guide_legend(order = 1, nrow = 2, byrow = TRUE, override.aes = list(shape = 21, size = 2.5)),
         shape = guide_legend(order = 2, nrow = 1, override.aes = list(fill = cook_grey, size = 2.5)))

# B: show every commonality component on a common scale, including the small centre increments.
joint <- read_csv(file.path(OUT, "rna_pc_confounder_joint.csv"), show_col_types = FALSE) %>% filter(PC == "PC1")
rep_pc <- read_csv(file.path(OUT, "sensitivity_patient_reps_pca.csv"), show_col_types = FALSE) %>%
  filter(PC == "PC1", set == "28 reps (HVG refit)")
stopifnot(nrow(joint) == 1L, nrow(rep_pc) == 1L)
commonality <- bind_rows(
  joint %>% transmute(analysis = "31 models", unique_subtype, unique_site, shared, unexplained = 1 - r2_joint),
  rep_pc %>% transmute(analysis = "28 patients", unique_subtype, unique_site, shared, unexplained = 1 - r2_joint)) %>%
  pivot_longer(-analysis, names_to = "component", values_to = "fraction") %>%
  mutate(analysis = factor(analysis, levels = c("31 models", "28 patients")),
         component = factor(component, levels = c("unexplained", "shared", "unique_site", "unique_subtype")),
         y = as.numeric(component) + ifelse(analysis == "31 models", 0.23, -0.23))
stopifnot(all(commonality$fraction >= 0), all(abs(tapply(commonality$fraction, commonality$analysis, sum) - 1) < 1e-8))
analysis_cols <- c(`31 models` = cook_rust, `28 patients` = cook_slate)
pB <- ggplot(commonality, aes(fraction * 100, as.numeric(component), colour = analysis)) +
  geom_vline(xintercept = 0, colour = cook_hair, linewidth = 0.25) +
  geom_point(aes(shape = analysis), size = 1.9) +
  geom_text(aes(label = sprintf("%.1f", 100 * fraction)), nudge_x = 5.4, hjust = 0,
            size = 2.45, colour = cook_ink) +
  scale_colour_manual(values = analysis_cols, guide = "none") +
  scale_shape_manual(values = c(`31 models` = 16, `28 patients` = 17), guide = "none") +
  scale_x_continuous(limits = c(-2, 82), breaks = c(0, 30, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 4.5), breaks = 1:4,
    labels = c("Unexplained", "Shared", "Unique centre", "Unique histotype"), expand = c(0, 0)) +
  facet_wrap(~ analysis, nrow = 1) +
  labs(title = "PC1 commonality", tag = "B", x = "PC1 variance (%)", y = NULL) + fig_theme() +
  theme(axis.text.y = element_text(size = 7), axis.ticks.y = element_blank(),
        strip.text = element_text(size = 7.2, face = "plain"), panel.spacing.x = unit(7, "pt"))

# C: descriptive gene-wise concordance across the 30 paired models.
pg <- read_csv(file.path(OUT, "integ_rnaprot_cor.csv"), show_col_types = FALSE) %>% filter(level == "per_gene")
cs <- read_csv(file.path(OUT, "integ_rnaprot_cor_summary.csv"), show_col_types = FALSE) %>% filter(metric == "per_gene_spearman")
stopifnot(nrow(cs) == 1L, nrow(pg) == cs$n_pergene_reported)
pC <- ggplot(pg, aes(spearman)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, boundary = 0,
    fill = cook_grey, colour = "white", linewidth = 0.15) +
  geom_density(colour = cook_slate, linewidth = 0.5) +
  geom_vline(xintercept = cs$median, colour = cook_rust, linewidth = 0.6) +
  annotate("text", x = -0.94, y = Inf, vjust = 1.15, hjust = 0, size = 2.5,
    colour = cook_ink, label = sprintf("Median %.3f", cs$median)) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1), expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0, 0.08))) +
  labs(title = sprintf("RNA-protein (%s genes)", format(nrow(pg), big.mark = ",")),
    tag = "C", x = "Per-gene Spearman r", y = "Density") + fig_theme()

ms <- read_csv(file.path(OUT, "rna_markers_summary.csv"), show_col_types = FALSE)
SUB6 <- c("HGS","CC","EC","MC","MMMT","SCCOHT")
sub_mean_cols <- paste0(SUB6, "_mean")
mmat <- as.matrix(ms[, sub_mean_cols]); rownames(mmat) <- ms$symbol
colnames(mmat) <- SUB6
# Match the inferential unit in the primary marker display. The individual-model
# panel remains available to inspect within-patient and between-model variation.
model_mmat <- mmat
tpm_f <- read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g_f <- read_csv(file.path(OUT, "tx2gene_matched.csv"), show_col_types = FALSE) %>%
  distinct(ensembl_gene_id, external_gene_name) %>%
  filter(!is.na(external_gene_name), external_gene_name != "")
m_id  <- as.matrix(tpm_f[, -1]); rownames(m_id) <- tpm_f[[1]]
i2s   <- setNames(t2g_f$external_gene_name, t2g_f$ensembl_gene_id)
common <- intersect(rownames(m_id), names(i2s))
logtpm <- log2(rowsum(m_id[common, , drop = FALSE], group = i2s[common]) + 1)
ann_sub <- ann %>% mutate(subtype = factor(subtype, levels = SUB6))
mk_mat  <- logtpm[ms$symbol, ann_sub$cell_line, drop = FALSE]
sm_chk  <- t(vapply(ms$symbol, function(g)
  tapply(mk_mat[g, ], ann_sub$subtype, mean)[SUB6], numeric(length(SUB6))))
stopifnot("recomputed subtype means must equal rna_markers_summary.csv (re-render, not recompute)" =
            max(abs(round(sm_chk, 2) - model_mmat)) < 1e-8)
message("Individual-model marker matrix agrees with the deposited subtype means.")

fmap <- read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)
rep_ids <- fmap$cell_line[fmap$patient_representative & fmap$has_rna]
rep_ann <- ann_sub %>% filter(cell_line %in% rep_ids)
stopifnot(nrow(rep_ann) == 28L)
mmat <- t(vapply(ms$symbol, function(g)
  tapply(mk_mat[g, rep_ann$cell_line], rep_ann$subtype, mean)[SUB6], numeric(length(SUB6))))
rownames(mmat) <- ms$symbol; colnames(mmat) <- SUB6
marker_effects <- read_csv(file.path(OUT, "rna_marker_effectsizes.csv"), show_col_types = FALSE) %>%
  filter(!is.na(marks_code))
expected_means <- mmat[cbind(match(marker_effects$symbol, rownames(mmat)),
                            match(marker_effects$marks_code, colnames(mmat)))]
stopifnot("Patient marker means must agree with the effect-size analysis" =
            max(abs(expected_means - marker_effects$mean_intended)) < 0.00051)
zmat <- t(scale(t(mmat)))                                   # z across the 6 subtype means, per marker
N_MARK <- nrow(ms)
# intended-subtype group per marker. MKI67 has no marks_code: it is a proliferation
# CONTROL, not an intended subtype marker, so it has NO expected cell and no box.
grp <- ifelse(is.na(ms$marks_code) | ms$marks_code == "NA", "Control", ms$marks_code)
sym_order <- ms$symbol                       # HGS..CC..EC..MC..MMMT..SCCOHT..MKI67
y_lvls <- rev(sym_order)                     # bottom -> top on a discrete y
# loss markers get a direction glyph on the axis label: without it, a DARK (low)
# boxed SCCOHT cell reads as a failed marker when it is the expected result.
dir_down  <- ms$direction == "down"
y_labels  <- setNames(ifelse(dir_down, paste0(ms$symbol, " ↓"), ms$symbol), ms$symbol)
y_labels[ms$symbol == "MKI67"] <- "MKI67 (control)"

mlong <- as_tibble(zmat, rownames = "symbol") %>%
  pivot_longer(-symbol, names_to = "subtype", values_to = "z") %>%
  mutate(subtype = factor(subtype, levels = SUB6),
         symbol  = factor(symbol, levels = y_lvls))
# expected cell = intended subtype column of each marker. Solid box = expected HIGH
# (up marker), dashed box = expected LOW (loss marker). The Control row has none.
expected_cell <- tibble(symbol = ms$symbol, grp_code = grp, direction = ms$direction) %>%
  filter(grp_code != "Control") %>%
  mutate(exp_sub = factor(grp_code, levels = SUB6),
         symbol  = factor(symbol, levels = y_lvls),
         expect  = ifelse(direction == "down", "expected LOW (loss)", "expected HIGH"))
stopifnot("every non-control marker needs an expected cell" =
            nrow(expected_cell) == sum(grp != "Control"))
# group boundaries + centres on the numeric y (1 = bottom); groups in sym_order top->bottom
grp_run <- rle(grp[match(sym_order, ms$symbol)])            # groups in table order
top_pos <- length(sym_order)                                # PAX8 sits at the top
bounds <- head(cumsum(grp_run$lengths), -1)                 # between-group indices (from top)
sep_y  <- top_pos - bounds + 0.5                            # y positions of separators
ctr_from_top <- cumsum(grp_run$lengths) - grp_run$lengths / 2 + 0.5
grp_ctr <- tibble(grp = grp_run$values, y = top_pos - ctr_from_top + 1) %>%
  mutate(n_mk = grp_run$lengths,
         lab  = grp)
EC_ROWS <- grp_ctr %>% filter(grp == "EC")

# D/E: two complementary views share marker row order and boundaries.
# Patient means for summary inference; all models for within-panel variation.
zmax <- max(abs(mlong$z), na.rm = TRUE)
patient_n <- rep_ann %>% count(subtype)
patient_labels <- setNames(sprintf("%s\n%d", patient_n$subtype, patient_n$n), patient_n$subtype)
pD <- ggplot(mlong, aes(subtype, symbol, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = expected_cell, aes(x = exp_sub, y = symbol, linetype = expect),
    inherit.aes = FALSE, fill = NA, colour = cook_ink, linewidth = 0.4, show.legend = FALSE) +
  scale_linetype_manual(values = c(`expected HIGH` = "solid", `expected LOW (loss)` = "22"), guide = "none") +
  geom_hline(yintercept = sep_y, colour = cook_hair, linewidth = 0.6) +
  scale_fill_cook_div(midpoint = 0, limits = c(-zmax, zmax), name = "Mean z",
    guide = guide_colourbar(direction = "horizontal", title.position = "left", barwidth = unit(38, "mm"),
                           barheight = unit(2.7, "mm"))) +
  scale_x_discrete(expand = c(0, 0), labels = patient_labels) +
  scale_y_discrete(expand = c(0, 0), labels = y_labels) +
  labs(title = "Patient means (28 patients)", tag = "D", x = "Histotype / patients", y = NULL) +
  fig_theme() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7, lineheight = 0.92), axis.text.y = element_text(size = 7.1),
    legend.position = "bottom", legend.justification = "left", legend.margin = margin(1, 0, 0, 0))

per_line <- as_tibble(mk_mat, rownames = "symbol") %>%
  pivot_longer(-symbol, names_to = "cell_line", values_to = "log2_tpm") %>%
  left_join(ann_sub %>% select(cell_line, subtype), by = "cell_line") %>%
  left_join(tibble(symbol = ms$symbol, intended = grp), by = "symbol") %>%
  mutate(is_intended = as.character(subtype) == intended, symbol = factor(symbol, levels = y_lvls))
pE <- ggplot(per_line, aes(log2_tpm, symbol)) +
  geom_vline(xintercept = 1, colour = cook_hair, linewidth = 0.45, linetype = "dashed") +
  geom_hline(yintercept = sep_y, colour = cook_hair, linewidth = 0.6) +
  geom_point(data = ~ filter(.x, !is_intended), colour = "#CBD5E1", size = 0.8, alpha = 0.9,
    position = position_jitter(width = 0, height = 0.13, seed = SEED)) +
  geom_point(data = ~ filter(.x, is_intended), aes(colour = subtype), size = 1.55, alpha = 0.92,
    position = position_jitter(width = 0, height = 0.13, seed = SEED)) +
  scale_colour_subtype(guide = "none", drop = TRUE) +
  scale_y_discrete(expand = expansion(add = 0.5), labels = NULL, breaks = NULL) +
  scale_x_continuous(breaks = seq(0, 12, 3), expand = expansion(mult = c(0.03, 0.03))) +
  labs(title = "Individual models (31 models)", tag = "E", x = "log2(TPM + 1)", y = NULL) +
  fig_theme() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "bottom")

# Match the heatmap and point row centres; reserve the same bottom-band height.
# The empty guide area is small and deliberate: it aligns shared marker rows.
marker_row <- pD + pE + plot_layout(widths = c(1.14, 1))
top <- pA + (pB / pC + plot_layout(heights = c(1.35, 1))) + plot_layout(widths = c(1.44, 1))
fig <- top / marker_row + plot_layout(heights = c(1, 1.26))
save_fig(fig, file.path(MSFIG, "fig3.pdf"), W2, 7.0)
save_fig(fig, file.path(MSFIG, "fig3.png"), W2, 7.0)
save_fig(pA, file.path(ASSETS, "f_rna_pca_subtype.png"), 4.4, 3.1)
save_fig(pA, file.path(ASSETS, "f_rna_pca_site.png"), 4.4, 3.1)
save_fig(pB, file.path(ASSETS, "f_pc1_commonality.png"), 3.0, 2.0)
save_fig(pC, file.path(ASSETS, "f_concordance.png"), 3.0, 2.0)
save_fig(marker_row, file.path(ASSETS, "f_rna_markers.png"), W2, 4.0)
cat(sprintf("Figure3: PCA31 models; concordance%d genes; markers%d rows; patient means%d patients.\n",
  nrow(pg), nrow(ms), nrow(rep_ann)))
message("35_fig3_biology.R complete.")
