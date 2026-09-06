# Figure 2: RNA/protein quality and recovered WES technical validation.
# WES values are independently parsed by scripts/24_wes_recovered_qc.py.
# External legend: reports/figure_redesign_2026-09-05/legends_expression.md
source("scripts/00_setup.R"); source("scripts/00b_figure_theme.R")
suppressPackageStartupMessages({library(tidyverse); library(patchwork); library(readxl); library(ggrastr)})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; count <- dplyr::count
set.seed(SEED)
ASSETS <- file.path(PROJ, "reports", "assets")
fig_theme <- function() theme_ovcan(base_size = 8.3) + theme(
  plot.title = element_text(size = 8.5, face = "bold", margin = margin(b = 5)),
  plot.tag = element_text(size = 11, face = "bold", family = FIG_FONT),
  plot.tag.position = "topleft", axis.text = element_text(size = 7.2),
  axis.title = element_text(size = 8), legend.text = element_text(size = 7),
  legend.title = element_text(size = 7), legend.key.size = unit(8, "pt"),
  legend.margin = margin(0, 0, 0, 0), plot.margin = margin(5, 5, 5, 5))
SITE_LVL <- c("Mes-Masson", "Huntsman", "Huntsman/Vanderhyden")
SITE_DISP <- c("Mes-Masson" = "CHUM", "Huntsman" = "BC Cancer", "Huntsman/Vanderhyden" = "OHRI")
site_col <- setNames(site_colours[SITE_LVL], SITE_DISP[SITE_LVL])
site_shp <- setNames(site_shapes[SITE_LVL], SITE_DISP[SITE_LVL])
qc <- read_csv(file.path(OUT, "rna_qc_metrics.csv"), show_col_types = FALSE) %>%
  mutate(site_f = factor(SITE_DISP[site], levels = SITE_DISP[SITE_LVL]),
         detected_k = detected / 1000, frag_m = n_processed_fragments / 1e6)
sc <- read_csv(file.path(OUT, "rna_qc_site_comparison.csv"), show_col_types = FALSE)
site_m <- sc %>% filter(block == "per-site medians") %>%
  transmute(site_f = factor(SITE_DISP[site], levels = SITE_DISP[SITE_LVL]),
            pseudoalign = pseudoalign_median, detected_k = detected_median / 1000,
            frag_m = n_processed_fragments_median / 1e6)
stopifnot(nrow(qc) == 31L, !anyNA(qc$site_f))
y_rng <- range(qc$detected_k) + c(-0.13, 0.15)
rna_panel <- function(xvar, xlab) {
  ggplot(qc, aes(.data[[xvar]], detected_k)) +
    geom_point(aes(fill = site_f, shape = site_f), colour = cook_ink, size = 2.1, stroke = 0.25) +
    geom_point(data = site_m, aes(shape = site_f), fill = NA, colour = cook_ink,
               size = 3.7, stroke = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = site_col, name = NULL) +
    scale_shape_manual(values = site_shp, name = NULL) +
    scale_y_continuous(limits = y_rng, breaks = 19:22, expand = c(0, 0)) +
    scale_x_continuous(n.breaks = 4, expand = expansion(mult = 0.05)) +
    labs(x = xlab, y = "Detected genes (thousands)") + fig_theme()
}
pA <- rna_panel("pseudoalign", "Pseudoalignment (%)") +
  labs(title = "RNA alignment (31 models)", tag = "A") + theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
pB <- rna_panel("frag_m", "Fragments (millions)") +
  labs(title = "Library depth", tag = "B", y = NULL) +
  guides(fill = "none", shape = "none") + theme(legend.position = "none")
# The shared RNA legend belongs to this row, outside the data region.
rna_row <- pA + pB + plot_layout(widths = c(1.06, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "left")

# C: paired spread. Retain all genes in summaries; visual range is the pooled 99th percentile.
dr <- read_csv(file.path(OUT, "prot_dynamic_range.csv"), show_col_types = FALSE)
long_iqr <- dr %>% select(gene, RNA = rna_iqr, Protein = prot_iqr) %>%
  pivot_longer(-gene, names_to = "assay", values_to = "iqr") %>%
  mutate(assay = factor(assay, levels = c("RNA", "Protein")))
med_iqr <- long_iqr %>% group_by(assay) %>% summarise(median = median(iqr), .groups = "drop")
IQR_CAP <- as.numeric(quantile(long_iqr$iqr, 0.99))
pC <- ggplot(long_iqr, aes(assay, iqr, fill = assay)) +
  geom_violin(colour = NA, alpha = 0.72, width = 0.8, scale = "width") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.3, colour = cook_ink) +
  geom_text(data = med_iqr, aes(y = IQR_CAP * 0.96, label = sprintf("%.2f", median)),
            size = 2.6, colour = cook_ink) +
  scale_fill_manual(values = c(RNA = cook_slate, Protein = cook_rust), guide = "none") +
  scale_y_continuous(breaks = 0:4, expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, IQR_CAP * 1.05)) +
  labs(title = "Cross-model spread", tag = "C", x = NULL, y = "IQR (log2 units)") + fig_theme()

# D: log-position encoding makes the small coverage classes legible without a broken axis.
bm <- read_csv(file.path(OUT, "prot_block_missingness.csv"), show_col_types = FALSE)
tier <- tibble(present_n_plex = 0:5) %>% left_join(bm %>% count(present_n_plex), by = "present_n_plex") %>%
  mutate(n = replace_na(n, 0L), key = ifelse(present_n_plex == 0, "No model channels", "Quantified"))
stopifnot(sum(tier$n) == 8427L, tier$n[tier$present_n_plex == 0] == 70L)
pD <- ggplot(tier, aes(n, present_n_plex)) +
  geom_segment(aes(x = 30, xend = n, yend = present_n_plex), colour = cook_grey, linewidth = 0.5) +
  geom_point(aes(colour = key), size = 2.4) +
  geom_text(aes(label = format(n, big.mark = ",", trim = TRUE)), hjust = -0.32,
            size = 2.5, colour = cook_ink) +
  scale_colour_manual(values = c(`No model channels` = cook_rust, Quantified = cook_slate), guide = "none") +
  scale_x_log10(limits = c(30, 35000), breaks = c(100, 1000, 10000), labels = c("100", "1,000", "10,000"),
                expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:5, expand = expansion(add = 0.5)) +
  labs(title = "Protein coverage", tag = "D", x = "Features (log scale)", y = "Quantified plexes") + fig_theme()

# E: reported and bridge-derived CV are distinct quantities; only the reported series has IQR.
cvd <- read_csv(file.path(OUT, "prot_cv_by_abundance.csv"), show_col_types = FALSE)
stopifnot(nrow(cvd) == 10L)
pE <- ggplot(cvd, aes(abundance_decile)) +
  geom_ribbon(aes(ymin = vendor_cv_q25, ymax = vendor_cv_q75), fill = cook_slate, alpha = 0.14) +
  geom_line(aes(y = vendor_cv_median, colour = "Reported"), linewidth = 0.6) +
  geom_line(aes(y = bridge_cv_pct_per_measurement, colour = "Bridge-derived"), linewidth = 0.6) +
  geom_point(aes(y = vendor_cv_median, colour = "Reported"), size = 1.7) +
  geom_point(aes(y = bridge_cv_pct_per_measurement, colour = "Bridge-derived"), size = 1.7) +
  scale_colour_manual(values = c(Reported = cook_slate, `Bridge-derived` = cook_rust), name = NULL) +
  scale_x_continuous(breaks = 1:10, expand = expansion(mult = 0.03)) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(0, 20, 5), expand = c(0, 0)) +
  labs(title = "Precision by abundance", tag = "E", x = "Protein abundance decile", y = "CV (%)") +
  fig_theme() + theme(legend.position = "top", legend.justification = "left")

bridge_ref <- read_csv(file.path(OUT, "prot_bridge_cor.csv"), show_col_types = FALSE)
bridge_agr <- read_csv(file.path(OUT, "prot_bridge_agreement.csv"), show_col_types = FALSE)
PROT_DIR <- file.path(DATA, "proteomics")
abund <- read_excel(file.path(PROT_DIR, "protein_relative_abundance.xlsx"))
lay   <- read_excel(file.path(PROT_DIR, "tmt.layout.xlsx")) %>%
  transmute(id, plex, tmt_label = `TMT.label`, name)
data_cols <- names(abund)[12:ncol(abund)]
cmap <- tibble(data_col = data_cols) %>% left_join(lay, by = c("data_col" = "id")) %>%
  mutate(role = case_when(grepl("^SM\\.iD", data_col) ~ "SM.iD",
                          tmt_label == 10 & plex >= 2 ~ "bridge",
                          TRUE ~ "sample"))
brs <- cmap %>% filter(role == "bridge") %>% transmute(name, bridge_col = data_col, bridge_plex = plex)
prm <- cmap %>% filter(role == "sample") %>% transmute(name, prim_col = data_col, prim_plex = plex)
bp  <- brs %>% left_join(prm, by = "name") %>% filter(!is.na(prim_col)) %>% arrange(prim_plex)

bridge_long <- purrr::map_dfr(seq_len(nrow(bp)), function(i) {
  x <- as.numeric(abund[[bp$prim_col[i]]]); y <- as.numeric(abund[[bp$bridge_col[i]]])
  ok <- is.finite(x) & is.finite(y)
  tibble(link = i, prim_plex = bp$prim_plex[i], bridge_plex = bp$bridge_plex[i],
         name = bp$name[i], prim = x[ok], bridge = y[ok])
})
# per-link stats + attach canonical labels from the reference table
bstat <- bridge_long %>% group_by(prim_plex, bridge_plex) %>%
  summarise(n = n(), pearson = cor(prim, bridge), .groups = "drop") %>%
  left_join(bridge_ref %>% select(cell_line, prim_plex, bridge_plex, pearson_ref = pearson, external),
            by = c("prim_plex","bridge_plex"))
stopifnot("bridge Pearson must match prot_bridge_cor.csv (re-render, not recompute)" =
            all(abs(bstat$pearson - bstat$pearson_ref) < 1e-3))
message("Bridge correlations validated against output/prot_bridge_cor.csv.")

# agreement statistics come from output/, and the n / bias / SD are re-derived from
# the same values plotted here so the panel and the table cannot disagree
agr <- bstat %>%
  left_join(bridge_agr %>% select(prim_plex, bridge_plex, n_proteins, bias, sd_diff,
                                  loa_lower, loa_upper, sd_single,
                                  cv_pct_per_measurement, median_abs_diff),
            by = c("prim_plex", "bridge_plex")) %>%
  mutate(loa_span = loa_upper - loa_lower)
stopifnot("agreement table must cover all four links" = !any(is.na(agr$sd_diff)))


# F: four wide facets form a single right-hand column with consistent axes.
bridge_long <- bridge_long %>% mutate(mean_ab = (prim + bridge) / 2, diff = prim - bridge)
chk <- bridge_long %>% group_by(prim_plex, bridge_plex) %>%
  summarise(n_check = n(), bias_check = mean(diff), sd_check = sd(diff), .groups = "drop")
agr <- agr %>% left_join(chk, by = c("prim_plex", "bridge_plex"))
stopifnot(all(agr$n_proteins == agr$n_check), max(abs(agr$bias - agr$bias_check)) < 1e-10,
          max(abs(agr$sd_diff - agr$sd_check)) < 1e-10)
agr <- agr %>% mutate(facet = sprintf("Plex %d-%d | %s%s\nn = %s",
  prim_plex, bridge_plex, cell_line, ifelse(external, "*", ""), format(n_proteins, big.mark = ",")))
bridge_long <- bridge_long %>% left_join(agr %>% select(prim_plex, bridge_plex, facet),
  by = c("prim_plex", "bridge_plex")) %>% mutate(facet = factor(facet, levels = agr$facet))
agr$facet <- factor(agr$facet, levels = agr$facet)
pF <- ggplot(bridge_long, aes(mean_ab, diff)) +
  geom_hline(yintercept = 0, colour = cook_hair, linewidth = 0.3) +
  rasterise(geom_point(size = 0.25, colour = cook_slate, alpha = 0.20), dpi = 500) +
  geom_hline(data = agr, aes(yintercept = bias), colour = cook_rust, linewidth = 0.45) +
  geom_hline(data = agr, aes(yintercept = loa_lower), colour = cook_rust, linewidth = 0.4, linetype = "dashed") +
  geom_hline(data = agr, aes(yintercept = loa_upper), colour = cook_rust, linewidth = 0.4, linetype = "dashed") +
  facet_wrap(~ facet, ncol = 1) +
  scale_x_continuous(breaks = c(6, 10, 14, 18), expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = c(-1, 0, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-1.6, 1.6)) +
  labs(title = "Bridge agreement", tag = "F", x = "Mean log2 abundance", y = "Primary - bridge (log2)") +
  fig_theme() + theme(strip.text = element_text(size = 7.2, hjust = 0, lineheight = 1.05),
                      panel.spacing.y = unit(8, "pt"))

left <- wrap_plots(rna_row, pC + pD + plot_layout(widths = c(1.02, 1.04)), pE,
  ncol = 1, heights = c(1.1, 1, 0.87))
expression_qc <- wrap_plots(left, pF, ncol = 2, widths = c(1.76, 1))

# G/H: the original mark-duplicate alignment stage, not interval-restricted
# recalibrated CRAMs. Coverage fractions are rounded to 0.01 in mosdepth reports.
wes_qc <- read_csv(file.path(OUT, "wes_qc_model_summary.csv"), show_col_types = FALSE)
stopifnot(nrow(wes_qc) == 23L, n_distinct(wes_qc$cell_line) == 23L)
pG <- ggplot(wes_qc, aes(mosdepth_md_mean_target_depth_x,
                        100 * mosdepth_md_fraction_target_ge_30x)) +
  geom_point(colour = cook_slate, size = 1.9, alpha = 0.85) +
  geom_text(data = wes_qc %>% filter(cell_line == "TOV2929D"),
            aes(label = cell_line), nudge_x = 0.6, nudge_y = 1.6,
            hjust = 0, size = 2.4, colour = cook_ink) +
  scale_x_continuous(limits = c(68, 91), breaks = seq(70, 90, 5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(63, 88), breaks = seq(65, 85, 5), expand = c(0, 0)) +
  labs(title = "Exome depth (23 models)", tag = "G", x = "Mean target depth (×)",
       y = "Target bases at ≥30× (%)") + fig_theme()

# Each line is a model; the median is calculated separately at every depth.
wes_profile <- read_csv(file.path(OUT, "wes_qc_coverage_profile.csv"), show_col_types = FALSE) %>%
  filter(stage == "md", depth_x <= 150)
stopifnot(n_distinct(wes_profile$cell_line) == 23L)
wes_median <- wes_profile %>% group_by(depth_x) %>%
  summarise(fraction_ge_depth = median(fraction_ge_depth), .groups = "drop")
pH <- ggplot(wes_profile, aes(depth_x, 100 * fraction_ge_depth)) +
  geom_line(aes(group = cell_line), colour = cook_slate, alpha = 0.25, linewidth = 0.35) +
  geom_line(data = wes_median, colour = cook_rust, linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 30, 60, 90, 120, 150), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(title = "Exome coverage", tag = "H", x = "Minimum depth (×)",
       y = "Target bases covered (%)") + fig_theme()

fig <- wrap_plots(expression_qc, pG + pH + plot_layout(widths = c(1, 1)),
                  ncol = 1, heights = c(6.15, 2.05))
save_fig(fig, file.path(MSFIG, "fig2.pdf"), W2, 8.2)
save_fig(fig, file.path(MSFIG, "fig2.png"), W2, 8.2)
save_fig(rna_row, file.path(ASSETS, "f_rna_qc.png"), 4.6, 2.5)
save_fig(pF, file.path(ASSETS, "f_prot_bridge.png"), 3.3, 6.0)
save_fig(pC + pD, file.path(ASSETS, "f_prot_compression.png"), 4.6, 2.5)
cat(sprintf("Figure2: %d RNA models; %d paired-spread genes; %d protein features; %d bridge points.\n",
  nrow(qc), nrow(dr), nrow(bm), nrow(bridge_long)))
cat(sprintf("Bridge display clips %.3f%% of points beyond +/-1.6 log2.\n", 100 * mean(abs(bridge_long$diff) > 1.6)))
message("34_fig2_qc.R complete.")
