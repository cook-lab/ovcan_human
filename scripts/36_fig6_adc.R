# Figure 6: ADC-target abundance atlas and descriptive distributions.
# Scientific interpretation is in the external figure legend.
source("scripts/00_setup.R"); source("scripts/00b_figure_theme.R")
suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggnewscale)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
set.seed(SEED)
ASSETS <- file.path(PROJ, "reports", "assets"); dir.create(ASSETS, showWarnings = FALSE, recursive = TRUE)

SUB_LVL   <- c("HGS","CC","EC","MC","MMMT","SCCOHT")
ROW_ORDER <- c("FOLR1","MSLN","ERBB2","SLC34A2","TACSTD2","CDH6","CD276","VTCN1","DPEP3")
ROW_LAB <- setNames(ROW_ORDER, ROW_ORDER)

# -----------------------------------------------------------------------------
# 1. Load; establish the drawn line set and the two per-modality line sets
# -----------------------------------------------------------------------------
adc  <- read_csv(file.path(OUT, "adc_expression.csv"), show_col_types = FALSE)
summ <- read_csv(file.path(OUT, "adc_subtype_summary.csv"), show_col_types = FALSE)
mods <- read_csv(file.path(OUT, "adc_modality_line_sets.csv"), show_col_types = FALSE)
bim  <- read_csv(file.path(OUT, "adc_folr1_bimodality.csv"), show_col_types = FALSE)
dyn  <- read_csv(file.path(OUT, "prot_dynamic_range.csv"), show_col_types = FALSE)
floorchk <- read_csv(file.path(OUT, "prot_compression_floor_check.csv"), show_col_types = FALSE)

rna_lines  <- adc %>% filter(assay == "RNA")     %>% distinct(cell_line) %>% pull()
prot_lines <- adc %>% filter(assay == "protein") %>% distinct(cell_line) %>% pull()
shared     <- intersect(rna_lines, prot_lines)
N_RNA <- length(rna_lines); N_PROT <- length(prot_lines); N_SHARED <- length(shared)
stopifnot("adc_modality_line_sets.csv disagrees with adc_expression.csv" =
            setequal(mods$cell_line, union(setdiff(rna_lines, prot_lines),
                                           setdiff(prot_lines, rna_lines))))

ann <- adc %>% distinct(cell_line, subtype) %>%
  filter(cell_line %in% shared) %>%
  mutate(subtype = factor(subtype, levels = SUB_LVL)) %>%
  arrange(subtype, cell_line)
check_palette_keys(levels(droplevels(ann$subtype)), subtype_colours, "subtype")
sub_n <- ann %>% count(subtype, name = "n", .drop = TRUE)
message("Fig6 drawn line set: ", N_SHARED, " lines (both modalities) — ",
        paste(sprintf("%s %d", sub_n$subtype, sub_n$n), collapse = ", "))
message("Per-modality: RNA ", N_RNA, " lines, protein ", N_PROT, " lines; modality-specific: ",
        paste(sprintf("%s (%s, %s only)", mods$cell_line, mods$subtype, mods$only_in), collapse = "; "))

# -----------------------------------------------------------------------------
# 2. Column x-positions: subtype blocks with a gap between them
# -----------------------------------------------------------------------------
GAP <- 0.9
xpos <- ann %>% group_by(subtype) %>% mutate(within = row_number()) %>% ungroup() %>%
  mutate(blk = as.integer(subtype)) %>%
  arrange(subtype, cell_line) %>%
  mutate(x = row_number() + GAP * (match(blk, sort(unique(blk))) - 1))
blk_range <- xpos %>% group_by(subtype) %>%
  summarise(x0 = min(x), x1 = max(x), n = n(), .groups = "drop") %>%
  mutate(xm = (x0 + x1) / 2)
XMAX <- max(xpos$x)

# -----------------------------------------------------------------------------
# 3. Row y-positions: protein block below, RNA block above, subtype strip on top
# -----------------------------------------------------------------------------
NR <- length(ROW_ORDER)
BLOCK_GAP <- 1.25
y_prot <- setNames(NR - seq_along(ROW_ORDER) + 1, ROW_ORDER)              # 9..1
y_rna  <- y_prot + NR + BLOCK_GAP                                        # 11.25..19.25
Y_STRIP <- max(y_rna) + 1.5                                              # gap before the strip
Y_BLAB  <- Y_STRIP + 0.95

# -----------------------------------------------------------------------------
# 4. ROW-SCALE: z per feature across the N_SHARED drawn lines (per modality)
# -----------------------------------------------------------------------------
hm <- adc %>% filter(cell_line %in% shared, symbol %in% ROW_ORDER) %>%
  select(cell_line, subtype, symbol, assay, log2_expr) %>%
  complete(nesting(cell_line, subtype), symbol = ROW_ORDER, assay = c("RNA","protein")) %>%
  group_by(assay, symbol) %>%
  mutate(z = if (all(is.na(log2_expr))) NA_real_
             else (log2_expr - mean(log2_expr, na.rm = TRUE)) / sd(log2_expr, na.rm = TRUE),
         iqr = if (all(is.na(log2_expr))) NA_real_ else IQR(log2_expr, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(xpos %>% select(cell_line, x), by = "cell_line") %>%
  mutate(y = ifelse(assay == "RNA", y_rna[symbol], y_prot[symbol]))
# Clamp display colours, preserving the underlying z scores and distinguishing NA.
# DPEP3 RNA has IQR zero; row scaling can emphasise small absolute differences.
ZOBS <- max(abs(hm$z), na.rm = TRUE)
ZMAX <- 3
nd <- hm %>% filter(is.na(z))                      # not detected (DPEP3 protein)
ND_ROWS <- unique(nd$symbol)
message("Fig6 row scaling: z per feature across the ", N_SHARED,
        " drawn lines; |z| max ", sprintf("%.2f", max(abs(hm$z), na.rm = TRUE)),
        "; not detected: ", paste(sprintf("%s (%s)", nd$symbol[1], nd$assay[1]), collapse = ", "))

# per-row cross-line IQR in log2 units — this is what one z unit is worth
iqr_lab <- hm %>% filter(!is.na(iqr)) %>% distinct(assay, symbol, iqr) %>%
  mutate(y = ifelse(assay == "RNA", y_rna[symbol], y_prot[symbol]))
X_IQR <- XMAX + 1.0
# Identify rows with negligible central spread for the external legend.
flat_row <- iqr_lab %>% filter(iqr < 1e-8)

pA <- ggplot() +
  # ---- data cells: z per feature, centered navy<->rust ----------------------
  geom_tile(data = filter(hm, !is.na(z)), aes(x, y, fill = z),
            width = 0.94, height = 0.9) +
  scale_fill_cook_div(midpoint = 0, oob = scales::squish,
                      name = "Within-target z score",
                      limits = c(-ZMAX, ZMAX),
                      breaks = c(-ZMAX, 0, ZMAX),
                      labels = c(sprintf("≤ -%g", ZMAX), "0", sprintf("≥ %g", ZMAX)),
                      guide = guide_colourbar(order = 1, barheight = unit(4.5, "pt"),
                                              barwidth = unit(42, "pt"),
                                              direction = "horizontal",
                                              title.position = "top")) +
  new_scale_fill() +
  # ---- not-detected cells: pale tile PLUS an x glyph, with its own legend.
  #      A pale tile on its own was indistinguishable from a z near 0.
  geom_tile(data = nd, aes(x, y), fill = "#F1F5F9", width = 0.94, height = 0.9) +
  geom_point(data = nd, aes(x, y, shape = "Not quantified"), size = 1.1,
             colour = cook_ink_muted, stroke = 0.4) +
  scale_shape_manual(values = c(`Not quantified` = 4), name = NULL,
                     guide = guide_legend(order = 2, keyheight = unit(6, "pt"))) +
  # ---- subtype strip, separated by a gap and drawn with a hairline border ---
  geom_tile(data = blk_range, aes(xm, Y_STRIP, fill = subtype, width = n),
            height = 0.85, colour = cook_ink, linewidth = 0.25) +
  scale_fill_subtype(guide = "none") +
  geom_text(data = blk_range, aes(xm, Y_BLAB, label = as.character(subtype)),
            size = 2.1, colour = cook_ink, family = FIG_FONT, fontface = "bold") +
  # ---- per-row absolute spread, so row scaling does not hide the units -----
  geom_text(data = iqr_lab, aes(X_IQR, y, label = sprintf("%.2f", iqr)),
            hjust = 0, size = 2.1, colour = cook_ink_muted, family = FIG_FONT) +
  annotate("text", x = X_IQR, y = c(max(y_rna) + 0.85, max(y_prot) + 0.85), hjust = 0,
           label = "IQR", size = 2.1, colour = cook_ink, family = FIG_FONT,
           fontface = "bold") +
  # ---- block labels --------------------------------------------------------
  annotate("text", x = 0.3, y = c(max(y_rna)+0.95, max(y_prot)+0.95),
           label = c("RNA", "Protein"), hjust = 1.1, size = 2.5, colour = cook_ink,
           fontface = "bold", family = FIG_FONT) +
  scale_x_continuous(breaks = xpos$x, labels = xpos$cell_line,
                     expand = expansion(add = c(0.6, 0.6))) +
  scale_y_continuous(breaks = c(y_rna, y_prot),
                     labels = c(ROW_LAB[names(y_rna)], ROW_LAB[names(y_prot)]),
                     expand = expansion(add = c(0.6, 0.5))) +
  coord_cartesian(xlim = c(0.5, X_IQR + 2.1), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_ovcan(base_size = 7.5) +
  theme(axis.text.x = element_text(size = 6.3, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.box = "horizontal",
        legend.justification = "left", legend.margin = margin(0, 6, 0, 0),
        legend.key.size = unit(6.5, "pt"),
        plot.margin = margin(3, 3, 2, 3))

# -----------------------------------------------------------------------------
# 5. PANEL B — FOLR1 within HGSC: show the distribution, not a bimodality claim
# -----------------------------------------------------------------------------
gv <- function(s) bim$value[bim$statistic == s]
BC <- gv("bimodality_coefficient")
BC_THR <- as.numeric(sub(".*BC <= ([0-9.]+).*", "\\1",
                         bim$interpretation[bim$statistic == "bimodality_coefficient"]))
DBIC <- gv("mclust_E_dBIC_G2_minus_G1"); N_LOW <- gv("mclust_E_G2_n_low")
N_HGS_B <- gv("n_HGS_lines")
folr1 <- adc %>% filter(symbol == "FOLR1", subtype == "HGS") %>%
  mutate(assay = factor(assay, levels = c("RNA","protein")),
         facet = ifelse(assay == "RNA", "RNA, log2(TPM + 1)", "Protein, log2 normalised abundance"),
         facet = factor(facet, levels = c("RNA, log2(TPM + 1)", "Protein, log2 normalised abundance")))
message(sprintf("Panel B: FOLR1 HGSC n=%g, range %.3f-%.3f, BC %.3f (threshold %.3f), dBIC %.2f -> %s, %g lines in the low mclust component",
                N_HGS_B, gv("min"), gv("max"), BC, BC_THR, DBIC,
                ifelse(DBIC < 0, "ONE component preferred", "two components preferred"), N_LOW))

pB <- ggplot(folr1, aes(log2_expr, 1)) +
  geom_point(position = position_jitter(height = 0.16, seed = SEED), size = 1.3,
             colour = cook_rust, alpha = 0.85) +
  facet_wrap(~ facet, ncol = 1, scales = "free_x") +
  scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL, expand = c(0, 0)) +
  labs(x = "FOLR1 abundance in HGSC models", y = NULL) +
  theme_ovcan(base_size = 7.5) +
  theme(strip.text = element_text(size = 7, hjust = 0),
        axis.line.y = element_blank(),
        plot.caption = element_text(size = 5.6, hjust = 0, colour = cook_ink_muted,
                                    lineheight = 1.15, family = FIG_FONT),
        panel.spacing = unit(6, "pt"), plot.margin = margin(3, 4, 2, 3))

# -----------------------------------------------------------------------------
# 6. PANEL C — protein-vs-RNA spread using range, IQR, and SD
#    These summaries differ in sensitivity to extremes. None separates biological
#    variability from detection floors, transformation, or assay processing.
# -----------------------------------------------------------------------------
STAT_LV <- c("Range", "IQR", "SD")
cmp <- dyn %>% filter(is_adc_target) %>%
  select(gene, n_paired, Range = range_ratio, IQR = iqr_ratio, SD = sd_ratio) %>%
  pivot_longer(all_of(STAT_LV), names_to = "stat", values_to = "ratio") %>%
  mutate(stat = factor(stat, levels = STAT_LV))
N_PAIRED <- unique(dyn$n_paired[dyn$is_adc_target])
stopifnot(length(N_PAIRED) == 1)
gene_ord <- cmp %>% filter(stat == "IQR") %>% arrange(ratio) %>% pull(gene)
cmp <- cmp %>% mutate(gene = factor(gene, levels = gene_ord))
tightest <- cmp %>% filter(stat == "IQR") %>% slice_min(ratio, n = 1) %>% pull(gene)
tight_row <- dyn %>% filter(gene == as.character(tightest))
prot_med_iqr <- median(dyn$prot_iqr, na.rm = TRUE)
adc_med_iqr  <- median(dyn$prot_iqr[dyn$is_adc_target], na.rm = TRUE)
stat_pal <- setNames(c(cook_ink_muted, cook_rust, cook_navy), STAT_LV)
message(sprintf("Panel C: %s range %.3f / IQR %.3f / SD %.3f (n=%d paired lines) | median protein IQR: whole panel %.3f, ADC targets %.3f",
                tightest, tight_row$range_ratio, tight_row$iqr_ratio, tight_row$sd_ratio,
                N_PAIRED, prot_med_iqr, adc_med_iqr))

pC <- ggplot(cmp, aes(ratio, gene, colour = stat, shape = stat)) +
  geom_line(aes(group = gene), colour = cook_hair, linewidth = 0.35, orientation = "y") +
  geom_point(size = 1.5) +
  scale_colour_manual(values = stat_pal, name = NULL,
                      guide = guide_legend(nrow = 1, override.aes = list(size = 2.1))) +
  scale_shape_manual(values = c(Range = 1, IQR = 16, SD = 17),
                     name = NULL, guide = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(0, 0.52), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.6, 8.4), clip = "off") +
  labs(x = "Protein / RNA spread", y = NULL) +
  theme_ovcan(base_size = 7.5) +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, -3, 0), legend.key.size = unit(6.5, "pt"),
        plot.margin = margin(3, 4, 2, 3))

# -----------------------------------------------------------------------------
# 7. Assemble + export
# -----------------------------------------------------------------------------
fig6 <- (pA / (pB | pC)) +
  plot_layout(heights = c(3.3, 1.8)) + plot_annotation(tag_levels = "A") &
  theme(text = element_text(family = FIG_FONT),
        plot.tag = element_text(size = 11, face = "bold", family = FIG_FONT))
save_fig(fig6, file.path(MSFIG, "fig6.pdf"), w = W2, h = 5.3)
save_fig(fig6, file.path(MSFIG, "fig6.png"), w = W2, h = 5.3)
save_fig(pA, file.path(ASSETS, "f_adc.png"), w = W2, h = 3.5)

# -----------------------------------------------------------------------------
# 8. Console check vs the recovery summary (per-modality n now carried)
# -----------------------------------------------------------------------------
cat(sprintf("\nFig6 written: %s  (%d targets x %d lines; RNA + protein, row-scaled)\n",
            file.path(MSFIG, "fig6.pdf"), length(ROW_ORDER), N_SHARED))
cat("Known-association recovery (from adc_subtype_summary.csv; note per-modality n):\n")
print(as.data.frame(summ %>% filter(!is.na(expected)) %>%
        select(symbol, common, assay, expected, top_mean_subtype, top_n_measured,
               top_n_lines_in_modality, top_mean_matches, top_line_in_expected)),
      row.names = FALSE)
message("36_fig6_adc.R complete.")
