#!/usr/bin/env Rscript
# Companion diagnostic; no numbered manuscript figure or clinical HRD classes.
ROOT <- normalizePath(Sys.getenv("OVCAN_PROJ", getwd()), mustWork = TRUE)
DEST <- file.path(ROOT, "reports/molecular_extension_2026-09-06/signatures")
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
source(file.path(ROOT, "scripts/00b_figure_theme.R"))
d <- read.csv(file.path(DEST, "sbs3_sensitivity_summary.csv"), check.names = FALSE)
burden <- read.csv(file.path(DEST, "substrate_burden.csv"), check.names = FALSE)
d <- d[d$opportunity_model == "target" & d$substrate %in% c("baseline", "rare_read_supported"), ]
stopifnot(nrow(d) == 92L, all(d$n_boots >= 100L))
ordering <- d[d$substrate == "baseline" & d$dictionary == "full", ]
ordering <- ordering[order(ordering$sbs3_fraction, ordering$cell_line), "cell_line"]
d$y <- match(d$cell_line, ordering) + ifelse(d$dictionary == "full", .13, -.13)
d$dictionary <- factor(d$dictionary, levels = c("full", "restricted"), labels = c("Full (60)", "Restricted (22)"))
palette <- c("Full (60)" = cook_rust, "Restricted (22)" = cook_slate)
base <- theme_classic(base_size = 9, base_family = "Arial") + theme(
  text = element_text(colour = cook_ink), axis.text = element_text(colour = cook_ink, size = 8),
  axis.title = element_text(size = 9), axis.line = element_line(linewidth = .3, colour = cook_ink_muted),
  axis.ticks = element_line(linewidth = .3, colour = cook_ink_muted),
  legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 9),
  plot.title = element_text(size = 10, face = "bold", hjust = 0),
  plot.margin = margin(5, 7, 4, 4), panel.grid.major.y = element_line(colour = "#F1F5F9", linewidth = .25))
one <- function(sub, title, labels = TRUE) {
  z <- d[d$substrate == sub, ]
  ggplot(z, aes(x = sbs3_fraction, y = y, colour = dictionary)) +
    geom_segment(aes(x = boot_fraction_lo95, xend = boot_fraction_hi95, yend = y), linewidth = .5, alpha = .8) +
    geom_point(size = 1.55) +
    scale_colour_manual(values = palette) +
    scale_x_continuous(limits = c(0, .8), breaks = c(0, .2, .4, .6, .8),
      labels = c("0", "20", "40", "60", "80"), expand = expansion(mult = c(.025, .025))) +
    scale_y_continuous(limits = c(.4, length(ordering) + .6), breaks = seq_along(ordering),
      labels = if (labels) ordering else NULL, expand = expansion(mult = 0)) +
    labs(x = "Fitted SBS3 contribution (%)", y = NULL, title = title) + base +
    theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
}
p1 <- one("baseline", "A  Baseline candidates")
p2 <- one("rare_read_supported", "B  Stricter candidates", labels = FALSE)
b <- burden[burden$substrate %in% c("baseline", "rare_read_supported"), ]
b$y <- match(b$cell_line, ordering)
b0 <- b[b$substrate == "baseline", ]; bs <- b[b$substrate == "rare_read_supported", ]
bs <- bs[match(b0$cell_line, bs$cell_line), ]
bc <- data.frame(y = b0$y, baseline = b0$n_snv, stricter = bs$n_snv,
                 label = paste0(b0$n_snv, " / ", bs$n_snv))
p3 <- ggplot(bc, aes(y = y)) + geom_text(aes(x = 0, label = label), hjust = 0,
  family = "Arial", size = 2.6, colour = cook_ink) +
  scale_y_continuous(limits = c(.4, length(ordering) + .6), expand = expansion(mult = 0)) +
  scale_x_continuous(limits = c(0, 1)) + labs(title = "C  SNVs", subtitle = "Baseline / stricter") +
  theme_void(base_size = 9, base_family = "Arial") + theme(
    plot.title = element_text(size = 10, face = "bold", colour = cook_ink),
    plot.subtitle = element_text(size = 7.5, colour = cook_ink_muted, margin = margin(b = 2)),
    plot.margin = margin(5, 0, 4, 4))
# Collect the two identical dictionary legends; counts are explicit beside fits.
p <- (p1 | p2 | p3) + plot_layout(widths = c(1, 1, .34), guides = "collect") &
  theme(legend.position = "bottom")
save_fig(p, file.path(DEST, "sbs3_sensitivity.pdf"), w = 7.2, h = 6.0)
save_fig(p, file.path(DEST, "sbs3_sensitivity.png"), w = 7.2, h = 6.0, dpi = 300)
message("Saved standalone SBS3 sensitivity plot in ", DEST)
