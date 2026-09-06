# Report diagnostic, not a numbered publication figure.
# Reproduces FGA(center) for the archived and corrected OV1369-R2 segments.
# Inputs and exact-center checks: script29c; no CNV fitting or cohort changes.
PROJ <- normalizePath(Sys.getenv("OVCAN_PROJ", unset = getwd()), mustWork = TRUE)
REPORT <- file.path(PROJ, "reports", "wes_cnv_coverage_2026-09-06")
source(file.path(PROJ, "scripts", "00b_figure_theme.R"))
suppressPackageStartupMessages(library(ggplot2))

read_csv_base <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
source_path <- function(path) {
  data_root <- Sys.getenv("OVCAN_DATA", unset = "")
  if (nzchar(data_root) && startsWith(path, "judy_archive/data/"))
    return(file.path(data_root, substring(path, nchar("judy_archive/data/")+1L)))
  file.path(PROJ, path)
}
check_sha <- function(path, expected) {
  stopifnot(file.exists(path), digest::digest(file = path, algo = "sha256") == expected)
}
map <- read_csv_base(file.path(PROJ, "output", "wes_recovered_provenance_cnv_support.csv"))
map <- map[map$cell_line == "OV1369-R2", , drop = FALSE]
new_map <- read_csv_base(file.path(PROJ, "output", "wes_cnv_target_only", "manifest.csv"))
new_map <- new_map[new_map$cell_line == "OV1369-R2", , drop = FALSE]
stopifnot(nrow(map) == 1L, nrow(new_map) == 1L)
old_path <- source_path(map$cns_source)
new_path <- file.path(PROJ, new_map$cns_path)
check_sha(old_path, map$cns_sha256); check_sha(new_path, new_map$cns_sha256)
auto <- paste0("chr", 1:22)
read_segments <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[x$chromosome %in% auto, , drop = FALSE]
  x$span <- x$end - x$start
  stopifnot(all(x$span > 0), all(is.finite(x$log2)))
  x
}
segments <- list(Archived = read_segments(old_path), `Target-only` = read_segments(new_path))
two_by_two <- read_csv_base(file.path(REPORT, "ov1369_centering_2x2.csv"))
alternatives <- read_csv_base(file.path(REPORT, "ov1369_alternative_centers.csv"))
old_center <- unique(two_by_two$center[two_by_two$center_source == "archived"])
new_center <- unique(two_by_two$center[two_by_two$center_source == "target_only"])
stopifnot(length(old_center) == 1L, length(new_center) == 1L)
bin_centers <- alternatives$center[grepl("^positive_target_bin_", alternatives$center_definition)]
stopifnot(length(bin_centers) == 2L)

fga <- function(x, center) sum(x$span[abs(x$log2-center) > .20]) / sum(x$span)
for (i in seq_len(nrow(two_by_two))) {
  row <- two_by_two[i, ]
  profile <- if (row$segment_source == "archived") "Archived" else "Target-only"
  stopifnot(abs(fga(segments[[profile]], row$center) - row$fga) < 1e-12)
}
for (i in seq_len(nrow(alternatives)))
  stopifnot(abs(fga(segments[["Target-only"]], alternatives$center[i]) -
                  alternatives$fga_on_target_only_segments[i]) < 1e-12)

# Dense evaluation grid plus exact declared/alternative centers. Curves are
# threshold sensitivity functions; no smooth model is fitted to these points.
centers <- sort(unique(c(seq(-.30, .10, by = .0001), old_center, new_center, bin_centers)))
curve <- do.call(rbind, lapply(names(segments), function(profile) {
  x <- segments[[profile]]
  data.frame(model = "OV1369-R2", profile = profile, center = centers,
    autosomal_span_bp = sum(x$span),
    fga = vapply(centers, function(center) fga(x, center), numeric(1)))
}))
write.csv(curve, file.path(REPORT, "ov1369_centering_curve.csv"), row.names = FALSE)
declared <- data.frame(profile = c("Archived", "Target-only"), center = c(old_center, new_center))
declared$fga <- mapply(function(profile, center) fga(segments[[profile]], center), declared$profile, declared$center)
declared$label <- c(sprintf("Archived median\nFGA %.3f", declared$fga[1]),
                    sprintf("Target-only median\nFGA %.3f", declared$fga[2]))
declared$label_x <- c(-.222, -.151)
declared$label_y <- c(.565, .93)
bin_points <- data.frame(center = bin_centers,
                        fga = vapply(bin_centers, function(x) fga(segments[["Target-only"]], x), numeric(1)))
colours <- c(Archived = cook_slate, `Target-only` = cook_rust)
p <- ggplot(curve, aes(center, fga, colour = profile)) +
  annotate("rect", xmin = min(bin_centers), xmax = max(bin_centers), ymin = .45, ymax = 1,
           fill = cook_grey, alpha = .65) +
  geom_hline(yintercept = seq(.5, 1, .1), colour = "#E8EDF2", linewidth = .3) +
  geom_line(linewidth = .8) +
  geom_segment(data = declared, aes(xend = center, y = .45, yend = fga),
               linewidth = .4, linetype = "dotted", show.legend = FALSE) +
  geom_point(data = declared, aes(center, fga), size = 2.7, show.legend = FALSE) +
  geom_point(data = bin_points, aes(center, fga), inherit.aes = FALSE,
             shape = 23, fill = "white", colour = cook_ink, size = 2.5, stroke = .65) +
  geom_segment(data = declared,
               aes(x = label_x, y = label_y, xend = center, yend = fga),
               linewidth = .35, show.legend = FALSE) +
  geom_label(data = declared, aes(x = label_x, y = label_y, label = label),
             fill = "white", linewidth = 0, label.padding = unit(2, "pt"),
             size = 3, lineheight = 1.05, show.legend = FALSE) +
  annotate("segment", x = min(bin_centers), xend = max(bin_centers),
           y = .507, yend = .507, colour = cook_ink_muted, linewidth = .4) +
  annotate("text", x = mean(bin_centers), y = .484, label = "Target-bin medians",
           size = 2.8, colour = cook_ink, vjust = .5) +
  scale_colour_manual(values = colours, name = NULL) +
  scale_x_continuous(breaks = seq(-.3, .1, .1), labels = function(x) sprintf("%.1f", x),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(.5, 1, .1), labels = function(x) sprintf("%.1f", x),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(-.3, .1), ylim = c(.45, 1.01), clip = "off") +
  labs(title = "OV1369-R2", x = expression("Center subtracted from segment log"[2]*" ratio"),
       y = "Fraction of autosomal segment spans altered") +
  theme_ovcan(10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        legend.position = "top", legend.justification = "left", legend.key.width = unit(22, "pt"),
        legend.margin = margin(0, 0, 2, 0),
        axis.title.x = element_text(margin = margin(t = 7)),
        axis.title.y = element_text(margin = margin(r = 7)),
        plot.margin = margin(7, 9, 7, 7))
save_fig(p, file.path(REPORT, "ov1369_centering_sensitivity.pdf"), w = 7.2, h = 3.8)
save_fig(p, file.path(REPORT, "ov1369_centering_sensitivity.png"), w = 7.2, h = 3.8, dpi = 300)
message("Report plot complete; all four declared-center values and alternative-center values match script29c.")
