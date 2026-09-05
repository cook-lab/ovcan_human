# Shared scientific-figure theme. Source after 00_setup.R.
# Brand reference: ~/Lab/Branding/docs/design-system.md and tokens/colors_and_type.css.
# Journal adaptation: Arial/Helvetica, rust with brand slate (#64748B), white pages.
# Signed quantities use a zero-centred slate-white-rust ramp; unsigned quantities
# use a sequential ramp. Semantic keys remain consistent across figures.
# Assay presence and confidence tiers use neutral fills. Variant classes use
# shapes in the redesigned oncoprint, and patient brackets supplement grouping.
# Explanatory narrative belongs in docs/manuscript/figures/figure_legends.md;
# preserve axes, units, panel tags and essential graphical keys in the exports.

suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(scales)
})

# ---- 1. Raw branding palettes (rustNavy) ------------------------------------
cook_diverging  <- c("#64748B","#94A3B8","#CBD5E1","#FFFFFF","#FABEA0","#EB6235","#C2410C") # slate <-> rust
cook_sequential <- c("#FEF3EE","#FDE2D1","#FABEA0","#F59065","#EB6235","#C2410C","#9C3409","#7A2A09") # rust ramp
cook_categorical<- c("#C2410C","#64748B","#0D9488","#A16207","#7C3AED","#94A3B8","#DC2626","#0369A1")

cook_rust <- "#C2410C"; cook_slate <- "#64748B"; cook_teal <- "#0D9488"
# cook_navy remains a compatibility name for older plotting code; marks now use slate.
cook_navy <- cook_slate
cook_ink <- "#334155"; cook_ink_muted <- "#64748B"; cook_grey <- "#E2E8F0"; cook_hair <- "#CBD5E1"

# ---- 2. Locked semantic mappings (SAME colour in every figure) --------------
# Subtypes (7; the panel's full set). HGS = rust (dominant/primary series).
subtype_colours <- c(
  HGS    = "#C2410C",  # rust
  LGS    = "#F59065",  # light rust (serous family, lighter)
  CC     = "#0D9488",  # teal
  EC     = "#7C3AED",  # violet
  MC     = "#A16207",  # ochre
  MMMT   = "#0369A1",  # blue
  SCCOHT = "#64748B"   # brand slate
)
# Contributing centre (source site).
# LOW-CHROMA WARM NEUTRAL, ordered by lightness — deliberately NOT a saturated
# palette. Two reasons: (i) site is a nuisance/design variable whose role in these
# figures is to show whether it structures the data, so it should not compete with
# the semantic subtype colours; (ii) the previous mapping WAS the subtype palette
# (Mes-Masson=HGS rust, Huntsman=MMMT blue, OHRI=CC teal), so a subtype-coloured
# panel beside a site-coloured panel invited direct misreading. Lightness-ordered
# neutrals are colour-vision-safe by construction and cannot be confused with the
# cool slate used for TMT plex.
#
# Keys cover BOTH the institutional names (preferred in text and captions, per
# figure review) and the exact strings in metadata/samples.csv. The metadata form
# "Huntsman/Vanderhyden" previously had NO entry here and would have silently
# rendered NA — the OHRI-derived line is that row.
site_colours <- c(
  # institutional (preferred for display)
  "CHUM"        = "#3F3A34", "BC Cancer" = "#8A8074", "OVCARE" = "#8A8074",
  "OHRI"        = "#C9C0B2",
  # metadata strings (must match samples.csv exactly)
  "Mes-Masson"  = "#3F3A34", "Huntsman"  = "#8A8074",
  "Huntsman/Vanderhyden" = "#C9C0B2", "Vanderhyden" = "#C9C0B2"
)
# NOTE: "Carey" is deliberately ABSENT from site_colours. samples.csv contains 13
# Carey/OVCARE LGS rows with provenance = external and analysis_include = N; they
# are not part of the 42-model resource (29 Mes-Masson + 12 Huntsman + 1
# Huntsman/Vanderhyden = 42). Leaving the key out means check_palette_keys() fails
# loudly if a figure ever plots unfiltered metadata, which is the behaviour we want.
# Double-encode site by shape wherever it is plotted as points, so the panel does
# not depend on low-chroma fills alone. Filled shapes (21/22/24) take a border.
site_shapes <- c(
  "CHUM" = 21, "BC Cancer" = 24, "OVCARE" = 24, "OHRI" = 22,
  "Mes-Masson" = 21, "Huntsman" = 24,
  "Huntsman/Vanderhyden" = 22, "Vanderhyden" = 22
)
# Somatic-confidence tiers (ordered: dark = highest confidence).
# NEUTRAL ordinal ramp by design. The previous rust ramp put Tier 2 (#EB6235)
# between the HGS (#C2410C) and LGS (#F59065) subtype swatches, which sit in the
# same oncoprint as an annotation strip — three variables, one hue family.
tier_colours <- c("Tier 1" = "#64748B", "Tier 2" = "#94A3B8", "Tier 3" = "#E2E8F0")

# Legacy class palette retained for reports; the manuscript oncoprint now uses
# class-specific glyphs, avoiding a second competing categorical fill scale.
variant_class_colours <- c(
  truncating   = "#64748B",  # slate — nonsense / frameshift / splice
  missense     = "#C2410C",  # rust
  inframe      = "#A16207",  # ochre — in-frame indel
  multi_hit    = "#7C3AED"   # violet
)

# Assay availability uses slate fill and white or pale unprofiled cells.
present_colours <- c("Yes" = "#64748B", "No" = "#F1F5F9", "TRUE" = "#64748B", "FALSE" = "#F1F5F9")

# TMT plex (batch; ordinal 1-5). Slate ramp — a nuisance variable must not borrow
# the subtype palette, which is what made the two protein-PCA panels misreadable.
plex_colours <- c("1" = "#C3D0DE", "2" = "#93A7BC", "3" = "#64798F",
                  "4" = "#41576E", "5" = "#22384F")

# Patient families (sublines of one patient). LOCKED, and deliberately a single-hue
# PLUM ramp rather than five distinct hues. Two problems are solved at once:
#   (i) the same five families were drawn in two different palettes across figures;
#  (ii) the old mapping reused four subtype hues (rust/teal/violet/blue), and the
#       family strip sits directly beside the subtype strip in the overview and
#       oncoprint figures — so a reader could read a family colour as a subtype.
# Plum is used by no other locked palette, so the family strip is unambiguous, and
# five ordered lightness steps still let same-family rows be matched at a glance.
# Where layout permits, ALSO bracket or label same-family rows: family is a grouping
# variable and a second channel reads better than colour alone.
family_colours <- c(
  "1369" = "#F0ABCB",  # light pink
  "2295" = "#DB2777",  # pink
  "3121" = "#BE185D",  # deep pink
  "3133" = "#831843",  # wine
  "3291" = "#9D174D"   # wine; avoid near-black block
)

# Within-HGSC descriptive pathway strata. LOCKED — these were defined per script
# and the same stratum name carried different colours in adjacent supplementary
# figures (and one colour mapped to two different strata). Semantics: warm =
# inflammatory, grey = low-signalling, teal = hypoxic.
# Keys cover BOTH the short display form and the exact string the data carries:
# hgs_heterogeneity.csv$cluster_label is "Inflammatory/NF-kB-EMT", so a bare
# "Inflammatory" key made check_palette_keys() fail on the raw label and forced two
# figure scripts to define the colour locally — the exact drift this palette exists to
# prevent. Same dual-spelling approach as site_colours.
stratum_colours <- c(
  "Inflammatory"           = "#C2410C",
  "Inflammatory/NF-kB-EMT" = "#C2410C",   # raw label in hgs_heterogeneity.csv
  "Low-signaling"          = "#94A3B8",
  "Low-signalling"         = "#94A3B8",   # tolerate both spellings
  "Hypoxic-glycolytic"     = "#0D9488"
)

# ConsensusOV HGSC class calls. IMR and MES are the two microenvironment-driven
# classes and were previously the closest pair in the palette; separated here.
consensusov_colours <- c(
  DIF = "#0369A1",  # differentiated  — blue
  PRO = "#7C3AED",  # proliferative   — violet
  IMR = "#A16207",  # immunoreactive  — ochre
  MES = "#64748B"   # mesenchymal     — slate
)

# ---- 3. Fonts — MANUSCRIPT figures use a standard sans (Arial/Helvetica) -----
# Deliberate deviation from the Inter brand font: journals expect Arial/Helvetica
# and standard fonts embed/render portably. (The HTML report keeps Inter via CSS.)
.have_font <- function(f) {
  if (!requireNamespace("systemfonts", quietly = TRUE)) return(FALSE)
  any(tolower(systemfonts::system_fonts()$family) == tolower(f))
}
.pick_font <- function(cands) { for (f in cands) if (.have_font(f)) return(f); "" }
FIG_FONT      <- .pick_font(c("Arial", "Helvetica", "Helvetica Neue", "Liberation Sans", "DejaVu Sans"))
FIG_FONT_DISP <- FIG_FONT   # no separate display face in manuscript figures

# ---- 4. Dense manuscript ggplot theme ---------------------------------------
# base_size 8 tuned for MULTI-PANEL composites; bump to 9-10 for standalone.
theme_ovcan <- function(base_size = 8, base_family = FIG_FONT) {
  h <- base_size / 2
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # family MUST be set on the root text element: %+replace% drops the
      # base_family that theme_classic() set, and cairo_pdf then falls back to
      # Helvetica. Setting it here propagates Arial to all inheriting text (PDF + PNG).
      text            = element_text(colour = cook_ink, family = base_family),
      # titles kept minimal by default — prefer labs(title=NULL) and let the caption narrate
      plot.title      = element_text(size = base_size, hjust = 0, colour = cook_ink,
                                     margin = margin(b = h), family = base_family),
      plot.subtitle   = element_text(size = base_size - 0.5, hjust = 0, colour = cook_ink_muted,
                                     margin = margin(b = h)),
      plot.tag        = element_text(size = base_size + 3, face = "bold", colour = cook_ink,
                                     family = FIG_FONT_DISP),
      axis.title      = element_text(size = base_size, colour = cook_ink),
      axis.text       = element_text(size = base_size - 1, colour = cook_ink),
      axis.line       = element_line(linewidth = 0.3, colour = cook_ink),
      axis.ticks      = element_line(linewidth = 0.3, colour = cook_ink),
      axis.ticks.length = unit(2, "pt"),
      legend.title    = element_text(size = base_size - 1, colour = cook_ink),
      legend.text     = element_text(size = base_size - 1, colour = cook_ink),
      legend.key.size = unit(9, "pt"),
      legend.margin   = margin(1, 1, 1, 1),
      legend.box.spacing = unit(3, "pt"),
      legend.background  = element_blank(), legend.key = element_blank(),
      strip.text      = element_text(size = base_size, colour = cook_ink, margin = margin(2,2,2,2)),
      strip.background= element_blank(),
      panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# ---- 4b. In-plot text geoms use the figure font too -------------------------
# Geom defaults are NOT controlled by the theme; without this, geom_text/label/
# annotate() and ggrepel fall back to the cairo_pdf default (Helvetica). Set once.
if (FIG_FONT != "") {
  try(update_geom_defaults("text",  list(family = FIG_FONT)), silent = TRUE)
  try(update_geom_defaults("label", list(family = FIG_FONT)), silent = TRUE)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    try(update_geom_defaults("text_repel",  list(family = FIG_FONT)), silent = TRUE)
    try(update_geom_defaults("label_repel", list(family = FIG_FONT)), silent = TRUE)
  }
}

# ---- 5. ggplot scale helpers -------------------------------------------------
scale_fill_subtype   <- function(...) scale_fill_manual(values = subtype_colours, ...)
scale_colour_subtype <- function(...) scale_colour_manual(values = subtype_colours, ...)
scale_color_subtype  <- scale_colour_subtype

# Locked categorical scales for the variables that previously drifted between
# figures. Use these rather than a local scale_*_manual() call.
scale_fill_plex      <- function(...) scale_fill_manual(values = plex_colours, ...)
scale_colour_plex    <- function(...) scale_colour_manual(values = plex_colours, ...)
scale_color_plex     <- scale_colour_plex
scale_fill_family    <- function(...) scale_fill_manual(values = family_colours, ...)
scale_colour_family  <- function(...) scale_colour_manual(values = family_colours, ...)
scale_color_family   <- scale_colour_family
scale_fill_stratum   <- function(...) scale_fill_manual(values = stratum_colours, ...)
scale_colour_stratum <- function(...) scale_colour_manual(values = stratum_colours, ...)
scale_color_stratum  <- scale_colour_stratum
scale_fill_tier      <- function(...) scale_fill_manual(values = tier_colours, ...)
scale_colour_tier    <- function(...) scale_colour_manual(values = tier_colours, ...)
scale_color_tier     <- scale_colour_tier
scale_fill_site      <- function(...) scale_fill_manual(values = site_colours, ...)
scale_colour_site    <- function(...) scale_colour_manual(values = site_colours, ...)
scale_color_site     <- scale_colour_site
scale_shape_site     <- function(...) scale_shape_manual(values = site_shapes, ...)

# Guard: catch a site/subtype/family label that has no locked colour BEFORE it
# renders as a silent NA grey. Call with the vector of labels you are about to plot.
check_palette_keys <- function(labels, palette, what = "value") {
  missing <- setdiff(unique(as.character(labels[!is.na(labels)])), names(palette))
  if (length(missing))
    stop(sprintf("no locked colour for %s: %s — add it to 00b_figure_theme.R",
                 what, paste(missing, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

# diverging: 0 (or `midpoint`) maps to the pale centre stop of the slate<->rust ramp
scale_fill_cook_div <- function(..., midpoint = 0, limits = NULL) {
  scale_fill_gradientn(colours = cook_diverging, limits = limits,
    rescaler = function(x, to = c(0,1), from = range(x, na.rm = TRUE))
      scales::rescale_mid(x, to, from, mid = midpoint), ...)
}
scale_colour_cook_div <- function(..., midpoint = 0, limits = NULL) {
  scale_colour_gradientn(colours = cook_diverging, limits = limits,
    rescaler = function(x, to = c(0,1), from = range(x, na.rm = TRUE))
      scales::rescale_mid(x, to, from, mid = midpoint), ...)
}
# sequential (non-centered continuous)
scale_fill_cook_seq   <- function(...) scale_fill_gradientn(colours = cook_sequential, ...)
scale_colour_cook_seq <- function(...) scale_colour_gradientn(colours = cook_sequential, ...)

# ---- 6. ComplexHeatmap / circlize colour functions --------------------------
# use for oncoprint fills, CNV heatmaps, z-scored expression heatmaps
cook_div_colfun <- function(max_abs = 2)
  circlize::colorRamp2(seq(-max_abs, max_abs, length.out = 7), cook_diverging)
cook_seq_colfun <- function(limits = c(0, 1))
  circlize::colorRamp2(seq(limits[1], limits[2], length.out = 8), cook_sequential)
# ComplexHeatmap global style: light hairline borders, standard sans serif.
ht_opt_cook <- function() {
  ComplexHeatmap::ht_opt(
    heatmap_border = TRUE,
    legend_border  = FALSE,
    simple_anno_size = grid::unit(3, "mm"),
    legend_title_gp = grid::gpar(fontsize = 8, fontfamily = FIG_FONT, col = cook_ink),
    legend_labels_gp= grid::gpar(fontsize = 7, fontfamily = FIG_FONT, col = cook_ink)
  )
}

# ---- 7. Export helpers (manuscript column widths, inches) --------------------
W1  <- 3.46   # single column  (~88 mm)
W15 <- 5.51   # 1.5 column     (~140 mm)
W2  <- 7.20   # double column  (~183 mm)
save_fig <- function(p, file, w, h, dpi = 400) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  if (grepl("\\.pdf$", file)) {
    ggplot2::ggsave(file, p, width = w, height = h, device = grDevices::cairo_pdf)
  } else {
    ggplot2::ggsave(file, p, width = w, height = h, dpi = dpi,
                    device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
  }
  invisible(file)
}

# Journal exports suppress captions only. Every remaining label is controlled
# explicitly in the source; there is no text-length heuristic that can silently
# hide a scientific label. The flag is kept for legacy diagnostic render paths.
FIG_PLAIN <- !identical(Sys.getenv("OVCAN_FIG_PLAIN", unset = "1"), "0")
if (FIG_PLAIN) {
  .plain_strip <- function(x) {
    if (inherits(x, "patchwork")) {
      if (length(x$patches$plots)) x$patches$plots <- lapply(x$patches$plots, .plain_strip)
      if (!is.null(x$patches$annotation)) x$patches$annotation$caption <- NULL
      x$labels$caption <- NULL
      x <- x & ggplot2::theme(plot.caption = ggplot2::element_blank())
    } else if (inherits(x, "ggplot")) {
      x$labels$caption <- NULL
      x <- x + ggplot2::theme(plot.caption = ggplot2::element_blank())
    }
    x
  }
  .save_fig_branded <- save_fig
  save_fig <- function(p, file, w, h, dpi = 400) {
    .save_fig_branded(.plain_strip(p), file, w, h, dpi = dpi)
  }
}
message("00b_figure_theme.R loaded | font: ", FIG_FONT,
        " | diverging palette: slate-white-rust | captions external")
