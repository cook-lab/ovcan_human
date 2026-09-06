# Figure 4: filtering, external expression comparison, driver candidates and CNV.
# Figure redesign 2026-09-05. Methods/caveats: reports/figure_redesign_2026-09-05/legends_genomics.md
# Read-only figure builder: all scientific quantities come from validated outputs.
source("scripts/00_setup.R"); source("scripts/00b_figure_theme.R")
suppressPackageStartupMessages({library(tidyverse); library(patchwork); library(ggnewscale)})
ASSET <- file.path(PROJ, "reports", "assets")
th <- function() theme_ovcan(8.5) + theme(axis.text = element_text(size = 7.5),
  legend.text = element_text(size = 7.5), legend.title = element_text(size = 7.5),
  legend.key.size = unit(9, "pt"), plot.margin = margin(5, 5, 5, 5))
tiers <- read_csv(file.path(OUT, "wes_driver_tiers.csv"), show_col_types = FALSE)
loads <- read_csv(file.path(OUT, "wes_mutation_load.csv"), show_col_types = FALSE)
fam <- read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)
models <- fam %>% filter(has_wes_maf) %>%
  mutate(subtype = factor(subtype, levels = c("HGS", "CC", "EC", "LGS", "MC"))) %>%
  arrange(subtype, patient_id, cell_line)
stopifnot(nrow(models) == 23, n_distinct(models$patient_id) == 16)
N <- nrow(models); NP <- n_distinct(models$patient_id)

# A: point trajectories preserve a meaningful logarithmic baseline.
cas <- read_csv(file.path(OUT, "wes_filter_cascade.csv"), show_col_types = FALSE)
stopifnot(setequal(cas$cell_line, models$cell_line), sum(cas$coding) == sum(loads$n_coding),
          all(cas$raw >= cas$pass), all(cas$pass >= cas$coding))
stages <- c("raw", "pass", "coding")
wf <- bind_rows(cas %>% summarise(across(all_of(stages), sum)) %>% mutate(scope = "All 23 models"),
                cas %>% filter(cell_line == "OV2295") %>% select(all_of(stages)) %>% mutate(scope = "OV2295")) %>%
  pivot_longer(all_of(stages), names_to = "stage", values_to = "n") %>%
  mutate(stage = factor(stage, levels = stages), scope = factor(scope, levels = c("All 23 models", "OV2295")))
pA <- ggplot(wf, aes(stage, n, group = scope, colour = scope)) +
  geom_line(linewidth = .65) + geom_point(size = 2.1) +
  geom_text(aes(label = scales::comma(n)), nudge_y = .12, vjust = -.65, size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = c("All 23 models" = cook_slate, OV2295 = cook_rust), name = NULL) +
  scale_x_discrete(labels = c(raw = "Annotated\nrecords", pass = "MAF PASS", coding = "Coding\ncandidates"),
                   expand = expansion(add = .27)) +
  scale_y_log10(breaks = c(100, 1000, 10000, 100000, 1000000),
               labels = c("100", "1,000", "10,000", "100,000", "1,000,000"),
               limits = c(100, 4000000), expand = c(0, 0)) +
  labs(x = NULL, y = "Variant count") + th() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0,0,0,0),
        legend.key.width = unit(13, "pt"))

# B: complete non-self comparison distribution, best non-self and self.
sm <- read_csv(file.path(OUT, "external_selfmatch_margin.csv"), show_col_types = FALSE) %>%
  arrange(margin_self_minus_best_nonself) %>% mutate(cell_line = factor(cell_line, levels = cell_line))
rho <- read_csv(file.path(OUT, "external_depmap_spearman_all.csv"), show_col_types = FALSE)
stopifnot(all(table(rho$our_line[rho$our_line %in% sm$cell_line]) == 67), nrow(sm) == 5)
nonself <- rho %>% filter(our_line %in% sm$cell_line, !is_selfpair) %>%
  mutate(cell_line = factor(our_line, levels = levels(sm$cell_line)))
mark <- bind_rows(sm %>% transmute(cell_line, value = best_nonself_spearman, comparison = "Best non-self"),
                 sm %>% transmute(cell_line, value = self_spearman, comparison = "Self"))
pB <- ggplot(sm, aes(y = cell_line)) +
  geom_segment(data = nonself, aes(x = spearman_hvg, xend = spearman_hvg,
               y = as.numeric(cell_line)-.10, yend = as.numeric(cell_line)+.10),
               linewidth = .35, colour = cook_slate, alpha = .3) +
  geom_segment(aes(x = best_nonself_spearman, xend = self_spearman, yend = cell_line),
               colour = cook_rust, linewidth = .65) +
  geom_point(data = mark, aes(x = value, shape = comparison, fill = comparison),
             colour = cook_rust, size = 2.5, stroke = .65) +
  scale_shape_manual(values = c("Best non-self" = 21, Self = 21), name = NULL) +
  scale_fill_manual(values = c("Best non-self" = "white", Self = cook_rust), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1), expand = c(.01, .01)) +
  scale_y_discrete(expand = expansion(add = .65)) +
  labs(x = "Spearman correlation", y = NULL) + th() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0,0,0,0))
save_fig(pB, file.path(ASSET, "f_external_concordance.png"), w = 4.5, h = 2.3)

# C: marker shape = variant class; marker fill = best priority tier.
trunc <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site")
ct <- tiers %>% mutate(tier_n = as.integer(sub("Tier", "", tier)),
  vc = case_when(variant_classification %in% trunc ~ "Truncating",
    variant_classification %in% c("In_Frame_Del", "In_Frame_Ins") ~ "In-frame",
    variant_classification == "Missense_Mutation" ~ "Missense", TRUE ~ NA_character_)) %>%
  group_by(gene, cell_line) %>% summarise(tier = paste("Tier", min(tier_n)),
    vc = if (n_distinct(vc) > 1) "Multiple classes" else unique(vc), .groups = "drop")
stopifnot(!anyNA(ct$vc))
gp <- tiers %>% distinct(gene) %>% left_join(tiers %>% filter(tier %in% c("Tier1", "Tier2")) %>%
  distinct(gene, patient_id) %>% count(gene, name = "np"), by = "gene") %>%
  mutate(np = replace_na(np, 0L)) %>% arrange(desc(np), gene)
G <- nrow(gp); genes <- gp$gene
xm <- setNames(seq_len(N), models$cell_line); ym <- setNames(G:1, genes)
bg <- expand_grid(cell_line = models$cell_line, gene = genes) %>% mutate(x = xm[cell_line], y = ym[gene])
ct <- ct %>% mutate(x = xm[cell_line], y = ym[gene])
gp <- gp %>% mutate(y = ym[gene])
subgroups <- models %>% mutate(x = row_number()) %>% group_by(subtype) %>%
  summarise(xmin = min(x)-.47, xmax = max(x)+.47, x = mean(range(x)), .groups = "drop") %>%
  mutate(label = as.character(subtype))
# Four groups with multiple profiled models are directly bracketed; remaining columns stand alone.
family_groups <- models %>% mutate(x = row_number()) %>% group_by(patient_id) %>% filter(n() > 1) %>%
  summarise(xmin = min(x)-.4, xmax = max(x)+.4, x = mean(range(x)), label = as.character(first(family)), .groups = "drop")
ld <- loads %>% filter(cell_line %in% models$cell_line) %>% mutate(x = xm[cell_line])
BASE <- G + 2.8; TOP <- 1.8
ly <- function(x) BASE + TOP * (log10(x)-2)/(log10(1500)-2)
XR <- N + 1.4; PW <- 2.4
pC <- ggplot() +
  geom_tile(data = bg, aes(x, y), fill = "#F8FAFC", colour = "white", linewidth = .25, width = .96, height = .96) +
  geom_point(data = ct, aes(x, y, fill = tier, shape = vc), size = 2.5, colour = cook_slate, stroke = .35) +
  scale_fill_tier(name = "Priority tier", guide = guide_legend(order = 1, nrow = 1,
    override.aes = list(shape = 22, size = 2.7))) +
  scale_shape_manual(values = c(Truncating = 22, Missense = 21, `In-frame` = 24, `Multiple classes` = 23),
    name = "Variant class", guide = guide_legend(order = 2, nrow = 1, override.aes = list(fill = cook_slate, size = 2.4))) +
  ggnewscale::new_scale_fill() +
  geom_rect(data = subgroups, aes(xmin = xmin, xmax = xmax, ymin = G+.65, ymax = G+1.25, fill = subtype)) +
  scale_fill_subtype(guide = "none") +
  geom_text(data = subgroups, aes(x, y = G+.96, label = label), colour = "white", size = 2.5) +
  geom_segment(data = family_groups, aes(x = xmin, xend = xmax, y = G+1.62, yend = G+1.62), colour = cook_slate, linewidth = .45) +
  geom_text(data = family_groups, aes(x, y = G+2.05, label = label), colour = cook_ink, size = 2.5) +
  annotate("segment", x = .52, xend = N+.48, y = ly(c(100,500,1500)), yend = ly(c(100,500,1500)),
           colour = cook_hair, linewidth = .3) +
  geom_segment(data = ld, aes(x, xend = x, y = BASE, yend = ly(n_coding)), colour = cook_slate, linewidth = 1.1) +
  geom_point(data = ld, aes(x, y = ly(n_coding)), colour = cook_slate, size = 1.1) +
  annotate("text", x = N+.68, y = ly(c(100,500,1500)), label = c("100","500","1,500"), size = 2.35, hjust = 0, colour = cook_ink) +
  geom_segment(data = gp, aes(x = XR, xend = XR + PW*np/NP, y = y, yend = y), colour = cook_slate, linewidth = 2.1) +
  geom_text(data = gp, aes(x = XR + PW + .3, y, label = np), hjust = 0, size = 2.5, colour = cook_ink) +
  annotate("text", x = XR+PW/2, y = G+1.1, label = "Patients\n(T1-2)", size = 2.5, lineheight = .95, colour = cook_ink) +
  annotate("segment", x = XR, xend = XR+PW, y = .2, yend = .2, colour = cook_hair, linewidth = .3) +
  annotate("text", x = c(XR, XR+PW), y = -.25, label = c("0",NP), size = 2.35, colour = cook_ink) +
  scale_x_continuous(breaks = seq_len(N), labels = models$cell_line, expand = c(0,0)) +
  scale_y_continuous(breaks = c(seq_len(G), G+.95, G+1.95, BASE+TOP/2),
    labels = c(paste0("<i>",names(sort(ym)),"</i>"), "Histotype", "Patient", "Coding<br>candidates"), expand = c(0,0)) +
  coord_cartesian(xlim = c(.45, XR+PW+.85), ylim = c(-.5, BASE+TOP+.25), clip = "off") +
  labs(x = NULL, y = NULL) + th() +
  theme(axis.text.x = element_text(size = 7.0, angle = 90, hjust = 1, vjust = .5),
    axis.text.y = ggtext::element_markdown(size = 7.5, family = FIG_FONT), axis.line = element_blank(), axis.ticks = element_blank(),
    legend.position = "bottom", legend.box = "vertical", legend.justification = "left",
    legend.spacing.y = unit(1, "pt"), legend.margin = margin(1,1,1,1))
save_fig(pC, file.path(ASSET, "f_wes_oncoplot.png"), w = W2, h = 4.7)

# D: all 39 assessed autosomal arms, with a patient-level denominator.
arm <- read_csv(file.path(OUT, "wes_cnv_arm_freq_patient.csv"), show_col_types = FALSE)
arm_order <- paste0(rep(1:22, each=2), rep(c("p","q"),22))
na_arms <- c("13p","14p","15p","21p","22p")
nd <- unique(arm$n_patients); stopifnot(length(nd)==1, nd==11)
ad <- arm %>% filter(!arm %in% na_arms) %>% mutate(arm=factor(arm, levels=intersect(arm_order,arm))) %>%
  select(arm, Gain=pct_patients_gain, Loss=pct_patients_loss) %>%
  pivot_longer(c(Gain,Loss), names_to="direction",values_to="pct") %>%
  mutate(signed=ifelse(direction=="Gain",pct,-pct))
pD <- ggplot(ad,aes(arm,signed,fill=direction)) + geom_col(width=.85) +
  geom_hline(yintercept=0,colour=cook_ink,linewidth=.3) +
  scale_fill_manual(values=c(Gain=cook_rust,Loss=cook_slate),name=NULL) +
  scale_y_continuous(limits=c(-100,100), breaks=c(-100,-50,0,50,100),labels=abs,expand=c(0,0)) +
  labs(x=NULL,y=sprintf("HGSC patients\n(%% of %d)",nd)) + th() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7.2),
        legend.position="top",legend.justification="left",legend.margin=margin(0,0,0,0))
fig4 <- ((pA | pB) / pC / pD) + plot_layout(heights=c(1.75,4.1,1.55)) +
  plot_annotation(tag_levels="A") & theme(text=element_text(family=FIG_FONT),
    plot.tag=element_text(size=11,face="bold",colour=cook_ink))
save_fig(fig4,file.path(MSFIG,"fig4.pdf"),w=W2,h=8.2)
save_fig(fig4,file.path(MSFIG,"fig4.png"),w=W2,h=8.2)
message("Figure 4 complete: 23 WES models, 16 mutation patients, 11 HGSC CNV patients.")
