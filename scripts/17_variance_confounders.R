# =============================================================================
# Script: 17_variance_confounders.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Revision workstream — quantify how much of the transcriptomic /
#          proteomic structure is subtype (biology) vs source site (batch) vs
#          patient family, and supplement the lenient marker "top-2" rule with
#          real effect sizes. Four deliverables:
#            A. JOINT / adjusted-R2 model for the RNA PCs (peer-review §3.5).
#               The published rna_pc_confounder.csv reports only UNIVARIATE
#               marginal R2 per PC, which cannot separate subtype from site
#               because the two are correlated (all HGS are Mes-Masson). Here we
#               regress each top PC on subtype + site JOINTLY and split the
#               variance into unique-subtype / unique-site / shared (commonality
#               analysis), and add univariate adjusted R2 to match
#               prot_pc_confounder.csv.
#            B. variancePartition genome-wide decomposition of per-gene /
#               per-protein variance into subtype + source_site + family
#               (+ TMT plex for protein) as random effects.
#            C. Passage-sensitivity check (§4 #8): do the top RNA PCs track
#               culture passage, and how badly do RNA vs WES passages mismatch
#               within a line?
#            D. Marker effect sizes (§3.6): Cohen's d + AUC of each canonical
#               marker in its intended subtype vs all other lines, computed on
#               patient-representative lines (no subline pseudoreplication), to
#               sit BESIDE the top-2-of-6 rule (not replace it).
#          Plus figure F: variance-partition violins (RNA + protein) and an
#          optional passage-check figure.
# Author:  Cook Lab (analyst: Claude — rna-variance workstream)
# Date:    2026-07-23
# Phase:   Revision (peer-review response §3.5 / §3.6 / §4 #8)
# -----------------------------------------------------------------------------
# CONFOUND STATED PLAINLY: HGS has no cross-centre replication
# (all 15 HGS lines are Mes-Masson/CHUM; all MMMT are Huntsman/BC). Clear cell
# (n=7: 2 Mes-Masson + 5 Huntsman) supplies the largest within-histotype
# centre comparison; EC and MC also span centres with fewer models. Every
# subtype-vs-site number below inherits this limit; the joint model can only
# separate the two using the mixed-site subtypes (CC/EC/MC).
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(tidyverse); library(SummarizedExperiment)
  library(lme4); library(parallel)            # variancePartition fallback (see Part B)
  library(scico)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
set.seed(SEED)

SUBTYPE_LEVELS <- c("HGS", "CC", "EC", "MC", "MMMT", "SCCOHT")

# --- lab viz helpers (visualization skill) -----------------------------------
okabe_ito <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000")
theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(text = element_text(colour = "black"),
          plot.title = element_text(size = 14, margin = margin(b = 8), hjust = 0),
          plot.subtitle = element_text(size = 9, colour = "grey30", margin = margin(b = 6), hjust = 0),
          axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10), legend.text = element_text(size = 9),
          strip.text = element_text(size = 11, face = "bold"),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          legend.key = element_blank(), panel.background = element_blank(),
          panel.grid = element_blank(), strip.background = element_blank(),
          plot.margin = margin(8, 8, 8, 8))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGS, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(PROJ, "reports", "assets", paste0(name, ".png")),
         p, width = w, height = h, dpi = 200)
  message("  saved figs/", name, ".pdf  +  reports/assets/", name, ".png")
}
term_cols <- c(subtype = "#C2410C", `Source site` = "#0F172A",
               `Patient family` = "#0072B2", `TMT plex` = "#009E73",
               Residuals = "#9CA3AF")

# =============================================================================
# LOAD
# =============================================================================
vsd <- readRDS(file.path(OUT, "rna_vst.rds"))
v   <- assay(vsd)                                   # 22,544 genes (ENSG) x 31 lines
ann <- readRDS(file.path(OUT, "rna_dds.rds")) |> colData() |> as.data.frame()
ann$subtype <- factor(ann$subtype, levels = SUBTYPE_LEVELS)
pca <- readRDS(file.path(OUT, "rna_pca.rds"))
stopifnot(identical(rownames(pca$x), ann$cell_line), identical(colnames(v), ann$cell_line))

ensure_family_map()   # [integration revision] 15 writes this map; see 00_setup.R.
                      # 07/08/16 already guarded it; 17 read it bare, so a clean
                      # checkout entered here failed on a file no script had run yet.
fam <- readr::read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)
samp <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE)

# family + a clean 2-level source (CHUM vs BC Cancer) for the RNA lines
site3 <- factor(ann$site)                                   # 3-level, matches published file
site2 <- factor(ifelse(ann$site == "Mes-Masson", "CHUM", "BC Cancer"),
                levels = c("CHUM", "BC Cancer"))            # institution-level contrast
# NB: line_family_map$family is NA for singleton lines; patient_id is populated
# for ALL lines (family number for sublings, cell_line for singletons) and is the
# correct per-patient grouping (collapses sublines, keeps singletons distinct).
rna_fam <- tibble(cell_line = ann$cell_line) |>
  left_join(fam |> select(cell_line, patient_id, patient_representative), by = "cell_line")
patient <- factor(rna_fam$patient_id)
stopifnot(!anyNA(patient))
cat(sprintf("RNA set: %d lines | %d subtypes | site3 {%s} | %d patients (sublines collapsed)\n",
            length(ann$cell_line), nlevels(ann$subtype),
            paste(levels(site3), collapse = ", "), nlevels(patient)))

# =============================================================================
# A. JOINT / ADJUSTED-R2 MODEL FOR THE RNA PCs  ->  rna_pc_confounder_joint.csv
# =============================================================================
# Commonality analysis on lm(PC ~ subtype + site): decompose the joint R2 into
# variance UNIQUE to subtype, UNIQUE to site, and SHARED (attributable to
# either because they are correlated). unique_site is the honest "how much does
# source site add BEYOND subtype" number a reviewer wants.
r2  <- function(y, f) summary(lm(y ~ f))$r.squared
ar2 <- function(y, f) summary(lm(y ~ f))$adj.r.squared
r2j <- function(y, a, b) summary(lm(y ~ a + b))$r.squared
ar2j<- function(y, a, b) summary(lm(y ~ a + b))$adj.r.squared

# WHICH SITE VARIABLE IS THE HEADLINE. Two exist and they differ materially
# (PC1: r2_site 0.313 with 3 levels vs 0.227 with 2). The 3-level variable is the
# one used for every *_site / unique_site / shared column below, because it is the
# metadata's own labelling; the 2-level institution contrast is carried alongside
# as *_2lvl. Any quoted number must say which. HEADLINE = 3-level.
HEADLINE_SITE <- "site3 (3 levels: Mes-Masson / Huntsman / Huntsman-Vanderhyden)"

# PASSAGE is added here as a third term (item: passage must be first-class).
# rna_passage_check.csv already showed passage adds ~15% of PC1 variance BEYOND
# site (a suppression pattern), which is larger than the genome-wide subtype
# component - so leaving it out of the joint model understated the confounding.
# Parsed here (again, cheaply) so section A does not depend on section C order.
parse_passage <- function(x) {                      # "p69" -> 69 ; "p39-40" -> 39.5
  x <- gsub("[Pp]", "", as.character(x))
  vapply(strsplit(x, "-"), function(z) mean(as.numeric(z)), numeric(1))
}
passage_vec <- samp |> filter(rna_seq == "Y", provenance == "generated") |>
  transmute(cell_line, rna_p = parse_passage(rna_passage)) |>
  (\(d) d$rna_p[match(ann$cell_line, d$cell_line)])()
stopifnot("passage missing for an RNA line" = !anyNA(passage_vec))

r2j3 <- function(y, a, b, c) summary(lm(y ~ a + b + c))$r.squared

joint <- map_dfr(1:5, function(k) {
  y <- pca$x[, k]
  rs <- r2(y, ann$subtype); rt3 <- r2(y, site3); rt2 <- r2(y, site2)
  rp <- r2(y, passage_vec)
  rj <- r2j(y, ann$subtype, site3)
  rj3 <- r2j3(y, ann$subtype, site3, passage_vec)     # + passage
  # aliasing guard: joint model must be full-rank (no dropped/NA coefficients)
  fit <- lm(y ~ ann$subtype + site3)
  aliased <- any(is.na(coef(fit)))
  arj <- ar2j(y, ann$subtype, site3)
  tibble(PC = paste0("PC", k),
         var_pct       = round(100 * pca$sdev[k]^2 / sum(pca$sdev^2), 1),
         r2_subtype    = rs,               r2_site = rt3,       # marginal (match published file)
         r2_site_2lvl  = rt2,              r2_passage = rp,
         adjr2_subtype = ar2(y, ann$subtype), adjr2_site = ar2(y, site3),  # match prot file
         adjr2_site_2lvl = ar2(y, site2),  adjr2_passage = ar2(y, passage_vec),
         r2_joint      = rj,               adjr2_joint = arj,
         # RAW commonality components
         unique_subtype = rj - rt3,        # variance subtype adds beyond site
         unique_site    = rj - rs,         # variance site adds beyond subtype  <-- key
         shared         = rs + rt3 - rj,
         # ADJUSTED commonality components. With a 6-level subtype + 3-level site
         # (~8 parameters) on 31 observations and four histotypes at n<=3, the raw
         # R2 is optimistic: fitting a group of n=2 costs one parameter and buys
         # near-perfect fit. Report BOTH; the adjusted joint R2 is ~8 points lower.
         adj_unique_subtype = arj - ar2(y, site3),
         adj_unique_site    = arj - ar2(y, ann$subtype),
         adj_shared         = ar2(y, ann$subtype) + ar2(y, site3) - arj,
         # 3-term model with passage
         r2_joint3_with_passage = rj3,
         unique_passage_beyond_subtype_site = rj3 - rj,
         headline_site_variable = HEADLINE_SITE,
         joint_aliased  = aliased)
})

# --- PERMUTATION NULL for the commonality components -------------------------
# The reader cannot tell whether unique-subtype = 42% is large without knowing
# what this n attains by chance. Two nulls, each respecting the confound:
#   unique-subtype null : shuffle SUBTYPE labels WITHIN site  (so the site
#                         composition of each subtype is preserved)
#   unique-site null    : shuffle SITE labels WITHIN subtype  (so single-site
#                         subtypes contribute nothing, which is the honest limit)
N_PERM <- 1000L
perm_within <- function(labels, by) {               # shuffle labels within groups of `by`
  out <- labels
  for (g in unique(by)) { i <- which(by == g); out[i] <- labels[i[sample.int(length(i))]] }
  out
}
perm_null <- map_dfr(1:5, function(k) {
  y  <- pca$x[, k]
  rs <- r2(y, ann$subtype); rt3 <- r2(y, site3)
  us_obs <- r2j(y, ann$subtype, site3) - rt3
  ut_obs <- r2j(y, ann$subtype, site3) - rs
  us_null <- replicate(N_PERM, r2j(y, perm_within(ann$subtype, site3), site3) - rt3)
  ut_null <- replicate(N_PERM, r2j(y, ann$subtype, perm_within(site3, ann$subtype)) - rs)
  tibble(PC = paste0("PC", k), n_perm = N_PERM,
         unique_subtype_obs = us_obs,
         unique_subtype_null_mean = mean(us_null),
         unique_subtype_null_median = median(us_null),
         unique_subtype_null_p95 = quantile(us_null, 0.95, names = FALSE),
         unique_subtype_perm_p = (1 + sum(us_null >= us_obs)) / (N_PERM + 1),
         unique_site_obs = ut_obs,
         unique_site_null_mean = mean(ut_null),
         unique_site_null_median = median(ut_null),
         unique_site_null_p95 = quantile(ut_null, 0.95, names = FALSE),
         unique_site_perm_p = (1 + sum(ut_null >= ut_obs)) / (N_PERM + 1))
})
perm_null <- perm_null %>% mutate(
  analysis_unit = "31 models; descriptive sensitivity with related sublines",
  primary_inference = "Use sensitivity_patient_reps_pc_permutation.csv: independent-patient reduced-model permutations")
readr::write_csv(perm_null, file.path(OUT, "rna_pc_confounder_permutation.csv"))
joint <- joint |> left_join(perm_null |> select(PC, unique_subtype_perm_p, unique_site_perm_p),
                            by = "PC")
readr::write_csv(joint, file.path(OUT, "rna_pc_confounder_joint.csv"))

# sanity: marginal R2 must reproduce the published rna_pc_confounder.csv
pub <- readr::read_csv(file.path(OUT, "rna_pc_confounder.csv"), show_col_types = FALSE)
stopifnot(all(abs(joint$r2_subtype - pub$r2_subtype) < 1e-6),
          all(abs(joint$r2_site    - pub$r2_site)    < 1e-6))

# -----------------------------------------------------------------------------
# A2. The within-clear-cell site control, ON THE RIGHT AXES
# -----------------------------------------------------------------------------
# The published version subset 7 clear-cell rows out of the PCA fitted to ALL 31
# lines and regressed those GLOBAL PC scores on site. Those PCs are defined by
# BETWEEN-histotype variance, so a within-CC slice of them is a residual and
# site's R2 on it is structurally small (0.044/0.058/0.057 -> the "4-6%" figure).
# That is not a within-group variance decomposition. Both are computed here:
#   global_pc_subset : the OLD statistic, kept so the change is auditable
#   within_cc_pca    : the CORRECT statistic - PCA refitted on the 7 CC lines only
# Neither is a licence to say site explains little, because with 2 vs 5 lines this
# test has essentially no power. The null expectation for a 1-df/5-df split is
# E[R2] = 1/(n-1) = 1/6 = 0.167, so a within-CC site R2 near 0.17 is exactly what
# NO site effect looks like. A within-group permutation p is therefore computed
# and is the number to quote. Conclusion the numbers support: UNDERPOWERED,
# indistinguishable from chance - not "site explains only 4-6%".
mixed <- ann$subtype %in% c("CC", "EC", "MC")
cc    <- ann$subtype == "CC"
NULL_R2_1DF <- 1 / (sum(cc) - 1)                    # E[R2] for a 1-df/(n-2)-df split

# refit the PCA within clear cell only (same recipe as script 02: top-2000 HVG on
# VST, centred, unscaled - but HVGs re-selected WITHIN the group, which is the
# point: the axes must describe within-CC variance)
v_cc   <- v[, cc, drop = FALSE]
rv_cc  <- apply(v_cc, 1, var)
top_cc <- head(order(rv_cc, decreasing = TRUE), 2000)
pca_cc <- prcomp(t(v_cc[top_cc, ]), center = TRUE, scale. = FALSE)
site_cc <- droplevels(site2[cc])
perm_r2 <- function(y, f, nperm = 2000L) {          # exact-ish permutation p on R2
  obs <- r2(y, f)
  null <- replicate(nperm, r2(y, sample(f)))
  (1 + sum(null >= obs)) / (nperm + 1)
}
within_cc <- map_dfr(1:3, function(k) tibble(
  approach   = c("global_pc_subset (OLD - wrong axes)", "within_cc_pca (CORRECT)"),
  PC         = paste0("PC", k),
  n          = sum(cc),
  n_CHUM     = sum(cc & site2 == "CHUM"),
  n_BC       = sum(cc & site2 == "BC Cancer"),
  var_pct    = c(round(100 * pca$sdev[k]^2 / sum(pca$sdev^2), 1),
                 round(100 * pca_cc$sdev[k]^2 / sum(pca_cc$sdev^2), 1)),
  site_r2    = c(r2(pca$x[cc, k], site_cc), r2(pca_cc$x[, k], site_cc)),
  site_p_anova = c(anova(lm(pca$x[cc, k] ~ site_cc))$`Pr(>F)`[1],
                   anova(lm(pca_cc$x[, k] ~ site_cc))$`Pr(>F)`[1]),
  site_p_perm  = c(perm_r2(pca$x[cc, k], site_cc), perm_r2(pca_cc$x[, k], site_cc)),
  null_expected_r2 = NULL_R2_1DF))
readr::write_csv(within_cc, file.path(OUT, "rna_within_cc_site.csv"))

# the mixed-site set (CC/EC/MC) kept as a wider, still-underpowered companion
site_within <- tibble(
  set = c("mixed-site subtypes (CC/EC/MC)", "clear cell only (CC)"),
  n   = c(sum(mixed), sum(cc)),
  n_CHUM = c(sum(mixed & site2 == "CHUM"), sum(cc & site2 == "CHUM")),
  n_BC   = c(sum(mixed & site2 == "BC Cancer"), sum(cc & site2 == "BC Cancer")),
  PC1_site_r2 = c(r2(pca$x[mixed, 1], droplevels(site2[mixed])),
                  r2(pca$x[cc, 1],    droplevels(site2[cc]))),
  PC2_site_r2 = c(r2(pca$x[mixed, 2], droplevels(site2[mixed])),
                  r2(pca$x[cc, 2],    droplevels(site2[cc]))))

cat("\n=== A. RNA PC confounder — joint / commonality decomposition ===\n")
cat("Headline site variable: ", HEADLINE_SITE, "\n", sep = "")
print(as.data.frame(joint |> select(PC, var_pct, r2_subtype, r2_site, r2_site_2lvl,
                                    r2_joint, adjr2_joint, unique_subtype, unique_site, shared) |>
                      mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat("\nADJUSTED commonality components (report these alongside the raw ones):\n")
print(as.data.frame(joint |> select(PC, adjr2_subtype, adjr2_site,
                                    adj_unique_subtype, adj_unique_site, adj_shared) |>
                      mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat("\nPASSAGE added to the joint model:\n")
print(as.data.frame(joint |> select(PC, r2_passage, r2_joint, r2_joint3_with_passage,
                                    unique_passage_beyond_subtype_site) |>
                      mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat(sprintf("\nPERMUTATION NULL for the commonality components (%d draws each):\n", N_PERM))
print(as.data.frame(perm_null |> select(PC, unique_subtype_obs, unique_subtype_null_mean,
                                        unique_subtype_null_p95, unique_subtype_perm_p,
                                        unique_site_obs, unique_site_null_mean,
                                        unique_site_null_p95, unique_site_perm_p) |>
                      mutate(across(where(is.numeric), ~signif(.x, 3)))), row.names = FALSE)
cat("\n!! CHECK the 'unique-site <=0.2% across PC1-PC3' claim against this table:\n")
for (k in 1:3)
  cat(sprintf("   %s unique_site = %.4f (%.2f%%)%s\n", joint$PC[k], joint$unique_site[k],
              100 * joint$unique_site[k],
              ifelse(joint$unique_site[k] > 0.002, "   <-- EXCEEDS 0.2%", "")))

cat("\n=== A2. Within-clear-cell site control (2 CHUM vs 5 BC Cancer) ===\n")
print(as.data.frame(within_cc |> mutate(across(where(is.numeric), ~signif(.x, 3)))), row.names = FALSE)
cat(sprintf("Null expectation for a 1-df/%d-df split: E[R2] = 1/%d = %.3f.\n",
            sum(cc) - 2, sum(cc) - 1, NULL_R2_1DF))
cat("The correct (within-CC) R2 values sit AT or near that null and no permutation p is\n")
cat("significant -> the honest conclusion is UNDERPOWERED / indistinguishable from chance,\n")
cat("NOT 'site explains only 4-6%'. Do not quote the global_pc_subset row.\n")
cat("\nWider mixed-site companion (also underpowered):\n")
print(as.data.frame(site_within |> mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)
if (any(joint$joint_aliased)) cat("!! WARNING: joint model aliased on some PC (site nested in subtype)\n")

# =============================================================================
# B. variancePartition  ->  rna_variancepartition.csv, prot_variancepartition.csv
# =============================================================================
# variancePartition IS installed but broken in this R: it calls lme4::findbars,
# which has moved to the 'reformulas' package -> "Initial model failed". We use
# the DOCUMENTED FALLBACK — the identical per-feature REML variance decomposition
# via lme4 directly: variance components (VarCorr) / total variance, parallelised.
# This reproduces variancePartition's fraction-of-variance-explained.
ctrl <- lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE, optimizer = "bobyqa")
# `fixed_label`: name to give the variance explained by the FIXED part of the model
# (used for continuous covariates such as passage, which cannot be a random effect
# with one observation per level). Computed as var(fitted values with the random
# effects zeroed) — the same quantity variancePartition attributes to fixed terms.
# `dropped` is returned as an attribute: features whose lmer fit errored were
# previously discarded silently, so n_features quietly fell short of the gene count
# the manuscript and figure captions quote.
vp_lmer <- function(mat, meta, form, ncores = 4, fixed_label = NULL) {
  fit_one <- function(y) {
    d <- data.frame(y = y, meta)
    out <- tryCatch({
      fit <- suppressMessages(suppressWarnings(lmer(form, data = d, REML = TRUE, control = ctrl)))
      vc  <- as.data.frame(VarCorr(fit)); v <- setNames(vc$vcov, vc$grp)
      if (!is.null(fixed_label)) {
        fx <- stats::var(as.numeric(stats::predict(fit, re.form = NA)))
        v  <- c(setNames(fx, fixed_label), v)
      }
      v
    }, error = function(e) NULL)
    if (is.null(out)) NULL else out / sum(out)
  }
  rows <- mclapply(seq_len(nrow(mat)), function(i) fit_one(mat[i, ]), mc.cores = ncores)
  keep <- !vapply(rows, is.null, logical(1))
  dropped <- rownames(mat)[!keep]
  M <- do.call(rbind, rows[keep]); rownames(M) <- rownames(mat)[keep]
  ord <- c("subtype", "source_site", "plex", "patient", "family", "passage", "Residual")
  M <- M[, intersect(c(ord, setdiff(colnames(M), ord)), colnames(M)), drop = FALSE]
  colnames(M)[colnames(M) == "Residual"] <- "Residuals"
  out <- as.data.frame(M)
  attr(out, "dropped")   <- dropped
  attr(out, "n_input")   <- nrow(mat)
  out
}
# MEAN AND IQR TRAVEL WITH THE MEDIAN. The per-gene distributions are strongly
# right-skewed and the median tells a different story from the mean (RNA subtype
# median 5.9% vs site 3.5%, but means 14.9% vs 14.5% with overlapping IQRs), so a
# median-only report is not defensible. All four are written for every term.
vp_summary <- function(df, assay_label, n_samples, model_label = "published",
                       design_note = NA_character_) {
  tibble(assay = assay_label, model = model_label, term = colnames(df),
         median_pct = 100 * apply(df, 2, median),
         mean_pct   = 100 * apply(df, 2, mean),
         q25_pct    = 100 * apply(df, 2, quantile, 0.25),
         q75_pct    = 100 * apply(df, 2, quantile, 0.75),
         n_features = nrow(df),
         n_features_input   = attr(df, "n_input"),
         n_features_dropped = length(attr(df, "dropped")),
         n_samples  = n_samples,
         design     = design_note)
}

# ---- DESIGN / IDENTIFIABILITY, recorded so it travels with the numbers -------
# `patient` is close to unidentifiable here: 31 lines, 28 patient levels, only 3
# patients replicated (one pair each), so the patient variance component and the
# residual are separable ONLY from those 3 pairs. A term with 28 levels on 31
# observations absorbs essentially all reproducible variance regardless of the
# biology, so "line/patient identity is the largest structured term" restates the
# design rather than reporting a finding. `source_site` has 3 levels and one of
# them is n=1 (BIN67), so its component also carries very wide implicit
# uncertainty. Two sensitivity models are therefore fitted and written alongside:
#   patient_dropped : the same model without `patient`
# A former sensitivity pooled unrelated singleton patients into one random-effect
# group; it has been retired because that invents a shared covariance. A second
# valid descriptive sensitivity adds passage as a fixed covariate (section C).
rep_patients <- names(which(table(patient) > 1))
site_counts  <- table(site3)
design_note_rna <- sprintf("31 lines, %d patient levels, %d replicated patients (%s); site levels %s",
                           nlevels(patient), length(rep_patients),
                           paste(rep_patients, collapse = "/"),
                           paste(sprintf("%s=%d", names(site_counts), as.integer(site_counts)),
                                 collapse = ", "))
cat("\n=== B. Design / identifiability ===\n"); cat(design_note_rna, "\n")

# ---- RNA: all expressed genes; subtype + source_site + patient -------------
rna_meta <- data.frame(row.names = ann$cell_line,
                       subtype = ann$subtype, source_site = site3, patient = patient,
                       passage_z = as.numeric(scale(passage_vec)))
cat(sprintf("\n=== B. lmer variance decomposition (RNA): %d genes x %d lines (4 cores) ===\n", nrow(v), ncol(v)))
vp_rna <- vp_lmer(v, rna_meta, y ~ (1 | subtype) + (1 | source_site) + (1 | patient))
rna_vp_tab <- vp_summary(vp_rna, "RNA", ncol(v), "published", design_note_rna)
print(as.data.frame(rna_vp_tab |> select(term, median_pct, mean_pct, q25_pct, q75_pct,
                                         n_features, n_features_dropped) |>
                      mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)
if (length(attr(vp_rna, "dropped"))) {
  readr::write_lines(attr(vp_rna, "dropped"),
                     file.path(OUT, "rna_variancepartition_dropped_genes.txt"))
  cat(sprintf("!! %d of %d genes were DROPPED (lmer fit failed) and are listed in\n",
              length(attr(vp_rna, "dropped")), attr(vp_rna, "n_input")))
  cat("   output/rna_variancepartition_dropped_genes.txt. n_features is the number to quote:\n")
  cat(sprintf("   %d, NOT %d.\n", nrow(vp_rna), attr(vp_rna, "n_input")))
  cat("   Dropped: ", paste(attr(vp_rna, "dropped"), collapse = ", "), "\n")
}

# ---- RNA sensitivity models -------------------------------------------------
cat("\n--- RNA sensitivity: patient dropped / + passage ---\n")
vp_rna_nopat <- vp_lmer(v, rna_meta, y ~ (1 | subtype) + (1 | source_site))
# Removed the former pooled-singleton random effect. Unrelated singleton donors
# cannot share one random intercept: that model invents correlation among them.
# Patient-dropped and passage-adjusted models remain descriptive sensitivities.
vp_rna_pass  <- vp_lmer(v, rna_meta,
                        y ~ passage_z + (1 | subtype) + (1 | source_site) + (1 | patient),
                        fixed_label = "passage")
rna_vp_sens <- bind_rows(
  vp_summary(vp_rna_nopat, "RNA", ncol(v), "patient_dropped", design_note_rna),
  vp_summary(vp_rna_pass,  "RNA", ncol(v), "published_plus_passage_fixed",
             paste0(design_note_rna, "; passage added as a scaled FIXED covariate")))
print(as.data.frame(rna_vp_sens |> select(model, term, median_pct, mean_pct, q25_pct, q75_pct) |>
                      mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)

# ---- Protein: complete-case proteins; + TMT plex axis ----------------------
pobj <- readRDS(file.path(OUT, "prot_matrix.rds"))
pmat <- pobj$mat; pann <- pobj$sample_ann
stopifnot(identical(colnames(pmat), pann$cell_line))
patient_prot <- factor(fam$patient_id[match(pann$cell_line, fam$cell_line)])
prot_meta <- data.frame(
  row.names   = pann$cell_line,
  subtype     = factor(pann$subtype_ss, levels = SUBTYPE_LEVELS),
  source_site = factor(pann$site_ss),
  plex        = factor(pann$plex),
  patient     = patient_prot)
cc_prot <- rowSums(is.na(pmat)) == 0
pmcc <- pmat[cc_prot, , drop = FALSE]
cat(sprintf("\n=== B. lmer variance decomposition (protein): %d complete-case proteins (of %d) x %d lines ===\n",
            nrow(pmcc), nrow(pmat), ncol(pmcc)))
cat(sprintf("   patients: %d | TMT plexes: %d\n", nlevels(prot_meta$patient), nlevels(prot_meta$plex)))
rep_patients_prot <- names(which(table(prot_meta$patient) > 1))
design_note_prot  <- sprintf("%d lines, %d patient levels, %d replicated patients (%s); site levels %s; %d TMT plexes",
                             ncol(pmcc), nlevels(prot_meta$patient), length(rep_patients_prot),
                             paste(rep_patients_prot, collapse = "/"),
                             paste(sprintf("%s=%d", names(table(prot_meta$source_site)),
                                           as.integer(table(prot_meta$source_site))), collapse = ", "),
                             nlevels(prot_meta$plex))
cat("   ", design_note_prot, "\n", sep = "")
vp_prot <- vp_lmer(pmcc, prot_meta, y ~ (1 | subtype) + (1 | source_site) + (1 | plex) + (1 | patient))
prot_vp_tab <- vp_summary(vp_prot, "protein", ncol(pmcc), "published", design_note_prot)
print(as.data.frame(prot_vp_tab |> select(term, median_pct, mean_pct, q25_pct, q75_pct,
                                          n_features, n_features_dropped) |>
                      mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)
if (length(attr(vp_prot, "dropped")))
  cat(sprintf("!! %d of %d proteins dropped (lmer fit failed); quote n_features = %d\n",
              length(attr(vp_prot, "dropped")), attr(vp_prot, "n_input"), nrow(vp_prot)))

cat("\n--- Protein sensitivity: patient dropped ---\n")
vp_prot_nopat <- vp_lmer(pmcc, prot_meta, y ~ (1 | subtype) + (1 | source_site) + (1 | plex))
prot_vp_sens <- vp_summary(vp_prot_nopat, "protein", ncol(pmcc), "patient_dropped", design_note_prot)
print(as.data.frame(prot_vp_sens |> select(model, term, median_pct, mean_pct, q25_pct, q75_pct) |>
                      mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)

# The PRIMARY model stays in the published filenames (Fig 3D reads these); the
# sensitivity fits go to a separate file so the primary table keeps its shape.
readr::write_csv(rna_vp_tab,  file.path(OUT, "rna_variancepartition.csv"))
readr::write_csv(prot_vp_tab, file.path(OUT, "prot_variancepartition.csv"))
readr::write_csv(bind_rows(rna_vp_sens, prot_vp_sens),
                 file.path(OUT, "variancepartition_sensitivity.csv"))

# =============================================================================
# C. PASSAGE-SENSITIVITY CHECK  ->  rna_passage_check.csv
# =============================================================================
rna_pass <- samp |> filter(rna_seq == "Y", provenance == "generated") |>
  transmute(cell_line, rna_passage, rna_p = parse_passage(rna_passage))
rna_pass <- tibble(cell_line = ann$cell_line, subtype = ann$subtype, site2 = site2) |>
  left_join(rna_pass, by = "cell_line")
stopifnot(!any(is.na(rna_pass$rna_p)))

# (a) does passage drive the top PCs? PC ~ passage R2/p. Also test the passage-
#     site confound: passage itself is almost perfectly split by site.
pass_pc <- map_dfr(1:5, function(k) {
  y <- pca$x[, k]; m <- lm(y ~ rna_pass$rna_p); f <- summary(m)$fstatistic
  # partial: does passage still explain the PC AFTER removing site?
  r_full <- summary(lm(y ~ site2 + rna_pass$rna_p))$r.squared
  r_site <- summary(lm(y ~ site2))$r.squared
  tibble(PC = paste0("PC", k),
         var_pct = round(100 * pca$sdev[k]^2 / sum(pca$sdev^2), 1),
         r2_passage = summary(m)$r.squared,
         p_passage  = pf(f[1], f[2], f[3], lower.tail = FALSE),
         r2_site    = r_site,
         partial_r2_passage_after_site = r_full - r_site)   # passage's unique add beyond site
})
# passage-vs-site confound magnitude
pass_site_r2 <- summary(lm(rna_pass$rna_p ~ site2))$r.squared
cat(sprintf("\n=== C. Passage check — passage is %.0f%% explained by site (CHUM high-passage vs BC low-passage) ===\n",
            100 * pass_site_r2))
print(as.data.frame(pass_pc |> mutate(across(where(is.numeric), ~signif(.x, 3)))), row.names = FALSE)

# (b) within-line cross-assay passage mismatch (RNA vs WES)
xassay <- samp |>
  filter(provenance == "generated",
         !is.na(rna_passage), rna_passage != "",
         !is.na(wes_passage), wes_passage != "") |>
  transmute(cell_line, subtype,
            rna_passage, wes_passage,
            rna_p = parse_passage(rna_passage),
            wes_p = parse_passage(wes_passage),
            passage_diff = wes_p - rna_p,
            abs_diff = abs(wes_p - rna_p)) |>
  arrange(desc(abs_diff))
cat(sprintf("\nCross-assay RNA-vs-WES passage mismatch (%d lines with both):\n", nrow(xassay)))
print(as.data.frame(xassay), row.names = FALSE)
cat(sprintf("  median |diff| = %.1f passages; max = %.0f (%s)\n",
            median(xassay$abs_diff), max(xassay$abs_diff),
            xassay$cell_line[which.max(xassay$abs_diff)]))

# Dedicated discordance table. Every cross-assay analysis in the resource (the
# RNA-protein concordance, the RNA-vs-WES identity checks, the ADC atlas) assumes
# the layers describe the same cells; for some lines they are ~20 passages apart,
# and passage explains a non-trivial share of PC1 (see pass_pc above). The range
# has to be quotable, so it is written out with the per-line detail.
disc <- xassay |>
  transmute(cell_line, subtype = as.character(subtype),
            rna_passage_label = rna_passage, wes_passage_label = wes_passage,
            rna_passage = rna_p, wes_passage = wes_p,
            passage_diff = passage_diff, abs_diff,
            wes_is_later = passage_diff > 0)
disc_summary <- tibble(
  statistic = c("n_lines_with_both", "median_abs_diff", "mean_abs_diff",
                "min_diff", "max_diff", "max_abs_diff", "n_lines_abs_diff_ge_10",
                "n_lines_abs_diff_ge_20", "spearman_rna_vs_wes_passage"),
  value = c(nrow(disc), median(disc$abs_diff), mean(disc$abs_diff),
            min(disc$passage_diff), max(disc$passage_diff), max(disc$abs_diff),
            sum(disc$abs_diff >= 10), sum(disc$abs_diff >= 20),
            cor(disc$rna_passage, disc$wes_passage, method = "spearman")))
readr::write_csv(bind_rows(
  disc |> mutate(block = "per line", .before = 1),
  disc_summary |> transmute(block = "summary", cell_line = statistic,
                            passage_diff = value)),
  file.path(OUT, "rna_passage_discordance.csv"))
cat("\nRNA-vs-WES passage discordance summary (output/rna_passage_discordance.csv):\n")
print(as.data.frame(disc_summary |> mutate(value = round(value, 3))), row.names = FALSE)
cat(sprintf("  QUOTABLE RANGE: RNA-vs-WES passage differs by %.0f to %+.0f passages across %d lines\n",
            min(disc$passage_diff), max(disc$passage_diff), nrow(disc)))

# one tidy CSV, two clearly-typed blocks (a `check` column disambiguates):
#   PC~passage regression  : id=PCk, r2_passage/p_passage/r2_site/partial add
#   cross-assay mismatch   : id=cell_line, rna_passage/wes_passage/passage_diff
pass_out <- bind_rows(
  pass_pc |> transmute(check = "PC~passage regression", id = PC, var_pct,
                       r2_passage, p_passage, r2_site, partial_r2_passage_after_site,
                       rna_passage = NA_real_, wes_passage = NA_real_, passage_diff = NA_real_),
  xassay |> transmute(check = "cross-assay passage mismatch", id = cell_line, var_pct = NA_real_,
                      r2_passage = NA_real_, p_passage = NA_real_, r2_site = NA_real_,
                      partial_r2_passage_after_site = NA_real_,
                      rna_passage = rna_p, wes_passage = wes_p, passage_diff))
readr::write_csv(pass_out, file.path(OUT, "rna_passage_check.csv"))

# =============================================================================
# D. MARKER EFFECT SIZES  ->  rna_marker_effectsizes.csv
# =============================================================================
# Collapse to ONE line per patient family (no subline pseudoreplication), then
# for each canonical marker compute Cohen's d and AUC of intended-subtype vs the
# rest, on VST. AUC is rank-based (robust for the n=2 rare subtypes); Cohen's d
# is on the VST scale. Direction-aware: 'up' markers should score d>0 / AUC>0.5,
# 'down' markers (SCCOHT SWI/SNF) d<0 / AUC<0.5 — we also report an ORIENTED AUC
# where 1 = perfect discrimination in the expected direction for every marker.
# For each patient, keep the flagged patient_representative among its RNA lines;
# if none of a patient's RNA lines is flagged, keep the first (deterministic).
rep_keep <- rna_fam |>
  group_by(patient_id) |>
  summarise(cell_line = if (any(patient_representative)) cell_line[patient_representative][1]
                        else cell_line[1], .groups = "drop") |>
  pull(cell_line)
dropped <- setdiff(ann$cell_line, rep_keep)
cat(sprintf("\n=== D. Marker effect sizes — collapsed to %d patient representatives (dropped sublines: %s) ===\n",
            length(rep_keep), paste(dropped, collapse = ", ")))
stopifnot(length(rep_keep) == 28)

# -----------------------------------------------------------------------------
# MARKER QUANTIFICATION IS NOW HARMONISED WITH SCRIPT 04 (was inconsistent).
# This script previously scored markers on VST, collapsing a multi-ENSG symbol by
# picking the single ENSG with the highest mean VST, while script 04 scored the
# SAME markers on symbol-SUMMED TPM. The two are different quantities, and yet the
# `lands_right` column from 04 was joined onto the effect sizes from here, so the
# published table mixed them. ONE RULE, used in both places:
#     collapse gene_id -> symbol by SUMMING TPM (TPM is additive, so summing
#     recovers the symbol-level abundance), then log2(TPM + 1).
# This is script 04's rule and also the rule used by scripts 12/13/19, so the
# whole resource now quantifies a gene symbol the same way. Consequence: Cohen's d
# changes (it is scale-dependent); AUC is rank-based and moves only where the two
# collapses reorder lines.
tpm <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"), show_col_types = FALSE)
id2sym_v <- t2g |> distinct(ensembl_gene_id, external_gene_name) |>
  filter(!is.na(external_gene_name), external_gene_name != "") |>
  (\(d) setNames(d$external_gene_name, d$ensembl_gene_id))()
m_id    <- as.matrix(tpm[, -1]); rownames(m_id) <- tpm$gene_id
common  <- intersect(rownames(m_id), names(id2sym_v))
logtpm  <- log2(rowsum(m_id[common, , drop = FALSE], group = id2sym_v[common]) + 1)
stopifnot("marker matrix lost the representative lines" = all(rep_keep %in% colnames(logtpm)))
x_rep   <- logtpm[, rep_keep, drop = FALSE]         # symbols x 28 patient reps
sub_rep <- ann$subtype[match(rep_keep, ann$cell_line)]

# EC markers (Hollis 2020) provide expression context for the two annotated
# EC models, including the discordant SWI/SNF evidence for TOV112D. Same set and
# same directions as script 04; VIM is scored once, under MMMT (see 04's note).
markers <- tribble(
  ~marks_code, ~symbol,   ~direction,
  "HGS","PAX8","up", "HGS","WT1","up", "HGS","MUC16","up", "HGS","MECOM","up", "HGS","SOX17","up",
  "CC","HNF1B","up", "CC","NAPSA","up", "CC","SPP1","up", "CC","GPX3","up", "CC","GCLC","up",
  "EC","ESR1","up", "EC","PGR","up", "EC","ARID1A","down",
  "MC","CDX2","up", "MC","MUC2","up", "MC","TFF1","up", "MC","TFF3","up", "MC","KRT20","up", "MC","MUC5AC","up",
  "MMMT","VIM","up", "MMMT","ZEB1","up", "MMMT","CDH2","up", "MMMT","SNAI2","up",
  "SCCOHT","SMARCA4","down", "SCCOHT","SMARCA2","down",
  NA_character_,"MKI67","up")
stopifnot("marker symbols must be unique" = !any(duplicated(markers$symbol)),
          "marker(s) absent from the symbol matrix" = all(markers$symbol %in% rownames(x_rep)))

get_expr <- function(sym) x_rep[sym, ]              # symbol-summed log2(TPM+1)
auc_fun <- function(x_in, x_rest) {                 # P(intended > rest), tie=0.5
  (sum(outer(x_in, x_rest, ">")) + 0.5 * sum(outer(x_in, x_rest, "=="))) /
    (length(x_in) * length(x_rest))
}
cohens_d_fun <- function(x_in, x_rest) {
  n_in <- length(x_in); n_rest <- length(x_rest)
  s_pool <- sqrt(((n_in - 1) * var(x_in) + (n_rest - 1) * var(x_rest)) / (n_in + n_rest - 2))
  if (!is.finite(s_pool) || s_pool == 0) return(NA_real_)
  (mean(x_in) - mean(x_rest)) / s_pool
}
# BOOTSTRAP CIs. The largest effects come from the smallest groups (SMARCA2 on
# n=2, KRT20 on n=3), so point estimates are not comparable to one another and a
# reader cannot see that without an interval. Stratified resampling: bootstrap
# WITHIN the intended group and WITHIN the rest, preserving both n's. With n=2 the
# intended group has only 3 distinct resamples, so many replicates are degenerate
# (zero within-group variance) — that is the honest signal that the estimate is
# unstable, and n_boot_valid records how often a finite value was obtainable.
N_BOOT <- 2000L
boot_ci <- function(x_in, x_rest, direction, stat = c("d", "auc"), conf = 0.95) {
  stat <- match.arg(stat)
  vals <- replicate(N_BOOT, {
    a <- sample(x_in, length(x_in), replace = TRUE)
    b <- sample(x_rest, length(x_rest), replace = TRUE)
    if (stat == "d") cohens_d_fun(a, b)
    else { u <- auc_fun(a, b); if (direction == "down") 1 - u else u }
  })
  ok <- is.finite(vals)
  q  <- if (any(ok)) quantile(vals[ok], c((1 - conf) / 2, 1 - (1 - conf) / 2), names = FALSE)
        else c(NA_real_, NA_real_)
  c(lo = q[1], hi = q[2], n_valid = sum(ok))
}
effsize_one <- function(marks_code, symbol, direction) {
  x <- get_expr(symbol)
  if (is.na(marks_code))                            # MKI67 proliferation control: unscored
    return(tibble(marks_code, symbol, direction,
                  mean_intended = NA_real_, mean_rest = NA_real_,
                  n_intended = NA_integer_, n_rest = NA_integer_,
                  cohens_d = NA_real_, cohens_d_lo = NA_real_, cohens_d_hi = NA_real_,
                  auc = NA_real_, auc_oriented = NA_real_,
                  auc_oriented_lo = NA_real_, auc_oriented_hi = NA_real_,
                  n_boot_valid = NA_integer_, wilcox_p = NA_real_, wilcox_p_min_attainable = NA_real_))
  inx <- sub_rep == marks_code
  x_in <- x[inx]; x_rest <- x[!inx]
  n_in <- sum(inx); n_rest <- sum(!inx)
  a  <- auc_fun(x_in, x_rest)
  cd <- boot_ci(x_in, x_rest, direction, "d")
  ca <- boot_ci(x_in, x_rest, direction, "auc")
  # smallest two-sided p the exact rank-sum test can return at these group sizes:
  # complete separation. For n_in=2 vs 26 that is ~0.0055, so several p-values are
  # FLOOR-LIMITED and cannot be compared with p-values from larger groups.
  p_min <- suppressWarnings(wilcox.test(seq_len(n_in), n_in + seq_len(n_rest))$p.value)
  tibble(marks_code, symbol, direction,
         mean_intended = mean(x_in), mean_rest = mean(x_rest),
         n_intended = n_in, n_rest = n_rest,
         cohens_d = cohens_d_fun(x_in, x_rest),
         cohens_d_lo = cd[["lo"]], cohens_d_hi = cd[["hi"]],
         auc = a, auc_oriented = if (direction == "down") 1 - a else a,
         auc_oriented_lo = ca[["lo"]], auc_oriented_hi = ca[["hi"]],
         n_boot_valid = as.integer(min(cd[["n_valid"]], ca[["n_valid"]])),
         wilcox_p = suppressWarnings(wilcox.test(x_in, x_rest)$p.value),
         wilcox_p_min_attainable = p_min)
}
eff <- pmap_dfr(markers, effsize_one)

# attach the published top-2 rule so effect sizes sit BESIDE it. Both sides are
# now computed on the SAME quantity (symbol-summed log2 TPM), so the join is valid.
# Recompute ranking on the SAME 28 representatives as the effect sizes. The
# previous join imported 31-model ranks and silently mixed the analysis units.
top2 <- purrr::pmap_dfr(markers, function(marks_code, symbol, direction) {
  vals <- tapply(x_rep[symbol, ], sub_rep, mean)
  rank_in <- if (is.na(marks_code)) NA_integer_ else
    as.integer(rank(-vals, ties.method = "min")[marks_code])
  recovered <- if (is.na(marks_code)) NA else if (direction == "up")
    rank_in <= 2L && vals[marks_code] > 1 else rank_in >= 5L
  tibble(symbol, rank_in, top_subtype = names(vals)[which.max(vals)], lands_right = recovered)
})
marker_eff <- eff |>
  left_join(top2, by = "symbol") |>
  # BH across the graded markers only (MKI67 carries no test). These were 22 (now
  # 25) unadjusted p-values in the published table.
  mutate(wilcox_p_bh = p.adjust(wilcox_p, method = "BH"),
         wilcox_p_floor_limited = is.finite(wilcox_p) & is.finite(wilcox_p_min_attainable) &
                                  wilcox_p <= wilcox_p_min_attainable * 1.001,
         quantification = "symbol-summed log2(TPM+1), 28 patient representatives") |>
  mutate(across(c(mean_intended, mean_rest, cohens_d, cohens_d_lo, cohens_d_hi,
                  auc, auc_oriented, auc_oriented_lo, auc_oriented_hi), ~round(.x, 3)),
         across(c(wilcox_p, wilcox_p_bh, wilcox_p_min_attainable), ~signif(.x, 3)))
stopifnot("BH adjustment must be computed on the graded markers" =
            sum(!is.na(marker_eff$wilcox_p_bh)) == sum(!is.na(marker_eff$marks_code)))
readr::write_csv(marker_eff, file.path(OUT, "rna_marker_effectsizes.csv"))

cat("\nMarker effect sizes (intended subtype vs rest, symbol-summed log2 TPM, 28 patient reps):\n")
print(as.data.frame(marker_eff |>
        select(marks_code, symbol, direction, n_intended, cohens_d, cohens_d_lo, cohens_d_hi,
               auc_oriented, auc_oriented_lo, auc_oriented_hi, rank_in, lands_right)),
      row.names = FALSE)
cat("\nTests: raw p, BH-adjusted p, and the smallest p attainable at that group size:\n")
print(as.data.frame(marker_eff |> filter(!is.na(marks_code)) |>
        select(symbol, n_intended, wilcox_p, wilcox_p_bh, wilcox_p_min_attainable,
               wilcox_p_floor_limited)), row.names = FALSE)
scored <- marker_eff |> filter(!is.na(marks_code))
cat(sprintf("\nMedian |Cohen's d| = %.2f | median oriented AUC = %.2f (n=%d scored markers; MKI67 excluded)\n",
            median(abs(scored$cohens_d), na.rm = TRUE),
            median(scored$auc_oriented, na.rm = TRUE), nrow(scored)))
cat(sprintf("Markers with oriented AUC >= 0.80: %d/%d ; large effect |d|>=0.8: %d/%d\n",
            sum(scored$auc_oriented >= 0.8, na.rm = TRUE), nrow(scored),
            sum(abs(scored$cohens_d) >= 0.8, na.rm = TRUE), nrow(scored)))
cat(sprintf("Significant after BH: %d/%d at 0.05. Floor-limited p-values (already at the exact\n",
            sum(scored$wilcox_p_bh < 0.05, na.rm = TRUE), nrow(scored)))
cat(sprintf("  test's minimum for their group size): %d/%d - these cannot be ranked against each other.\n",
            sum(scored$wilcox_p_floor_limited, na.rm = TRUE), nrow(scored)))
cat(sprintf("Bootstrap CIs (%d stratified replicates) are in cohens_d_lo/hi and auc_oriented_lo/hi.\n", N_BOOT))
cat("  With n=2, percentile bootstrap intervals are poorly calibrated even when all resamples yield finite statistics; the\n")
cat("  finite-resample count does not validate interval coverage or generalisability.\n")

# =============================================================================
# E. FIGURES
# =============================================================================
# ---- F: variance-partition violins (RNA + protein) --------------------------
pretty_term <- c(subtype = "subtype", source_site = "Source site",
                 patient = "Patient family", plex = "TMT plex", Residuals = "Residuals")
# Label the facets with the number of features the model ACTUALLY fitted, not the
# input gene count: vp_lmer drops features whose fit fails, so "22,544 genes" was
# wrong by the number of dropped genes.
rna_lab  <- sprintf("RNA (%s genes)", format(nrow(vp_rna), big.mark = ","))
prot_lab <- sprintf("Protein (%s proteins)", format(nrow(vp_prot), big.mark = ","))
vp_long <- bind_rows(
  as.data.frame(vp_rna)  |> mutate(feature = rownames(as.data.frame(vp_rna))) |>
    pivot_longer(-feature, names_to = "term", values_to = "frac") |> mutate(assay = rna_lab),
  as.data.frame(vp_prot) |> mutate(feature = rownames(as.data.frame(vp_prot))) |>
    pivot_longer(-feature, names_to = "term", values_to = "frac") |>
    mutate(assay = prot_lab)) |>
  mutate(pct = 100 * frac,
         term = recode(term, !!!pretty_term),
         term = factor(term, levels = c("subtype","Source site","Patient family","TMT plex","Residuals")),
         assay = factor(assay, levels = c(rna_lab, prot_lab)))
med_lab <- vp_long |> group_by(assay, term) |> summarise(med = median(pct), .groups = "drop")

pF <- ggplot(vp_long, aes(term, pct, fill = term)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.85) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.35) +
  geom_text(data = med_lab, aes(term, 102, label = sprintf("%.1f%%", med)),
            inherit.aes = FALSE, size = 3, vjust = 0) +
  facet_wrap(~ assay, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = term_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25)) +
  labs(title = "Variance partition: what drives per-feature variance?",
       subtitle = "Per-feature REML decomposition (lme4). Box = IQR; label = genome-wide median. Subtype and site are partly confounded (all HGS share one site).",
       x = NULL, y = "% variance explained") +
  theme_lab() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(pF, "f_variance_partition", 9.5, 5)

# ---- optional: passage-check figure (confound + cross-assay mismatch) --------
pcd <- tibble(PC1 = pca$x[, 1], passage = rna_pass$rna_p, site = rna_pass$site2,
              subtype = rna_pass$subtype, cell_line = ann$cell_line)
pP1 <- ggplot(pcd, aes(passage, PC1, colour = site)) +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_colour_manual(values = c(CHUM = "#C2410C", `BC Cancer` = "#0F172A")) +
  labs(title = "Culture passage is confounded with source site",
       subtitle = sprintf("Passage is %.0f%% explained by site (CHUM high-passage vs BC low-passage); passage alone explains only R²=%.2f of PC1 (n.s.)",
                          100 * pass_site_r2, pass_pc$r2_passage[1]),
       x = "RNA passage number", y = sprintf("PC1 (%.1f%%)", pass_pc$var_pct[1]), colour = "Source") +
  theme_lab()
xa <- xassay |> mutate(cell_line = fct_reorder(cell_line, abs_diff))
pP2 <- ggplot(xa) +
  geom_segment(aes(x = rna_p, xend = wes_p, y = cell_line, yend = cell_line), colour = "grey70", linewidth = 0.8) +
  geom_point(aes(rna_p, cell_line, colour = "RNA"), size = 2.4) +
  geom_point(aes(wes_p, cell_line, colour = "WES"), size = 2.4) +
  scale_colour_manual(values = c(RNA = "#C2410C", WES = "#0F172A"), name = NULL) +
  labs(title = "Within-line RNA vs WES passage mismatch",
       subtitle = "Same line, different passage per assay (e.g. TOV112D RNA p63 vs WES P83)",
       x = "Passage number", y = NULL) +
  theme_lab() + theme(legend.position = c(0.9, 0.2))
save_fig(pP1, "f_passage_check", 6.5, 4.8)
save_fig(pP2, "f_passage_check_crossassay", 6.5, 4.8)

cat("\n--- 17_variance_confounders.R complete ---\n")
cat("Tables: output/{rna_pc_confounder_joint, rna_pc_confounder_permutation, rna_within_cc_site,\n")
cat("        rna_variancepartition, prot_variancepartition, variancepartition_sensitivity,\n")
cat("        rna_passage_check, rna_passage_discordance, rna_marker_effectsizes}.csv\n")
cat("Figures: figs/f_variance_partition.pdf (+ f_passage_check*), reports/assets/*.png\n")

# =============================================================================
# G. Environment record
# =============================================================================
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
write_session_info("17_variance_confounders")

message("\n17_variance_confounders.R complete.")
