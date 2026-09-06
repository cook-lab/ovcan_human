# =============================================================================
# Script: 22_wes_signature_refit.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Mutational-signature refitting and uncertainty, complementing the
#          cosine screen in 16_wes_signatures_msi.R. Exploratory fits do not
#          establish MMR deficiency, MSI or HRD in a model.
# Author:  Cook Lab (analyst: Claude)  |  Date: 2026-07-24
# =============================================================================
#
# WHY THIS SCRIPT EXISTS (peer review M7 / referee minor 12)
#   Script 16 computes cos_sim_matrix() against COSMIC and thresholds a GROUP
#   MAXIMUM. That is a similarity to one reference spectrum, not an exposure
#   estimate: it has no refit, no uncertainty, and no reconstruction check, and
#   SBS1/5/6/15/44 are mutually similar, so TOV21G's cos(MMR-d)=0.877 against a
#   0.651 clock group is a MODEST separation in cosine space. A refit is what can
#   compare attribution under specified reference sets; it does not by itself
#   discriminate clinical MMR status or eliminate tumour-only confounding.
#
#   This script therefore reports, per line:
#     - sparse exposures from fit_to_signatures_strict() (backward selection)
#     - bootstrap distributions from fit_to_signatures_bootstrapped() with 95%
#       intervals AND a per-signature selection frequency (how often across
#       bootstraps the signature survives selection at all) — for a 300-2,400
#       variant spectrum the selection frequency is the honest stability metric
#     - the RECONSTRUCTION cosine (fitted vs observed), which the screen never gave
#     - how many variants entered each spectrum
#     - the exact COSMIC release and package versions
#
# *** SUBSTRATE CAVEATS — read before interpreting any exposure ***
#   1. TUMOUR-ONLY, NO MATCHED NORMAL. The substrate is exome-wide MAF PASS
#      variants with population AF <= 0.001 (the same filter as 07/16): 303-2,417
#      SNVs per line. The retained coding count is not an estimate of residual germline.
#      Unknown germline contamination can distort any fitted component; neither
#      its magnitude nor the direction of signature bias is identified here.
#   2. EXOME SUBSTRATE vs GENOME-DERIVED REFERENCE. COSMIC SBS signatures are
#      defined on genome-wide trinucleotide frequencies; these spectra come from a
#      recovered SeqCap EZ Exome v3 design and derived GRCh38 target bins.
#      The primary historical refit retains genome-reference opportunities; the
#      helper sourced below deposits target-opportunity sensitivity separately.
#      Script44 extends SBS3 filter/dictionary/bootstrap sensitivity. Relative
#      attribution remains sensitive to opportunities, candidate origin and fit design.
#   3. SMALL COUNTS. 300-2,400 variants is at or below the usual floor for stable
#      fitting of 60 reference signatures. This is why (a) selection is strict/
#      backward, (b) a hypothesis-restricted reference set is fitted alongside the
#      full one, and (c) intervals — not point estimates — are reported.
#
# INPUT : output/wes_sbs_context.csv  (96 x 23 SBS matrix written by script 16)
#         output/wes_msi_mmr.csv      (per-line cosine screen + load, for comparison)
# OUTPUT: output/wes_signature_refit_exposures.csv
#         output/wes_signature_refit_summary.csv
#         figs/f_wes_signature_refit.pdf
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})
# MutationalPatterns pulls in S4Vectors, which exports rename/first/second and masks
# the dplyr verbs. Pin the dplyr ones explicitly (same idiom as scripts 10/11/18).
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count  <- dplyr::count;  slice  <- dplyr::slice
set.seed(SEED)

SBS_CTX  <- file.path(OUT, "wes_sbs_context.csv")
stopifnot("output/wes_sbs_context.csv missing — run 16_wes_signatures_msi.R first" =
            file.exists(SBS_CTX))

# Reference release must match 16 exactly, or the refit and the screen are not
# comparable. COSMIC_v3.2 is the newest COSMIC bundled with MutationalPatterns.
COSMIC_SOURCE <- "COSMIC_v3.2"
COSMIC_GENOME <- "GRCh38"
N_BOOT        <- 200L      # bootstrap replicates per line
MAX_DELTA     <- 0.004     # fit_to_signatures_strict default: cosine cost of dropping a sig
MMR_D_SIGS <- c("SBS6","SBS14","SBS15","SBS20","SBS21","SBS26","SBS44")
POLE_SIGS  <- c("SBS10a","SBS10b","SBS10c","SBS10d","SBS28")
CLOCK_SIGS <- c("SBS1","SBS5","SBS40")
# Hypothesis-restricted reference retained as a sensitivity prior. The omitted
# signatures are not established to be irrelevant to these models. Compare with
# the full 60-signature reference; a restricted-only exposure is not confirmation
# of a pathway defect. The set excludes alternatives including platinum signatures.
RESTRICTED_SIGS <- unique(c(CLOCK_SIGS, MMR_D_SIGS, POLE_SIGS,
                            "SBS2","SBS13",   # APOBEC
                            "SBS3",           # HR deficiency
                            "SBS8","SBS17a","SBS17b","SBS18"))  # common in HGSC series

theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace% theme(
    text = element_text(colour = "black"),
    plot.title = element_text(size = rel(1.1), hjust = 0, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = rel(0.82), hjust = 0, colour = "grey30", margin = margin(b = 6)),
    axis.title = element_text(size = rel(0.95)), axis.text = element_text(size = rel(0.8), colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    legend.title = element_text(size = rel(0.8)), legend.text = element_text(size = rel(0.75)),
    legend.key = element_blank(), panel.background = element_blank(),
    panel.grid = element_blank(), strip.background = element_blank(),
    strip.text = element_text(size = rel(0.85), face = "bold"), plot.margin = margin(6, 8, 6, 6))
}

# =============================================================================
# 1. DEPENDENCIES — report, never silently skip
# =============================================================================
need <- c("MutationalPatterns")
miss <- setdiff(need, rownames(installed.packages()))
if (length(miss)) {
  stop("Signature refit requires: ", paste(miss, collapse = ", "),
       "\n  Install: BiocManager::install(c(", paste0('"', miss, '"', collapse = ", "), "))",
       "\n  NOT skipping silently: the manuscript's signature claim depends on this refit.")
}
suppressPackageStartupMessages(library(MutationalPatterns))
MP_VER <- as.character(packageVersion("MutationalPatterns"))
message("MutationalPatterns ", MP_VER, " | COSMIC ", COSMIC_SOURCE, " / ", COSMIC_GENOME,
        " | ", N_BOOT, " bootstraps | max_delta ", MAX_DELTA)

# =============================================================================
# 2. SUBSTRATE — the 96-context matrix written by 16 (identical variant set)
# =============================================================================
ctx <- readr::read_csv(SBS_CTX, show_col_types = FALSE)
mut_mat <- as.matrix(ctx[, -1]); rownames(mut_mat) <- ctx$context      # 96 x 23
n_used  <- colSums(mut_mat)
message(sprintf("SBS-96 matrix: %d contexts x %d lines | variants per spectrum: median %.0f (range %d-%d)",
                nrow(mut_mat), ncol(mut_mat), median(n_used), min(n_used), max(n_used)))
stopifnot("SBS matrix is not 96 contexts" = nrow(mut_mat) == 96L)

screen <- readr::read_csv(file.path(OUT, "wes_msi_mmr.csv"), show_col_types = FALSE) %>%
  select(cell_line, subtype, n_coding, is_hypermutator, msi_mmr_call,
         any_of(c("cos_mmr_d","cos_clock","cos_pole","cos_margin_mmrd_vs_clock",
                  "sbs_call","best_cosmic")))

cosmic <- get_known_signatures(muttype = "snv", source = COSMIC_SOURCE, genome = COSMIC_GENOME)
message(sprintf("COSMIC reference: %d signatures (%s)", ncol(cosmic),
                paste(head(colnames(cosmic), 6), collapse = ", ")))
restricted <- cosmic[, intersect(RESTRICTED_SIGS, colnames(cosmic)), drop = FALSE]
message(sprintf("Restricted reference: %d signatures", ncol(restricted)))

# =============================================================================
# 3. STRICT REFIT (sparse, backward selection) — full and restricted references
# =============================================================================
refit_one <- function(sigs, tag) {
  fit <- fit_to_signatures_strict(mut_mat, sigs, max_delta = MAX_DELTA, method = "backwards")
  contrib <- fit$fit_res$contribution           # signatures x lines (absolute variant counts)
  recon   <- fit$fit_res$reconstructed          # 96 x lines
  cosines <- vapply(colnames(mut_mat), function(s) cos_sim(mut_mat[, s], recon[, s]), numeric(1))
  list(tag = tag, contribution = contrib, reconstruction_cosine = cosines,
       n_selected = colSums(contrib > 0))
}
full_fit <- refit_one(cosmic,     "full COSMIC v3.2 (60 sigs)")
rest_fit <- refit_one(restricted, sprintf("restricted (%d sigs)", ncol(restricted)))

cat("\n=== Reconstruction cosine (fitted vs observed spectrum) — the check the screen never made ===\n")
print(as.data.frame(tibble(
  cell_line   = colnames(mut_mat),
  n_snv_used  = as.integer(n_used),
  recon_cos_full = round(full_fit$reconstruction_cosine, 3),
  n_sigs_full    = as.integer(full_fit$n_selected),
  recon_cos_restricted = round(rest_fit$reconstruction_cosine, 3),
  n_sigs_restricted    = as.integer(rest_fit$n_selected)) %>%
  arrange(desc(n_snv_used))), row.names = FALSE)

# =============================================================================
# 4. BOOTSTRAP — intervals + selection frequency (the honest stability metric)
# -----------------------------------------------------------------------------
# This is the expensive step (200 replicates x 23 lines x 2 reference sets, ~40 min
# with backward selection over 60 signatures), so the SUMMARY is cached to
# output/wes_signature_refit_bootstrap.csv. The cache is reused only when the
# recorded parameters, reference-set sizes and line set all match; otherwise it is
# recomputed. Delete the file (or set OVCAN_REFIT_FORCE=1) to force a fresh run.
# =============================================================================
BOOT_CACHE <- file.path(OUT, "wes_signature_refit_bootstrap.csv")
FORCE_BOOT <- Sys.getenv("OVCAN_REFIT_FORCE", unset = "0") == "1"

# Cache identity covers the actual counts/reference spectra, software and RNG
# policy. Parameter counts alone cannot detect a changed spectrum or reference.
BOOT_CACHE_SHA <- digest::digest(list(mut_mat = mut_mat, cosmic = cosmic,
  restricted = restricted, version = MP_VER, seed = SEED, n_boot = N_BOOT,
  max_delta = MAX_DELTA, method = "strict_backwards", seed_policy = "sample-and-reference-v1"),
  algo = "sha256")
BOOT_WORKERS <- max(1L, as.integer(Sys.getenv("OVCAN_BOOT_WORKERS", "4")))
boot_summary <- function(sigs, tag) {
  one <- function(sample_name) {
    sample_seed <- (strtoi(substr(digest::digest(paste(tag, sample_name), serialize = FALSE), 1, 7), 16L) + SEED)
    set.seed(sample_seed)
    bt <- fit_to_signatures_bootstrapped(mut_mat[, sample_name, drop = FALSE], sigs,
      n_boots = N_BOOT, max_delta = MAX_DELTA, method = "strict_backwards", verbose = FALSE)
    # MutationalPatterns drops never-selected signatures; restore explicit zeros
    # before summarising so absence cannot masquerade as missing cache coverage.
    aligned <- matrix(0, nrow(bt), ncol(sigs), dimnames = list(rownames(bt), colnames(sigs)))
    aligned[, colnames(bt)] <- bt
    as_tibble(as.data.frame(aligned), rownames = "row") %>%
      pivot_longer(-row, names_to = "signature", values_to = "contribution") %>%
      group_by(signature) %>%
      summarise(boot_mean = mean(contribution), boot_median = median(contribution),
        boot_lo95 = quantile(contribution, 0.025), boot_hi95 = quantile(contribution, 0.975),
        boot_selected_frac = mean(contribution > 0), n_boots = dplyr::n(), .groups = "drop") %>%
      mutate(cell_line = sample_name, reference_set = tag, bootstrap_seed = sample_seed)
  }
  res <- parallel::mclapply(colnames(mut_mat), one, mc.cores = BOOT_WORKERS,
                            mc.set.seed = FALSE, mc.preschedule = FALSE)
  stopifnot("A signature bootstrap worker failed" = !any(vapply(res, inherits, FALSE, "try-error")))
  bind_rows(res)
}

cache_ok <- FALSE
if (file.exists(BOOT_CACHE) && !FORCE_BOOT) {
  cached <- readr::read_csv(BOOT_CACHE, show_col_types = FALSE)
  cache_ok <- all(c("cell_line","signature","reference_set","boot_selected_frac",
                    "n_boots","cache_max_delta","cache_n_ref_full","cache_n_ref_restricted","cache_sha256") %in%
                  names(cached)) &&
    setequal(unique(cached$cell_line), colnames(mut_mat)) &&
    all(!is.na(cached$cache_sha256) & cached$cache_sha256 == BOOT_CACHE_SHA) &&
    all(cached$n_boots == N_BOOT, na.rm = TRUE) &&
    all(cached$cache_max_delta == MAX_DELTA, na.rm = TRUE) &&
    all(cached$cache_n_ref_full == ncol(cosmic), na.rm = TRUE) &&
    all(cached$cache_n_ref_restricted == ncol(restricted), na.rm = TRUE) &&
    # COVERAGE [integration revision]. The checks above are all on cache KEYS; none
    # of them notices a cache that is keyed correctly but INCOMPLETE. That is not
    # hypothetical: the pre-integration cache on disk held only the ever-selected
    # (line, signature) pairs (916 full + 334 restricted = 1,250) rather than the
    # full grid for that historical 22-model run (22x60 + 22x22 = 1,804), and it
    # satisfied every key check. Reusing such a cache is silently destructive,
    # because the join at "5. TIDY EXPOSURE TABLE" treats an absent row as
    # boot_selected_frac = 0 via replace_na() and then DROPS it from
    # wes_signature_refit_exposures.csv. Require the complete grid instead.
    sum(cached$reference_set == "full")       == ncol(mut_mat) * ncol(cosmic) &&
    sum(cached$reference_set == "restricted") == ncol(mut_mat) * ncol(restricted)
  message(if (cache_ok) "Reusing cached bootstrap summary: " else
          "Cached bootstrap does not match the current parameters, recomputing: ",
          BOOT_CACHE)
}
if (cache_ok) {
  boot_full <- cached %>% filter(reference_set == "full")
  boot_rest <- cached %>% filter(reference_set == "restricted")
} else {
  message("Bootstrapping the strict refit (", N_BOOT, " replicates x ", ncol(mut_mat),
          " lines x 2 reference sets) — this takes tens of minutes...")
  boot_full <- boot_summary(cosmic,     "full")
  boot_rest <- boot_summary(restricted, "restricted")
  readr::write_csv(bind_rows(boot_full, boot_rest) %>%
                     mutate(cache_max_delta = MAX_DELTA, cache_sha256 = BOOT_CACHE_SHA,
                            cache_n_ref_full = ncol(cosmic),
                            cache_n_ref_restricted = ncol(restricted),
                            cosmic_source = COSMIC_SOURCE,
                            mutationalpatterns_version = MP_VER),
                   BOOT_CACHE)
  message("Wrote bootstrap cache -> ", BOOT_CACHE)
}

# =============================================================================
# 5. TIDY EXPOSURE TABLE
# =============================================================================
tidy_contrib <- function(fit, tag) {
  # NB use if_else on a per-row `total` column, NOT ifelse(sum(x) > 0, ...): with a
  # length-1 condition ifelse() returns a length-1 result, which silently recycles the
  # FIRST signature's fraction across every row of the group.
  as_tibble(as.data.frame(fit$contribution), rownames = "signature") %>%
    pivot_longer(-signature, names_to = "cell_line", values_to = "contribution") %>%
    group_by(cell_line) %>%
    mutate(total_fitted = sum(contribution),
           rel_contribution = if_else(total_fitted > 0, contribution / total_fitted,
                                      NA_real_)) %>%
    ungroup() %>% mutate(reference_set = tag)
}
sig_group <- function(s) case_when(s %in% MMR_D_SIGS ~ "MMR-d / MSI",
                                  s %in% POLE_SIGS  ~ "POLE proofreading",
                                  s %in% CLOCK_SIGS ~ "clock-like (age)",
                                  TRUE               ~ "other")

boot_all <- bind_rows(boot_full, boot_rest) %>%
  select(cell_line, signature, reference_set, boot_mean, boot_median,
         boot_lo95, boot_hi95, boot_selected_frac, n_boots)
# The filter at the end of this pipeline uses replace_na(boot_selected_frac, 0), so a
# gap in boot_all silently removes a signature from the exposure table instead of
# erroring. Assert the bootstrap covers the whole grid before that can happen.
# [integration revision]
stopifnot("bootstrap summary does not cover every (line, signature, reference_set) — cache is incomplete" =
            nrow(boot_all) == ncol(mut_mat) * (ncol(cosmic) + ncol(restricted)),
          "bootstrap summary has duplicate (line, signature, reference_set) rows" =
            !any(duplicated(boot_all[, c("cell_line", "signature", "reference_set")])))
exposures <- bind_rows(tidy_contrib(full_fit, "full"),
                       tidy_contrib(rest_fit, "restricted")) %>%
  left_join(boot_all, by = c("cell_line", "signature", "reference_set")) %>%
  mutate(signature_group = sig_group(signature),
         n_snv_used = as.integer(n_used[cell_line]),
         cosmic_source = COSMIC_SOURCE, cosmic_genome = COSMIC_GENOME,
         mutationalpatterns_version = MP_VER,
         fit_method = "fit_to_signatures_strict(method='backwards')",
         max_delta = MAX_DELTA, n_boots_requested = N_BOOT, bootstrap_input_sha256 = BOOT_CACHE_SHA) %>%
  left_join(screen %>% select(cell_line, subtype, is_hypermutator), by = "cell_line") %>%
  # keep every signature the fit or the bootstrap ever touched; drop the always-zero
  # ones so the table is readable rather than 23 x 60 x 2 mostly-zero rows
  filter(contribution > 0 | replace_na(boot_selected_frac, 0) > 0) %>%
  arrange(reference_set, cell_line, desc(rel_contribution))
readr::write_csv(exposures, file.path(OUT, "wes_signature_refit_exposures.csv"))
message("Wrote output/wes_signature_refit_exposures.csv (", nrow(exposures), " rows)")

# =============================================================================
# 6. PER-LINE SUMMARY — refit vs screen, side by side
# =============================================================================
grp_rel <- function(df, g) df %>% filter(signature_group == g) %>%
  group_by(cell_line) %>% summarise(v = sum(rel_contribution, na.rm = TRUE), .groups = "drop")
ex_full <- exposures %>% filter(reference_set == "full")
top_sig <- ex_full %>% group_by(cell_line) %>%
  slice_max(rel_contribution, n = 1, with_ties = FALSE) %>%
  transmute(cell_line, top_signature = signature,
            top_rel = round(rel_contribution, 3),
            top_boot_selected_frac = round(boot_selected_frac, 3))
mmrd_rel <- grp_rel(ex_full, "MMR-d / MSI") %>% rename(rel_mmr_d = v)
clock_rel <- grp_rel(ex_full, "clock-like (age)") %>% rename(rel_clock = v)
pole_rel <- grp_rel(ex_full, "POLE proofreading") %>% rename(rel_pole = v)
# strongest MMR-d signature per line, with its bootstrap interval — the number that
# replaces "SBS6 cosine 0.88" as the evidence for the MMR-deficiency claim
mmrd_best <- ex_full %>% filter(signature_group == "MMR-d / MSI") %>%
  group_by(cell_line) %>% slice_max(rel_contribution, n = 1, with_ties = FALSE) %>%
  transmute(cell_line, mmr_d_top_signature = signature,
            mmr_d_top_rel = round(rel_contribution, 3),
            mmr_d_top_boot_lo95 = round(boot_lo95, 1),
            mmr_d_top_boot_hi95 = round(boot_hi95, 1),
            mmr_d_top_boot_selected_frac = round(boot_selected_frac, 3))

rest_mmrd <- grp_rel(exposures %>% filter(reference_set == "restricted"), "MMR-d / MSI") %>%
  rename(rel_mmr_d_restricted = v)
refit_summary <- tibble(
    cell_line = colnames(mut_mat),
    n_snv_used = as.integer(n_used),
    n_sigs_selected_full = as.integer(full_fit$n_selected),
    recon_cosine_full = round(full_fit$reconstruction_cosine, 4),
    n_sigs_selected_restricted = as.integer(rest_fit$n_selected),
    recon_cosine_restricted = round(rest_fit$reconstruction_cosine, 4)) %>%
  left_join(screen, by = "cell_line") %>%
  left_join(rest_mmrd, by = "cell_line") %>%
  left_join(top_sig,  by = "cell_line") %>%
  left_join(mmrd_best, by = "cell_line") %>%
  left_join(mmrd_rel,  by = "cell_line") %>%
  left_join(clock_rel, by = "cell_line") %>%
  left_join(pole_rel,  by = "cell_line") %>%
  mutate(across(c(rel_mmr_d, rel_clock, rel_pole), ~round(replace_na(., 0), 3)),
         cosmic_source = COSMIC_SOURCE, cosmic_genome = COSMIC_GENOME,
         mutationalpatterns_version = MP_VER, n_boots = N_BOOT, max_delta = MAX_DELTA,
         substrate = "exome-wide Mutect2 PASS & popAF<=0.001 SNVs (tumour-only, no matched normal)",
         exome_renormalisation = "not applied to primary screen/refit; archived target-bin opportunity sensitivity deposited separately") %>%
  arrange(desc(rel_mmr_d))
readr::write_csv(refit_summary, file.path(OUT, "wes_signature_refit_summary.csv"))
message("Wrote output/wes_signature_refit_summary.csv (", nrow(refit_summary), " lines)")

cat("\n=== REFIT vs COSINE SCREEN — does the MMR-d component survive a fit? ===\n")
print(as.data.frame(refit_summary %>%
  select(cell_line, subtype, n_snv_used, cos_mmr_d, cos_margin_mmrd_vs_clock,
         rel_mmr_d, mmr_d_top_signature, mmr_d_top_boot_selected_frac,
         rel_clock, recon_cosine_full, sbs_call) %>%
  mutate(across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

hyper <- refit_summary %>% filter(replace_na(is_hypermutator, FALSE))
if (nrow(hyper)) {
  others <- refit_summary %>% filter(!replace_na(is_hypermutator, FALSE))
  cat(sprintf("\n=> %s: MMR-d relative exposure %.3f (top %s, selected in %.1f%% of %d bootstraps);\n   the other %d lines span %.3f-%.3f. Reconstruction cosine %.3f vs %.3f-%.3f.\n",
              hyper$cell_line[1], hyper$rel_mmr_d[1], hyper$mmr_d_top_signature[1],
              100 * hyper$mmr_d_top_boot_selected_frac[1], N_BOOT,
              nrow(others), min(others$rel_mmr_d), max(others$rel_mmr_d),
              hyper$recon_cosine_full[1],
              min(others$recon_cosine_full), max(others$recon_cosine_full)))
  cat("   Interpretation: a refit estimates contributions conditional on a chosen reference set,\n")
  cat("   which a cosine to one reference cannot. Read the bootstrap interval and the\n")
  cat("   selection frequency, not the point exposure — n_snv_used is small.\n")
}

# =============================================================================
# 7. FIGURE — exposures with bootstrap intervals (diagnostic; for a figure agent)
# =============================================================================
grp_cols <- c("MMR-d / MSI" = COOK_RUST, "POLE proofreading" = COOK_NAVY,
              "clock-like (age)" = "grey55", "other" = "grey80")
lvl <- refit_summary$cell_line

pA <- ex_full %>%
  mutate(cell_line = factor(cell_line, levels = rev(lvl))) %>%
  ggplot(aes(rel_contribution, cell_line, fill = signature_group)) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)), labels = scales::percent) +
  labs(title = "A  Strict refit exposures (full COSMIC v3.2)",
       subtitle = sprintf("fit_to_signatures_strict, backward selection, max_delta=%.3f; ordered by MMR-d exposure",
                          MAX_DELTA),
       x = "relative contribution", y = NULL) +
  theme_lab() + theme(legend.position = "top")

top_line <- lvl[1]
pB <- ex_full %>% filter(cell_line == top_line) %>%
  mutate(signature = reorder(signature, boot_median)) %>%
  ggplot(aes(boot_median, signature, colour = signature_group)) +
  geom_errorbarh(aes(xmin = boot_lo95, xmax = boot_hi95), height = 0, linewidth = 0.6) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%.0f%% boots", 100 * boot_selected_frac)),
            hjust = -0.2, size = 2.2, colour = "grey25") +
  scale_colour_manual(values = grp_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.30))) +
  labs(title = sprintf("B  %s exposures with bootstrap 95%% intervals", top_line),
       subtitle = sprintf("%d bootstraps of the strict refit; %s SNVs entered the spectrum",
                          N_BOOT, scales::comma(n_used[top_line])),
       x = "fitted variant count", y = NULL) +
  theme_lab() + theme(legend.position = "none")

pC <- refit_summary %>%
  mutate(cell_line = factor(cell_line, levels = rev(lvl))) %>%
  ggplot(aes(recon_cosine_full, cell_line)) +
  geom_col(width = 0.74, fill = COOK_NAVY) +
  geom_vline(xintercept = 0.90, linetype = 2, colour = "grey45", linewidth = 0.4) +
  coord_cartesian(xlim = c(0.5, 1)) +
  labs(title = "C  Reconstruction cosine",
       subtitle = "fitted vs observed 96-context spectrum (dashed 0.90)",
       x = "cosine(observed, reconstructed)", y = NULL) +
  theme_lab()

fig <- (pA | (pB / pC)) + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    title = "Exploratory mutational-signature refit with uncertainty",
    subtitle = paste0("Tumour-only exome-wide MAF PASS & popAF<=0.001 SNVs (", min(n_used), "-",
                      max(n_used), " per line); COSMIC ", COSMIC_SOURCE, "/", COSMIC_GENOME,
                      "; MutationalPatterns ", MP_VER, ". Primary fit uses genome opportunities; ",
                      "recovered-target sensitivity is reported separately. Residual germline can distort fitted components."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 8.2, colour = "grey30")))
ggsave(file.path(FIGS, "f_wes_signature_refit.pdf"), fig, width = 13.5, height = 8.6)
message("Wrote figs/f_wes_signature_refit.pdf")

message("\n22_wes_signature_refit.R complete. Outputs: ",
        "output/wes_signature_refit_exposures.csv, output/wes_signature_refit_summary.csv, ",
        "figs/f_wes_signature_refit.pdf")

# The archived derived target bins permit a bounded opportunity sensitivity.
source("scripts/lib/wes_target_signature_sensitivity.R")
if (!is.null(warnings())) writeLines(capture.output(warnings()),
  file.path(OUT, "wes_signature_run_warnings.txt"))
