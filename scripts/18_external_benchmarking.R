# =============================================================================
# Script: 18_external_benchmarking.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Revision workstream — EXTERNAL identity / benchmarking against public
#          references, addressing peer-review §4 (essential #5) and §3.2:
#
#          TASK 1  CCLE / DepMap expression concordance (the strongest external
#                  identity check). (a) enumerate which of our 42 lines exist in
#                  DepMap (via Model.csv + Cellosaurus RRIDs); (b) for the
#                  overlap, correlate THIS project's RNA against DepMap 24Q4
#                  expression and show a high on-diagonal self-correlation that
#                  stands out from all other DepMap ovarian lines; plus a driver-
#                  mutation cross-check (TP53 / SMARCA4 / clear-cell drivers,
#                  TOV21G hypermutation) against DepMap somatic-mutation matrices.
#          TASK 2  COMPUTE the ConsensusOV subtype calls for the HGSC RNA lines with
#                  the consensusOV package (version recorded), report the per-class
#                  random-forest PROBABILITIES and the top-vs-second MARGIN, and
#                  reconcile them against the labels inherited from samples.csv
#                  `notes`. The caveat stands and the margins strengthen it:
#                  TCGA/ConsensusOV subtypes are substantially TME/stroma-driven and
#                  thus of questionable validity on pure tumor-cell cultures.
#                  [review revision — the calls were previously STRING-EXTRACTED from
#                  free-text notes with no classifier version recorded, while the
#                  manuscript claimed they "were generated ... (script 10)".]
#          TASK 3  Cellosaurus accessions + STR / mycoplasma documentation status
#                  for all 42 lines (a Scientific Data authentication requirement;
#                  current-stock STR documentation is unconfirmed). Records the Cellosaurus
#                  accession, whether the originator has a documented STR profile,
#                  and any "problematic / misidentified" Cellosaurus flag.
#
# Author:  Cook Lab (analyst: Claude — external-benchmark workstream)
# Date:    2026-07-23
# Phase:   Revision
# -----------------------------------------------------------------------------
# EXTERNAL DATA (downloaded to output/external/ by this workstream; see delta):
#   DepMap Public 24Q4 (Figshare article 27993248):
#     Model.csv                                  (cell-model metadata + RRID)
#     OmicsExpressionProteinCodingGenesTPMLogp1  (log2(TPM+1); rows=models)
#     OmicsSomaticMutationsMatrixHotspot.csv     (hotspot count per gene)
#     OmicsSomaticMutationsMatrixDamaging.csv    (damaging count per gene)
#   Cellosaurus API (api.cellosaurus.org) JSON cached in output/external/cellosaurus/
# IDENTITY RULE: DepMap IDs (ACH-*) and Cellosaurus accessions (CVCL_*) are never
#   invented — every one is verified against Model.csv / the Cellosaurus record.
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(tidyverse); library(matrixStats); library(data.table)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(jsonlite); library(patchwork)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
set.seed(SEED)

EXT  <- file.path(OUT, "external")
CELLO <- file.path(EXT, "cellosaurus")
DEPMAP_RELEASE <- "DepMap Public 24Q4"

# --- lab viz helpers (visualization skill; identical to scripts 14/17) --------
theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(text = element_text(colour = "black"),
          plot.title = element_text(size = 13, margin = margin(b = 6), hjust = 0, face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey30", margin = margin(b = 6), hjust = 0),
          axis.title = element_text(size = 11), axis.text = element_text(size = 9, colour = "black"),
          legend.title = element_text(size = 9), legend.text = element_text(size = 8),
          strip.text = element_text(size = 10, face = "bold"),
          axis.line = element_line(colour = "black", linewidth = 0.4),
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

# =============================================================================
# LOAD our data + external metadata
# =============================================================================
samp <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE)
ours <- samp %>% filter(provenance == "generated", analysis_include == "Y")
stopifnot(nrow(ours) == 42)

model <- readr::read_csv(file.path(EXT, "Model.csv"), show_col_types = FALSE)
ov_models <- readr::read_csv(file.path(EXT, "depmap_ovarian_models.csv"), show_col_types = FALSE)

norm_nm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

# =============================================================================
# PART 1a — DepMap OVERLAP (recomputed, not hardcoded)
# =============================================================================
model <- model %>% mutate(stripped_dm = norm_nm(StrippedCellLineName))
ours  <- ours  %>% mutate(stripped = norm_nm(cell_line))

overlap <- ours %>%
  inner_join(model %>% select(ModelID, StrippedCellLineName, RRID,
                              OncotreeLineage, OncotreePrimaryDisease, OncotreeSubtype,
                              stripped_dm),
             by = c("stripped" = "stripped_dm")) %>%
  select(cell_line, subtype, source_site, rna_seq, wes_mut,
         depmap_ACH = ModelID, depmap_name = StrippedCellLineName, depmap_RRID = RRID,
         depmap_oncotree_lineage = OncotreeLineage,
         depmap_oncotree_subtype = OncotreeSubtype)

message("\n== DepMap overlap: ", nrow(overlap), " of 42 lines ==")
print(as.data.frame(overlap %>% select(cell_line, subtype, depmap_ACH, depmap_RRID, depmap_oncotree_subtype)))
stopifnot(all(c("OV90","TOV21G","TOV112D","BIN67","COV434") %in% overlap$cell_line))

# =============================================================================
# PART 1b — EXPRESSION CONCORDANCE
#   ours (linear TPM) -> collapse to Entrez, log2(TPM+1)   [matches DepMap units]
#   DepMap (already log2(TPM+1)) -> collapse to Entrez
#   Spearman (primary; rank-based, robust to cross-platform scale) + Pearson,
#   over the top-2000 most variable shared genes.
#
# *** TWO ASYMMETRIES THAT INFLATE THE ABSOLUTE CORRELATIONS *** [review revision]
#   (1) The two sides are collapsed to Entrez DIFFERENTLY, and this cannot be fixed
#       without DepMap's pre-collapse matrix: OUR side sums LINEAR TPM across the
#       ENSGs mapping to one Entrez ID and then logs it, whereas DepMap ships
#       already-logged values, so its duplicate Entrez columns are AVERAGED IN LOG
#       SPACE (collapse_mean below). sum-then-log != mean-of-logs, so a systematic
#       offset exists for every multi-ENSG gene.
#   (2) HVGs are selected on the COMBINED 31+67 matrix, so genes that differ most
#       between the two PLATFORMS are preferentially selected.
#   Both inflate the absolute Spearman (0.74-0.88) for SELF and NON-SELF pairs
#   ALIKE, which is exactly why SPECIFICITY — the self-vs-best-other margin and the
#   rank of the self match among all 67 DepMap ovarian lines — IS THE SIGNAL, and the
#   absolute coefficient is not. The margin table written below is the reportable
#   statistic; a sensitivity correlation over ALL shared genes (S_all, no HVG
#   selection) is carried alongside so the HVG effect is visible.
# =============================================================================
collapse_mean <- function(mat) {                 # average rows sharing a rowname
  if (!any(duplicated(rownames(mat)))) return(mat)
  ag  <- rowsum(mat, group = rownames(mat))
  cnt <- as.integer(table(rownames(mat))[rownames(ag)])
  sweep(ag, 1, cnt, "/")
}

## our TPM -> Entrez log2(TPM+1)
tpm <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
our_ensg <- as.matrix(tpm[, -1]); rownames(our_ensg) <- tpm$gene_id     # linear TPM, ENSG x 31
map <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(our_ensg),
                             keytype = "ENSEMBL", columns = "ENTREZID") %>%
  filter(!is.na(ENTREZID))
map <- map[!duplicated(map$ENSEMBL), ]           # 1 Entrez per ENSG (drop 1:many ambiguity)
ent <- setNames(map$ENTREZID, map$ENSEMBL)[rownames(our_ensg)]
keep <- !is.na(ent)
our_ent <- rowsum(our_ensg[keep, ], group = ent[keep])   # sum linear TPM per Entrez
our_log <- log2(our_ent + 1)                             # Entrez x 31

## DepMap ovarian expression -> Entrez (already log2(TPM+1))
dm_raw <- data.table::fread(file.path(EXT, "depmap_expr_ovarian.csv"),
                            data.table = FALSE, check.names = FALSE)
colnames(dm_raw)[1] <- "ModelID"
gcols  <- colnames(dm_raw)[-1]
gentrez <- sub(".*\\(([0-9]+)\\)$", "\\1", gcols)
dm_mat <- t(as.matrix(dm_raw[, -1])); rownames(dm_mat) <- gentrez; colnames(dm_mat) <- dm_raw$ModelID
dm_log <- collapse_mean(dm_mat)                          # Entrez x 67

## common genes + HVG
common <- intersect(rownames(our_log), rownames(dm_log))
comb   <- cbind(our_log[common, ], dm_log[common, ])
ord    <- order(matrixStats::rowVars(comb), decreasing = TRUE)
n_hvg  <- 2000L
hvg    <- common[ord][seq_len(n_hvg)]
message(sprintf("  shared Entrez genes: %d ; HVG used: %d", length(common), n_hvg))

## correlation matrices: our 31 (rows) x DepMap 67 (cols)
S_hvg <- cor(our_log[hvg, ], dm_log[hvg, ], method = "spearman")   # 31 x 67
P_hvg <- cor(our_log[hvg, ], dm_log[hvg, ], method = "pearson")
S_all <- cor(our_log[common, ], dm_log[common, ], method = "spearman")  # sensitivity

ach2name <- setNames(ov_models$StrippedCellLineName, ov_models$ModelID)
ach2sub  <- setNames(ov_models$OncotreeSubtype,     ov_models$ModelID)

## per-overlap-line identity statistics (self vs all 67 DepMap ovarian)
concordance <- overlap %>%
  rowwise() %>%
  mutate(
    self_spearman_hvg  = S_hvg[cell_line, depmap_ACH],
    self_pearson_hvg   = P_hvg[cell_line, depmap_ACH],
    self_spearman_all  = S_all[cell_line, depmap_ACH],
    # rank of the correct DepMap line among all 67 DepMap ovarian lines (1 = best)
    self_rank_of67     = rank(-S_hvg[cell_line, ])[depmap_ACH],
    n_depmap_ovarian   = ncol(S_hvg),
    # 2nd-best DepMap ovarian match for our line + margin
    runnerup_ACH       = names(sort(S_hvg[cell_line, ], decreasing = TRUE))[2],
    runnerup_spearman  = sort(S_hvg[cell_line, ], decreasing = TRUE)[2],
    margin_self_vs_next= self_spearman_hvg - runnerup_spearman,
    # reciprocal check: among our 31 lines, is this our-line the best match for the DepMap line?
    depmap_best_our    = rownames(S_hvg)[which.max(S_hvg[, depmap_ACH])],
    reciprocal_best    = depmap_best_our == cell_line
  ) %>% ungroup() %>%
  mutate(runnerup_name = ach2name[runnerup_ACH]) %>%
  select(cell_line, subtype, depmap_ACH, depmap_name, depmap_RRID,
         depmap_oncotree_subtype, self_spearman_hvg, self_pearson_hvg,
         self_spearman_all, self_rank_of67, n_depmap_ovarian,
         runnerup_name, runnerup_spearman, margin_self_vs_next,
         depmap_best_our, reciprocal_best)

message("\n== Expression identity (self vs 67 DepMap ovarian) ==")
print(as.data.frame(concordance %>% select(cell_line, depmap_name, self_spearman_hvg,
                                            self_rank_of67, runnerup_name, runnerup_spearman,
                                            margin_self_vs_next, reciprocal_best)))

## long-form 31 x 5 table for the CSV + heatmap (our lines x the 5 overlap DepMap lines)
long <- expand_grid(cell_line = rownames(S_hvg), depmap_ACH = overlap$depmap_ACH) %>%
  mutate(spearman_hvg  = S_hvg[cbind(cell_line, depmap_ACH)],
         pearson_hvg   = P_hvg[cbind(cell_line, depmap_ACH)],
         depmap_name   = ach2name[depmap_ACH],
         our_subtype   = ours$subtype[match(cell_line, ours$cell_line)],
         is_selfpair   = norm_nm(cell_line) == norm_nm(depmap_name)) %>%
  left_join(overlap %>% select(depmap_ACH, depmap_RRID, depmap_oncotree_subtype), by = "depmap_ACH") %>%
  left_join(concordance %>% select(cell_line, self_rank_of67, margin_self_vs_next), by = "cell_line") %>%
  mutate(self_rank_of67 = ifelse(is_selfpair, self_rank_of67, NA_integer_),
         margin_self_vs_next = ifelse(is_selfpair, margin_self_vs_next, NA_real_)) %>%
  select(our_line = cell_line, our_subtype, depmap_name, depmap_ACH, depmap_RRID,
         depmap_oncotree_subtype, spearman_hvg, pearson_hvg, is_selfpair,
         self_rank_of67, margin_self_vs_next)

readr::write_csv(long, file.path(OUT, "external_ccle_concordance.csv"))
message("  wrote output/external_ccle_concordance.csv (", nrow(long), " rows)")

# --- SELF-MATCH MARGIN AGAINST ALL 67 DepMap OVARIAN LINES [review revision] ---
# The manuscript claims "rank 1 of 67, reciprocal-best in both directions", but the
# only figure is a 5x5 submatrix, which cannot show a rank of 1 in 67 and cannot show
# the MARGIN — and with the collapse/HVG asymmetries above, the margin is the claim
# that carries. Two tables are emitted so a figure agent can draw it directly:
#   external_selfmatch_margin.csv    one row per overlap line: self, best non-self,
#                                    margin, rank of 67, reciprocal-best, z of self
#                                    against that line's own non-self distribution
#   external_depmap_spearman_all.csv the full 31 x 67 long matrix behind it
selfmatch <- overlap %>% rowwise() %>%
  mutate(
    self_spearman       = S_hvg[cell_line, depmap_ACH],
    self_rank_of67      = as.integer(rank(-S_hvg[cell_line, ])[depmap_ACH]),
    n_depmap_ovarian    = ncol(S_hvg),
    best_nonself_ACH    = names(which.max(S_hvg[cell_line, setdiff(colnames(S_hvg), depmap_ACH)])),
    best_nonself_spearman = max(S_hvg[cell_line, setdiff(colnames(S_hvg), depmap_ACH)]),
    margin_self_minus_best_nonself = self_spearman - best_nonself_spearman,
    nonself_median      = median(S_hvg[cell_line, setdiff(colnames(S_hvg), depmap_ACH)]),
    nonself_sd          = sd(S_hvg[cell_line, setdiff(colnames(S_hvg), depmap_ACH)]),
    self_z_vs_nonself   = (self_spearman - nonself_median) / nonself_sd,
    # reciprocal direction: among OUR 31 lines, is this line the DepMap line's best?
    depmap_best_our     = rownames(S_hvg)[which.max(S_hvg[, depmap_ACH])],
    depmap_best_our_2nd = rownames(S_hvg)[order(-S_hvg[, depmap_ACH])[2]],
    reciprocal_best     = depmap_best_our == cell_line,
    reciprocal_margin   = S_hvg[cell_line, depmap_ACH] -
                          S_hvg[depmap_best_our_2nd, depmap_ACH],
    # sensitivity: same statistics over ALL shared genes (no HVG selection)
    self_spearman_allgenes = S_all[cell_line, depmap_ACH],
    self_rank_of67_allgenes = as.integer(rank(-S_all[cell_line, ])[depmap_ACH])) %>%
  ungroup() %>%
  mutate(best_nonself_name = ach2name[best_nonself_ACH],
         best_nonself_subtype = ach2sub[best_nonself_ACH],
         n_hvg = n_hvg, n_shared_entrez = length(common),
         depmap_release = DEPMAP_RELEASE,
         collapse_note = "ours: sum linear TPM per Entrez then log2(x+1); DepMap: mean of log2(TPM+1) — asymmetric by construction",
         hvg_note = "HVGs selected on the COMBINED 31+67 matrix; inflates absolute rho for self AND non-self alike") %>%
  select(cell_line, subtype, depmap_ACH, depmap_name, depmap_RRID,
         depmap_oncotree_subtype, self_spearman, self_rank_of67, n_depmap_ovarian,
         best_nonself_name, best_nonself_ACH, best_nonself_subtype,
         best_nonself_spearman, margin_self_minus_best_nonself,
         nonself_median, nonself_sd, self_z_vs_nonself,
         reciprocal_best, depmap_best_our, depmap_best_our_2nd, reciprocal_margin,
         self_spearman_allgenes, self_rank_of67_allgenes,
         n_hvg, n_shared_entrez, depmap_release, collapse_note, hvg_note) %>%
  arrange(desc(margin_self_minus_best_nonself))
readr::write_csv(selfmatch, file.path(OUT, "external_selfmatch_margin.csv"))

all67 <- as.data.frame(S_hvg) %>% rownames_to_column("our_line") %>%
  pivot_longer(-our_line, names_to = "depmap_ACH", values_to = "spearman_hvg") %>%
  mutate(spearman_allgenes = S_all[cbind(our_line, depmap_ACH)],
         depmap_name    = ach2name[depmap_ACH],
         depmap_subtype = ach2sub[depmap_ACH],
         our_subtype    = ours$subtype[match(our_line, ours$cell_line)],
         self_ACH       = setNames(overlap$depmap_ACH, overlap$cell_line)[our_line],
         is_selfpair    = !is.na(self_ACH) & depmap_ACH == self_ACH,
         our_line_has_depmap_namesake = our_line %in% overlap$cell_line) %>%
  select(-self_ACH) %>%
  arrange(our_line, desc(spearman_hvg))
readr::write_csv(all67, file.path(OUT, "external_depmap_spearman_all.csv"))

message("\n== SELF-MATCH MARGIN vs all ", ncol(S_hvg), " DepMap ovarian lines ==")
print(as.data.frame(selfmatch %>%
  select(cell_line, depmap_name, self_spearman, self_rank_of67, best_nonself_name,
         best_nonself_spearman, margin_self_minus_best_nonself, self_z_vs_nonself,
         reciprocal_best, reciprocal_margin) %>%
  mutate(across(where(is.numeric), ~round(., 3)))))
message("  wrote output/external_selfmatch_margin.csv (", nrow(selfmatch), " overlap lines)")
message("  wrote output/external_depmap_spearman_all.csv (", nrow(all67), " our-line x DepMap-line pairs)")
stopifnot("a self-match is no longer rank 1 of 67" = all(selfmatch$self_rank_of67 == 1L),
          "a self-match is no longer reciprocal-best" = all(selfmatch$reciprocal_best),
          "a self-match margin is non-positive" =
            all(selfmatch$margin_self_minus_best_nonself > 0))

# =============================================================================
# PART 1c — DRIVER MUTATION CROSS-CHECK (DepMap somatic mutation matrices)
# =============================================================================
read_mut_matrix <- function(path) {
  m <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  colnames(m)[1] <- "ModelID"
  rownames(m) <- m$ModelID; m$ModelID <- NULL
  colnames(m) <- sub(" \\([0-9]+\\)$", "", colnames(m))     # strip "(entrez)"
  m
}
hot <- read_mut_matrix(file.path(EXT, "depmap_hotspot_overlap5.csv"))
dam <- read_mut_matrix(file.path(EXT, "depmap_damaging_overlap5.csv"))
key_genes <- c("TP53","SMARCA4","SMARCA2","ARID1A","PIK3CA","KRAS","CTNNB1","PTEN","SMAD4","BRAF","PPP2R1A")

drv <- overlap %>% select(cell_line, depmap_ACH) %>%
  rowwise() %>%
  mutate(depmap_total_damaging = sum(as.numeric(unlist(dam[depmap_ACH, ])), na.rm = TRUE),
         depmap_total_hotspot  = sum(as.numeric(unlist(hot[depmap_ACH, ])), na.rm = TRUE)) %>%
  ungroup()
message("\n== DepMap mutation burden (damaging count per line; TOV21G expected high = MSI/hypermutation) ==")
print(as.data.frame(drv))

gene_status <- function(ach, gene) {
  h <- if (gene %in% colnames(hot)) as.numeric(hot[ach, gene]) else NA
  d <- if (gene %in% colnames(dam)) as.numeric(dam[ach, gene]) else NA
  paste0("hot=", ifelse(is.na(h),".",h), "/dam=", ifelse(is.na(d),".",d))
}
drv_tab <- expand_grid(cell_line = overlap$cell_line, gene = key_genes) %>%
  left_join(overlap %>% select(cell_line, depmap_ACH), by = "cell_line") %>%
  rowwise() %>% mutate(depmap = gene_status(depmap_ACH, gene)) %>% ungroup() %>%
  select(cell_line, gene, depmap) %>%
  pivot_wider(names_from = cell_line, values_from = depmap)
message("\n== DepMap driver status (hotspot/damaging counts) for overlap lines ==")
print(as.data.frame(drv_tab))

## our WES driver calls for the 3 overlap lines that have WES
stopifnot("output/wes_mutations_filtered.csv missing — run 07_wes_mutations.R first" =
            file.exists(file.path(OUT, "wes_mutations_filtered.csv")))
w <- readr::read_csv(file.path(OUT, "wes_mutations_filtered.csv"), show_col_types = FALSE)
our_drv <- w %>% filter(cell_line %in% overlap$cell_line, Hugo_Symbol %in% key_genes) %>%
  select(cell_line, Hugo_Symbol, HGVSp_Short, Variant_Classification, vaf, is_driver)
message("\n== OUR WES driver calls (overlap lines with WES: OV90/TOV21G/TOV112D) ==")
print(as.data.frame(our_drv))
our_burden <- as.data.frame(table(cell_line = w$cell_line[w$cell_line %in% overlap$cell_line]),
                             responseName = "our_PASS_variants", stringsAsFactors = FALSE)
message("\n== OUR WES PASS-variant burden (overlap lines) ==")
print(as.data.frame(our_burden))

# --- WRITE THE CROSS-CHECK TO A FILE [review revision] -----------------------
# The manuscript cites specific results from this comparison (TOV112D TP53 R175H +
# SMARCA4 damaging, OV90 SMAD4, all five TOV21G clear-cell drivers corroborated;
# TOV112D KRAS-A59T and OV90's BRAF in-frame indel NOT corroborated) but the
# comparison existed only as console output — untraceable in the deposited record.
# One row per (line, gene): our call, DepMap's hotspot/damaging counts, and an
# explicit agreement verdict.
tiers_f <- file.path(OUT, "wes_driver_tiers.csv")
tiers <- if (file.exists(tiers_f)) {
  readr::read_csv(tiers_f, show_col_types = FALSE) %>%
    group_by(cell_line, gene) %>%
    summarise(our_tier = paste(sort(unique(tier)), collapse = ";"), .groups = "drop")
} else {
  tibble(cell_line = character(), gene = character(), our_tier = character())
}

xcheck <- expand_grid(cell_line = overlap$cell_line, gene = key_genes) %>%
  left_join(overlap %>% select(cell_line, depmap_ACH, depmap_name), by = "cell_line") %>%
  left_join(w %>% filter(Hugo_Symbol %in% key_genes) %>%
              group_by(cell_line, gene = Hugo_Symbol) %>%
              summarise(our_variants = paste(paste0(HGVSp_Short, " (", Variant_Classification, ")"),
                                             collapse = " ; "),
                        our_n_calls = dplyr::n(), .groups = "drop"),
            by = c("cell_line","gene")) %>%
  left_join(tiers, by = c("cell_line","gene")) %>%
  rowwise() %>%
  mutate(depmap_hotspot  = if (gene %in% colnames(hot)) as.numeric(hot[depmap_ACH, gene]) else NA_real_,
         depmap_damaging = if (gene %in% colnames(dam)) as.numeric(dam[depmap_ACH, gene]) else NA_real_) %>%
  ungroup() %>%
  mutate(our_has_wes = cell_line %in% w$cell_line,
         our_call    = replace_na(our_n_calls, 0L) > 0,
         depmap_call = replace_na(depmap_hotspot, 0) > 0 | replace_na(depmap_damaging, 0) > 0,
         agreement = case_when(
           !our_has_wes              ~ "no WES for this line (RNA/protein only)",
           our_call &  depmap_call   ~ "corroborated (both)",
           our_call & !depmap_call   ~ "OURS ONLY — not corroborated by DepMap",
           !our_call &  depmap_call  ~ "DEPMAP ONLY — we do not call it",
           TRUE                      ~ "neither")) %>%
  select(cell_line, depmap_name, depmap_ACH, gene, our_has_wes, our_variants,
         our_n_calls, our_tier, depmap_hotspot, depmap_damaging, agreement) %>%
  arrange(cell_line, match(gene, key_genes))
xcheck$depmap_release <- DEPMAP_RELEASE
burden_x <- overlap %>% select(cell_line, depmap_ACH) %>%
  left_join(drv %>% select(cell_line, depmap_total_damaging, depmap_total_hotspot),
            by = "cell_line") %>%
  left_join(our_burden %>% rename(our_PASS_coding_variants = our_PASS_variants),
            by = "cell_line") %>%
  mutate(depmap_release = DEPMAP_RELEASE)
readr::write_csv(xcheck,   file.path(OUT, "external_depmap_driver_crosscheck.csv"))
readr::write_csv(burden_x, file.path(OUT, "external_depmap_burden.csv"))
message("\n== DepMap driver cross-check (written, no longer console-only) ==")
print(as.data.frame(xcheck %>% filter(agreement != "no WES for this line (RNA/protein only)",
                                      agreement != "neither") %>%
  select(cell_line, gene, our_variants, our_tier, depmap_hotspot, depmap_damaging, agreement)),
  row.names = FALSE)
message("  wrote output/external_depmap_driver_crosscheck.csv (", nrow(xcheck), " line x gene rows)")
message("  wrote output/external_depmap_burden.csv (", nrow(burden_x), " overlap lines)")

# =============================================================================
# PART 2 — ConsensusOV subtype calls: COMPUTED, with probabilities and margins
# -----------------------------------------------------------------------------
# [review revision] These calls were previously STRING-EXTRACTED from samples.csv
# free-text `notes` — inherited from an archived analysis with no recorded classifier
# version or settings — while the manuscript claimed they "were generated ... (script
# 10)". Neither part was true. We now RUN consensusOV::get.consensus.subtypes() and
# record the package version, the settings, and the per-class random-forest
# probabilities. The inherited labels are retained as `consensusov_call_inherited`
# and reconciled against the computed ones.
#
# WHY THE MARGIN MATTERS MORE THAN THE LABEL: the paper's own argument is that
# TCGA/ConsensusOV subtypes are substantially TME/stroma-driven and therefore of
# questionable validity on pure tumour-cell cultures. A line called MES with a
# 0.34/0.33 top-vs-second margin makes that point far better than a bare label, so
# the margin and all four class probabilities are emitted per line.
#
# INPUT-SET SENSITIVITY: get.consensus.subtypes() rescales the supplied expression
# matrix internally, so a line's call can depend on WHICH samples are in the matrix.
# We therefore run it twice — on the 15 HGSC RNA lines (the set the manuscript
# reports) and on all 31 RNA lines — and report whether any of the 15 calls move.
# For a claim about label reliability that instability is itself a result.
# =============================================================================
CONSOV_OK <- requireNamespace("consensusOV", quietly = TRUE)
CONSOV_SETTINGS <- "concordant.tumors.only = TRUE; remove.using.cutoff = FALSE (package defaults)"
inherited <- ours %>%
  transmute(cell_line, subtype, source_site,
            has_rna = rna_seq == "Y", applicable = subtype == "HGS",
            consensusov_call_inherited = str_match(notes, "ConsensusOV:([A-Za-z]+)")[, 2])

if (!CONSOV_OK) {
  # Documented fallback: keep the parse, but NAME it as inherited and unreproducible.
  consov <- inherited %>%
    mutate(consensusov_call = NA_character_,
           consensusov_provenance = paste0(
             "INHERITED from metadata/samples.csv `notes`; consensusOV NOT installed in ",
             "this environment, so the calls could NOT be recomputed. Classifier version ",
             "and settings are UNKNOWN. Must be described in the manuscript as inherited ",
             "and unreproducible, not as generated by this pipeline."),
           consensusov_version = NA_character_, consensusov_settings = NA_character_)
  warning("consensusOV is not installed: ConsensusOV calls remain INHERITED and ",
          "unreproducible. Install with BiocManager::install('consensusOV').",
          call. = FALSE)
} else {
  # NB: consensusOV exports a function named `margin`, which would mask
  # ggplot2::margin and break theme_lab() further down. Call it namespaced only —
  # do NOT library(consensusOV).
  CONSOV_VER <- as.character(packageVersion("consensusOV"))
  message("\n== ConsensusOV ", CONSOV_VER, " — computing calls (", CONSOV_SETTINGS, ") ==")

  run_consov <- function(lines, tag) {
    em <- our_log[, intersect(lines, colnames(our_log)), drop = FALSE]
    stopifnot("no RNA columns for the requested ConsensusOV input set" = ncol(em) > 0)
    set.seed(SEED)                       # random forest -> seed matters
    res <- consensusOV::get.consensus.subtypes(em, rownames(em))
    pr  <- as.data.frame(res$rf.probs)
    colnames(pr) <- sub("_consensus$", "", colnames(pr))
    srt <- t(apply(as.matrix(pr), 1, function(r) sort(r, decreasing = TRUE)[1:2]))
    tibble(cell_line = colnames(em),
           consensusov_call = sub("_consensus$", "", as.character(res$consensusOV.subtypes)),
           prob_top = srt[, 1], prob_second = srt[, 2],
           margin_top_vs_second = srt[, 1] - srt[, 2],
           second_call = colnames(pr)[apply(as.matrix(pr), 1,
                                            function(r) order(-r)[2])],
           input_set = tag, n_input_lines = ncol(em)) %>%
      bind_cols(pr %>% rename_with(~paste0("prob_", .x)))
  }
  hgs_rna <- ours %>% filter(subtype == "HGS", rna_seq == "Y") %>% pull(cell_line)
  all_rna <- intersect(ours$cell_line, colnames(our_log))
  co_hgs <- run_consov(hgs_rna, sprintf("HGSC RNA lines (n=%d)", length(hgs_rna)))
  co_all <- run_consov(all_rna, sprintf("all RNA lines (n=%d)", length(all_rna)))

  # input-set stability for the 15 reported lines
  stab <- co_hgs %>% select(cell_line, call_hgs_set = consensusov_call,
                            margin_hgs_set = margin_top_vs_second) %>%
    left_join(co_all %>% select(cell_line, call_all_set = consensusov_call,
                                margin_all_set = margin_top_vs_second), by = "cell_line") %>%
    mutate(call_stable_across_input_sets = call_hgs_set == call_all_set)

  consov <- inherited %>%
    left_join(co_hgs %>% select(-input_set, -n_input_lines), by = "cell_line") %>%
    left_join(stab %>% select(cell_line, call_all_set, margin_all_set,
                              call_stable_across_input_sets), by = "cell_line") %>%
    mutate(agrees_with_inherited = ifelse(is.na(consensusov_call) | is.na(consensusov_call_inherited),
                                          NA, consensusov_call == consensusov_call_inherited),
           consensusov_version = CONSOV_VER,
           consensusov_settings = CONSOV_SETTINGS,
           consensusov_provenance = paste0(
             "COMPUTED by scripts/18_external_benchmarking.R with consensusOV ", CONSOV_VER,
             " on log2(TPM+1) Entrez-collapsed expression of the ", length(hgs_rna),
             " HGSC RNA lines; seed ", SEED, "; ", CONSOV_SETTINGS,
             ". Inherited labels from samples.csv `notes` retained for comparison."))

  cat("\n=== ConsensusOV COMPUTED calls, probabilities and margins (HGSC RNA lines) ===\n")
  print(as.data.frame(consov %>% filter(!is.na(consensusov_call)) %>%
    select(cell_line, consensusov_call, prob_top, second_call, prob_second,
           margin_top_vs_second, consensusov_call_inherited, agrees_with_inherited,
           call_all_set, call_stable_across_input_sets) %>%
    arrange(margin_top_vs_second) %>%
    mutate(across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

  cmp <- consov %>% filter(!is.na(consensusov_call))
  cat(sprintf("\n  computed:  %s\n  inherited: %s\n",
              paste(sprintf("%s %d", names(table(cmp$consensusov_call)),
                            as.integer(table(cmp$consensusov_call))), collapse = " / "),
              paste(sprintf("%s %d", names(table(cmp$consensusov_call_inherited)),
                            as.integer(table(cmp$consensusov_call_inherited))), collapse = " / ")))
  cat(sprintf("  agreement with the inherited labels: %d of %d lines\n",
              sum(cmp$agrees_with_inherited, na.rm = TRUE), nrow(cmp)))
  cat(sprintf("  calls stable when the input set changes from %d HGSC to all %d RNA lines: %d of %d\n",
              length(hgs_rna), length(all_rna),
              sum(cmp$call_stable_across_input_sets, na.rm = TRUE), nrow(cmp)))
  cat(sprintf("  top-vs-second margin: median %.3f, min %.3f (%s) — a small margin means the\n    label is a coin-flip between two classes, which is the paper's point.\n",
              median(cmp$margin_top_vs_second), min(cmp$margin_top_vs_second),
              cmp$cell_line[which.min(cmp$margin_top_vs_second)]))
}
readr::write_csv(consov, file.path(OUT, "consensusov_calls.csv"))
message("  wrote output/consensusov_calls.csv (", nrow(consov), " lines, ",
        sum(!is.na(consov$consensusov_call)), " with a computed call)")

# =============================================================================
# PART 3 — Cellosaurus accessions + STR / problematic-flag status (all 42)
# =============================================================================
parse_cello <- function(name, rrid_hint = NA) {
  safe <- gsub("[/ ]", "_", name)
  f <- file.path(CELLO, paste0("search_", safe, ".json"))
  base <- tibble(cell_line = name, in_cellosaurus = FALSE,
                 cellosaurus_accession = NA_character_, cellosaurus_identifier = NA_character_,
                 str_profile_documented = FALSE, n_str_markers = 0L, str_source = NA_character_,
                 problematic_flag = NA_character_, match_confidence = "not found",
                 rrid_matches_depmap = NA)
  if (!file.exists(f)) return(base)
  j <- tryCatch(jsonlite::read_json(f), error = function(e) NULL)
  cl <- j$Cellosaurus$`cell-line-list`
  if (is.null(cl) || length(cl) == 0) return(base)

  # EXACT NAME MATCH ONLY [review revision]. This previously fell back to the search
  # TOP HIT when no name matched, which can attach a WRONG Cellosaurus accession
  # (and therefore a wrong STR/problematic-flag status) to a line. All 30 current
  # hits are exact matches, so nothing is lost by refusing the fallback; a future
  # rename or a new line now returns "not found" rather than a plausible-looking
  # wrong accession, and the assertion after this function fails loudly.
  target <- norm_nm(name)
  pick <- NA; conf <- NA
  for (i in seq_along(cl)) {
    nms <- vapply(cl[[i]]$`name-list`, function(z) z$value, character(1))
    if (target %in% norm_nm(nms)) { pick <- i; conf <- "high (exact name match)"; break }
  }
  if (is.na(pick)) {
    base$match_confidence <- sprintf("not found (no exact name match among %d search hits; top hit NOT used)",
                                     length(cl))
    return(base)
  }
  e <- cl[[pick]]

  acc <- vapply(e$`accession-list`, function(z) z$value,
                character(1))[vapply(e$`accession-list`, function(z) z$type, character(1)) == "primary"][1]
  ident <- vapply(e$`name-list`, function(z) z$value,
                  character(1))[vapply(e$`name-list`, function(z) z$type, character(1)) == "identifier"][1]
  markers <- e$`str-list`$`marker-list`
  nmark <- if (is.null(markers)) 0L else length(markers)
  comments <- e$`comment-list`
  ccat <- vapply(comments, function(z) z$category %||% NA_character_, character(1))
  cval <- vapply(comments, function(z) z$value    %||% NA_character_, character(1))
  prob <- if (any(ccat == "Problematic cell line", na.rm = TRUE))
    paste(cval[ccat == "Problematic cell line"], collapse = " | ") else NA_character_
  strsrc <- if (any(grepl("STR profile", cval), na.rm = TRUE))
    cval[grepl("STR profile", cval)][1] else NA_character_

  tibble(cell_line = name, in_cellosaurus = TRUE,
         cellosaurus_accession = acc, cellosaurus_identifier = ident,
         str_profile_documented = nmark > 0, n_str_markers = as.integer(nmark),
         str_source = strsrc, problematic_flag = prob,
         match_confidence = conf,
         rrid_matches_depmap = if (!is.na(rrid_hint)) isTRUE(acc == unname(rrid_hint)) else NA)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

rrid_hint <- setNames(overlap$depmap_RRID, overlap$cell_line)
cello <- purrr::map_dfr(ours$cell_line, ~ parse_cello(.x, rrid_hint[.x]))
cello <- ours %>% select(cell_line, subtype, source_site) %>% left_join(cello, by = "cell_line") %>%
  mutate(current_stock_str_status = "documentation pending author confirmation",
         current_stock_mycoplasma_status = "documentation pending author confirmation")
readr::write_csv(cello, file.path(OUT, "cellosaurus_str_status.csv"))

# Fail loudly rather than silently carrying a fuzzy accession [review revision].
bad <- cello %>% filter(in_cellosaurus, match_confidence != "high (exact name match)")
if (nrow(bad)) print(as.data.frame(bad %>% select(cell_line, cellosaurus_accession, match_confidence)))
stopifnot("a Cellosaurus accession was attached WITHOUT an exact name match — verify it by hand before using this table" =
            nrow(bad) == 0,
          "an RRID cross-check against DepMap Model.csv failed" =
            all(cello$rrid_matches_depmap[!is.na(cello$rrid_matches_depmap)]))

message("\n== Cellosaurus / STR summary ==")
message("  in Cellosaurus: ", sum(cello$in_cellosaurus), " / 42 ; absent: ",
        paste(cello$cell_line[!cello$in_cellosaurus], collapse = ", "))
message("  with documented STR profile: ", sum(cello$str_profile_documented), " / 42")
message("  flagged 'Problematic cell line':")
print(as.data.frame(cello %>% filter(!is.na(problematic_flag)) %>%
                      select(cell_line, cellosaurus_accession, problematic_flag)))
message("  RRID cross-check vs DepMap (should all be TRUE):")
print(as.data.frame(cello %>% filter(!is.na(rrid_matches_depmap)) %>%
                      select(cell_line, cellosaurus_accession, rrid_matches_depmap)))

# =============================================================================
# FIGURE — f_external_concordance
#   A: 5 DepMap overlap lines (rows) x 31 our lines (cols) Spearman heatmap;
#      the correct (self) cell outlined in rust — each DepMap line's brightest
#      column is its namesake.
#   B: each overlap line's Spearman to all 67 DepMap ovarian lines; self in rust.
# =============================================================================
sub_order <- c("HGS","CC","EC","MC","MMMT","SCCOHT")
col_ord <- ours %>% filter(cell_line %in% rownames(S_hvg)) %>%
  mutate(subtype = factor(subtype, levels = sub_order)) %>%
  arrange(subtype, cell_line) %>% pull(cell_line)

hmA <- long %>%
  mutate(our_line = factor(our_line, levels = col_ord),
         depmap_lab = paste0(depmap_name, "\n(", depmap_ACH, ")"),
         depmap_lab = factor(depmap_lab,
                             levels = paste0(overlap$depmap_name, "\n(", overlap$depmap_ACH, ")")))
selfcells <- hmA %>% filter(is_selfpair)

pA <- ggplot(hmA, aes(our_line, depmap_lab, fill = spearman_hvg)) +
  geom_tile(colour = "grey92", linewidth = 0.2) +
  geom_tile(data = selfcells, colour = COOK_RUST, linewidth = 1.1, fill = NA, width = 0.98, height = 0.98) +
  geom_text(data = selfcells, aes(label = sprintf("%.2f", spearman_hvg)),
            colour = COOK_RUST, fontface = "bold", size = 3) +
  scale_fill_gradient(low = "white", high = COOK_NAVY, limits = c(0, 1),
                      name = "Spearman\n(2000 HVG)") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(title = "A  Cross-dataset expression identity: our RNA vs DepMap 24Q4",
       subtitle = "Each DepMap overlap line (row) matches its namesake among our 31 lines (rust outline); columns grouped by our subtype",
       x = "This resource: RNA-seq lines (grouped by subtype)", y = "DepMap ovarian line") +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.position = "right", panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4))

# Panel B: self vs all 67 DepMap ovarian
allpairs <- as.data.frame(S_hvg) %>% rownames_to_column("our_line") %>%
  pivot_longer(-our_line, names_to = "depmap_ACH", values_to = "spearman") %>%
  filter(our_line %in% overlap$cell_line) %>%
  left_join(overlap %>% select(our_line = cell_line, self_ACH = depmap_ACH), by = "our_line") %>%
  mutate(is_self = depmap_ACH == self_ACH,
         our_line = factor(our_line, levels = overlap$cell_line))
pB <- ggplot(allpairs, aes(our_line, spearman)) +
  geom_jitter(data = filter(allpairs, !is_self), width = 0.18, height = 0,
              colour = "grey70", size = 1, alpha = 0.7) +
  geom_point(data = filter(allpairs, is_self), colour = COOK_RUST, size = 3.2) +
  labs(title = "B  Specificity: each line vs all 67 DepMap ovarian lines",
       subtitle = "Rust = correct DepMap line (self); grey = the other 66 ovarian lines",
       x = NULL, y = "Spearman correlation (2000 HVG)") +
  theme_lab() + theme(axis.text.x = element_text(face = "bold"))

fig <- pA / pB + patchwork::plot_layout(heights = c(1.15, 1))
save_fig(fig, "f_external_concordance", w = 9.5, h = 8.8)

message("\nDONE 18_external_benchmarking.R")
