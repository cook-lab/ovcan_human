# =============================================================================
# Script: 10_authentication.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Cell-line IDENTITY / histotype authentication (Phase 4), integrating
#          the completed RNA (Phase 1), proteomics (Phase 2) and WES (Phase 3)
#          layers. Two deliverables:
#            (1) a MULTI-OMIC SWI/SNF-deficiency panel (SMARCA4/SMARCA2/ARID1A/
#                SMARCB1 across RNA expression, protein abundance, WES mutation)
#                that independently tests the Karnezis-2021 reclassifications
#                COV434 -> SCCOHT and TOV112D -> dedifferentiated carcinoma, and
#            (2) a per-line histotype-consistency table (does each line's
#                molecular profile agree with its labeled subtype?) carrying a
#                populated STR_status column: the Cellosaurus reference profile
#                where one exists, plus the in-house profile once available.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   4 (authentication) — script 1 of 2 (mucinous authenticity = script 11)
# =============================================================================
# Assay-aware notes (these drive how each layer is weighted):
#   - SWI/SNF loss is NOT read the same way in every layer:
#       * SMARCA4 (BRG1) loss in SCCOHT is largely POST-TRANSCRIPTIONAL
#         (inactivating mutation / protein loss with RETAINED mRNA). BIN67 is the
#         textbook case: its SMARCA4 mRNA is normal, so RNA alone MISSES it;
#         protein and/or WES are the definitive readouts.
#       * SMARCA2 (BRM) loss is EPIGENETIC (promoter silencing), so it shows up
#         at the mRNA level (very low TPM). Here RNA is the sensitive layer.
#     Therefore we require BOTH RNA and protein/WES to authenticate the SWI/SNF
#     group, and interpret each gene in the layer where its loss is detectable.
#   - TMT proteomics reports RELATIVE abundance and is subject to ratio
#     compression / co-isolation interference, which BLUNTS the apparent
#     magnitude of a true protein loss (a null protein rarely reads as -Inf).
#     So we judge protein loss by RANK/z among lines and direction, not absolute
#     fold-change, and say so.
#   - WES is tumor-only (no matched normal); mutation CALLS are trustworthy for
#     canonical, well-supported drivers (used here) but somatic FREQUENCIES are
#     not — we use presence of a canonical driver as supportive evidence, with
#     the limitation stated. Only 22-23 CHUM lines have WES; VOA/BIN67/COV434 do
#     not, so their genomics columns are NA and rely on RNA+protein.
#   - Expression "consistency" for well-sampled subtypes with POSITIVE markers
#     (HGS/CC/MC/MMMT) is scored quantitatively; EC has no robust positive panel
#     and SCCOHT is loss-defined, so those are called from the SWI/SNF panel +
#     genomics and flagged as such rather than force-fit to a marker score.
# =============================================================================

source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(matrixStats)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count  <- dplyr::count; slice <- dplyr::slice
set.seed(SEED)

# CB-safe status palette (Okabe-Ito subset) for the authentication tiles
STATUS_COLS <- c(consistent = "#0072B2", discordant = "#D55E00",
                 partial = "#E69F00", "not assessed" = "#BDBDBD")
subtype_lvls <- c("HGS","LGS","CC","EC","MC","MMMT","SCCOHT")
sub_cols <- setNames(RColorBrewer::brewer.pal(8, "Dark2")[c(1,7,2,3,4,5,6)], subtype_lvls)
# (HGS/CC/EC/MC/MMMT/SCCOHT keep their script-01/02/04 Dark2 colours; LGS = 7th)

# -----------------------------------------------------------------------------
# 0. Master line roster (generated lines only) + per-assay availability
# -----------------------------------------------------------------------------
samples <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE)
gen <- samples %>%
  filter(provenance == "generated") %>%
  transmute(
    cell_line, labeled_subtype = subtype, subtype_status, source_site,
    has_rna  = rna_seq == "Y",
    has_prot = proteomics == "Y",
    has_wes_cnv = wes_cnv == "Y",
    has_wes_mut = wes_mut == "Y",
    notes) %>%
  mutate(labeled_subtype = factor(labeled_subtype, levels = subtype_lvls))
message(sprintf("Generated lines: %d (RNA %d, protein %d, WES-CNV %d, WES-mut %d)",
                nrow(gen), sum(gen$has_rna), sum(gen$has_prot),
                sum(gen$has_wes_cnv), sum(gen$has_wes_mut)))

# -----------------------------------------------------------------------------
# 1. Load the three omic layers
# -----------------------------------------------------------------------------
## RNA: TPM -> symbol-collapsed log2 matrix (same recipe as script 04) --------
tpm  <- readr::read_csv(file.path(OUT, "rna_tpm.csv"), show_col_types = FALSE)
t2g  <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"), show_col_types = FALSE)
g2e  <- t2g %>% distinct(ensembl_gene_id, external_gene_name) %>%
  filter(!is.na(external_gene_name), external_gene_name != "")
m_id <- as.matrix(tpm[, -1]); rownames(m_id) <- tpm$gene_id
id2sym <- setNames(g2e$external_gene_name, g2e$ensembl_gene_id)
common <- intersect(rownames(m_id), names(id2sym))
tpm_sym <- rowsum(m_id[common, , drop = FALSE], group = id2sym[common])  # symbols x lines
logtpm  <- log2(tpm_sym + 1)
rna_lines <- colnames(logtpm)
zscore <- function(mat) t(scale(t(mat)))                 # z across lines, per gene

## Protein: log2 relative abundance, symbols x lines --------------------------
# [integration revision] Read through the shared loader in 00_setup.R. This script
# previously carried its own `filter(!is.na(protein), !duplicated(protein))` guard
# written against the OLD matrix shape (8,430 rows, no-symbol rows named "NA"/"NA.1",
# duplicate symbols named SYMBOL.1). Script 05 now handles all of that at source —
# 8,427 rows, unique names, non-representatives named SYMBOL|UNIPROT — so the local
# guard was dead code that hid the change. read_prot_matrix() also derives the 70
# zero-plex rows (NA in all 31 lines) and asserts them against prot_qc.csv.
prot <- read_prot_matrix()
prot_zero_plex <- attr(prot, "zero_plex")
prot_lines <- colnames(prot)

## WES mutations (filtered, canonical-driver-aware) ---------------------------
stopifnot("output/wes_mutations_filtered.csv missing — run 07_wes_mutations.R first" =
            file.exists(file.path(OUT, "wes_mutations_filtered.csv")))
wm <- readr::read_csv(file.path(OUT, "wes_mutations_filtered.csv"), show_col_types = FALSE)
wes_lines <- sort(unique(wm$cell_line))
trunc_classes <- c("Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation",
                   "Splice_Site","Nonstop_Mutation","Translation_Start_Site")

## Somatic-confidence tiers per (line, gene) — needed so a Tier-3 "cannot exclude
## germline" call cannot silently drive an authentication verdict. [review revision]
tiers <- readr::read_csv(file.path(OUT, "wes_driver_tiers.csv"), show_col_types = FALSE)
tier_min <- tiers %>% group_by(cell_line, gene) %>%
  summarise(tier_best = min(as.integer(sub("Tier", "", tier))),
            tiers_all = paste(sort(unique(tier)), collapse = ";"),
            .groups = "drop")

# Variant-level tiers prevent an unqualified truncation inheriting the tier of
# another candidate in the same gene/model.
tier_key <- c("cell_line", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
tier_variant <- tiers %>% transmute(across(all_of(tier_key)), Hugo_Symbol = gene,
  tier_best = as.integer(sub("Tier", "", tier)), tiers_all = tier)
stopifnot(!anyDuplicated(tier_variant[, c(tier_key, "Hugo_Symbol")]))

## CNV fraction-genome-altered
fga <- readr::read_csv(file.path(OUT, "wes_cnv_fga.csv"), show_col_types = FALSE)

# =============================================================================
# SECTION 1 — MULTI-OMIC SWI/SNF-DEFICIENCY PANEL
# =============================================================================
swi_genes <- c("SMARCA4","SMARCA2","ARID1A","SMARCB1")
key_lines <- c("TOV112D","COV434","BIN67")
stopifnot(all(swi_genes %in% rownames(logtpm)), all(swi_genes %in% rownames(prot)))
# A zero-plex protein is NA in all 31 lines, so a z-score or rank built from it would
# be an artefact of missingness, not a measurement. None of the SWI/SNF panel is
# zero-plex; assert it rather than trust it. [integration revision]
stopifnot("a SWI/SNF panel protein is zero-plex (NA in all 31 lines) — its protein z/rank would be meaningless" =
            !any(swi_genes %in% prot_zero_plex))

## --- RNA z + rank (rank 1 = lowest of the RNA lines) ------------------------
rna_z <- zscore(logtpm[swi_genes, , drop = FALSE])
rna_rank <- t(apply(logtpm[swi_genes, , drop = FALSE], 1,
                    function(x) rank(x, ties.method = "min")))
n_rna <- length(rna_lines)

## --- Protein z + rank -------------------------------------------------------
prot_z <- zscore(prot[swi_genes, , drop = FALSE])
prot_rank <- t(apply(prot[swi_genes, , drop = FALSE], 1,
                     function(x) rank(x, ties.method = "min", na.last = "keep")))
n_prot <- length(prot_lines)

## --- WES mutation status per gene per line ----------------------------------
# [review revision] The somatic-confidence TIER now travels with each call. The
# deficiency logic below previously used truncating status irrespective of tier, so
# Tier-3 "cannot exclude germline" SMARCA2 calls (OV2295, TOV21G) fed
# swisnf_deficient — i.e. an authentication verdict rested on evidence the pipeline
# itself declares not defensibly somatic.
swi_mut <- wm %>%
  filter(Hugo_Symbol %in% swi_genes) %>%
  mutate(mclass = ifelse(Variant_Classification %in% trunc_classes,
                         "truncating",
                         ifelse(Variant_Classification == "Missense_Mutation",
                                "missense", "inframe/other"))) %>%
  left_join(tier_variant, by = c(tier_key, "Hugo_Symbol")) %>%
  group_by(cell_line, Hugo_Symbol) %>%
  # if multiple hits, keep the most damaging (truncating > missense > other)
  summarise(trunc_defensible = any(mclass == "truncating" & !is.na(tier_best) & tier_best <= 2L),
    mclass = c("truncating","missense","inframe/other")[
    min(match(mclass, c("truncating","missense","inframe/other")))],
    hgvsp = paste(unique(gsub("^p\\.", "", HGVSp_Short)), collapse = ";"),
    tier_best = suppressWarnings(min(tier_best, na.rm = TRUE)),
    tiers_all = paste(sort(unique(na.omit(tiers_all))), collapse = ";"),
    .groups = "drop") %>%
  mutate(tier_best = ifelse(is.finite(tier_best), tier_best, NA_integer_))

wes_status <- function(cl, g) {
  if (!cl %in% wes_lines) return(NA_character_)            # no WES for this line
  hit <- swi_mut %>% filter(cell_line == cl, Hugo_Symbol == g)
  if (nrow(hit) == 0) return("no retained candidate")
  paste0(hit$mclass[1], " (", hit$hgvsp[1], ")",
         if (!is.na(hit$tier_best[1])) sprintf(" [Tier%d]", hit$tier_best[1]) else "")
}
# tier-aware truncation lookup used by the deficiency calls
wes_trunc_t12 <- function(cl, g) {
  if (!cl %in% wes_lines) return(NA)
  hit <- swi_mut %>% filter(cell_line == cl, Hugo_Symbol == g)
  if (nrow(hit) == 0) return(FALSE)
  isTRUE(hit$trunc_defensible[1])
}

## --- Assemble the wide per-line SWI/SNF panel (all generated lines) ----------
gv <- function(mat, g, cl) if (cl %in% colnames(mat)) mat[g, cl] else NA_real_
swisnf_panel <- gen %>%
  select(cell_line, subtype = labeled_subtype, source_site,
         has_rna, has_prot) %>%
  mutate(has_wes = cell_line %in% wes_lines) %>%
  rowwise() %>%
  mutate(
    SMARCA4_rna_z  = round(gv(rna_z,  "SMARCA4", cell_line), 2),
    SMARCA4_rna_rank = gv(rna_rank, "SMARCA4", cell_line),
    SMARCA4_prot_z = round(gv(prot_z, "SMARCA4", cell_line), 2),
    SMARCA4_prot_rank = gv(prot_rank, "SMARCA4", cell_line),
    SMARCA4_wes    = wes_status(cell_line, "SMARCA4"),
    SMARCA4_prot_rank_of  = n_prot,
    SMARCA4_rna_rank_of   = n_rna,
    SMARCA2_rna_z  = round(gv(rna_z,  "SMARCA2", cell_line), 2),
    SMARCA2_rna_rank = gv(rna_rank, "SMARCA2", cell_line),
    SMARCA2_prot_z = round(gv(prot_z, "SMARCA2", cell_line), 2),
    SMARCA2_prot_rank = gv(prot_rank, "SMARCA2", cell_line),
    SMARCA2_wes    = wes_status(cell_line, "SMARCA2"),
    ARID1A_rna_z   = round(gv(rna_z, "ARID1A", cell_line), 2),
    ARID1A_rna_rank  = gv(rna_rank,  "ARID1A", cell_line),
    ARID1A_prot_z  = round(gv(prot_z, "ARID1A", cell_line), 2),
    ARID1A_prot_rank = gv(prot_rank, "ARID1A", cell_line),
    ARID1A_wes     = wes_status(cell_line, "ARID1A"),
    SMARCB1_rna_z  = round(gv(rna_z, "SMARCB1", cell_line), 2),
    SMARCB1_rna_rank  = gv(rna_rank,  "SMARCB1", cell_line),
    SMARCB1_prot_z = round(gv(prot_z, "SMARCB1", cell_line), 2),
    SMARCB1_prot_rank = gv(prot_rank, "SMARCB1", cell_line),
    SMARCB1_wes    = wes_status(cell_line, "SMARCB1"),
    # raw (un-z'd) values so the panel is quotable without the z reference set
    SMARCA4_rna_log2tpm  = round(gv(logtpm, "SMARCA4", cell_line), 2),
    SMARCA4_prot_log2    = round(gv(prot,   "SMARCA4", cell_line), 2),
    SMARCA2_rna_log2tpm  = round(gv(logtpm, "SMARCA2", cell_line), 2),
    SMARCA2_prot_log2    = round(gv(prot,   "SMARCA2", cell_line), 2),
    # tier-aware truncation evidence (see swi_mut above)
    SMARCA4_trunc_tier12 = wes_trunc_t12(cell_line, "SMARCA4"),
    SMARCA2_trunc_tier12 = wes_trunc_t12(cell_line, "SMARCA2"),
    ARID1A_trunc_tier12  = wes_trunc_t12(cell_line, "ARID1A"),
    SMARCB1_trunc_tier12 = wes_trunc_t12(cell_line, "SMARCB1")
  ) %>% ungroup()

## --- Deficiency call (evidence-based, per gene, layer-appropriate) ----------
# SMARCA4 loss: TIER1-2 truncating OR protein z <= -1 (post-transcriptional; RNA not required)
# SMARCA2 loss: RNA very low (rank <= 4 of 31 OR z <= -1.5; epigenetic silencing) OR TIER1-2 truncating
# ARID1A/SMARCB1 loss: TIER1-2 truncating OR (protein z <= -1.5 AND rna z <= -1) [rarely the driver here]
# [review revision] "truncating" is now TIER-RESTRICTED (see swi_mut). A Tier-3 call
# means germline cannot be excluded in tumour-only WES, so it is recorded in the
# panel (*_wes carries the tier) but is NOT admitted as deficiency evidence. Tier-3
# truncating calls that were previously admitted: OV2295 SMARCA2, TOV21G SMARCA2.
is_trunc <- function(s) !is.na(s) & grepl("^truncating", s)
isTRUE_v <- function(x) !is.na(x) & x          # NA (no WES) is not evidence
swisnf_panel <- swisnf_panel %>%
  mutate(
    SMARCA4_loss = isTRUE_v(SMARCA4_trunc_tier12) |
      (!is.na(SMARCA4_prot_z) & SMARCA4_prot_z <= -1),
    SMARCA2_loss = isTRUE_v(SMARCA2_trunc_tier12) |
      (!is.na(SMARCA2_rna_rank) & SMARCA2_rna_rank <= 4) |
      (!is.na(SMARCA2_rna_z) & SMARCA2_rna_z <= -1.5),
    ARID1A_loss  = isTRUE_v(ARID1A_trunc_tier12),
    SMARCB1_loss = isTRUE_v(SMARCB1_trunc_tier12),
    swisnf_deficient = SMARCA4_loss | SMARCA2_loss | ARID1A_loss | SMARCB1_loss,
    # transparency: which lines carry a truncating call that is TIER-3 ONLY and is
    # therefore recorded but NOT counted as evidence
    swisnf_tier3_only_calls = {
      f <- function(g, w, t) if (is_trunc(w) && !isTRUE(t)) g else NULL
      purrr::pmap_chr(list(SMARCA4_wes, SMARCA4_trunc_tier12,
                           SMARCA2_wes, SMARCA2_trunc_tier12,
                           ARID1A_wes,  ARID1A_trunc_tier12,
                           SMARCB1_wes, SMARCB1_trunc_tier12),
        function(a4w, a4t, a2w, a2t, arw, art, b1w, b1t) {
          v <- c(f("SMARCA4", a4w, a4t), f("SMARCA2", a2w, a2t),
                 f("ARID1A", arw, art), f("SMARCB1", b1w, b1t))
          if (length(v)) paste(v, collapse = ", ") else NA_character_
        })
    })

# human-readable evidence string
ev <- function(r) {
  e <- c()
  if (r$SMARCA4_loss) e <- c(e, paste0("SMARCA4 loss[",
      paste(c(if (isTRUE(r$SMARCA4_trunc_tier12)) "WES-trunc(Tier1-2)",
              if (!is.na(r$SMARCA4_prot_z) && r$SMARCA4_prot_z <= -1)
                sprintf("prot z=%.1f (rank %s/%s)", r$SMARCA4_prot_z,
                        r$SMARCA4_prot_rank, r$SMARCA4_prot_rank_of)), collapse = "+"), "]"))
  if (r$SMARCA2_loss) e <- c(e, paste0("SMARCA2 loss[",
      paste(c(if (isTRUE(r$SMARCA2_trunc_tier12)) "WES-trunc(Tier1-2)",
              if (!is.na(r$SMARCA2_rna_rank) && r$SMARCA2_rna_rank <= 4)
                sprintf("RNA rank=%d/%d", r$SMARCA2_rna_rank, r$SMARCA4_rna_rank_of)),
            collapse = "+"), "]"))
  if (r$ARID1A_loss)  e <- c(e, "ARID1A truncating(WES, Tier1-2)")
  if (r$SMARCB1_loss) e <- c(e, "SMARCB1 truncating(WES, Tier1-2)")
  if (!is.na(r$swisnf_tier3_only_calls))
    e <- c(e, paste0("[NOT counted: Tier-3 truncating call in ",
                     r$swisnf_tier3_only_calls, " — germline not excludable]"))
  if (length(e) == 0) "SWI/SNF intact (by available layers)" else paste(e, collapse = "; ")
}
swisnf_panel$swisnf_evidence <- vapply(seq_len(nrow(swisnf_panel)),
                                      function(i) ev(swisnf_panel[i, ]), character(1))
swisnf_panel <- swisnf_panel %>% arrange(subtype, cell_line)
readr::write_csv(swisnf_panel, file.path(OUT, "auth_swisnf_panel.csv"))

cat("\n=== SWI/SNF panel: key lines (TOV112D, COV434, BIN67) ===\n")
print(as.data.frame(swisnf_panel %>% filter(cell_line %in% key_lines) %>%
  select(cell_line, subtype,
         SMARCA4_rna_rank, SMARCA4_prot_z, SMARCA4_wes,
         SMARCA2_rna_rank, SMARCA2_prot_z, SMARCA2_wes,
         swisnf_deficient, swisnf_evidence)), row.names = FALSE)
cat("\n=== All SWI/SNF-deficient lines (any layer) ===\n")
print(as.data.frame(swisnf_panel %>% filter(swisnf_deficient) %>%
  select(cell_line, subtype, swisnf_evidence)), row.names = FALSE)
if (any(!is.na(swisnf_panel$swisnf_tier3_only_calls))) {
  cat("\n=== Tier-3-ONLY truncating SWI/SNF calls: recorded but NOT admitted as evidence ===\n")
  print(as.data.frame(swisnf_panel %>% filter(!is.na(swisnf_tier3_only_calls)) %>%
    select(cell_line, subtype, swisnf_tier3_only_calls, swisnf_deficient)), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 1a. PER-LINE SWI/SNF TABLE FOR A FIGURE  [review revision]
# -----------------------------------------------------------------------------
# The manuscript claims two SWI/SNF-based reclassifications (COV434 -> SCCOHT,
# TOV112D -> dedifferentiated carcinoma) and cites Fig. 5 for them, but NO figure in
# the set contains a SWI/SNF panel — SMARCA4/SMARCA2 per-line RNA and protein appear
# nowhere. This long-format table carries everything such a panel needs: per line and
# gene, the raw value, the z, the rank and its denominator, the WES call with its
# somatic-confidence tier, and whether that call was admitted as evidence.
swisnf_long <- bind_rows(
  as_tibble(logtpm[swi_genes, , drop = FALSE], rownames = "gene") %>%
    pivot_longer(-gene, names_to = "cell_line", values_to = "value") %>%
    mutate(assay = "RNA", unit = "log2(TPM+1)"),
  as_tibble(prot[swi_genes, , drop = FALSE], rownames = "gene") %>%
    pivot_longer(-gene, names_to = "cell_line", values_to = "value") %>%
    mutate(assay = "protein", unit = "log2 relative abundance (TMT, PIS ratio)")) %>%
  left_join(bind_rows(
      as_tibble(rna_z[swi_genes, , drop = FALSE], rownames = "gene") %>%
        pivot_longer(-gene, names_to = "cell_line", values_to = "z") %>% mutate(assay = "RNA"),
      as_tibble(prot_z[swi_genes, , drop = FALSE], rownames = "gene") %>%
        pivot_longer(-gene, names_to = "cell_line", values_to = "z") %>% mutate(assay = "protein")),
    by = c("gene","cell_line","assay")) %>%
  left_join(bind_rows(
      as_tibble(rna_rank[swi_genes, , drop = FALSE], rownames = "gene") %>%
        pivot_longer(-gene, names_to = "cell_line", values_to = "rank_low_is_lowest") %>%
        mutate(assay = "RNA", rank_of = n_rna),
      as_tibble(prot_rank[swi_genes, , drop = FALSE], rownames = "gene") %>%
        pivot_longer(-gene, names_to = "cell_line", values_to = "rank_low_is_lowest") %>%
        mutate(assay = "protein", rank_of = n_prot)),
    by = c("gene","cell_line","assay")) %>%
  left_join(swi_mut %>% select(cell_line, gene = Hugo_Symbol, wes_class = mclass,
                               wes_hgvsp = hgvsp, wes_tier_best = tier_best,
                               wes_trunc_admitted = trunc_defensible),
            by = c("cell_line","gene")) %>%
  left_join(gen %>% select(cell_line, subtype = labeled_subtype, source_site), by = "cell_line") %>%
  mutate(has_wes = cell_line %in% wes_lines,
         z_reference_set = ifelse(assay == "RNA",
                                  sprintf("z across all %d RNA lines", n_rna),
                                  sprintf("z across all %d protein lines", n_prot))) %>%
  select(cell_line, subtype, source_site, gene, assay, unit, value, z,
         rank_low_is_lowest, rank_of, z_reference_set, has_wes, wes_class, wes_hgvsp,
         wes_tier_best, wes_trunc_admitted) %>%
  arrange(gene, assay, rank_low_is_lowest)
readr::write_csv(swisnf_long, file.path(OUT, "auth_swisnf_long.csv"))

# The BIN67 case with its actual values and rank margin — the single most quotable
# authentication lesson in the paper ("mRNA retained, protein second-lowest") and
# currently supported only by the word "second-lowest".
cat("\n=== SWI/SNF ranked values for a figure (SMARCA4/SMARCA2; rank 1 = LOWEST) ===\n")
for (g in c("SMARCA4","SMARCA2")) for (a in c("RNA","protein")) {
  s <- swisnf_long %>% filter(gene == g, assay == a, !is.na(value)) %>%
    arrange(rank_low_is_lowest)
  cat(sprintf("  %s %s (n=%d): lowest 4 = %s | %s rank %s (%.2f, z=%.2f); margin to next = %.3f\n",
              g, a, nrow(s), paste(head(s$cell_line, 4), collapse = ", "),
              "BIN67",
              ifelse("BIN67" %in% s$cell_line,
                     as.character(s$rank_low_is_lowest[s$cell_line == "BIN67"]), "NA"),
              ifelse("BIN67" %in% s$cell_line, s$value[s$cell_line == "BIN67"], NA_real_),
              ifelse("BIN67" %in% s$cell_line, s$z[s$cell_line == "BIN67"], NA_real_),
              if ("BIN67" %in% s$cell_line) {
                i <- which(s$cell_line == "BIN67")
                if (i < nrow(s)) s$value[i + 1L] - s$value[i] else NA_real_
              } else NA_real_))
}
message("Wrote output/auth_swisnf_long.csv (", nrow(swisnf_long),
        " line x gene x assay rows) — SWI/SNF figure input")

# -----------------------------------------------------------------------------
# 1b. SWI/SNF multi-omic figure (RNA z + protein z tracks; WES mutations overlaid)
# -----------------------------------------------------------------------------
# lines shown = those with RNA and/or protein (the SWI/SNF story lives here);
# WES-only lines are captured in the CSV. Columns ordered by subtype; the three
# reclassification-candidate lines are highlighted on the x-axis.
show_lines <- swisnf_panel %>% filter(has_rna | has_prot)
line_order <- show_lines %>% arrange(subtype, cell_line) %>% pull(cell_line)

long_z <- bind_rows(
  as_tibble(rna_z,  rownames = "gene") %>% pivot_longer(-gene, names_to = "cell_line",
             values_to = "z") %>% mutate(assay = "RNA (z log2 TPM)"),
  as_tibble(prot_z, rownames = "gene") %>% pivot_longer(-gene, names_to = "cell_line",
             values_to = "z") %>% mutate(assay = "Protein (z log2 abundance)")
) %>%
  filter(cell_line %in% line_order) %>%
  mutate(gene = factor(gene, levels = rev(swi_genes)),
         cell_line = factor(cell_line, levels = line_order),
         assay = factor(assay, levels = c("RNA (z log2 TPM)","Protein (z log2 abundance)")),
         z = pmax(pmin(z, 2), -2))

# WES mutation points (only SWI/SNF hits), placed on both assay facets
mut_pts <- swi_mut %>%
  filter(cell_line %in% line_order) %>%
  transmute(cell_line = factor(cell_line, levels = line_order),
            gene = factor(Hugo_Symbol, levels = rev(swi_genes)),
            mclass = factor(mclass, levels = c("truncating","missense","inframe/other"))) %>%
  tidyr::crossing(assay = factor(c("RNA (z log2 TPM)","Protein (z log2 abundance)"),
                                 levels = c("RNA (z log2 TPM)","Protein (z log2 abundance)")))

xcol <- ifelse(line_order %in% key_lines, COOK_RUST, "grey30")
xface <- ifelse(line_order %in% key_lines, "bold", "plain")

p_swi <- ggplot(long_z, aes(cell_line, gene)) +
  geom_tile(aes(fill = z), colour = "white", linewidth = 0.3) +
  geom_point(data = mut_pts, aes(shape = mclass), size = 2.4, stroke = 0.9,
             colour = "black", fill = "black") +
  scale_fill_gradient2(low = COOK_NAVY, mid = "white", high = COOK_RUST,
                       midpoint = 0, limits = c(-2, 2), name = "z") +
  scale_shape_manual(values = c(truncating = 4, missense = 21, "inframe/other" = 24),
                     name = "WES mutation", drop = FALSE) +
  facet_wrap(~ assay, ncol = 1) +
  labs(title = "Multi-omic SWI/SNF status across OvCAN cell lines",
       subtitle = paste0("SMARCA4 loss = post-transcriptional (protein/WES); SMARCA2 loss = epigenetic (mRNA).  ",
                         "Rust labels = published dedifferentiated/SCCOHT models (TOV112D, COV434, BIN67).\n",
                         "RNA n=", n_rna, " lines; protein n=", n_prot, " lines; WES mutations from ",
                         length(wes_lines), " lines (X = truncating)."),
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   colour = xcol, face = xface, size = 7),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right")
ggsave(file.path(FIGS, "10_swisnf_panel.pdf"), p_swi, width = 11, height = 5.2)

# =============================================================================
# SECTION 2 — PER-LINE HISTOTYPE-CONSISTENCY CHECK
# =============================================================================
# Expression scoring -----------------------------------------------------------
# Per-line expression authentication in 2D CULTURE is intrinsically limited:
# cultured lines erode in-vivo markers (WT1/MUC16 down in HGS) and broadly
# express secretory/mesenchymal genes (TFF1/TFF3, VIM), so a naive "which subtype
# marker panel scores highest" (argmax) spuriously reclassifies bona fide HGS
# lines. We therefore DON'T argmax. Instead we score the two subtypes that have
# truly SPECIFIC, lineage-defining programs — clear cell (HNF1B/NAPSA) and
# mucinous/intestinal (CDX2/MUC5AC/MUC2/KRT20) — plus a serous/Mullerian
# retention score (PAX8/SOX17/WT1). A line is expression-authenticated POSITIVELY
# for CC/MC if it expresses that specific program; HGS (the serous default) is
# authenticated by ABSENCE of a competing specific program (+ its genomics);
# SCCOHT by SWI/SNF loss; EC has no specific positive panel (-> genomics/STR).
# Multivariate subtype SEPARATION (scripts 02/04) is the population-level
# expression validation; this per-line layer is supportive, not definitive.
logz <- zscore(logtpm)                                    # all symbols, z across lines
spec_panels <- list(CC = c("HNF1B","NAPSA"),
                    MC = c("CDX2","MUC5AC","MUC2","KRT20"))
spec_panels <- lapply(spec_panels, function(g) intersect(g, rownames(logz)))
ser_genes <- intersect(c("PAX8","SOX17","WT1"), rownames(logz))
prog_score <- function(cl, genes) if (cl %in% colnames(logz))
  mean(logz[genes, cl]) else NA_real_
STRONG <- 1.0     # z threshold for a "strong" specific program (competing signal)
PRESENT <- 0.5    # z threshold for a line to be expressing its own specific program

# TP53 / driver / CNV summaries from WES --------------------------------------
tp53 <- wm %>% filter(Hugo_Symbol == "TP53") %>%
  group_by(cell_line) %>% summarise(TP53 = paste0("mut(",
    paste(unique(gsub("^p\\.", "", HGVSp_Short)), collapse = ";"), ")"), .groups = "drop")
# Driver calls, SPLIT BY SOMATIC-CONFIDENCE TIER [review revision]. `key_drivers`
# previously pooled every retained driver call with no tier qualification and was
# carried into Table S1 with no tier column at all, so TOV81D read as
# key_drivers = "BRCA2" — one of exactly the two BRCA2 calls the manuscript states
# are Tier 3 and excluded ("defensible somatic BRCA1/2 is zero").
drv_struct <- wm %>% filter(is_driver) %>% distinct(cell_line, gene = Hugo_Symbol) %>%
  left_join(tier_min, by = c("cell_line", "gene"))
drivers <- drv_struct %>%
  group_by(cell_line) %>%
  summarise(key_drivers = paste(sort(unique(gene)), collapse = ", "),
            key_drivers_tier12 = paste(sort(unique(gene[!is.na(tier_best) & tier_best <= 2L])),
                                       collapse = ", "),
            key_drivers_tier3  = paste(sort(unique(gene[!is.na(tier_best) & tier_best == 3L])),
                                       collapse = ", "),
            key_drivers_annotated = paste(sort(unique(ifelse(
              !is.na(tier_best) & tier_best == 3L, paste0(gene, "(T3)"), gene))),
              collapse = ", "),
            .groups = "drop") %>%
  mutate(across(c(key_drivers_tier12, key_drivers_tier3), ~ifelse(.x == "", NA_character_, .x)))
# *** AUTOSOME-RESTRICTED FGA (C2 fix) ***
# This previously read `fga_0.2`, the LEGACY chrX-inclusive column. 08_wes_cnv.R was
# explicitly revised to make `fga_auto_0.2` the headline metric because chrX carries a
# pooled-normal sex-composition artifact, and the manuscript reports autosome-
# restricted FGA throughout — so the authentication table was scoring genome
# instability on a quantity the pipeline itself declares an artifact, and Table S1
# consequently showed OV90 with fga_autosome = 0.269 beside genomics_consistent =
# "consistent" under a rule stated as FGA > 0.30. CONSEQUENCE OF THE FIX: OV90 is the
# only line that crosses the 0.30 threshold (0.3084 with chrX vs 0.2689 autosome), so
# its genomics_consistent moves consistent -> partial and cnv_instability high ->
# intermediate. No other line changes category (asserted below).
fga_l <- fga %>% transmute(cell_line,
  FGA = round(fga_auto_0.2, 3),
  FGA_metric = "fga_auto_0.2 (autosome-restricted; chrX excluded as a pooled-normal artifact)",
  FGA_withX_legacy = round(fga_0.2, 3),
  cnv_instability = cut(fga_auto_0.2, c(-Inf, 0.15, 0.30, Inf),
                        labels = c("low","intermediate","high")),
  cnv_instability_withX_legacy = cut(fga_0.2, c(-Inf, 0.15, 0.30, Inf),
                        labels = c("low","intermediate","high")))
chg <- fga_l %>% filter(as.character(cnv_instability) != as.character(cnv_instability_withX_legacy))
cat("\n=== C2: lines whose CNV-instability category changes when chrX is excluded ===\n")
if (nrow(chg)) print(as.data.frame(chg %>% select(cell_line, FGA_withX_legacy, FGA,
  cnv_instability_withX_legacy, cnv_instability)), row.names = FALSE) else cat("  none\n")
stopifnot("more lines changed CNV-instability category than the expected single line (OV90) — investigate before accepting" =
            nrow(chg) <= 1)

# canonical-driver expectations per subtype (for genomics_consistent)
canon <- list(
  HGS  = "TP53",
  LGS  = "MAPK(KRAS/BRAF/NRAS) or quiet genome; TP53-wt",
  CC   = "ARID1A/PIK3CA/PTEN",
  EC   = "CTNNB1/PIK3CA/ARID1A/PTEN; TP53-wt (typical EC)",
  MC   = "KRAS (+/-CDKN2A/ERBB2); no APC",
  MMMT = "TP53 (carcinosarcoma)",
  SCCOHT = "SMARCA4 (usually protein-level)"
)

perline <- gen %>%
  left_join(tp53, by = "cell_line") %>%
  left_join(drivers, by = "cell_line") %>%
  left_join(fga_l, by = "cell_line") %>%
  left_join(swisnf_panel %>% select(cell_line, swisnf_deficient, swisnf_evidence,
                                    swisnf_tier3_only_calls,
                                    SMARCA4_loss, SMARCA2_loss, ARID1A_loss, SMARCB1_loss),
            by = "cell_line") %>%
  rowwise() %>%
  mutate(
    has_wes = cell_line %in% wes_lines,
    CC_prog_z     = round(prog_score(cell_line, spec_panels$CC), 2),
    MC_prog_z     = round(prog_score(cell_line, spec_panels$MC), 2),
    serous_prog_z = round(prog_score(cell_line, ser_genes), 2),
    TP53 = ifelse(has_wes & is.na(TP53), "no retained candidate", TP53)
  ) %>% ungroup()

# expression_consistent (consistent / partial / discordant / not assessed) ----
# "discordant" (a strong claim) is reserved for a POSITIVE competing specific
# program (e.g. an HGS/MMMT line strongly expressing the CC or MC/intestinal
# program). A line that merely fails to express its OWN eroded marker panel is
# "partial", not discordant (absence of a positive marker in 2D culture is not
# evidence of mislabel). EC has no specific panel -> "not assessed".
perline <- perline %>%
  rowwise() %>%
  mutate(
    expression_competing = {
      comp <- c()
      if (as.character(labeled_subtype) != "CC" && !is.na(CC_prog_z) && CC_prog_z > STRONG)
        comp <- c(comp, sprintf("clear-cell program (HNF1B/NAPSA z=%.1f)", CC_prog_z))
      if (as.character(labeled_subtype) != "MC" && !is.na(MC_prog_z) && MC_prog_z > STRONG)
        comp <- c(comp, sprintf("mucinous/intestinal program (CDX2/MUC5AC z=%.1f)", MC_prog_z))
      if (length(comp) == 0) NA_character_ else paste(comp, collapse = "; ")
    },
    expression_consistent = {
      lab <- as.character(labeled_subtype)
      if (!has_rna) "not assessed"
      else if (lab == "MC") {
        if (!is.na(MC_prog_z) && MC_prog_z > PRESENT) "consistent"
        else if (!is.na(MC_prog_z) && MC_prog_z > 0) "partial" else "partial"
      }
      else if (lab == "CC") {
        if (!is.na(CC_prog_z) && CC_prog_z > PRESENT) "consistent" else "partial"
      }
      else if (lab %in% c("HGS","MMMT")) {
        if (is.na(expression_competing)) "consistent" else "discordant"
      }
      else if (lab == "SCCOHT") {
        if (isTRUE(swisnf_deficient)) "consistent" else "partial"
      }
      else "not assessed"     # EC: no specific positive panel -> genomics/STR
    },
    # WHAT KIND OF EVIDENCE the call rests on [review revision]. "26 expression-
    # consistent" cannot be quoted without this: for HGS/MMMT "consistent" means only
    # the ABSENCE of a competing clear-cell or mucinous program (culture-eroded
    # markers make a positive serous call impossible), and for SCCOHT it reuses the
    # very SWI/SNF evidence that produced the reclassification, which is circular.
    # Only CC and MC are authenticated on a POSITIVE lineage program.
    expression_basis = {
      lab <- as.character(labeled_subtype)
      if (!has_rna) "not assessed (no RNA)"
      else if (lab %in% c("CC","MC")) "positive lineage program (specific markers expressed)"
      else if (lab %in% c("HGS","MMMT")) "absence of a competing CC/MC program (NOT a positive serous call)"
      else if (lab == "SCCOHT") "SWI/SNF loss — same evidence as the reclassification (circular)"
      else "not assessed (no specific positive panel for this subtype)"
    }
  ) %>% ungroup()

# genomics_consistent (Y/partial/N/not assessed) ------------------------------
# [review revision] Look the driver up in the STRUCTURED table rather than regex-
# matching a comma-joined string — the string form was a needless failure mode
# (a gene name that is a substring of another, a changed separator) when the tidy
# table is right here.
drv_by_line <- split(drv_struct$gene, drv_struct$cell_line)
has_drv <- function(cl, genes) length(intersect(genes, drv_by_line[[cl]])) > 0
perline <- perline %>%
  rowwise() %>%
  mutate(genomics_consistent = {
    lab <- as.character(labeled_subtype)
    if (!has_wes && is.na(FGA)) "not assessed"   # neither mutations nor CNV
    else if (lab == "HGS") {
      if (!is.na(TP53) && grepl("^mut", TP53) && !is.na(FGA) && FGA > 0.30) "consistent"
      else if (!is.na(TP53) && grepl("^mut", TP53)) "partial"
      else if (!has_wes && !is.na(FGA) && FGA > 0.30) "partial"  # CNV-only (e.g. TOV3121D, no MAF): high instability supports HGS
      else "discordant"
    }
    else if (!has_wes) "not assessed"   # non-HGS lines need mutation calls (driver logic below)
    else if (lab == "CC") {
      if (has_drv(cell_line, c("ARID1A","PIK3CA","PTEN"))) "consistent"
      else if (has_drv(cell_line, c("KRAS"))) "partial" else "discordant"
    }
    else if (lab == "EC") {
      # typical EC = TP53-wt + CTNNB1/PIK3CA/ARID1A; TP53-mut + SMARCA4 loss = dedifferentiated
      if (!is.na(TP53) && grepl("^mut", TP53) && isTRUE(swisnf_deficient)) "discordant"
      else if (has_drv(cell_line, c("CTNNB1","PIK3CA","ARID1A","PTEN"))) "consistent"
      else "partial"
    }
    else if (lab == "MC") {
      if (has_drv(cell_line, c("KRAS")) && !has_drv(cell_line, c("APC"))) "consistent"
      else if (has_drv(cell_line, c("APC"))) "discordant" else "partial"
    }
    else if (lab == "LGS") {
      # quiet genome + TP53-wt is the LGS signature; MAPK not always recovered
      if ((is.na(TP53) || TP53 == "no retained candidate") && !is.na(FGA) && FGA < 0.15) "consistent"
      else if (!is.na(TP53) && grepl("^mut", TP53)) "discordant" else "partial"
    }
    else "not assessed"
  }) %>% ungroup()

# flags + provisional molecular call (curated from the evidence above) --------
perline <- perline %>%
  rowwise() %>%
  mutate(
    flags = {
      lab <- as.character(labeled_subtype)
      f <- c()
      if (isTRUE(grepl("CONFLICT", subtype_status)))
        f <- c(f, "label conflict in source metadata")
      # --- independently recovered reclassifications (Karnezis 2021) ---
      if (lab == "EC" && !is.na(TP53) && grepl("^mut", TP53) && isTRUE(swisnf_deficient))
        f <- c(f, "TP53-mut + SWI/SNF (SMARCA4/SMARCA2) loss -> dedifferentiated carcinoma, NOT typical EC (Karnezis 2021)")
      if (lab == "SCCOHT" && isTRUE(swisnf_deficient))
        f <- c(f, "SWI/SNF-deficient (SMARCA4 protein/WES +/- SMARCA2 mRNA) supports SCCOHT (Karnezis 2021)")
      if (cell_line %in% c("OV90","OV3331") && !is.na(TP53) && grepl("^mut", TP53))
        f <- c(f, "Adeno-vs-HGS label conflict: TP53-mut + high CNV support HGS-family")
      if (lab == "HGS" && !is.na(serous_prog_z) && serous_prog_z < -1)
        f <- c(f, paste0("serous/Mullerian markers very low (PAX8/WT1/SOX17 z=", serous_prog_z,
                        ") — serous identity NOT confirmed by expression (culture erosion vs atypia); HGS call rests on genomics"))
      # --- expression competing program (the genuine expression discordance) ---
      if (!is.na(expression_competing))
        f <- c(f, paste0("expression resembles ", expression_competing,
                        " (review vs labeled ", lab, ")"))
      # --- SWI/SNF signals interpreted in subtype context ---
      if (isTRUE(ARID1A_loss) && lab == "CC")
        f <- c(f, "ARID1A truncating mutation = canonical clear-cell driver (supportive)")
      if (isTRUE(swisnf_deficient) && lab == "CC" && !isTRUE(ARID1A_loss))
        f <- c(f, paste0("atypical SWI/SNF signal for CC: ", swisnf_evidence, " — verify"))
      if (isTRUE(swisnf_deficient) && lab == "HGS")
        f <- c(f, paste0("secondary SWI/SNF subunit variant: ", swisnf_evidence))
      # --- mucinous origin caveat (script 11 does the full work-up) ---
      if (cell_line == "TOV2414")
        f <- c(f, "ovarian mucinous — published derivation (Sauriol 2020: KRAS G12A); see script 11")
      if (lab == "MC" && cell_line != "TOV2414")
        f <- c(f, "ovarian-vs-GI origin unresolved (no external provenance) — see script 11 + STR/IHC")
      # --- EC after reclassification ---
      if (cell_line == "VOA4395")
        f <- c(f, "becomes sole EC (n=1) if TOV112D reclassified; no WES — STR pending")
      if (length(f) == 0) "none" else paste(f, collapse = "; ")
    }) %>% ungroup()

perline <- perline %>%
  mutate(provisional_molecular_call = dplyr::case_when(
    cell_line == "TOV112D" ~ "dedifferentiated carcinoma (candidate SMARCA4 loss; historical label 'EC')",
    cell_line == "COV434"  ~ "SCCOHT (candidate SMARCA4/SMARCA2 loss; historical label 'granulosa')",
    cell_line == "BIN67"   ~ "SCCOHT (low SMARCA4 protein + low SMARCA2 RNA)",
    cell_line %in% c("OV90","OV3331") & serous_prog_z < -1 ~
      "HGS-family carcinoma; serous identity NOT confirmed by expression (Adeno label unresolved)",
    cell_line %in% c("OV90","OV3331") ~ "HGS (supported over 'Adenocarcinoma')",
    genomics_consistent == "consistent" & expression_consistent %in% c("consistent","partial") ~
      paste0(as.character(labeled_subtype), " (concordant)"),
    genomics_consistent == "consistent" ~ paste0(as.character(labeled_subtype), " (genomics-supported)"),
    expression_consistent == "consistent" ~ paste0(as.character(labeled_subtype), " (expression-supported)"),
    genomics_consistent == "partial" ~ paste0(as.character(labeled_subtype), " (partial genomic support)"),
    expression_consistent == "partial" ~ paste0(as.character(labeled_subtype), " (partial expression support)"),
    TRUE ~ paste0(as.character(labeled_subtype), " (as labeled; STR pending)")))

# --- STR status, populated rather than left as an empty placeholder ----------
# This column previously shipped as all-NA, which reads as "not checked" when in
# fact the public reference status IS known for every line. Script 18 parses the
# Cellosaurus records; run_all.sh therefore runs 18 before 10 (neither reads the
# other's outputs, so the ordering is free). In-house profiles generated for this
# resource are appended by the same join once str_inhouse.csv exists.
str_f <- file.path(OUT, "cellosaurus_str_status.csv")
stopifnot("output/cellosaurus_str_status.csv missing — run 18_external_benchmarking.R first" =
            file.exists(str_f))
str_pub <- readr::read_csv(str_f, show_col_types = FALSE) %>%
  transmute(cell_line, in_cellosaurus, cellosaurus_accession,
            str_profile_documented, n_str_markers)

inh_f <- file.path(OUT, "str_inhouse.csv")
str_inh <- if (file.exists(inh_f)) {
  readr::read_csv(inh_f, show_col_types = FALSE) %>%
    transmute(cell_line, inhouse_str_result)          # e.g. "match to reference (16/16)"
} else {
  tibble(cell_line = character(), inhouse_str_result = character())
}

perline <- perline %>%
  left_join(str_pub, by = "cell_line") %>%
  left_join(str_inh, by = "cell_line") %>%
  mutate(STR_status = dplyr::case_when(
    !is.na(inhouse_str_result) & isTRUE(str_profile_documented) ~
      paste0("in-house profile, ", inhouse_str_result,
             "; reference profile in Cellosaurus ", cellosaurus_accession,
             " (", n_str_markers, " markers)"),
    !is.na(inhouse_str_result) ~
      paste0("in-house profile, ", inhouse_str_result, "; no Cellosaurus reference profile"),
    str_profile_documented ~
      paste0("reference profile in Cellosaurus ", cellosaurus_accession,
             " (", n_str_markers, " markers); current-stock STR documentation unconfirmed"),
    in_cellosaurus ~
      paste0("Cellosaurus record ", cellosaurus_accession,
             " carries no STR profile; current-stock STR documentation unconfirmed"),
    TRUE ~ "no Cellosaurus record and current-stock STR documentation unconfirmed")) %>%
  select(-in_cellosaurus, -cellosaurus_accession, -str_profile_documented,
         -n_str_markers, -inhouse_str_result)

stopifnot("STR_status must be populated for every line" = !any(is.na(perline$STR_status)))
cat("\n=== STR status ===\n"); print(table(sub(";.*", "", perline$STR_status)))

auth_table <- perline %>%
  select(cell_line, labeled_subtype, source_site, has_rna, has_prot, has_wes,
         CC_prog_z, MC_prog_z, serous_prog_z,
         expression_consistent, expression_basis, expression_competing,
         TP53, key_drivers, key_drivers_tier12, key_drivers_tier3,
         key_drivers_annotated,
         FGA, FGA_metric, FGA_withX_legacy, cnv_instability, genomics_consistent,
         swisnf_deficient, swisnf_tier3_only_calls,
         flags, provisional_molecular_call, STR_status) %>%
  arrange(labeled_subtype, cell_line)
readr::write_csv(auth_table, file.path(OUT, "auth_perline_table.csv"))

cat("\n=== Per-line authentication: consistency tally ===\n")
print(auth_table %>% count(expression_consistent, name = "n_expr"))
print(auth_table %>% count(genomics_consistent, name = "n_geno"))

# THE BREAKDOWN THAT MUST TRAVEL WITH THE HEADLINE COUNT [review revision]
cat("\n=== 'expression_consistent' BY EVIDENCE TYPE — quote this, not the bare count ===\n")
print(as.data.frame(auth_table %>% filter(expression_consistent == "consistent") %>%
  count(labeled_subtype, expression_basis, name = "n")), row.names = FALSE)
eb <- auth_table %>% filter(expression_consistent == "consistent") %>% count(expression_basis)
cat(sprintf("  total 'consistent' = %d, of which:\n%s\n",
            sum(eb$n), paste(sprintf("    %d  %s", eb$n, eb$expression_basis), collapse = "\n")))
cat("  => only the 'positive lineage program' rows are positive authentication; the\n",
    "     HGS/MMMT rows are absence-of-competing-program and the SCCOHT rows reuse the\n",
    "     reclassification evidence. Do not report a single number without this split.\n")
cat("\n=== Tier-3-only driver calls that must NOT be quoted as somatic ===\n")
print(as.data.frame(auth_table %>% filter(!is.na(key_drivers_tier3)) %>%
  select(cell_line, labeled_subtype, key_drivers_tier12, key_drivers_tier3)), row.names = FALSE)
cat("\n=== Lines with flags (discordances / conflicts surfaced) ===\n")
print(as.data.frame(auth_table %>% filter(flags != "none") %>%
  select(cell_line, labeled_subtype, expression_consistent, genomics_consistent,
         provisional_molecular_call, flags)), row.names = FALSE)

# -----------------------------------------------------------------------------
# 2b. Histotype-consistency figure (status tiles + CNV instability bar)
# -----------------------------------------------------------------------------
# All four tracks are encoded as CONSISTENCY WITH THE LABELED SUBTYPE'S
# EXPECTATION (not raw status), so the semantics are uniform: blue = as expected,
# orange = contradicts the label (reclassification signal), gold = weak/atypical,
# grey = not assessable. Under this framing TP53-mutated reads as EXPECTED (blue)
# for HGS but as a NOT-typical-EC signal (orange) for TOV112D; TP53-wildtype reads
# as EXPECTED (blue) for LGS/CC/MC; SWI/SNF-deficiency reads as EXPECTED (blue) for
# SCCOHT but as a reclassification signal (orange) for the EC label.
expect_tp53 <- function(st, tp53, has_wes) {
  if (!isTRUE(has_wes) || is.na(tp53)) return("not assessed")
  mut <- grepl("^mut", tp53)
  if (st %in% c("HGS","MMMT")) if (mut) "consistent" else "discordant"
  else if (st == "EC")         if (mut) "discordant" else "consistent"
  else if (st %in% c("LGS","CC","MC")) if (mut) "partial" else "consistent"
  else "not assessed"
}
expect_swisnf <- function(st, def, arid, assessable) {
  if (!isTRUE(assessable) || is.na(def)) return("not assessed")
  if (st == "SCCOHT") { if (def) "consistent" else "partial" }
  else if (st == "EC") { if (def) "discordant" else "consistent" }
  else if (st == "CC") { if (!def) "consistent" else if (isTRUE(arid)) "consistent" else "partial" }
  else { if (def) "partial" else "consistent" }        # HGS / MC / LGS / MMMT
}
tile_df <- perline %>%
  rowwise() %>%
  transmute(cell_line, labeled_subtype,
            Expression = expression_consistent,
            Genomics   = genomics_consistent,
            `TP53`     = expect_tp53(as.character(labeled_subtype), TP53, has_wes),
            `SWI/SNF`  = expect_swisnf(as.character(labeled_subtype), swisnf_deficient,
                                       ARID1A_loss, has_rna || has_prot)) %>%
  ungroup() %>%
  pivot_longer(c(Expression, Genomics, `TP53`, `SWI/SNF`),
               names_to = "metric", values_to = "status") %>%
  mutate(metric = factor(metric, levels = c("Expression","Genomics","TP53","SWI/SNF")),
         status = factor(status, levels = names(STATUS_COLS)))
lo <- auth_table %>% arrange(labeled_subtype, cell_line) %>% pull(cell_line)
tile_df$cell_line <- factor(tile_df$cell_line, levels = rev(lo))

p_tile <- ggplot(tile_df, aes(metric, cell_line, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = STATUS_COLS, name = NULL, drop = FALSE,
                    labels = c(consistent = "consistent w/ label",
                               discordant = "contradicts label (reclass. signal)",
                               partial = "weak / atypical",
                               "not assessed" = "not assessable (data / STR pending)")) +
  facet_grid(labeled_subtype ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Per-line histotype consistency with labeled subtype (OvCAN generated lines)",
       subtitle = paste0("Each cell = agreement with what the LABELED subtype predicts.  ",
                        "TOV112D & COV434 are consistent with Karnezis-2021 reclassifications;\n",
                        "OV90/OV3331 Adeno-vs-HGS conflict resolves to HGS.  ",
                        "STR profiles pending (BC Cancer/OVCARE)."),
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        panel.grid = element_blank(), legend.position = "top",
        legend.text = element_text(size = 7.5))
ggsave(file.path(FIGS, "10_histotype_consistency.pdf"), p_tile, width = 8, height = 9)

cat("\nOutputs: output/auth_swisnf_panel.csv, output/auth_perline_table.csv\n")
cat("Figures: figs/10_swisnf_panel.pdf, figs/10_histotype_consistency.pdf\n")
message("\n10_authentication.R complete.")
