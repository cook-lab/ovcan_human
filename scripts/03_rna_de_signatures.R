# =============================================================================
# Script: 03_rna_de_signatures.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: One-vs-rest differential expression per subtype, up-signatures, and
#          GO-BP pathway recovery (fgsea) — the "recapitulates known biology"
#          evidence for the Data Descriptor.
#          All contrasts describe this selected panel. Sample size alone does
#          not establish generalisable histotype signatures, particularly for
#          the rare groups (EC/MMMT/SCCOHT n=2, MC n=3).
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   1 (RNA re-processing) — step 3 of 4 (this = DE/signatures/GO)
# =============================================================================
# CONFOUNDER NOTE: one-versus-rest contrasts describe this panel. All HGS models
# share one source centre, and related sublines are not independent replicates.
# The additive subtype+site design is full-rank because CC spans centres, but
# missing histotype-centre combinations prevent validation of that additive
# assumption. Marker/pathway agreement does not eliminate source-aligned effects.
# Patient-representative comparisons are provided by script 21; no sample-count
# threshold turns these selected models into generalisable histotype signatures.
#
# BUG FIX (archived RNA notebook, lines 609/703): CC and MC signatures were
#   built by reading the DE CSV WITHOUT row.names=1, then `Gene <- rownames(df)`
#   overwrote symbols with "1","2",... (integer strings), which then failed the
#   downstream `Gene %in% rownames(tpm)` filter and silently emptied those two
#   signatures. Here we never round-trip through a CSV: we map Ensembl gene_id
#   -> symbol via the pinned tx2gene map for ALL subtypes and assert the
#   signature symbols are real (Section 2/2b).
# =============================================================================

source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(DESeq2); library(fgsea)
  library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi); library(ggrepel)
})
set.seed(SEED)

# -----------------------------------------------------------------------------
# 1. Inputs: DESeq2 object (counts + length offset from tximport) + gene map
# -----------------------------------------------------------------------------
dds <- readRDS(file.path(OUT, "rna_dds.rds"))
dds$subtype <- factor(dds$subtype, levels = c("HGS","CC","EC","MC","MMMT","SCCOHT"))
subtypes <- levels(dds$subtype)
n_by     <- table(dds$subtype)
message("Subtype n (RNA analysis set): ",
        paste(names(n_by), as.integer(n_by), sep = "=", collapse = ", "))

# pinned Ensembl gene_id -> symbol (distinct; the fix relies on this map)
t2g <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"),
                       show_col_types = FALSE)
gene_meta <- t2g %>% distinct(ensembl_gene_id, external_gene_name, gene_biotype)
id2sym <- gene_meta %>% distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  { setNames(.$external_gene_name, .$ensembl_gene_id) }

# All one-versus-rest model comparisons describe the sampled panel.
formal_ok <- character() # All model-level contrasts are descriptive.
message("Descriptive panel contrasts: ",
        paste(setdiff(subtypes, formal_ok), collapse = ", "))

# -----------------------------------------------------------------------------
# 2. One-vs-rest DE per subtype (binary target-vs-Other refit)
#    Reuses the tximport normalizationFactors (length offset) stored in dds;
#    changing the design does not clear them, so the counts path stays correct.
# -----------------------------------------------------------------------------
run_ovr <- function(st) {
  d <- dds
  d$grp <- factor(ifelse(d$subtype == st, st, "Other"), levels = c("Other", st))
  design(d) <- ~ grp
  d <- DESeq(d, quiet = TRUE)
  res <- results(d, contrast = c("grp", st, "Other"), alpha = 0.05)
  as_tibble(as.data.frame(res), rownames = "gene_id") %>%
    mutate(symbol  = unname(id2sym[gene_id]),
           subtype = st,
           n_group = as.integer(n_by[st]),
           evidence = "descriptive; model-level, source-confounded") %>%
    arrange(padj)
}

de_list <- setNames(lapply(subtypes, run_ovr), subtypes)

# per-subtype DE tables + one combined table
for (st in subtypes)
  readr::write_csv(de_list[[st]], file.path(OUT, sprintf("rna_de_%s.csv", st)))
de_all <- bind_rows(de_list)
readr::write_csv(de_all, file.path(OUT, "rna_de_all.csv"))

# -----------------------------------------------------------------------------
# 2b. Up-signatures (padj<0.05 & log2FC>1), symbol-mapped for ALL subtypes.
#     Keep only real symbols; collapse duplicate symbols to the best padj.
# -----------------------------------------------------------------------------
make_sig <- function(tab) {
  tab %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange > 1,
           !is.na(symbol), symbol != "") %>%
    arrange(padj) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    transmute(subtype, symbol, gene_id, log2FoldChange, padj,
              n_group, evidence)
}
sig_list <- lapply(de_list, make_sig)
for (st in subtypes)
  readr::write_csv(sig_list[[st]], file.path(OUT, sprintf("rna_signatures_%s.csv", st)))
sig_all <- bind_rows(sig_list)
readr::write_csv(sig_all, file.path(OUT, "rna_signatures_all.csv"))

# --- ASSERT the bug is fixed: signature symbols are real gene names, not indices
sig_sizes <- sig_all %>% dplyr::count(subtype, name = "n_sig") %>%
  right_join(tibble(subtype = subtypes), by = "subtype") %>%
  mutate(n_sig = replace_na(n_sig, 0L),
         subtype = factor(subtype, levels = subtypes)) %>%
  left_join(tibble(subtype = factor(names(n_by), levels = subtypes),
                   n_lines = as.integer(n_by)), by = "subtype") %>%
  arrange(subtype)
bad <- sig_all %>% filter(grepl("^[0-9]+$", symbol))
stopifnot("BUG NOT FIXED: integer-string symbols present in signatures" = nrow(bad) == 0)
stopifnot("CC signature empty — the archived Gene bug may have recurred" =
            sig_sizes$n_sig[sig_sizes$subtype == "CC"] > 0)
stopifnot("MC signature empty — the archived Gene bug may have recurred" =
            sig_sizes$n_sig[sig_sizes$subtype == "MC"] > 0)
cat("\n=== Up-signature sizes (padj<0.05 & log2FC>1; real symbols) ===\n")
print(as.data.frame(sig_sizes))
cat("Bug-fix check: CC & MC signatures non-empty with real symbols. Examples —\n")
cat("  CC:", paste(head(sig_list[["CC"]]$symbol, 8), collapse = ", "), "\n")
cat("  MC:", paste(head(sig_list[["MC"]]$symbol, 8), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 3. Figure: descriptive HGS volcano (BH-adjusted values on y-axis)
# -----------------------------------------------------------------------------
hgs <- de_list[["HGS"]] %>%
  filter(!is.na(padj)) %>%
  mutate(padj_plot = pmax(padj, min(padj[padj > 0], na.rm = TRUE) / 10),
         sig = case_when(padj < 0.05 & log2FoldChange >  1 ~ "up in HGS",
                         padj < 0.05 & log2FoldChange < -1 ~ "down in HGS",
                         TRUE ~ "ns"))
lab <- hgs %>% filter(sig != "ns", symbol != "", !is.na(symbol)) %>%
  group_by(sig) %>% slice_min(padj, n = 12) %>% ungroup()
volc_cols <- c("up in HGS" = COOK_RUST, "down in HGS" = COOK_NAVY, "ns" = "grey80")
p_volc <- ggplot(hgs, aes(log2FoldChange, -log10(padj_plot), colour = sig)) +
  geom_point(size = 1.1, alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") +
  ggrepel::geom_text_repel(data = lab, aes(label = symbol), size = 2.4,
                           max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = volc_cols) +
  labs(x = "log2 fold-change (HGS vs rest)", y = "-log10 adjusted p",
       colour = NULL, title = "HGS one-vs-rest differential expression",
       subtitle = sprintf("Descriptive panel contrast; related models and source confounding limit inference.")) +
  theme_minimal(base_size = 11)
ggsave(file.path(FIGS, "03_rna_hgs_volcano.pdf"), p_volc, width = 7.5, height = 6)

# -----------------------------------------------------------------------------
# 4. Figure: signature sizes by subtype (n annotated; evidence flagged)
# -----------------------------------------------------------------------------
sz <- sig_sizes %>%
  mutate(evidence = "descriptive panel comparison",
         lab = sprintf("%d lines", n_lines))
p_sz <- ggplot(sz, aes(subtype, n_sig, fill = evidence)) +
  geom_col() +
  geom_text(aes(label = lab), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("descriptive panel comparison" = COOK_NAVY)) +
  labs(x = NULL, y = "up-signature genes (padj<0.05 & log2FC>1)", fill = NULL,
       title = "Subtype up-signature sizes",
       subtitle = "All contrasts are descriptive; patient-representative results are in script 21") +
  theme_minimal(base_size = 11)
ggsave(file.path(FIGS, "03_rna_signature_sizes.pdf"), p_sz, width = 7, height = 5)

# -----------------------------------------------------------------------------
# 5. GO-BP recovery via fgsea (threshold-free; ranked by Wald stat)
#    fgsea on the full ranked list is more robust than ORA on a tiny n<=3
#    signature, and does not depend on the arbitrary signature cutoff.
# -----------------------------------------------------------------------------
message("Building GO-BP gene sets from pinned org.Hs.eg.db ...")
go2eg  <- as.list(org.Hs.egGO2ALLEGS)                    # Entrez per GO (w/ descendants)
onto   <- AnnotationDbi::Ontology(names(go2eg))
bp_ids <- names(go2eg)[!is.na(onto) & onto == "BP"]
eg2sym <- unlist(as.list(org.Hs.egSYMBOL))               # Entrez -> symbol
go_bp  <- lapply(go2eg[bp_ids], function(e)
  unique(na.omit(unname(eg2sym[unique(e)]))))
keep   <- lengths(go_bp) >= 10 & lengths(go_bp) <= 500    # drop tiny/huge terms
go_bp  <- go_bp[keep]
names(go_bp) <- paste0(names(go_bp), " ", AnnotationDbi::Term(names(go_bp)))
message(sprintf("GO-BP sets used (size 10-500): %d", length(go_bp)))

rank_syms <- function(tab) {
  v <- tab %>% filter(!is.na(stat), !is.na(symbol), symbol != "") %>%
    group_by(symbol) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  setNames(v$stat, v$symbol)
}
# fgsea uses Monte Carlo sampling. Set the seed within each call so results are
# reproducible independently of call order, and use 50,000 initial permutations.
# Section 5b records the actual adjusted-p ranges across three seeds; precision
# checks do not remove patient dependence or source confounding.
GSEA_NPERM_SIMPLE <- 50000L

run_gsea <- function(st, seed = SEED) {
  set.seed(seed)
  ranks <- rank_syms(de_list[[st]])
  fg <- suppressWarnings(fgsea(pathways = go_bp, stats = ranks,
                               minSize = 10, maxSize = 500, eps = 0,
                               nPermSimple = GSEA_NPERM_SIMPLE))
  fg %>% as_tibble() %>%
    mutate(subtype = st, n_group = as.integer(n_by[st]),
           evidence = "descriptive; model-level, source-confounded",
           go_id = sub(" .*", "", pathway),
           term  = sub("^GO:[0-9]+ ", "", pathway),
           leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) %>%
    arrange(padj)
}
gsea_all <- bind_rows(lapply(subtypes, run_gsea))
readr::write_csv(gsea_all, file.path(OUT, "rna_de_gsea_go.csv"))

# --- did the EXPECTED known biology come back? (targeted, honest probe) -------
#   One-vs-rest among (uniformly proliferating) cancer cell lines is a stringent
#   contrast: shared programs (proliferation) do not differentiate, and 6369 GO
#   terms carry a heavy multiple-testing penalty. We therefore report the BEST
#   positive-NES hit per expected program with BOTH nominal p and padj, and grade
#   recovery: recovered (padj<0.05) / suggestive (nominal p<0.05) / not recovered.
expected <- tribble(
  ~st,      ~label,                       ~pattern,
  "HGS",    "DNA repair (HR/DSB)",        "double-strand break repair|homologous recombination|mismatch repair|DNA repair",
  "CC",     "oxidative/glutathione/detox","glutathione|xenobiotic|sulfur amino acid|reactive oxygen|oxidative stress",
  "MC",     "glycan/oligosaccharide",     "oligosaccharide|glycosylation|O-glycan|mucin|glycoprotein metabolic",
  "MMMT",   "EMT / migration / ECM",      "epithelial to mesenchymal|extracellular matrix|cell migration|cell motility",
  "SCCOHT", "cell cycle / mitosis",       "sister chromatid|mitotic|centromere|replication fork|chromosome segregation"
)
grade_recovery <- function(gsea_tbl) {
  purrr::map_dfr(seq_len(nrow(expected)), function(i) {
    sub_hits <- gsea_tbl %>%
      filter(subtype == expected$st[i], grepl(expected$pattern[i], term, ignore.case = TRUE))
    pos <- sub_hits %>% filter(NES > 0) %>% arrange(pval)
    best <- if (nrow(pos) > 0) pos[1, ] else sub_hits %>% arrange(pval) %>% slice_head(n = 1)
    n_sig <- sum(pos$padj < 0.05, na.rm = TRUE)
    status <- if (nrow(pos) == 0) "not recovered (no positive-NES term)"
              else if (isTRUE(best$padj < 0.05)) "recovered (padj<0.05)"
              else if (isTRUE(best$pval < 0.05)) "suggestive (nominal p<0.05)"
              else "not recovered"
    tibble(subtype = expected$st[i], program = expected$label[i],
           n_pos_sig = n_sig, best_term = best$term,
           NES = round(best$NES, 2), pval = signif(best$pval, 3),
           padj = signif(best$padj, 3), status = status,
           n_group = as.integer(n_by[expected$st[i]]))
  })
}
recovery <- grade_recovery(gsea_all)
readr::write_csv(recovery, file.path(OUT, "rna_de_gsea_recovery.csv"))
cat("\n=== GO recovery of EXPECTED subtype biology (best positive-NES hit) ===\n")
for (i in seq_len(nrow(recovery)))
  cat(sprintf("[%s n=%d] %-28s -> %s\n     best: %-46s NES=%+.2f p=%.1e padj=%.1e (%d terms padj<0.05)\n",
              recovery$subtype[i], recovery$n_group[i], recovery$program[i], recovery$status[i],
              substr(recovery$best_term[i], 1, 46), recovery$NES[i],
              recovery$pval[i], recovery$padj[i], recovery$n_pos_sig[i]))

# --- 5b. Seed stability of the recovery GRADES --------------------------------
#     The grade is a threshold call on a Monte-Carlo p-value, so it must be shown
#     not to move with the RNG. Re-run the 5 graded subtypes at two further seeds
#     and record the padj range and whether the grade ever changes. If `stable` is
#     FALSE for a program, that program's grade must NOT be reported as a result.
GSEA_SEEDS <- c(SEED, SEED + 1L, SEED + 2L)
graded_sts <- unique(expected$st)
stab <- purrr::map_dfr(GSEA_SEEDS, function(sd) {
  tb <- if (identical(sd, SEED)) gsea_all
        else bind_rows(lapply(graded_sts, run_gsea, seed = sd))
  grade_recovery(tb) %>% mutate(seed = sd)
})
stab_summary <- stab %>% group_by(subtype, program) %>%
  summarise(n_seeds = dplyr::n(),
            padj_min = min(padj), padj_max = max(padj),
            pval_min = min(pval), pval_max = max(pval),
            n_distinct_status = dplyr::n_distinct(status),
            statuses = paste(sort(unique(status)), collapse = " | "),
            stable = dplyr::n_distinct(status) == 1L, .groups = "drop")
readr::write_csv(stab, file.path(OUT, "rna_de_gsea_recovery_seeds.csv"))
readr::write_csv(stab_summary, file.path(OUT, "rna_de_gsea_recovery_stability.csv"))
cat(sprintf("\n--- Recovery-grade stability across %d fgsea seeds (nPermSimple=%d) ---\n",
            length(GSEA_SEEDS), GSEA_NPERM_SIMPLE))
print(as.data.frame(stab_summary %>% dplyr::select(subtype, program, padj_min,
                                                   padj_max, stable, statuses)),
      row.names = FALSE)
if (any(!stab_summary$stable))
  cat("!! FLAG: the grade for the program(s) above is RNG-dependent even at this precision;\n",
      "   report the padj range, not the grade.\n")

# -----------------------------------------------------------------------------
# 6. Figure: top recovered GO-BP terms per subtype (positive NES)
# -----------------------------------------------------------------------------
# unique per-row labels so ggplot orders terms independently within each facet
# (avoids a tidytext::reorder_within dependency); the zero-width marker after
# the term is stripped from the axis labels below.
top_go <- gsea_all %>%
  filter(NES > 0, padj < 0.05) %>%
  group_by(subtype) %>% slice_min(padj, n = 8, with_ties = FALSE) %>% ungroup() %>%
  mutate(subtype = factor(subtype, levels = subtypes),
         term_short = ifelse(nchar(term) > 45, paste0(substr(term, 1, 44), "..."), term),
         row = row_number(),
         lab = paste0(term_short, "@@@", row))
p_go <- ggplot(top_go, aes(NES, reorder(lab, NES))) +
  geom_segment(aes(x = 0, xend = NES, yend = reorder(lab, NES)), colour = "grey80") +
  geom_point(aes(colour = -log10(padj), size = size)) +
  scale_y_discrete(labels = function(x) sub("@@@.*", "", x)) +
  scale_colour_viridis_c(option = "mako", direction = -1, end = 0.9) +
  facet_wrap(~ subtype, scales = "free_y", ncol = 2) +
  labs(x = "NES (enriched in subtype)", y = NULL,
       colour = "-log10 padj", size = "set size",
       title = "GO-BP recovery per subtype (fgsea, top 8 by padj)",
       subtitle = "Descriptive model-level contrasts; related models and source confounding limit inference") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(FIGS, "03_rna_gsea_go_dotplot.pdf"), p_go, width = 11, height = 8)

# -----------------------------------------------------------------------------
# 7. Console report
# -----------------------------------------------------------------------------
cat("\n=== DE summary (one-vs-rest; genes with padj<0.05) ===\n")
print(de_all %>% group_by(subtype) %>%
        summarise(n_group = dplyr::first(n_group),
                  n_padj05 = sum(padj < 0.05, na.rm = TRUE),
                  n_up_sig = sum(padj < 0.05 & log2FoldChange > 1, na.rm = TRUE),
                  evidence = dplyr::first(evidence), .groups = "drop") %>%
        mutate(subtype = factor(subtype, levels = subtypes)) %>% arrange(subtype) %>%
        as.data.frame())
cat("\nPSEUDOREPLICATION NOTE: this contrast is PER LINE. The 15 HGS lines come from 12\n")
cat("  patients (1369, 2295, 3133 each contribute two models), so these observations are dependent.\n")
cat("  The model-level Wald/BH statistics do not account for that dependence.\n")
cat("  scripts/21_rna_sensitivity.R re-runs the contrasts and this\n")
cat("  GO recovery on the 28 patient representatives; report the two together.\n")

cat("\nOutputs: output/rna_de_*.csv, rna_signatures_*.csv, rna_de_gsea_go.csv\n")
cat("Figures: figs/03_rna_hgs_volcano.pdf, 03_rna_signature_sizes.pdf, 03_rna_gsea_go_dotplot.pdf\n")

# -----------------------------------------------------------------------------
# 8. Environment record
# -----------------------------------------------------------------------------
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
write_session_info("03_rna_de_signatures")

message("\n03_rna_de_signatures.R complete.")
