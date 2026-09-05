# =============================================================================
# Script: 21_rna_sensitivity.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Revision workstream — PATIENT-REPRESENTATIVE sensitivity analysis for
#          the expression arm. Patient-level collapse is applied to the genomic
#          frequencies (07/08) and to the marker effect sizes (17D), but NOT to
#          the PCA / t-SNE / silhouette (02), the one-vs-rest DE and GSEA (03),
#          the Hallmark scoring (04), the RNA-protein concordance (12) or the ADC
#          means (13) — and that switch of analysis unit was never flagged. This
#          script re-runs the load-bearing parts on ONE LINE PER PATIENT so the
#          two can be reported side by side.
#          Structure: 31 RNA lines from 28 patients (3 patients contribute 2 lines
#          each: families 1369, 2295, 3133). Within HGS it is 15 lines from 12
#          patients, i.e. 3 of the 15 DESeq2 "replicates" are additional models from already-represented patients,
#          which makes the HGS Wald/BH statistics anti-conservative by an amount
#          nobody had quantified. That is the open question here.
#          Deliverables:
#            A. PCA + silhouette + PC confounder + commonality on the 28 reps
#            B. One-vs-rest DESeq2 (all 6 subtypes) + GO recovery on the 28 reps
#            C. Side-by-side comparison tables
# Author:  Cook Lab (analyst: Claude — rna-variance workstream)
# Date:    2026-07-24
# Phase:   Revision (peer-review response: pseudoreplication in the expression arm)
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES NOT DO: it does not claim the 28-rep numbers are the
# "correct" ones and the 31-line numbers wrong. Both are legitimate — 31 lines is
# the resource as deposited, 28 patients is the unit for any inference about
# patients. The defect being fixed is that only one was reported, without saying
# which. Report both.
# =============================================================================
source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(SummarizedExperiment); library(DESeq2)
  library(matrixStats); library(cluster)
  library(fgsea); library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
set.seed(SEED)

SUBTYPE_LEVELS <- c("HGS", "CC", "EC", "MC", "MMMT", "SCCOHT")
N_HVG          <- 2000L    # matches script 02
N_PC_SIL       <- 10L      # matches script 02 (silhouette on PC1-10)
GSEA_NPERM_SIMPLE <- 50000L  # matches script 03 (see the precision note there)

# =============================================================================
# LOAD
# =============================================================================
vsd <- readRDS(file.path(OUT, "rna_vst.rds"))
v   <- assay(vsd)
dds <- readRDS(file.path(OUT, "rna_dds.rds"))
ann <- as.data.frame(colData(dds))
ann$subtype <- factor(ann$subtype, levels = SUBTYPE_LEVELS)
stopifnot(identical(colnames(v), ann$cell_line))

ensure_family_map()   # [integration revision] see 00_setup.R; same gap as 17.
fam <- readr::read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)
rna_fam <- tibble(cell_line = ann$cell_line) %>%
  left_join(fam %>% select(cell_line, patient_id, patient_representative), by = "cell_line")
stopifnot("family map does not cover every RNA line" = !anyNA(rna_fam$patient_id))

# One line per patient. Same deterministic rule as 17D: keep the flagged
# patient_representative among a patient's RNA lines; if none of them is flagged,
# keep the first (which is deterministic because colData order is fixed).
rep_keep <- rna_fam %>% group_by(patient_id) %>%
  summarise(cell_line = if (any(patient_representative)) cell_line[patient_representative][1]
                        else cell_line[1], .groups = "drop") %>% pull(cell_line)
dropped_lines <- setdiff(ann$cell_line, rep_keep)
stopifnot("expected 28 patient representatives" = length(rep_keep) == 28L)

site3_all <- factor(ann$site)
sub_all   <- ann$subtype
keep_i    <- match(rep_keep, ann$cell_line)
sub_rep   <- droplevels(sub_all[keep_i])
site3_rep <- droplevels(site3_all[keep_i])

cat(sprintf("\n=== Patient-representative set: %d lines -> %d patients ===\n",
            ncol(v), length(rep_keep)))
cat("Dropped sublines: ", paste(dropped_lines, collapse = ", "), "\n", sep = "")
cat("Subtype n, 31 lines vs 28 reps:\n")
print(as.data.frame(full_join(
  as_tibble(table(subtype = sub_all), .name_repair = "minimal") %>% rename(n_31_lines = n),
  as_tibble(table(subtype = sub_rep), .name_repair = "minimal") %>% rename(n_28_reps = n),
  by = "subtype")), row.names = FALSE)

# =============================================================================
# A. PCA / silhouette / PC confounder / commonality on the 28 representatives
# =============================================================================
# TWO variants of "recompute the PCA", because they are not the same thing and the
# difference is exactly the kind of choice that should not be silent:
#   refit_hvg : subset the VST to 28 lines and re-select the top-2000 HVG WITHIN
#               those lines, then PCA. This is script 02's recipe applied to the
#               28-line matrix and is the primary sensitivity result.
#   fixed_hvg : keep the HVGs chosen on all 31 lines and only drop the 3 columns.
#               Isolates the effect of removing the sublines from the effect of
#               re-selecting features.
pca_variant <- function(mat, refit_hvg) {
  if (refit_hvg) {
    rv <- rowVars(mat); top <- head(order(rv, decreasing = TRUE), N_HVG)
  } else {
    rv <- rowVars(v); top <- head(order(rv, decreasing = TRUE), N_HVG)
  }
  prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
}
r2  <- function(y, f) summary(lm(y ~ f))$r.squared
ar2 <- function(y, f) summary(lm(y ~ f))$adj.r.squared
r2j <- function(y, a, b) summary(lm(y ~ a + b))$r.squared
ar2j<- function(y, a, b) summary(lm(y ~ a + b))$adj.r.squared

sil_by_subtype <- function(pcx, subf) {
  d   <- dist(pcx[, 1:N_PC_SIL])
  sil <- silhouette(as.integer(subf), d)
  tibble(subtype = levels(subf)[sil[, "cluster"]], width = sil[, "sil_width"]) %>%
    group_by(subtype) %>% summarise(n = dplyr::n(), mean_sil = mean(width), .groups = "drop") %>%
    bind_rows(tibble(subtype = "ALL", n = nrow(pcx), mean_sil = mean(sil[, "sil_width"])))
}

pca31 <- prcomp(t(v[head(order(rowVars(v), decreasing = TRUE), N_HVG), ]),
                center = TRUE, scale. = FALSE)
v_rep <- v[, rep_keep, drop = FALSE]
pca28_refit <- pca_variant(v_rep, refit_hvg = TRUE)
pca28_fixed <- pca_variant(v_rep, refit_hvg = FALSE)

pc_stats <- function(pcx, subf, sitef, label) {
  map_dfr(1:3, function(k) {
    y <- pcx$x[, k]
    rj <- r2j(y, subf, sitef)
    tibble(set = label, PC = paste0("PC", k),
           var_pct        = 100 * pcx$sdev[k]^2 / sum(pcx$sdev^2),
           r2_subtype     = r2(y, subf), r2_site = r2(y, sitef),
           r2_joint       = rj, adjr2_joint = ar2j(y, subf, sitef),
           adjr2_subtype = ar2(y, subf), adjr2_site = ar2(y, sitef),
           adj_unique_subtype = ar2j(y, subf, sitef) - ar2(y, sitef),
           adj_unique_site = ar2j(y, subf, sitef) - ar2(y, subf),
           adj_shared = ar2(y, subf) + ar2(y, sitef) - ar2j(y, subf, sitef),
           unique_subtype = rj - r2(y, sitef),
           unique_site    = rj - r2(y, subf),
           shared         = r2(y, subf) + r2(y, sitef) - rj)
  })
}
pc_cmp <- bind_rows(
  pc_stats(pca31,       sub_all, site3_all, "31 lines"),
  pc_stats(pca28_refit, sub_rep, site3_rep, "28 reps (HVG refit)"),
  pc_stats(pca28_fixed, sub_rep, site3_rep, "28 reps (HVG fixed from 31)"))

sil_cmp <- bind_rows(
  sil_by_subtype(pca31$x,       sub_all) %>% mutate(set = "31 lines"),
  sil_by_subtype(pca28_refit$x, sub_rep) %>% mutate(set = "28 reps (HVG refit)"),
  sil_by_subtype(pca28_fixed$x, sub_rep) %>% mutate(set = "28 reps (HVG fixed from 31)"))

cat("\n=== A. PCA / commonality: 31 lines vs 28 patient representatives ===\n")
print(as.data.frame(pc_cmp %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat("\nSilhouette by subtype (Euclidean on PC1-10; ALL = overall mean):\n")
print(as.data.frame(sil_cmp %>% select(set, subtype, n, mean_sil) %>%
                      mutate(mean_sil = round(mean_sil, 4))), row.names = FALSE)

# Patient-aware reduced-model permutation tests. The full design is fixed;
# residuals from the reduced model are permuted within its categorical groups,
# then added back to fitted values (Freedman-Lane). One model per patient avoids
# exchanging related sublines as independent observations. Tests remain conditional
# on exchangeability within centre/histotype and cannot resolve empty design cells.
N_PC_PERM <- 9999L
pc_permutation <- function(term) {
  reduced_group <- if (term == "subtype") site3_rep else sub_rep
  full_q <- qr(model.matrix(~ sub_rep + site3_rep))
  red_q <- qr(model.matrix(~ reduced_group))
  Y <- pca28_refit$x[, 1:3, drop = FALSE]
  residual <- qr.resid(red_q, Y); fitted <- Y - residual
  df_num <- full_q$rank - red_q$rank; df_den <- nrow(Y) - full_q$rank
  f_stat <- function(y) {
    rss_full <- colSums(qr.resid(full_q, y)^2)
    rss_red <- colSums(qr.resid(red_q, y)^2)
    ((rss_red - rss_full) / df_num) / (rss_full / df_den)
  }
  observed <- f_stat(Y)
  set.seed(SEED)
  null <- replicate(N_PC_PERM, {
    ii <- seq_len(nrow(Y))
    for (g in levels(reduced_group)) {
      jj <- which(reduced_group == g); ii[jj] <- jj[sample.int(length(jj))]
    }
    f_stat(fitted + residual[ii, , drop = FALSE])
  })
  tibble(PC = paste0("PC", 1:3), tested_term = term,
         unit = "28 patient representatives; HVGs and PCA refitted", n_patients = nrow(Y),
         observed_F = observed, df_numerator = df_num, df_denominator = df_den,
         permutation_p = (1 + rowSums(null >= observed - 1e-12)) / (N_PC_PERM + 1),
         p_bh_3_PCs = p.adjust(permutation_p, "BH"), n_perm = N_PC_PERM, seed = SEED,
         method = "reduced-model residual permutation within reduced-factor groups; fixed full design")
}
readr::write_csv(bind_rows(pc_permutation("subtype"), pc_permutation("site")),
                 file.path(OUT, "sensitivity_patient_reps_pc_permutation.csv"))

# =============================================================================
# B. One-vs-rest DESeq2 + GO recovery on the 28 representatives
# =============================================================================
# This is the part nobody had computed. The published HGS contrast gives DESeq2 15
# "replicates" of which 3 are additional models from already-represented patients; here it gets 12 independent
# patients. Everything else (design, alpha, the up-signature cutoff, the GO-BP set
# construction, the recovery grading) is identical to script 03, so the ONLY thing
# that changes is the sample set.
dds_rep <- dds[, rep_keep]
dds_rep$subtype <- droplevels(factor(dds_rep$subtype, levels = SUBTYPE_LEVELS))
# re-derive the length-corrected normalization factors on the 28-column matrix
# (the stored ones were centred over 31 samples)
dds_rep <- estimateSizeFactors(dds_rep)
n_by_rep <- table(dds_rep$subtype)
n_by_all <- table(dds$subtype)
formal_ok_rep <- character() # Patient representatives remove pseudoreplication, not source confounding
cat(sprintf("\n=== B. DE on 28 reps | subtype n: %s | exploratory panel comparisons: %s ===\n",
            paste(names(n_by_rep), as.integer(n_by_rep), sep = "=", collapse = ", "),
            paste(names(n_by_rep), collapse = ", ")))

t2g <- readr::read_csv(file.path(OUT, "tx2gene_matched.csv"), show_col_types = FALSE)
id2sym <- t2g %>% distinct(ensembl_gene_id, external_gene_name) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  (\(d) setNames(d$external_gene_name, d$ensembl_gene_id))()

run_ovr_rep <- function(st) {
  d <- dds_rep
  d$grp <- factor(ifelse(d$subtype == st, st, "Other"), levels = c("Other", st))
  design(d) <- ~ grp
  d <- DESeq(d, quiet = TRUE)
  res <- results(d, contrast = c("grp", st, "Other"), alpha = 0.05)
  as_tibble(as.data.frame(res), rownames = "gene_id") %>%
    mutate(symbol = unname(id2sym[gene_id]), subtype = st,
           n_group = as.integer(n_by_rep[st]),
           evidence = "exploratory patient-representative association; source-confounded") %>% arrange(padj)
}
subtypes_rep <- levels(dds_rep$subtype)
de_rep <- setNames(lapply(subtypes_rep, run_ovr_rep), subtypes_rep)
readr::write_csv(bind_rows(de_rep), file.path(OUT, "sensitivity_patient_reps_de.csv"))

# --- GO-BP sets, built exactly as in script 03 -------------------------------
message("Building GO-BP gene sets from pinned org.Hs.eg.db ...")
go2eg  <- as.list(org.Hs.egGO2ALLEGS)
onto   <- AnnotationDbi::Ontology(names(go2eg))
bp_ids <- names(go2eg)[!is.na(onto) & onto == "BP"]
eg2sym <- unlist(as.list(org.Hs.egSYMBOL))
go_bp  <- lapply(go2eg[bp_ids], function(e) unique(na.omit(unname(eg2sym[unique(e)]))))
go_bp  <- go_bp[lengths(go_bp) >= 10 & lengths(go_bp) <= 500]
names(go_bp) <- paste0(names(go_bp), " ", AnnotationDbi::Term(names(go_bp)))
message(sprintf("GO-BP sets used (size 10-500): %d", length(go_bp)))

rank_syms <- function(tab) {
  x <- tab %>% filter(!is.na(stat), !is.na(symbol), symbol != "") %>%
    group_by(symbol) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  setNames(x$stat, x$symbol)
}
run_gsea_rep <- function(st) {
  set.seed(SEED)                                   # per-call, as in script 03
  fg <- suppressWarnings(fgsea(pathways = go_bp, stats = rank_syms(de_rep[[st]]),
                               minSize = 10, maxSize = 500, eps = 0,
                               nPermSimple = GSEA_NPERM_SIMPLE))
  fg %>% as_tibble() %>%
    mutate(subtype = st, go_id = sub(" .*", "", pathway),
           term = sub("^GO:[0-9]+ ", "", pathway),
           leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) %>%
    arrange(padj)
}
gsea_rep <- bind_rows(lapply(subtypes_rep, run_gsea_rep))
readr::write_csv(gsea_rep, file.path(OUT, "sensitivity_patient_reps_gsea_go.csv"))

# the same expected-program patterns as script 03, so the grades are comparable
expected <- tribble(
  ~st,      ~label,                       ~pattern,
  "HGS",    "DNA repair (HR/DSB)",        "double-strand break repair|homologous recombination|mismatch repair|DNA repair",
  "CC",     "oxidative/glutathione/detox","glutathione|xenobiotic|sulfur amino acid|reactive oxygen|oxidative stress",
  "MC",     "glycan/oligosaccharide",     "oligosaccharide|glycosylation|O-glycan|mucin|glycoprotein metabolic",
  "MMMT",   "EMT / migration / ECM",      "epithelial to mesenchymal|extracellular matrix|cell migration|cell motility",
  "SCCOHT", "cell cycle / mitosis",       "sister chromatid|mitotic|centromere|replication fork|chromosome segregation")
grade_recovery <- function(gsea_tbl, n_lookup) {
  map_dfr(seq_len(nrow(expected)), function(i) {
    hits <- gsea_tbl %>% filter(subtype == expected$st[i],
                                grepl(expected$pattern[i], term, ignore.case = TRUE))
    pos  <- hits %>% filter(NES > 0) %>% arrange(pval)
    best <- if (nrow(pos) > 0) pos[1, ] else hits %>% arrange(pval) %>% slice_head(n = 1)
    status <- if (nrow(pos) == 0) "not recovered (no positive-NES term)"
              else if (isTRUE(best$padj < 0.05)) "recovered (padj<0.05)"
              else if (isTRUE(best$pval < 0.05)) "suggestive (nominal p<0.05)"
              else "not recovered"
    tibble(subtype = expected$st[i], program = expected$label[i],
           n_group = as.integer(n_lookup[expected$st[i]]),
           n_pos_sig = sum(pos$padj < 0.05, na.rm = TRUE),
           best_term = best$term, NES = round(best$NES, 2),
           pval = signif(best$pval, 3), padj = signif(best$padj, 3), status = status)
  })
}
rec_rep <- grade_recovery(gsea_rep, n_by_rep)
rec_pub <- readr::read_csv(file.path(OUT, "rna_de_gsea_recovery.csv"), show_col_types = FALSE)

go_cmp <- rec_pub %>%
  transmute(subtype, program, n_group_31 = n_group, best_term_31 = best_term,
            NES_31 = NES, pval_31 = pval, padj_31 = padj, status_31 = status) %>%
  full_join(rec_rep %>%
              transmute(subtype, program, n_group_28 = n_group, best_term_28 = best_term,
                        NES_28 = NES, pval_28 = pval, padj_28 = padj, status_28 = status),
            by = c("subtype", "program")) %>%
  mutate(status_changed = status_31 != status_28,
         same_best_term = best_term_31 == best_term_28)
readr::write_csv(go_cmp, file.path(OUT, "sensitivity_patient_reps_go_comparison.csv"))

cat("\n=== B. GO recovery: 31 lines vs 28 patient representatives ===\n")
print(as.data.frame(go_cmp %>% select(subtype, program, n_group_31, n_group_28,
                                      padj_31, padj_28, status_31, status_28,
                                      status_changed)), row.names = FALSE)
if (any(go_cmp$status_changed))
  cat("!! The recovery grade CHANGES under patient collapse for the program(s) above.\n",
      "   Report both, or report the padj values rather than the grade.\n")

de_cmp <- bind_rows(
  bind_rows(de_rep) %>% group_by(subtype) %>%
    summarise(set = "28 reps", n_group = dplyr::first(n_group),
              n_padj05 = sum(padj < 0.05, na.rm = TRUE),
              n_up_sig = sum(padj < 0.05 & log2FoldChange > 1, na.rm = TRUE), .groups = "drop"),
  readr::read_csv(file.path(OUT, "rna_de_all.csv"), show_col_types = FALSE) %>%
    group_by(subtype) %>%
    summarise(set = "31 lines", n_group = dplyr::first(n_group),
              n_padj05 = sum(padj < 0.05, na.rm = TRUE),
              n_up_sig = sum(padj < 0.05 & log2FoldChange > 1, na.rm = TRUE), .groups = "drop")) %>%
  mutate(subtype = factor(subtype, levels = SUBTYPE_LEVELS)) %>% arrange(subtype, set)
cat("\nDE counts, 31 lines vs 28 reps (a drop is the expected direction: the published\n")
cat("  contrast counted 3 additional models from already-represented patients as independent replicates):\n")
print(as.data.frame(de_cmp), row.names = FALSE)

# =============================================================================
# C. One comparison table for the manuscript
# =============================================================================
cmp <- function(metric, v31, v28, note = NA_character_)
  tibble(analysis = NA_character_, metric = metric,
         value_31_lines = v31, value_28_reps = v28, delta = v28 - v31, note = note)
g <- function(tb, set, st, col) tb[[col]][tb$set == set & tb$subtype == st]
gp <- function(tb, set, pc, col) tb[[col]][tb$set == set & tb$PC == pc]

sens <- bind_rows(
  # --- PCA / silhouette
  cmp("PC1 variance explained (%)",
      gp(pc_cmp, "31 lines", "PC1", "var_pct"),
      gp(pc_cmp, "28 reps (HVG refit)", "PC1", "var_pct")) %>%
    mutate(analysis = "02 PCA"),
  cmp("PC2 variance explained (%)",
      gp(pc_cmp, "31 lines", "PC2", "var_pct"),
      gp(pc_cmp, "28 reps (HVG refit)", "PC2", "var_pct")) %>%
    mutate(analysis = "02 PCA"),
  cmp("PC1 R2 subtype",
      gp(pc_cmp, "31 lines", "PC1", "r2_subtype"),
      gp(pc_cmp, "28 reps (HVG refit)", "PC1", "r2_subtype")) %>%
    mutate(analysis = "02/17 PC confounder"),
  cmp("PC1 R2 site (3-level)",
      gp(pc_cmp, "31 lines", "PC1", "r2_site"),
      gp(pc_cmp, "28 reps (HVG refit)", "PC1", "r2_site")) %>%
    mutate(analysis = "02/17 PC confounder"),
  cmp("PC1 unique-subtype (commonality)",
      gp(pc_cmp, "31 lines", "PC1", "unique_subtype"),
      gp(pc_cmp, "28 reps (HVG refit)", "PC1", "unique_subtype")) %>%
    mutate(analysis = "17A commonality"),
  cmp("PC1 unique-site (commonality)",
      gp(pc_cmp, "31 lines", "PC1", "unique_site"),
      gp(pc_cmp, "28 reps (HVG refit)", "PC1", "unique_site")) %>%
    mutate(analysis = "17A commonality"),
  cmp("PC1 shared (commonality)",
      gp(pc_cmp, "31 lines", "PC1", "shared"),
      gp(pc_cmp, "28 reps (HVG refit)", "PC1", "shared")) %>%
    mutate(analysis = "17A commonality"),
  cmp("PC1 adjusted joint R2",
      gp(pc_cmp, "31 lines", "PC1", "adjr2_joint"),
      gp(pc_cmp, "28 reps (HVG refit)", "PC1", "adjr2_joint")) %>%
    mutate(analysis = "17A commonality"),
  bind_rows(lapply(c(SUBTYPE_LEVELS, "ALL"), function(st)
    cmp(sprintf("silhouette %s", st),
        g(sil_cmp, "31 lines", st, "mean_sil"),
        g(sil_cmp, "28 reps (HVG refit)", st, "mean_sil"),
        sprintf("n: %s -> %s", g(sil_cmp, "31 lines", st, "n"),
                g(sil_cmp, "28 reps (HVG refit)", st, "n"))) %>%
      mutate(analysis = "02 silhouette"))),
  # --- DE
  bind_rows(lapply(subtypes_rep, function(st)
    cmp(sprintf("DE genes padj<0.05 (%s one-vs-rest)", st),
        de_cmp$n_padj05[de_cmp$subtype == st & de_cmp$set == "31 lines"],
        de_cmp$n_padj05[de_cmp$subtype == st & de_cmp$set == "28 reps"],
        sprintf("n_group: %d -> %d", as.integer(n_by_all[st]), as.integer(n_by_rep[st]))) %>%
      mutate(analysis = "03 DE"))),
  bind_rows(lapply(subtypes_rep, function(st)
    cmp(sprintf("up-signature genes padj<0.05 & log2FC>1 (%s)", st),
        de_cmp$n_up_sig[de_cmp$subtype == st & de_cmp$set == "31 lines"],
        de_cmp$n_up_sig[de_cmp$subtype == st & de_cmp$set == "28 reps"]) %>%
      mutate(analysis = "03 signatures"))),
  # --- GO
  bind_rows(lapply(seq_len(nrow(go_cmp)), function(i)
    cmp(sprintf("GO recovery padj (%s: %s)", go_cmp$subtype[i], go_cmp$program[i]),
        go_cmp$padj_31[i], go_cmp$padj_28[i],
        sprintf("%s -> %s", go_cmp$status_31[i], go_cmp$status_28[i])) %>%
      mutate(analysis = "03 GO recovery"))))
sens <- sens %>%
  mutate(unit_31 = "cell line (n=31)", unit_28 = "patient representative (n=28)",
         across(c(value_31_lines, value_28_reps, delta), ~round(.x, 5))) %>%
  relocate(analysis, metric)
readr::write_csv(sens, file.path(OUT, "sensitivity_patient_reps.csv"))
readr::write_csv(pc_cmp,  file.path(OUT, "sensitivity_patient_reps_pca.csv"))
readr::write_csv(sil_cmp, file.path(OUT, "sensitivity_patient_reps_silhouette.csv"))

cat("\n=== C. Sensitivity comparison table (output/sensitivity_patient_reps.csv) ===\n")
print(as.data.frame(sens %>% select(analysis, metric, value_31_lines, value_28_reps, delta)),
      row.names = FALSE)

cat("\nOutputs: output/sensitivity_patient_reps{,_pca,_silhouette,_de,_gsea_go,_go_comparison}.csv\n")

# =============================================================================
# D. Environment record
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
write_session_info("21_rna_sensitivity")

message("\n21_rna_sensitivity.R complete.")
