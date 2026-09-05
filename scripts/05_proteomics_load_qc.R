# =============================================================================
# Script: 05_proteomics_load_qc.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Load the Morin-lab TMT proteomics (processed relative-abundance +
#          peptide tables), map channels -> cell lines via the TMT layout +
#          samples.csv, build a proteins x 31-line matrix, and assemble the
#          Technical-Validation QC we can from the available data:
#            - protein / peptide identifications
#            - missingness pattern (which is PER-PLEX in TMT, not per-sample)
#            - CV distribution (Morin-provided) + peptide coverage
#            - inter-plex daisy-chain BRIDGE replicate correlation (batch QC)
#            - principled missingness handling (presence threshold, not na.omit)
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   2 (Proteomics) — step 1 of 2 (this = load + QC; 06 = separation)
# -----------------------------------------------------------------------------
# CONSTRAINT (PI, 2026-07-23): TMT data were processed by the Morin lab. We hold
# the processed relative-abundance + peptide tables and their CV / peptide-
# coverage columns, but NOT the raw MS or a reproducible upstream pipeline.
# So we DOCUMENT the processing and assemble QC from what is available; we do
# NOT reconstruct the pipeline. Provenance of each QC metric is annotated below.
#
# WHAT THE MORIN LAB PROVIDED (documented, not recomputed here):
#   - Relative abundances are log2 reporter intensities normalized to a Pooled
#     Internal Standard (PIS, TMT ch1 in every plex) -> "relative abundance".
#   - qvalue column: proteins already filtered to <=0.01 (1% FDR).
#   - "Number of peptides" / "Number of unique peptides" / "Npeptides_quant":
#     peptide coverage per protein (min 2 peptides -> already peptide-filtered).
#   - "CV replicates": Morin's coefficient-of-variation across replicate
#     measurements (percent). This is the vendor CV distribution.
# WHAT WE COMPUTE HERE: protein/peptide counts, per-plex missingness, presence
#   threshold impact, and the bridge-replicate correlation (batch reproducibility).
# =============================================================================
source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(matrixStats)
})

PROT_DIR <- file.path(DATA, "proteomics")
f_abund  <- file.path(PROT_DIR, "protein_relative_abundance.xlsx")
f_pep    <- file.path(PROT_DIR, "peptide_ratio.xlsx")
f_layout <- file.path(PROT_DIR, "tmt.layout.xlsx")
stopifnot("Missing proteomics input(s)" = all(file.exists(f_abund, f_pep, f_layout)))

SUBTYPE_LEVELS <- c("HGS", "CC", "EC", "MC", "MMMT", "SCCOHT")

# =============================================================================
# 1. Sample sheet -> the 31 generated proteomics analysis lines
#    (EXCLUDE 5 external Carey LGS; bridges/SM.iD handled separately below)
# =============================================================================
ss <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE)
prot_ss <- ss %>%
  filter(proteomics == "Y", provenance == "generated", analysis_include == "Y") %>%
  mutate(subtype = factor(subtype, levels = SUBTYPE_LEVELS),
         site    = factor(source_site),
         plex    = factor(tmt_plex))
stopifnot("Expected 31 generated proteomics lines" = nrow(prot_ss) == 31L)
gen_lines <- prot_ss$cell_line
message(sprintf("Proteomics analysis set: %d generated lines | subtypes: %s",
                nrow(prot_ss),
                paste(names(table(prot_ss$subtype)), table(prot_ss$subtype),
                      sep = "=", collapse = ", ")))

# =============================================================================
# 2. TMT layout -> data-column map
#    Data-column headers in the abundance table are the layout `id` values
#    (verified: e.g. 'OV1369_R2', bridges 'VOA10816.1', spike-ins 'SM.iD.N').
#    Layout name uses '_' where samples.csv uses '-' (normalize on join).
#    Channel roles: TMT ch1 = PIS (denominator, absent from data), ch2-8 =
#    samples, ch9 = empty (absent), ch10 = inter-plex bridge (daisy-chain), and
#    ch11 = SM+iD spike-in. COV434 uniquely occupies plex-1 ch10 as a real
#    biological sample (no prior plex to bridge from), so it is a "sample".
# =============================================================================
lay <- read_excel(f_layout) %>%
  transmute(id, plex, tmt_label = `TMT.label`, name, type = Type, site_layout = Site,
            cell_line = gsub("_", "-", name))

abund <- read_excel(f_abund)
meta_cols <- names(abund)[1:11]                    # Symbol ... Npeptides_quant
data_cols <- names(abund)[12:ncol(abund)]          # 45 data columns == layout id
stopifnot("Expected 45 TMT data columns" = length(data_cols) == 45L)

colmap <- tibble(data_col = data_cols) %>%
  left_join(lay, by = c("data_col" = "id")) %>%
  mutate(role = case_when(
    grepl("^SM\\.iD", data_col)        ~ "SM.iD",   # spike-in standard (exclude)
    tmt_label == 10 & plex >= 2        ~ "bridge",  # daisy-chain technical rep
    TRUE                               ~ "sample")) # ch2-8 + COV434 (plex1 ch10)
stopifnot("Unmapped TMT data column(s)" = all(!is.na(colmap$plex)))

# attach provenance / canonical annotation from samples.csv
colmap <- colmap %>%
  left_join(ss %>% select(cell_line, provenance, analysis_include,
                          subtype_ss = subtype, site_ss = source_site,
                          tmt_channel),
            by = "cell_line")

biol_gen <- colmap %>%
  filter(role == "sample", provenance == "generated", analysis_include == "Y")
stopifnot("Biological generated columns != 31" = nrow(biol_gen) == 31L,
          "Column set mismatch vs samples.csv" =
            setequal(biol_gen$cell_line, gen_lines))
ext_lines <- colmap %>% filter(role == "sample", provenance == "external") %>% pull(cell_line)
message(sprintf("Columns classified: %d biological-generated, %d external(Carey LGS), %d daisy-chain bridges, %d SM.iD spike-in",
                nrow(biol_gen), length(ext_lines),
                sum(colmap$role == "bridge"), sum(colmap$role == "SM.iD")))

# =============================================================================
# 3. Build the proteins x 31 biological matrix
# -----------------------------------------------------------------------------
# FEATURE-COUNT HYGIENE (revision). The archived code set rownames with
# make.unique(abund$Symbol), which had two silent consequences:
#   (a) the 3 search rows with NO gene symbol became rows literally named
#       "NA" / "NA.1" / "NA.2". Downstream readers then dropped only the first
#       (readr parses the bare "NA" back to a missing value), so the deposited
#       denominator silently became 8,429 while the text still said 8,430 —
#       which is why "8,430 = 6,856 + 1,573" is off by one.
#   (b) the extra rows of a duplicated symbol became SYMBOL.1 / SYMBOL.2, which
#       every symbol-keyed analysis downstream orphans, while match() silently
#       kept whichever row happened to come first.
# Both are fixed here, and output/prot_feature_accounting.csv is written so every
# protein count quoted anywhere traces back to ONE denominator.
#   - NA-symbol rows are DROPPED with a logged count. They carry a Uniprot
#     accession but no gene symbol, so they cannot enter a symbol-keyed analysis.
#   - Duplicated symbols KEEP all their rows (distinct accessions, real
#     measurements), but exactly one is the SYMBOL REPRESENTATIVE and carries the
#     bare symbol as its row name. Deterministic rule, applied in order:
#       1. most non-missing values across the 31 lines
#       2. tie -> most quantified peptides (Npeptides_quant)
#       3. tie -> lowest q-value
#       4. tie -> first Uniprot accession alphabetically  (fully deterministic)
#     Non-representative rows are named SYMBOL|UNIPROT, so they are visibly a
#     second accession rather than a look-alike "SYMBOL.1", and are flagged
#     symbol_representative = FALSE in prot_qc.csv.
# =============================================================================
n_rows_search <- nrow(abund)
na_symbol     <- is.na(abund$Symbol) | abund$Symbol == ""
message(sprintf("Feature hygiene: %d search rows | %d with no gene symbol -> DROPPED (accessions: %s)",
                n_rows_search, sum(na_symbol),
                paste(abund$Uniprot[na_symbol], collapse = ", ")))
abund <- abund[!na_symbol, ]
stopifnot("Uniprot accession is not a unique key" = !any(duplicated(abund$Uniprot)))

feat <- tibble(symbol    = abund$Symbol,
               uniprot   = abund$Uniprot,
               n_peptides        = as.integer(abund$`Number of peptides`),
               n_unique_peptides = as.integer(abund$`Number of unique peptides`),
               npeptides_quant   = as.integer(abund$Npeptides_quant),
               qvalue            = as.numeric(abund$qvalue),
               cv_replicates     = as.numeric(abund$`CV replicates`),
               isoDoping         = as.logical(abund$isoDoping))

# order biological columns by subtype then cell line for readability
biol_gen <- biol_gen %>%
  mutate(subtype = factor(subtype_ss, levels = SUBTYPE_LEVELS)) %>%
  arrange(subtype, cell_line)
mat <- as.matrix(abund[, biol_gen$data_col])
colnames(mat) <- biol_gen$cell_line
mode(mat) <- "numeric"

# --- deterministic symbol representative (rule documented above) --------------
# Resolve on the FULL ordered key so the winner does not depend on table order.
feat$.n_present <- rowSums(!is.na(mat))
rep_idx <- feat %>% mutate(.i = row_number()) %>%
  arrange(symbol, dplyr::desc(.n_present), dplyr::desc(npeptides_quant),
          qvalue, uniprot) %>%
  group_by(symbol) %>% slice_head(n = 1) %>% ungroup() %>% pull(.i)
feat$symbol_representative <- seq_len(nrow(feat)) %in% rep_idx
feat$row <- ifelse(feat$symbol_representative, feat$symbol,
                   paste(feat$symbol, feat$uniprot, sep = "|"))
stopifnot("row names not unique after representative assignment" =
            !any(duplicated(feat$row)),
          "every symbol must have exactly one representative" =
            sum(feat$symbol_representative) == dplyr::n_distinct(feat$symbol))
rownames(mat) <- feat$row
feat <- feat %>% select(row, symbol, uniprot, symbol_representative,
                        everything(), -.n_present)
n_dup_symbols <- sum(duplicated(feat$symbol))
message(sprintf("Feature hygiene: %d rows retained | %d distinct symbols | %d non-representative duplicate-symbol rows (renamed SYMBOL|UNIPROT)",
                nrow(feat), dplyr::n_distinct(feat$symbol), n_dup_symbols))

# =============================================================================
# 4. Missingness — TMT missingness is PER-PLEX (all channels in a plex share the
#    same quantified/absent protein set). So per-SAMPLE missingness is degenerate
#    (identical within a plex); the meaningful unit is the plex (the MS run).
#    The archive notebook used na.omit (listwise deletion), which drops any
#    protein absent from ANY plex. We instead threshold on presence and DOCUMENT.
# =============================================================================
plex_of_col   <- setNames(as.integer(biol_gen$plex), biol_gen$cell_line)
plexes        <- sort(unique(plex_of_col))
present_line  <- rowSums(!is.na(mat))                      # per protein, 0..31
present_plex  <- sapply(plexes, function(p)
  rowSums(!is.na(mat[, names(plex_of_col)[plex_of_col == p], drop = FALSE])) > 0)
present_n_plex <- rowSums(present_plex)                    # 0..5

n_total    <- nrow(mat)
n_complete <- sum(present_line == 31L)                     # == na.omit (all 5 plexes)
n_pass50   <- sum(present_line >= 16L)                     # >= 50% of 31 lines
n_any      <- sum(present_line >= 1L)                      # quantified SOMEWHERE
n_zero     <- sum(present_line == 0L)                      # NA in all 31 lines
message(sprintf("Missingness: %d proteins total | complete/na.omit=%d | >=50%% lines=%d (+%d vs na.omit) | >=1 line=%d | 0 lines=%d",
                n_total, n_complete, n_pass50, n_pass50 - n_complete, n_any, n_zero))

feat <- feat %>%
  mutate(present_n_lines = present_line,
         present_n_plex  = present_n_plex,
         pass_presence50 = present_line >= 16L,
         complete_case   = present_line == 31L,
         # ZERO-PLEX proteins: in the search output but NA in ALL 31 analysis
         # lines, so they carry no quantification anywhere. They are RETAINED in
         # the deposited matrix (deleting rows from a deposited resource loses the
         # identification), but they must be flagged, because they sit inside the
         # "absent from at least one whole set" figure while being absent from all
         # five. Flag travels here and in prot_zero_plex_proteins.csv.
         zero_plex       = present_n_plex == 0L)

# --- 4b. Zero-plex characterisation -> prot_zero_plex_proteins.csv ------------
zero_tbl <- feat %>% filter(zero_plex) %>%
  select(row, symbol, uniprot, symbol_representative,
         n_peptides, n_unique_peptides, npeptides_quant, qvalue,
         cv_replicates, isoDoping, present_n_lines, present_n_plex)
readr::write_csv(zero_tbl, file.path(OUT, "prot_zero_plex_proteins.csv"))
cat(sprintf("\nZero-plex proteins (identified but quantified in 0 of 31 lines): %d\n", nrow(zero_tbl)))
if (nrow(zero_tbl))
  cat(sprintf("  peptides/protein median %g (range %g-%g); q-value median %.3g; all flagged zero_plex=TRUE in prot_qc.csv\n",
              median(zero_tbl$n_peptides, na.rm = TRUE),
              min(zero_tbl$n_peptides, na.rm = TRUE), max(zero_tbl$n_peptides, na.rm = TRUE),
              median(zero_tbl$qvalue, na.rm = TRUE)))

# --- 4c. Feature accounting -> prot_feature_accounting.csv --------------------
#     ONE table so every protein count in the manuscript is traceable. Any
#     "N proteins" statement anywhere should be quotable from a row of this file.
accounting <- tibble::tribble(
  ~level, ~n, ~definition,
  "search output rows",              n_rows_search,
    "rows in protein_relative_abundance.xlsx",
  "dropped: no gene symbol",         as.integer(sum(na_symbol)),
    "Uniprot accession present, Symbol missing; cannot enter a symbol-keyed analysis",
  "analysis rows (matrix)",          n_total,
    "rows of prot_abundance_matrix.csv / prot_matrix.rds",
  "distinct gene symbols",           as.integer(dplyr::n_distinct(feat$symbol)),
    "unique symbols among the analysis rows",
  "non-representative duplicate rows", as.integer(n_dup_symbols),
    "second+ accession for a symbol; row named SYMBOL|UNIPROT, symbol_representative=FALSE",
  "quantified in >=1 line",          n_any,
    "present_n_lines >= 1",
  "quantified in 0 lines (zero-plex)", n_zero,
    "present_n_plex == 0; retained in the matrix, flagged zero_plex=TRUE",
  "complete case (all 5 plexes)",    n_complete,
    "present in all 31 lines; == the archived na.omit set; PCA/variance-partition input",
  "absent from >=1 whole plex",       as.integer(n_total - n_complete),
    "present_n_plex < 5; INCLUDES the zero-plex proteins",
  "absent from >=1 plex, excl. zero-plex", as.integer(n_total - n_complete - n_zero),
    "present_n_plex in 1..4",
  "presence filter >=50% of lines",  n_pass50,
    "present_n_lines >= 16; the recommended reuse filter",
  "peptide-table accessions",        NA_integer_,
    "filled in section 5 below")
cat("\n=== Feature accounting (single denominator for every protein count) ===\n")
print(as.data.frame(accounting), row.names = FALSE)
stopifnot("complete + absent-from-any must equal the analysis row count" =
            n_complete + (n_total - n_complete) == n_total)

# =============================================================================
# 5. Peptide table — count + q-value distribution + independent coverage.
#    Read only the 4 annotation columns (skip the 45 data cols) for speed.
# =============================================================================
pep <- read_excel(f_pep,
                  col_types = c("text", "text", "numeric", "numeric",
                                rep("skip", 45)))
n_peptides_total <- nrow(pep)
pep_per_prot     <- dplyr::count(pep, Accession, name = "n_pep")
message(sprintf("Peptides quantified: %d across %d protein accessions (median %g peptides/protein)",
                n_peptides_total, nrow(pep_per_prot), median(pep_per_prot$n_pep)))

# Reconcile the peptide table against the abundance table — the two counts differ
# and the archived script printed both without ever relating them.
accounting$n[accounting$level == "peptide-table accessions"] <- nrow(pep_per_prot)
accounting <- accounting %>%
  add_row(level = "peptide accessions absent from the analysis matrix",
          n = as.integer(sum(!pep_per_prot$Accession %in% feat$uniprot)),
          definition = "in peptide_ratio.xlsx but not in the abundance table (incl. the dropped no-symbol rows)") %>%
  add_row(level = "peptides quantified", n = as.integer(n_peptides_total),
          definition = "rows in peptide_ratio.xlsx")
readr::write_csv(accounting, file.path(OUT, "prot_feature_accounting.csv"))
cat(sprintf("\nPeptide/abundance reconciliation: %d peptide accessions vs %d matrix rows; %d peptide accessions have no abundance row\n",
            nrow(pep_per_prot), n_total,
            sum(!pep_per_prot$Accession %in% feat$uniprot)))

# =============================================================================
# 6. Bridge-replicate correlation — the batch-reproducibility metric.
#    Each daisy-chain bridge (plex k ch10) re-runs one sample from plex k-1.
#    Correlate the sample's primary column vs its bridge column across proteins.
#    NOTE: VOA3993 is an external Carey LGS line; its bridge pair is still a
#    valid TECHNICAL replicate for reproducibility QC (flagged 'external').
# =============================================================================
bridges <- colmap %>% filter(role == "bridge") %>%
  transmute(name, bridge_col = data_col, bridge_plex = plex)
prims   <- colmap %>% filter(role == "sample") %>%
  transmute(name, prim_col = data_col, prim_plex = plex)
bpairs  <- bridges %>% left_join(prims, by = "name")
stopifnot("Bridge without matching primary" = all(!is.na(bpairs$prim_col)))

bridge_cor <- lapply(seq_len(nrow(bpairs)), function(i) {
  x  <- as.numeric(abund[[bpairs$prim_col[i]]])
  y  <- as.numeric(abund[[bpairs$bridge_col[i]]])
  ok <- is.finite(x) & is.finite(y)
  tibble(cell_line   = gsub("_", "-", bpairs$name[i]),
         prim_plex   = bpairs$prim_plex[i],
         bridge_plex = bpairs$bridge_plex[i],
         n_proteins  = sum(ok),
         pearson     = cor(x[ok], y[ok], method = "pearson"),
         spearman    = cor(x[ok], y[ok], method = "spearman"),
         external    = gsub("_", "-", bpairs$name[i]) %in% ext_lines)
}) %>% bind_rows() %>% arrange(bridge_plex)
readr::write_csv(bridge_cor, file.path(OUT, "prot_bridge_cor.csv"))

# =============================================================================
# 7. Per-line QC (note the per-plex degeneracy) and per-protein QC table
# =============================================================================
sample_qc <- biol_gen %>%
  transmute(cell_line, subtype, site = site_ss, plex, channel = tmt_channel) %>%
  mutate(n_detected = colSums(!is.na(mat))[cell_line],
         pct_missing = 100 * (1 - n_detected / n_total))
readr::write_csv(sample_qc, file.path(OUT, "prot_sample_qc.csv"))

readr::write_csv(feat, file.path(OUT, "prot_qc.csv"))

# =============================================================================
# 8. Save matrix + annotation. The deposited/documented resource matrix keeps
#    every symbol-carrying search row with NAs preserved (see the feature-hygiene
#    note in section 3: the 3 no-symbol rows are dropped, duplicate-symbol rows
#    are kept but renamed SYMBOL|UNIPROT). Presence, representative and zero-plex
#    flags live in prot_qc.csv; the count reconciliation lives in
#    prot_feature_accounting.csv. Quote counts from those, not from memory.
# =============================================================================
readr::write_csv(as_tibble(mat, rownames = "protein"),
                 file.path(OUT, "prot_abundance_matrix.csv"))
readr::write_csv(
  biol_gen %>% transmute(cell_line, subtype, site = site_ss, plex,
                         channel = tmt_channel, data_col, provenance),
  file.path(OUT, "prot_sample_annotation.csv"))
saveRDS(list(mat = mat, feat = feat, sample_ann = biol_gen),
        file.path(OUT, "prot_matrix.rds"))

# =============================================================================
# 9. QC figures (colorblind-safe; show distributions, not just summaries)
# =============================================================================
sub_cols  <- setNames(RColorBrewer::brewer.pal(6, "Dark2"), SUBTYPE_LEVELS)
plex_cols <- setNames(viridisLite::viridis(5), as.character(1:5))

# 9a. Missingness / per-plex structure ---------------------------------------
plex_prot <- tibble(plex = factor(plexes),
                    n_proteins = colSums(present_plex),
                    n_lines = as.integer(table(plex_of_col)[as.character(plexes)]))
p_miss <- ggplot(plex_prot, aes(plex, n_proteins, fill = plex)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = n_complete, linetype = "dashed", colour = COOK_NAVY) +
  geom_hline(yintercept = n_total, linetype = "dotted", colour = "grey50") +
  geom_text(aes(label = n_proteins), vjust = -0.4, size = 3) +
  scale_fill_manual(values = plex_cols, guide = "none") +
  annotate("text", x = 0.7, y = n_complete, label = sprintf("complete (all plexes) = %d", n_complete),
           vjust = -0.4, hjust = 0, size = 2.8, colour = COOK_NAVY) +
  annotate("text", x = 0.7, y = n_total, label = sprintf("total = %d", n_total),
           vjust = -0.4, hjust = 0, size = 2.8, colour = "grey40") +
  labs(x = "TMT plex (MS run)", y = "Proteins quantified",
       title = "Proteomics missingness is per-plex, not per-sample",
       subtitle = sprintf("Within a plex all channels share the quantified set; %d lines total",
                          nrow(biol_gen))) +
  theme_minimal(base_size = 11)

tier <- tibble(present_n_plex = factor(0:5, levels = 0:5),
               n_proteins = as.integer(table(factor(feat$present_n_plex, levels = 0:5))))
p_tier <- ggplot(tier, aes(present_n_plex, n_proteins)) +
  geom_col(width = 0.7, fill = COOK_RUST) +
  geom_text(aes(label = n_proteins), vjust = -0.4, size = 3) +
  labs(x = "# plexes a protein is quantified in (0-5)", y = "Proteins",
       title = "Protein detection across plexes",
       subtitle = "Peak at 5/5 = the complete-case set; presence threshold recovers the mid-tiers") +
  theme_minimal(base_size = 11)
ggsave(file.path(FIGS, "05_prot_missingness.pdf"), p_miss,  width = 7.5, height = 5)
ggsave(file.path(FIGS, "05_prot_presence_tiers.pdf"), p_tier, width = 7, height = 4.5)

# 9b. CV distribution (Morin) + peptide coverage ------------------------------
cv_df  <- feat %>% filter(is.finite(cv_replicates))
pep_df <- feat %>% filter(is.finite(n_peptides))
cv_med  <- median(cv_df$cv_replicates)
pep_med <- median(pep_df$n_peptides)
p_cv <- ggplot(cv_df, aes(cv_replicates)) +
  geom_histogram(bins = 60, fill = COOK_NAVY, colour = NA) +
  geom_vline(xintercept = cv_med, colour = COOK_RUST, linetype = "dashed") +
  scale_x_continuous(limits = c(0, quantile(cv_df$cv_replicates, 0.99))) +
  labs(x = "CV across replicates (%) - Morin-provided", y = "Proteins",
       title = "Replicate CV distribution (vendor QC)",
       subtitle = sprintf("median = %.1f%% (n=%d proteins with CV; x truncated at 99th pct)",
                          cv_med, nrow(cv_df))) +
  theme_minimal(base_size = 11)
p_pep <- ggplot(pep_df, aes(n_peptides)) +
  geom_histogram(bins = 60, fill = COOK_NAVY, colour = NA) +
  geom_vline(xintercept = pep_med, colour = COOK_RUST, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Peptides per protein (log10) - Morin-provided", y = "Proteins",
       title = "Peptide coverage per protein",
       subtitle = sprintf("median = %g peptides/protein; min 2 (already peptide-filtered)", pep_med)) +
  theme_minimal(base_size = 11)
ggsave(file.path(FIGS, "05_prot_cv_distribution.pdf"), p_cv,  width = 6.5, height = 4.5)
ggsave(file.path(FIGS, "05_prot_peptide_coverage.pdf"), p_pep, width = 6.5, height = 4.5)

# 9c. Bridge-replicate scatter (batch reproducibility) ------------------------
bridge_long <- lapply(seq_len(nrow(bpairs)), function(i) {
  cl <- gsub("_", "-", bpairs$name[i])
  tibble(cell_line = cl,
         primary   = as.numeric(abund[[bpairs$prim_col[i]]]),
         bridge    = as.numeric(abund[[bpairs$bridge_col[i]]]))
}) %>% bind_rows() %>% filter(is.finite(primary), is.finite(bridge))
lab_df <- bridge_cor %>%
  mutate(cell_line = factor(cell_line, levels = bridge_cor$cell_line),
         lab = sprintf("r=%.3f, rho=%.3f\nn=%d%s",
                       pearson, spearman, n_proteins,
                       ifelse(external, "  (external)", "")))
bridge_long$cell_line <- factor(bridge_long$cell_line, levels = bridge_cor$cell_line)
p_bridge <- ggplot(bridge_long, aes(primary, bridge)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
  geom_point(size = 0.4, alpha = 0.25, colour = COOK_NAVY) +
  geom_text(data = lab_df, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.05, vjust = 1.2, size = 2.9, colour = COOK_RUST, inherit.aes = FALSE) +
  facet_wrap(~ cell_line, nrow = 1) +
  labs(x = "log2 relative abundance - primary channel",
       y = "bridge channel (next plex, ch10)",
       title = "Inter-plex daisy-chain bridge reproducibility",
       subtitle = "Same physical sample re-run in the following plex; each point = one protein") +
  coord_equal() + theme_minimal(base_size = 10)
ggsave(file.path(FIGS, "05_prot_bridge_cor.pdf"), p_bridge, width = 12, height = 3.6)

# =============================================================================
# 10. Console report
# =============================================================================
cat("\n=== PROTEOMICS QC SUMMARY (script 05) ===\n")
cat(sprintf("Proteins quantified (>=1 line): %d / %d in table\n", n_any, n_total))
cat(sprintf("Peptides quantified: %d\n", n_peptides_total))
cat(sprintf("Complete-case (all 5 plexes; == na.omit): %d\n", n_complete))
cat(sprintf(">=50%% of 31 lines (presence threshold): %d  (+%d vs na.omit, +%.1f%%)\n",
            n_pass50, n_pass50 - n_complete, 100 * (n_pass50 - n_complete) / n_complete))
cat(sprintf("Replicate CV (Morin): median %.2f%%, IQR %.2f-%.2f%%\n",
            median(cv_df$cv_replicates),
            quantile(cv_df$cv_replicates, .25), quantile(cv_df$cv_replicates, .75)))
cat("\nPer-plex proteins quantified (missingness unit):\n"); print(as.data.frame(plex_prot))
cat("\nBridge-replicate correlation (inter-plex technical reproducibility):\n")
print(as.data.frame(bridge_cor), digits = 4)
cat(sprintf("\nBridge summary: Pearson %.3f-%.3f, Spearman %.3f-%.3f\n",
            min(bridge_cor$pearson), max(bridge_cor$pearson),
            min(bridge_cor$spearman), max(bridge_cor$spearman)))
cat("\nNote: per-SAMPLE missingness is degenerate (identical within a plex); the\n",
    "meaningful missingness unit is the plex. isoDoping flags", sum(feat$isoDoping),
    "spike-in-doped proteins (retained; flagged in prot_qc.csv).\n")

# =============================================================================
# 11. Environment record (see the note in 01_rna_load_qc.R section 8)
# =============================================================================
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
write_session_info("05_proteomics_load_qc")

message("\n05_proteomics_load_qc.R complete. Outputs in output/prot_* ; QC figures in figs/05_*.")
