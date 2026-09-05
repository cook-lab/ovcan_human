# =============================================================================
# Script: 01_rna_load_qc.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Load kallisto RNA-seq for the 31-line GENERATED analysis set,
#          summarize to gene level (pinned tx2gene, no live biomaRt), build a
#          DESeq2 object (proper counts + length offset), and report QC.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# Phase:   1 (RNA re-processing) — step 1 of 2 (this = load+QC; 02 = separation)
# =============================================================================
source("scripts/00_setup.R")
check_pkgs()
suppressPackageStartupMessages({
  library(tidyverse); library(tximport); library(DESeq2)
  library(matrixStats); library(jsonlite)   # jsonlite: parse kallisto run_info.json
  # biomaRt is loaded lazily only on a tx2gene cache miss (see section 2), so the
  # pipeline runs on machines without biomaRt whenever the cached map is present.
})

# 1. RNA analysis set (generated, retained) -----------------------------------
ss  <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE)
rna <- ss %>%
  filter(rna_seq == "Y", provenance == "generated", analysis_include == "Y") %>%
  mutate(
    subtype = factor(subtype, levels = c("HGS", "CC", "EC", "MC", "MMMT", "SCCOHT")),
    site    = factor(source_site),
    file    = file.path(DATA, "rna_seq", rna_sample_id, "abundance.tsv")
  )
stopifnot("Expected 31 RNA lines" = nrow(rna) == 31,
          "Missing abundance.tsv"  = all(file.exists(rna$file)))
files <- setNames(rna$file, rna$cell_line)
message(sprintf("RNA analysis set: %d lines | subtypes: %s",
                nrow(rna), paste(names(table(rna$subtype)), table(rna$subtype),
                                 sep = "=", collapse = ", ")))

# 2. Exact transcript reference recovered from the Ensembl release-93 cDNA FASTA
# Every versioned transcript ID and sequence length matches all 185,299 archived
# kallisto targets. The quantification date establishes no particular Ensembl
# release; the former "104-era" inference was incorrect. Re-summarising the
# existing transcript estimates with this exact map repairs the 3,529 targets
# lost by the release-105 map; raw-read re-quantification is unnecessary for this
# transcript-to-gene assignment correction. See data/reference/rna_reference_provenance.json.
T2G_MD5 <- "da6b6edd4d01d8d462802b7fe060875b"
T2G_NROW <- 185299L
t2g_ref <- file.path(PROJ, "data", "reference", "tx2gene_ensembl_rel93.csv")
t2g_file <- file.path(OUT, "tx2gene_matched.csv")
stopifnot("Restore the pinned exact transcript map in data/reference" = file.exists(t2g_ref),
          "Matched map checksum changed" = unname(tools::md5sum(t2g_ref)) == T2G_MD5)
t2g <- readr::read_csv(t2g_ref, show_col_types = FALSE)
stopifnot(nrow(t2g) == T2G_NROW,
          !anyDuplicated(t2g$transcript_id_versioned),
          !anyNA(t2g$ensembl_gene_id))
file.copy(t2g_ref, t2g_file, overwrite = TRUE)
tx2gene <- t2g %>% transmute(TXNAME = transcript_id_versioned, GENEID = ensembl_gene_id)
gene_meta <- t2g %>% distinct(ensembl_gene_id, external_gene_name, gene_biotype,
                              seq_region, primary_assembly)
readr::write_csv(gene_meta, file.path(OUT, "rna_gene_annotation.csv"))

# 3. Validate the full versioned target set and length in EVERY library before
# summarisation. Preserve primary and alternative reference loci as separately
# annotated Ensembl genes; do not silently drop counts assigned to either.
recon <- purrr::imap_dfr(files, function(f, cl) {
  ab <- readr::read_tsv(f, show_col_types = FALSE)
  idx <- match(ab$target_id, t2g$transcript_id_versioned)
  stopifnot("Versioned kallisto targets must match the complete reference" =
              nrow(ab) == nrow(t2g) && !anyDuplicated(ab$target_id) && !anyNA(idx),
            "Target lengths must match the reference FASTA" =
              all(ab$length == t2g$transcript_length[idx]))
  alt <- !t2g$primary_assembly[idx]
  tibble(cell_line = cl, n_index_targets = nrow(ab), n_unmapped = 0L,
         frac_targets_unmapped = 0, frac_tpm_unmapped = 0,
         frac_counts_unmapped = 0, versioned_ids_match = TRUE, lengths_match = TRUE,
         n_alternative_locus_targets = sum(alt),
         frac_tpm_alternative_loci = sum(ab$tpm[alt]) / sum(ab$tpm),
         frac_counts_alternative_loci = sum(ab$est_counts[alt]) / sum(ab$est_counts),
         transcript_estimated_count_total = sum(ab$est_counts),
         transcript_tpm_total = sum(ab$tpm),
         index_provenance = "All versioned IDs and sequence lengths match Ensembl release-93 GRCh38 cDNA FASTA",
         tx2gene_release = "Ensembl 93; exact versioned mapping; all reference loci retained")
})
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,
                ignoreTxVersion = FALSE, countsFromAbundance = "no")
recon <- recon %>% mutate(
  gene_estimated_count_total = colSums(txi$counts)[cell_line],
  gene_tpm_total = colSums(txi$abundance)[cell_line],
  count_fraction_retained = gene_estimated_count_total / transcript_estimated_count_total,
  tpm_fraction_retained = gene_tpm_total / transcript_tpm_total)
stopifnot("Gene summarisation must conserve all transcript estimated counts" =
            all(abs(recon$count_fraction_retained - 1) < 1e-8),
          "Gene summarisation must conserve transcript TPM" =
            all(abs(recon$tpm_fraction_retained - 1) < 1e-8))
readr::write_csv(recon, file.path(OUT, "rna_reference_reconciliation.csv"))
# Retain the old diagnostic filename with an explicit resolution rather than
# leaving stale tests of the now-repaired mismatch in the active output tree.
readr::write_csv(tibble(status = "resolved by exact versioned transcript mapping",
                        n_models = nrow(recon), n_targets = nrow(t2g),
                        n_genes = nrow(txi$counts), max_frac_unmapped = 0,
                        reference_release = 93L,
                        note = "All target IDs and lengths verified; all estimated counts and TPM retained"),
                 file.path(OUT, "rna_reference_sensitivity.csv"))
message(sprintf("Exact-reference gene matrix: %d genes x %d models; 100%% transcript estimates retained",
                nrow(txi$counts), ncol(txi$counts)))

# 4. DESeq2 object (counts + avgTxLength offset — the correct counts path) -----
dds <- DESeqDataSetFromTximport(txi, colData = as.data.frame(rna), design = ~ subtype)
n0  <- nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 2          # expressed in >=2 lines
dds  <- dds[keep, ]
dds  <- estimateSizeFactors(dds)
vsd  <- vst(dds, blind = TRUE)                    # for PCA/clustering (script 02)
message(sprintf("Gene filter: %d -> %d genes (>=10 counts in >=2 lines)", n0, nrow(dds)))

# 5. Save matrices/objects -----------------------------------------------------
tpm    <- txi$abundance
counts <- counts(dds, normalized = FALSE)
saveRDS(txi,  file.path(OUT, "rna_txi.rds"))
saveRDS(dds,  file.path(OUT, "rna_dds.rds"))
saveRDS(vsd,  file.path(OUT, "rna_vst.rds"))
readr::write_csv(as_tibble(tpm, rownames = "gene_id"),    file.path(OUT, "rna_tpm.csv"))
readr::write_csv(as_tibble(counts, rownames = "gene_id"), file.path(OUT, "rna_counts.csv"))
readr::write_csv(rna, file.path(OUT, "rna_sample_annotation.csv"))

# 6. QC metrics ----------------------------------------------------------------
#    TWO DEPTH COLUMNS, DO NOT CONFLATE THEM:
#      assigned_gene_counts  = colSums of the POST-FILTER gene-level tximport
#          estimated counts. This is an assignment total, NOT a sequencing depth:
#          it excludes unassigned/unmapped fragments and the genes dropped by the
#          >=10-in->=2 filter. (It was previously mislabelled `lib_size`, and the
#          manuscript read it as "median library size in reads".)
#      n_processed_fragments = kallisto's `n_processed` from run_info.json, i.e.
#          the true number of sequenced FRAGMENTS entering pseudoalignment. The
#          libraries are PAIRED-END, so these are fragments (read pairs), not
#          reads; a "reads" figure would be ~2x this.
run_info <- purrr::map_dfr(seq_len(nrow(rna)), function(i) {
  ri <- jsonlite::fromJSON(file.path(DATA, "rna_seq", rna$rna_sample_id[i], "run_info.json"))
  tibble(cell_line = rna$cell_line[i],
         n_processed_fragments = as.numeric(ri$n_processed),
         n_pseudoaligned       = as.numeric(ri$n_pseudoaligned),
         kallisto_version      = ri$kallisto_version,
         index_version         = as.integer(ri$index_version))
})
stopifnot("run_info missing for some RNA line" = nrow(run_info) == nrow(rna),
          "run_info index_version is not the archived index" =
            all(run_info$index_version == 10L))

qc <- tibble(
  cell_line   = rna$cell_line,
  subtype     = rna$subtype,
  site        = rna$site,
  pseudoalign = rna$rna_pseudoalign_pct,
  assigned_gene_counts = colSums(counts),
  detected    = colSums(counts >= 1)
) %>% left_join(run_info %>% select(cell_line, n_processed_fragments), by = "cell_line")
readr::write_csv(qc, file.path(OUT, "rna_qc_metrics.csv"))
cat("\n=== RNA QC summary ===\n")
print(qc %>% summarise(
  n = n(),
  pseudo_med = median(pseudoalign), pseudo_min = min(pseudoalign),
  assignedM_med = median(assigned_gene_counts)/1e6,
  fragM_med = median(n_processed_fragments)/1e6,
  fragM_min = min(n_processed_fragments)/1e6,
  fragM_max = max(n_processed_fragments)/1e6,
  det_med = median(detected)))

# 6b. Per-site QC comparison — TESTED, not asserted ----------------------------
#     The archived framing was that the site difference in pseudoalignment rate
#     "does not translate into a difference in genes detected". That is an
#     absence-of-effect claim drawn from two medians, and it does not hold: the
#     LOWER-alignment site detects MORE genes (Wilcoxon p ~ 0.008), and across all
#     31 lines pseudoalignment rate is NEGATIVELY correlated with gene detection.
#     Depth is reported here as a first-class metric because it does predict gene
#     detection line-by-line, but note what the tests actually support: the SITE
#     difference in depth is NOT significant (medians differ, distributions
#     overlap), so depth is a mechanism for the detection spread, not a
#     demonstrated site effect. Tests are Wilcoxon rank-sum on the two sites with
#     n>1 (Mes-Masson vs Huntsman); the n=1 site cannot enter a test. Both the
#     exact and the normal-approximation p are written for genes detected, because
#     the two differ slightly (0.0082 vs 0.0098) and both have been quoted.
site_tab <- qc %>% group_by(site) %>%
  summarise(n = n(),
            pseudoalign_median = median(pseudoalign),
            detected_median    = median(detected),
            assigned_gene_counts_median = median(assigned_gene_counts),
            n_processed_fragments_median = median(n_processed_fragments),
            .groups = "drop") %>% arrange(desc(n))
testable <- site_tab$site[site_tab$n > 1]
stopifnot("expected exactly 2 sites with n>1 for the rank-sum test" = length(testable) == 2L)
wx <- function(col, exact = TRUE) {
  a <- qc[[col]][qc$site == testable[1]]; b <- qc[[col]][qc$site == testable[2]]
  suppressWarnings(wilcox.test(a, b, exact = exact, correct = !exact)$p.value)
}
# pseudoalignment rate vs genes detected across ALL 31 lines: if alignment rate
# drove detection this would be positive. It is negative.
pa_det_r <- cor.test(qc$pseudoalign, qc$detected, method = "pearson")
pa_dep_r <- cor.test(qc$pseudoalign, qc$n_processed_fragments, method = "pearson")
dep_det_r<- cor.test(qc$n_processed_fragments, qc$detected, method = "pearson")

site_cmp <- bind_rows(
  site_tab %>% mutate(block = "per-site medians", statistic = NA_character_,
                      estimate = NA_real_, p_value = NA_real_, .before = 1),
  tibble(block = "site test (Wilcoxon rank-sum)",
         site = paste(testable, collapse = " vs "),
         n = sum(site_tab$n[site_tab$site %in% testable]),
         statistic = c("genes detected (exact)",
                       "genes detected (normal approx, continuity-corrected)",
                       "n_processed_fragments", "pseudoalignment %",
                       "assigned_gene_counts"),
         estimate = NA_real_,
         p_value = c(wx("detected"), wx("detected", exact = FALSE),
                     wx("n_processed_fragments"),
                     wx("pseudoalign"), wx("assigned_gene_counts"))),
  tibble(block = "correlation across all 31 lines", site = "all", n = nrow(qc),
         statistic = c("pearson(pseudoalign, detected)",
                       "pearson(pseudoalign, n_processed_fragments)",
                       "pearson(n_processed_fragments, detected)"),
         estimate = c(pa_det_r$estimate, pa_dep_r$estimate, dep_det_r$estimate),
         p_value  = c(pa_det_r$p.value,  pa_dep_r$p.value,  dep_det_r$p.value)))
readr::write_csv(site_cmp, file.path(OUT, "rna_qc_site_comparison.csv"))

cat("\nBy site (batch check — medians, then tests):\n")
print(as.data.frame(site_tab), row.names = FALSE)
cat(sprintf("\nWilcoxon %s vs %s: genes detected p=%.4g (normal approx %.4g) | fragments p=%.4g | assigned counts p=%.4g | pseudoalign%% p=%.4g\n",
            testable[1], testable[2], wx("detected"), wx("detected", exact = FALSE),
            wx("n_processed_fragments"), wx("assigned_gene_counts"), wx("pseudoalign")))
cat(sprintf("Across all %d lines: pseudoalign vs genes detected r=%.3f (p=%.3g) -> the lower-alignment\n",
            nrow(qc), pa_det_r$estimate, pa_det_r$p.value))
cat(sprintf("  site detects MORE genes, not fewer. Depth predicts detection line-by-line (r=%.3f, p=%.3g),\n",
            dep_det_r$estimate, dep_det_r$p.value))
cat(sprintf("  but the SITE depth difference is not significant (p=%.2f) - report depth as the mechanism for\n",
            wx("n_processed_fragments")))
cat("  the detection spread, NOT as an established site effect.\n")

# 7. QC figure -----------------------------------------------------------------
qc_long <- qc %>% mutate(cell_line = fct_reorder(cell_line, detected))
p <- ggplot(qc_long, aes(detected, cell_line, colour = subtype, shape = site)) +
  geom_point(size = 2.5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Genes detected (counts >= 1)", y = NULL,
       title = "RNA-seq QC — genes detected per line",
       subtitle = sprintf("n=%d generated lines; pseudoalignment median %.1f%%",
                           nrow(qc), median(qc$pseudoalign))) +
  theme_minimal(base_size = 11)
ggsave(file.path(FIGS, "01_rna_qc_detected.pdf"), p, width = 7, height = 7)

# 8. Environment record --------------------------------------------------------
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
write_session_info("01_rna_load_qc")

message("\n01_rna_load_qc.R complete. Outputs in output/ ; QC figure in figs/.")
