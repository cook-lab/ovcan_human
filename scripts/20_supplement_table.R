#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 20_supplement_table.R
#
# Consolidated per-line QC / metadata / identity SUPPLEMENT (peer-review add #11).
# One row per generated model (n=42), merging every per-line field a reuser wants:
# provenance + patient/family, per-assay availability, RNA/protein/WES QC, the
# corrected genomics (autosome FGA, hypermutator/MSI, TP53, tiered drivers),
# external identity (CCLE/DepMap-derived ConsensusOV, Cellosaurus STR), and the
# authentication calls. This is the reusable Table S1 backbone for the descriptor.
#
# INPUT : metadata/{samples.csv, line_family_map.csv} + output/{rna_qc_metrics,
#         prot_sample_qc, wes_mutation_load, wes_cnv_fga, wes_msi_mmr,
#         wes_driver_tiers, consensusov_calls, cellosaurus_str_status,
#         auth_perline_table}.csv
# OUTPUT: output/supplement_per_line.csv
#
# Pure join of already-materialized results (no recomputation); safe to run last.
#
# *** REVISED (review): DRIVER CALLS ARE NOW TIER-QUALIFIED ***
#   Usage Notes call this table "the recommended entry point for reuse", so it is the
#   artifact people will actually use — and it previously carried `key_drivers` with
#   NO TIER COLUMN AT ALL. That directly contradicted Technical Validation §3
#   ("defensible somatic BRCA1/2 is zero: both candidate BRCA2 calls are Tier 3 and
#   excluded"): TOV81D read `key_drivers = "BRCA2"` as its ONLY driver, and TOV3133D
#   read "BRCA2, CDK12, TP53" — exactly those two calls. Nine (line, gene) pairs were
#   Tier-3-only and appeared unqualified. We now join wes_driver_tiers.csv and emit
#   drivers_tier12 / drivers_tier3 / drivers_annotated, plus the hypermutator
#   passenger-risk `context` column that existed in the tier table but never reached
#   here. `key_drivers` (all tiers pooled) is retained for continuity but is no
#   longer the only driver column.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })

PROJ <- Sys.getenv("OVCAN_PROJ", unset = getwd())
op   <- function(...) file.path(PROJ, "output", ...)
mp   <- function(...) file.path(PROJ, "metadata", ...)
rd   <- function(f) readr::read_csv(f, show_col_types = FALSE)

# [integration revision] line_family_map.csv is GENERATED (by 15), not a source input,
# so a clean checkout that starts here has no map and this script used to abort on it.
# 07/08/16/17/21 call ensure_family_map() from 00_setup.R; this script deliberately
# sources nothing (it is a pure join and must stay cheap), so it regenerates the map
# the same way — by sourcing 15 — rather than duplicating the rule.
if (!file.exists(mp("line_family_map.csv"))) {
  message("metadata/line_family_map.csv absent -> running 15_patient_family_map.R")
  local(source(file.path(PROJ, "scripts", "15_patient_family_map.R"), local = TRUE))
}

need <- c(mp("line_family_map.csv"), mp("samples.csv"),
          op("rna_qc_metrics.csv"), op("prot_sample_qc.csv"),
          op("wes_mutation_load.csv"), op("wes_cnv_fga.csv"), op("wes_msi_mmr.csv"),
          op("wes_driver_tiers.csv"), op("consensusov_calls.csv"),
          op("cellosaurus_str_status.csv"), op("auth_perline_table.csv"))
missing <- need[!file.exists(need)]
if (length(missing))
  stop("Table S1 inputs missing (run the upstream scripts first):\n  ",
       paste(basename(missing), collapse = "\n  "))

fam <- rd(mp("line_family_map.csv"))              # 42 generated lines = the base
ss  <- rd(mp("samples.csv"))
rnaqc <- rd(op("rna_qc_metrics.csv"))
prqc  <- rd(op("prot_sample_qc.csv"))
load  <- rd(op("wes_mutation_load.csv"))
fga   <- rd(op("wes_cnv_fga.csv"))
msi   <- rd(op("wes_msi_mmr.csv"))
cons  <- rd(op("consensusov_calls.csv"))
str_  <- rd(op("cellosaurus_str_status.csv"))
auth  <- rd(op("auth_perline_table.csv"))
tiers <- rd(op("wes_driver_tiers.csv"))

# --- tier-resolved driver strings, per line ---------------------------------
# A (line, gene) pair takes its BEST (lowest) tier across that pair's calls. Genes
# with any Tier1-2 call are prioritized candidates; somatic/germline status
# remains unresolved at every tier in tumour-only WES.
drv <- tiers %>%
  mutate(tier_n = as.integer(str_replace(tier, "Tier", ""))) %>%
  group_by(cell_line, gene) %>%
  summarise(tier_best = min(tier_n), .groups = "drop") %>%
  group_by(cell_line) %>%
  summarise(
    drivers_tier12 = paste(sort(unique(gene[tier_best <= 2L])), collapse = ", "),
    drivers_tier3  = paste(sort(unique(gene[tier_best == 3L])), collapse = ", "),
    drivers_annotated = paste(sort(unique(ifelse(tier_best == 3L,
                                                 paste0(gene, "(T3)"), gene))), collapse = ", "),
    n_drivers_tier12 = sum(tier_best <= 2L),
    n_drivers_tier3  = sum(tier_best == 3L),
    .groups = "drop") %>%
  mutate(across(c(drivers_tier12, drivers_tier3),
                ~ifelse(.x == "", NA_character_, .x)))

# hypermutator passenger-risk context — one value per line in the tier table
drv_ctx <- tiers %>% filter(!is.na(context)) %>%
  group_by(cell_line) %>%
  summarise(driver_context = paste(unique(context), collapse = " | "), .groups = "drop")

sup <- fam %>%
  transmute(cell_line, subtype, source_site, patient_id, family,
            n_lines_in_family, patient_representative,
            has_rna, has_prot, has_wes_cnv, has_wes_maf) %>%
  # --- sample-sheet metadata -------------------------------------------------
  left_join(ss %>% transmute(cell_line,
                             subtype_conflict = subtype_status != "consensus",
                             all_three, tmt_plex, tmt_channel,
                             rna_passage, wes_passage), by = "cell_line") %>%
  # --- RNA QC ----------------------------------------------------------------
  # [review revision] 01_rna_load_qc.R renamed its `lib_size` column to
  # `assigned_gene_counts` and added the true sequenced-fragment count
  # `n_processed_fragments` (the old column was a SUM OF GENE-LEVEL ESTIMATED COUNTS
  # after filtering, not a read count, and was reported as "median library size of
  # 56.5M reads"). Both are carried here under names that say what they are. The
  # fallback keeps this join working against an older rna_qc_metrics.csv.
  left_join(rnaqc %>% transmute(cell_line,
                                rna_pseudoalign_pct = pseudoalign,
                                rna_assigned_gene_counts_M = round(
                                  (if ("assigned_gene_counts" %in% names(rnaqc)) assigned_gene_counts
                                   else lib_size) / 1e6, 1),
                                rna_sequenced_fragments_M = if ("n_processed_fragments" %in% names(rnaqc))
                                  round(n_processed_fragments / 1e6, 1) else NA_real_,
                                rna_genes_detected = detected), by = "cell_line") %>%
  # --- Protein QC ------------------------------------------------------------
  left_join(prqc %>% transmute(cell_line,
                               prot_n_detected = n_detected,
                               prot_pct_missing = round(pct_missing, 1)), by = "cell_line") %>%
  # --- WES: load / hypermutator / MSI ---------------------------------------
  left_join(load %>% transmute(cell_line,
                               wes_n_coding = n_coding,
                               wes_indel_frac = round(indel_frac, 3),
                               wes_tstv = round(tstv, 2),
                               wes_is_hypermutator = is_hypermutator), by = "cell_line") %>%
  left_join(msi %>% transmute(cell_line,
                              msi_mmr_call, sbs_best_cosmic = best_cosmic), by = "cell_line") %>%
  # --- WES: copy number (autosome-restricted) --------------------------------
  # NB fga_autosome has always been the autosome-restricted metric here; the C2 fix
  # was in script 10, which previously scored cnv_instability / genomics_consistent
  # on the chrX-inclusive column — that is why this table used to show OV90 with
  # fga_autosome = 0.269 beside genomics_consistent = "consistent" under a rule
  # stated as FGA > 0.30. Both columns now derive from the same autosome metric.
  left_join(fga %>% transmute(cell_line,
                              fga_autosome = round(fga_auto_0.2, 3),
                              fga_withX_legacy = round(fga_0.2, 3),
                              high_fga_flag), by = "cell_line") %>%
  # --- WES: drivers / authentication -----------------------------------------
  left_join(auth %>% transmute(cell_line,
                               TP53 = TP53, key_drivers,
                               expression_consistent, expression_basis,
                               genomics_consistent,
                               swisnf_deficient, swisnf_tier3_only_calls,
                               provisional_molecular_call, auth_flags = flags), by = "cell_line") %>%
  # --- WES: TIER-QUALIFIED driver calls (C3 fix) ------------------------------
  left_join(drv,     by = "cell_line") %>%
  left_join(drv_ctx, by = "cell_line") %>%
  # --- external identity -----------------------------------------------------
  left_join(cons %>% transmute(cell_line, consensusov_call,
                               consensusov_prob_top = if ("prob_top" %in% names(cons)) prob_top else NA_real_,
                               consensusov_margin = if ("margin_top_vs_second" %in% names(cons)) margin_top_vs_second else NA_real_,
                               consensusov_call_inherited = if ("consensusov_call_inherited" %in% names(cons)) consensusov_call_inherited else NA_character_,
                               consensusov_provenance = if ("consensusov_provenance" %in% names(cons)) consensusov_provenance else NA_character_),
            by = "cell_line") %>%
  left_join(str_ %>% transmute(cell_line, cellosaurus_accession,
                               str_profile_documented, problematic_flag, current_stock_str_status,
                               current_stock_mycoplasma_status), by = "cell_line") %>%
  arrange(match(subtype, c("HGS","LGS","CC","EC","MC","MMMT","SCCOHT")),
          patient_id, cell_line)

outfile <- op("supplement_per_line.csv")
write_csv(sup, outfile)

cat(sprintf("Wrote %s : %d lines x %d columns\n", outfile, nrow(sup), ncol(sup)))
cat(sprintf("  coverage: RNA %d | protein %d | WES-MAF %d | tri-omic %d\n",
            sum(sup$has_rna), sum(sup$has_prot), sum(sup$has_wes_maf),
            sum(sup$all_three == "Y", na.rm = TRUE)))
cat(sprintf("  STR documented (Cellosaurus): %d/42 | no public STR: %d\n",
            sum(sup$str_profile_documented, na.rm = TRUE),
            sum(!sup$str_profile_documented | is.na(sup$str_profile_documented))))
cat(sprintf("  hypermutator/MSI: %s\n",
            paste(sup$cell_line[which(sup$wes_is_hypermutator)], collapse = ", ")))
cat(sprintf("  high-FGA (median-centering caution): %s\n",
            paste(sup$cell_line[which(sup$high_fga_flag)], collapse = ", ")))

# --- C3 verification: no Tier-3 call is presented unqualified ----------------
cat("\n=== Tier-qualified driver calls (C3): lines with Tier-3-only genes ===\n")
t3 <- sup %>% filter(!is.na(drivers_tier3)) %>%
  select(cell_line, subtype, key_drivers, drivers_tier12, drivers_tier3, drivers_annotated)
print(as.data.frame(t3), row.names = FALSE)
cat(sprintf("  %d lines carry a Tier-3-only driver gene (%d such (line, gene) pairs).\n",
            nrow(t3), sum(sup$n_drivers_tier3, na.rm = TRUE)))
brca <- sup %>% filter(grepl("BRCA", key_drivers)) %>%
  select(cell_line, key_drivers, drivers_tier12, drivers_tier3)
cat("  BRCA1/2-bearing lines — Tier1-2 column must be BRCA-free for the paper's\n  'defensible somatic BRCA1/2 = 0' statement to hold:\n")
print(as.data.frame(brca), row.names = FALSE)
stopifnot("a BRCA1/2 call is now Tier1-2, contradicting 'defensible somatic BRCA1/2 = 0'" =
            !any(grepl("BRCA", brca$drivers_tier12)))
cat(sprintf("  driver_context (hypermutator passenger risk) carried for: %s\n",
            paste(sup$cell_line[!is.na(sup$driver_context)], collapse = ", ")))
cat(sprintf("  genomics_consistent tally: %s\n",
            paste(sprintf("%s %d", names(table(sup$genomics_consistent)),
                          as.integer(table(sup$genomics_consistent))), collapse = " / ")))
