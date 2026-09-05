# =============================================================================
# Script: 16_wes_signatures_msi.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Panel-wide mutation LOAD, MMR/MSI assessment, and mutational-context
#          (SBS) profiling — driven by peer-review 3.8: TOV21G carries a mutation
#          load >3x the next line and was previously unflagged. Decide whether
#          TOV21G is best described as a candidate MMR-deficient / MSI-high (or
#          POLE-mutant) clear-cell model vs a tumor-only artifact, on converging
#          evidence: load, indel burden, MMR-gene status, and the SBS-96 spectrum.
# Author:  Cook Lab (analyst: Claude)  |  Date: 2026-07-23
# =============================================================================
#
# *** ASSAY CAVEAT — TUMOR-ONLY WES, NO MATCHED NORMAL ***
#   These are Mutect2 tumor-only calls (no matched normal). Absolute mutation
#   burden / TMB and definitive MSI status CANNOT be called from these data
#   (residual germline persists in ~100%-pure lines; Halperin 2017; Little 2021).
#   We therefore treat load as a RELATIVE within-panel metric and MSI/MMR status
#   as QUALITATIVE, built from CONVERGING proxies (load + indel fraction +
#   MMR-gene status + SBS-96 context), not a clinical MSI assay. Calibrated
#   language throughout; the SBS spectrum is the most discriminating proxy.
#
# *** GENOME BUILD — VERIFIED GRCh38/hg38 (do not trust the MAF header) ***
#   The archived MAF `NCBI_Build` column literally reads "GRCh37", but that is a
#   spurious vcf2maf default: the calls are GRCh38. Evidence (this project):
#     - Mutect2 PoN = 1000g_pon.hg38.vcf.gz; VCF ##contig chr1 length 248,956,422
#       (= GRCh38; GRCh37 chr1 = 249,250,621).
#     - Driver hotspots sit at GRCh38 coordinates: KRAS G12 chr12:25,245,350
#       (GRCh37 = 25,398,284); PIK3CA chr3:179.2 Mb (GRCh37 ~178.9 Mb);
#       TP53 chr17:7.67 Mb (GRCh37 ~7.57 Mb).
#   => trinucleotide context uses BSgenome.Hsapiens.UCSC.hg38. This build
#      mislabel is a deposition item (flag for the descriptor).
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(patchwork); library(ggrepel)
})
set.seed(SEED)

# ---- Lab theme + palettes (match scripts 07/12/14) --------------------------
theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace% theme(
    text = element_text(colour = "black"),
    plot.title = element_text(size = rel(1.1), hjust = 0, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = rel(0.82), hjust = 0, colour = "grey30", margin = margin(b = 6)),
    axis.title = element_text(size = rel(0.95)), axis.text = element_text(size = rel(0.8), colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    legend.title = element_text(size = rel(0.8)), legend.text = element_text(size = rel(0.75)),
    legend.key = element_blank(), legend.background = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid = element_blank(), strip.background = element_blank(),
    strip.text = element_text(size = rel(0.85), face = "bold"), plot.margin = margin(6, 8, 6, 6))
}
# Subtype palette — identical to the 07 oncoplot (Okabe-Ito; no red-green pair).
sub_cols <- c(HGS = "#0072B2", CC = "#E69F00", EC = "#009E73",
              MC = "#CC79A7", LGS = "#56B4E9", MMMT = "#D55E00", SCCOHT = "#000000")

# ---- Config -----------------------------------------------------------------
MMR_ENZYME <- c("MLH1","MSH2","MSH6","PMS2")             # direct MMR enzymes (Lynch core)
MMR_INDIRECT <- c("EPCAM")                               # 3' deletions silence MSH2 (indirect)
MMR_OTHER <- c("PMS1","MSH3")                            # secondary MMR
PROOF     <- c("POLE","POLD1")                           # replicative proofreading
POLE_EXO  <- c(268L, 471L)                               # POLE exonuclease domain (codons)
TRUNCATING <- c("Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation",
                "Splice_Site","Translation_Start_Site","Nonstop_Mutation")
POP_AF_MAX <- 0.001                                      # match script 07 germline filter
# COSMIC diagnostic signature groups (SBS).
# *** EXACT REFERENCE RELEASE (review revision) *** — the reference version was
# never recorded in any output. COSMIC_v3.2 is the newest COSMIC release bundled
# with MutationalPatterns (3.18.0 offers COSMIC / SIGNAL / SPARSE / COSMIC_v3.1 /
# COSMIC_v3.2 only), so v3.3+ signatures are NOT available in this environment.
# Both strings are written into wes_msi_mmr.csv and wes_sbs_cosine.csv.
COSMIC_SOURCE <- "COSMIC_v3.2"
COSMIC_GENOME <- "GRCh38"
MMR_D_SIGS <- c("SBS6","SBS14","SBS15","SBS20","SBS21","SBS26","SBS44")  # defective MMR / MSI
POLE_SIGS  <- c("SBS10a","SBS10b","SBS10c","SBS10d","SBS28")             # POLE/POLD1 proofreading
CLOCK_SIGS <- c("SBS1","SBS5","SBS40")                                   # clock-like (age)
# COSINE CALL RULE (review revision). The previous rule thresholded the MMR-d GROUP
# MAXIMUM at an absolute 0.75 with no justification. SBS1/5/6/15/44 are mutually
# similar, so an absolute cosine cut carries little information; what the data
# actually support is a PANEL-RELATIVE separation — TOV21G is the only line whose
# MMR-d similarity exceeds its clock-like similarity at all (0.877 vs 0.651; the
# other 21 lines sit at MMR-d 0.51-0.70 with clock > MMR-d). The call is therefore
# a MARGIN over the clock group, and the margin itself is reported per line. The old
# absolute screen is retained as a reported column (cos_mmr_d_ge_075), not a gate.
SBS_MARGIN_MIN <- 0.10   # cos(group) - cos(clock) required for a candidate call
SBS_ABS_LEGACY <- 0.75   # legacy absolute screen — reported, no longer decisive

# =============================================================================
# 1. LOAD — filtered MAF (coding non-synonymous PASS candidates) + metadata
#    (read once into memory; wes_mutations_filtered.csv is owned by another agent)
# =============================================================================
stopifnot("output/wes_mutations_filtered.csv missing — run 07_wes_mutations.R first" =
            file.exists(file.path(OUT, "wes_mutations_filtered.csv")))
maf <- readr::read_csv(file.path(OUT, "wes_mutations_filtered.csv"), show_col_types = FALSE)
ensure_family_map()          # [review revision] 15 writes this map; see 00_setup.R
fam <- readr::read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)
message(sprintf("Filtered MAF: %d coding-nonsyn candidates across %d lines",
                nrow(maf), dplyr::n_distinct(maf$cell_line)))

# =============================================================================
# 2. PER-LINE MUTATION LOAD — SNV / indel / MNV, indel fraction, Ts/Tv, flag
# =============================================================================
is_transition <- function(ref, alt) {
  pair <- paste0(pmin(ref, alt), pmax(ref, alt))
  pair %in% c("AG", "CT")           # A<->G, C<->T
}
# Ts/Tv on single-base substitutions only
snv <- maf %>%
  filter(Variant_Type == "SNP",
         Reference_Allele  %in% c("A","C","G","T"),
         Tumor_Seq_Allele2 %in% c("A","C","G","T")) %>%
  mutate(is_ts = is_transition(Reference_Allele, Tumor_Seq_Allele2))
tstv <- snv %>% group_by(cell_line) %>%
  summarise(n_ts = sum(is_ts), n_tv = sum(!is_ts), .groups = "drop") %>%
  mutate(tstv = ifelse(n_tv > 0, n_ts / n_tv, NA_real_),
         ts_frac = n_ts / (n_ts + n_tv))

load_line <- maf %>%
  group_by(cell_line, subtype) %>%
  summarise(
    n_coding   = n(),
    n_snv      = sum(Variant_Type == "SNP"),
    n_ins      = sum(Variant_Type == "INS"),
    n_del      = sum(Variant_Type == "DEL"),
    n_indel    = sum(Variant_Type %in% c("INS","DEL")),
    n_mnv      = sum(Variant_Type %in% c("DNP","TNP","ONP")),
    .groups = "drop") %>%
  mutate(indel_frac = n_indel / (n_snv + n_indel)) %>%
  left_join(tstv %>% select(cell_line, n_ts, n_tv, tstv, ts_frac), by = "cell_line")

# Robust within-panel outlier statistics (coding-nonsyn load is the review metric).
# [review revision] The hypermutator flag now comes from the SHARED rule in
# 00_setup.R (hypermutator_stats), which 07 also calls — previously 07 used
# ">3x median" and this script used "robust_z > 5 & fold > 3", i.e. one
# manuscript-facing flag with two implementations.
hm_stats <- hypermutator_stats(load_line$n_coding, load_line$cell_line)
med_load <- hm_stats$median_load[1]
mad_load <- hm_stats$mad_load[1]                    # median absolute deviation (x1.4826)
load_line <- load_line %>%
  left_join(hm_stats %>% select(cell_line, fold_over_median, robust_z, is_hypermutator),
            by = "cell_line") %>%
  arrange(desc(n_coding)) %>%
  mutate(load_rank = row_number())

# =============================================================================
# 3. PER-PATIENT COLLAPSE (avoid pseudoreplication of same-patient sublines)
#    Rule: within a patient, the MAF-bearing "patient_representative" if it has a
#    MAF, else the alphabetically-first MAF-bearing subline. TOV21G is a singleton.
# =============================================================================
fam_maf <- fam %>% filter(cell_line %in% load_line$cell_line)      # 22 MAF lines
wes_rep <- fam_maf %>%
  group_by(patient_id) %>%
  mutate(is_wes_patient_rep =
           if (any(patient_representative)) cell_line == cell_line[patient_representative][1]
           else cell_line == sort(cell_line)[1],
         n_family_maf_lines = n()) %>%
  ungroup() %>%
  select(cell_line, patient_id, family, n_family_maf_lines, is_wes_patient_rep)

load_line <- load_line %>% left_join(wes_rep, by = "cell_line")

# Per-patient table (one row / patient), with within-family load range shown.
load_patient <- load_line %>%
  group_by(patient_id, family, n_family_maf_lines) %>%
  summarise(subtype        = subtype[is_wes_patient_rep][1],
            rep_line        = cell_line[is_wes_patient_rep][1],
            n_coding_rep    = n_coding[is_wes_patient_rep][1],
            indel_frac_rep  = indel_frac[is_wes_patient_rep][1],
            tstv_rep        = tstv[is_wes_patient_rep][1],
            family_load_min = min(n_coding), family_load_max = max(n_coding),
            is_hypermutator = any(is_hypermutator),
            .groups = "drop") %>%
  mutate(fold_over_median = n_coding_rep / med_load) %>%
  arrange(desc(n_coding_rep))

cat("\n=== PER-LINE mutation load (ranked) ===\n")
print(as.data.frame(load_line %>%
  select(load_rank, cell_line, subtype, n_coding, n_snv, n_indel, indel_frac,
         tstv, fold_over_median, robust_z, is_hypermutator) %>%
  mutate(across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
cat(sprintf("\nMedian coding-nonsyn load = %.0f | MAD = %.1f | top = %s (%d) | 2nd = %s (%d) | fold(1st/2nd)=%.2f\n",
            med_load, mad_load, load_line$cell_line[1], load_line$n_coding[1],
            load_line$cell_line[2], load_line$n_coding[2],
            load_line$n_coding[1] / load_line$n_coding[2]))
cat("\n=== PER-PATIENT mutation load (collapsed; n patients =", nrow(load_patient), ") ===\n")
print(as.data.frame(load_patient %>%
  mutate(across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

# =============================================================================
# 4. MMR / PROOFREADING GENE STATUS (coding non-syn hits in the filtered MAF)
# =============================================================================
mmr_genes <- c(MMR_ENZYME, MMR_INDIRECT, MMR_OTHER, PROOF)
mmr_hits <- maf %>%
  filter(Hugo_Symbol %in% mmr_genes) %>%
  mutate(class = ifelse(Variant_Classification %in% TRUNCATING, "truncating",
                 ifelse(Variant_Classification == "Missense_Mutation", "missense", "other")),
         # First integer in the reconstructed HGVSp string. [review revision] The
         # previous lookbehind required a letter immediately before the digits, so it
         # returned NA for the position-only forms 07 now emits when the affected
         # residue identity is not recoverable (e.g. "p.2_3fs") — which silently
         # dropped TOV21G's EPCAM codon and any POLE indel from the domain check.
         codon = suppressWarnings(as.integer(str_extract(HGVSp_Short, "[0-9]+")))) %>%
  select(cell_line, subtype, Hugo_Symbol, Variant_Classification, class, HGVSp_Short,
         codon, Consequence, IMPACT, CLIN_SIG, pop_af_max, vaf, germline_like_vaf)

cat("\n=== MMR / proofreading gene hits (panel-wide, coding non-syn) ===\n")
print(as.data.frame(mmr_hits), row.names = FALSE)
if (nrow(mmr_hits) == 0) cat("  (none)\n")

# POLE exonuclease-domain check (canonical ultramutator hotspots ~codon 268-471)
pole_exo <- mmr_hits %>%
  filter(Hugo_Symbol == "POLE", !is.na(codon), codon >= POLE_EXO[1], codon <= POLE_EXO[2])
cat(sprintf("\nPOLE exonuclease-domain (codon %d-%d) damaging hits: %d\n",
            POLE_EXO[1], POLE_EXO[2], nrow(pole_exo)))

# Per-line MMR/proofreading summary
# *** EVERY MMR-PANEL HIT IS SURFACED EXPLICITLY (review revision) ***
#   The manuscript states "no MMR-enzyme coding mutation is present" for TOV21G,
#   yet this panel records an EPCAM truncating call for exactly that line — the only
#   truncating hit anywhere in the panel. Whether real or artefactual, the pipeline
#   cannot surface a hit and have the text deny it. The summary below therefore
#   carries (a) `mmr_panel_hit`, true for ANY hit in the whole panel, (b) the literal
#   variant string per gene, and (c) `mmr_note`, a per-line caveat sentence that
#   travels with the row into wes_msi_mmr.csv and Table S1's source, so the text is
#   forced to address it.
mmr_line <- load_line %>%
  transmute(cell_line, subtype, load_rank, n_coding, n_indel, indel_frac, tstv,
            is_hypermutator) %>%
  left_join(
    mmr_hits %>% group_by(cell_line) %>%
      summarise(mmr_enzyme_hit = any(Hugo_Symbol %in% MMR_ENZYME),   # MLH1/MSH2/MSH6/PMS2
                epcam_hit      = any(Hugo_Symbol == "EPCAM"),         # indirect (MSH2 silencing)
                mmr_genes_hit = paste(sort(unique(paste0(Hugo_Symbol, "(", class, ")"))), collapse = "; "),
                mmr_variants  = paste(paste0(Hugo_Symbol, " ", HGVSp_Short, " [",
                                             Variant_Classification, "; VAF ",
                                             ifelse(is.na(vaf), "NA (indel)", sprintf("%.3f", vaf)),
                                             "; ClinVar ", ifelse(is.na(CLIN_SIG), "-", CLIN_SIG), "]"),
                                      collapse = " ; "),
                mmr_any_truncating = any(class == "truncating"),
                pole_hit      = any(Hugo_Symbol == "POLE"),
                pold1_hit     = any(Hugo_Symbol == "POLD1"),
                .groups = "drop"),
    by = "cell_line") %>%
  mutate(across(c(mmr_enzyme_hit, epcam_hit, pole_hit, pold1_hit, mmr_any_truncating),
                ~replace_na(., FALSE)),
         mmr_genes_hit = replace_na(mmr_genes_hit, ""),
         mmr_variants  = replace_na(mmr_variants, ""),
         mmr_panel_hit = nzchar(mmr_genes_hit),
         has_pole_exo  = cell_line %in% pole_exo$cell_line)

# Per-line MMR caveat sentence. TOV21G's EPCAM call gets the full statement,
# including the mechanism mismatch: the panel flags EPCAM because 3'-END deletions
# silence MSH2 in cis, whereas THIS variant is a 5' frameshift at codon 2-3, which
# truncates EPCAM itself and does NOT invoke the indirect-MSH2 mechanism.
epcam_detail <- mmr_hits %>% filter(Hugo_Symbol == "EPCAM") %>%
  transmute(cell_line, txt = sprintf("%s (%s)", HGVSp_Short, Variant_Classification),
            codon5 = !is.na(codon) & codon <= 50)
mmr_line <- mmr_line %>%
  left_join(epcam_detail, by = "cell_line") %>%
  mutate(mmr_note = dplyr::case_when(
    epcam_hit ~ paste0(
      "MMR-PANEL HIT — REPORT, DO NOT OMIT: EPCAM ", txt,
      ". EPCAM is on this panel as an INDIRECT MMR mechanism (3'-end deletions silence MSH2 in cis); ",
      if_else(replace_na(codon5, FALSE),
              "this call is a 5' frameshift near the start codon, which truncates EPCAM itself and does NOT invoke that mechanism. ",
              "the position is consistent with the 3' mechanism. "),
      "TUMOUR-ONLY WES, no matched normal => GERMLINE CANNOT BE EXCLUDED (EPCAM/Lynch germline deletions exist). ",
      "The variant is novel in dbSNP, carries no VAF support (indel representations do not join the VCF), ",
      "and sits in the panel's one hypermutator, where passenger risk is elevated. ",
      "Requires MMR IHC / MSI-PCR and germline testing before any mechanistic claim."),
    mmr_enzyme_hit ~ paste0(
      "MMR-PANEL HIT — MMR-ENZYME (Lynch core) coding variant: ", mmr_variants,
      ". Tumour-only WES => germline cannot be excluded; requires IHC/MSI-PCR."),
    mmr_panel_hit ~ paste0(
      "MMR-PANEL HIT — secondary/proofreading gene: ", mmr_variants,
      ". Tumour-only WES => germline cannot be excluded; not an established MMR-loss mechanism on its own."),
    TRUE ~ "no retained coding candidate in the MMR / proofreading panel")) %>%
  select(-txt, -codon5)

cat("\n*** MMR-PANEL HITS — these MUST appear in the manuscript text ***\n")
if (any(mmr_line$mmr_panel_hit)) {
  print(as.data.frame(mmr_line %>% filter(mmr_panel_hit) %>%
    select(cell_line, subtype, n_coding, is_hypermutator, mmr_enzyme_hit, epcam_hit,
           mmr_any_truncating, mmr_variants)), row.names = FALSE)
  for (i in which(mmr_line$mmr_panel_hit))
    cat(sprintf("  - %s: %s\n", mmr_line$cell_line[i], mmr_line$mmr_note[i]))
} else cat("  none in the 22-line panel\n")

# =============================================================================
# 5. MUTATIONAL CONTEXT — SBS-96 profile (GRCh38) + COSMIC v3.2 cosine
#    Substrate: exome-wide PASS somatic-candidate SNVs from the RAW MAFs (a
#    superset of the coding-nonsyn filtered MAF), re-derived with the SAME
#    PASS + population-AF filter as script 07 — coding-only SNVs are too few /
#    selection-biased for a stable 96-context. Read-only on the archived source.
# =============================================================================
SBS_OK <- all(c("MutationalPatterns","BSgenome.Hsapiens.UCSC.hg38","GenomicRanges") %in%
              rownames(installed.packages()))
sbs_line <- NULL; mut_mat <- NULL; sbs_note <- ""

if (SBS_OK) {
  suppressPackageStartupMessages({
    library(GenomicRanges); library(MutationalPatterns); library(BSgenome.Hsapiens.UCSC.hg38)
  })
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  std_chr    <- paste0("chr", c(1:22, "X", "Y"))
  MUT_DIR    <- file.path(DATA, "wes - old", "mutect2")
  canon_key  <- function(x) toupper(gsub("[-_]", "", sub("_[Pp][0-9]+$", "", x)))

  # map raw MAF dirs -> cell_line (same logic as 07), keep the 22 filtered lines
  source("scripts/lib/wes_maf_inputs.R")
  maf_map <- wes_maf_inputs(PROJ, DATA, OUT) %>%
    transmute(path = source_maf, cell_line) %>%
    left_join(load_line %>% select(cell_line, subtype), by = "cell_line")
  stopifnot("MAFs did not map to all filtered models" = setequal(maf_map$cell_line, load_line$cell_line))

  read_pass <- function(path, cell_line) {
    r <- fread(path, sep = "\t", skip = "Hugo_Symbol", quote = "",
               na.strings = c("", "NA", "."), showProgress = FALSE)
    afc <- intersect(c("gnomADe_AF","AF","AA_AF","EA_AF"), names(r))
    n0  <- function(v) { x <- suppressWarnings(as.numeric(v)); ifelse(is.na(x), 0, x) }
    r[, pop_af_max := do.call(pmax, lapply(.SD, n0)), .SDcols = afc]
    r <- r[FILTER == "PASS" & pop_af_max <= POP_AF_MAX]                 # exome-wide candidates
    snv <- r[Variant_Type == "SNP" &
               Reference_Allele  %in% c("A","C","G","T") &
               Tumor_Seq_Allele2 %in% c("A","C","G","T") &
               Chromosome %in% std_chr]
    gr <- GRanges(seqnames = snv$Chromosome,
                  ranges  = IRanges(start = snv$Start_Position, width = 1),
                  REF = Biostrings::DNAStringSet(snv$Reference_Allele),
                  ALT = Biostrings::DNAStringSetList(as.list(snv$Tumor_Seq_Allele2)))
    GenomeInfoDb::seqlevels(gr) <- std_chr
    GenomeInfoDb::genome(gr)    <- "hg38"
    list(gr = gr,
         ew_snv   = sum(r$Variant_Type == "SNP"),
         ew_indel = sum(r$Variant_Type %in% c("INS","DEL")))
  }
  message("Reading raw MAFs for exome-wide PASS SNVs (SBS substrate)...")
  parsed <- Map(read_pass, maf_map$path, maf_map$cell_line)
  names(parsed) <- maf_map$cell_line
  grl <- GRangesList(lapply(parsed, `[[`, "gr"))
  ew  <- tibble(cell_line = names(parsed),
                ew_pass_snv   = vapply(parsed, `[[`, numeric(1), "ew_snv"),
                ew_pass_indel = vapply(parsed, `[[`, numeric(1), "ew_indel")) %>%
         mutate(ew_pass_indel_frac = ew_pass_indel / (ew_pass_snv + ew_pass_indel))

  # 96-context matrix + COSMIC cosine similarity (reference release pinned above)
  mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)      # 96 x 22
  cosmic  <- tryCatch(get_known_signatures(muttype = "snv", source = COSMIC_SOURCE,
                                           genome = COSMIC_GENOME),
                      error = function(e) { message("COSMIC fetch failed: ", conditionMessage(e)); NULL })
  if (!is.null(cosmic))
    message(sprintf("COSMIC reference: %s / %s (%d signatures) | MutationalPatterns %s",
                    COSMIC_SOURCE, COSMIC_GENOME, ncol(cosmic),
                    as.character(packageVersion("MutationalPatterns"))))

  six_frac <- function(v) {
    ty <- sub("^.\\[(.*)\\].$", "\\1", rownames(mut_mat))
    tapply(v, ty, sum) / sum(v)
  }
  sbs_line <- tibble(cell_line = colnames(mut_mat),
                     n_snv_used = colSums(mut_mat)) %>%
    left_join(ew, by = "cell_line")
  six <- t(apply(mut_mat, 2, six_frac)); colnames(six) <- paste0("frac_", colnames(six))
  sbs_line <- bind_cols(sbs_line, as_tibble(six))

  if (!is.null(cosmic)) {
    cs <- cos_sim_matrix(mut_mat, cosmic)                              # lines x signatures
    grp_max <- function(sigs) {
      sigs <- intersect(sigs, colnames(cs))
      if (!length(sigs)) return(rep(NA_real_, nrow(cs)))
      apply(cs[, sigs, drop = FALSE], 1, max)
    }
    grp_which <- function(sigs) {                    # WHICH signature carries the group max
      sigs <- intersect(sigs, colnames(cs))
      if (!length(sigs)) return(rep(NA_character_, nrow(cs)))
      setNames(sigs[apply(cs[, sigs, drop = FALSE], 1, which.max)], rownames(cs))
    }
    best_sig  <- setNames(colnames(cs)[apply(cs, 1, which.max)], rownames(cs))
    top3      <- apply(cs, 1, function(r) paste(names(sort(r, decreasing = TRUE))[1:3],
                                                collapse = ";"))  # named by rownames(cs)
    sbs_line <- sbs_line %>%
      mutate(cos_mmr_d  = grp_max(MMR_D_SIGS)[cell_line],
             cos_pole   = grp_max(POLE_SIGS)[cell_line],
             cos_clock  = grp_max(CLOCK_SIGS)[cell_line],
             # WHICH signature each group max came from — "cosine 0.88" is otherwise
             # unattributable, and the text/figure disagreed on SBS20 vs SBS6/44/15.
             cos_mmr_d_sig = grp_which(MMR_D_SIGS)[cell_line],
             cos_pole_sig  = grp_which(POLE_SIGS)[cell_line],
             cos_clock_sig = grp_which(CLOCK_SIGS)[cell_line],
             best_cosmic = best_sig[cell_line],
             top3_cosmic = top3[cell_line],
             # MARGINS (review revision): the discriminating quantity, reported per
             # line, replacing an unjustified absolute 0.75 cut on the group maximum.
             cos_margin_mmrd_vs_clock = cos_mmr_d - cos_clock,
             cos_margin_mmrd_vs_pole  = cos_mmr_d - cos_pole,
             cos_margin_pole_vs_clock = cos_pole  - cos_clock,
             cos_mmr_d_ge_075 = cos_mmr_d >= SBS_ABS_LEGACY,   # legacy screen, reported only
             sbs_call = case_when(
               cos_margin_mmrd_vs_clock >= SBS_MARGIN_MIN & cos_mmr_d > cos_pole ~ "MMR-d/MSI-like",
               cos_margin_pole_vs_clock >= SBS_MARGIN_MIN & cos_pole  > cos_mmr_d ~ "POLE-like",
               TRUE                                                              ~ "clock/indeterminate"),
             cosmic_source = COSMIC_SOURCE, cosmic_genome = COSMIC_GENOME,
             mutationalpatterns_version = as.character(packageVersion("MutationalPatterns")))
    cat(sprintf("\n=== SBS-96 cosine vs %s (%s) — group maxima AND margins, per line ===\n",
                COSMIC_SOURCE, COSMIC_GENOME))
    cat("    n_snv_used = exome-wide PASS & popAF<=0.001 SNVs entering that line's spectrum.\n")
    cat(sprintf("    CALL RULE: cos(group) - cos(clock) >= %.2f (panel-relative). Legacy absolute\n    screen cos_mmr_d >= %.2f is reported (cos_mmr_d_ge_075) but no longer decisive.\n",
                SBS_MARGIN_MIN, SBS_ABS_LEGACY))
    print(as.data.frame(sbs_line %>%
      select(cell_line, n_snv_used, ew_pass_indel_frac, cos_mmr_d, cos_mmr_d_sig,
             cos_clock, cos_pole, cos_margin_mmrd_vs_clock, cos_mmr_d_ge_075,
             best_cosmic, sbs_call) %>%
      arrange(desc(cos_mmr_d)) %>% mutate(across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
    cat(sprintf("  => MMR-d-vs-clock margin: %s %+.3f ; the other %d lines %+.3f to %+.3f (all negative)\n",
                sbs_line$cell_line[which.max(sbs_line$cos_margin_mmrd_vs_clock)],
                max(sbs_line$cos_margin_mmrd_vs_clock),
                nrow(sbs_line) - 1L,
                min(sort(sbs_line$cos_margin_mmrd_vs_clock, decreasing = TRUE)[-1]),
                max(sort(sbs_line$cos_margin_mmrd_vs_clock, decreasing = TRUE)[-1])))
    cat("  CAVEAT: SBS1/5/6/15/44 are mutually similar, so 0.88 vs a 0.65 clock group is a\n",
        "  MODEST separation in cosine space. Cosine screening is not a fit — see\n",
        "  22_wes_signature_refit.R for a bootstrapped fit_to_signatures refit.\n")

    # FULL cosine profile (not just the three group maxima) — review revision.
    cos_out <- as.data.frame(cs) %>% rownames_to_column("cell_line") %>%
      pivot_longer(-cell_line, names_to = "signature", values_to = "cosine") %>%
      mutate(group = case_when(signature %in% MMR_D_SIGS ~ "MMR-d / MSI",
                               signature %in% POLE_SIGS  ~ "POLE proofreading",
                               signature %in% CLOCK_SIGS ~ "clock-like (age)",
                               TRUE                       ~ "other"),
             cosmic_source = COSMIC_SOURCE, cosmic_genome = COSMIC_GENOME) %>%
      left_join(sbs_line %>% select(cell_line, n_snv_used), by = "cell_line") %>%
      arrange(cell_line, desc(cosine))
    readr::write_csv(cos_out, file.path(OUT, "wes_sbs_cosine.csv"))
    message("Wrote output/wes_sbs_cosine.csv (", nrow(cos_out), " line x signature cosines)")
  }
  # store the 96-context matrix (rows = contexts, cols = lines)
  sbs96 <- as.data.frame(mut_mat) %>% rownames_to_column("context")
  readr::write_csv(sbs96, file.path(OUT, "wes_sbs_context.csv"))
  message("Wrote output/wes_sbs_context.csv (96 x ", ncol(mut_mat), ")")
  sbs_note <- sprintf("SBS-96 computed on GRCh38 for %d lines (exome-wide PASS SNVs).",
                      ncol(mut_mat))
} else {
  sbs_note <- "SBS-96 SKIPPED (MutationalPatterns / BSgenome.Hsapiens.UCSC.hg38 not installed)."
  message(sbs_note)
}

# =============================================================================
# 6. INTEGRATIVE MSI / MMR ASSESSMENT TABLE
# =============================================================================
have_sbs <- !is.null(sbs_line) && "cos_mmr_d" %in% names(sbs_line)

msi_tbl <- mmr_line %>%
  left_join(
    if (have_sbs)
      sbs_line %>% select(cell_line, ew_pass_snv, ew_pass_indel, ew_pass_indel_frac,
                          n_snv_used, cos_mmr_d, cos_mmr_d_sig, cos_pole, cos_pole_sig,
                          cos_clock, cos_clock_sig, cos_margin_mmrd_vs_clock,
                          cos_margin_mmrd_vs_pole, cos_mmr_d_ge_075,
                          best_cosmic, top3_cosmic, sbs_call,
                          cosmic_source, cosmic_genome, mutationalpatterns_version)
    else tibble(cell_line = mmr_line$cell_line),
    by = "cell_line")

# Precompute SBS confirmation vectors (vectorised; scalar-safe when no SBS).
# [review revision] mmrd_confirm now uses the SAME margin rule as sbs_call, so the
# integrative call and the spectrum call cannot diverge.
pole_confirm <- if (have_sbs) msi_tbl$sbs_call == "POLE-like"      else rep(TRUE, nrow(msi_tbl))
mmrd_confirm <- if (have_sbs) msi_tbl$sbs_call == "MMR-d/MSI-like" else rep(TRUE, nrow(msi_tbl))

msi_tbl <- msi_tbl %>%
  mutate(
    # Integrative qualitative call — tumor-only => CANDIDATE, not a definitive
    # clinical MSI/POLE call. Requires converging load + indel + (SBS) evidence.
    msi_mmr_call = dplyr::case_when(
      has_pole_exo & pole_confirm                                 ~ "candidate POLE-ultramutator",
      is_hypermutator & indel_frac >= 0.20 & mmrd_confirm         ~ "candidate MMR-d / MSI-high",
      is_hypermutator                                             ~ "hypermutated, mechanism unresolved",
      TRUE                                                        ~ "no hypermutation signal")) %>%
  arrange(load_rank)

cat("\n=== INTEGRATIVE MSI/MMR assessment (per line) ===\n")
print(as.data.frame(msi_tbl %>%
  select(load_rank, cell_line, subtype, n_coding, indel_frac, tstv, mmr_panel_hit,
         mmr_enzyme_hit, epcam_hit, has_pole_exo,
         any_of(c("cos_mmr_d","cos_margin_mmrd_vs_clock","sbs_call")), msi_mmr_call) %>%
  mutate(across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

# =============================================================================
# 7. WRITE OUTPUTS
# =============================================================================
# Per-line load (per-patient is derivable: filter is_wes_patient_rep == TRUE).
load_out <- load_line %>%
  select(load_rank, cell_line, subtype, patient_id, family, n_family_maf_lines,
         is_wes_patient_rep, n_coding, n_snv, n_ins, n_del, n_indel, indel_frac,
         n_mnv, n_ts, n_tv, tstv, ts_frac, fold_over_median, robust_z, is_hypermutator)
readr::write_csv(load_out, file.path(OUT, "wes_mutation_load.csv"))
readr::write_csv(msi_tbl,  file.path(OUT, "wes_msi_mmr.csv"))
message("Wrote output/wes_mutation_load.csv (", nrow(load_out), " lines) and output/wes_msi_mmr.csv")

# =============================================================================
# 8. FIGURE — f_wes_hypermutation (4 panels)
#    A load bar (indel fraction) | B MMR-gene status | C SBS-96 | D COSMIC cosine
# =============================================================================
lvl <- load_line$cell_line          # ranked, highest load first

# Panel A — ranked coding-load bar, fill = indel fraction
pA <- load_line %>%
  mutate(cell_line = factor(cell_line, levels = rev(lvl)),
         lab = scales::comma(n_coding)) %>%
  ggplot(aes(n_coding, cell_line, fill = indel_frac)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.2) +
  geom_vline(xintercept = med_load, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_text(aes(label = lab), hjust = -0.12, size = 2.1, colour = "grey15") +
  geom_text(data = ~dplyr::filter(., is_hypermutator),
            aes(x = n_coding * 0.5, label = "hypermutator"),
            colour = "white", fontface = "bold", size = 2.6) +
  scale_fill_viridis_c(option = "magma", end = 0.93, name = "indel\nfraction") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "A  Coding mutation load (tumor-only; relative)",
       subtitle = sprintf("TOV21G %s vs panel median %.0f (%.1f×); dashed = median",
                          scales::comma(load_line$n_coding[1]), med_load, load_line$fold_over_median[1]),
       x = "coding non-synonymous variants", y = NULL) +
  theme_lab()

# Panel B — MMR / proofreading gene status (top-load lines x genes)
topN   <- head(lvl, 10)
tile <- tidyr::expand_grid(cell_line = topN, Hugo_Symbol = mmr_genes) %>%
  left_join(mmr_hits %>% group_by(cell_line, Hugo_Symbol) %>%
              summarise(class = dplyr::case_when(any(class == "truncating") ~ "truncating",
                                                 any(class == "missense")  ~ "missense",
                                                 TRUE ~ "other"), .groups = "drop"),
            by = c("cell_line","Hugo_Symbol")) %>%
  mutate(class = tidyr::replace_na(class, "none"),
         cell_line = factor(cell_line, levels = rev(topN)),
         Hugo_Symbol = factor(Hugo_Symbol, levels = mmr_genes))
pB <- ggplot(tile, aes(Hugo_Symbol, cell_line, fill = class)) +
  geom_tile(colour = "grey80", linewidth = 0.4) +
  scale_fill_manual(values = c(truncating = COOK_RUST, missense = COOK_NAVY,
                               other = "grey55", none = "grey95"), name = "coding\nalteration") +
  labs(title = "B  MMR / proofreading gene status",
       subtitle = paste0("top-10 lines by load; TOV21G carries an EPCAM truncating call ",
                         "(the only truncating MMR-panel hit) — no MMR-enzyme or POLE hit.\n",
                         "Tumour-only WES: germline cannot be excluded (see mmr_note in wes_msi_mmr.csv)"),
       x = NULL, y = NULL) +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

# Panels C/D — SBS (only if computed)
if (have_sbs) {
  cbf6 <- c("#56B4E9","#000000","#D55E00","#999999","#009E73","#CC79A7")  # CVD-safe 6-class
  pC <- plot_96_profile(mut_mat[, "TOV21G", drop = FALSE], condensed = TRUE, colors = cbf6) +
    labs(title = "C  TOV21G SBS-96 spectrum",
         subtitle = sprintf("%s exome-wide PASS SNVs (GRCh38)",
                            scales::comma(sbs_line$n_snv_used[sbs_line$cell_line == "TOV21G"]))) +
    theme_lab() +
    theme(strip.text = element_text(size = rel(0.8)),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none")

  diag_sigs <- intersect(c(MMR_D_SIGS, POLE_SIGS, CLOCK_SIGS), colnames(cs))
  grp_of <- function(s) dplyr::case_when(s %in% MMR_D_SIGS ~ "MMR-d / MSI",
                                         s %in% POLE_SIGS  ~ "POLE",
                                         TRUE ~ "clock (age)")
  dD <- tibble(sig = diag_sigs, cos = cs["TOV21G", diag_sigs], grp = grp_of(diag_sigs))
  pD <- ggplot(dD, aes(cos, reorder(sig, cos), fill = grp)) +
    geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", cos)), hjust = -0.15, size = 2.2, colour = "grey15") +
    scale_fill_manual(values = c("MMR-d / MSI" = COOK_RUST, "POLE" = COOK_NAVY,
                                 "clock (age)" = "grey60"), name = NULL) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.14))) +
    labs(title = "D  TOV21G spectrum vs COSMIC v3.2",
         subtitle = "cosine similarity to diagnostic SBS signatures",
         x = "cosine similarity", y = NULL) +
    theme_lab() + theme(legend.position = "top")
  fig <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 0.95))
  fig_h <- 9.6
} else {
  fig <- (pA | pB); fig_h <- 5.2
}

fig <- fig + plot_annotation(
  title = "WES hypermutation, MMR/MSI status, and mutational context",
  subtitle = paste("Tumor-only WES (no matched normal): load is a RELATIVE within-panel metric;",
                    "MSI/MMR status is QUALITATIVE (converging proxies), not a clinical assay.", sbs_note),
  theme = theme(plot.title = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(size = 8.4, colour = "grey30")))

ggsave(file.path(FIGS, "f_wes_hypermutation.pdf"), fig, width = 12.5, height = fig_h)
ggsave(file.path(PROJ, "reports", "assets", "f_wes_hypermutation.png"), fig,
       width = 12.5, height = fig_h, dpi = 200)
message("Wrote figs/f_wes_hypermutation.pdf and reports/assets/f_wes_hypermutation.png")

message("\n16_wes_signatures_msi.R complete.")
message("Verdict scaffold: TOV21G -> ",
        msi_tbl$msi_mmr_call[msi_tbl$cell_line == "TOV21G"])
message("MMR-panel note for the text: ", msi_tbl$mmr_note[msi_tbl$cell_line == "TOV21G"])
message("Cosine screening is NOT a signature fit — run 22_wes_signature_refit.R for ",
        "bootstrapped exposures, reconstruction cosine and the pinned COSMIC release.")
