# =============================================================================
# Script: 07_wes_mutations.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Re-analyse the WES SNV/indel calls (Mutect2, "old" module) for the
#          23 generated WES lines. RE-FILTER the tumor-only calls and report a
#          CANONICAL-DRIVER landscape only, with TP53 in HGS as a positive
#          control. Phase 3 (WES), step 1 of 3.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# =============================================================================
#
# *** CRITICAL ASSAY CAVEAT — TUMOR-ONLY WES, NO MATCHED NORMAL ***
#   These Mutect2 calls were generated WITHOUT a matched normal (the MAF's
#   n_depth/n_ref_count/n_alt_count columns are empty for every variant).
#   Tumor-only somatic calling is germline-contaminated: ~224 private germline
#   variants/sample and ~50-70% FDR after standard filtering, and germline is
#   HARDEST to separate in ~100%-pure samples like cell lines (Halperin 2017,
#   10.1186/s12920-017-0296-8; Little 2021 UNMASC, 10.1093/narcan/zcab040).
#   => Per-gene "mutation frequencies" from these data are NOT reliable somatic
#      rates. We therefore (a) apply the Mutect2 FILTER column that the archived
#      analysis evidently ignored [it is the mechanism behind the implausible
#      "ATM 100% / ATR 75% / BRCA2 majority in HGSC" claims], (b) add a
#      population-AF germline filter, and (c) restrict the reported landscape to
#      CANONICAL DRIVERS + report TP53 in HGS as a positive control.
#   The residual FDR cannot be removed without a matched normal; treat all retained
#   calls as candidates; canonical TP53 alterations serve as a positive control.
#
# WHY re-filtering works here (verified on disk):
#   Each MAF holds ALL Mutect2 records, not just PASS (e.g. OV2295: 25,914 rows,
#   only 493 PASS). Mutect2 was run WITH a germline resource + panel-of-normals:
#   the FILTER field already carries "germline", "panel_of_normals",
#   "common_variant" flags. Applying FILTER==PASS removes the germline-driven
#   artefact calls (ATR: 10 calls -> 0 PASS; PTEN 8 -> 0; BRCA2 3 -> 0), while
#   the real somatic TP53 hit (private, gnomAD AF 1.4e-6) is retained.
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(maftools)
  library(ComplexHeatmap); library(circlize); library(grid)
})
set.seed(SEED)

MUT_DIR <- file.path(DATA, "wes - old", "mutect2")
stopifnot("Mutect2 dir not found" = dir.exists(MUT_DIR))

# Canonical drivers for ovarian carcinoma (the reported landscape). ------------
# Rationale: TP53 (HGSC ~universal), homologous-recombination (BRCA1/2, CDK12),
# RB1/CCNE axis, PI3K/RAS (PTEN/PIK3CA/KRAS/NRAS/BRAF), NF1, SWI/SNF
# (ARID1A/ARID1B/SMARCA4/SMARCA2 — clear-cell/endometrioid/dedifferentiated),
# CTNNB1, PPP2R1A, CDKN2A, ERBB2. This is the set the descriptor reports.
CANONICAL_DRIVERS <- c("TP53","BRCA1","BRCA2","RB1","NF1","PTEN","PIK3CA",
                       "KRAS","NRAS","BRAF","CTNNB1","ARID1A","ARID1B",
                       "PPP2R1A","CDK12","CDKN2A","SMARCA4","SMARCA2","ERBB2")
# Genes to check the artefact is gone (were implausibly frequent pre-filter):
ARTEFACT_WATCH <- c("ATM","ATR","BRCA2","FANCM","POLE","POLQ")
POP_AF_MAX <- 0.001   # common-germline threshold (gnomAD/1000G/ESP)

# maftools' non-synonymous set (what counts as a driver "hit"):
NONSYN <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
            "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site",
            "Nonstop_Mutation","Translation_Start_Site")

# Canonical-key normaliser: strip passage suffix, drop separators, upcase.
canon_key <- function(x) toupper(gsub("[-_]", "",
                          sub("_[Pp][0-9]+$", "", x)))

# 1. Sample sheet -> the 23 WES-mutation lines + subtype ----------------------
ss  <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE)
wes <- ss %>%
  filter(wes_mut == "Y", provenance == "generated") %>%
  transmute(cell_line, subtype,
            key = canon_key(cell_line))
stopifnot("Expected 23 WES-mutation lines" = nrow(wes) == 23)
message(sprintf("WES-mutation lines in sample sheet: %d", nrow(wes)))

# Resolve archived MAFs and annotation-only recoveries through one inventory.
source("scripts/lib/wes_maf_inputs.R")
maf_inputs <- wes_maf_inputs(PROJ, DATA, OUT)
readr::write_csv(maf_inputs %>% mutate(across(c(source_maf, source_vcf),
  ~ sub(paste0("^", PROJ, "/"), "", .x))), file.path(OUT, "wes_input_manifest.csv"))
maf_map <- maf_inputs %>% transmute(path = source_maf, source_vcf, cell_line) %>%
  left_join(wes %>% select(cell_line, subtype), by = "cell_line")
stopifnot("Every declared WES model must have a MAF" = nrow(maf_map) == nrow(wes),
          !anyNA(maf_map$subtype))
message(sprintf("Usable MAFs: %d/%d; recovered from archived annotated VCF: %s",
  nrow(maf_map), nrow(wes), paste(maf_inputs$cell_line[grepl("recovered", maf_inputs$source_kind)], collapse = ", ")))

# 2. Load all MAFs (raw; keep every record so we can show the filter effect) --
read_one_maf <- function(path, cell_line, subtype) {
  dt <- fread(path, sep = "\t", skip = "Hugo_Symbol", quote = "",
              na.strings = c("", "NA", "."), showProgress = FALSE)
  dt[, `:=`(cell_line = cell_line, subtype = subtype)]
  dt[]
}
raw <- pmap(list(maf_map$path, maf_map$cell_line, maf_map$subtype), read_one_maf) %>%
  rbindlist(fill = TRUE)
message(sprintf("Loaded %d raw Mutect2 records across %d lines",
                nrow(raw), length(unique(raw$cell_line))))

# 3. Derive population AF, VAF (from VCF), reconstruct protein change ----------
# Population AF across the annotation columns present in the MAF; NA -> 0.
af_cols <- c("gnomADe_AF","AF","AA_AF","EA_AF")     # gnomAD exome, 1000G, ESP
stopifnot("Expected population-AF columns absent" = all(af_cols %in% names(raw)))
# NA -> 0 (variant absent from that db => not common there). Rare multiallelic
# "a,b" strings coerce to NA->0 here, but Mutect2's own common_variant/
# multiallelic FILTER flags already remove those (belt-and-suspenders w/ PASS).
num_or0 <- function(v) { x <- suppressWarnings(as.numeric(v)); ifelse(is.na(x), 0, x) }
raw[, pop_af_max := do.call(pmax, lapply(.SD, num_or0)), .SDcols = af_cols]

# HGVSp_Short is EMPTY in every one of these MAFs (the column is all-NA), so the
# protein change is RECONSTRUCTED from VEP's Amino_acids + Protein_position.
# *** REVISED (review): the previous reconstruction emitted non-HGVS strings ***
#   It pasted "p.<ref><pos><alt>" verbatim, which yields identifiers that will not
#   match ClinVar/COSMIC: "p.L639X" for a frameshift (VEP writes the unknown new
#   reading frame as "X"), "p.-156-157P" for an insertion, "p.KAL623-625X" for a
#   multi-residue frameshift. The manuscript quoted one of them (SMARCA4 p.L639X)
#   as a variant identifier. We now emit standard short notation wherever these two
#   fields determine it, and FLAG the rest:
#     missense/nonsense/nonstop/start-loss   Q/*   192       -> p.Q192*        canonical
#     frameshift (VEP alt carries "X")       L/X   639       -> p.L639fs       canonical
#     in-frame deletion (VEP alt "-")        CP/-  1255-1256 -> p.C1255_P1256del canonical
#     in-frame insertion (VEP ref "-")       -/P   156-157   -> p.156_157insP  NOT canonical
#         (HGVS requires the flanking residue identities; Amino_acids does not carry them)
#     multi-residue substitution             KA/NS 512-513   -> p.K512_A513delinsNS
#         (a valid FORM but not the minimally-reduced HGVS description) NOT canonical
#   Two provenance columns travel with the table: `hgvsp_reconstructed` (TRUE for
#   every row — none of these strings came from the MAF) and `hgvsp_canonical`
#   (which strings are safe to quote as variant identifiers).
AA1 <- "^[ACDEFGHIKLMNPQRSTVWY*]$"
mk_hgvsp <- function(aa, pos, vclass) {
  n   <- length(aa)
  aa  <- ifelse(is.na(aa), "", aa)
  pp  <- sub("/.*$", "", ifelse(is.na(pos), "", pos))   # keep protein position only
  p1  <- suppressWarnings(as.integer(sub("-.*$", "", pp)))          # first residue
  p2  <- suppressWarnings(as.integer(sub("^[0-9]*-", "", pp)))      # last residue
  p2  <- ifelse(is.na(p2), p1, p2)
  ref <- ifelse(grepl("/", aa), sub("/.*$", "", aa), aa)
  alt <- ifelse(grepl("/", aa), sub("^.*/", "", aa), aa)            # synon: single letter
  r1  <- substr(ref, 1, 1); rn <- substr(ref, nchar(ref), nchar(ref))
  hg  <- rep(NA_character_, n); canon <- rep(NA, n)
  ok  <- aa != "" & !is.na(p1)
  # frameshift first: VEP encodes the new (unknown-length) frame as "X". Reference
  # allele "-" means an insertion-frameshift, so the first CHANGED residue is not
  # identified by these fields -> position-only string, marked non-canonical.
  fs  <- ok & (grepl("X", alt, fixed = TRUE) |
               vclass %in% c("Frame_Shift_Del", "Frame_Shift_Ins"))
  fsk <- fs & grepl(AA1, r1)                      # first affected residue known
  hg[fsk] <- paste0("p.", r1[fsk], p1[fsk], "fs"); canon[fsk] <- TRUE
  fsu <- fs & !fsk
  hg[fsu] <- ifelse(p1[fsu] == p2[fsu], paste0("p.", p1[fsu], "fs"),
                    paste0("p.", p1[fsu], "_", p2[fsu], "fs")); canon[fsu] <- FALSE
  # single-residue substitution: missense, nonsense (*), nonstop, translation start
  s1 <- ok & !fs & grepl(AA1, ref) & grepl(AA1, alt)
  hg[s1] <- paste0("p.", ref[s1], p1[s1], alt[s1]); canon[s1] <- TRUE
  # in-frame deletion (VEP alt "-")
  dl <- ok & !fs & !s1 & alt == "-" & ref != "-"
  hg[dl] <- ifelse(nchar(ref[dl]) == 1L,
                   paste0("p.", ref[dl], p1[dl], "del"),
                   paste0("p.", r1[dl], p1[dl], "_", rn[dl], p2[dl], "del"))
  canon[dl] <- nchar(ref[dl]) == (p2[dl] - p1[dl] + 1L)
  # in-frame insertion (VEP ref "-") — flanking residues unrecoverable here
  is_ <- ok & !fs & !s1 & ref == "-" & alt != "-"
  hg[is_] <- paste0("p.", p1[is_], "_", p2[is_], "ins", alt[is_]); canon[is_] <- FALSE
  # remainder: multi-residue substitution / complex delins. Trim the shared prefix
  # and suffix so VEP's "CP/C" at 1255-1256 reduces to the HGVS form p.P1256del
  # rather than an unreduced p.C1255_P1256delinsC (CDK12 in these data).
  ot <- which(ok & !fs & !s1 & !dl & !is_)
  if (length(ot)) {
    red <- vapply(ot, function(i) {
      r <- strsplit(ref[i], "")[[1]]; a <- strsplit(alt[i], "")[[1]]
      lo <- 0L
      while (lo < min(length(r), length(a)) && r[lo + 1L] == a[lo + 1L]) lo <- lo + 1L
      hi <- 0L
      while (hi < min(length(r), length(a)) - lo &&
             r[length(r) - hi] == a[length(a) - hi]) hi <- hi + 1L
      rr <- r[seq_len(max(length(r) - hi - lo, 0L)) + lo]
      av <- paste(a[seq_len(max(length(a) - hi - lo, 0L)) + lo], collapse = "")
      q1 <- p1[i] + lo; q2 <- q1 + max(length(rr), 1L) - 1L
      if (!nzchar(av) && length(rr) == 1L)  c(paste0("p.", rr, q1, "del"), "TRUE")
      else if (!nzchar(av) && length(rr) > 1L)
        c(paste0("p.", rr[1], q1, "_", rr[length(rr)], q2, "del"), "TRUE")
      else if (!length(rr))                                   # pure insertion
        c(paste0("p.", q1 - 1L, "_", q1, "ins", av), "FALSE")
      else if (length(rr) == 1L && nchar(av) == 1L) c(paste0("p.", rr, q1, av), "TRUE")
      else c(paste0("p.", rr[1], q1, "_", rr[length(rr)], q2, "delins", av), "FALSE")
    }, character(2))
    hg[ot] <- red[1, ]; canon[ot] <- red[2, ] == "TRUE"
  }
  list(hgvsp = hg, canonical = canon)
}
.hgv <- mk_hgvsp(raw$Amino_acids, raw$Protein_position, raw$Variant_Classification)
raw[, `:=`(HGVSp_Short = .hgv$hgvsp, hgvsp_canonical = .hgv$canonical,
           hgvsp_reconstructed = TRUE)]      # TRUE by construction; see note above
rm(.hgv)

# VAF: not carried in the MAF (t_depth empty). Pull the Mutect2 FORMAT/AF from
# each VCF and join on Tumor_Sample_Barcode + Chromosome + vcf_pos (+ alt for
# SNVs). Exact for SNVs; indel representations differ -> VAF stays NA there.
read_vcf_vaf <- function(vcf) {
  if (is.na(vcf)) return(NULL)
  if (grepl("\\.gz$", vcf)) {
    con <- gzfile(vcf, "rt"); lines <- readLines(con); close(con)
    v <- fread(text = paste(lines[!grepl("^##", lines)], collapse = "\n"), sep = "\t", header = TRUE, showProgress = FALSE)
  } else v <- fread(vcf, skip = "#CHROM", sep = "\t", header = TRUE, showProgress = FALSE)
  setnames(v, 1:10, c("CHROM","POS","ID","REF","ALT","QUAL","FILT","INFO","FORMAT","SMP"))
  fmt <- strsplit(v$FORMAT[1], ":", fixed = TRUE)[[1]]
  ai  <- match("AF", fmt)
  if (is.na(ai)) return(NULL)
  af  <- as.numeric(vapply(strsplit(v$SMP, ":", fixed = TRUE),
                           function(z) z[ai], character(1)))
  data.table(CHROM = v$CHROM, vcf_pos = v$POS, ALT = v$ALT, vaf = af)
}
vaf_dt <- map(maf_map$source_vcf, read_vcf_vaf) |>
  set_names(unique(maf_map$path))
vaf_dt <- imap(vaf_dt, function(d, p) {
  if (is.null(d)) return(NULL)
  d[, cell_line := maf_map$cell_line[maf_map$path == p]][]
}) %>% rbindlist()
raw <- merge(raw, vaf_dt,
             by.x = c("cell_line","Chromosome","vcf_pos","Tumor_Seq_Allele2"),
             by.y = c("cell_line","CHROM","vcf_pos","ALT"),
             all.x = TRUE, sort = FALSE)
message(sprintf("VAF joined for %.0f%% of records (SNVs; indels differ in rep.)",
                100 * mean(!is.na(raw$vaf))))

# 4. Re-filtering -------------------------------------------------------------
raw[, is_pass         := FILTER == "PASS"]
raw[, is_common_germ  := pop_af_max > POP_AF_MAX]
raw[, is_coding_nonsyn:= Variant_Classification %in% NONSYN]
# Retained somatic candidates: PASS AND not common germline.
raw[, retained := is_pass & !is_common_germ]
# Soft germline-like VAF flag (informational only — CANNOT remove germline in
# ~100%-pure lines, where het germline ~0.5 overlaps subclonal somatic and LOH
# drives real somatic drivers toward ~1.0; per Little 2021 UNMASC).
raw[, germline_like_vaf := !is.na(vaf) & vaf >= 0.90]
message(sprintf("VAF available for %.0f%% of retained coding candidates",
                100 * mean(!is.na(raw$vaf[raw$retained & raw$is_coding_nonsyn]))))

cat("\n=== Filtering cascade (all records -> retained candidates) ===\n")
print(data.frame(
  step = c("raw records","FILTER==PASS","  & not common germline (pop AF>0.001)",
           "  & coding non-synonymous"),
  n    = c(nrow(raw), sum(raw$is_pass),
           sum(raw$retained), sum(raw$retained & raw$is_coding_nonsyn))))
# Deposit the cascade from the complete input inventory (including recovered MAFs).
cascade <- raw[, .(raw = .N, pass = sum(is_pass), rare = sum(retained),
                   coding = sum(retained & is_coding_nonsyn)), by = .(cell_line, subtype)]
readr::write_csv(as_tibble(cascade), file.path(OUT, "wes_filter_cascade.csv"))
# Residual-burden reality check: even after filtering, per-line coding count far
# exceeds real HGSC exonic somatic burden (~50-80) => residual germline persists
# (tumor-only, pure lines; Halperin ~224/sample). This is WHY we report canonical
# drivers only and do NOT report mutation burden or non-driver frequencies.
per_line <- raw[retained & is_coding_nonsyn, .N, by = cell_line]
cat(sprintf("Retained coding-nonsyn per line: median %.0f (range %d-%d) -- far above real\n  HGSC exonic burden (~50-80) => residual germline remains; driver-only reporting.\n",
            median(per_line$N), min(per_line$N), max(per_line$N)))

# 5. Filtering VALIDATION — is the germline artefact gone? --------------------
watch_tbl <- map_dfr(ARTEFACT_WATCH, function(g) {
  pre  <- raw[Hugo_Symbol == g & is_coding_nonsyn,
              .(n_lines = uniqueN(cell_line))]$n_lines
  post <- raw[Hugo_Symbol == g & retained & is_coding_nonsyn,
              .(n_lines = uniqueN(cell_line))]$n_lines
  tibble(gene = g,
         lines_prefilter  = ifelse(length(pre)  == 0, 0L, pre),
         lines_retained   = ifelse(length(post) == 0, 0L, post),
         pct_prefilter    = 100 * ifelse(length(pre)  == 0, 0, pre)  / nrow(maf_map),
         pct_retained     = 100 * ifelse(length(post) == 0, 0, post) / nrow(maf_map))
})
cat("\n=== Artefact watch: coding-nonsyn calls, pre-filter vs retained (of 22 lines) ===\n")
print(as.data.frame(watch_tbl), row.names = FALSE, digits = 3)

# 6. POSITIVE CONTROL — TP53 in HGS -------------------------------------------
hgs_lines <- wes %>% filter(subtype == "HGS", cell_line %in% maf_map$cell_line) %>% pull(cell_line)
tp53_hgs  <- raw[Hugo_Symbol == "TP53" & retained & is_coding_nonsyn &
                   cell_line %in% hgs_lines, uniqueN(cell_line)]
cat(sprintf("\n=== POSITIVE CONTROL: TP53 in HGS ===\nTP53 (retained, coding) in %d/%d HGS WES lines = %.1f%% (expect ~universal; Ahmed 2010 HGS 96.7%%)\n",
            tp53_hgs, length(hgs_lines), 100 * tp53_hgs / length(hgs_lines)))
tp53_missing <- setdiff(hgs_lines,
                        raw[Hugo_Symbol=="TP53" & retained & is_coding_nonsyn, unique(cell_line)])
if (length(tp53_missing))
  cat("  HGS lines WITHOUT a retained TP53 hit:", paste(tp53_missing, collapse=", "), "\n")

# 7. Per-driver frequency BY SUBTYPE (n annotated) ----------------------------
drv <- raw[retained & is_coding_nonsyn & Hugo_Symbol %in% CANONICAL_DRIVERS]
subs <- wes %>% filter(cell_line %in% maf_map$cell_line)
freq_by_subtype <- expand_grid(Hugo_Symbol = CANONICAL_DRIVERS,
                               subtype = unique(subs$subtype)) %>%
  left_join(subs %>% count(subtype, name = "n_lines"), by = "subtype") %>%
  left_join(
    drv %>% as_tibble() %>% distinct(Hugo_Symbol, subtype, cell_line) %>%
      count(Hugo_Symbol, subtype, name = "n_mutated"),
    by = c("Hugo_Symbol","subtype")) %>%
  mutate(n_mutated = replace_na(n_mutated, 0L),
         pct = round(100 * n_mutated / n_lines, 1)) %>%
  arrange(match(Hugo_Symbol, CANONICAL_DRIVERS), subtype)
cat("\n=== Canonical-driver frequency by subtype (retained, coding; PER LINE) — head ===\n")
print(freq_by_subtype %>% filter(n_mutated > 0) %>% as.data.frame(), row.names = FALSE)

# =============================================================================
# 7b. PATIENT-LEVEL RE-TABULATION + SOMATIC-CONFIDENCE TIERING  [review revision]
# -----------------------------------------------------------------------------
# Peer review (§3.1, §3.3): per-LINE frequencies PSEUDOREPLICATE patients — the
# Mes-Masson OV/TOV sublines and -R re-derivations are the SAME tumour genome
# (confirmed at the variant level: all four 3133 lines share TP53 p.Q192* AND
# CDK12 p.L123fs; all three 2295 lines share TP53 p.I195T; the 1369 pair shares
# p.G244C). We therefore collapse to INDEPENDENT PATIENTS via
# metadata/line_family_map.csv (a driver counts for a patient if >=1 of that
# patient's lines carries it) and report n_patients ALONGSIDE n_lines.
#
# Separately, the `germline_like_vaf` flag does NOT add somatic confidence: in
# ~100%-pure lines a clonal somatic driver with LOH sits at VAF->1, so the flag
# labels bona fide somatic TP53 hotspots (R175H, Q192*) "germline-like". We
# replace it with a per-CALL candidate-prioritization TIER and EXCLUDE Tier3 (lowest-priority candidates) from headline frequencies.
# Every tier remains potentially germline; prioritization is not a somatic classifier.
# =============================================================================
ensure_family_map()          # [review revision] 15 writes this map; see 00_setup.R
fam <- readr::read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)
famm <- fam %>% filter(cell_line %in% maf_map$cell_line) %>%
  transmute(cell_line, patient_id, family, n_lines_in_family, patient_representative)
stopifnot("family map does not cover the MAF line set" =
            setequal(famm$cell_line, maf_map$cell_line))
message(sprintf("Patient collapse (MAF set): %d lines -> %d patients (%d multi-line families)",
                nrow(famm), n_distinct(famm$patient_id),
                sum(famm$n_lines_in_family > 1 & famm$patient_representative)))

# ---- (i) Candidate-prioritization tiering of every canonical-driver call ----------
TRUNC   <- c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
             "Splice_Site","Nonstop_Mutation","Translation_Start_Site")
INFRAME <- c("In_Frame_Del","In_Frame_Ins")
STRONG_TSG <- c("PTEN","RB1","NF1","ARID1A","SMARCA4")   # truncating LOF -> Tier1 (TP53 handled separately)
codon_of <- function(h) suppressWarnings(as.integer(stringr::str_extract(h, "[0-9]+")))
has_path <- function(cs) {
  # Match positive ClinVar terms explicitly; conflict/uncertain text does not
  # become pathogenic merely because it contains the substring "pathogenic".
  if (is.na(cs) || grepl("conflict|uncertain|benign", cs, ignore.case = TRUE)) return(FALSE)
  any(strsplit(cs, "[,/|&]", perl = TRUE)[[1]] %in% c("pathogenic", "likely_pathogenic"))
}

assign_tier <- function(gene, hgvsp, vclass, vaf, clin) {
  cod   <- codon_of(hgvsp)
  high_vaf <- !is.na(vaf) & vaf >= 0.85       # descriptive; does not establish LOH
  trunc <- vclass %in% TRUNC; infr <- vclass %in% INFRAME; mis <- vclass == "Missense_Mutation"
  path  <- has_path(clin)
  incod <- function(v) !is.na(cod) && cod %in% v
  # BRCA1/2: rare germline cannot be classified somatic in tumor-only WES.
  if (gene %in% c("BRCA1","BRCA2"))
    return(c("Tier3","BRCA1/2 in tumor-only WES: rare germline cannot be excluded (no matched normal); not defensibly somatic"))
  if (gene == "TP53") {
    if (trunc) return(c("Tier1","truncating TP53 (HGSC-typical LOF)"))
    if (mis && (path || incod(95:289)))
      return(c("Tier1", paste0("TP53 DBD missense", if (high_vaf) "+high VAF (LOH unassessed)" else "", if (path) "+ClinVar-path" else "")))
    return(c("Tier2","TP53 non-canonical candidate (outside DBD / no positive ClinVar annotation)"))
  }
  if (gene == "KRAS") {
    if (incod(c(12,13,61)))  return(c("Tier1", sprintf("KRAS codon %d hotspot", cod)))
    if (incod(c(59,117,146))) return(c("Tier2", sprintf("KRAS codon %d activating (non-canonical)", cod)))
    return(c("Tier3","KRAS non-hotspot"))
  }
  if (gene == "NRAS") {
    if (incod(c(12,13,61))) return(c("Tier1", sprintf("NRAS codon %d hotspot", cod)))
    return(c("Tier3","NRAS non-hotspot"))
  }
  if (gene == "BRAF") {
    if (grepl("600", hgvsp)) return(c("Tier1","BRAF V600 hotspot"))
    if (infr && incod(480:500)) return(c("Tier2","BRAF beta3-alphaC in-frame indel, plausibly activating (class II/III)"))
    return(c("Tier3","BRAF non-V600"))
  }
  if (gene == "PIK3CA") {
    if (incod(c(88,111,118,344,345,420,453,542,545,546,726,1043,1047,1049)) || path)
      return(c("Tier1","PIK3CA known activating hotspot"))
    return(c("Tier2","PIK3CA missense, non-canonical"))
  }
  if (gene == "CTNNB1") {
    if (incod(c(32,33,34,37,41,45))) return(c("Tier1","CTNNB1 degron (exon 3) hotspot"))
    return(c("Tier3","CTNNB1 non-degron missense (uncertain)"))
  }
  if (gene %in% STRONG_TSG) {
    if (trunc) return(c("Tier1", sprintf("truncating LOF in %s (canonical tumour-suppressor)", gene)))
    if (mis && path) return(c("Tier2", sprintf("%s ClinVar-pathogenic missense", gene)))
    return(c("Tier3", sprintf("%s missense/in-frame, uncharacterised", gene)))
  }
  if (gene == "CDK12") {
    if (trunc) return(c("Tier2","CDK12 truncating (HGSC driver; biallelic LOF not confirmable tumor-only)"))
    return(c("Tier3","CDK12 in-frame / low-VAF / non-hotspot missense"))
  }
  if (gene == "ARID1B") {
    if (trunc) return(c("Tier2","ARID1B truncating (SWI/SNF; less established in HGSC)"))
    return(c("Tier3","ARID1B non-truncating"))
  }
  if (gene == "CDKN2A") {
    if (trunc)        return(c("Tier2","CDKN2A truncating"))
    if (mis && path)  return(c("Tier2","CDKN2A ClinVar-pathogenic missense candidate"))
    return(c("Tier3","CDKN2A uncharacterised missense; high VAF alone does not establish pathogenicity or LOH"))
  }
  # everything else (e.g. SMARCA2 — not an established driver in these histotypes)
  return(c("Tier3", sprintf("%s: not an established driver here / uncharacterised", gene)))
}

drv_calls <- raw[retained & is_coding_nonsyn & Hugo_Symbol %in% CANONICAL_DRIVERS] %>%
  as_tibble() %>%
  left_join(famm %>% select(cell_line, patient_id, family), by = "cell_line") %>%
  rowwise() %>%
  mutate(.t = list(assign_tier(Hugo_Symbol, HGVSp_Short, Variant_Classification, vaf, CLIN_SIG)),
         tier = .t[[1]], rationale = .t[[2]]) %>%
  ungroup() %>% select(-.t)

# Hypermutator context (cross-note from wes-signatures): TOV21G is a candidate
# MSI-high / MMR-deficient clear-cell line (1,416 coding candidates = 6.9x panel
# median, 3.4x the next line; SBS6/15/44). In a hypermutator, passenger SNVs land
# in driver genes by chance, so its long driver list is inflated. Its CANONICAL CC
# drivers (ARID1A truncating, PIK3CA, KRAS) remain real; secondary/SWI-SNF calls
# should be read in that context (see f_wes_hypermutation). Flag it in the table.
# [review revision] ONE hypermutator rule for the whole pipeline. 07 previously
# used ">3x median" while 16 used "robust_z > 5 & fold > 3" — two implementations of
# one manuscript-facing flag. Both now call hypermutator_stats() from 00_setup.R.
tmb_line <- raw[retained & is_coding_nonsyn, .(n_coding = .N), by = cell_line]
hm_stats <- hypermutator_stats(tmb_line$n_coding, tmb_line$cell_line)
hypermut <- hm_stats$cell_line[hm_stats$is_hypermutator]
message(sprintf("Hypermutator (shared rule: robust_z > %g AND load > %gx median): %s",
                HYPERMUT_ROBUST_Z, HYPERMUT_FOLD, paste(hypermut, collapse = ", ")))

tiers_out <- drv_calls %>%
  mutate(context = if_else(cell_line %in% hypermut,
           "hypermutator/MSI-high candidate (see f_wes_hypermutation): non-canonical calls carry elevated passenger risk",
           NA_character_)) %>%
  transmute(cell_line, patient_id, family, subtype, gene = Hugo_Symbol,
            protein_change = HGVSp_Short, variant_classification = Variant_Classification,
            vaf, germline_like_vaf, clin_sig = CLIN_SIG, tier, rationale, context,
            Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2,
            hgvsp_canonical, hgvsp_reconstructed,
            interpretation = "heuristic candidate priority; somatic/germline origin and LOH are unassessed") %>%
  arrange(tier, gene, patient_id, cell_line)
readr::write_csv(tiers_out, file.path(OUT, "wes_driver_tiers.csv"))
cat("\n=== Candidate-prioritization tier distribution (of", nrow(tiers_out), "driver calls) ===\n")
print(tiers_out %>% count(tier), row.names = FALSE)
cat("\n=== Tier3 (cannot-exclude-germline) driver calls — EXCLUDED from headline freqs ===\n")
print(tiers_out %>% filter(tier == "Tier3") %>%
        select(cell_line, subtype, gene, protein_change, vaf, rationale) %>% as.data.frame(),
      row.names = FALSE)

# ---- (ii) Per-gene x subtype frequency: LINE vs PATIENT, all vs Tier1-2 -----
subtype_lines <- famm %>% left_join(distinct(subs, cell_line, subtype), by = "cell_line") %>%
  count(subtype, name = "n_lines")
subtype_pats  <- famm %>% left_join(distinct(subs, cell_line, subtype), by = "cell_line") %>%
  distinct(subtype, patient_id) %>% count(subtype, name = "n_patients")

count_muts <- function(calls) {
  list(
    line = calls %>% distinct(Hugo_Symbol, subtype, cell_line) %>%
            count(Hugo_Symbol, subtype, name = "n_lines_mut"),
    pat  = calls %>% distinct(Hugo_Symbol, subtype, patient_id) %>%
            count(Hugo_Symbol, subtype, name = "n_patients_mut"))
}
all_c  <- count_muts(drv_calls)
t12_c  <- count_muts(drv_calls %>% filter(tier %in% c("Tier1","Tier2")))

freq_patient <- expand_grid(Hugo_Symbol = CANONICAL_DRIVERS,
                            subtype = unique(subs$subtype)) %>%
  left_join(subtype_lines, by = "subtype") %>%
  left_join(subtype_pats,  by = "subtype") %>%
  left_join(all_c$line, by = c("Hugo_Symbol","subtype")) %>%
  left_join(all_c$pat,  by = c("Hugo_Symbol","subtype")) %>%
  left_join(t12_c$pat %>% rename(n_patients_mut_tier12 = n_patients_mut),
            by = c("Hugo_Symbol","subtype")) %>%
  mutate(across(c(n_lines_mut, n_patients_mut, n_patients_mut_tier12), ~replace_na(., 0L)),
         pct_lines             = round(100 * n_lines_mut            / n_lines,    1),
         pct_patients          = round(100 * n_patients_mut         / n_patients, 1),
         pct_patients_tier12   = round(100 * n_patients_mut_tier12  / n_patients, 1)) %>%
  arrange(match(Hugo_Symbol, CANONICAL_DRIVERS), subtype)
readr::write_csv(freq_patient, file.path(OUT, "wes_driver_freq_patient.csv"))

# Augment the existing by-subtype table with patient-level columns (non-destructive).
freq_by_subtype_aug <- freq_by_subtype %>%
  left_join(freq_patient %>% select(Hugo_Symbol, subtype, n_patients,
                                    n_patients_mut, pct_patients,
                                    n_patients_mut_tier12, pct_patients_tier12),
            by = c("Hugo_Symbol","subtype"))
readr::write_csv(freq_by_subtype_aug, file.path(OUT, "wes_driver_freq_by_subtype.csv"))

# ---- (iii) Headline reconciliations printed for the report ------------------
gp <- function(g, s) freq_patient %>% filter(Hugo_Symbol == g, subtype == s)
tp53_hgs_p  <- gp("TP53","HGS"); cdk12_hgs_p <- gp("CDK12","HGS")
tp53_all_lines <- drv_calls %>% filter(Hugo_Symbol=="TP53") %>% distinct(cell_line) %>% nrow()
tp53_all_pats  <- drv_calls %>% filter(Hugo_Symbol=="TP53") %>% distinct(patient_id) %>% nrow()
n_all_lines <- nrow(famm); n_all_pats <- n_distinct(famm$patient_id)
cat("\n=== HEADLINE PATIENT-LEVEL RECONCILIATIONS ===\n")
cat(sprintf("TP53: all-subtype %d/%d MAF lines (%.0f%%) = %d/%d patients (%.0f%%); HGS %d/%d lines (100%%) = %d/%d patients (%.0f%%)\n",
            tp53_all_lines, n_all_lines, 100*tp53_all_lines/n_all_lines,
            tp53_all_pats, n_all_pats, 100*tp53_all_pats/n_all_pats,
            tp53_hgs_p$n_lines_mut, tp53_hgs_p$n_lines,
            tp53_hgs_p$n_patients_mut, tp53_hgs_p$n_patients, tp53_hgs_p$pct_patients))
cat(sprintf("CDK12 (HGS): %d/%d lines (%.1f%%)  ->  %d/%d patients (%.1f%%)  ->  Tier1-2 only %d/%d patients (%.1f%%)\n",
            cdk12_hgs_p$n_lines_mut, cdk12_hgs_p$n_lines, cdk12_hgs_p$pct_lines,
            cdk12_hgs_p$n_patients_mut, cdk12_hgs_p$n_patients, cdk12_hgs_p$pct_patients,
            cdk12_hgs_p$n_patients_mut_tier12, cdk12_hgs_p$n_patients, cdk12_hgs_p$pct_patients_tier12))
brca2_defensible <- freq_patient %>% filter(Hugo_Symbol=="BRCA2") %>%
  summarise(tier12_patients = sum(n_patients_mut_tier12)) %>% pull()
cat(sprintf("BRCA2: per-line HGS %d/%d (5.9%%) + LGS 1/1; ALL-subtype 2/%d lines (9%%);  prioritized candidates (Tier1-2) = %d patients across ALL subtypes\n",
            gp("BRCA2","HGS")$n_lines_mut, gp("BRCA2","HGS")$n_lines, n_all_lines, brca2_defensible))
cat("\n=== Patient-level driver frequency (genes with >=1 patient) — head ===\n")
print(freq_patient %>% filter(n_patients_mut > 0) %>%
        select(Hugo_Symbol, subtype, n_lines, n_lines_mut, pct_lines,
               n_patients, n_patients_mut, pct_patients, n_patients_mut_tier12) %>%
        as.data.frame(), row.names = FALSE)

# 8. Filtered variant table (retained coding/splice candidates) ---------------
# *** THIS IS THE PIPELINE'S CENTRAL MUTATION TABLE *** — 10, 11, 16 and 18 all read
# output/wes_mutations_filtered.csv. The write below was previously COMMENTED OUT and
# the script instead READ the file in order to assert its own row count, so from a
# clean output/ nothing produced it, 07 aborted, and four downstream scripts failed
# (every mutation-derived number in the manuscript was non-reproducible). The write
# is restored and the drift guard is now a CONTENT HASH: comparing nrow() alone let
# any logic change that preserved 6,036 rows pass silently.
keep_cols <- c("cell_line","subtype","Hugo_Symbol","Chromosome","Start_Position",
               "End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification",
               "Variant_Type","HGVSp_Short","hgvsp_canonical","hgvsp_reconstructed",
               "Transcript_ID","Consequence","IMPACT","dbSNP_RS",
               "Existing_variation","CLIN_SIG","pop_af_max","vaf","germline_like_vaf",
               "is_driver","FILTER")
raw[, is_driver := Hugo_Symbol %in% CANONICAL_DRIVERS]
filtered_out <- raw[retained & is_coding_nonsyn, ..keep_cols][order(cell_line, Chromosome, Start_Position)]
readr::write_csv(as_tibble(filtered_out), file.path(OUT, "wes_mutations_filtered.csv"))

# ---- Content-hash drift guard ------------------------------------------------
# sha256 over the canonically-sorted identity columns of the retained variant set:
# any change to WHICH variants are retained, their coordinates, their classification
# or their reconstructed protein change trips this assertion. If a filtering change
# is intentional, review the new hash and re-pin FILTERED_EXPECT$sha256 deliberately.
filtered_sha256 <- function(d) {
  k <- c("cell_line","Hugo_Symbol","Chromosome","Start_Position","End_Position",
         "Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","HGVSp_Short")
  x <- as.data.frame(d)[, k]
  x[] <- lapply(x, function(v) ifelse(is.na(v), "", as.character(v)))
  x <- x[order(x$cell_line, x$Chromosome, x$Start_Position, x$End_Position,
               x$Reference_Allele, x$Tumor_Seq_Allele2, x$Hugo_Symbol), ]
  digest::digest(paste(do.call(paste, c(unname(x), list(sep = "\t"))), collapse = "\n"),
                 algo = "sha256")
}
FILTERED_EXPECT <- list(
  sha256   = "4febdd74482a09ef5818be7f6c7ffe2cd53f09e5e750bfcb4d86bacb1967e6b4",
  n_rows   = 6194L,     # 582,474 raw -> 16,081 PASS -> 15,995 rare -> 6,194 coding nonsyn
  n_lines  = 23L,
  med_line = 206)     # per-line median (range 133-1,416; TOV21G is the 1,416)
got <- list(sha256   = filtered_sha256(filtered_out),
            n_rows   = nrow(filtered_out),
            n_lines  = uniqueN(filtered_out$cell_line),
            med_line = median(filtered_out[, .N, by = cell_line]$N))
cat(sprintf("\n=== wes_mutations_filtered.csv WRITTEN: %d rows / %d lines / median %.1f per line ===\n    sha256(identity cols) = %s\n",
            got$n_rows, got$n_lines, got$med_line, got$sha256))
stopifnot(
  "retained-variant COUNT drifted from the pinned value" =
    got$n_rows == FILTERED_EXPECT$n_rows,
  "retained-variant LINE COUNT drifted from the pinned value" =
    got$n_lines == FILTERED_EXPECT$n_lines,
  "per-line median retained-variant count drifted from the pinned value" =
    isTRUE(all.equal(got$med_line, FILTERED_EXPECT$med_line)),
  "retained-variant CONTENT HASH drifted: the filtering logic or the variant annotation changed. Review the diff, then re-pin FILTERED_EXPECT$sha256 deliberately." =
    identical(got$sha256, FILTERED_EXPECT$sha256))
cat(sprintf("    HGVSp_Short: all %d reconstructed from Amino_acids/Protein_position; %d canonical, %d NOT canonical (never quote as identifiers), %d not derivable\n",
            got$n_rows, sum(filtered_out$hgvsp_canonical, na.rm = TRUE),
            sum(!filtered_out$hgvsp_canonical, na.rm = TRUE),
            sum(is.na(filtered_out$hgvsp_canonical))))

# 9. Canonical-driver ONCOPRINT — patient-family track + somatic-tier encoding
# -----------------------------------------------------------------------------
# Rebuilt with ComplexHeatmap::oncoPrint (was maftools) so we can (a) show a
# PATIENT-FAMILY annotation track — the four identical 3133 columns read as ONE
# patient — and (b) encode the candidate-prioritization TIER per cell (fill = tier;
# box height = truncating-LOF vs missense/in-frame). Gene labels on the right
# report n_PATIENTS (not lines), and per-line % is suppressed (pseudoreplicated).
n_pat_total <- n_distinct(famm$patient_id)                       # 16 patients (MAF set)

# Per (gene,line) cell: top (most-confident) tier + whether it is a truncating LOF.
cell_tab <- drv_calls %>%
  mutate(tier_n = as.integer(sub("Tier", "", tier)),
         is_lof = Variant_Classification %in% TRUNC) %>%
  group_by(Hugo_Symbol, cell_line) %>%
  summarise(top = min(tier_n),
            lof = any(is_lof[tier_n == min(tier_n)]), .groups = "drop") %>%
  mutate(cat = paste0("T", top, ifelse(lof, "_LOF", "_MIS")))

gene_order <- freq_patient %>% group_by(Hugo_Symbol) %>%
  summarise(np = sum(n_patients_mut), nl = sum(n_lines_mut), .groups = "drop") %>%
  filter(np > 0) %>% arrange(desc(np), desc(nl)) %>% pull(Hugo_Symbol)

col_meta <- subs %>% left_join(famm, by = "cell_line") %>%
  mutate(subtype = factor(subtype, levels = c("HGS","CC","EC","LGS","MC")),
         fam_grp = ifelse(n_lines_in_family > 1, family, "singleton")) %>%
  arrange(subtype, patient_id, cell_line)
col_order <- col_meta$cell_line

mat <- matrix("", length(gene_order), nrow(col_meta),
              dimnames = list(gene_order, col_order))
for (i in seq_len(nrow(cell_tab))) {
  g <- cell_tab$Hugo_Symbol[i]; s <- cell_tab$cell_line[i]
  if (g %in% gene_order) mat[g, s] <- cell_tab$cat[i]
}

tier_fill <- c(`1` = COOK_RUST, `2` = "#E8A13A", `3` = "#BEBEBE")   # rust / amber / grey
draw_cell <- function(x, y, w, h, t, lof)
  grid.rect(x, y, w * 0.88, h * ifelse(lof, 0.94, 0.40),
            gp = gpar(fill = tier_fill[t], col = "white", lwd = 0.4))
alter_fun <- list(
  background = function(x, y, w, h) grid.rect(x, y, w*0.95, h*0.95, gp = gpar(fill = "#F2F2F2", col = NA)),
  T1_LOF = function(x,y,w,h) draw_cell(x,y,w,h,"1",TRUE),
  T1_MIS = function(x,y,w,h) draw_cell(x,y,w,h,"1",FALSE),
  T2_LOF = function(x,y,w,h) draw_cell(x,y,w,h,"2",TRUE),
  T2_MIS = function(x,y,w,h) draw_cell(x,y,w,h,"2",FALSE),
  T3_LOF = function(x,y,w,h) draw_cell(x,y,w,h,"3",TRUE),
  T3_MIS = function(x,y,w,h) draw_cell(x,y,w,h,"3",FALSE))
alter_col <- c(T1_LOF=COOK_RUST, T1_MIS=COOK_RUST, T2_LOF="#E8A13A",
               T2_MIS="#E8A13A", T3_LOF="#BEBEBE", T3_MIS="#BEBEBE")

# Annotation tracks: patient family (multi-line families coloured; singletons grey)
fam_lv   <- c("1369","2295","3121","3133","3291")
fam_cols <- c(setNames(c("#1B9E77","#7570B3","#E7298A","#66A61E","#A6761D"), fam_lv),
              singleton = "grey88")
sub_cols <- c(HGS="#0072B2", CC="#E69F00", EC="#009E73", LGS="#56B4E9", MC="#CC79A7")
# TMB track (retained coding candidates/line): flags the TOV21G hypermutator spike
# so the figure is self-consistent (a reviewer sees WHY TOV21G carries many calls).
tmb_vec  <- setNames(tmb_line$n_coding, tmb_line$cell_line)[col_order]
tmb_fill <- ifelse(col_order %in% hypermut, COOK_RUST, COOK_NAVY)
top_anno <- HeatmapAnnotation(
  `coding\ncandidates` = anno_barplot(tmb_vec, border = FALSE,
       gp = gpar(fill = tmb_fill, col = NA), height = unit(1.3, "cm"),
       axis_param = list(gp = gpar(fontsize = 6))),
  `Patient family` = col_meta$fam_grp,
  Subtype          = as.character(col_meta$subtype),
  col = list(`Patient family` = fam_cols, Subtype = sub_cols),
  annotation_name_gp = gpar(fontsize = 7.5), simple_anno_size = unit(4, "mm"),
  gap = unit(1, "mm"),
  annotation_legend_param = list(
    `Patient family` = list(at = c(fam_lv, "singleton"),
                            labels = c(paste0("family ", fam_lv), "single-line patient"))))

np_tab <- freq_patient %>% group_by(Hugo_Symbol) %>%
  summarise(np = sum(n_patients_mut), .groups = "drop")
np_vec <- setNames(np_tab$np, np_tab$Hugo_Symbol)[gene_order]
right_anno <- rowAnnotation(
  `patients (of 16)` = anno_barplot(np_vec, border = FALSE, gp = gpar(fill = COOK_NAVY),
                                    width = unit(2.4, "cm"), axis_param = list(gp = gpar(fontsize = 7))),
  n = anno_text(sprintf("%d", np_vec), gp = gpar(fontsize = 7.5), just = "left"),
  annotation_name_gp = gpar(fontsize = 8), annotation_name_rot = 0)

lgd_tier <- Legend(title = "Candidate priority (fill)",
  labels = c("Tier1 canonical", "Tier2 plausible", "Tier3 uncertain candidate"),
  legend_gp = gpar(fill = c(COOK_RUST, "#E8A13A", "#BEBEBE")))
lgd_class <- Legend(title = "Variant class (box height)",
  labels = c("Truncating / splice (LOF)", "Missense / in-frame"),
  graphics = list(
    function(x,y,w,h) grid.rect(x,y,w,h*0.94, gp = gpar(fill="grey40", col=NA)),
    function(x,y,w,h) grid.rect(x,y,w,h*0.40, gp = gpar(fill="grey40", col=NA))))

op <- oncoPrint(mat, alter_fun = alter_fun, col = alter_col,
  row_order = gene_order, column_order = col_order,
  top_annotation = top_anno, right_annotation = right_anno,
  show_pct = FALSE, show_heatmap_legend = FALSE,
  remove_empty_columns = FALSE, remove_empty_rows = FALSE,
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  column_names_gp = gpar(fontsize = 7), pct_gp = gpar(fontsize = 8),
  column_title = "OvCAN WES canonical drivers - patient-collapsed, candidate-tiered (re-filtered, tumor-only)",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"))

caveat <- paste0(
  "TUMOR-ONLY WES (no matched normal): coding candidates after FILTER==PASS + gnomAD/1000G/ESP removal. ",
  "Columns are cell LINES grouped by patient family (top track); gene bars count PATIENTS (n=", n_pat_total,
  "), not lines. Tier3 excluded from headline frequencies. n(HGS)=", length(hgs_lines),
  " lines / ", tp53_hgs_p$n_patients, " patients. ",
  hypermut[1], " is a hypermutator (MSI-high candidate; TMB bar) - its secondary calls carry passenger risk (see f_wes_hypermutation).")

render_onco <- function(open_dev) {
  open_dev()
  draw(op, annotation_legend_list = list(lgd_tier, lgd_class),
       annotation_legend_side = "right", padding = unit(c(12, 2, 2, 2), "mm"))
  grid.text(caveat, x = 0.5, y = 0.015,
            gp = gpar(fontsize = 6.5, fontface = "italic", col = "grey25"))
  dev.off()
}
render_onco(function() pdf(file.path(FIGS, "f_wes_oncoplot.pdf"), width = 12, height = 7))
render_onco(function() ragg::agg_png(file.path(PROJ, "reports", "assets", "f_wes_oncoplot.png"),
                                     width = 12, height = 7, units = "in", res = 200))
render_onco(function() pdf(file.path(FIGS, "07_oncoplot_drivers.pdf"), width = 12, height = 7))  # legacy name
message("Wrote figs/f_wes_oncoplot.pdf + reports/assets/f_wes_oncoplot.png (and legacy figs/07_oncoplot_drivers.pdf)")

message("\n07_wes_mutations.R complete. Outputs: output/wes_mutations_filtered.csv (NOW WRITTEN ",
        "by this script — content-hash guarded), output/wes_driver_freq_by_subtype.csv (augmented), ",
        "output/wes_driver_freq_patient.csv, output/wes_driver_tiers.csv ; ",
        "figure figs/f_wes_oncoplot.pdf (+ reports/assets PNG).")
