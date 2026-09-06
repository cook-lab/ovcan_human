# =============================================================================
# Script: 08_wes_cnv.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Genome-wide copy-number landscape for the 23 generated WES lines from
#          target-only CNVkit reconstruction. Per-sample median-centre the log2
#          segment means, build a binned genome heatmap highlighting canonical
#          HGSC events, and compute CN QC (n segments, fraction genome altered).
#          Phase 3 (WES), step 2 of 3.
#          *** REVISED (peer review, §3.1/§3.4) ***
#          - AUTOSOME-RESTRICTED FGA (chrX is a pooled-normal sex-composition
#            artifact — see below — and was previously INCLUDED in FGA).
#          - chrX dropped from the heatmap (chrY already was).
#          - PATIENT-LEVEL arm-event frequencies (collapse sublines to patients
#            via metadata/line_family_map.csv); patient-family annotation track.
#          - Explicit caveats: median-centering bias in high-FGA genomes;
#            pooled UNMATCHED public normals are a VALIDITY concern, not hygiene.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# =============================================================================
#
# ASSAY NOTES (verified on disk; see reports/assessment/02_data_inventory.md §4):
#  - CNVkit v0.9.10, GRCh38 (GATK assembly38), WES target/antitarget intervals,
#    `--drop-low-coverage`. One pooled normal reference (byte-identical across all
#    23 tumors), built from 5 PUBLIC healthy exomes (SRR4039087-97, PRJNA339046,
#    a skin-cancer study). The author confirms the same capture kit was used
#    for model and public-reference exomes (v5 comment31). Provider records
#    identify SeqCap EZ Exome v3; the original and lifted target BEDs are recovered.
#    These are unmatched normals; relative ratios and median centring do not
#    establish absolute ploidy or a diploid tumour baseline.
#  - CHR X ARTIFACT: the public normals' sex composition and CNVkit sex-reference configuration
#    are unverified. A systematic chrX shift is observed but its origin cannot
#    be assigned to a sex-reference effect rather than tumour copy number. chrY was already dropped; we now also drop chrX
#    from the heatmap AND restrict FGA / gain-loss / arm calls to AUTOSOMES.
#  - MEDIAN-CENTERING CAVEAT: per-sample probe-weighted median-centering assumes
#    the modal copy state is neutral. In genomes where >50% is altered (several
#    HGSC lines here have autosome FGA > 0.7) the median may sit on an altered state, so the
#    neutral baseline and gain/loss fractions can be biased. Flagged
#    per-line below (fga_auto_0.2 > 0.7).
#  - The September coverage audit identified four positive-depth antitarget bins
#    that overlap capture targets and carry extreme ratios. Script 29 removes
#    all antitarget rows and reruns native CNVkit 0.9.10 CBS on the unchanged
#    archived target ratios. Use its .cns log2 profiles, not integer .call.cns.
#    Original inputs remain unchanged. chrY and alt/random/Un contigs are dropped.
#  - No allele-specific / BAF output exists (no `--vcf` was passed) — relevant to
#    the HRD feasibility question handled in 09_wes_hrd.R.
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(matrixStats)
  library(GenomicRanges); library(ComplexHeatmap); library(circlize)
})
set.seed(SEED)

CNV_MANIFEST <- file.path(OUT, "wes_cnv_target_only", "manifest.csv")
stopifnot("Run scripts/29_wes_cnv_target_only.py before script 08" = file.exists(CNV_MANIFEST))

MAIN_CHR   <- paste0("chr", c(1:22, "X"))     # loaded (for the chrX diagnostic); chrY/alt dropped
AUTOSOMES  <- paste0("chr", 1:22)             # headline analysis is autosome-restricted (chrX = artifact)
FGA_THRESH <- 0.20                            # |log2| gain/loss call for FGA
BIN_WIDTH  <- 1e7                             # 10 Mb heatmap bins
# ARM-LEVEL CALL PARAMETERS [review revision] — the manuscript quotes six arm
# frequencies but never stated the threshold that produced them. Both parameters are
# now named here, written into wes_cnv_arm_freq_patient.csv, and swept in a
# sensitivity table (wes_cnv_arm_freq_sensitivity.csv) so the referee's question
# ("what log2 threshold defines a gain or loss?") is answerable from the outputs.
ARM_LOG2_THRESH <- FGA_THRESH   # |log2c_auto| for a segment to count as gained/lost
ARM_MAJORITY    <- 0.50         # fraction of arm length that must be gained/lost
ARM_THRESH_SWEEP   <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.40)
ARM_MAJORITY_SWEEP <- c(0.30, 0.40, 0.50, 0.60, 0.70)

# GRCh38 (hg38) main-chromosome lengths (bp) for a fixed bin layout.
HG38 <- c(chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
          chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
          chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
          chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
          chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
          chr21=46709983, chr22=50818468, chrX=156040895)

# UCSC hg38 cytobands, downloaded 2026-09-05; 0-based, half-open intervals.
# Exclude the full acen interval from both arms, rather than approximating a
# centromere by a midpoint. The pinned local reference makes this reproducible.
CYTO_FILE <- file.path(PROJ, "data", "reference", "hg38_cytoBand.tsv")
stopifnot("hg38 cytoband reference missing" = file.exists(CYTO_FILE),
          "hg38 cytoband reference content changed" =
            digest::digest(file = CYTO_FILE, algo = "sha256") ==
            "ce1b6033a5243e7c5022660b952d2ec33243e307e909afcaeec1894641a5208f")
cyto <- fread(CYTO_FILE, col.names = c("chromosome", "start", "end", "band", "stain"))
centromeres <- cyto[chromosome %in% AUTOSOMES & stain == "acen",
                    .(cen_start = min(start), cen_end = max(end)), by = chromosome]
stopifnot(nrow(centromeres) == 22L)
arm_bounds <- rbind(
  centromeres[, .(chromosome, arm = paste0(sub("chr", "", chromosome), "p"),
                  arm_start = 0, arm_end = cen_start)],
  centromeres[, .(chromosome, arm = paste0(sub("chr", "", chromosome), "q"),
                  arm_start = cen_end, arm_end = unname(HG38[chromosome]))])
arm_bounds[, arm_bp := arm_end - arm_start]
arm_bounds[, assessable := !arm %in% paste0(c(13, 14, 15, 21, 22), "p")]
readr::write_csv(as_tibble(arm_bounds), file.path(OUT, "wes_cnv_arm_boundaries.csv"))

# CNVkit dirs are "<line>_P<passage>_new" (passage mid-string): strip _new, then
# the passage suffix, then separators. Harmless on bare cell_line inputs.
canon_key <- function(x)
  toupper(gsub("[-_]", "", sub("_[Pp][0-9]+$", "", sub("_new$", "", x))))

# Canonical HGSC copy-number events to highlight (gene ~ hg38 position, Mb).
HGSC_EVENTS <- tribble(
  ~label,           ~dir,   ~chrom, ~pos_mb,
  "MECOM/PRKCI 3q26","gain","chr3", 169.5,
  "SOX2 3q26",       "gain","chr3", 181.7,
  "MYC 8q24",        "gain","chr8", 127.7,
  "CCNE1 19q12",     "gain","chr19",29.8,
  "20q",             "gain","chr20",50.0,
  "TP53 17p",        "loss","chr17",7.7,
  "RB1 13q",         "loss","chr13",48.3,
  "PTEN 10q",        "loss","chr10",87.9)

# 1. Sample sheet -> 23 WES-CNV lines + subtype -------------------------------
ss  <- readr::read_csv(SAMPLE_SHEET, show_col_types = FALSE)
wes <- ss %>% filter(wes_cnv == "Y", provenance == "generated") %>%
  transmute(cell_line, subtype, key = canon_key(cell_line))
stopifnot("Expected 23 WES-CNV lines" = nrow(wes) == 23)

# 1b. Patient / subline-family map (peer review §3.1) -------------------------
ensure_family_map()          # [review revision] 15 writes this map; see 00_setup.R
fam <- readr::read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)
famm <- as.data.table(fam %>% transmute(cell_line, patient_id, family, n_lines_in_family,
                                        is_multiline_family, patient_representative))

# 2. Read the hash-pinned target-only segment manifest -----------------------
cns_tbl <- readr::read_csv(CNV_MANIFEST, show_col_types = FALSE) %>%
  transmute(cell_line, path = file.path(PROJ, cns_path), cns_sha256) %>%
  left_join(wes, by = "cell_line")
stopifnot("Expected all 23 unique CNV lines" = nrow(cns_tbl) == 23L &&
            !anyDuplicated(cns_tbl$cell_line) && setequal(cns_tbl$cell_line, wes$cell_line),
          "Missing target-only CNS file" = all(file.exists(cns_tbl$path)))
actual_hash <- vapply(cns_tbl$path, function(p) digest::digest(file = p, algo = "sha256"), character(1))
stopifnot("Target-only CNS hash differs from script 29 manifest" = all(actual_hash == cns_tbl$cns_sha256))
message("Loaded 23 hash-verified target-only CNVkit segment files from script 29")
# Attach patient/family; order rows subtype -> patient -> line so families are adjacent.
cns_tbl <- cns_tbl %>% left_join(famm, by = "cell_line") %>%
  mutate(subtype = factor(subtype, levels = c("HGS","CC","EC","LGS","MC"))) %>%
  arrange(subtype, patient_id, cell_line)
message(sprintf("CNV patient collapse: %d lines -> %d patients (HGS: %d lines -> %d patients)",
                nrow(cns_tbl), n_distinct(cns_tbl$patient_id),
                sum(cns_tbl$subtype == "HGS"),
                n_distinct(cns_tbl$patient_id[cns_tbl$subtype == "HGS"])))

# 3. Read segments, keep main chroms, per-sample median-centre ----------------
# Two centrings: log2c = centred on ALL main chroms (legacy, incl chrX); log2c_auto
# = centred on AUTOSOMES only (headline). chrX is retained in the table for the
# diagnostic but excluded from every headline metric below.
read_cns <- function(path, cell_line, subtype) {
  dt <- fread(path, sep = "\t",
              select = c("chromosome","start","end","log2","depth","probes","weight"),
              showProgress = FALSE)
  dt <- dt[chromosome %in% MAIN_CHR]
  wmed  <- matrixStats::weightedMedian(dt$log2, w = dt$probes, na.rm = TRUE)
  auto  <- dt$chromosome %in% AUTOSOMES
  wmeda <- matrixStats::weightedMedian(dt$log2[auto], w = dt$probes[auto], na.rm = TRUE)
  dt[, `:=`(cell_line = cell_line, subtype = as.character(subtype),
            log2_raw = log2, log2c = log2 - wmed, log2c_auto = log2 - wmeda)]
  dt[]
}
segs <- pmap(list(cns_tbl$path, cns_tbl$cell_line, cns_tbl$subtype), read_cns) %>%
  rbindlist()
segs[, chromosome := factor(chromosome, levels = MAIN_CHR)]
message(sprintf("Loaded %d main-chrom segments across %d lines (median %.0f/line)",
                nrow(segs), uniqueN(segs$cell_line),
                median(segs[, .N, by = cell_line]$N)))

# Sanity: autosome vs all-chrom centring differ negligibly (median is autosome-
# dominated), so switching FGA to autosomes is a chrX-EXCLUSION effect, not a
# re-centring artefact.
center_delta <- segs[, .(d = abs(unique(log2c - log2c_auto))), by = cell_line]
message(sprintf("Median-centre shift (all-chrom vs autosome-only): max %.4f, median %.4f log2 across lines",
                max(center_delta$d), median(center_delta$d)))

readr::write_csv(
  segs[, .(cell_line, subtype, chromosome, start, end, probes, weight,
           log2_raw, log2c, log2c_auto)] %>% as_tibble(),
  file.path(OUT, "wes_cnv_segments.csv"))

# 3b. chrX artifact diagnostic ------------------------------------------------
segs[, seg_len := end - start]
xdiag <- segs[chromosome == "chrX",
              .(chrX_median_log2c   = matrixStats::weightedMedian(log2c_auto, probes, na.rm = TRUE),
                chrX_frac_altered   = sum(seg_len[abs(log2c_auto) > 0.20]) / sum(seg_len)),
              by = cell_line]
cat("\n=== chrX diagnostic (log2c relative to AUTOSOME baseline) — systematic shift; origin unconfirmed ===\n")
print(as.data.frame(xdiag[order(chrX_median_log2c)]), row.names = FALSE, digits = 3)
cat(sprintf("chrX median-log2c across lines: median %.2f (IQR %.2f to %.2f); median chrX frac-altered %.2f\n  => a systematic shift of uncertain origin; sex-reference configuration and tumour CN cannot be separated here. EXCLUDED from FGA/heatmap.\n",
            median(xdiag$chrX_median_log2c),
            quantile(xdiag$chrX_median_log2c, .25), quantile(xdiag$chrX_median_log2c, .75),
            median(xdiag$chrX_frac_altered)))

# 4. CN QC — fraction genome altered: AUTOSOME-RESTRICTED (headline) + legacy --
fga_calc <- function(d, l2c) {
  d[, .(n_segments        = .N,
        total_mb          = sum(seg_len) / 1e6,
        fga_0.2           = sum(seg_len[abs(get(l2c)) > 0.20]) / sum(seg_len),
        fga_0.3           = sum(seg_len[abs(get(l2c)) > 0.30]) / sum(seg_len),
        frac_gain         = sum(seg_len[get(l2c) >  0.20]) / sum(seg_len),
        frac_loss         = sum(seg_len[get(l2c) < -0.20]) / sum(seg_len)),
    by = .(cell_line, subtype)]
}
qc_legacy <- fga_calc(segs, "log2c")                       # incl chrX (reproduces prior file)
qc_auto   <- fga_calc(segs[chromosome %in% AUTOSOMES], "log2c_auto")   # AUTOSOME-restricted (headline)

cnv_qc <- qc_legacy[, .(cell_line, subtype,
                        n_segments, total_assessed_mb = total_mb,
                        fga_0.2, fga_0.3, frac_gain, frac_loss)] %>%
  merge(qc_auto[, .(cell_line,
                    n_segments_auto = n_segments,
                    total_assessed_autosome_mb = total_mb,
                    fga_auto_0.2 = fga_0.2, fga_auto_0.3 = fga_0.3,
                    frac_gain_auto = frac_gain, frac_loss_auto = frac_loss)],
        by = "cell_line") %>%
  merge(xdiag, by = "cell_line") %>%
  merge(famm[, .(cell_line, patient_id, n_lines_in_family, patient_representative)],
        by = "cell_line") %>%
  mutate(high_fga_flag = fga_auto_0.2 > 0.70) %>%           # median-centering caveat applies
  as.data.table()
setorder(cnv_qc, subtype, -fga_auto_0.2)
readr::write_csv(as_tibble(cnv_qc), file.path(OUT, "wes_cnv_fga.csv"))

cat("\n=== CN QC — AUTOSOME-restricted FGA (headline) vs legacy with-chrX ===\n")
print(as.data.frame(cnv_qc[, .(cell_line, subtype, patient_id,
        fga_auto_0.2 = round(fga_auto_0.2, 3), fga_withX_0.2 = round(fga_0.2, 3),
        chrX_med = round(chrX_median_log2c, 2), high_fga_flag)]), row.names = FALSE)
# Subtype medians at PATIENT level.
# NB (review revision): a patient's value is the MEAN of that patient's lines, so a
# "patient-level median" over a 2-patient subtype is a median of two means. For
# clear cell specifically, this is the midpoint of two different patient values,
# not an estimate of a typical clear-cell FGA. n_patients is printed alongside for this
# reason; rare-subtype rows are 1-2 patients and must be read as such.
pat_fga <- cnv_qc[, .(fga_auto_0.2 = mean(fga_auto_0.2), n_lines_for_patient = .N),
                  by = .(patient_id, subtype)]
# Label self-consistency (review revision): the subtype column is the LABELLED
# histotype from samples.csv. For n<=2 subtypes we print the constituent lines so a
# reader cannot mistake e.g. "EC (n=1)" for a histotype estimate — that row is
# TOV112D, the line scripts 10/31 reclassify AWAY from endometrioid (TP53-mut +
# SMARCA4 truncation -> dedifferentiated carcinoma, Karnezis 2021). Same for LGS
# (TOV81D) and MC (TOV2414).
sub_lines <- cnv_qc[, .(lines = paste(cell_line, collapse = ", ")), by = subtype]
fga_by_sub <- cnv_qc[, .(n_lines = .N, fga_line_med = round(median(fga_auto_0.2), 3)), by = subtype][
      pat_fga[, .(n_patients = .N, fga_patient_med = round(median(fga_auto_0.2), 3)), by = subtype],
      on = "subtype"][sub_lines, on = "subtype"]
fga_by_sub[, label_note := fifelse(
  n_lines <= 2,
  paste0("labelled subtype, n<=2 lines (", lines, ") — see 10 for reclassification"),
  "")]
cat("\n=== Autosome FGA by LABELLED subtype: line-level vs patient-level median ===\n")
print(as.data.frame(fga_by_sub[, .(subtype, n_lines, fga_line_med, n_patients,
                                   fga_patient_med, label_note)]), row.names = FALSE)

# chrX-exclusion effect — SCOPE IT EXPLICITLY (review revision, audit #55).
# The manuscript quoted "chrX exclusion barely changes high-FGA HGSC genomes
# (0.61->0.62)", but 0.61->0.62 is the median across ALL 23 lines, not the HGSC
# subset. Both are computed and labelled here so they cannot be confused again.
delta_all <- cnv_qc[, .(n = .N,
                        withX   = median(fga_0.2),
                        autosome = median(fga_auto_0.2),
                        delta   = median(fga_0.2 - fga_auto_0.2))]
delta_hgs <- cnv_qc[subtype == "HGS", .(n = .N,
                        withX   = median(fga_0.2),
                        autosome = median(fga_auto_0.2),
                        delta   = median(fga_0.2 - fga_auto_0.2))]
cat("\n=== chrX-exclusion effect on FGA — median, BY SCOPE (do not conflate) ===\n")
print(as.data.frame(rbind(
  data.frame(scope = "ALL WES-CNV lines", delta_all),
  data.frame(scope = "HGS lines only",    delta_hgs))), row.names = FALSE, digits = 3)
cat(sprintf("  all %d lines: with-chrX %.3f -> autosome %.3f | HGS only (n=%d): %.3f -> %.3f\n",
            delta_all$n, delta_all$withX, delta_all$autosome,
            delta_hgs$n, delta_hgs$withX, delta_hgs$autosome))
cat(sprintf("High-FGA lines (>0.70, median-centering caveat): %s\n",
            paste(cnv_qc$cell_line[cnv_qc$high_fga_flag], collapse = ", ")))

# 5. ARM-LEVEL events: per line, then collapsed to PATIENTS (peer review §3.1) -
# Intersect each segment with both chromosome arms. A segment spanning the
# centromere contributes only its overlap to each arm; centromeric sequence is
# excluded. Acrocentric p arms are not assessed. Fractions use the full annotated
# non-centromeric arm length, so >50% means >50% of the arm (not >50% of a subset
# of segments that happened to have their midpoint on that arm).
seg_a <- merge(segs[chromosome %in% AUTOSOMES], arm_bounds,
               by = "chromosome", allow.cartesian = TRUE)
seg_a[, overlap_bp := pmax(0, pmin(end, arm_end) - pmax(start, arm_start))]
seg_a <- seg_a[overlap_bp > 0 & assessable]
stopifnot("segments overlap beyond chromosome arm length" =
  seg_a[, .(covered = sum(overlap_bp), arm_bp = arm_bp[1]),
        by = .(cell_line, arm)][, all(covered <= arm_bp)])

# Parameterised arm caller so the headline table and the threshold-sensitivity
# sweep share ONE implementation (review revision: the threshold was previously
# a bare literal inside this expression and appeared nowhere in any output).
arm_freq_at <- function(thr = ARM_LOG2_THRESH, maj = ARM_MAJORITY) {
  al <- seg_a[, .(gain_frac = sum(overlap_bp[log2c_auto >  thr]) / arm_bp[1],
                  loss_frac = sum(overlap_bp[log2c_auto < -thr]) / arm_bp[1],
                  assessed_fraction = sum(overlap_bp) / arm_bp[1]),
              by = .(cell_line, subtype, arm)]
  al[, call := fifelse(assessed_fraction <= maj, "not assessed",
                       fifelse(gain_frac > maj, "gain", fifelse(loss_frac > maj, "loss", "neutral")))]
  al <- merge(al, famm[, .(cell_line, patient_id)], by = "cell_line")
  hgs <- al[subtype == "HGS"]
  nl <- uniqueN(hgs$cell_line); np <- uniqueN(hgs$patient_id)
  af <- hgs[, .(n_lines = nl, n_lines_gain = sum(call == "gain"),
                n_lines_loss = sum(call == "loss"), n_patients = np,
                n_patients_gain = uniqueN(patient_id[call == "gain"]),
                n_patients_loss = uniqueN(patient_id[call == "loss"])), by = arm]
  af[, `:=`(pct_lines_gain    = round(100 * n_lines_gain    / n_lines,    0),
            pct_lines_loss    = round(100 * n_lines_loss    / n_lines,    0),
            pct_patients_gain = round(100 * n_patients_gain / n_patients, 0),
            pct_patients_loss = round(100 * n_patients_loss / n_patients, 0),
            log2_threshold    = thr,
            arm_majority_frac = maj)]
  list(arm_line = al, arm_freq = af, n_lines = nl, n_patients = np)
}

# HGS focus (canonical arm landscape); denominators = 18 lines / 11 patients.
main <- arm_freq_at()
arm_line <- main$arm_line
readr::write_csv(as_tibble(arm_line), file.path(OUT, "wes_cnv_arm_calls.csv"))
arm_freq <- main$arm_freq
n_hgs_lines <- main$n_lines; n_hgs_pat <- main$n_patients
# Record the call parameters IN the deposited table (review revision) so the six
# quoted frequencies are never separable from the threshold that produced them.
arm_freq[, `:=`(centring = "per-sample autosome probe-weighted median (log2c_auto)",
                pct_per_patient = round(100 / n_patients, 1),
                arm_reference = "UCSC hg38 cytoBand acen; 0-based half-open; centromeres excluded",
                fraction_denominator = "full annotated non-centromeric arm length",
                patient_rule = "at least one derived line with the event")]
# order by the stronger of gain/loss at patient level
arm_freq[, top_pct := pmax(pct_patients_gain, pct_patients_loss)]
setorder(arm_freq, -top_pct)
readr::write_csv(as_tibble(arm_freq[, !"top_pct"]), file.path(OUT, "wes_cnv_arm_freq_patient.csv"))

cat(sprintf("\n=== HGS arm-event frequency: LINE (n=%d) vs PATIENT (n=%d) — top arms ===\n",
            n_hgs_lines, n_hgs_pat))
cat(sprintf("    CALL RULE: |log2c_auto| > %.2f over > %.0f%% of arm length; each patient = %.1f%% (n=%d)\n",
            ARM_LOG2_THRESH, 100 * ARM_MAJORITY, 100 / n_hgs_pat, n_hgs_pat))
show_arms <- c("3q","8q","19q","20q","17p","13q","10q","5q","4q","8p","22q","6q","16q")
print(as.data.frame(arm_freq[arm %in% show_arms][order(match(arm, show_arms)),
        .(arm, gain_lines = pct_lines_gain, gain_pat = pct_patients_gain,
          loss_lines = pct_lines_loss, loss_pat = pct_patients_loss)]), row.names = FALSE)
cat("=> Patient frequencies count at least one derived model with the event; this does not establish a shared patient-trunk event.\n")

# 5a. THRESHOLD SENSITIVITY (review revision) ---------------------------------
# Sweep the two call parameters over a reasonable range and report the six arms the
# manuscript quotes, so a referee can see how much of each frequency is the rule.
MS_ARMS <- c("3q","20q","17p","8q","13q","19q")
MS_DIR  <- c("3q"="gain","20q"="gain","17p"="loss","8q"="gain","13q"="loss","19q"="gain")
sens <- rbindlist(c(
  lapply(ARM_THRESH_SWEEP,   function(t) arm_freq_at(thr = t, maj = ARM_MAJORITY)$arm_freq),
  lapply(setdiff(ARM_MAJORITY_SWEEP, ARM_MAJORITY),
                             function(m) arm_freq_at(thr = ARM_LOG2_THRESH, maj = m)$arm_freq)))
sens <- sens[arm %in% MS_ARMS]
sens[, reported_direction := MS_DIR[arm]]
sens[, pct_patients_reported := fifelse(reported_direction == "gain",
                                        pct_patients_gain, pct_patients_loss)]
sens[, is_headline := log2_threshold == ARM_LOG2_THRESH & arm_majority_frac == ARM_MAJORITY]
setorder(sens, arm_majority_frac, log2_threshold, arm)
readr::write_csv(as_tibble(sens[, .(arm, reported_direction, log2_threshold, arm_majority_frac,
                                    n_patients, n_patients_gain, n_patients_loss,
                                    pct_patients_gain, pct_patients_loss,
                                    pct_patients_reported, n_lines, pct_lines_gain,
                                    pct_lines_loss, is_headline)]),
                file.path(OUT, "wes_cnv_arm_freq_sensitivity.csv"))
cat("\n=== Arm-call SENSITIVITY: the six quoted frequencies vs the call parameters ===\n")
print(as.data.frame(dcast(sens[arm_majority_frac == ARM_MAJORITY],
                          log2_threshold ~ arm, value.var = "pct_patients_reported")[
                          , c("log2_threshold", MS_ARMS), with = FALSE]), row.names = FALSE)
cat(sprintf("  (rows = |log2| threshold at majority %.0f%%; headline row = %.2f)\n",
            100 * ARM_MAJORITY, ARM_LOG2_THRESH))
print(as.data.frame(dcast(sens[log2_threshold == ARM_LOG2_THRESH],
                          arm_majority_frac ~ arm, value.var = "pct_patients_reported")[
                          , c("arm_majority_frac", MS_ARMS), with = FALSE]), row.names = FALSE)
cat(sprintf("  (rows = arm-majority fraction at |log2| %.2f; headline row = %.2f)\n",
            ARM_LOG2_THRESH, ARM_MAJORITY))

# 5b. Point-locus canonical-event recovery: LINE vs PATIENT -------------------
hgs_lines <- cns_tbl$cell_line[cns_tbl$subtype == "HGS"]
hgs_pat_map <- famm[cell_line %in% hgs_lines, .(cell_line, patient_id)]
event_recovery <- HGSC_EVENTS %>% rowwise() %>%
  mutate(lines_hit = list({
    reg <- segs[cell_line %in% hgs_lines & chromosome == chrom &
                start < pos_mb*1e6 & end >= pos_mb*1e6]
    hit <- if (dir == "gain") reg[log2c_auto >  0.20] else reg[log2c_auto < -0.20]
    unique(hit$cell_line)
  }),
  n_lines    = length(lines_hit),
  n_patients = n_distinct(hgs_pat_map$patient_id[hgs_pat_map$cell_line %in% lines_hit])) %>%
  ungroup() %>%
  mutate(pct_lines    = round(100 * n_lines    / length(hgs_lines), 0),
         pct_patients = round(100 * n_patients / n_distinct(hgs_pat_map$patient_id), 0)) %>%
  select(-lines_hit)
cat(sprintf("\n=== Canonical HGSC point-locus recovery: LINE (n=%d) vs PATIENT (n=%d) ===\n",
            length(hgs_lines), n_distinct(hgs_pat_map$patient_id)))
print(as.data.frame(event_recovery %>% select(label, dir, pct_lines, pct_patients,
                                               n_lines, n_patients)), row.names = FALSE)

# 6. Genome-wide CN landscape heatmap (AUTOSOMES; + patient-family track) ------
bins <- suppressWarnings(tileGenome(HG38[AUTOSOMES], tilewidth = BIN_WIDTH,
                                    cut.last.tile.in.chrom = TRUE))
bin_chr <- factor(as.character(seqnames(bins)), levels = AUTOSOMES)

bin_sample <- function(scell) {
  d  <- segs[cell_line == scell & chromosome %in% AUTOSOMES]
  gr <- GRanges(as.character(d$chromosome), IRanges(d$start + 1L, d$end))
  h  <- findOverlaps(bins, gr)
  ov <- pintersect(bins[queryHits(h)], gr[subjectHits(h)])
  w  <- width(ov); val <- d$log2c_auto[subjectHits(h)]; qh <- queryHits(h)
  num <- tapply(val * w, qh, sum); den <- tapply(w, qh, sum)
  out <- rep(NA_real_, length(bins))
  out[as.integer(names(num))] <- num / den
  out
}
mat <- suppressWarnings(
  vapply(cns_tbl$cell_line, bin_sample, numeric(length(bins))))     # bins x samples
mat <- t(mat)                                                        # samples x bins
rownames(mat) <- cns_tbl$cell_line

sub_cols <- c(HGS="#0072B2", CC="#E69F00", EC="#009E73", LGS="#56B4E9", MC="#CC79A7")
sub_cols <- sub_cols[names(sub_cols) %in% as.character(unique(cns_tbl$subtype))]
row_sub  <- factor(as.character(cns_tbl$subtype), levels = names(sub_cols))
n_lab    <- setNames(paste0(names(sub_cols), " (n=", table(row_sub)[names(sub_cols)], ")"),
                     names(sub_cols))

# Patient-family track: multi-line families coloured, single-line patients grey.
fam_lv   <- c("1369","2295","3121","3133","3291")
fam_cols <- c(setNames(c("#1B9E77","#7570B3","#E7298A","#66A61E","#A6761D"), fam_lv),
              singleton = "grey88")
fam_grp  <- ifelse(cns_tbl$n_lines_in_family > 1, cns_tbl$family, "singleton")

# blue(loss)-white-gain(red) diverging, colourblind-safe (ColorBrewer RdBu ends).
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#2166AC", "#F7F7F7", "#B2182B"))

# Map canonical-event loci to bin columns for labelled marks along the top.
ev <- HGSC_EVENTS %>%
  mutate(bin_idx = map2_int(chrom, pos_mb, function(cc, pp) {
    idx <- which(as.character(seqnames(bins)) == cc &
                 start(bins) <= pp*1e6 & end(bins) >= pp*1e6)
    if (length(idx)) idx[1] else NA_integer_ }))
ev <- ev %>% filter(!is.na(bin_idx))
top_anno <- HeatmapAnnotation(
  events = anno_mark(at = ev$bin_idx, labels = ev$label, side = "top",
                     labels_gp = gpar(fontsize = 7),
                     link_gp = gpar(col = ifelse(ev$dir == "gain", "#B2182B", "#2166AC"))),
  show_annotation_name = FALSE)
left_anno <- rowAnnotation(
  `Patient family` = fam_grp, subtype = row_sub,
  col = list(`Patient family` = fam_cols, subtype = sub_cols),
  annotation_legend_param = list(
    subtype = list(labels = n_lab),
    `Patient family` = list(at = c(fam_lv, "singleton"),
                            labels = c(paste0("family ", fam_lv), "single-line patient"))),
  annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(4, "mm"),
  show_annotation_name = TRUE)

ht <- Heatmap(mat, name = "log2\n(centred,\nautosome)", col = col_fun, na_col = "grey85",
   cluster_rows = FALSE, cluster_columns = FALSE,
   row_split = row_sub, row_title_rot = 0, row_title_gp = gpar(fontsize = 9, fontface = "bold"),
   column_split = bin_chr, cluster_column_slices = FALSE,
   column_title_gp = gpar(fontsize = 7), column_title_rot = 90,
   column_gap = unit(0.3, "mm"), border = TRUE,
   show_column_names = FALSE, row_names_gp = gpar(fontsize = 6.5),
   left_annotation = left_anno, top_annotation = top_anno,
   heatmap_legend_param = list(direction = "horizontal"))

cap <- paste0("CNVkit copy number (segment log2, per-sample AUTOSOME-median-centred; ",
              "chrX AND chrY dropped - chrX reference-sex effect cannot be excluded). ",
              "Rows grouped by patient family (left track); the four 3133 lines are ONE patient. ",
              "Reference = 5 UNMATCHED PUBLIC exomes (PRJNA339046) - same capture kit (author-confirmed); ",
              "SeqCap EZ Exome v3; target-only resegmentation. Median centring may bias the baseline in highly altered genomes. ",
              "n=", nrow(mat), " lines / ", n_distinct(cns_tbl$patient_id), " patients.")

render_cnv <- function(open_dev) {
  open_dev()
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       column_title = "OvCAN WES - genome-wide copy-number landscape (CNVkit, autosomes)",
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       padding = unit(c(16, 2, 2, 2), "mm"))
  grid::grid.text(cap, x = 0.5, y = 0.02,
                  gp = grid::gpar(fontsize = 6, fontface = "italic", col = "grey25"))
  dev.off()
}
render_cnv(function() pdf(file.path(FIGS, "f_wes_cnv.pdf"), width = 13, height = 8))
render_cnv(function() ragg::agg_png(file.path(PROJ, "reports", "assets", "f_wes_cnv.png"),
                                    width = 13, height = 8, units = "in", res = 200))
render_cnv(function() pdf(file.path(FIGS, "08_cnv_landscape.pdf"), width = 13, height = 8))  # legacy name
message("Wrote figs/f_wes_cnv.pdf + reports/assets/f_wes_cnv.png (and legacy figs/08_cnv_landscape.pdf)")

message("\n08_wes_cnv.R complete. Outputs: output/wes_cnv_segments.csv (adds log2c_auto), ",
        "output/wes_cnv_fga.csv (augmented: autosome FGA + patient_id + chrX diagnostic), ",
        "output/wes_cnv_arm_freq_patient.csv (now carries the call threshold), ",
        "output/wes_cnv_arm_freq_sensitivity.csv (NEW: threshold sweep) ; ",
        "figure figs/f_wes_cnv.pdf (+ reports/assets PNG).")
