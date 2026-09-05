# =============================================================================
# Script: 09_wes_hrd.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: FEASIBILITY GATE for genomic HRD (scarHRD) from the archived WES.
#          Determine — programmatically, from the files — whether genuine genomic
#          HRD (HRD-LOH + TAI + LST) can be computed. It CANNOT: scarHRD needs
#          allele-specific copy number, which the archived CNVkit outputs do not
#          contain. We DO NOT fabricate a score; we document the blocker and the
#          path to a real HRD score. Phase 3 (WES), step 3 of 3.
# Author:  Cook Lab (analyst: Claude)
# Date:    2026-07-23
# =============================================================================
#
# WHY genomic HRD needs allele-specific CN:
#   The genomic-scar HRD score = HRD-LOH (Abkevich 2012) + telomeric allelic
#   imbalance (Birkbak 2012) + large-scale state transitions (Popova 2012). The composite requires allele-specific
#   information for LOH and TAI; an LST-like count alone is not the composite HRD score. scarHRD (Sztupinszki 2018) therefore requires
#   per-segment ALLELE-SPECIFIC copy number (major "A_cn" + minor "B_cn", i.e. a
#   B-allele-frequency-derived split), typically from Sequenza or ASCAT.
#   TOTAL copy number / log2 (what plain CNVkit produces without a SNP VCF) is
#   NOT sufficient — you cannot see LOH in total CN alone.
#
#   NB: HRDetect (Davies 2017) and CHORD (Nguyen 2020) are WGS-only (they use
#   structural-variant signatures) and do not apply to these exomes. The Peng
#   expression signature is a transcriptional state, NOT a genomic scar — it is
#   dropped (category error), per the analysis plan and literature review Theme 8.
# =============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })

CNV_DIR <- file.path(DATA, "cnvkit wes - new")
HRD_MD  <- file.path(OUT, "wes_hrd_feasibility.md")
PARAM_CSV <- file.path(OUT, "wes_pipeline_parameters.csv")

# 1. Data check: do the CNVkit outputs carry allele-specific / BAF info? -------
# Allele-specific evidence would be columns like baf / cn1 / cn2 / A_cn / B_cn in
# the .cnr/.cns/.call.cns, and a `--vcf`/`-v` SNP VCF in the cnvkit batch call.
allele_cols <- c("baf", "cn1", "cn2", "A_cn", "B_cn", "minor_cn", "major_cn")
hdr_of <- function(f) if (length(f) && !is.na(f[1]))
  strsplit(readLines(f[1], n = 1), "\t")[[1]] else character(0)
all_cns <- list.files(CNV_DIR, pattern = "\\.cns$", recursive = TRUE, full.names = TRUE)
hdr_cnr  <- hdr_of(list.files(CNV_DIR, pattern = "\\.cnr$", recursive = TRUE, full.names = TRUE))
hdr_cns  <- hdr_of(all_cns[!grepl("\\.(call|bintest)\\.cns$", all_cns)])
hdr_call <- hdr_of(all_cns[grepl("\\.call\\.cns$", all_cns)])
found_allele <- intersect(allele_cols, unique(c(hdr_cnr, hdr_cns, hdr_call)))

# CNVkit batch calls: was a SNP VCF (`--vcf`/`-v`) passed (the BAF input)?
cmds <- readLines(file.path(CNV_DIR, "commands.txt"))
batch_calls <- grep("cnvkit\\.py batch", cmds, value = TRUE)
n_with_vcf  <- sum(grepl("(--vcf|(^|\\s)-v(\\s|$))", batch_calls))

# Tumor recal BAMs archived? (a Sequenza/ASCAT re-run would need them)
tumor_bams  <- list.files(DATA, pattern = "recal\\.bam$",
                          recursive = TRUE, full.names = TRUE)
any_bams    <- list.files(DATA, pattern = "\\.bam$",
                          recursive = TRUE, full.names = TRUE)

cat("\n=== HRD feasibility — data checks ===\n")
cat("Allele-specific columns found in CNVkit .cnr/.cns/.call.cns:",
    if (length(found_allele)) paste(found_allele, collapse = ", ") else "NONE", "\n")
cat(sprintf(".call.cns columns: %s\n", paste(hdr_call, collapse = ", ")))
cat(sprintf("CNVkit batch calls: %d total, %d passed a SNP --vcf (BAF input)\n",
            length(batch_calls), n_with_vcf))
cat(sprintf("Tumor recal BAMs archived: %d (any .bam: %d -> %s)\n",
            length(tumor_bams), length(any_bams),
            if (length(any_bams)) paste(basename(any_bams), collapse = ", ") else "none"))

allele_specific_available <- any(vapply(list(c("cn1", "cn2"), c("A_cn", "B_cn"),
  c("major_cn", "minor_cn")), function(pair) all(pair %in% found_allele), logical(1)))
# A --vcf argument or BAF field alone is not a derived major/minor CN solution.
FEASIBLE <- allele_specific_available   # dispositive data gate

# 1b. RECOVERABLE PIPELINE PARAMETERS -> machine-readable table  [review revision]
# -----------------------------------------------------------------------------
# The referee's M2 asks for versions and key parameters for every tool named. The
# CNVkit invocation IS fully recoverable — it is in the archived `commands.txt`
# this script already reads — but it existed only as a comment in 08's header, so
# Methods could not cite it. Parse it here and write output/wes_pipeline_parameters.csv
# so every value is traceable to the file it came from, and so what is NOT
# recoverable is recorded as explicitly NOT recoverable rather than left blank.
grab1 <- function(pattern, x, group = 1L) {
  m <- regmatches(x, regexec(pattern, x))
  v <- unique(vapply(m, function(z) if (length(z) >= group + 1L) z[group + 1L] else NA_character_,
                     character(1)))
  v <- v[!is.na(v)]
  if (length(v)) paste(v, collapse = " | ") else NA_character_
}
cnvkit_ver  <- grab1("cnvkit-([0-9.]+)--", batch_calls)
cnvkit_img  <- grab1("(cnvkit-[^ ]+\\.img)", batch_calls)
ref_fasta   <- grab1("-f ([^ ]+)", batch_calls)
target_bed  <- grab1("-t ([^ ]+)", batch_calls)
antitgt_bed <- grab1("-a ([^ ]+)", batch_calls)
pon_glob    <- grab1("-n ([^ ]+)", batch_calls)
# analysis-affecting flags only: drop --output-dir (bookkeeping) and the "--pyhdfd"
# fragment the container image name contributes to a naive `--[a-z-]+` match.
extra_flags <- paste(setdiff(sort(unique(unlist(regmatches(batch_calls,
                       gregexpr("--[a-z][a-z-]+", batch_calls))))),
                     c("--output-dir", "--pyhdfd")), collapse = " ")
n_batch     <- length(batch_calls)

params <- data.frame(
  step = c(rep("CNV calling (CNVkit)", 8), rep("WES capture", 2),
           rep("CNV reference (panel of normals)", 2), "SNV calling (Mutect2)",
           "Genome build", "Proteomics search"),
  parameter = c("tool", "version", "container image", "reference genome FASTA",
                "target intervals", "antitarget intervals", "flags", "n samples processed",
                "capture kit identity", "capture interval file",
                "normals used", "kit-matched to tumours?",
                "matched normal", "build actually used", "search parameters"),
  value = c("cnvkit.py batch", cnvkit_ver, cnvkit_img, ref_fasta,
            target_bed, antitgt_bed, extra_flags, as.character(n_batch),
            "NOT RECOVERABLE", basename(ifelse(is.na(target_bed), "", target_bed)),
            pon_glob, "SAME CAPTURE KIT (author-confirmed, v5 comment 31); kit name/design version pending",
            "NONE (tumour-only; n_ref/n_alt columns empty in every MAF)",
            "GRCh38/hg38 (MAF NCBI_Build field reads 'GRCh37' — a spurious vcf2maf default; see 16 header)",
            "NOT RECOVERABLE"),
  recoverable = c(rep(TRUE, 8), FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
  evidence = c(rep("judy_archive/data/cnvkit wes - new/commands.txt", 8),
               "no capture-kit record in the archive; intervals_sorted.bed itself is not archived",
               "judy_archive/data/cnvkit wes - new/commands.txt",
               "judy_archive/data/cnvkit wes - new/commands.txt (5 public exomes, PRJNA339046)",
               "author v5 comment 31 confirms the same capture kit; documentary kit name/design version pending",
               "MAF n_depth/n_ref_count/n_alt_count empty for every variant",
               "Mutect2 PoN 1000g_pon.hg38.vcf.gz; VCF ##contig chr1 length 248,956,422; hotspot coordinates",
               "judy_archive/data/proteomics/Readme.md is a 1-byte EMPTY file"),
  action = c(rep("cite in Methods", 8),
             "REQUIRES PI: ask which exome capture kit was used",
             "cite in Methods", "cite in Methods",
             "REQUIRES PI: supply the shared capture-kit name and design version",
             "state as a limitation", "state in Methods (and flag the MAF header as a deposition item)",
             "REQUIRES PI / Morin lab: instrument, MS2-vs-MS3, search engine + version, database, enzyme, modifications, tolerances, FDR, interference filter, peptide-to-protein rollup"),
  stringsAsFactors = FALSE)
# Archived per-normal coverage names and derived target bins provide more
# provenance than the batch command's glob and unarchived vendor BED path.
normal_files <- list.files(CNV_DIR, pattern = "^SRR[0-9]+.*targetcoverage\\.cnn$",
                           recursive = TRUE, full.names = TRUE)
normal_ids <- sort(unique(sub("\\..*$", "", basename(normal_files))))
target_files <- list.files(CNV_DIR, pattern = "^intervals_sorted\\.target\\.bed$",
                           recursive = TRUE, full.names = TRUE)
target_hash <- unique(vapply(target_files, function(p) digest::digest(file = p, algo = "sha256"), ""))
stopifnot("Archived target-bin files disagree" = length(target_hash) == 1L)
file.copy(target_files[1], file.path(OUT, "wes_cnvkit_target_intervals.bed"), overwrite = TRUE)
params$value[params$parameter == "normals used"] <- paste(normal_ids, collapse = "; ")
params$evidence[params$parameter == "normals used"] <- "archived per-normal *.targetcoverage.cnn / *.antitargetcoverage.cnn filenames"
params$recoverable[params$parameter == "kit-matched to tumours?"] <- TRUE
params <- bind_rows(params, tibble(
  step = "CNV target footprint", parameter = "derived target bins",
  value = paste0("wes_cnvkit_target_intervals.bed; ", length(target_files), " identical archived copies; sha256 ", target_hash),
  recoverable = TRUE, evidence = sub(paste0("^", PROJ, "/"), "", target_files[1]),
  action = "deposit derived bins; recover original vendor BED and kit design/version separately"))

# GATK and VEP versions are present in primary VCF headers for every model.
source("scripts/lib/wes_maf_inputs.R")
inputs <- wes_maf_inputs(PROJ, DATA, OUT, recover = FALSE)
header_provenance <- map_dfr(seq_len(nrow(inputs)), function(i) {
  f <- inputs$source_vcf[i]
  con <- if (grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "rt")
  h <- readLines(con); close(con); h <- h[grepl("^##", h)]
  gatk <- h[grepl("^##GATKCommandLine=<ID=(Mutect2|FilterMutectCalls),", h)]
  vep <- h[grepl("^##VEP=", h)]
  tibble(cell_line = inputs$cell_line[i], source_vcf = sub(paste0("^", PROJ, "/"), "", f),
    mutect2_version = grab1('Version="([^"]+)"', gatk[grepl("ID=Mutect2,", gatk)]),
    filtermutectcalls_version = grab1('Version="([^"]+)"', gatk[grepl("ID=FilterMutectCalls,", gatk)]),
    vep_version = grab1('^##VEP="([^"]+)"', vep),
    vep_assembly = grab1('assembly="([^"]+)"', vep),
    gnomad_exomes = grab1('gnomADe="([^"]+)"', vep),
    gnomad_genomes = grab1('gnomADg="([^"]+)"', vep),
    clinvar_release = grab1('ClinVar="([^"]+)"', vep),
    dbsnp_release = grab1('dbSNP="([^"]+)"', vep))
})
readr::write_csv(header_provenance, file.path(OUT, "wes_vcf_header_provenance.csv"))
params <- bind_rows(params, tibble(step = c("SNV calling (Mutect2)", "Variant filtering", "Variant annotation"),
  parameter = c("Mutect2 version", "FilterMutectCalls version", "Ensembl VEP version"),
  value = c(paste(unique(header_provenance$mutect2_version), collapse = ";"),
            paste(unique(header_provenance$filtermutectcalls_version), collapse = ";"),
            paste(unique(header_provenance$vep_version), collapse = ";")),
  recoverable = TRUE, evidence = "all primary VCF headers; wes_vcf_header_provenance.csv",
  action = "cite in Methods"))
readr::write_csv(params, PARAM_CSV)
cat("\n=== Recoverable WES pipeline parameters (-> wes_pipeline_parameters.csv) ===\n")
print(as.data.frame(params[, c("step","parameter","value","recoverable")]), row.names = FALSE)
cat(sprintf("\n%d of %d parameters recoverable from the archive; %d REQUIRE the PI: %s\n",
            sum(params$recoverable), nrow(params), sum(!params$recoverable),
            paste(params$parameter[!params$recoverable], collapse = ", ")))

# 2. Read-only tooling inventory. Installing software cannot remedy absent
# allele-specific data; workflow execution must not mutate the user's library.
scarhrd_installed <- requireNamespace("scarHRD", quietly = TRUE)
copynumber_avail <- requireNamespace("copynumber", quietly = TRUE)
install_note <- "not attempted: software installation cannot resolve the missing allele-specific input"

# 3. Verdict + feasibility report (NO fabricated score) ------------------------
verdict <- if (FEASIBLE) "FEASIBLE" else "NOT FEASIBLE FROM ARCHIVED DATA"
cat(sprintf("\n=== GENOMIC HRD FEASIBILITY VERDICT: %s ===\n", verdict))
if (!FEASIBLE)
  cat("Reason: no allele-specific copy number (no BAF split); scarHRD LOH/TAI/LST\n",
      "cannot be computed from total-CN CNVkit output. No score is produced.\n")

md <- c(
"# Genomic HRD — feasibility assessment (WES)",
"",
paste0("**Date:** 2026-07-23  |  **Script:** `scripts/09_wes_hrd.R`  |  **Project:** OvCAN human ovarian cancer cell-line resource"),
"",
paste0("## Verdict: **", verdict, "**"),
"",
"Genuine genomic HRD (the genomic-scar score = HRD-LOH + telomeric allelic imbalance + large-scale state transitions) **cannot be computed from the archived WES data**, because the copy-number calls are **total copy number only — there is no allele-specific / B-allele-frequency information**. We therefore report **no HRD score** rather than a fabricated one.",
"",
"## What genomic HRD requires",
"- The composite HRD-scar score combines LOH, TAI, and LST; LOH and TAI require allele-specific information, so an LST-like count from total-copy profiles cannot substitute for the composite score.",
"- `scarHRD` (Sztupinszki 2018, npj Breast Cancer; HRD-sum r=0.87 vs SNP array, robust at 30x) consumes **per-segment allele-specific copy number** (major `A_cn` + minor `B_cn`), normally produced by **Sequenza** or **ASCAT** from a BAF-bearing SNP VCF.",
"- **Total copy number / log2 (plain CNVkit) is insufficient:** LOH is invisible in total CN. HRDetect (Davies 2017) and CHORD (Nguyen 2020) are **WGS-only** (SV-signature based) and do not apply to exomes.",
"",
"## Evidence from the archived files (this script)",
paste0("- Allele-specific columns (`baf`/`cn1`/`cn2`/`A_cn`/`B_cn`) in CNVkit `.cnr`/`.cns`/`.call.cns`: **",
       if (length(found_allele)) paste(found_allele, collapse=", ") else "NONE", "**."),
paste0("- `.call.cns` columns are total-CN only: `", paste(hdr_call, collapse="`, `"), "` (has `cn`/`log2`, no allele split)."),
paste0("- CNVkit `batch` calls passing a SNP `--vcf` (the BAF input that would enable allele-specific segmentation): **",
       n_with_vcf, " of ", length(batch_calls), "** — none. CNVkit was run with `--diagram --scatter` but **no `--vcf`**."),
paste0("- Tumor `recal.bam`s archived (needed to re-run Sequenza/ASCAT): **", length(tumor_bams),
       "**. The only archived BAM is a single *public normal* (",
       if (length(any_bams)) paste(basename(any_bams), collapse=", ") else "none",
       "; PRJNA339046). The tumors' recal BAMs live on HPC scratch (`/scratch/asmab/...`), not in the archive."),
"",
"## scarHRD tooling note (secondary)",
paste0("- `scarHRD` installed in this environment: **", scarhrd_installed, "**."),
paste0("- Its dependency **`copynumber`** is installed locally (Bioconductor ", as.character(BiocManager::version()),
       "): **", copynumber_avail, "**."),
if (!is.na(install_note)) paste0("- Install attempt outcome: `", install_note, "`.") else "- scarHRD already present.",
"- This is a *secondary* obstacle only; even with scarHRD installed, the **data blocker above is dispositive**.",
"",
"## Recommended path to a genuine genomic HRD score",
"1. Recover the tumors' recalibrated BAMs (nf-core/Sarek `recal.bam`, HPC scratch) — the CNV/SNV inputs.",
"2. Generate **allele-specific CN** per line with **Sequenza** (or ASCAT/FACETS): call heterozygous SNP BAFs (a population SNP panel suffices for tumor-only; a matched normal is better) + depth ratio, fit cellularity/ploidy, emit segments with `A_cn`/`B_cn`.",
"3. Run **`scarHRD`** on the Sequenza segments to obtain HRD-LOH + TAI + LST and the HRD-sum.",
"4. Interpretation caveat (Takamatsu 2024, *Sci Data*): in cell lines, HRD scars **persist** and do **not** predict in-vitro platinum/PARP sensitivity — report HRD as a genomic-scar descriptor, not a drug-response predictor.",
"",
"## Note on the dropped expression 'HRD'",
"The Peng et al. 2014 expression signature previously labelled 'HRD' is a **transcriptional state, not a genomic scar** (published transcriptomic HRD signatures barely overlap and none is guideline-endorsed). It is **dropped**; 'HRD' in this resource must mean the genomic-scar score above.",
"")
writeLines(md, HRD_MD)
message("\n09_wes_hrd.R complete. Wrote feasibility report -> ", HRD_MD,
        "\n  and machine-readable pipeline parameters -> ", PARAM_CSV,
        "\n  Verdict: ", verdict, " (no HRD score fabricated).")
