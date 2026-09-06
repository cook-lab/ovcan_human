# =============================================================================
# Script: 09_wes_hrd.R
# Project: OvCAN human ovarian cancer cell-line multi-omic resource
# Purpose: Gate direct genomic-scar scoring from existing CNV columns.
#          Current relative CNV profiles lack a validated major/minor solution.
#          VCF allele counts exist, but require a separate input-eligibility
#          audit and allele-specific fit before scarHRD. No score is produced.
#          See reports/clinical_classification_2026-09-06/HRD_FEASIBILITY.md.
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
#   Standard HRDetect and CHORD models require genome-wide mutation/SV features
#   absent here. HRDetect also has a separately evaluated WES model; the original
#   WGS model cannot be applied with absent features replaced by zero.
#   Expression signatures describe a transcriptional association, not a scar or
#   current functional HR assay. This script does not compute either.
#
#   Historical provenance parser below predates the September recovery. Its
#   original archive scan is not a complete current inventory. Scripts 23-26 and
#   the dated completion/coverage reports supersede old unrecoverable labels.
#   Text revised 2026-09-06 without executing this provenance-writing script.
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

# Historical original-archive BAM inventory; not a cluster existence check.
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
FEASIBLE <- allele_specific_available   # column-presence gate only; any fit still needs validation

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

# 2. Read-only tooling inventory; do not install packages as part of this gate.
# Future allele-specific fitting requires separately verified input eligibility.
scarhrd_installed <- requireNamespace("scarHRD", quietly = TRUE)
copynumber_avail <- requireNamespace("copynumber", quietly = TRUE)
install_note <- "not attempted: allele-specific input preparation/fitting is a separate analysis"

# 3. Verdict + feasibility report (NO fabricated score) ------------------------
verdict <- if (FEASIBLE) "ALLELE-SPECIFIC COLUMNS PRESENT; FIT VALIDATION REQUIRED" else "NO VALIDATED ALLELE-SPECIFIC SOLUTION FOR DIRECT SCAR SCORING"
cat(sprintf("\n=== GENOMIC HRD FEASIBILITY VERDICT: %s ===\n", verdict))
if (!FEASIBLE)
  cat("Reason: no allele-specific copy number (no BAF split); scarHRD LOH/TAI/LST\n",
      "cannot be computed from total-CN CNVkit output. No score is produced.\n")

md <- c(
"# Genomic-scar HRD — feasibility assessment (WES)",
"",
"**Assessment revised:** 2026-09-06. Original archive inventory: 2026-07-23; script `scripts/09_wes_hrd.R`.",
"",
paste0("## Current gate: **", verdict, "**"),
"",
"No composite genomic-scar score is produced. Current CNV profiles contain relative total-copy ratios, without a validated major/minor copy-number and ploidy solution. Those profiles cannot directly supply LOH/TAI/LST-based scarHRD input. This is not a claim that WES-based HRD analysis is inherently impossible or that the VCFs lack allele counts.",
"",
"The [September HRD review](../reports/clinical_classification_2026-09-06/HRD_FEASIBILITY.md) supersedes the earlier blanket feasibility verdict. It verifies AD/AF in all 23 VCFs, audits germline-site ascertainment, identifies an exact pathogenic BRCA2 allele, and separates genotype, scars, mutational signatures and current function. The [targeted next-step request](../docs/cluster/CLINICAL_CLASSIFICATION_NEXT_STEPS.md) defines the proposed cluster work.",
"",
"## What the current evidence supports",
"",
"- scarHRD needs allele-specific segments and ploidy. Total-copy transitions or FGA alone are insufficient; neither equals the composite scar score.",
"- The original CNVkit column scan below is historical. Current target-only CNS files and their source CNR hashes are in `output/wes_cnv_target_only/manifest.csv`; the old antitarget-containing profiles are superseded.",
"- Existing VCF FORMAT AD/AF can support variant-level review. All 23 caller headers nevertheless record `--genotype-germline-sites false`, `--genotype-pon-sites false`, and zero interval padding. PASS-only filtering removes most retained common SNPs. Neither the coding MAF nor the present VCF is an automatically eligible germline SNP dataset for PureCN.",
"- PureCN supports tumour-only exome and cell-line allele-specific inference, including CNVkit inputs, once an appropriate SNP set is supplied. Standard Sequenza/FACETS workflows use matched tumour-normal allelic evidence; unrelated reference exomes are not patient-matched normals.",
"- Standard HRDetect and CHORD models require mutation/SV features not supplied here. A separately trained WES HRDetect model was evaluated in the original study; it is incorrect to call every HRDetect application WGS-only. No HRDetect/CHORD result is generated here.",
"",
"## Historical original-archive scan (scope limited to this script)",
"",
paste0("- Allele-specific columns in scanned original CNVkit files: ",
       if (length(found_allele)) paste(found_allele, collapse=", ") else "NONE", "."),
paste0("- Original `.call.cns` columns: `", paste(hdr_call, collapse="`, `"), "`. The `cn` column alone does not establish ploidy or an allele split."),
paste0("- Recorded original batch calls with a SNP `--vcf`: ", n_with_vcf, " of ", length(batch_calls), "."),
paste0("- Files ending `recal.bam` in the original DATA tree: ", length(tumor_bams),
       "; any BAM: ", length(any_bams), ". This does not establish current cluster availability or identity of later CRAM files."),
"- Upstream parameters generated by this historical parser are not the authoritative September provenance assessment; consult scripts 23-26 and the dated WES completion/coverage reports.",
"",
"## Tooling and proposed path",
"",
paste0("- scarHRD installed at this scan: ", scarhrd_installed, "; copynumber available: ", copynumber_avail, "."),
paste0("- Installation: ", install_note, "."),
"- First retrieve existing germline-inclusive allele counts and normal SNP mapping-bias evidence. If absent, obtain verified alignment/reference identities and authorize a bounded input-generation and PureCN pilot separately. Five-normal coverage pairs are already recovered; do not request those again.",
"- Compare plausible purity/ploidy solutions, inspect BAF/log-ratio fits and gene-locus allelic states, and retain unresolved fits as unassessed. Only accepted allele-specific segments can proceed to exploratory LOH/TAI/LST scoring. Do not transfer a commercial clinical cutoff to cell-line scar sums without validation.",
"- Genomic scars can persist after restoration of repair. They are not direct measurements of current HR function or guaranteed drug response. RAD51-foci assays and measured drug sensitivity address different dimensions; expression abundance is not their substitute.",
"",
"The previously dropped expression signature remains excluded from the genomic-scar output. HRD is a broader biological concept: an expression association, a genomic scar and a functional assay should be named and validated according to what each measures.",
"")
writeLines(md, HRD_MD)
message("\n09_wes_hrd.R complete. Wrote feasibility report -> ", HRD_MD,
        "\n  and machine-readable pipeline parameters -> ", PARAM_CSV,
        "\n  Verdict: ", verdict, " (no HRD score fabricated).")
