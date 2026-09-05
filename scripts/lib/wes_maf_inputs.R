# Canonical WES MAF inventory. Recover annotation-only omissions from existing
# VEP-annotated VCFs; preserve the raw archive and record source checksums.
wes_maf_inputs <- function(PROJ, DATA, OUT, recover = TRUE) {
  mut_dir <- file.path(DATA, "wes - old", "mutect2")
  key <- function(x) toupper(gsub("[-_]", "", sub("_[Pp][0-9]+$", "", x)))
  samples <- read.csv(file.path(PROJ, "metadata", "samples.csv"), check.names = FALSE)
  wanted <- samples[samples$provenance == "generated" & samples$wes_mut == "Y", ]
  maf <- list.files(mut_dir, pattern = "\\.maf$", recursive = TRUE, full.names = TRUE)
  archive_key <- key(basename(dirname(maf)))
  out <- data.frame(cell_line = wanted$cell_line, source_maf = maf[match(key(wanted$cell_line), archive_key)],
                    source_kind = "archived MAF", source_vcf = NA_character_)
  for (i in seq_len(nrow(out))) {
    dirs <- list.dirs(mut_dir, recursive = FALSE, full.names = TRUE)
    d <- dirs[key(basename(dirs)) == key(out$cell_line[i])]
    stopifnot(length(d) == 1L)
    vcf <- list.files(d, pattern = "filtered_VEP\\.ann\\.vcf(\\.gz)?$", full.names = TRUE)
    vcf <- vcf[order(grepl("\\.gz$", vcf))]
    if (length(vcf)) out$source_vcf[i] <- vcf[1]
    if (!is.na(out$source_maf[i])) next
    stopifnot("A missing MAF has no annotated source VCF" = length(vcf) >= 1L)
    recovered <- file.path(OUT, "recovered_maf", paste0(basename(d), ".maf"))
    if (!file.exists(recovered)) {
      stopifnot("Recover the missing MAF first" = recover)
      dir.create(dirname(recovered), recursive = TRUE, showWarnings = FALSE)
      tmp_vcf <- tempfile(pattern = "ovcan_recovery_", fileext = ".vcf")
      con <- if (grepl("\\.gz$", vcf[1])) gzfile(vcf[1], "rt") else file(vcf[1], "rt")
      writeLines(readLines(con), tmp_vcf); close(con)
      ref <- system.file("extdata", "single_sequences.2bit", package = "BSgenome.Hsapiens.UCSC.hg38")
      stopifnot("Installed hg38 TwoBit reference required for recovery" = file.exists(ref))
      adapter <- file.path(PROJ, "scripts", "lib", "twobit_faidx.R")
      # No VEP rerun: preserve the archived consequence annotations. The MAF
      # barcode uses the historic convention, matching the other archived MAFs.
      args <- c(file.path(PROJ, "scripts", "vendor", "vcf2maf-1.6.22.pl"),
        "--input-vcf", tmp_vcf, "--output-maf", recovered,
        "--tumor-id", basename(d), "--inhibit-vep", "--ncbi-build", "GRCh38",
        "--ref-fasta", ref, "--samtools-exec", adapter, "--tabix-exec", adapter)
      status <- system2("perl", shQuote(args), stdout = paste0(recovered, ".conversion.log"),
                        stderr = paste0(recovered, ".conversion.log"))
      unlink(tmp_vcf)
      stopifnot("vcf2maf recovery failed" = status == 0L, file.exists(recovered))
    }
    out$source_maf[i] <- recovered
    out$source_kind[i] <- "recovered MAF from archived VEP VCF; vcf2maf 1.6.22 --inhibit-vep"
  }
  stopifnot(!anyNA(out$source_maf), !anyDuplicated(out$cell_line))
  out$maf_sha256 <- vapply(out$source_maf, function(f) digest::digest(file = f, algo = "sha256"), "")
  out$vcf_sha256 <- vapply(out$source_vcf, function(f) if (is.na(f)) NA_character_ else
                          digest::digest(file = f, algo = "sha256"), "")
  out
}
