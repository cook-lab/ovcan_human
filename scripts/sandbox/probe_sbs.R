suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(MutationalPatterns)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
ref <- "BSgenome.Hsapiens.UCSC.hg38"
std <- paste0("chr", c(1:22, "X", "Y"))
read_pass_snv <- function(path) {
  r <- fread(path, sep = "\t", skip = "Hugo_Symbol", quote = "",
             na.strings = c("", "NA", "."), showProgress = FALSE)
  afc <- intersect(c("gnomADe_AF", "AF", "AA_AF", "EA_AF"), names(r))
  n0 <- function(v) { x <- suppressWarnings(as.numeric(v)); ifelse(is.na(x), 0, x) }
  r[, pop_af_max := do.call(pmax, lapply(.SD, n0)), .SDcols = afc]
  r <- r[FILTER == "PASS" & pop_af_max <= 0.001 & Variant_Type == "SNP" &
           Reference_Allele %in% c("A", "C", "G", "T") &
           Tumor_Seq_Allele2 %in% c("A", "C", "G", "T") &
           Chromosome %in% std]
  gr <- GRanges(seqnames = r$Chromosome,
                ranges = IRanges(start = r$Start_Position, width = 1),
                REF = Biostrings::DNAStringSet(r$Reference_Allele),
                ALT = Biostrings::DNAStringSetList(as.list(r$Tumor_Seq_Allele2)))
  GenomeInfoDb::seqlevels(gr) <- std
  GenomeInfoDb::genome(gr) <- "hg38"
  gr
}
f21 <- list.files("judy_archive/data/wes - old/mutect2", pattern = "TOV21G.*\\.maf$",
                  recursive = TRUE, full.names = TRUE)
gr <- read_pass_snv(f21)
cat("TOV21G PASS SNVs:", length(gr), "\n")
grl <- GRangesList(TOV21G = gr)
mm <- mut_matrix(vcf_list = grl, ref_genome = ref)
cat("mut_matrix dim:", paste(dim(mm), collapse = " x "), " colsum:", sum(mm), "\n")
types <- sub("^.\\[(.*)\\].$", "\\1", rownames(mm))
cat("6-class collapse:\n"); print(tapply(mm[, 1], types, sum))
sig <- tryCatch(get_known_signatures(muttype = "snv", source = "COSMIC_v3.2", genome = "GRCh38"),
                error = function(e) { cat("get_known_signatures ERR:", conditionMessage(e), "\n"); NULL })
if (!is.null(sig)) {
  cat("known sig matrix dim:", paste(dim(sig), collapse = " x "), "\n")
  cs <- cos_sim_matrix(mm, sig)
  cat("Top 10 COSMIC cosine matches for TOV21G:\n")
  print(round(sort(cs[1, ], decreasing = TRUE)[1:10], 3))
}
