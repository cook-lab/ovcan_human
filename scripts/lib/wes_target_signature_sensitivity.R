# Exploratory capture-footprint sensitivity, sourced by script22 after the main
# unnormalised refit. Archived CNVkit-derived target bins are a recoverable proxy
# for the sequenced territory, not a per-sample callable mask or vendor manifest.
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
source("scripts/lib/wes_maf_inputs.R")
bed_files <- list.files(file.path(DATA, "cnvkit wes - new"),
  pattern = "^intervals_sorted\\.target\\.bed$", recursive = TRUE, full.names = TRUE)
bed_sha <- vapply(bed_files, function(p) digest::digest(file = p, algo = "sha256"), "")
stopifnot(length(bed_files) > 0L, length(unique(bed_sha)) == 1L,
  unique(bed_sha) == "2177970ffdbc933e488067f812e3ae4760310bd9a40fadf4e2fa5389480d955d")
bed <- fread(bed_files[1], select = 1:3, col.names = c("chr", "start", "end"))
std_chr <- paste0("chr", c(1:22, "X", "Y"))
bed <- bed[chr %in% std_chr]
territory <- reduce(GRanges(bed$chr, IRanges(bed$start + 1L, bed$end)))
ref <- BSgenome.Hsapiens.UCSC.hg38
# Context opportunities include one flank either side of each target union, so
# every targeted central base is counted once, including interval endpoints.
tri_cache <- file.path(OUT, "wes_target_trinucleotide_opportunities.csv")
if (file.exists(tri_cache)) {
  tri <- read.csv(tri_cache)
  tri_ok <- all(tri$target_bed_sha256 == unique(bed_sha)) && nrow(tri) == 32L
} else tri_ok <- FALSE
if (!tri_ok) {
  message("Counting hg38 and archived target-footprint trinucleotide opportunities...")
  genome64 <- Reduce(`+`, lapply(std_chr, function(ch) oligonucleotideFrequency(ref[[ch]], width = 3L)))
  flank <- resize(territory, width(territory) + 2L, fix = "center")
  # Target intervals are internal to the reference chromosomes in this archive.
  stopifnot(all(start(flank) > 0L), all(end(flank) <= seqlengths(ref)[as.character(seqnames(flank))]))
  target64 <- colSums(oligonucleotideFrequency(getSeq(ref, flank), width = 3L))
  pyr <- names(genome64)[substr(names(genome64), 2, 2) %in% c("C", "T")]
  rc <- as.character(reverseComplement(DNAStringSet(pyr)))
  g <- unname(genome64[pyr] + genome64[rc]); t <- unname(target64[pyr] + target64[rc])
  tri <- data.frame(trinucleotide = pyr, genome_opportunities = g,
    target_opportunities = t, genome_fraction = g / sum(g), target_fraction = t / sum(t),
    target_bed_sha256 = unique(bed_sha), target_union_bp = sum(width(territory)))
  readr::write_csv(tri, tri_cache)
}

inputs <- wes_maf_inputs(PROJ, DATA, OUT, recover = FALSE)
one_target_sbs <- function(path) {
  d <- fread(path, skip = "Hugo_Symbol", quote = "", na.strings = c("", "NA", "."),
    select = c("FILTER", "gnomADe_AF", "AF", "AA_AF", "EA_AF", "Variant_Type",
               "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
  numeric0 <- function(z) { z <- suppressWarnings(as.numeric(z)); z[is.na(z)] <- 0; z }
  af <- do.call(pmax, lapply(d[, c("gnomADe_AF", "AF", "AA_AF", "EA_AF"), with = FALSE], numeric0))
  d <- d[FILTER == "PASS" & af <= 0.001 & Variant_Type == "SNP" & Chromosome %in% std_chr &
           Reference_Allele %in% c("A", "C", "G", "T") & Tumor_Seq_Allele2 %in% c("A", "C", "G", "T")]
  g <- GRanges(d$Chromosome, IRanges(d$Start_Position, width = 1L),
       REF = DNAStringSet(d$Reference_Allele), ALT = DNAStringSetList(as.list(d$Tumor_Seq_Allele2)))
  g <- subsetByOverlaps(g, territory)
  seqlevels(g) <- std_chr; genome(g) <- "hg38"
  g
}
target_gr <- GRangesList(lapply(inputs$source_maf, one_target_sbs))
names(target_gr) <- inputs$cell_line
target_mat <- mut_matrix(target_gr, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
readr::write_csv(as.data.frame(target_mat) %>% tibble::rownames_to_column("context"),
                file.path(OUT, "wes_sbs_context_target_restricted.csv"))
ctx_tri <- paste0(substr(rownames(target_mat), 1, 1), substr(rownames(target_mat), 3, 3),
                  substr(rownames(target_mat), 7, 7))
mult <- (tri$target_fraction / tri$genome_fraction)[match(ctx_tri, tri$trinucleotide)]
stopifnot(length(mult) == 96L, all(is.finite(mult)), all(mult > 0))
# Transform P(context | signature, genome) to target opportunities, then renormalise
# each signature to sum1. Restrict observed variants to the same target territory.
normalized_reference <- cosmic * mult
normalized_reference <- sweep(normalized_reference, 2, colSums(normalized_reference), "/")
readr::write_csv(as.data.frame(normalized_reference) %>%
  mutate(context = rownames(target_mat), .before = 1), file.path(OUT, "wes_cosmic_target_normalized.csv"))
fit_target <- function(reference, mode) {
  f <- fit_to_signatures_strict(target_mat, reference, max_delta = MAX_DELTA, method = "backwards")$fit_res
  rel <- sweep(f$contribution, 2, colSums(f$contribution), "/")
  data.frame(cell_line = colnames(target_mat), mode = mode,
    n_snv_target = as.integer(colSums(target_mat)),
    n_snv_unrestricted = as.integer(n_used[colnames(target_mat)]),
    mmr_relative_exposure = colSums(rel[intersect(MMR_D_SIGS, rownames(rel)), , drop = FALSE]),
    reconstruction_cosine = vapply(seq_len(ncol(target_mat)), function(i)
      cos_sim(target_mat[, i], f$reconstructed[, i]), numeric(1)),
    target_bed_sha256 = unique(bed_sha), target_union_bp = sum(width(territory)),
    caveat = "exploratory point refit; derived CNVkit target footprint, not sample-callable mask; tumour-only candidates")
}
target_sensitivity <- bind_rows(
  fit_target(cosmic, "target-restricted variants; genome reference"),
  fit_target(normalized_reference, "target-restricted variants; target-opportunity-adjusted reference"),
  fit_target(normalized_reference[, intersect(RESTRICTED_SIGS, colnames(normalized_reference)), drop = FALSE],
             "target-restricted variants; target-adjusted restricted reference"))
readr::write_csv(target_sensitivity, file.path(OUT, "wes_signature_target_sensitivity.csv"))
message("Target-opportunity sensitivity: ", sum(width(territory)), "bp; ", nrow(inputs), "models; TOV21G:")
print(target_sensitivity[target_sensitivity$cell_line == "TOV21G", 1:6], row.names = FALSE)
