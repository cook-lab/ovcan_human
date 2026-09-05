#!/usr/bin/env Rscript
# Limited faidx-compatible adapter for vcf2maf reference checks, using the
# already-installed hg38 TwoBit reference. No alignment/indexing functionality.
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) >= 3L, a[1] == "faidx", file.exists(a[2]))
suppressPackageStartupMessages({library(rtracklayer); library(GenomicRanges)})
regions <- a[-c(1, 2)]
chr <- sub(":.*$", "", regions)
lo <- as.integer(sub(".*:([0-9]+)-[0-9]+$", "\\1", regions))
hi <- as.integer(sub(".*-([0-9]+)$", "\\1", regions))
stopifnot(all(!is.na(lo)), all(!is.na(hi)), all(lo > 0), all(hi >= lo))
seqs <- import(TwoBitFile(a[2]), which = GRanges(chr, IRanges(lo, hi)))
stopifnot(length(seqs) == length(regions))
for (i in seq_along(regions)) cat(">", regions[i], "\n", as.character(seqs[[i]]), "\n", sep = "")
