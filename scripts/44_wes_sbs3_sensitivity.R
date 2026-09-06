#!/usr/bin/env Rscript
# Exploratory SBS3 identifiability/sensitivity, separate from canonical outputs.
# Does not estimate HRD probability, clinical status, somatic burden or HR function.
# No setup.R: it refreshes canonical provenance; this script writes only DEST.
ROOT <- normalizePath(Sys.getenv("OVCAN_PROJ", getwd()), mustWork = TRUE)
DEST <- file.path(ROOT, "reports/molecular_extension_2026-09-06/signatures")
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)
for (p in c("MutationalPatterns", "data.table", "digest", "jsonlite", "GenomicRanges",
            "Biostrings", "BSgenome.Hsapiens.UCSC.hg38")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Required installed package: ", p)
}
# Use namespaces: attaching MutationalPatterns/NMF can trigger unrelated CPU
# autodetection failure under the desktop sandbox. No package code is patched.
N_BOOT <- as.integer(Sys.getenv("OVCAN_SBS3_BOOTS", "200"))
WORKERS <- as.integer(Sys.getenv("OVCAN_SBS3_WORKERS", "3"))
stopifnot(N_BOOT >= 100L, WORKERS >= 1L)
SEED <- 440609L
MAX_DELTA <- 0.004
COSMIC_SOURCE <- "COSMIC_v3.2"
GENOME <- "GRCh38"
std_chr <- paste0("chr", c(1:22, "X", "Y"))
restricted_names <- c("SBS1", "SBS5", "SBS40", "SBS6", "SBS14", "SBS15", "SBS20",
  "SBS21", "SBS26", "SBS44", "SBS10a", "SBS10b", "SBS10c", "SBS10d", "SBS28",
  "SBS2", "SBS13", "SBS3", "SBS8", "SBS17a", "SBS17b", "SBS18")
sha <- function(p) digest::digest(file = p, algo = "sha256")
read_csv <- function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
write_csv <- function(x, name) utils::write.csv(x, file.path(DEST, name), row.names = FALSE, na = "")
resolve <- function(p) {
  if (startsWith(p, "/")) return(p)
  if (nzchar(Sys.getenv("OVCAN_DATA")) && startsWith(p, "judy_archive/data/"))
    return(file.path(Sys.getenv("OVCAN_DATA"), sub("^judy_archive/data/", "", p)))
  file.path(ROOT, p)
}
matrix_csv <- function(p) {
  x <- read_csv(p)
  m <- as.matrix(x[, -1L, drop = FALSE]); rownames(m) <- x[[1L]]; m
}
inputs <- read_csv(file.path(ROOT, "output/wes_input_manifest.csv"))
family <- read_csv(file.path(ROOT, "metadata/line_family_map.csv"))
family <- family[match(inputs$cell_line, family$cell_line), ]
stopifnot(nrow(inputs) == 23L, !anyNA(family$patient_id), length(unique(family$patient_id)) == 16L)
target_mat_old <- matrix_csv(file.path(ROOT, "output/wes_sbs_context_target_restricted.csv"))
cosmic <- MutationalPatterns::get_known_signatures("snv", source = COSMIC_SOURCE, genome = GENOME)
stopifnot(nrow(cosmic) == 96L, ncol(cosmic) == 60L)
# This package's bundled reference has no row names; its SBS96 order is the
# mut_matrix order used by script22. The independent transform equality below
# verifies this mapping against the already deposited normalized reference.
rownames(cosmic) <- rownames(target_mat_old)
target_ref_old <- matrix_csv(file.path(ROOT, "output/wes_cosmic_target_normalized.csv"))
tri <- read_csv(file.path(ROOT, "output/wes_target_trinucleotide_opportunities.csv"))
bed_path <- file.path(ROOT, "output/wes_cnvkit_target_intervals.bed")
bed_hash <- sha(bed_path)
stopifnot(all(tri$target_bed_sha256 == bed_hash), nrow(tri) == 32L,
          bed_hash == "2177970ffdbc933e488067f812e3ae4760310bd9a40fadf4e2fa5389480d955d")
bed <- data.table::fread(bed_path, select = 1:3, col.names = c("chr", "start", "end"))
bed <- bed[chr %in% std_chr]
territory <- GenomicRanges::reduce(GenomicRanges::GRanges(bed$chr, IRanges::IRanges(bed$start + 1L, bed$end)))
stopifnot(sum(BiocGenerics::width(territory)) == 63465385L)
trinuc <- paste0(substr(rownames(cosmic), 1, 1), substr(rownames(cosmic), 3, 3), substr(rownames(cosmic), 7, 7))
mult <- (tri$target_fraction / tri$genome_fraction)[match(trinuc, tri$trinucleotide)]
target_ref <- sweep(cosmic * mult, 2L, colSums(cosmic * mult), "/")
stopifnot(max(abs(target_ref - target_ref_old[rownames(cosmic), colnames(cosmic)])) < 1e-12)
refs <- list(genome_full = cosmic, genome_restricted = cosmic[, restricted_names],
             target_full = target_ref, target_restricted = target_ref[, restricted_names])
dict <- do.call(rbind, lapply(names(refs), function(nm) data.frame(reference = nm,
  signature = colnames(refs[[nm]]), sbs3_reference_cosine = as.numeric(
    MutationalPatterns::cos_sim_matrix(refs[[nm]][, "SBS3", drop = FALSE], refs[[nm]])))))
write_csv(dict, "reference_dictionaries.csv")

# Build exactly the canonical target-restricted baseline from MAF PASS and the
# four population fields; VCF AD/AF adds read-support screens, not somatic proof.
af_cols <- c("gnomADe_AF", "AF", "AA_AF", "EA_AF")
numeric0 <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  stopifnot(all(is.na(x) | !is.na(z)))
  z[is.na(z)] <- 0; z
}
ref_genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
substrates <- c("baseline", "rare", "rare_read_supported")
matrices <- lapply(substrates, function(x) matrix(0L, 96L, 23L,
                    dimnames = list(rownames(cosmic), inputs$cell_line)))
names(matrices) <- substrates
burden <- vector("list", nrow(inputs))
for (i in seq_len(nrow(inputs))) {
  id <- inputs$cell_line[i]
  maf_path <- resolve(inputs$source_maf[i]); vcf_path <- resolve(inputs$source_vcf[i])
  stopifnot(sha(maf_path) == inputs$maf_sha256[i], sha(vcf_path) == inputs$vcf_sha256[i])
  d <- data.table::fread(maf_path, skip = "Hugo_Symbol", quote = "", na.strings = c("", "NA", "."),
    select = c("FILTER", af_cols, "Variant_Type", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
  d <- d[FILTER == "PASS" & Variant_Type == "SNP" & Chromosome %in% std_chr &
         Reference_Allele %in% c("A", "C", "G", "T") & Tumor_Seq_Allele2 %in% c("A", "C", "G", "T")]
  d$pop_af_max <- do.call(pmax, lapply(d[, af_cols, with = FALSE], numeric0))
  d <- d[pop_af_max <= 0.001]
  g <- GenomicRanges::GRanges(d$Chromosome, IRanges::IRanges(d$Start_Position, width = 1L))
  keep <- GenomicRanges::countOverlaps(g, territory) > 0
  d <- d[keep]; g <- g[keep]
  stopifnot(!anyDuplicated(paste(d$Chromosome, d$Start_Position, d$Reference_Allele, d$Tumor_Seq_Allele2)))
  con <- if (grepl("\\.gz$", vcf_path)) gzfile(vcf_path, "rt") else file(vcf_path, "rt")
  lines <- readLines(con); close(con)
  z <- data.table::tstrsplit(lines[!startsWith(lines, "#")], "\t", fixed = TRUE,
                           keep = c(1L, 2L, 4L, 5L, 9L, 10L))
  key <- paste(z[[1]], z[[2]], z[[3]], z[[4]])
  m <- match(paste(d$Chromosome, d$Start_Position, d$Reference_Allele, d$Tumor_Seq_Allele2), key)
  stopifnot(!anyNA(m))
  geno <- lapply(seq_along(m), function(j) setNames(strsplit(z[[6]][m[j]], ":", fixed = TRUE)[[1]],
                                                       strsplit(z[[5]][m[j]], ":", fixed = TRUE)[[1]]))
  ad <- lapply(geno, function(x) suppressWarnings(as.integer(strsplit(x[["AD"]], ",", fixed = TRUE)[[1]])))
  stopifnot(all(lengths(ad) == 2L), !anyNA(unlist(ad)))
  depth <- vapply(ad, sum, 0L); alt_count <- vapply(ad, `[`, 0L, 2L)
  caller_af <- vapply(geno, function(x) as.numeric(x[["AF"]]), 0)
  stopifnot(all(is.finite(caller_af)), all(caller_af >= 0 & caller_af <= 1))
  triseq <- BSgenome::getSeq(ref_genome, GenomicRanges::resize(g, width = 3L, fix = "center"))
  stopifnot(all(substr(as.character(triseq), 2L, 2L) == d$Reference_Allele))
  purine <- d$Reference_Allele %in% c("A", "G")
  alt <- d$Tumor_Seq_Allele2
  triseq[purine] <- Biostrings::reverseComplement(triseq[purine])
  alt[purine] <- chartr("ACGT", "TGCA", alt[purine])
  ctx <- paste0(substr(as.character(triseq), 1L, 1L), "[", substr(as.character(triseq), 2L, 2L),
                ">", alt, "]", substr(as.character(triseq), 3L, 3L))
  stopifnot(all(ctx %in% rownames(cosmic)))
  masks <- list(baseline = rep(TRUE, nrow(d)), rare = d$pop_af_max <= 1e-5,
    rare_read_supported = d$pop_af_max <= 1e-5 & depth >= 20 & alt_count >= 5 & caller_af >= 0.05)
  for (nm in substrates) matrices[[nm]][, id] <- tabulate(match(ctx[masks[[nm]]], rownames(cosmic)), nbins = 96L)
  burden[[i]] <- data.frame(cell_line = id, patient_id = family$patient_id[i], subtype = family$subtype[i],
    substrate = substrates, n_snv = vapply(masks, sum, 0L), baseline_n_snv = nrow(d),
    baseline_read_supported = sum(depth >= 20 & alt_count >= 5 & caller_af >= .05),
    baseline_missing_population_annotation = sum(rowSums(is.na(d[, af_cols, with = FALSE])) == 4L),
    population_af_max = c(.001, 1e-5, 1e-5), min_AD_sum = c(NA, NA, 20),
    min_AD_alt = c(NA, NA, 5), min_caller_AF = c(NA, NA, .05),
    n_nonzero_contexts = vapply(masks, function(v) length(unique(ctx[v])), 0L))
  message("Substrate verified: ", id, " ", paste(vapply(masks, sum, 0L), collapse = " / "))
}
stopifnot(identical(unname(matrices$baseline), unname(target_mat_old[rownames(cosmic), inputs$cell_line])))
burden <- do.call(rbind, burden)
burden$fraction_retained <- burden$n_snv / burden$baseline_n_snv
burden$count_band <- ifelse(burden$n_snv < 200, "<200 SNVs", ifelse(burden$n_snv < 500, "200-499 SNVs", ">=500 SNVs"))
write_csv(burden, "substrate_burden.csv")
for (nm in substrates) write_csv(data.frame(context = rownames(cosmic), matrices[[nm]], check.names = FALSE), paste0("sbs96_", nm, ".csv"))

# Each job uses a sample/substrate seed that is identical across dictionaries and
# opportunity models, giving matched multinomial resamples for sensitivity.
jobs <- expand.grid(cell_line = inputs$cell_line, substrate = substrates,
                    reference = names(refs), stringsAsFactors = FALSE)
jobs$bootstrap <- jobs$substrate != "rare" # Intermediate AF-only screen is a point-fit diagnostic.
cache_dir <- file.path(DEST, "cache")
dir.create(cache_dir, showWarnings = FALSE)
script_hash <- sha(file.path(ROOT, "scripts/44_wes_sbs3_sensitivity.R"))
run_one <- function(i) {
  job <- jobs[i, ]; id <- job$cell_line; substrate <- job$substrate; ref_name <- job$reference
  mat <- matrices[[substrate]][, id, drop = FALSE]; sigs <- refs[[ref_name]]
  seed <- SEED + strtoi(substr(digest::digest(paste(id, substrate), serialize = FALSE), 1L, 7L), 16L)
  key <- digest::digest(list(mat, sigs, seed, N_BOOT, MAX_DELTA, job$bootstrap,
    as.character(utils::packageVersion("MutationalPatterns")), "strict-backwards-v1"), algo = "sha256")
  cache <- file.path(cache_dir, paste0(id, "_", substrate, "_", ref_name, "_", key, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  stopifnot(sum(mat) > 0L)
  fit <- MutationalPatterns::fit_to_signatures_strict(mat, sigs, max_delta = MAX_DELTA,
                                                    method = "backwards")$fit_res
  contrib <- setNames(rep(0, ncol(sigs)), colnames(sigs)); contrib[rownames(fit$contribution)] <- fit$contribution[, 1L]
  frac <- contrib / sum(contrib)
  nnls <- MutationalPatterns::fit_to_signatures(mat, sigs)
  nnls_no3 <- MutationalPatterns::fit_to_signatures(mat, sigs[, colnames(sigs) != "SBS3", drop = FALSE])
  sse <- sum((mat - nnls$reconstructed)^2); sse_no3 <- sum((mat - nnls_no3$reconstructed)^2)
  stopifnot(sse_no3 >= sse - 1e-6)
  summary <- data.frame(cell_line = id, substrate = substrate, reference = ref_name,
    dictionary = ifelse(grepl("restricted", ref_name), "restricted", "full"),
    opportunity_model = ifelse(startsWith(ref_name, "target"), "target", "genome"),
    n_snv = sum(mat), n_reference_signatures = ncol(sigs), n_selected = sum(contrib > 0),
    sbs3_count = contrib["SBS3"], total_fitted = sum(contrib), sbs3_fraction = frac["SBS3"],
    reconstruction_cosine = as.numeric(MutationalPatterns::cos_sim(as.numeric(mat), as.numeric(fit$reconstructed))),
    nnls_without_sbs3_delta_sse_per_n_squared = (sse_no3 - sse) / sum(mat)^2,
    nnls_with_sbs3_cosine = as.numeric(MutationalPatterns::cos_sim(as.numeric(mat), as.numeric(nnls$reconstructed))),
    nnls_without_sbs3_cosine = as.numeric(MutationalPatterns::cos_sim(as.numeric(mat), as.numeric(nnls_no3$reconstructed))),
    n_boots = if (job$bootstrap) N_BOOT else 0L, seed = seed, cache_sha256 = key,
    boot_fraction_lo95 = NA_real_, boot_fraction_hi95 = NA_real_, boot_fraction_median = NA_real_,
    boot_selected_fraction = NA_real_, boot_count_lo95 = NA_real_, boot_count_hi95 = NA_real_)
  boot <- NULL; competitors <- NULL
  if (job$bootstrap) {
    set.seed(seed)
    bt <- MutationalPatterns::fit_to_signatures_bootstrapped(mat, sigs, n_boots = N_BOOT,
      max_delta = MAX_DELTA, method = "strict_backwards", verbose = FALSE)
    aligned <- matrix(0, N_BOOT, ncol(sigs), dimnames = list(NULL, colnames(sigs)))
    aligned[, colnames(bt)] <- bt
    stopifnot(all(rowSums(aligned) > 0))
    bf <- aligned[, "SBS3"] / rowSums(aligned)
    boot <- data.frame(cell_line = id, substrate = substrate, reference = ref_name,
      replicate = seq_len(N_BOOT), sbs3_count = aligned[, "SBS3"], total_fitted = rowSums(aligned), sbs3_fraction = bf)
    summary$boot_fraction_lo95 <- quantile(bf, .025)
    summary$boot_fraction_hi95 <- quantile(bf, .975)
    summary$boot_fraction_median <- median(bf)
    summary$boot_selected_fraction <- mean(bf > 0)
    summary$boot_count_lo95 <- quantile(aligned[, "SBS3"], .025)
    summary$boot_count_hi95 <- quantile(aligned[, "SBS3"], .975)
    bc <- suppressWarnings(vapply(colnames(aligned), function(s) stats::cor(aligned[, "SBS3"], aligned[, s]), 0))
    competitors <- data.frame(cell_line = id, substrate = substrate, reference = ref_name,
      signature = names(bc), bootstrap_count_correlation_with_sbs3 = bc)
  }
  result <- list(summary = summary, boot = boot, competitors = competitors,
    exposure = data.frame(cell_line = id, substrate = substrate, reference = ref_name,
                         signature = names(contrib), count = contrib, fraction = frac))
  saveRDS(result, cache)
  message("Fit complete: ", id, " / ", substrate, " / ", ref_name)
  result
}
message("Running ", nrow(jobs), " fits; ", sum(jobs$bootstrap), " with ", N_BOOT, " bootstrap resamples; workers=", WORKERS)
results <- parallel::mclapply(seq_len(nrow(jobs)), run_one, mc.cores = WORKERS,
                             mc.set.seed = FALSE, mc.preschedule = FALSE)
failed <- which(vapply(results, inherits, FALSE, "try-error"))
if (length(failed)) stop("Signature job failed: ", as.character(results[[failed[1L]]]))
summaries <- do.call(rbind, lapply(results, `[[`, "summary"))
summaries <- merge(summaries, family[, c("cell_line", "patient_id", "subtype", "patient_representative")], by = "cell_line", sort = FALSE)
summaries$interpretation <- "conditional SBS3 fit; not an HRD classification; tumour-only, no sample-callable mask"
write_csv(summaries, "sbs3_sensitivity_summary.csv")
boots <- do.call(rbind, lapply(results, `[[`, "boot"))
write_csv(boots, "sbs3_bootstrap_replicates.csv")
write_csv(do.call(rbind, lapply(results, `[[`, "competitors")), "sbs3_bootstrap_competitors.csv")
write_csv(do.call(rbind, lapply(results, `[[`, "exposure")), "all_signature_point_exposures.csv")
contrasts <- list(); k <- 0L
for (id in inputs$cell_line) for (sub in c("baseline", "rare_read_supported")) {
  b <- boots[boots$cell_line == id & boots$substrate == sub, ]
  for (pair in list(c("genome_full", "target_full"), c("genome_restricted", "target_restricted"),
                    c("target_full", "target_restricted"))) {
    a <- b[b$reference == pair[1L], ]; z <- b[b$reference == pair[2L], ]
    stopifnot(identical(a$replicate, z$replicate))
    delta <- z$sbs3_fraction - a$sbs3_fraction
    k <- k + 1L
    contrasts[[k]] <- data.frame(cell_line = id, substrate = sub, from_reference = pair[1L], to_reference = pair[2L],
      paired_boot_delta_median = median(delta), paired_boot_delta_lo95 = quantile(delta, .025),
      paired_boot_delta_hi95 = quantile(delta, .975), n_boots = N_BOOT)
  }
}
write_csv(do.call(rbind, contrasts), "paired_reference_sensitivity.csv")
patient <- aggregate(sbs3_fraction ~ patient_id + substrate + reference, summaries, mean)
names(patient)[4L] <- "mean_model_sbs3_fraction"
patient$n_models <- vapply(seq_len(nrow(patient)), function(i) sum(summaries$patient_id == patient$patient_id[i] &
  summaries$substrate == patient$substrate[i] & summaries$reference == patient$reference[i]), 0L)
write_csv(patient, "patient_descriptive_sensitivity.csv")
source_paths <- c("output/wes_input_manifest.csv", "output/wes_sbs_context_target_restricted.csv",
  "output/wes_cosmic_target_normalized.csv", "output/wes_target_trinucleotide_opportunities.csv",
  "output/wes_cnvkit_target_intervals.bed", "metadata/line_family_map.csv", "scripts/44_wes_sbs3_sensitivity.R")
write_csv(data.frame(path = source_paths, sha256 = vapply(file.path(ROOT, source_paths), sha, "")), "source_manifest.csv")
jsonlite::write_json(list(script_sha256 = script_hash, models = nrow(inputs), patients = length(unique(family$patient_id)),
  n_boot = N_BOOT, workers = WORKERS, seed_base = SEED, strict_max_delta = MAX_DELTA,
  seed_policy = "same sample/substrate seed across references; paired multinomial resamples",
  cosmic_source = COSMIC_SOURCE, genome = GENOME, target_union_bp = sum(BiocGenerics::width(territory)),
  reference_signatures = vapply(refs, ncol, 0L), baseline_exactly_reproduces_canonical_target_matrix = TRUE,
  source_maf_and_vcf_hashes_verified = TRUE, samples = inputs,
  versions = setNames(lapply(c("MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg38", "data.table", "digest"),
    function(p) as.character(utils::packageVersion(p))), c("MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg38", "data.table", "digest")),
  R = as.character(getRversion()), completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  caveat = "Bootstrap is conditional sampling variability, not clinical confidence or uncertainty in germline contamination, capture callability, reference choice, or calling."),
  file.path(DEST, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
writeLines(capture.output(sessionInfo()), file.path(DEST, "session_info.txt"))
message("Completed standalone SBS3 sensitivity; all outputs in ", DEST)
