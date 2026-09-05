# Exploration: hypermutation / MSI / mutational-context (wes-signatures workstream)
# Goal: confirm TOV21G load, indel fraction, Ts/Tv, MMR/POLE hits; size raw PASS set.
source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })

maf <- readr::read_csv(file.path(OUT, "wes_mutations_filtered.csv"), show_col_types = FALSE)
fam <- readr::read_csv(file.path(META, "line_family_map.csv"), show_col_types = FALSE)

cat("=== Filtered MAF overview ===\n")
cat("rows:", nrow(maf), " lines:", dplyr::n_distinct(maf$cell_line), "\n")
cat("\nVariant_Type levels:\n"); print(table(maf$Variant_Type))
cat("\nVariant_Classification levels:\n"); print(sort(table(maf$Variant_Classification), decreasing = TRUE))

# --- Per-line load: SNV / indel / MNV, indel fraction, Ts/Tv ---
is_ts <- function(ref, alt) {
  p <- paste0(pmin(ref, alt), pmax(ref, alt))
  p %in% c("AG", "CT")   # transitions A<->G, C<->T
}
snv <- maf %>% filter(Variant_Type == "SNP",
                      Reference_Allele %in% c("A","C","G","T"),
                      Tumor_Seq_Allele2 %in% c("A","C","G","T")) %>%
  mutate(ts = is_ts(Reference_Allele, Tumor_Seq_Allele2))

load_tbl <- maf %>%
  group_by(cell_line, subtype) %>%
  summarise(
    n_total   = n(),
    n_snv     = sum(Variant_Type == "SNP"),
    n_ins     = sum(Variant_Type == "INS"),
    n_del     = sum(Variant_Type == "DEL"),
    n_indel   = sum(Variant_Type %in% c("INS","DEL")),
    n_mnv     = sum(Variant_Type %in% c("DNP","TNP","ONP")),
    .groups = "drop") %>%
  mutate(indel_frac = n_indel / (n_snv + n_indel))

tstv <- snv %>% group_by(cell_line) %>%
  summarise(ts = sum(ts), tv = sum(!ts), tstv = ts/tv, .groups="drop")
load_tbl <- load_tbl %>% left_join(tstv, by = "cell_line") %>%
  arrange(desc(n_total))

med <- median(load_tbl$n_total); madv <- mad(load_tbl$n_total)
load_tbl <- load_tbl %>%
  mutate(fold_over_median = n_total / med,
         robust_z = (n_total - med) / madv)

cat("\n=== Per-line load (ranked) ===\n")
print(as.data.frame(load_tbl %>% mutate(across(where(is.numeric), ~round(.,3)))), row.names = FALSE)

cat(sprintf("\nMedian total coding-nonsyn load = %.0f ; MAD = %.1f\n", med, madv))
top2 <- load_tbl$n_total[1:2]
cat(sprintf("TOP line: %s = %d ; 2nd = %s = %d ; fold(1st/2nd) = %.2f ; fold(1st/median) = %.2f\n",
            load_tbl$cell_line[1], top2[1], load_tbl$cell_line[2], top2[2],
            top2[1]/top2[2], top2[1]/med))

# --- MMR / proofreading gene scan ---
MMR   <- c("MLH1","MSH2","MSH6","PMS2","PMS1","MSH3","EPCAM")
PROOF <- c("POLE","POLD1")
cat("\n=== MMR / proofreading gene hits across panel (coding-nonsyn) ===\n")
mmr_hits <- maf %>% filter(Hugo_Symbol %in% c(MMR, PROOF)) %>%
  select(cell_line, subtype, Hugo_Symbol, Variant_Classification, HGVSp_Short,
         Consequence, IMPACT, CLIN_SIG, pop_af_max, vaf, germline_like_vaf)
print(as.data.frame(mmr_hits), row.names = FALSE)

cat("\n=== TOV21G: all MMR/proofreading hits ===\n")
print(as.data.frame(mmr_hits %>% filter(cell_line == "TOV21G")), row.names = FALSE)

# POLE exonuclease domain check (~codons 268-471)
pole <- maf %>% filter(Hugo_Symbol == "POLE") %>%
  mutate(codon = as.integer(str_extract(HGVSp_Short, "(?<=[A-Z])[0-9]+")))
cat("\n=== POLE variants (codon parsed; exonuclease ~268-471) ===\n")
if (nrow(pole)) print(as.data.frame(pole %>% select(cell_line, HGVSp_Short, codon,
      Variant_Classification, IMPACT, CLIN_SIG, vaf)), row.names = FALSE) else cat("none\n")

# --- Size the raw PASS SNV set (for possible SBS 96-context) ---
MUT_DIR <- file.path(DATA, "wes - old", "mutect2")
canon_key <- function(x) toupper(gsub("[-_]", "", sub("_[Pp][0-9]+$", "", x)))
maf_files <- list.files(MUT_DIR, pattern = "\\.maf$", recursive = TRUE, full.names = TRUE)
# TOV21G raw
f21 <- maf_files[grepl("TOV21G", maf_files, ignore.case = TRUE)]
cat("\n=== Raw MAF PASS sizing (TOV21G) ===\nfile:", f21, "\n")
r <- fread(f21, sep = "\t", skip = "Hugo_Symbol", quote = "",
           na.strings = c("","NA","."), showProgress = FALSE)
af_cols <- intersect(c("gnomADe_AF","AF","AA_AF","EA_AF"), names(r))
num_or0 <- function(v){ x <- suppressWarnings(as.numeric(v)); ifelse(is.na(x),0,x) }
r[, pop_af_max := do.call(pmax, lapply(.SD, num_or0)), .SDcols = af_cols]
r[, retained := FILTER == "PASS" & pop_af_max <= 0.001]
cat("raw rows:", nrow(r), " PASS:", sum(r$FILTER=="PASS"),
    " retained(PASS&popAF):", sum(r$retained), "\n")
cat("retained SNV (single-base):",
    sum(r$retained & r$Variant_Type=="SNP" &
        r$Reference_Allele %in% c("A","C","G","T") &
        r$Tumor_Seq_Allele2 %in% c("A","C","G","T")), "\n")
cat("retained by Variant_Type:\n"); print(table(r$Variant_Type[r$retained]))
cat("\nNCBI_Build (raw MAF):", unique(as.character(r$NCBI_Build)), "\n")
