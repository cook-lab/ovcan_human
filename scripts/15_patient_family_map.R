#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 15_patient_family_map.R
#
# Canonical line -> patient -> subline-family map for the OvCAN panel.
#
# WHY: several "cell lines" are sublines/isolates from the SAME patient (the
# Mes-Masson OV/TOV isolates and -R re-derivations). Peer review flagged that
# genomic frequencies computed per-LINE pseudoreplicate these patients — e.g.
# the 3133 family's TP53 p.Q192* + CDK12 frameshift is one event counted 4x.
# This produces the single source of truth so every downstream tabulation can
# collapse to INDEPENDENT PATIENTS and every per-line figure can carry a family
# annotation track.
#
# INPUT : metadata/samples.csv        (pure metadata -> NO pipeline dependency)
# OUTPUT: metadata/line_family_map.csv
#
# Depends only on samples.csv, so it may be run at any time (it is consumed by
# the re-tabulations in 07/08 and by 16-20; it does not read pipeline outputs).
# Family membership is from CHUM provenance (samples.csv `notes`) and was
# confirmed at the variant level during review (shared somatic TP53 alleles:
# 1369->G244C, 2295->I195T, 3133->Q192*).
#
# *** EXECUTION ORDER — NUMERIC ORDER IS NOT DEPENDENCY ORDER *** [review revision]
#   This script writes metadata/line_family_map.csv, which 07, 08, 16 and 20
#   CONSUME. A numeric top-to-bottom run therefore only worked because the map was
#   already committed; from a clean checkout 07 would fail. Rather than renumber
#   this script (the manuscript, reports/ and docs/ all cite it as "script 15"),
#   00_setup.R now provides ensure_family_map(), and 07/08/16 call it before
#   reading the map: it re-runs this file if the map is absent. Run order is
#   therefore either `15 -> 07 -> 08 -> 16 -> ... -> 20` (preferred, explicit) or
#   plain numeric order (safe, because the guard regenerates the map on demand).
# ---------------------------------------------------------------------------

PROJ <- Sys.getenv("OVCAN_PROJ", unset = getwd())

samples <- read.csv(file.path(PROJ, "metadata", "samples.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)

# Resource = generated, analysis-included lines only (external/Carey removed).
gen <- samples[samples$provenance == "generated" & samples$analysis_include == "Y", ]

# --- Explicit same-patient families ----------------------------------------
families <- list(
  "1369" = c("OV1369-R2", "TOV1369"),
  "2295" = c("OV2295", "OV2295-R2", "TOV2295-R"),
  "3133" = c("OV3133-R", "OV3133-R2", "TOV3133D", "TOV3133G"),
  "3291" = c("OV3291", "TOV3291G"),
  "3121" = c("TOV3121D", "TOV3121EP")
)

# sanity: every declared family member must exist in the generated set
declared <- unlist(families, use.names = FALSE)
missing  <- setdiff(declared, gen$cell_line)
if (length(missing)) stop("Family members not found in generated lines: ",
                          paste(missing, collapse = ", "))

patient_id <- setNames(gen$cell_line, gen$cell_line)   # default: singleton = own name
family     <- setNames(rep(NA_character_, nrow(gen)), gen$cell_line)
for (fam in names(families)) {
  patient_id[families[[fam]]] <- fam
  family[families[[fam]]]     <- fam
}

gen$patient_id <- patient_id[gen$cell_line]
gen$family     <- family[gen$cell_line]

# --- Per-assay availability -------------------------------------------------
gen$has_rna     <- gen$rna_seq    == "Y"
gen$has_prot    <- gen$proteomics == "Y"
gen$has_wes_cnv <- gen$wes_cnv    == "Y"
gen$has_wes_maf <- gen$wes_mut    == "Y" & !grepl("no MAF", gen$wes_mut_note)

# --- Family size + a deterministic representative ---------------------------
# Representative = most-complete-omics line per patient (tri-omic preferred),
# ties broken alphabetically. Used where a single line per patient is wanted
# (e.g. expression PCA); frequency collapses instead use ">=1 line per patient".
gen$omics_score <- gen$has_rna + gen$has_prot + gen$has_wes_maf + 0.1 * gen$has_wes_cnv

tab <- table(gen$patient_id)
gen$n_lines_in_family   <- as.integer(tab[gen$patient_id])
gen$is_multiline_family <- gen$n_lines_in_family > 1L

ord <- order(gen$patient_id, -gen$omics_score, gen$cell_line)
g   <- gen[ord, ]
g$patient_representative <- !duplicated(g$patient_id)

out <- g[, c("cell_line", "subtype", "source_site", "patient_id", "family",
             "n_lines_in_family", "is_multiline_family", "patient_representative",
             "has_rna", "has_prot", "has_wes_cnv", "has_wes_maf")]
out <- out[order(out$patient_id, out$cell_line), ]

outfile <- file.path(PROJ, "metadata", "line_family_map.csv")
write.csv(out, outfile, row.names = FALSE)

# --- Console summary --------------------------------------------------------
n_lines    <- nrow(out)
n_patients <- length(unique(out$patient_id))
hgs <- out[out$subtype == "HGS", ]
hgs_maf         <- hgs[hgs$has_wes_maf, ]
hgs_maf_patients <- length(unique(hgs_maf$patient_id))

cat(sprintf("Wrote %s\n", outfile))
cat(sprintf("Generated lines: %d  ->  independent patients: %d  (collapsed %d lines)\n",
            n_lines, n_patients, n_lines - n_patients))
cat(sprintf("Multi-line families: %d  (%s)\n",
            sum(out$is_multiline_family[out$patient_representative]),
            paste(names(families), collapse = ", ")))
cat(sprintf("HGS with WES-MAF: %d lines -> %d patients\n",
            nrow(hgs_maf), hgs_maf_patients))
cat("\nMulti-line families:\n")
print(out[out$is_multiline_family, c("cell_line","patient_id","subtype",
                                    "patient_representative","has_wes_maf")],
      row.names = FALSE)
