# Project Spec: OvCAN Human Ovarian Cancer Cell-Line Multi-Omic Resource

**Last updated:** 2026-07-23
**Status:** Planning → re-analysis
**Target venue:** *Scientific Data* (Data Descriptor)
**Detailed synthesis:** `reports/00_synthesis_and_recommendations.md`

## Overview
A uniform multi-omic characterization of a curated panel of human ovarian cancer cell-line models (the OvCAN / Ovarian Cancer Canada collection), generated across three matched layers — bulk RNA-seq (kallisto), TMT proteomics (Morin lab, GSC), and whole-exome sequencing (Sarek: Mutect2 + CNVkit). Sources: CHUM (Mes-Masson/Provencher; TOV/OV lines), BC Cancer/OVCARE (Huntsman; VOA lines), and OHRI (Vanderhyden; BIN67). Repackaged from an MSc thesis (J. Sobh) into a data descriptor.

## Framing (decided)
A **validation-focused Data Descriptor**, not a discovery paper. The paper argues that the data recapitulate known subtype biology (so they can be trusted and reused) and provides the QC, authentication, and Data-Records machinery for reuse. One featured **usage example** is retained: a subtype-resolved **ADC-target expression atlas** (framed as model-selection utility, not clinical target discovery). Discovery biology (heterogeneity, SWI/SNF convergence, single-line stories) is deferred to a possible companion paper.

## Questions the resource answers (descriptor-appropriate)
1. **Quality:** Are the RNA-seq, proteomics, and WES data technically sound (depth/alignment, protein/peptide identifications + batch reproducibility, coverage/on-target)?
2. **Identity:** Are the lines authentic and correctly annotated (STR; histotype from CNV+mutation+expression), and does the data independently recover published reclassifications (COV434→SCCOHT; TOV112D→dedifferentiated)?
3. **Biological fidelity:** Do the models recapitulate canonical subtype biology — subtype separation (RNA + protein), known markers/GO terms, canonical genomics (TP53 ubiquity in HGSC, HGSC CNV landscape), and expected RNA–protein concordance (benchmarked to CPTAC ovarian, 0.38/0.45)?
4. **Reuse:** How can others use it to select and interrogate models (by subtype, genotype, ADC-target expression, genomic HRD, within-subtype heterogeneity)?

## Resource definition
- **Generated data** (source of truth: `metadata/samples.csv`, `provenance=generated`): **42 cell-line models** with in-house multi-omic data — RNA-seq **31** · TMT proteomics **31** · WES (CNV + SNV) **23** · tri-omic backbone **13** (all CHUM). Provenance confirmed empirically via FASTQ paths (in-house run NS.1676.003 vs external Carey). Subtypes: HGSC, clear cell, mucinous, endometrioid, MMMT/carcinosarcoma, SCCOHT; plus LGS (TOV81D only — annotated as the sole LGSC in the collection, WES-only).
- **Framing nuance:** we generated the *omic data*; the cell lines themselves are pre-existing/published — the descriptor documents our data on them and cites original derivation (Provencher/Létourneau/Sauriol for CHUM TOV/OV; Huntsman for VOA; Gamwell for BIN67; van den Berg-Bakker for COV434).
- **Excluded as external:** the LGSC/Carey models (data pulled from published Carey/OVCARE work, not generated here) — cited (e.g., Shrestha 2021), not deposited.
- **Deposition:** full *generated* dataset with QC-based include/exclude flags — GEO/SRA (RNA-seq, WES), PRIDE/MassIVE via ProteomeXchange (proteomics).

## Novelty (calibrated)
Breadth across subtypes + scale (~40 lines) + uniform matched WES/RNA-seq/TMT + curated/corrected annotations; anchored most strongly by **carcinosarcoma** (no modern omics-characterized 2D lines) and **SCCOHT** (~3–4 lines worldwide). Not "first multi-omic ovarian cell-line data" (cf. Shrestha 2021, Coscia 2016) and not "CCLE lacks rare subtypes."

## Deliverables
1. *Scientific Data* Data Descriptor (5 main figures; see ANALYSIS_PLAN.md).
2. Deposited raw + processed data with a documented sample sheet.
3. A reproducible, version-pinned analysis pipeline (`scripts/`, `renv`).

## Authorship (from latest draft; to confirm)
J. Sobh · G. Negri · E. Carmona · B. C. Vanderhyden · G. Morin · D. Huntsman · A.-M. Mes-Masson · D. P. Cook (corresponding).

## Out of scope
LGSC/Carey external data; hypothesis-driven discovery biology; drug-response/functional assays beyond any existing metadata.
