<!--
FIRST FULL DRAFT — Scientific Data Data Descriptor — OvCAN human ovarian cancer cell-line multi-omic resource
Prepared: 2026-07-24 by descriptor-draft (Cook Lab / Claude), from reports/01–03, PROJECT_SPEC, and reports/lit_review/.
Structure follows Scientific Data Data Descriptor policy: Title · Abstract · Background & Summary · Methods · Data Records · Technical Validation · Usage Notes · Code Availability. NO Results/Discussion (by editorial policy).
Figure/table numbering follows reports/03_figures_and_webapp_plan.md §B (the FINAL numbering): Figs 1–6 + S1–S8; Table 1 = Data Records; Table S1 = supplement_per_line.csv.
Citations are given inline as author-year for ease of PI verification against reports/lit_review/; convert to Scientific Data numbered style at submission (see References note at end).
Deposition specifics are marked [PLACEHOLDER: …]; genuine gaps are marked [NOTE: …]; numbers to re-confirm are marked [CHECK].

STATUS 2026-07-27: SUPERSEDED — HISTORICAL. The live manuscript is docs/manuscript/v4/OvCAN_data_descriptor_v4.md (this file is the v1 draft; v2, v3 and v4 followed). Do not quote statistics from this file. Claims corrected by the methodological audit of 2026-07-27 (reports/06_methodological_validity_and_submission_readiness_review.md) are annotated inline as [CORRECTED 2026-07-27]; the authoritative versions are in reports/01_multiomic_characterization_results.md and the v4 draft.
-->

# Title

**A uniform multi-omic resource of 42 human ovarian cancer cell-line models spanning common and rare histotypes**

<!-- Recommended title (leads with scope + the rare-subtype differentiator, states the number). Alternates for PI consideration:
 (2) Matched RNA-seq, quantitative proteomics, and whole-exome sequencing of the OvCAN ovarian cancer cell-line panel
 (3) Multi-omic characterization of ovarian cancer cell-line models across seven histotypes, including SCCOHT and carcinosarcoma -->

---

# Abstract

<!-- Target ≤170 words, unreferenced (per Scientific Data policy). Current count ≈ 168. -->

Epithelial ovarian cancer comprises several histotypes that differ in genetics, cell of origin, and treatment response, yet existing cell-line omics resources are weighted toward high-grade serous and other common subtypes. We describe OvCAN, a uniformly processed multi-omic resource for 42 human ovarian cancer cell-line models derived from 34 patients and spanning seven histotypes: high-grade serous, clear cell, mucinous, endometrioid, low-grade serous, small-cell carcinoma of the ovary hypercalcemic type (SCCOHT), and carcinosarcoma. The resource comprises bulk RNA-seq (31 lines), tandem-mass-tag quantitative proteomics (31 lines), and whole-exome sequencing (23 lines), with 13 lines profiled by all three assays. The lines originate from three centres and carry curated, literature-corrected histotype annotations. Raw and processed data are deposited in public repositories. Technical validation shows that the data are of high quality, recover canonical subtype biology and genetics, including TP53 mutation in all 11 high-grade serous patients, and molecularly corroborate line identity against DepMap and Cellosaurus. A consolidated per-line table and worked reuse examples support model selection.

---

# Background & Summary

Epithelial ovarian cancer is not a single disease. Its major histotypes, high-grade serous carcinoma (HGSC), clear cell, endometrioid, mucinous, and low-grade serous carcinoma, differ in cell of origin, driver genetics, and response to therapy (Cancer Genome Atlas Research Network 2011; Wiegand et al. 2010; Cheasley et al. 2019; Hollis et al. 2020). Rarer entities extend this diversity further: SCCOHT is defined by SMARCA4 (BRG1) inactivation (Ramos et al. 2014; Witkowski et al. 2014; Jelinic et al. 2014; Karnezis et al. 2016), and ovarian carcinosarcoma remains among the most aggressive and least modelled presentations (Wong et al. 2025). Preclinical work depends on cell-line models that faithfully represent these distinct histotypes, but the field's most-used "HGSC" lines are a recognised weakness: genomic benchmarking against tumours placed SKOV3 and A2780 among the least tumour-like HGSC models (Domcke et al. 2013; Anglesio et al. 2013), and cell-line misidentification is widespread and systematically catalogued (Korch et al. 2012; Huang et al. 2017; ICLAC 2026).

Large multi-omic cell-line resources have set both the technical standard and the biological expectations for reuse. The Cancer Cell Line Encyclopedia established genomic and transcriptomic profiling at scale (Barretina et al. 2012; Ghandi et al. 2019), later extended by tandem-mass-tag (TMT) proteomics on 375 lines (Nusinow et al. 2020) and the ProCan proteomic map of 949 lines (Gonçalves et al. 2022). These pan-cancer resources are weighted toward common cancers and toward HGSC within ovarian cancer, and rare ovarian histotypes are thinly represented. Ovarian-specific efforts have characterised smaller panels: label-free proteomics of 26 lines (Coscia et al. 2016), the OCCP transcriptomic panel of 39 lines (Beaufort et al. 2014), and, most closely related to the present work, matched WES, RNA-seq, and mass-spectrometry proteomics of 14 low-grade serous lines from the group that supplies several of the lines described here (Shrestha et al. 2021). No pan-histotype ovarian cell-line panel with uniform matched WES, RNA-seq, and TMT proteomics has been described, and no ovarian-cell-line multi-omic *Scientific Data* Data Descriptor was found (literature review, this project).

The OvCAN collection (Ovarian Cancer Canada) assembles curated ovarian cancer cell-line models from three centres: the Centre hospitalier de l'Université de Montréal (CHUM; the Mes-Masson and Provencher TOV/OV lines), BC Cancer / OVCARE (the Huntsman VOA lines), and the Ottawa Hospital Research Institute (OHRI; the Vanderhyden BIN67 line). Here we describe uniform multi-omic data generated on 42 of these models, drawn from 34 patients and spanning the seven histotypes above (Fig. 1). We generated the omic data; the cell lines themselves are pre-existing and previously published, and we cite their original derivation throughout (Provencher et al. 2000; Ouellet et al. 2008; Létourneau et al. 2012; Fleury et al. 2015, 2016; Sauriol et al. 2020; Gamwell et al. 2013; van den Berg-Bakker et al. 1993). The resource's contribution is best stated with calibration: its value lies in breadth across histotypes, scale (~40 lines), uniform matched processing of three omic layers, and curated annotations that incorporate two published reclassifications, rather than in being the first ovarian cell-line omics dataset. Its strongest gap-filling coverage is for histotypes that are scarce elsewhere, namely SCCOHT (of which only a handful of lines exist worldwide) and carcinosarcoma (which lacks modern omics-characterised 2D models; Wong et al. 2025).

Consistent with the aims of a Data Descriptor, this paper documents how the resource was assembled and processed and validates its technical quality for reuse; it does not test biological hypotheses or draw new biological conclusions. Technical Validation proceeds in three steps that together support reliability: sequencing and proteome data quality (Fig. 2); recovery of canonical, histotype-specific biology and genetics, which is evidence that the data are trustworthy rather than a discovery (Figs 3, 4); and molecular corroboration of line identity, including independent recovery of two published reclassifications, against DepMap and Cellosaurus (Figs 4, 5). Usage Notes then provide a consolidated per-line table (Table S1) and worked reuse examples for model selection, including a subtype-resolved antibody-drug-conjugate (ADC) target atlas (Fig. 6) and a candidate rare-subtype model (Fig. 5).

---

# Methods

> This section documents how the data and the validation analyses were produced, in enough detail to reproduce. Parameters are taken from the numbered analysis scripts (`scripts/00`–`scripts/20`).

## Cell lines, provenance, and annotation

The resource comprises 42 human ovarian cancer cell-line models from three centres: CHUM (Mes-Masson/Provencher TOV and OV lines), BC Cancer / OVCARE (Huntsman VOA lines), and OHRI (Vanderhyden BIN67). Original derivation and prior characterisation are cited per line: the founding CHUM lines including TOV-21G, TOV-81D, OV-90, and TOV-112D (Provencher et al. 2000); additional serous lines (Ouellet et al. 2008); matched primary/recurrent HGSC lines that are the source of the OV2295 and OV3133 lineages (Létourneau et al. 2012); HGSC lines including the germline BRCA1/BRCA2-mutant OV4453/OV4485 (Fleury et al. 2015, 2016); ten newer lines including one mucinous and one clear cell line (Sauriol et al. 2020); the SCCOHT line BIN-67 (Gamwell et al. 2013); and COV434, originally described as a granulosa-tumour line (van den Berg-Bakker et al. 1993). We generated the omic data described here; we did not derive the lines.

Provenance of the generated data was verified empirically from raw-file paths and run identifiers: all deposited RNA-seq derives from a single in-house sequencing run (`NS.1676.003`), and all proteomics derives from a single five-set TMT design, distinguishing these data from externally sourced low-grade serous (Carey/OVCARE) data that were excluded from the resource and are cited rather than deposited (Shrestha et al. 2021).

Several lines are sublines or primary/recurrent isolates from a shared patient, so the 42 lines represent 34 independent patients. Patient-family structure is recorded in `metadata/line_family_map.csv` (script `15`): five multi-line families, 1369 (2 lines), 2295 (3), 3121 (2), 3133 (4), and 3291 (2), were collapsed to one representative per patient for all genomic frequency statistics, with the representative chosen by omic coverage. Family membership was confirmed at the variant level (for example, all four family-3133 lines share the TP53 p.Q192\* allele). Genomic frequencies throughout are reported per patient, with `n_patients` stated alongside `n_lines`; per-line values are retained in the supplement for transparency.

Histotype annotation is recorded in `metadata/samples.csv` (the single source of truth for sample identity, assay availability, passage, and QC flags). Annotations are curated and, where the literature supports it, corrected (see Technical Validation §4). Across the 34 patients, the histotype composition is HGSC (16 patients / 24 lines), clear cell (8 / 8), mucinous (3 / 3), endometrioid (2 / 2), carcinosarcoma (MMMT; 2 / 2), SCCOHT (2 / 2), and low-grade serous (1 / 1, WES-only).

## RNA sequencing and processing

Libraries were prepared and sequenced in a single run (`NS.1676.003`) [CHECK: read configuration and library-prep kit to be confirmed for Methods]. Reads were pseudoaligned with kallisto v0.46 against a transcriptome index built from Ensembl release 105. Transcript-to-gene mapping used a release-105 `tx2gene` table that is pinned and cached in the repository; the biomaRt query that generated it is invoked lazily only on a cache miss, which makes the mapping reproducible and removes a dependence on a live database call. Transcript abundances were summarised to the gene level with tximport (`ignoreTxVersion = TRUE`), yielding 39,568 genes. A DESeq2 object was constructed from gene-level counts; genes expressed at least 10 counts in at least two lines were retained (22,544 genes), and a variance-stabilising transformation (VST) was applied for distance-based and dimensionality-reduction analyses. Thirty-one lines passed RNA-seq QC and constitute the RNA analysis set.

## Proteomics and processing

Quantitative proteomics was generated by the Morin laboratory (Canada's Michael Smith Genome Sciences Centre) using isobaric tandem-mass-tag (TMT) labelling across five 11-plex sets. Two features of the design are load-bearing for reuse. First, inter-set normalisation is a common-reference hub: channel 1 of every set carries a Pooled Internal Standard (PIS), a common reference to which each channel is ratio-normalised; scripts `05`/`06` apply no further inter-set centring. Second, channel 10 carries a separate technical-replicate "bridge" sample that forms a daisy-chain of four consecutive links across the five sets (set 1↔2↔3↔4↔5), used for reproducibility QC rather than for normalisation. [NOTE: raw mass-spectrometry files are required for ProteomeXchange/PRIDE deposition and are being obtained from the Morin laboratory; processed abundance matrices are in hand.]

The search results comprise 8,430 proteins across 146,830 peptides. Missingness is structural rather than random: it is per-set and all-or-nothing, because a protein quantified in one TMT set may be absent from another. Of 8,430 proteins, 6,856 (81.3%) are quantified across all five sets, and 1,573 (18.7%) are absent from at least one whole set (median 13 of 31 lines lost for an affected protein; per-protein set coverage is tabulated in `output/prot_block_missingness.csv`). Because set composition is correlated with histotype, present/absent status for the affected proteins is set-conditional and not biologically interpretable. A presence filter requiring quantification in at least 50% of lines retains 7,734 proteins for analysis. Thirty-one lines constitute the proteomics analysis set.

## Whole-exome sequencing and processing

Whole-exome sequencing was processed with the nf-core/Sarek pipeline, calling somatic single-nucleotide variants and indels with Mutect2 in tumour-only mode (no matched normal) and copy number with CNVkit against a pooled diploid-normal reference [CHECK: capture kit identity to be stated in Methods]. Twenty-three lines have CNVkit copy-number profiles; 22 have Mutect2 mutation calls (MAF).

Variant filtering corrects a root-cause error in the archived analysis, which ignored the Mutect2 `FILTER` field so that MAF files contained all candidate calls rather than passing calls (for example, the OV2295 MAF contained 25,914 rows, of which 493 are `PASS`). We retained `FILTER == PASS` calls and removed variants present in gnomAD, 1000 Genomes, or ESP population-allele-frequency databases. Because the lines are near-100% pure and unmatched by a normal, residual germline variation persists after filtering (~206 coding variants per line). Accordingly, drivers are reported per patient with an explicit somatic-confidence tier rather than as a mutation-burden metric. Each of 50 canonical driver calls was assigned to one of three tiers (`output/wes_driver_tiers.csv`): Tier 1 (n=27), a hotspot or clear loss-of-function in a tumour-suppressor gene; Tier 2 (n=11), plausible somatic; Tier 3 (n=12), cannot exclude germline (excluded from headline frequencies). The archived `germline_like_vaf` flag was retired as a somatic filter, because in near-pure lines it flags bona fide somatic hotspots (for example TP53 R175H and Q192\*) whose variant-allele fraction approaches 1 through loss of heterozygosity rather than germline origin.

Copy number was analysed autosome-restricted with a per-sample median-centred baseline. The chromosome-X signal was excluded as a pooled-normal sex-composition artifact (per-line chrX median log2 copy-ratio swings from −1.23 to +0.83 with the sign inverting between lines); chromosome Y was already excluded. Fraction of the genome altered (FGA) is defined on autosomes. Three caveats are documented and carried into interpretation: the five public normal exomes in the CNVkit reference (PRJNA339046) are unmatched and their capture-kit concordance is unconfirmed (Talevich et al. 2016); median-centring is flagged for the two genomes with autosome FGA > 0.7 (TOV3121D, TOV2929D); and the archived MAF files are aligned to GRCh38 despite an erroneous `NCBI_Build = GRCh37` header (verified from the panel-of-normals, contig lengths, and hotspot coordinates), to be corrected on deposition.

## Computational analysis supporting validation

Analyses were run under system R v4.5.2 with a fixed seed (1234); `scripts/00_setup.R` records paths, package versions, and the seed.

*Dimensionality reduction and confounder modelling.* Subtype separation was assessed by principal-component analysis on VST expression (RNA) and on the protein abundance matrix, with silhouette widths computed per annotated histotype. Because all HGSC lines derive from one centre, histotype is partially confounded with source site; rather than univariate R², we fit a joint model `PC ~ subtype + site` and decomposed the variance by commonality analysis into unique-subtype, unique-site, and shared components (script `17`). The one histotype present at both centres, clear cell (n=7: 2 CHUM + 5 BC Cancer), provides the sole clean cross-site control. Genome-wide variance was decomposed with a per-gene linear mixed model fit by REML (an equivalent of `variancePartition`, which is unavailable in this R build; see Code Availability).

*Marker and pathway recovery.* Canonical histotype markers were scored by a top-2-of-6 rule and quantified by effect size (Cohen's d and marker-vs-rest AUC) across 28 patient-level RNA representatives (`output/rna_marker_effectsizes.csv`, script `04`). Gene-set and GO recovery used ranked differential-expression signatures; formal differential expression was restricted to HGSC (adequately sampled) and treated descriptively for histotypes with n≤3.

*Cross-assay concordance.* RNA-protein concordance was computed as per-gene and per-line Spearman correlation on shared genes, and stratified by protein dynamic-range terciles (scripts `12`, `19`).

*Legacy subtype classification.* ConsensusOV HGSC subtype calls were generated for the 15 HGSC RNA lines and are surfaced as an orthogonal legacy annotation with an explicit microenvironment caveat (script `10`; see Technical Validation §2).

*Mutational signatures.* Single-base-substitution spectra were extracted with MutationalPatterns on the GRCh38 BSgenome and compared to COSMIC SBS reference signatures by cosine similarity (script `16`).

*External benchmarking.* Cross-platform identity and driver concordance were assessed against DepMap Public 24Q4 (Figshare deposition 27993248) for the five overlapping lines, and STR/record status was queried from the Cellosaurus API (script `18`). External data are cached in `output/external/`.

---

# Data Records

All raw and processed data are being deposited in community repositories, and accessions will be finalised before publication. RNA-seq and WES raw reads and processed matrices are deposited under a GEO SuperSeries linked to SRA; proteomics raw and search-result files are deposited via ProteomeXchange/PRIDE; and processed multi-omic matrices, metadata, and derived tables are deposited to a figshare/Zenodo record. The sample sheet `metadata/samples.csv` is the master join key across all records; every processed matrix is keyed to the `cell_line` identifier.

The data records are summarised in Table 1.

**Table 1. Data Records.**

| Record | Assay | Repository | Accession | Format | n | File(s) |
|---|---|---|---|---|---|---|
| Raw RNA-seq reads | Bulk RNA-seq | SRA (via GEO SuperSeries) | [PLACEHOLDER: SRP/GEO] | FASTQ | 31 | per-line paired-end FASTQ |
| Processed RNA-seq | Bulk RNA-seq | GEO | [PLACEHOLDER: GSE] | CSV/TSV | 31 | gene × line matrices (TPM, counts, VST) |
| Raw mass spectrometry | TMT proteomics | ProteomeXchange / PRIDE | [PLACEHOLDER: PXD] | RAW + mzML + search | 31 | 5 TMT 11-plex sets |
| Processed proteomics | TMT proteomics | figshare / PRIDE | [PLACEHOLDER: DOI/PXD] | CSV | 31 | protein abundance matrix; per-protein set coverage |
| WES mutation calls | WES | GEO/SRA or figshare | [PLACEHOLDER] | VCF / MAF | 22 | per-line filtered MAF |
| WES copy number | WES | figshare / Zenodo | [PLACEHOLDER: DOI] | CSV / SEG | 23 | CNVkit segment tables |
| Raw WES reads | WES | SRA | [PLACEHOLDER — pending, see note] | FASTQ | 23 | per-line reads |
| Metadata & derived tables | All | figshare / Zenodo | [PLACEHOLDER: DOI] | CSV | 42 | `samples.csv`, `line_family_map.csv`, `supplement_per_line.csv` (Table S1), driver tiers, marker effect sizes, ConsensusOV calls, block-missingness |

File organisation, naming keyed to `cell_line`, and per-matrix column and unit definitions are provided in the repository README and in `metadata/samples.csv`.

[NOTE: three deposition-completeness items gate what Data Records can claim and require PI resolution.
 (1) WES raw reads: the analysis notes indicate tumour BAMs are not archived; confirm that FASTQs exist for SRA deposit, otherwise the WES record is processed-only (VCF/MAF/CNV) with a stated limitation.
 (2) Raw MS availability from the Morin laboratory for ProteomeXchange.
 (3) Deposition hygiene: correct the MAF `NCBI_Build` header (GRCh37→GRCh38) and strip the bundled third-party normal exomes from the CNVkit outputs before deposit.]

---

# Technical Validation

## 1. Sequencing and proteome data quality (Fig. 2)

*RNA-seq.* Pseudoalignment rates are high and uniform across the 31 lines (median 91.1%, range 85.8–93.1%), with a median library size of 56.5M reads and a median of ~20,119 genes detected per line. A modest source-site difference in pseudoalignment rate (Mes-Masson 92.2% vs Huntsman 88.1%) does not translate into a difference in the number of genes detected, so it does not affect downstream detection. Pipeline identity is further supported externally for the five lines catalogued in DepMap (§4).

*Proteomics.* Bridge technical-replicate reproducibility is high: per-link Pearson correlation is 0.991–0.994 across the four daisy-chain links, with a median coefficient of variation of 5.3%. This level of reproducibility is sufficient for technical validation independent of the raw search pipeline (which is nonetheless required for deposition). The structural missingness described in Methods is quantified per protein (`output/prot_block_missingness.csv`; Fig. S4): 81.3% of proteins are complete across all five sets, and the remaining 18.7% are set-conditional. Reporting present/absent status for the affected proteins as biology would be an artifact, and the resource is documented accordingly.

*WES.* Applying the `FILTER == PASS` and population-allele-frequency filters restores a credible variant landscape from MAF files that had contained all candidate calls (the OV2295 example, 25,914 rows → 493 PASS, is shown as a filtering waterfall in Fig. 4). [NOTE: on-target and coverage metrics can be reported only if the raw reads or alignment metrics are recoverable; if not, they should be stated as not computable from the archived outputs (Fig. S8 is contingent on this).]

## 2. The data recover canonical subtype biology (Fig. 3)

Recovering known, histotype-specific biology from independent assays is evidence that the data are trustworthy, and this is the intended reading of this section rather than a biological finding.

*Subtype separation is driven by biology, not source site.* Histotypes separate on the leading RNA principal components (PC1 20.7%, PC2 10.4%) and, more weakly, in protein. The confounding concern (all HGSC lines derive from one centre) is addressed by the joint model: on `PC ~ subtype + site` with commonality decomposition, source site adds ≤0.2% of PC1 variance beyond histotype (unique-subtype 42%, unique-site 0.2%, shared 31%), and the same holds on PC2–PC3. The clean cross-site control agrees: within clear cell (2 CHUM + 5 BC Cancer), source lab explains only 4–6% of the leading components, so clear-cell lines from the two centres co-cluster by biology. Genome-wide variance decomposition is concordant: subtype ≥ site in RNA (median 5.9% vs 3.5%) and subtype ≫ site ≈ set in protein (8.6% vs 0.0% vs 0.9%), with line/patient identity the largest structured term, as expected for distinct cell lines. Silhouette widths are reported with sample sizes and are modest for the well-sampled histotypes (HGSC 0.16, clear cell 0.12, mucinous 0.15); the high values for MMMT (0.74) and SCCOHT (0.82) rest on n=2 and are not over-interpreted, and the two endometrioid lines do not co-cluster (see §4).

*Marker and pathway recovery.* Of 22 canonical histotype markers, 16 land in the expected histotype by a top-2-of-6 rule (`output/rna_marker_effectsizes.csv`; Fig. 3, Fig. S7). Effect sizes quantify the rule: median oriented marker-vs-rest AUC is 0.69 and median |Cohen's d| is 0.72, with 8 of 22 markers reaching AUC ≥ 0.80 (KRT20 d = 3.1, CDX2 d = 2.8, SPP1 d = 1.5, HNF1B d = 1.4). The six markers that miss the top-2 cutoff are small-effect or shared-lineage rather than failures (for example PAX8 is pan-Müllerian, d = 0.53; MUC2 is silenced in 2D culture, d = −0.29; MECOM d ≈ 0). GO recovery aligns with known biology: HGSC DNA-repair, clear-cell oxidative/glutathione, and SCCOHT cell-cycle signatures are recovered; mucinous glycan biology is suggestive; MMMT EMT is not recovered at n=2. One caveat is on record: the strongest HGSC GO axis (OXPHOS/ribosome biogenesis) is confounded with site and is not correctable, so lineage-specific marker recovery, which is not a plausible batch artifact, is the load-bearing evidence here.

*RNA-protein concordance sits where this assay pair lands.* [CORRECTED 2026-07-27: the heading previously read "at the expected, compression-limited ceiling"; the narrower protein spread is consistent with TMT ratio compression but these data cannot isolate a compression factor or establish a causal concordance ceiling, and the sentence beginning "The moderate value reflects an intrinsic TMT ratio-compression ceiling" below is retired for the same reason.] On 30 lines and 8,212 shared genes, per-gene median Spearman correlation is 0.40 (per-line median 0.41, IQR 0.36–0.44). This sits on the same-disease CPTAC ovarian proteogenomic benchmark (mean 0.38 / median 0.45; Zhang et al. 2016) and within the cell-line range (ProCan 0.42, CCLE 0.48, integrated landscape 0.58; Gonçalves et al. 2022; Nusinow et al. 2020; Jarnuczak et al. 2021), against a measurement-reproducibility ceiling of ~0.72 (Upadhya & Ryan 2022). The moderate value reflects an intrinsic TMT ratio-compression ceiling rather than weak data: protein cross-line spread is ~0.30× the RNA spread, and concordance rises monotonically with protein dynamic range (protein-IQR terciles: median Spearman 0.30 / 0.40 / 0.53), which is why resources sharing the TMT limitation land at the same value.

*Legacy HGSC subtype calls (reported with caveat).* ConsensusOV calls for the 15 HGSC RNA lines (differentiated 6, mesenchymal 4, immunoreactive 3, proliferative 2; `output/consensusov_calls.csv`, Table S5, Fig. S6) are provided as an orthogonal legacy annotation, not a validated axis. The TCGA/ConsensusOV mesenchymal and immunoreactive subtypes are substantially microenvironment-driven, being defined by stromal and immune compartments (Zhang et al. 2019; Schwede et al. 2020; Olbrecht et al. 2021; Chen et al. 2018). Pure tumour-cell cultures lack these compartments, so the fact that 7 of 15 lines are called mesenchymal or immunoreactive is itself evidence that the labels are not ground truth for cell lines. They are reported as provenance only; the tumour-cell-intrinsic within-HGSC strata (Usage Notes; Fig. S8) are the appropriate resolution.

## 3. Genomic fidelity (Fig. 4)

The WES data recover canonical ovarian cancer genetics, which validates the sequencing and the filtering strategy.

*Driver mutations.* TP53 is mutated in all 11 HGSC patients (17/17 HGSC lines), the expected near-universal positive control for HGSC (Ahmed et al. 2010), and it survives patient collapse. Patient-level tabulation also corrects pseudoreplication that inflates per-line frequencies: the apparent CDK12 enrichment (6/17 lines, 35%) collapses to 2–3 of 11 patients (18–27%), matching the ~3% TCGA-HGSC rate, because one family's single frameshift was counted four times. Reporting is restricted to robust events (TP53; KRAS G12/G13 hotspots in clear cell and mucinous; CTNNB1/PIK3CA/PTEN in clear cell; SMARCA4 truncation in TOV112D), each carrying a somatic-confidence tier (27 Tier 1, 11 Tier 2, 12 Tier 3; Table S2). Defensible somatic BRCA1/2 is zero: both candidate BRCA2 calls are Tier 3 and excluded (one is incoherent with an otherwise quiet, HR-proficient genome; the other is a noise-floor N-terminal missense at VAF 0.038), and tumour-only WES with population-AF filtering cannot classify a rare germline BRCA1/2 variant as somatic (Halperin et al. 2017; Little et al. 2021).

*Copy number.* CNV recovers the textbook HGSC landscape, and unlike point mutations, arm-level frequencies are near-invariant to patient collapse because they are trunk events shared across sublines (per patient: 3q gain 82%, 20q gain 91%, 17p loss 82%, 8q gain 73%, 13q loss 64%, 19q gain 55%). This asymmetry is itself informative: pseudoreplication distorts mutation frequencies badly but CNV frequencies barely. Autosome-restricted FGA orders the histotypes HGSC 0.62 (n=11) > clear cell 0.37 (n=2) > mucinous 0.32 (n=1) > endometrioid 0.23 (n=1) > low-grade serous 0.02 (n=1); the non-HGSC values rest on n=1 or on two dissimilar clear-cell lines and are reported without smoothing. Excluding chromosome X (the pooled-normal artifact) barely changes high-FGA HGSC genomes (median 0.61→0.62) but corrects quiet genomes (TOV81D 0.073→0.021).

## 4. Cell-line identity and authentication (Figs 4, 5)

*External identity validation.* Five of the 42 lines are catalogued in DepMap (OV90, TOV21G, TOV112D, BIN67, COV434; verified by DepMap `Model.csv` and Cellosaurus RRID; no BC Cancer VOA line exists publicly). Each self-matches its DepMap namesake at Spearman 0.74–0.88, ranking first of 67 DepMap ovarian lines and reciprocal-best in both directions; because the comparison is cross-platform (kallisto vs RSEM), specificity is the signal (Fig. 4). Driver calls corroborate identity: DepMap independently recovers TOV21G hypermutation (568 damaging variants vs 7–18 in other lines), TOV112D TP53-R175H and SMARCA4 damaging, OV90 SMAD4, and all five TOV21G clear-cell drivers. Two non-canonical tumour-only variants (TOV112D KRAS-A59T, OV90 BRAF in-frame indel) are not corroborated in DepMap and are flagged.

*Independent recovery of two published reclassifications.* The resource's own multi-omic data recover two literature reclassifications, reinforced by DepMap (Fig. 5; Karnezis et al. 2021). COV434, historically a granulosa-tumour line, and BIN67 are consistent with SCCOHT: both show SMARCA4/SMARCA2 loss, DepMap independently shows SMARCA4 damaging in both, and Cellosaurus flags COV434 as reclassified to SCCOHT. BIN67 is the instructive case: SMARCA4 mRNA is retained but protein is second-lowest (post-transcriptional loss), so SMARCA4 IHC/protein is confirmatory and RNA-only authentication would miss it. TOV112D, historically endometrioid, is consistent with dedifferentiated carcinoma: SMARCA4 truncation (p.L639X), low protein, silenced SMARCA2 mRNA, and TP53 R175H, an SWI/SNF-null pattern rather than typical endometrioid.

*STR and mycoplasma statement.* No STR profiling or mycoplasma testing was performed in-house. For 30 of the 42 models a documented STR profile deposited by the originating laboratories exists in Cellosaurus, and each accession is listed in the supplement (Table S1, Table S4). The remaining 12 models (the 11 BC Cancer/OVCARE VOA lines and TOV3121D) have no public STR record, and this set contains the identity-doubt lines below, so in-house STR and IHC are specifically requested for them. For the five DepMap lines, identity is further supported by the expression concordance and driver cross-check above.

*Histotype consistency and flagged discordances.* Most lines are histotype-concordant across expression and genomics (26 expression-consistent, 20 genomics-consistent; `output/auth_perline_table.csv`). Documented flags for reuse: OV90 is best read as an HGSC-family carcinoma whose serous identity is not confirmed by expression (PAX8/WT1/SOX17/CK7 all ≈ 0, the lowest serous profile among HGSC lines); VOA5436 (annotated MMMT) expresses a strong clear-cell program; VOA4841 (clear cell) has the lowest SMARCA4 of all lines; VOA4395 becomes the sole endometrioid line once TOV112D is reclassified. Mucinous authenticity is the highest-risk histotype (Stordal et al. 2024; Meagher et al. 2025): TOV2414 is a well-authenticated ovarian mucinous line (KRT7+/PAX8+/MUC5AC+/SATB2−, KRAS G12A; Sauriol et al. 2020; Cellosaurus CVCL_A1SR), whereas VOA8762 and VOA8771 show an intestinal/GI-leaning profile (CK7-low, PAX8-low, CDX2-high), have no external provenance and no public STR, and their ovarian origin is not supported by expression (Fig. 5). Genuine ovarian mucinous coverage may therefore be n=1 (TOV2414) pending STR and IHC.

[NOTE: genomic HRD is not computable from the archived data. CNVkit produced total copy number only (no allele-specific/BAF; the tumour BAMs are not archived), so scarHRD cannot run and no HRD score was computed (Sztupinszki et al. 2018). A genuine score would require recovering the WES BAMs (Sequenza → scarHRD); this is a PI decision. Expression-based "HRD" proxies are not genomic scars and are not reported as HRD (Peng et al. 2014; Takamatsu et al. 2024).]

---

# Usage Notes

## Model selection by histotype and genotype

The consolidated per-line supplement (`output/supplement_per_line.csv`; Table S1) is the recommended entry point for reuse. It provides one row per line with assay availability, QC metrics, patient-family membership, tiered driver calls, autosome FGA, MSI/MMR status, histotype-consistency and authentication flags, ConsensusOV call, and Cellosaurus/STR status, so that a reuser can select models by histotype, genotype, data availability, or authentication confidence in a single table.

## ADC-target expression atlas (Fig. 6)

Antibody-drug conjugates are a fast-moving therapeutic class in ovarian cancer, and patient selection is expression-dependent (Matulonis et al. 2023; Moore et al. 2023; Meric-Bernstam et al. 2024), which motivates a subtype-resolved target-expression atlas in models (`output/adc_expression.csv`; Fig. 6). Known associations recover in both RNA and protein: mesothelin is highest in HGSC, and HER2/ERBB2 is higher in clear cell and mucinous. FOLR1 is strongly bimodal within HGSC (RNA range 0.06–9.5), and the highest FRα expressers panel-wide are HGSC lines (top: TOV3133G), which mirrors the clinical FRα heterogeneity that mandates IHC selection for mirvetuximab and is itself a model-selection feature. Two reuse caveats apply. First, TMT ratio compression narrows protein target ranges 3–5× relative to RNA (for example MSLN protein range 1.8 vs RNA 9.3), so shortlists should be led with RNA and confirmed with protein rank-consistency; cross-assay-consistent examples include HER2→TOV3392D, NaPi2b→VOA12539, and CDH6→OV3331. Second, target expression is necessary but not sufficient for ADC benefit: several agents with high target expression had negative pivotal trials (NaPi2b, mesothelin, DPEP3; Hamilton et al. 2020), so the atlas is a preclinical model-selection aid, not clinical target validation.

## Rare-subtype and candidate models

The resource is enriched for histotypes that are scarce in existing resources. It provides two SCCOHT models (COV434, BIN67) and two carcinosarcoma models (VOA5217, VOA5436) with matched omics. It also nominates a candidate rare-biology model for reuse: TOV21G (clear cell) is a clear hypermutation outlier (1,416 coding candidate variants, 6.9× the panel median and 3.4× the next line), enriched for indels, with a GRCh38 SBS-96 spectrum matching mismatch-repair-deficient/MSI signatures (SBS6 cosine 0.88, plus SBS44/SBS15/SBS20) and no POLE signal, alongside biallelic ARID1A truncation, PIK3CA H1047Y, and KRAS G13C (Fig. 5). No MMR-enzyme coding mutation is present, consistent with the epigenetic (MLH1-methylation) mechanism typical of MSI-high clear-cell carcinoma, which WES cannot detect. Because the WES is tumour-only, this is a converging-evidence call: TOV21G is best described as a candidate MMR-deficient/MSI-high ovarian clear-cell model for immunotherapy or MMR biology, pending MMR IHC/MSI-PCR confirmation. The hypermutation is independently recovered in DepMap (§4), and the other clear-cell line (TOV3392D) is not hypermutated, so the signal is line-specific.

## Within-HGSC heterogeneity (example only; Fig. S8)

As a model-selection example, the 15 HGSC lines resolve into three descriptive pathway-activity strata (inflammatory/NF-κB-EMT; low-signaling; hypoxic-glycolytic), corroborated orthogonally by PROGENy, with same-patient families mostly co-clustering (`output/hgs_heterogeneity.csv`; Fig. S8). This is provided only to illustrate within-histotype model selection (n=15); no survival or discovery claims are made.

## Reuse caveats (essential)

- **Mutations are tumour-only.** Use the per-patient, tiered driver calls (Table S2), not a mutation-burden metric; treat Tier 3 calls as germline-indeterminate.
- **Genomic HRD is not computable** from the archived data; do not infer HRD from these files (see Technical Validation §4 note).
- **Proteomics has TMT ratio compression and structural per-set missingness;** lead quantitative comparisons with RNA and treat set-conditional present/absent calls as non-biological (Table S6, Fig. S4).
- **Copy number uses unmatched public normals of unconfirmed capture-kit concordance;** confirm concordance before quantitative CNV reuse, and note the two high-FGA genomes flagged for median-centring.
- **Obtain authenticated stocks** with STR and mycoplasma clearance for any critical reuse, particularly the VOA lines and the flagged identity-doubt lines. [CORRECTED 2026-07-27: no in-house STR profiling or mycoplasma testing was performed in this project, and a Cellosaurus record documents a third-party reference profile rather than a match to these stocks — so the panel must never be described as authenticated.]

---

# Code Availability

The analysis is implemented as numbered R scripts (`scripts/00`–`scripts/20`) that run top-to-bottom under system R v4.5.2, with `scripts/00_setup.R` recording paths, package versions, and a fixed random seed (1234). Transcript-to-gene mapping is pinned to Ensembl release 105 via a cached `tx2gene` table (biomaRt is invoked lazily only on a cache miss). External datasets (DepMap Public 24Q4, Figshare deposition 27993248; Cellosaurus API) are cached under `output/external/` for reproducibility. One documented substitution: `variancePartition` is unavailable in this R build (a `lme4`/`reformulas` incompatibility), so an equivalent per-gene `lme4` REML variance decomposition is used and documented in `scripts/17`. [NOTE: the `renv` lockfile will be finalised last, once the environment settles, and referenced here; the public repository host (GitHub with a Zenodo DOI) is to be confirmed by the PI.]

---

<!-- ============================================================ -->
# References

*Draft citation format note: references are given inline as author-year for ease of verification against `reports/lit_review/`. Convert to Scientific Data's numbered/superscript style at submission. DOIs below are copied from the verified literature review (`reports/lit_review/ovarian-cell-line-multiomic-resource-literature-review.md`); all were existence-verified there (OpenAlex/DOI), and none is retracted.*

- Ahmed AA, et al. (2010). Driver mutations in TP53 are ubiquitous in high grade serous carcinoma of the ovary. *J Pathol* 221:49–56. 10.1002/path.2696
- Anglesio MS, et al. (2013). Type-specific cell line models for type-specific ovarian cancer research. *PLoS One* 8:e72162. 10.1371/journal.pone.0072162
- Barretina J, et al. (2012). The Cancer Cell Line Encyclopedia. *Nature* 483:603–607. 10.1038/nature11003
- Beaufort CM, et al. (2014). Ovarian cancer cell line panel (OCCP). *PLoS One* 9:e103988. 10.1371/journal.pone.0103988
- Cancer Genome Atlas Research Network (2011). Integrated genomic analyses of ovarian carcinoma. *Nature* 474:609–615. 10.1038/nature10166
- Cheasley D, et al. (2019). The molecular origin and taxonomy of mucinous ovarian carcinoma. *Nat Commun* 10:3935. 10.1038/s41467-019-11862-x
- Chen GM, et al. (2018). Consensus on molecular subtypes of high-grade serous ovarian carcinoma. *Clin Cancer Res* 24:5037–5047. 10.1158/1078-0432.CCR-18-0784
- Coscia F, et al. (2016). Integrative proteomic profiling of ovarian cancer cell lines. *Nat Commun* 7:12645. 10.1038/ncomms12645
- Domcke S, et al. (2013). Evaluating cell lines as tumour models by comparison of genomic profiles. *Nat Commun* 4:2126. 10.1038/ncomms3126
- Fleury H, et al. (2015). Novel high-grade serous epithelial ovarian cancer cell lines (BRCA1/BRCA2). *Genes Cancer* 6:378–398. 10.18632/genesandcancer.76
- Fleury H, et al. (2016). Cumulative defects in DNA repair pathways and olaparib. *Oncotarget* 7:40152–40168. 10.18632/oncotarget.10308
- Gamwell LF, et al. (2013). The SCCOHT cell line BIN-67. *Orphanet J Rare Dis* 8:33. PMC3635907
- Ghandi M, et al. (2019). Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature* 569:503–508. 10.1038/s41586-019-1186-3
- Gonçalves E, et al. (2022). Pan-cancer proteomic map of 949 human cell lines. *Cancer Cell* 40:835–849. 10.1016/j.ccell.2022.06.010
- Halperin RF, et al. (2017). A method to reduce ancestry related germline false positives in tumor only somatic variant calling. *BMC Med Genomics* 10:61. 10.1186/s12920-017-0296-8
- Hamilton E, et al. (2020). Tamrintamab pamozirine (SC-003, anti-DPEP3). *Gynecol Oncol* 158:640–645. 10.1016/j.ygyno.2020.05.038
- Hollis RL, et al. (2020). Molecular stratification of endometrioid ovarian carcinoma. *Nat Commun* 11:4995. 10.1038/s41467-020-18819-5
- Huang Y, et al. (2017). Authentication of 278 cancer cell lines. *PLoS One* 12:e0170384. 10.1371/journal.pone.0170384
- ICLAC (2026). Register of Misidentified Cell Lines, v14. iclac.org
- Jarnuczak AF, et al. (2021). An integrated landscape of protein expression in human cancer. *Sci Data* 8:115. 10.1038/s41597-021-00890-2
- Jelinic P, et al. (2014). Recurrent SMARCA4 mutations in SCCOHT. *Nat Genet* 46:424–426. 10.1038/ng.2922
- Karnezis AN, et al. (2016). Dual loss of SMARCA4/SMARCA2 in SCCOHT. *J Pathol* 238:389–400. 10.1002/path.4633
- Karnezis AN, et al. (2021). Re-assigning the histologic identities of COV434 and TOV-112D. *Gynecol Oncol* 160:568–578. 10.1016/j.ygyno.2020.12.004
- Korch C, et al. (2012). DNA profiling analysis of endometrial and ovarian cell lines. *Gynecol Oncol* 127:241–248. 10.1016/j.ygyno.2012.06.017
- Kuo KT, et al. (2009). Frequent activating mutations of PIK3CA in ovarian clear cell carcinoma. *Am J Pathol* 174:1597–1601. 10.2353/ajpath.2009.081000
- Létourneau IJ, et al. (2012). Derivation and characterization of matched cell lines from primary and recurrent serous ovarian cancer. *BMC Cancer* 12:379. 10.1186/1471-2407-12-379
- Little P, et al. (2021). UNMASC: tumor-only variant calling with unmatched normal controls. *NAR Cancer* 3:zcab040. 10.1093/narcan/zcab040
- Matulonis UA, et al. (2023). Mirvetuximab soravtansine in FRα-high platinum-resistant ovarian cancer (SORAYA). *J Clin Oncol* 41. 10.1200/JCO.22.01900
- Meagher NS, et al. (2025). Cellular origins of mucinous ovarian carcinoma. *J Pathol*. 10.1002/path.6407
- Meric-Bernstam F, et al. (2024). Trastuzumab deruxtecan in HER2-expressing solid tumors (DESTINY-PanTumor02). *J Clin Oncol* 42:47–58. 10.1200/JCO.23.02005
- Moore KN, et al. (2023). Mirvetuximab soravtansine in FRα-positive platinum-resistant ovarian cancer (MIRASOL). *N Engl J Med* 389:2162–2174. 10.1056/NEJMoa2309169
- Nusinow DP, et al. (2020). Quantitative proteomics of the Cancer Cell Line Encyclopedia. *Cell* 180:387–402. 10.1016/j.cell.2019.12.023
- Olbrecht S, et al. (2021). HGSC refined with single-cell RNA sequencing. *Genome Med* 13:111. 10.1186/s13073-021-00922-x
- Ouellet V, et al. (2008). Characterization of three new serous epithelial ovarian cancer cell lines. *BMC Cancer* 8:152. 10.1186/1471-2407-8-152
- Peng G, et al. (2014). Genome-wide transcriptome profiling of homologous recombination DNA repair. *Nat Commun* 5:3361. 10.1038/ncomms4361
- Provencher DM, et al. (2000). Characterization of four novel epithelial ovarian cancer cell lines. *In Vitro Cell Dev Biol Anim* 36:357–361. 10.1290/1071-2690(2000)036<0357:COFNEO>2.0.CO;2
- Ramos P, et al. (2014). SMARCA4 inactivating mutations in SCCOHT. *Nat Genet* 46:427–429. 10.1038/ng.2928
- Sauriol A, et al. (2020). Modeling the diversity of epithelial ovarian cancer through ten novel cell lines. *Cancers* 12:2222. 10.3390/cancers12082222
- Schwede M, et al. (2020). The impact of stroma admixture on molecular subtypes and prognostic gene signatures in serous ovarian cancer. *Cancer Epidemiol Biomarkers Prev* 29:509–519. 10.1158/1055-9965.EPI-18-1359
- Shrestha R, et al. (2021). Multiomic characterization of low-grade serous ovarian carcinoma cell lines. *Cancer Res* 81:1681–1694. 10.1158/0008-5472.CAN-20-2222
- Stordal B, et al. (2024). Authenticating ovarian cancer cell lines by histotype. *Mol Biol Rep* 51:784. 10.1007/s11033-024-09747-4
- Sztupinszki Z, et al. (2018). Migrating SNP-array HRD measures to NGS (scarHRD). *npj Breast Cancer* 4:16. 10.1038/s41523-018-0066-6
- Takamatsu S, et al. (2024). Homologous recombination deficiency unrelated to platinum and PARP inhibitor response in cell line libraries. *Sci Data* 11:171. 10.1038/s41597-024-03018-4
- Talevich E, et al. (2016). CNVkit: genome-wide copy number detection from targeted sequencing. *PLoS Comput Biol* 12:e1004873. 10.1371/journal.pcbi.1004873
- Upadhya SR, Ryan CJ (2022). Experimental reproducibility limits mRNA–protein correlations. *Cell Rep Methods* 2:100288. 10.1016/j.crmeth.2022.100288
- van den Berg-Bakker CAM, et al. (1993). Establishment and characterization of 7 ovarian carcinoma cell lines and one granulosa tumor cell line (incl. COV434). *Int J Cancer* 53:613–620. 10.1002/ijc.2910530415
- Wiegand KC, et al. (2010). ARID1A mutations in endometriosis-associated ovarian carcinomas. *N Engl J Med* 363:1532–1543. 10.1056/NEJMoa1008433
- Witkowski L, et al. (2014). Germline and somatic SMARCA4 mutations characterize SCCOHT. *Nat Genet* 46:438–443. 10.1038/ng.2931
- Wong NKY, et al. (2025). Modeling gynecologic carcinosarcoma (PDX + cell line). *Transl Oncol* 63:102591.
- Zhang H, et al. (2016). Integrated proteogenomic characterization of human high-grade serous ovarian cancer (CPTAC). *Cell* 166:755–765. 10.1016/j.cell.2016.05.069
- Zhang Q, Wang C, Cliby WA (2019). Cancer-associated stroma contributes to the mesenchymal subtype signature. *Gynecol Oncol* 152:368–374. 10.1016/j.ygyno.2018.11.014

<!-- ============================================================
DRAFT NOTES FOR PI (not part of the manuscript body)

PLACEHOLDERS to resolve (deposition; PI handles):
- All repository accessions in Table 1 (GEO/SRA, ProteomeXchange/PRIDE, figshare/Zenodo DOI).
- WES raw-reads deposit is contingent on FASTQ/BAM availability (see Data Records note).
- Raw MS from Morin lab for ProteomeXchange (see Methods/Proteomics note).

CHECK items (numbers/parameters to confirm at Methods stage):
- RNA-seq read configuration and library-prep kit (run NS.1676.003).
- WES capture-kit identity.
(All quantitative results in the body are sourced from reports/01_multiomic_characterization_results.md and output/ CSVs.)

NOTE items (genuine gaps / decisions flagged in body):
- Table numbering: the cohort "Table 1" from the outline is replaced by Fig. 1 per the figures plan (§B), so the Data Records table is Table 1 (renumbered from the outline's Table 2). Confirm.
- Genomic HRD not computable from archived data (§4 note) — PI decision on BAM recovery.
- On-target/coverage WES QC (Fig. S8) contingent on recoverable raw data.
- Deposition hygiene: MAF NCBI_Build header correction; strip third-party normals from CNVkit outputs.
- renv lockfile + repository host to finalize (Code Availability).
- Authorship list per PROJECT_SPEC to confirm.
============================================================ -->
