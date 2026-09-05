# Lit Review — Cluster Q3/Q4: Pan-cancer cell-line multi-omic resources & Scientific Data Descriptor norms

Specialist: lit-search (CCLE / proteomics / Data Descriptor norms)
Date: 2026-07-23
Scope: (Q3) large-scale cancer cell-line multi-omic resources + RNA–protein concordance benchmarks; (Q4) Scientific Data Data Descriptor precedents + editorial expectations.
Purpose: benchmark an OvCAN ovarian cell-line resource (bulk RNA-seq + TMT proteomics + WES, ~31 lines, 7 subtypes) against precedent.

> **IMPORTANT SEED CORRECTION (read first).** The seed listed *"Frejno et al. 2021 Scientific Data … 'integrated landscape of protein expression in human cancer' (191 lines + 246 tumours)."* The paper with that title, DOI, and those exact numbers is **Jarnuczak et al. 2021** (a Vizcaíno/EMBL-EBI group meta-analysis), **not Frejno**. **Frejno et al.** is a *different, real* paper — "Proteome activity landscapes of tumor cell lines determine drug responses," **2020 Nature Communications** (125 cell lines). Both are reported below and disambiguated. Attribute the "integrated landscape" 0.58 correlation to **Jarnuczak**, not Frejno.

---

## Findings

### Q3 — Cancer cell-line multi-omic resources & RNA–protein correlation benchmarks

### CCLE genomics — original (Barretina 2012)
- **Citation:** Barretina J, Caponigro G, Stransky N, et al. (2012). "The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity." *Nature* 483:603–607. DOI: 10.1038/nature11003. https://www.nature.com/articles/nature11003
- **Key finding:** Foundational CCLE resource: **947 human cancer cell lines** characterized for gene expression, chromosomal copy number, and massively parallel sequencing, screened against 24 anticancer compounds to build predictive models of drug sensitivity. Establishes the "large cell-line panel + multi-omic characterization" template the field now expects.
- **Relevance:** Medium — historical anchor / lineage citation for any cell-line multi-omic resource; not proteomics.
- **Source type:** Primary

### CCLE genomics — next-generation characterization (Ghandi 2019)
- **Citation:** Ghandi M, Huang FW, Jané-Valbuena J, et al. (2019). "Next-generation characterization of the Cancer Cell Line Encyclopedia." *Nature* 569:503–508. DOI: 10.1038/s41586-019-1186-3. (search-confirmed title/journal/year/DOI)
- **Key finding:** Expanded multi-omic CCLE; per search summary, **RNA-Seq profiling of 1,076 cell lines** plus additional layers (WES/WGS, RRBS methylation, miRNA, RPPA, metabolomics, etc.). Frames the current multi-omic CCLE that downstream proteomics (Nusinow) was layered onto. **Note:** exact per-assay line counts were taken from the search summary, not the primary PDF (see Gaps).
- **Relevance:** High — the canonical multi-omic cell-line reference; direct comparator for "how many lines / how many data layers."
- **Source type:** Primary

### CCLE proteomics — Nusinow 2020 (the key mRNA–protein benchmark)
- **Citation:** Nusinow DP, Szpyt J, Ghandi M, et al. (2020). "Quantitative Proteomics of the Cancer Cell Line Encyclopedia." *Cell* 180(2):387–402.e16. DOI: 10.1016/j.cell.2019.12.023. https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6
- **Key finding:** TMT (10-plex, SPS-MS3 on Orbitrap Fusion Lumos) quantitative proteomics of **375 cell lines across 22 lineages**; **>9,000 proteins quantified on average per experiment** (42 multiplex experiments, 504 MS runs, >1,500 instrument-hours). Original paper reports **average per-gene mRNA–protein correlation ≈ 0.5**. Gene sets at plasma membrane / ECM correlate higher; many protein complexes correlate lower.
- **Data accession:** MassIVE **MSV000085836** (ProteomeXchange partner).
- **Relevance:** High — the single most-cited RNA-vs-protein concordance benchmark for cell lines; TMT chemistry matches the OvCAN proteomics assay.
- **Source type:** Primary

### CCLE mRNA–protein correlation — precise re-analysis (Upadhya & Ryan 2022)
- **Citation:** Upadhya SR, Ryan CJ (2022). "Experimental reproducibility limits the correlation between mRNA and protein abundances in tumor proteomic profiles." *Cell Reports Methods* 2(9):100288. DOI: 10.1016/j.crmeth.2022.100288. https://pmc.ncbi.nlm.nih.gov/articles/PMC9499981/
- **Key finding:** Re-analysis reports the **CCLE median mRNA–protein correlation ≈ 0.48**, versus a **replicate-proteome reproducibility ≈ 0.72** — i.e., measurement reproducibility (not just post-transcriptional regulation) caps the achievable correlation. Proteins with more reproducible measurements have higher mRNA–protein correlation. Gives a defensible **~0.5 target and a reproducibility ceiling** to benchmark OvCAN against.
- **Relevance:** High — exact number + the correct interpretive framing (don't over-read a "modest" 0.5; it's partly a measurement ceiling).
- **Source type:** Primary (methods/re-analysis)

### ProCan pan-cancer proteomic map (Gonçalves 2022)
- **Citation:** Gonçalves E, Poulos RC, Cai Z, et al. (2022). "Pan-cancer proteomic map of 949 human cell lines." *Cancer Cell* 40(8):835–849.e8. DOI: 10.1016/j.ccell.2022.06.010. https://pmc.ncbi.nlm.nih.gov/articles/PMC9387775/
- **Key finding:** **949 cell lines across 28 tissues / >40 cancer types**; **8,498 proteins quantified (6,692 supported by >1 peptide)** by **DIA-SWATH-MS** (SCIEX 6600 TripleTOF, 100 variable windows, 90-min gradients). **Median protein-wise mRNA–protein Pearson r = 0.42** (cell-type-enriched proteins higher). Protein networks are more strongly co-regulated than transcriptomics; proteome predicts drug response ~ as well as transcriptome.
- **Data accession:** PRIDE **PXD030304** (raw DIA-MS + spectral library); portal https://cellmodelpassports.sanger.ac.uk; drug data at cancerRxgene.org.
- **Relevance:** High — largest cell-line proteomic map; DIA benchmark (contrast with OvCAN's TMT); provides the lower end (~0.42) of the correlation range.
- **Source type:** Primary

### Integrated landscape of protein expression — Jarnuczak 2021 (NOT Frejno) — Scientific Data descriptor
- **Citation:** Jarnuczak AF, Najgebauer H, Barzine M, Kundu DJ, Ghavidel F, Perez-Riverol Y, Papatheodorou I, Brazma A, Vizcaíno JA (2021). "An integrated landscape of protein expression in human cancer." *Scientific Data* 8:115. DOI: 10.1038/s41597-021-00890-2. https://pmc.ncbi.nlm.nih.gov/articles/PMC8065022/
- **Key finding:** Meta-analysis reprocessing **11 large-scale quantitative proteomics datasets (7,171 MS runs)** into a reference map of **191 cancer cell lines + 246 clinical tumours across 13 lineages**; **15,443 gene products with ≥1 unique peptide (67.8% of the human reference proteome)**. **mRNA–protein correlation: median 0.58 (range 0.43–0.66) across cell lines.** Technical validation included correlating normalized quantities against independent absolute protein concentrations (25 proteins × 6 lines; Spearman ≈ 0.82).
- **Data accession:** PRIDE **PXD013455**; Expression Atlas **E-PROT-18…E-PROT-28**; transcriptomics **E-MTAB-2706 / -2770 / -3983** (ArrayExpress).
- **Relevance:** High — doubles as a Q3 correlation benchmark (the 0.58/0.43–0.66 the seed cited) AND a Q4 *Scientific Data proteomics Data Descriptor* precedent (Technical Validation + Data Records structure).
- **Source type:** Primary / Data Descriptor

### Proteome activity landscapes — the actual Frejno et al. paper (2020)
- **Citation:** Frejno M, Meng C, Ruprecht B, et al. (2020). "Proteome activity landscapes of tumor cell lines determine drug responses." *Nature Communications* 11:3639. DOI: 10.1038/s41467-020-17336-9. https://www.nature.com/articles/s41467-020-17336-9
- **Key finding:** Integrated quantitative data for **~10,000 proteins and ~55,000 phosphosites across 125 cancer cell lines** ("proteome activity landscapes") and linked them to drug response. This is the real "Frejno" paper — related to but distinct from the Jarnuczak "integrated landscape" descriptor above.
- **Relevance:** Medium — supporting cell-line proteomics resource; disambiguates the seed. Did not retrieve its specific mRNA–protein correlation value (see Gaps).
- **Source type:** Primary

### Q4 — Data Descriptor precedents & editorial norms

### Breast-cancer cell-line TMT proteomics Data Descriptor (Kalocsay 2023) — closest precedent
- **Citation:** Kalocsay M, Berberich MJ, Everley RA, et al. (2023). "Proteomic profiling across breast cancer cell lines and models." *Scientific Data* 10:514. DOI: 10.1038/s41597-023-02355-0. https://pmc.ncbi.nlm.nih.gov/articles/PMC10403526/
- **Key finding (QC / normalization detail — the requested specifics):**
  - **60 breast-cancer cell lines/models**, quantified to a depth of **~13,000 proteins**; plus 5 biological replicates of MCF 10A.
  - **TMT10plex and TMT11plex.** **Bridge channel = a single mixed sample derived from six cell lines**, included in every TMT set. Two-step normalization: **within-set normalization to the bridge sample**, then **across-set normalization to a reference set**; protein-level values by peptide averaging; intensities scaled 0–100.
  - **QC metrics:** missed-cleavage threshold **<15%**; **TMT labeling efficiency >95%**; peptide filtered to **FDR <1%**; **protein-level FDR <1%**; **isolation specificity (purity) filter >0.7**; **MCF 10A replicate correlation = 0.719** at protein level.
  - **Data Records / accessions:** PRIDE **PXD026581**; Synapse **syn32672593** (Level 1–3 data); Figshare 10.6084/m9.figshare.c.6443633.v2; distributed via the **LINCS** portal.
- **Relevance:** High — the direct template for the OvCAN TMT descriptor: bridge-channel design, IRS-style within→across-set normalization, and the exact QC panel reviewers will expect (labeling efficiency, FDR, isolation purity, replicate correlation).
- **Source type:** Data Descriptor (Primary)

### LL-100 cell-line panel — WES + RNA-seq Data Descriptor precedent (Quentmeier 2019)
- **Citation:** Quentmeier H, Pommerenke C, Dirks WG, Eberth S, Koeppel M, MacLeod RAF, Nagel S, Steube K, Uphoff CC, Drexler HG (2019). "The LL-100 panel: 100 cell lines for blood cancer studies." *Scientific Reports* 9:8218. DOI: 10.1038/s41598-019-44491-x. https://pmc.ncbi.nlm.nih.gov/articles/PMC6547646/
- **Key finding (technical-validation template for WES+RNA-seq cell-line panels):** **100 leukemia/lymphoma lines across 22 entities.** Assays: **WES (>10M reads, 2×151 bp, ≥50× coverage)**, **RNA-seq (>29M reads, 2×151 bp)**, plus microarray. QC/validation: **cell-line authentication by DNA (STR) profiling + cytogenetics**; mycoplasma/virus screening; trimming (fastq-mcf) + FastQC; **STAR v2.5.3a** alignment (Gencode v26); variant calling with **GATK HaplotypeCaller + VarScan** (quality ≥20, depth ≥10, VAF ≥0.2); CNAs via **control-FREEC v11.0**. Accessions: ENA **PRJEB30297** (WES) / **PRJEB30312** (RNA-seq); ArrayExpress **E-MTAB-7722** (WES) / **E-MTAB-7721** (RNA-seq).
- **Relevance:** High — nearest analog to OvCAN's RNA-seq + WES arms; shows the expected authentication + read-depth + alignment + variant-calling validation stack and repository choices (ENA/ArrayExpress). (Scientific Reports, not Scientific Data, but same publisher family / descriptor style.)
- **Source type:** Primary (resource paper)

### WES/WGS reference-dataset Data Descriptor (SEQC2, Scientific Data 2021)
- **Citation:** *Scientific Data* (2021). "Whole genome and exome sequencing reference datasets from a multi-center and cross-platform benchmark study." DOI: 10.1038/s41597-021-01077-5. https://www.nature.com/articles/s41597-021-01077-5
- **Key finding:** SEQC2-consortium paired tumor–normal reference WGS/WES generated across **16 library protocols, 7 platforms, 6 centers** — a *Scientific Data* precedent for how sequencing (incl. WES) Data Descriptors present cross-protocol technical validation and reference-material QC.
- **Relevance:** Medium — precedent for WES technical-validation framing / benchmarking language (relevant to OvCAN's tumor-only WES arm).
- **Source type:** Data Descriptor
- **Note:** citation captured from search result; full author list not fetched (see Gaps).

### Scientific Data — editorial expectations for a Data Descriptor
- **Citation:** Scientific Data — Submission Guidelines, "For Referees," and "The Data Descriptor – making your data reusable." https://www.nature.com/sdata/submission-guidelines ; https://www.nature.com/sdata/policies/for-referees ; https://blogs.nature.com/scientificdata/2013/09/19/the-data-descriptor-making-your-data-reusable/
- **Key finding:**
  - Data Descriptor has defined sections: **Background & Summary, Methods, Data Records, Technical Validation, Usage Notes, Code Availability** (plus standard front/back matter). **Methods has no length limit** — full experimental design, acquisition, and computational processing so others can reproduce.
  - **Acceptance is NOT based on perceived impact or novelty.** Descriptors "will not be expected to contain in-depth analyses or new scientific conclusions"; authors must instead **support the rigour and technical quality** of the data. Peer review stays focused on **data quality and reusability**, not interpretation.
  - Data must be **deposited in a public, community-endorsed repository**: sequencing → GEO / SRA / ArrayExpress (now BioStudies) / ENA (raw reads → ENA); proteomics → **ProteomeXchange / PRIDE (or MassIVE)** with identifications, PTMs, and supporting spectra.
- **Relevance:** High — defines exactly the machinery reviewers expect and confirms the "minimal interpretation, maximal validation/reusability" posture for the OvCAN descriptor.
- **Source type:** Editorial policy (Primary)

---

## Benchmark summary tables (for quick reference)

**RNA–protein correlation benchmarks (cell lines):**
| Resource | Method | Lines | Proteins | mRNA–protein correlation |
|---|---|---|---|---|
| CCLE / Nusinow 2020 | TMT-MS3 | 375 | ~9,000 avg/expt | ~0.5 (median ≈ 0.48 per Upadhya & Ryan 2022; ceiling ~0.72 reproducibility) |
| ProCan / Gonçalves 2022 | DIA-SWATH | 949 | 8,498 (6,692 >1 pep) | median protein-wise Pearson r = 0.42 |
| Integrated landscape / Jarnuczak 2021 | meta (11 datasets) | 191 (+246 tumours) | 15,443 gene products | median 0.58 (range 0.43–0.66) |
| Frejno 2020 | proteome activity | 125 | ~10,000 (+55k phosphosites) | not retrieved |

**Takeaway for OvCAN:** a cell-line RNA–protein correlation in the **~0.4–0.6** band is the published norm; ~0.5 is squarely expected and partly reproducibility-limited.

---

## Search Log
- Query 1: "Ghandi 2019 Nature next-generation characterization CCLE proteomics mRNA" — 8 examined, 3 relevant
- Query 2: "Nusinow 2020 Cell Quantitative Proteomics CCLE mRNA protein correlation number proteins" — 8 examined, 4 relevant
- Query 3: "Gonçalves 2022 Cancer Cell pan-cancer proteomic map 949 human cell lines ProCan mRNA protein correlation" — 10 examined, 5 relevant
- Query 4: "Frejno 2021 integrated landscape of protein expression in human cancer 191 cell lines mRNA protein correlation" — 6 examined, 3 relevant (surfaced the Jarnuczak attribution)
- Query 5: "breast cancer cell line TMT proteomics Scientific Data 2023 s41597-023-02355-0 data descriptor" — 10 examined, 4 relevant
- Query 6: "Barretina 2012 Nature CCLE predictive modeling drug sensitivity number cell lines" — 9 examined, 2 relevant
- Query 7: "Scientific Data data descriptor technical validation reusability editorial requirements author guidelines" — 8 examined, 5 relevant
- Query 8: "Nusinow 2020 … DOI 10.1016 PMC accession ProteomeXchange PXD" — 9 examined, 3 relevant (MassIVE task URL)
- Query 9: "Frejno 'proteome activity landscapes' vs 'integrated landscape' 2021 Jarnuczak Vizcaino distinguish" — 9 examined, 3 relevant (confirmed disambiguation)
- Query 10: "CCLE proteomics Nusinow mRNA protein correlation median Pearson 0.48" — 9 examined, 3 relevant (0.48 vs 0.72)
- Query 11: "Scientific Data data descriptor cell line panel WES / multi-omics technical validation accession" — 8 examined, 4 relevant (LL-100, SEQC2)
- Query 12: "'Data Descriptor' cell lines RNA-seq proteomics GEO ArrayExpress PRIDE reusability" — 8 examined, 3 relevant
- Query 13: "'Experimental reproducibility limits the correlation…' Cell Reports Methods 2022 Upadhya Ryan DOI" — 8 examined, 2 relevant

WebFetch (numbers/verification): ProCan PMC9387775 ✓; Jarnuczak PMC8065022 ✓; Kalocsay breast-cancer descriptor PMC10403526 ✓; LL-100 PMC6547646 ✓; MassIVE dataset page (Nusinow MSV000085836) ✓. Blocked/failed: cell.com fulltext (403), nature.com article/policy pages (redirect to idp.nature.com login), digitalcommons PDF (403), bioRxiv fulltext (403), web.archive.org (unavailable) — routed around all via PMC/search.

## Gaps
- **Frejno 2020 mRNA–protein correlation value** not retrieved (only protein/phosphosite counts and line count). If needed, fetch Nature Communications 11:3639 directly.
- **Ghandi 2019 exact per-assay line counts** (RNA-seq vs WES/WGS/RPPA/etc.) taken from search summary, not the primary Nature PDF; the "1,076 RNA-seq" figure should be confirmed against the paper before quoting in the descriptor.
- **SEQC2 WES/WGS descriptor**: only title/DOI captured from search; full author list and specific QC metrics not fetched.
- **Scientific Data "For Referees" exact checklist wording** not fetched (page requires login; general criteria captured from the public submission-guidelines/blog snippets). The substantive requirements (data quality + reusability, not novelty; repository deposition; Technical Validation section) are well-supported.
- Did not verify a separate **Frejno et al. 2017 Mol Syst Biol** ("Pharmacoproteomic characterisation of colon/rectal cancer") that the seed's "Mol Syst Biol" mention may point to — flag if the lead wants it traced.

## Breadth Flag
All primary seed papers (Barretina, Ghandi, Nusinow, ProCan/Gonçalves, "integrated landscape") located and verified, plus the requested breast-cancer TMT descriptor with full QC. Beyond the seeds I added: the Upadhya & Ryan (2022) reproducibility re-analysis (gives the precise CCLE 0.48 vs 0.72 ceiling), the LL-100 WES+RNA-seq cell-line panel (direct analog to OvCAN's sequencing arms), and the SEQC2 WES/WGS reference descriptor. Coverage of Q3 correlation benchmarks and Q4 descriptor/editorial norms is strong. Main residual: the two Frejno-authored papers' internal correlation numbers were not the focus and remain partially open (see Gaps). Recommend one follow-up cluster if the lead wants an *ovarian-specific* cell-line RNA/proteomics descriptor precedent (none surfaced; likely a genuine gap that strengthens the OvCAN paper's novelty-of-resource framing).
