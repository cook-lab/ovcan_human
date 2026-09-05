# Co-citation Chain — Explore/Verify Pass (09)

**Task:** Verification + web-grounding of 9 papers surfaced by OpenAlex co-citation analysis as foundational works missed by keyword search. Every finding below is drawn from a WebSearch/WebFetch result (not memory). Each paper reported individually.

**Date:** 2026-07-23
**Result:** 9/9 verified. 0 dropped.

---

## Findings

### Paper 1 — TP53 ubiquity in HGSC (genomics validation)
- **Citation:** Ahmed AA, Etemadmoghadam D, Temple J, Lynch AG, Riad M, Sharma R, Stewart C, Fereday S, Caldas C, Defazio A, Bowtell D, Brenton JD. "Driver mutations in TP53 are ubiquitous in high grade serous carcinoma of the ovary." *J Pathol.* 2010 May;221(1):49–56. DOI 10.1002/path.2696. PMID 20229506.
- **Verification:** Verified via DOI + title + PMID (Wiley/PubMed/PMC3262968 all consistent).
- **Key finding (web-grounded):** Sequenced exons 2–11 and intron–exon boundaries of *TP53* in tumour DNA from 145 patients. Pathogenic/driver *TP53* mutations identified in **96.7% (119/123)** of high-grade (poorly differentiated) serous carcinoma cases — establishing near-ubiquity of *TP53* mutation as the defining genomic event of HGSC.
- **Target section:** Genomics validation: TP53 ubiquity.

### Paper 2 — CPTAC/TCGA proteogenomics benchmark
- **Citation:** Zhang H, Liu T, Zhang Z, et al. (CPTAC Investigators). "Integrated Proteogenomic Characterization of Human High-Grade Serous Ovarian Cancer." *Cell.* 2016 Jul 28;166(3):755–765. DOI 10.1016/j.cell.2016.05.069. PMID 27372738.
- **Verification:** Verified via DOI + title + PMID (Cell/ScienceDirect S0092867416306730 + PubMed abstract).
- **Key finding (web-grounded):** Mass-spectrometry proteomic characterization of **174 ovarian tumours previously profiled by TCGA, of which 169 were HGSC**. mRNA–protein integration: **90.6% of gene–protein pairs showed positive correlation and 79.4% were significantly correlated** (mean Spearman ρ = 0.38, median 0.45); metabolic pathways and interferon response correlated strongly, ribosome and mRNA-splicing poorly. Copy-number–protein integration mapped **29,393 CNAs against 3,202 proteins**, revealing which CNAs propagate to protein abundance. Derived protein signatures (from CNA-affected proteins) that predicted overall survival in independent datasets. Linked **histone H4 acetylation (K12, K16) to homologous-recombination deficiency**, proposing acetylation state as a therapy-stratification marker.
- **Target section:** Proteogenomics benchmark (CPTAC/TCGA ovarian).

### Paper 3 — PIK3CA in ovarian clear cell carcinoma
- **Citation:** Kuo KT, Mao TL, Jones S, et al. "Frequent Activating Mutations of PIK3CA in Ovarian Clear Cell Carcinoma." *Am J Pathol.* 2009 May;174(5):1597–1601. DOI 10.2353/ajpath.2009.081000. PMID 19349352.
- **Verification:** Verified via DOI + title + PMID (Am J Pathol / Johns Hopkins Pure / PubMed).
- **Key finding (web-grounded):** Analyzed **97 ovarian clear cell carcinomas (CCCs)** (18 affinity-purified fresh tumours, 69 microdissected FFPE, 10 cell lines) for KRAS, BRAF, PIK3CA, TP53, PTEN, CTNNB1. Mutation frequencies: **PIK3CA 33%**, TP53 15%, KRAS 7%, PTEN 5%, CTNNB1 3%, BRAF 1%. In the subset of 28 affinity-purified CCCs + CCC cell lines (higher tumour purity), **PIK3CA mutation frequency rose to 46%** — establishing PIK3CA/PI3K-pathway activation as a hallmark of clear cell carcinoma.
- **Target section:** Clear cell biology — PIK3CA mutation frequency.

### Paper 4 — Mucinous ovarian cancer as a separate entity
- **Citation:** Hess V, A'Hern R, et al. "Mucinous Epithelial Ovarian Cancer: A Separate Entity Requiring Specific Treatment." *J Clin Oncol.* 2004. DOI 10.1200/JCO.2004.08.078. (Verified landing page ascopubs.org/doi/10.1200/JCO.2004.08.078; PubMed record present.)
- **Verification:** Verified via DOI + title (ASCO Publications + PubMed). Exact volume/page/PMID not independently confirmed in this pass — flagged for the citation manager; the substantive finding below is grounded.
- **Key finding (web-grounded):** Case–control comparison of **27 mucinous epithelial ovarian cancer cases vs 54 non-mucinous controls**, all treated with platinum-based regimens. Objective response rate was **26.3% (95% CI 9.2–51.2) for mucinous vs 64.9% (95% CI 47.5–79.8) for controls (P = .01)** — significantly poorer platinum chemoresponse, motivating the argument that mucinous EOC is a distinct entity needing histology-specific treatment.
- **Target section:** Mucinous as a distinct entity — chemoresponse vs serous.

### Paper 5 — Konecny molecular subtype classifier
- **Citation:** Konecny GE, Wang C, Hamidi H, et al. "Prognostic and Therapeutic Relevance of Molecular Subtypes in High-Grade Serous Ovarian Cancer." *J Natl Cancer Inst.* 2014 Oct;106(10):dju249. DOI 10.1093/jnci/dju249. PMID 25269487 (PMC4271115).
- **Verification:** Verified via DOI + title + PMID (Oxford Academic / PubMed / PMC).
- **Key finding (web-grounded):** Gene-expression profiling of **174 HGSOC (Mayo Clinic)** confirmed four transcriptional subtypes de novo: **immunoreactive-like, differentiated-like, proliferative-like, mesenchymal-like**. Overall survival differed significantly (log-rank P = .006): best for immunoreactive-like; significantly worse for **proliferative-like (adjusted HR 1.89, 95% CI 1.18–3.02, P = .008)** and **mesenchymal-like (adjusted HR 2.45, 95% CI 1.43–4.18, P = .001)** vs immunoreactive-like. Replicated in an external cohort of **185 HGSOC (Bonome)**. This is the "Konecny classifier" — one of the four subtype schemes reconciled by consensusOV.
- **Target section:** Molecular subtyping — Konecny classifier.

### Paper 6 — Helland C5/MYCN–LIN28B–let-7 subtype
- **Citation:** Helland Å, Anglesio MS, George J, Cowin PA, Johnstone CN, House CM, Sheppard KE, Etemadmoghadam D, Melnyk N, Rustgi AK, Phillips WA, Johnsen H, Holm R, Kristensen GB, Birrer MJ, Pearson RB, Børresen-Dale AL, Huntsman DG, DeFazio A, Creighton CJ, Smyth GK, Bowtell DDL. "Deregulation of MYCN, LIN28B and LET7 in a Molecular Subtype of Aggressive High-Grade Serous Ovarian Cancers." *PLoS One.* 2011 Apr;6(4):e18064. DOI 10.1371/journal.pone.0018064. (PMC3076323.)
- **Verification:** Verified via DOI + title (PLoS ONE + PMC + Harvard DASH + Monash Research author list).
- **Key finding (web-grounded):** Across **>900 primary serous ovarian tumour samples**, deregulation of the MYCN–LIN28B–let-7 oncogenic stem-cell-renewal axis was found to be **specifically associated with the C5 molecular subtype**, the most aggressive HGSC subtype. First demonstration that a defined oncogenic pathway is selectively altered in one molecular subtype of serous ovarian cancer, nominating it for targeted intervention. (This work underpins the "Helland" subtype scheme referenced alongside Konecny in consensusOV reconciliation.)
- **Target section:** Molecular subtyping — Helland classifier.

### Paper 7 — Langdon PEO/PEA cell lines (historical landmark)
- **Citation:** Langdon SP, Lawrie SS, Hay FG, Hawkes MM, McDonald A, Hayward IP, Schol DJ, Hilgers J, Leonard RC, Smyth JF. "Characterization and Properties of Nine Human Ovarian Adenocarcinoma Cell Lines." *Cancer Res.* 1988 Nov 1;48(21):6166–6172. PMID 3167863. URL https://aacrjournals.org/cancerres/article/48/21/6166/493169
- **Verification:** Verified via PMID + URL (AACR Cancer Research + PubMed).
- **Key finding (web-grounded):** Established and characterized **nine human ovarian adenocarcinoma cell lines**. From poorly differentiated adenocarcinoma (ascites/pleural effusion): **PEO1, PEO4, PEO6 (one patient); PEA1, PEA2 (second patient); PEO16 (third patient)**. From a well-differentiated serous adenocarcinoma patient: **PEO14, PEO23 (ascites) and TO14 (solid metastasis)**. Foundational, still widely used HGSC/serous line panel (PEO1/PEO4 are canonical BRCA2-mutant platinum-sensitive/resistant models).
- **Target section:** Cell-line resources — historical landmark.

### Paper 8 — van den Berg-Bakker COV series (incl. COV434)
- **Citation:** van den Berg-Bakker CAM, Hagemeijer A, Franken-Postma EM, Smit VTHBM, Kuppen PJK, van Ravenswaay Claasen HH, Cornelisse CJ, Schrier PI. "Establishment and characterization of 7 ovarian carcinoma cell lines and one granulosa tumor cell line: growth features and cytogenetics." *Int J Cancer.* 1993 Feb 20;53(4):613–620. DOI 10.1002/ijc.2910530415. PMID 8436435.
- **Verification:** Verified via DOI + title + PMID (culturecollections.org.uk/ECACC profile cites it as the COV434 origin paper; PubMed record consistent).
- **Key finding (web-grounded):** Establishment paper for the Leiden **COV ("Cancer OVary") series** — seven ovarian carcinoma lines plus **one granulosa tumour cell line, COV434**. Per the ECACC cell-line profile (which cites this paper as the origin reference), **COV434 was isolated in 1984 from a primary granulosa cell tumour of a 27-year-old woman with metastatic granulosa cell carcinoma and was originally described as a granulosa cell line.**
- **Tie-in confirmed (Karnezis reclassification):** Karnezis AN, et al. "Re-assigning the histologic identities of COV434 and TOV-112D ovarian cancer cell lines." *Gynecol Oncol.* 2021 Feb. PMID 33328126 (PMC10039450). This paper concludes **COV434 is a bona fide SCCOHT (small cell carcinoma of the ovary, hypercalcemic type) line** based on morphology, IHC, genetics (SMARCA4 loss), and clinical features — COV434 xenografts raise serum calcium, and SMARCA4 re-expression dramatically suppresses growth (TOV-112D reassigned to dedifferentiated carcinoma). So the original granulosa designation from van den Berg-Bakker 1993 has been overturned. Directly supports the review's SCCOHT/COV434 narrative.
- **Target section:** Cell-line resources — COV series; COV434 → SCCOHT reclassification.

### Paper 9 — curatedOvarianData resource
- **Citation:** Ganzfried BF, Riester M, Haibe-Kains B, Risch T, Tyekucheva S, Jazic I, Wang XV, Ahmadifar M, Birrer MJ, Parmigiani G, Huttenhower C, Waldron L. "curatedOvarianData: clinically annotated data for the ovarian cancer transcriptome." *Database (Oxford).* 2013 Apr 2;2013:bat013. DOI 10.1093/database/bat013. PMID 23550061 (PMC3625954).
- **Verification:** Verified via DOI + title (Oxford Academic + PMC + Bioconductor + Waldron Lab).
- **Key finding (web-grounded):** R/Bioconductor experiment-data package providing **uniformly preprocessed microarray gene-expression data + manually curated, documented clinical metadata for 2,970 ovarian cancer patients across 23 studies spanning 11 measurement platforms**, delivered as documented ExpressionSet objects. Purpose-built for reproducible development/validation of prognostic models, investigation of molecular subtypes, and benchmarking of machine-learning algorithms on gene-expression data.
- **Target section:** Ovarian transcriptome data resource.

---

## Search Log
- Query 1: `Ahmed 2010 "Driver mutations in TP53 are ubiquitous in high grade serous carcinoma of the ovary" Journal of Pathology` — ~9 examined, 5 relevant (Wiley, PubMed 20229506, PMC3262968, Semantic Scholar, Google Scholar).
- Query 2: `Zhang 2016 "Integrated Proteogenomic Characterization of Human High-Grade Serous Ovarian Cancer" Cell CPTAC` — ~9 examined, 4 relevant (Cell/ScienceDirect S0092867416306730, PubMed 27372738, PMC).
- Query 3: `Kuo 2009 "Frequent Activating Mutations of PIK3CA in Ovarian Clear Cell Carcinoma" American Journal of Pathology` — ~7 examined, 4 relevant (Am J Pathol, PubMed 19349352, Johns Hopkins Pure, SciRP).
- Query 4: `Hess 2004 "Mucinous Epithelial Ovarian Cancer: A Separate Entity Requiring Specific Treatment" JCO response rate` — ~10 examined, 2 directly relevant (ascopubs DOI landing, PubMed).
- Query 5: `Konecny 2014 "Prognostic and Therapeutic Relevance of Molecular Subtypes in High-Grade Serous Ovarian Cancer" JNCI four subtypes survival` — ~8 examined, 3 relevant (Oxford Academic, PubMed 25269487, PMC4271115).
- Query 6: `Helland 2011 "Deregulation of MYCN, LIN28B and LET7 ..." PLoS One` — ~10 examined, 4 relevant (PLoS ONE, PMC3076323, Harvard DASH, Monash Research).
- Query 7: `Langdon 1988 "Characterization and properties of nine human ovarian adenocarcinoma cell lines" Cancer Research PEO PEA` — ~8 examined, 3 relevant (AACR 48(21):6166, PubMed 3167863, ResearchGate).
- Query 8: `van den Berg-Bakker 1993 "Establishment and characterization of 7 ovarian carcinoma cell lines and one granulosa tumor cell line" IJC COV434` — ~9 examined, 5 relevant (ECACC profile, PubMed 8436435, Karnezis reclassification PMID 33328126).
- Query 9: `Ganzfried 2013 curatedOvarianData "clinically annotated data for the ovarian cancer transcriptome" Database Oxford` — ~9 examined, 5 relevant (Oxford Academic, PMC3625954, dblp, Bioconductor, Waldron Lab).
- Query 10 (fetch): PubMed 27372738 abstract — Zhang 2016 quantitative integration numbers.
- Query 11 (fetch): PMC4271115 — Konecny 2014 four subtype names + hazard ratios.
- Query 12 (fetch): ECACC COV434 profile — origin paper + granulosa designation.
- Query 13: `"Re-assigning the histologic identities of COV434 and TOV-112D" Karnezis 2021 SCCOHT` — ~7 examined, 3 relevant (ScienceDirect, PubMed 33328126, PMC10039450).

## Unverifiable / dropped
- None. All 9 papers verified.

## Flags for citation manager
- Paper 4 (Hess 2004): DOI 10.1200/JCO.2004.08.078 and the substantive response-rate finding are grounded, but exact **volume/issue/page and PMID were not independently confirmed** in this pass — verify before final reference list. (Commonly cited as *J Clin Oncol.* 2004;22(6):1040–1044, but confirm.)
- Paper 6 (Helland 2011) and Paper 8 (van den Berg-Bakker 1993) PMIDs (21533229 / 8436435) — PMID for Helland was not directly displayed in the search snippet; confirm from PubMed. van den Berg-Bakker PMID 8436435 was displayed and is confirmed.
