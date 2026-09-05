# Cluster 01 — Ovarian cancer cell-line resources & authentication

Literature search specialist report. Every finding below was verified via WebSearch/WebFetch
(July 2026). Findings are reported individually, not synthesized. Seed papers were confirmed and
expanded from. A curated grey-literature resource (the OvCAN Collection catalogue) was also read
directly because it *is* the resource the Data Descriptor situates.

---

## Findings

### Q1 — Established OC cell-line panels & how they were molecularly characterized

### Domcke et al. 2013 — genomic benchmarking of 47 OC lines vs TCGA HGSOC (seed, confirmed)
- **Citation:** Domcke S, Sinha R, Levine DA, Sander C, Schultz N (2013). "Evaluating cell lines as tumour models by comparison of genomic profiles." *Nature Communications* 4:2126. DOI: 10.1038/ncomms3126 · https://www.nature.com/articles/ncomms3126 · PMID 23839242
- **Key finding:** Compared copy-number, mutation and mRNA profiles of 47 ovarian cancer cell lines against TCGA HGSOC tumours and ranked "suitability." The most commonly used/most-cited lines (notably **SKOV3 and A2780**) had among the *lowest* HGSOC-suitability scores; rarely used lines (**Kuramochi, OVSAHO**, plus SNU-119/COV362/OVCAR4, and possibly CAOV3/OVCAR3/OVCAR8) most closely resembled HGSOC tumours. Foundational argument that the field's workhorse lines are poor HGSOC models.
- **Relevance:** High — the canonical "which lines are good vs poor HGSOC models" reference; directly frames why a rare-subtype-aware, authenticated panel is needed.
- **Source type:** Primary

### Beaufort et al. 2014 — Ovarian Cancer Cell Line Panel (OCCP) (seed, confirmed)
- **Citation:** Beaufort CM, Helmijr JCA, Piskorz AM, Hoogstraat M, Ruigrok-Ritstier K, Besselink N, Murtaza M, van IJcken WFJ, Heine AAJ, Smid M, Koudijs MJ, Brenton JD, Berns EMJJ, Helleman J (2014). "Ovarian Cancer Cell Line Panel (OCCP): Clinical Importance of In Vitro Morphological Subtypes." *PLoS ONE* 9(9):e103988. DOI: 10.1371/journal.pone.0103988 · PMID 25230021
- **Key finding:** Uniformly characterized **39 OC cell lines** (33 from unique patients) for growth, mRNA/miRNA expression, exon sequencing, and drug response. Assigned **14 as high-grade serous, 4 serous-type, 1 low-grade serous, 20 non-serous**. Defined three in-vitro morphologies — Epithelial (n=21), Round (n=7), Spindle (n=12) — that mapped to the TCGA/Tothill C4/C5/C1 expression subtypes; spindle morphology tracked with advanced stage and poor prognosis.
- **Relevance:** High — one of the two most-cited uniform OC cell-line panels; benchmark for panel scope and characterization depth. (Note a 2015 correction/erratum exists: DOI 10.1371/journal.pone.0122284.)
- **Source type:** Primary

### Anglesio et al. 2013 — type-specific (histotype) cell-line models (seed, confirmed)
- **Citation:** Anglesio MS, Wiegand KC, Melnyk N, Chow C, Salamanca C, Prentice LM, Senz J, Yang W, Spillman MA, Cochrane DR, Shumansky K, Shah SP, Kalloger SE, Huntsman DG (2013). "Type-Specific Cell Line Models for Type-Specific Ovarian Cancer Research." *PLoS ONE* 8(9):e72162. DOI: 10.1371/journal.pone.0072162 · PMID 24023729
- **Key finding:** Classified **32 "ovarian cancer" cell lines** into histotypes using mutation profiles, IHC mutation-surrogates and validated immunomarkers. Confirmed **clear-cell** models (e.g., TOV21G, JHOC-5) carry characteristic **ARID1A and PIK3CA** mutations and clear-cell immunoprofiles; explicitly **questioned the use of SKOV3 and A2780 as HGSOC models**. Argues histotype-matched lines must be used for histotype-specific research.
- **Relevance:** High — establishes the rare-subtype (clear-cell) authentication logic and independently flags SKOV3/A2780; directly relevant to a multi-histotype resource.
- **Source type:** Primary

### Ince et al. 2015 — 25 new OCI lines that phenocopy primary tumours (seed, confirmed)
- **Citation:** Ince TA, Sousa AD, Jones MA, Harrell JC, Agoston ES, Krohn M, Selfors LM, Liu W, Chen K, Yong M, Buchwald P, Wang B, Hale KS, Cohick E, Sergent P, Witt A, Kozhekbaeva Z, Gao S, Agoston AT, Merritt MA, Foster R, Rueda BR, Crum CP, Brugge JS, Mills GB (2015). "Characterization of twenty-five ovarian tumour cell lines that phenocopy primary tumours." *Nature Communications* 6:7419. DOI: 10.1038/ncomms8419 · PMID 26080861 · PMC4473807
- **Key finding:** Introduced a defined medium (**"Ovarian Carcinoma Modified Ince medium," OCMI**; OCMIe variant with 17β-estradiol) enabling establishment of **25 new lines from 26 attempts (>95% efficiency**, vs <1% with standard media). Lines retain the genomic landscape, histopathology and molecular features of the parental tumours and stratify with prognostically distinct primary-tumour groups. Subtypes span **papillary serous (most), clear cell, endometrioid, mucinous**, plus rare **carcinosarcoma and dysgerminoma**; exact per-subtype counts not tabulated in the abstract/body text retrieved.
- **Relevance:** High — the "phenocopy" gold-standard resource and a rare direct precedent for deliberately including rare histotypes (carcinosarcoma) in a new panel.
- **Source type:** Primary

### Coscia et al. 2016 — deep proteomic map of OC cell lines (proteomics precedent)
- **Citation:** Coscia F, Watters KM, Curtis M, Eckert MA, Chiang CY, Tyanova S, Montag A, Lastra RR, Lengyel E, Mann M (2016). "Integrative proteomic profiling of ovarian cancer cell lines reveals precursor cell associated proteins and functional status." *Nature Communications* 7:12645. DOI: 10.1038/ncomms12645 · PMID 27561551 · PMC5007461
- **Key finding:** Single-run MS proteomics of **26 OC cell lines** plus HGSOC tumours, immortalized OSE and fallopian-tube epithelium; **>10,000 proteins** quantified. Resolved three cell-line groups (epithelial / clear-cell / mesenchymal) and derived a **67-protein signature** that split both cell lines and the CPTAC/TCGA tumour proteome into epithelial vs mesenchymal HGSOC clusters, suggesting distinct cells of origin.
- **Relevance:** High — the closest published *proteomic* cell-line-map precedent; a direct benchmark for the TMT-proteomics arm of the resource (n, protein depth, integration with CPTAC/TCGA).
- **Source type:** Primary

### Thu et al. 2017 — comprehensively characterized HGSOC-representative panel
- **Citation:** Thu KL, Papari-Zareei M, Stastny V, Song K, Peyton M, Martinez VD, Zhang YA, Castro IB, Varella-Garcia M, Liang H, Xing C, Kittler R, Milchgrub S, Castrillon DH, Davidson HL, Reynolds CP, Lam WL, Lea J, Gazdar AF (2017). "A comprehensively characterized cell line panel highly representative of clinical ovarian high-grade serous carcinomas." *Oncotarget* 8(31):50489–50505. DOI: 10.18632/oncotarget.9929 · PMC5584155
- **Key finding:** Assembled and deeply profiled (mutation, copy-number, expression) a panel selected to match clinical HGSOC molecular features, reinforcing the Domcke-style argument that a specific subset of lines (not the most-cited ones) faithfully represents HGSOC.
- **Relevance:** High — a second independent "HGSOC-representative panel" benchmark for panel-selection criteria.
- **Source type:** Primary

### Haley et al. 2016 — functional characterization of an HGSOC line panel
- **Citation:** Haley J, Tomar S, Pulliam N, Xiong S, Perkins SM, Karpf AR, Mitra S, Nephew KP, Mitra AK (2016). "Functional characterization of a panel of high-grade serous ovarian cancer cell lines as representative experimental models of the disease." *Oncotarget* 7(22):32810–32820. DOI: 10.18632/oncotarget.9053 · PMID 27147568 · PMC5078053
- **Key finding:** Functionally profiled genomically "validated" HGSOC lines — **CAOV3, COV362, Kuramochi, OVCAR4, OVCAR5, OVCAR8, OVSAHO, SNU119** — for migration, invasion, proliferation, clonogenicity, EMT and cisplatin resistance. **OVCAR5, OVCAR8 and Kuramochi** were most aggressive; **SNU119 and OVSAHO** least active. Provides functional annotation on top of genomic suitability.
- **Relevance:** Medium — useful functional companion to Domcke/Thu for HGSOC line selection.
- **Source type:** Primary

### Barnes et al. 2021 — transcriptional stratification of OC lines into 5 histotypes
- **Citation:** Barnes BM, Nelson L, Tighe A, Burghel GJ, Lin IH, Desai S, McGrail JC, Morgan RD, Taylor SS (2021). "Distinct transcriptional programs stratify ovarian cancer cell lines into the five major histological subtypes." *Genome Medicine* 13(1):140. DOI: 10.1186/s13073-021-00952-5 · PMC8408985
- **Key finding:** Used transcriptional profiles to assign OC cell lines to the five major histotypes (HGSOC, LGSOC, ENOC, CCOC, MOC), providing a reference for histotype membership of commonly used lines. (Subtyping methodology itself is another cluster's remit; listed here as a characterization/resource paper.)
- **Relevance:** Medium — histotype-annotation reference for cell lines; complements Anglesio 2013.
- **Source type:** Primary

### CHUM (Mes-Masson / Provencher) TOV/OV line-establishment series

### Provencher et al. 2000 — first TOV/OV lines (seed, confirmed)
- **Citation:** Provencher DM, Lounis H, Champoux L, Tétrault M, Manderson EN, Wang JC, Eydoux P, Savoie R, Tonin PN, Mes-Masson AM (2000). "Characterization of four novel epithelial ovarian cancer cell lines." *In Vitro Cellular & Developmental Biology – Animal* 36(6):357–361. DOI: 10.1290/1071-2690(2000)036<0357:COFNEO>2.0.CO;2 · PMID 10949993
- **Key finding:** Established and characterized **4 lines — TOV-21G (clear cell), TOV-81D, OV-90, TOV-112D** — from chemo/radiotherapy-naive patients, with mutation spectra across BRCA2, TGFβ-RII, KRAS2, TP53, CDKN2A. Founding paper of the CHUM TOV(solid tumour)/OV(ascites) naming convention.
- **Relevance:** High — origin of several lines in the resource (TOV-21G, OV-90, TOV-112D).
- **Source type:** Primary

### Ouellet et al. 2008 — three new serous CHUM lines
- **Citation:** Ouellet V, Zietarska M, Portelance L, Lafontaine J, Madore J, Puiffe ML, Arcand SL, Shen Z, Hébert J, Tonin PN, Provencher DM, Mes-Masson AM (2008). "Characterization of three new serous epithelial ovarian cancer cell lines." *BMC Cancer* 8:152. DOI: 10.1186/1471-2407-8-152 · PMID 18507860 · PMC2467432
- **Key finding:** Derived **3 lines from poorly differentiated serous tumours — TOV-1946, TOV-2223G, and matched ascites OV-1946**; all carried somatic **TP53** mutations (TOV-1946/OV-1946 shared the same mutation); only TOV-1946/OV-1946 formed spheroids and xenograft tumours.
- **Relevance:** High — lines present in the resource (TOV-1946, OV-1946, TOV-2223G per OvCAN catalogue).
- **Source type:** Primary

### Létourneau et al. 2012 — matched primary/recurrent HGSOC lines (seed, confirmed)
- **Citation:** Létourneau IJ, Quinn MCJ, Wang LL, Portelance L, Caceres KY, Cyr L, Delvoye N, Meunier L, de Ladurantaye M, Shen Z, Arcand SL, Tonin PN, Provencher DM, Mes-Masson AM (2012). "Derivation and characterization of matched cell lines from primary and recurrent serous ovarian cancer." *BMC Cancer* 12:379. DOI: 10.1186/1471-2407-12-379 · PMID 22931248
- **Key finding:** Established **9 lines from 3 HGSOC patients** sampled at initial diagnosis and at post-chemo relapse (solid-tumour TOV + ascites OV), characterized for spheroid growth, migration, anchorage independence, in-vivo tumorigenicity and carboplatin/paclitaxel response — a matched pre/post-treatment model set.
- **Relevance:** High — supplies the OV2295/TOV2295(R), OV1369/OV1369(R2), TOV/OV3133 matched series in the resource.
- **Source type:** Primary

### Fleury et al. 2015 — six new HGSOC lines spanning sporadic + hereditary disease
- **Citation:** Fleury H, Communal L, Carmona E, Portelance L, Arcand SL, Rahimi K, Tonin PN, Provencher D, Mes-Masson AM (2015). "Novel high-grade serous epithelial ovarian cancer cell lines that reflect the molecular diversity of both the sporadic and hereditary disease." *Genes & Cancer* 6(9-10):378–398. DOI: 10.18632/genesandcancer.76 · PMID 26622941 · PMC4633166
- **Key finding:** Established **6 lines — TOV2978G, TOV3041G, TOV3291G (tumour) and OV866(2), OV4453, OV4485 (ascites)**. Exome sequencing: somatic **TP53** in five lines; one line with a germline **BRCA1** splice-site mutation and another with a recurrent germline **BRCA2** nonsense mutation — i.e., BRCA-mutant HGSOC models.
- **Relevance:** High — supplies BRCA1/2-mutant HGSOC lines (OV4453, OV4485) in the resource.
- **Source type:** Primary

### Sauriol et al. 2020 — ten novel lines including rarer histotypes
- **Citation:** Sauriol SA, Simeone K, Portelance L, Meunier L, Leclerc-Desaulniers K, de Ladurantaye M, Chergui M, Kendall-Dupont J, Rahimi K, Carmona E, Provencher DM, Mes-Masson AM (2020). "Modeling the Diversity of Epithelial Ovarian Cancer through Ten Novel Well Characterized Cell Lines Covering Multiple Subtypes of the Disease." *Cancers* 12(8):2222. DOI: 10.3390/cancers12082222 · PMID 32784519 · PMC7465288
- **Key finding:** Ten new spontaneously immortalized patient-derived lines with mutational, biomarker, and in-vitro/in-vivo growth data: **8 HGSC, 1 mucinous, 1 clear cell** — an explicit CHUM effort to extend beyond HGSOC into rarer subtypes.
- **Relevance:** High — directly extends CHUM rare-subtype coverage (mucinous TOV2414, clear-cell TOV3392D appear in the resource).
- **Source type:** Primary

### Fleury et al. 2016 — DNA-repair/PARP functional characterization of 18 CHUM HGSOC lines
- **Citation:** Fleury H, Carmona E, Morin VG, Meunier L, Masson JY, Tonin PN, Provencher D, Mes-Masson AM (2016). "Cumulative defects in DNA repair pathways drive the PARP inhibitor response in high-grade serous epithelial ovarian cancer cell lines." *Oncotarget* 7(26):40152–40168. DOI: 10.18632/oncotarget.10308 · PMID 27374179
- **Key finding:** Profiled **18 HGSOC lines** for gene expression, homologous-recombination functionality and olaparib sensitivity; PARPi response was driven by cumulative HR-pathway defects, not HR status alone. (This is the recurring "27374179" reference against many CHUM HGSC lines in the OvCAN catalogue.)
- **Relevance:** Medium — functional/drug-response annotation layer for the CHUM HGSOC lines in the resource.
- **Source type:** Primary

### Rare-subtype model authentication & establishment

### Karnezis et al. 2020/2021 — re-assigning COV434 and TOV-112D (spans Q1 + Q2)
- **Citation:** Karnezis AN, Wang Y, Keul J, Tessier-Cloutier B, Magrill J, Kommoss S, Senz J, Yang W, Proctor L, Schmidt D, Clement PB, Gilks CB, Huntsman DG, Kommoss F (2021). "Re-assigning the histologic identities of COV434 and TOV-112D ovarian cancer cell lines." *Gynecologic Oncology* 160(2):568–578. DOI: 10.1016/j.ygyno.2020.12.004 · PMID 33328126
- **Key finding:** Using the original tumour blocks, **COV434** (long labelled a granulosa-cell-tumour line) is re-diagnosed as **SCCOHT** (morphology, IHC, SMARCA4/genetics), and **TOV-112D** (long labelled grade-3 endometrioid) as **dedifferentiated ovarian carcinoma**. Demonstrates that even "known" rare-subtype lines have been mis-histotyped.
- **Relevance:** High — directly affects two lines the resource carries (COV434 as SCCOHT; TOV-112D listed as endometrioid in the OvCAN catalogue — a discrepancy worth reconciling in the descriptor).
- **Source type:** Primary

### Gamwell et al. 2013 — BIN-67 SCCOHT model
- **Citation:** Gamwell LF, Gambaro K, Merziotis M, Crane C, Arcand SL, Bourada V, Davis C, Squire JA, Huntsman DG, Tonin PN, Vanderhyden BC (2013). "Small cell ovarian carcinoma: genomic stability and responsiveness to therapeutics." *Orphanet Journal of Rare Diseases* 8:33. PMID 23433318 · PMC3635907 · https://pmc.ncbi.nlm.nih.gov/articles/PMC3635907/
- **Key finding:** Characterized the **BIN-67 SCCOHT line** — genomic stability (SNP array, spectral karyotyping, sequencing), xenografts recapitulating human SCCOHT, chemoresistance but sensitivity to oncolytic viruses (VSV, JX-594). One of only ~3 SCCOHT cell models in existence.
- **Relevance:** High — BIN-67 is in the resource (Vanderhyden lab, per OvCAN catalogue) and is a scarce rare-subtype model. (First-author name to re-confirm at proof stage; PMID/PMC verified.)
- **Source type:** Primary

### The resource itself — OvCAN Collection catalogue (grey literature, read directly)
- **Citation:** Tone A (Ovarian Cancer Canada) (2023). "Building your scientific toolbox: OCC's OvCAN Collection of ovarian cancer research models." Catalogue v4, 1-Dec-2023. https://ovariancanada.org/wp-content/uploads/2023/12/1-Dec-2023_OvCAN-Collection-v4.pdf
- **Key finding:** A **virtual collection** of OCC/OvCAN-funded high-fidelity OC models. Patient-derived **cell lines** are organized by histotype and lab: **HGSC** (Mes-Masson/Provencher CHUM TOV/OV lines — dozens); **LGSC** (BC Cancer, Carey/Lee lab — ~20 **VOA** lines, e.g., VOA-1312, VOA-1056, VOA-3723, VOA-4627, iOv241Ca; plus CHUM TOV81D); **clear cell** (TOV21G, TOV3392D; Huntsman VOA4841/6861/10816/12539, XVOA295/867); **endometrioid** (TOV112D; VOA4395, VOA5596); **mucinous** (TOV2414; VOA8762, VOA8771); **adenocarcinoma** (OV90, OV3331); **carcinosarcoma/MMMT** (Huntsman VOA5217, VOA5436); **SCCOHT** (Huntsman VOA12721, COV434; Vanderhyden BIN67). Data types listed per line include **Clariom S microarray, WES, SNP arrays (Affy/Illumina), RNA-Seq, and proteomics/RPPA/MS** — i.e., the multi-omic layers of the resource. Catalogue also covers organoids/PDX/syngeneic models (out of scope for this cluster).
- **Relevance:** High (context) — this is the pre-publication catalogue of the exact collection the Data Descriptor formalizes; establishes provenance (CHUM + BC Cancer VOA + Vanderhyden), histotype breadth, and existing omics. **No peer-reviewed primary descriptor of the integrated OvCAN multi-omic cell-line resource yet exists** — that is the novelty gap the descriptor fills.
- **Source type:** Grey literature / resource catalogue

---

### Q2 — Authentication norms & misidentification of OC cell lines

### ICLAC Register of Misidentified Cell Lines (authoritative register)
- **Citation:** International Cell Line Authentication Committee (ICLAC). "Register of Misidentified Cell Lines," version 14, released 15-Feb-2026 (originated by Capes-Davis A & Freshney RI, 2010). https://iclac.org/databases/cross-contaminations/ · Guide: https://iclac.org/wp-content/uploads/ICLAC_Guide-to-Human-Cell-Line-Authentication_02-Mar-2023.pdf
- **Key finding:** Curated register of cell lines known to be misidentified by cross-contamination or mislabelling. **v14 lists 608 cell lines**, of which **560 are misidentified with no known authentic stock** (Table 1). The Guide sets the community norm that any line not entering a lab via a public bank must be STR-authenticated. Example ovarian misidentifications include HeLa-derived lines.
- **Relevance:** High — the definitive reference list to cross-check every line in the resource against; a Data Descriptor reviewer will expect this check.
- **Source type:** Register / reference (primary curation)

### ANSI/ATCC ASN-0002 — STR profiling standard (methodological norm)
- **Citation:** ATCC Standards Development Organization Workgroup ASN-0002. "Authentication of Human Cell Lines: Standardization of Short Tandem Repeat (STR) Profiling." ANSI/ATCC ASN-0002-2011; **Revised ASN-0002-2022**. https://webstore.ansi.org/standards/atcc/ansiatccasn00022022 · NIST program: https://www.nist.gov/programs-projects/cell-line-authentication/cell-line-id-and-authentication-human-cell-lines
- **Key finding:** The consensus standard specifying STR-profiling methodology, data analysis, QC, match interpretation and a searchable public database for unambiguous human cell-line authentication; recommends a core panel of autosomal STR loci (expanded in the 2022 revision). Widely required by journals/funders.
- **Relevance:** High — the method the resource's authentication section should cite/conform to.
- **Source type:** Standard

### Korch et al. 2012 — DNA profiling of endometrial & ovarian lines (seed, confirmed)
- **Citation:** Korch C, Spillman MA, Jackson TA, Jacobsen BM, Murphy SK, Lessey BA, Jordan VC, Bradford AP (2012). "DNA profiling analysis of endometrial and ovarian cell lines reveals misidentification, redundancy and contamination." *Gynecologic Oncology* 127(1):241–248. DOI: 10.1016/j.ygyno.2012.06.017 · PMID 22710073
- **Key finding:** Profiled **51 ovarian cancer lines** by STR, TP53 SNPs and microsatellite instability. Found **10 redundant lines** and **5 lines (A2008, OV2008, C13, SK-OV-4, SK-OV-6) that are actually cervical (HeLa-lineage) cells**, plus widespread misidentification/loss of integrity across endometrial+ovarian lines. The core ovarian-specific authentication study.
- **Relevance:** High — names specific ovarian lines to exclude/flag; direct precedent for an ovarian-line authentication audit.
- **Source type:** Primary

### Stordal et al. 2024 — CNV + mutation needed to authenticate OC lines
- **Citation:** Stordal B, Farrelly AM, Hennessy BT (2024). "Chromosomal copy number and mutational status are required to authenticate ovarian cancer cell lines as appropriate cell models." *Molecular Biology Reports* 51(1):784. DOI: 10.1007/s11033-024-09747-4 · PMC11213756
- **Key finding:** Argues STR identity alone is insufficient for OC lines: **copy-number profile and mutational status** (e.g., near-universal TP53 in HGSOC) are also required to confirm a line is an appropriate model of its assigned histotype. Extends authentication from "is it the right line" to "is it the right disease model."
- **Relevance:** High — motivates using the resource's own WES/CNV to authenticate histotype assignment, not just STR.
- **Source type:** Primary (analysis/review hybrid)

### Huang et al. 2017 — cross-contamination survey of 278 lines (general-rate context)
- **Citation:** Huang Y, Liu Y, Zheng C, Shen C (2017). "Investigation of Cross-Contamination and Misidentification of 278 Widely Used Tumor Cell Lines." *PLoS ONE* 12(1):e0170384. DOI: 10.1371/journal.pone.0170384 · PMC5249119
- **Key finding:** STR-profiled 278 lines from 28 institutes; **46.0% (128/278) were cross-contaminated or misidentified**, many HeLa-derived. Establishes the alarming base rate of misidentification in routine tumour-line stocks (includes ovarian lines among those with MSI, though not itemized).
- **Relevance:** Medium — quantifies why authentication is non-optional; general (not ovarian-specific).
- **Source type:** Primary

### Dirks/Yu-type large cross-contamination survey (482 lines) — supporting context
- **Citation:** Yu M, Selvaraj SK, Liang-Chu MMY, et al. (2017). "A Combination of Species Identification and STR Profiling Identifies Cross-contaminated Cells from 482 Human Tumor Cell Lines." *Scientific Reports* 7:9774. DOI: 10.1038/s41598-017-09660-w · PMID 28851942 · PMC5575032
- **Key finding:** Large-scale STR + species-ID survey of 482 human tumour cell lines identifying cross-contaminated/misidentified stocks; reinforces the field-wide prevalence documented by Huang 2017. (Author list to confirm at proof stage; identifiers verified.)
- **Relevance:** Medium — corroborating prevalence data.
- **Source type:** Primary

### Korch et al. 2026 — practical STR genotype-analysis guide (current norm)
- **Citation:** Korch C, et al. (2026). "A Biomedical Researcher's Guide for Analyzing Short Tandem Repeat (STR) Genotypes of Human Cell Lines and in Vitro Tissue Samples Using Three Standard Authentication Algorithms." *Current Protocols*. DOI: 10.1002/cpz1.70387
- **Key finding:** Step-by-step protocol for interpreting STR genotypes against the three standard match algorithms — the how-to companion to ASN-0002 for a methods/authentication section.
- **Relevance:** Medium — practical citation for the descriptor's authentication methods.
- **Source type:** Protocol

### SKOV3 / A2780 as poor HGSOC models — covered above by Domcke 2013 and Anglesio 2013
- Both seed papers (reported under Q1) are the primary evidence that **SKOV3 and A2780 are poor/likely-mis-histotyped HGSOC models** (low genomic suitability; classified non-serous). Beaufort 2014 and Haley 2016 similarly place OVCAR3/CAOV3/Kuramochi/OVSAHO — not SKOV3/A2780 — as HGSOC-representative. No separate citation added to avoid duplication.

---

## Search Log
- Query 1: "Ince 2015 Nature Communications 25 ovarian tumour cell lines phenocopy primary tumours" — 10 examined, 3 relevant
- Query 2: "Domcke 2013 Nature Communications evaluating cell lines as tumour models genomic profiles ovarian" — 10 examined, 4 relevant
- Query 3: "Beaufort 2014 PLoS One Ovarian Cancer Cell Line Panel OCCP characterization" — 9 examined, 3 relevant
- Query 4: "Anglesio 2013 PLoS One type-specific ovarian cancer cell line models clear cell" — 9 examined, 2 relevant
- Query 5: "Korch 2012 Gynecologic Oncology ovarian cancer cell line authentication STR profiling misidentification" — 8 examined, 3 relevant
- Query 6: "Provencher 2000 In Vitro Cellular Developmental Biology TOV OV ... CHUM" — 8 examined, 2 relevant
- Query 7: "Létourneau 2012 BMC Cancer novel ... ovarian cancer cell lines TOV" — 10 examined, 2 relevant
- Query 8: "Fleury Tonin transcriptomic profiling ovarian cancer cell lines TOV" — 10 examined, 3 relevant
- Query 9: "ICLAC register misidentified cross-contaminated cell lines ovarian STR authentication" — 9 examined, 4 relevant
- Query 10: "SKOV3 A2780 OVCAR3 not representative high-grade serous ovarian carcinoma model" — 8 examined, 4 relevant
- Query 11: "ANSI ATCC ASN-0002 standard STR profiling authentication human cell lines consensus" — 10 examined, 4 relevant
- Query 12: "Kuramochi OVSAHO OVCAR4 OVCAR8 most suitable HGSOC cell line models Domcke Beaufort" — 9 examined, 4 relevant
- Query 13: "SCCOHT cell lines BIN67 COV434 SMARCA4 ... hypercalcemic type model" — 7 examined, 3 relevant
- Query 14: "mucinous ovarian cancer cell lines misidentified colorectal gastric contamination authentic" — 8 examined, 3 relevant
- Query 15: "Coscia 2016 Nature Communications integrative proteomic profiling ovarian cancer cell lines precursor" — 9 examined, 2 relevant
- Query 16: "BC Cancer Agency VOA ovarian clear cell endometrioid cell lines Anglesio Huntsman" — 9 examined, 3 relevant
- Query 17: "ovarian carcinosarcoma cell line malignant mixed Mullerian tumor model established" — 8 examined, 0 cell-line-specific (reviews only)
- Query 18: "comprehensively characterized cell line panel highly representative ... HGSC Oncotarget authors" — 5 examined, 2 relevant
- Query 19: "Fleury 2015 six new ... ovarian cancer cell lines ... DOI" — 8 examined, 2 relevant
- Query 20: "Re-assigning the histologic identities of COV434 and TOV-112D ..." — 6 examined, 2 relevant
- Query 21: "OvCAN collection Ovarian Cancer Canada VOA TOV OV cell lines multi-omic" — 7 examined, 0 (no primary pub; catalogue found via Q16)
- Query 22: "Korch DNA profiling endometrial ovarian ... Gynecologic Oncology 2012 DOI" — 7 examined, 2 relevant
- Query 23: "Chromosomal copy number and mutational status ... authenticate ovarian cancer cell lines" — 5 examined, 2 relevant
- Query 24: "Distinct transcriptional programs stratify ovarian cancer cell lines five histological subtypes" — 7 examined, 2 relevant
- Query 25: "Modeling the Diversity of EOC through Ten Novel ... Cell Lines authors DOI" — 5 examined, 1 relevant
- Query 26: "Ouellet Characterization of three new serous ... BMC Cancer 2008 DOI" — 8 examined, 2 relevant
- Query 27: "investigation cross-contamination misidentification 278 widely used tumor cell lines STR" — 6 examined, 2 relevant
- Query 28: "Ince 2015 ovarian cell lines OCI panel histological subtypes ..." — 8 examined, 1 relevant
- Query 29: "PMID 27374179 ovarian cancer cell lines proteomic ..." — 2 rounds, identified as Fleury 2016 Oncotarget
- Query 30: "Functional characterization of a panel of HGSOC cell lines ... Oncotarget authors" — 6 examined, 1 relevant
- Query 31: "BIN-67 SCCOHT ... 2013 Gamwell Vanderhyden PMID 23433318" — 8 examined, 2 relevant
- WebFetch A: PMC4473807 (Ince 2015) — subtype composition + OCMI medium + efficiency
- WebFetch B: journals.plos.org 0170384 (Huang 2017) — authors + rate
- WebFetch C + Read: OvCAN Collection catalogue PDF (20 pp; read pp.1–16) — full line/subtype/omics inventory

## Gaps
- **Ovarian carcinosarcoma (MMMT) established 2D cell lines — literature gap (genuine).** Repeated searches returned only clinical/pathology reviews and a patient-derived *organoid* model; no well-known, authenticated ovarian-carcinosarcoma *cell line* is described in the primary literature. The OvCAN catalogue lists two unpublished ones (VOA5217, VOA5436, Huntsman lab). This is the resource's clearest rare-subtype novelty and a real hole in the field.
- **A dedicated primary publication describing the BC Cancer "VOA" cell-line collection — search-limitation likely + partial literature-gap.** VOA lines appear across papers (Anglesio 2013; LGSC papers cited in the catalogue: PMIDs 27822414, 28545687, 30636931, 26076164) but I did not confirm a single "derivation/characterization of the VOA panel" descriptor. Worth a targeted follow-up on those LGSC PMIDs (esp. the Carey/Lee low-grade serous line series).
- **Exact per-histotype counts in Ince 2015 (OCI) — search-limitation.** The paper confirms the subtypes present (incl. carcinosarcoma, dysgerminoma) but the retrieved text did not tabulate n per subtype; would need full-text supplementary tables.
- **No peer-reviewed integrated multi-omic descriptor of the OvCAN cell-line resource — literature gap (this is the novelty).** Only the OCC grey-literature catalogue exists.
- **Mucinous-line authenticity specifics — partially covered.** Cheasley et al. 2019 *Nat Commun* "The molecular origin and taxonomy of mucinous ovarian carcinoma" (DOI 10.1038/s41467-019-11862-x) and Meagher et al. 2025 *J Pathol* "Cellular origins of mucinous ovarian carcinoma" (DOI 10.1002/path.6407) bear on distinguishing authentic mucinous OC from GI contaminants, but their cell-line-authentication content was not deeply verified — flagged as leads for the rare-subtype-biology cluster rather than reported as findings here.

## Breadth Flag
- **Individual rare-subtype line-establishment papers** (single clear-cell, single LGSC lines — e.g., FDOV1, CAISMOV24, OVPA8, iOvCa/ iOv241Ca origin paper) were surfaced but not individually verified; a deeper pass could enumerate every rare-subtype line and its origin paper.
- **BC Cancer LGSC (Carey/Lee) line series (PMIDs 27822414, 28545687, 30636931)** deserves its own focused search — likely the primary source(s) for the ~20 VOA LGSC lines in the resource.
- **The two large cross-contamination surveys** (Huang 2017; the 482-line Sci Rep paper) — author lists should be re-confirmed at proof stage; identifiers are verified.
- Overlap note: SKOV3/A2780 "poor HGSOC model" evidence and the COV434/TOV-112D reassignment sit at the Q1/Q2 boundary and are also relevant to the molecular-subtyping and rare-subtype-biology clusters.
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
# Lit Review — Cluster 03: Ovarian Cancer Molecular Subtyping & Its Validity on Cell Lines

**Agent:** lit-subtyping
**Date:** 2026-07-23
**Questions:** Q5 (how HGSC transcriptomic subtypes are defined / classifiers built & validated) and Q6 (the critique that subtypes are driven by the tumour microenvironment, not tumour-cell-intrinsic state, and implications for applying classifiers to pure tumour-cell cultures / cell lines / organoids).

All claims below are grounded in WebSearch/WebFetch performed this session. Where a fast-model fetch produced a number, I note the source page. Preprints are flagged. Report is per-paper; no synthesis.

---

## Findings

### Q5 — The classifier / subtyping papers

### Tothill et al. 2008 — first molecular subtypes (C1–C6), stromal C1 = worst prognosis
- **Citation:** Tothill RW, Tinker AV, George J, Brown R, Fox SB, Lade S, Johnson DS, Trivett MK, Etemadmoghadam D, Locandro B, Traficante N, Fereday S, Hung JA, Chiew Y-E, Haviv I, Australian Ovarian Cancer Study Group, Gertig D, DeFazio A, Bowtell DDL (2008). "Novel Molecular Subtypes of Serous and Endometrioid Ovarian Cancer Linked to Clinical Outcome." *Clinical Cancer Research* 14(16):5198–5208. PMID 18698038. DOI 10.1158/1078-0432.CCR-08-0196. https://aacrjournals.org/clincancerres/article/14/16/5198/72783
- **Key finding:** K-means clustering of microarray profiles from 285 serous/endometrioid ovarian, peritoneal and fallopian-tube tumours defined six subtypes (C1–C6). Four (C1, C2, C4, C5) captured high-grade serous disease. Critically, the subtypes are defined in large part by non-tumour compartments: **C1 = high desmoplasia / reactive stroma; C2 and C4 = high intratumoral CD3+ T-cell content; C5 = mesenchymal.** The **C1 (reactive-stroma/desmoplastic) subtype had the worst overall survival**, and a later search source noted ~40% of C1 tumours had >50% stromal component. This is the foundational observation that "subtype" tracks stromal/immune admixture.
- **Relevance:** High — the historical origin of ovarian subtyping; already shows subtypes are keyed to stroma (C1) and T cells (C2/C4), directly seeding Q6.
- **Source type:** Primary

### TCGA 2011 — the canonical four subtypes (immunoreactive / differentiated / proliferative / mesenchymal); NOT prognostic
- **Citation:** Cancer Genome Atlas Research Network (2011). "Integrated genomic analyses of ovarian carcinoma." *Nature* 474(7353):609–615. DOI 10.1038/nature10166. https://www.nature.com/articles/nature10166
- **Key finding:** NMF consensus clustering on ~1,500 variable genes across 489 high-grade serous tumours yielded four transcriptional subtypes, each named for its dominant expression program: **Immunoreactive** = T-cell chemokines CXCL11/CXCL10 and receptor CXCR3; **Proliferative** = HMGA2, SOX11, plus proliferation markers MCM2/PCNA (and low differentiation markers); **Differentiated** = MUC16/MUC1 and secretory fallopian-tube marker SLPI; **Mesenchymal** = HOX genes and "markers suggestive of increased stromal components" including myofibroblast and pericyte markers. Per the WebFetch of the paper, **"survival duration did not differ significantly for transcriptional subtypes"** — i.e., the four subtypes were NOT significantly prognostic. Prognostic stratification instead came from a **separate 193-gene survival signature** (108 poor- / 85 favourable-prognosis genes) trained on 215 samples and validated in 4 independent datasets. TP53 mutated in ~96% of tumours.
- **Relevance:** High — the reference classifier and subtype nomenclature. Note that two of the four subtype names (immunoreactive, mesenchymal) are, by TCGA's own gene descriptions, immune- and stroma-defined.
- **Source type:** Primary

### Verhaak et al. 2013 — CLOVAR classifier + prognostic signature
- **Citation:** Verhaak RGW, Tamayo P, Yang J-Y, Hubbard D, Zhang H, Creighton CJ, Fereday S, Lawrence M, Carter SL, Mermel CH, Kostic AD, Etemadmoghadam D, Saksena G, Cibulskis K, Duraisamy S, Levanon K, Sougnez C, Tsherniak A, Gomez S, Onofrio R, Gabriel S, Chin L, Zhang N, Spellman PT, Zhang Y, Akbani R, Hoadley KA, Kahn A, Köbel M, Huntsman D, Soslow RA, Defazio A, Birrer MJ, Gray JW, Weinstein JN, Bowtell DD, Drapkin R, Mesirov JP, Getz G, Levine DA, Meyerson M (2013). "Prognostically relevant gene signatures of high-grade serous ovarian carcinoma." *Journal of Clinical Investigation* 123(1):517–525. DOI 10.1172/JCI65833. https://www.jci.org/articles/view/65833
- **Key finding:** Integrated the TCGA subtype and survival classifiers into a framework called **CLOVAR (Classification of Ovarian Cancer)**. Uses the same four subtype names (differentiated, immunoreactive, mesenchymal, proliferative). A supervised **193-gene prognostic signature** was built from 215 TCGA profiles and, together with the subtypes, **validated on an independent dataset of 879 HGS-OvCa expression profiles.** The worst-outcome group (**23% of cases**) had **median survival 23 months and 63% platinum-resistance**, versus **46 months and 23%** in others; combining with BRCA1/2 status, residual disease and stage improved stratification. Per WebFetch, the paper notes some tumours are "characterized by high numbers of infiltrating T lymphocytes or stromal cells" (immunoreactive = T-cell infiltration) but did NOT, in the extracted text, explicitly equate the mesenchymal subtype with Tothill C1 stroma.
- **Relevance:** High — the classifier + prognostic signature most often applied downstream; the standard bulk-trained tool. n and HR-equivalent effect sizes captured.
- **Source type:** Primary

### Chen et al. 2018 — consensusOV (the classifier the prior analysis used)
- **Citation:** Chen GM, Kannan L, Geistlinger L, Kofia V, Safikhani Z, Gendoo DMA, Parmigiani G, Birrer M, Haibe-Kains B, Waldron L (2018). "Consensus on Molecular Subtypes of High-Grade Serous Ovarian Carcinoma." *Clinical Cancer Research* 24(20):5037–5047. DOI 10.1158/1078-0432.CCR-18-0784. PMID 30084834. https://pmc.ncbi.nlm.nih.gov/articles/PMC6207081/
- **Key finding:** Reconciled prior HGSOC classifiers (**Helland/Tothill 2008–2011, TCGA/Verhaak 2011–2013, Konecny 2014** — all of which nominally produce the same four subtypes). Built **consensusOV**, a **random-forest classifier trained on concordantly-subtyped tumours across 15 datasets (1,770 HGSOC tumours)**, using **100 genes from Verhaak et al.** encoded as **binary gene-pair (rank) relationships** so it applies across expression platforms without cross-dataset normalization; evaluated by leave-one-dataset-out cross-validation. Introduced a **"margin" score** (difference between top and second subtype scores) to flag robustly-classifiable vs intermediate tumours. Per WebFetch, the subtypes retain their microenvironment definitions: **Immunoreactive = elevated lymphocyte (p<0.05) and neutrophil (p<0.10) infiltration; Mesenchymal = desmoplasia associated with infiltrating stromal cells.**
- **Relevance:** High — this is the exact classifier applied to the 15 HGSC lines. Note the classifier's own paper defines two of four subtypes by immune/stromal infiltrate — compartments absent from pure cell-line cultures.
- **Source type:** Primary

---

### Q6 — The tumour-microenvironment critique + applicability to cell lines / organoids

### Zhang, Wang & Cliby 2019 — mesenchymal markers are expressed by stroma, not cancer cells (IHC)
- **Citation:** Zhang Q, Wang C, Cliby WA (2019). "Cancer-associated stroma significantly contributes to the mesenchymal subtype signature of serous ovarian cancer." *Gynecologic Oncology* 152(2):368–374. DOI 10.1016/j.ygyno.2018.11.014. PMID 30448260. https://pubmed.ncbi.nlm.nih.gov/30448260/
- **Key finding:** IHC for eight mesenchymal-subtype signature genes (ACTA2, COL5A1, COL11A1, FAP, POSTN, VCAN, ZEB1, p-SMAD2) in 15 mesenchymal-subtype HGSOC cases. **Four proteins (COL5A1, VCAN, FAP, ZEB1) were almost exclusively expressed by stroma, not cancer cells**, and stromal expression dominated for the other four (ACTA2, COL11A1, POSTN, p-SMAD2). Conclusion: "the existing molecular classification reflects a significant stromal contribution" — the mesenchymal subtype signature is substantially cancer-associated-fibroblast/stroma-derived rather than tumour-cell-intrinsic. (This is the "Zhang et al. 2019" from the seed list.)
- **Relevance:** High — direct protein-level evidence that mesenchymal-subtype genes live in the stromal compartment; a pure tumour-cell culture would lack the cells that generate this signature.
- **Source type:** Primary

### Schwede et al. 2020 — stromal admixture reclassifies subtypes and abolishes prognostic signatures (THE key critique)
- **Citation:** Schwede M, Waldron L, Mok SC, Wei W, Basunia A, Merritt MA, Mitsiades CS, Parmigiani G, Harrington DP, Quackenbush J, Birrer MJ, Culhane AC (2020; online 2019). "The Impact of Stroma Admixture on Molecular Subtypes and Prognostic Gene Signatures in Serous Ovarian Cancer." *Cancer Epidemiology, Biomarkers & Prevention* 29(2):509–519. DOI 10.1158/1055-9965.EPI-18-1359. PMID 31871106. https://pmc.ncbi.nlm.nih.gov/articles/PMC7448721/ (bioRxiv preprint: 10.1101/496406)
- **Key finding:** Built a **688-gene tumour–stroma signature (461 stroma-overexpressed / 227 epithelial-overexpressed)** from paired microdissected epithelium and stroma. In computational mixing simulations (n=38 microdissected pairs), **small stromal increases destabilize subtype calls: at 10% added stroma 6/38 tumours changed subtype (e.g., C4→C1/C2); at 30% stroma ~15/38 (>one-third) reclassified. The TCGA classifier was similarly vulnerable (5/38 reclassified at 20%, 8/38 at 30%).** Furthermore, **five published prognostic gene signatures lost independent prognostic value after adjusting for pathologist-scored stromal content** in multivariate models. Authors conclude **"single-cell analyses may be required to refine the molecular subtypes of HGSOC,"** because bulk classifications conflate tumour-intrinsic biology with stromal admixture.
- **Relevance:** High — the single most quantitative demonstration that ovarian "molecular subtype" is a function of stromal fraction. Strong basis for caveating/dropping subtype calls on samples whose stromal fraction differs radically from the training tumours (e.g., pure cell lines = ~0% stroma). NOTE: the seed's "Schwede ~2013" appears to conflate this paper with a *different* Schwede paper (see next entry).
- **Source type:** Primary

### Schwede et al. 2013 — (seed disambiguation) different paper: stem-cell-like signature / Type II subtype
- **Citation:** Schwede M, Spentzos D, Bentink S, Hofmann O, Haibe-Kains B, Harrington D, Quackenbush J, Culhane AC (2013). "Stem Cell-Like Gene Expression in Ovarian Cancer Predicts Type II Subtype and Prognosis." *PLoS ONE* 8(3):e57799. DOI 10.1371/journal.pone.0057799. https://pmc.ncbi.nlm.nih.gov/articles/PMC3594231/
- **Key finding:** Applied embryonic/adult/cancer stem-cell signatures to 145 serous ovarian tumours; a stem-cell-like signature was prognostic within high-stage serous disease and helped define Type I/II subtypes. This is NOT the stroma-admixture critique — it is a separate, earlier Schwede paper. Flagged only to resolve the seed's "Schwede ~2013 and later 2020" (the 2013 item is this stem-cell paper; the stroma-admixture critique is the 2020 CEBP paper above).
- **Relevance:** Medium (disambiguation) — prevents mis-citation; not itself a TME-confounding critique.
- **Source type:** Primary

### Olbrecht et al. 2021 — scRNA-seq: subtype label is set by which stromal/immune cells are present
- **Citation:** Olbrecht S, Busschaert P, Qian J, Vanderstichele A, Loverix L, Van Gorp T, Van Nieuwenhuysen E, Han S, Van den Broeck A, Coosemans A, Van Rompuy A-S, Lambrechts D, Vergote I (2021). "High-grade serous tubo-ovarian cancer refined with single-cell RNA sequencing: specific cell subtypes influence survival and determine molecular subtype classification." *Genome Medicine* 13:111. DOI 10.1186/s13073-021-00922-x. PMID 34238352. https://pmc.ncbi.nlm.nih.gov/articles/PMC8268616/
- **Key finding:** scRNA-seq of 18,403 cells from 7 treatment-naïve HGSTOC patients defined 43 phenotypes (11 cancer, 32 stromal). Mapping TCGA subtype signatures onto single cells showed the subtypes are compartment-specific: **fibroblast subclusters score highest for the mesenchymal subtype; immune cells for immunoreactive; myofibroblasts (FB_MYH11) and mesothelial cells (FB_CALB2) for the differentiated subtype.** They demonstrate "an important contribution of genes expressed by fibroblasts to identify the mesenchymal subtype" — i.e., **subtype assignment is driven by stromal/immune cell composition rather than malignant-cell state.** Relative frequencies of myofibroblasts, TGF-β CAFs, mesothelial and lymphatic endothelial cells predicted poor outcome; plasma cells favourable.
- **Relevance:** High — single-cell proof that three of the four subtype signatures are emitted by non-malignant cells. Implies a pure tumour-cell culture cannot legitimately be "mesenchymal" or "immunoreactive."
- **Source type:** Primary

### Hippen et al. 2023 (preprint) — deconvolution: ConsensusOV subtype probabilities track cell-type proportions
- **Citation:** Hippen AA, Davidson NR, Barnard ME, Weber LM, Gertz J, Doherty JA, Hicks SC, Greene CS (2023). "Deconvolution reveals compositional differences in high-grade serous ovarian cancer subtypes." *bioRxiv* preprint. DOI 10.1101/2023.06.14.544991. https://www.biorxiv.org/content/10.1101/2023.06.14.544991v1
- **Key finding:** Deconvolved bulk HGSOC transcriptomes using single-cell references (164 tumour reference samples across two datasets) and compared inferred cell-type proportions to TCGA subtypes. Reported that **ConsensusOV subtype probabilities reflect tumour cellular composition**, i.e., subtype differences correspond to differences in stromal/immune fractions rather than distinct malignant programs.
- **Relevance:** High — population-scale deconvolution linking the exact consensusOV tool to composition. (Preprint — treat as suggestive/not peer-reviewed; note it is the predecessor work to Tanis et al. 2026 from the same group.)
- **Source type:** Primary (preprint)

### Tanis et al. 2026 (preprint) — most direct statement: subtypes are composition, not malignant state
- **Citation:** Tanis S, Lixandrão M, Ivich A, Grieshober L, Lawson-Michod KA, Collin LJ, Peres LC, Salas LA, Marks JR, Bitler BG, Greene CS, Schildkraut JM, Doherty JA, Davidson NR (2026). "Transcriptomic subtypes in high-grade serous ovarian cancer are driven by tumor cellular composition." *bioRxiv* preprint, posted 2026-04-21. https://www.biorxiv.org/content/10.64898/2026.04.16.719000v1 (PMID listed as 42079185; DOI as shown 10.64898/2026.04.16.719000 — unusual prefix, confirm final identifier before citing in manuscript).
- **Key finding:** Integrated single-cell-derived pseudobulk simulations with deconvolution of **1,834 primary HGSC tumours** (RNA-seq + microarray). **Cellular composition alone predicted subtype labels with ROC-AUC 0.81–0.95; the mesenchymal subtype showed the strongest composition-driven signal.** A secondary, composition-independent expression signal existed but did not define the dominant subtype structure. Conclusion: HGSC transcriptomic subtypes are "**features of the tumor ecosystem rather than discrete malignant states.**"
- **Relevance:** High — the strongest, most recent, most quantitative statement of the exact thesis behind Q6. Directly supports the position that consensusOV/TCGA subtypes are not tumour-cell-intrinsic. (Preprint — flag as not yet peer-reviewed; verify DOI.)
- **Source type:** Primary (preprint)

### Supporting context (from search snippets; secondary)
- **TCGA subtype prognosis is contested:** search sources (e.g., a systematic review/meta-analysis "Better or worse? The prognostic role of the mesenchymal subtype…", PMC9582683, and a review PMC10740166) indicate the original TCGA subtypes were not significantly prognostic, while later re-analyses report proliferative and mesenchymal as worst OS. Reported here as context; full citations not individually verified this session (see Gaps).

---

## Search Log
- Query 1: "Tothill 2008 novel molecular subtypes serous ovarian cancer Clinical Cancer Research" — 8 examined, 1 primary relevant
- Query 2: "TCGA 2011 Integrated genomic analyses ovarian carcinoma four transcriptional subtypes…" — 8 examined, 1 primary relevant
- Query 3: "Verhaak 2013 prognostic gene expression signature HGSC CLOVAR classifier" — 8 examined, 1 primary relevant
- Query 4: "Chen 2018 Consensus Molecular Subtypes HGSOC consensusOV Clinical Cancer Research" — 7 examined, 1 primary relevant
- Query 5: "Schwede stromal tumor content confounds ovarian cancer molecular subtype prognosis" — 7 examined, 1 primary relevant (2020 CEBP)
- Query 6: "Zhang 2019 tumor epithelial vs stroma contribution ovarian subtype single-cell" — 8 examined, 2 relevant (Zhang/Cliff Gynecol Oncol; Olbrecht)
- Query 7: "single-cell RNA-seq HGSC mesenchymal subtype stromal fibroblast immunoreactive immune deconvolution" — 8 examined, 3 relevant (Olbrecht, Hippen, Tanis)
- Query 8: "TCGA molecular subtype ovarian cancer cell lines classification applicability caveat pure tumor cells" — 7 examined, 1 relevant + cell-line-histotype leads
- Query 9: "Schwede 2013 microdissection ovarian cancer stroma epithelium…" — 8 examined, disambiguation
- Query 10: "mesenchymal subtype ovarian CAF desmoplasia stromal reactive not tumor intrinsic" — 10 examined, 1 relevant (Zhang)
- Query 11: "Schwede 2013 PLoS One ovarian cancer stem cell gene expression subtype prognosis" — 7 examined, 1 (disambiguation)
- Query 12: "TCGA ovarian subtypes not significantly associated OS; Tothill C1 stromal worst; Verhaak mesenchymal" — 9 examined, context
- Query 13: consensusOV/TCGA subtype applied to cell lines/organoids — 10 examined, gap-confirming
- Query 14: exact-title search for Zhang/Cliby Gynecol Oncol paper — authorship confirmed
- Query 15: exact-title search for Hippen deconvolution preprint — authorship confirmed
- Query 16: Tanis/Davidson 2026 verification — venue/authors confirmed
- WebFetch: Chen 2018 (PMC6207081); Verhaak 2013 (jci.org/65833); Schwede 2020 (PMC7448721); Olbrecht 2021 (PMC8268616); Zhang 2019 (PubMed 30448260); TCGA 2011 (PMC3163504); Tanis 2026 (PubMed page). AACR (aacrjournals) and ScienceDirect returned 403; used PMC/PubMed mirrors instead.

## Gaps
1. **No paper directly applying consensusOV/TCGA molecular subtypes to a panel of HGSC CELL LINES with a methodological critique.** The critique literature is all on bulk *tumours*. The inference "cell lines lack the stroma/immune cells that define mesenchymal & immunoreactive subtypes, so those calls are uninterpretable" is strongly supported by Zhang 2019, Schwede 2020, Olbrecht 2021, Hippen 2023, and Tanis 2026 — but is not, to my search, directly demonstrated in a cell-line study. This is an important gap for the Data Descriptor to state explicitly.
2. **Cell-line subtyping papers found are HISTOTYPE, not TCGA molecular subtype:** e.g., "Classification of ovarian cancer cell lines using transcriptional profiles defines the five major pathological subtypes" (bioRxiv 2020) and "Subtype Characterization of Ovarian Cancer Cell Lines Using Machine Learning" (Cancers 2025, doi 10.3390/cancers17213509) classify HGSC vs clear-cell/endometrioid/etc., not immunoreactive/mesenchymal/etc. Not verified in depth — flagged as leads for the cell-line-panel cluster.
3. **Organoids:** no organoid-specific paper on TCGA/consensusOV subtype applicability was surfaced. Gap.
4. **TCGA-subtype prognosis meta-analysis** (PMC9582683) and review (PMC10740166) noted only from snippets; full citations not individually verified.
5. **Tothill 2008 DOI** (10.1158/1078-0432.CCR-08-0196) follows AACR convention and matches the verified PMID/URL, but the DOI string itself was not independently fetched. **Tanis 2026 DOI prefix (10.64898)** is unusual for bioRxiv and should be re-verified before manuscript citation.

## Breadth Flag
Coverage of Q5 (four landmark classifier papers) and Q6 (the TME-confounding critique) is strong and multi-source, spanning IHC (Zhang), microdissection+simulation (Schwede), single-cell (Olbrecht), and population-scale deconvolution (Hippen, Tanis) — a consistent, convergent body of evidence from 2008→2026. The main breadth gap is the absence of a *direct* cell-line/organoid subtype-application study; the parent team should treat the cell-line applicability conclusion as an evidence-based inference, not a directly-cited result. Adjacent cell-line *histotype* classifiers overlap the "cell-line panels/authentication" cluster and were not pursued in depth here.
# Lit Review — Cluster 04: Rare ovarian cancer subtype biology, defining markers, and SWI/SNF convergence

**Scope:** Q7 (SCCOHT genetics), Q8 (defining markers of clear cell / mucinous / low-grade serous / endometrioid), Q9 (shared SWI/SNF disruption and the EC–SCCOHT convergence question).
**Rule applied:** every claim below is grounded in a WebSearch or WebFetch performed this session (see Search Log). Papers reported individually; no cross-paper synthesis.

---

## Findings

### Q7 — SCCOHT: the defining genetics (SMARCA4)

#### Ramos et al. 2014 — SMARCA4 germline + somatic mutations (one of three simultaneous 2014 discovery papers)
- **Citation:** Ramos P, Karnezis AN, Craig DW, et al. (2014). "Small cell carcinoma of the ovary, hypercalcemic type, displays frequent inactivating germline and somatic mutations in SMARCA4." *Nature Genetics* 46:427–429. https://www.nature.com/articles/ng.2928 (Corrigendum: https://www.nature.com/articles/ng0714-759a)
- **Key finding:** Inactivating germline and somatic mutations in the SWI/SNF chromatin-remodeling gene **SMARCA4** in **69% (9/13)** of SCCOHT cases; SMARCA4 (BRG1) protein loss by IHC in **82% (14/17)** of SCCOHT tumors versus only **0.4% (2/485)** of other primary ovarian tumors. Implicates SMARCA4 as the crucial oncogenic driver of SCCOHT.
- **Relevance:** High — one of the two SCCOHT lines in the resource; establishes the defining lesion and the diagnostic IHC contrast.
- **Source type:** Primary

#### Witkowski et al. 2014 — germline SMARCA4 in familial SCCOHT
- **Citation:** Witkowski L, Carrot-Zhang J, Albrecht S, et al. (2014). "Germline and somatic SMARCA4 mutations characterize small cell carcinoma of the ovary, hypercalcemic type." *Nature Genetics* 46:438–443. https://www.nature.com/articles/ng.2931 (PubMed: https://pubmed.ncbi.nlm.nih.gov/24658002/)
- **Key finding:** Exome sequencing of six individuals from three SCCOHT families found segregating deleterious **germline SMARCA4 mutations** in all three families (confirmed in a fourth family). Across cases with available DNA, **at least one germline or somatic deleterious SMARCA4 mutation in 30 of 32 cases**. Establishes SCCOHT as a hereditary as well as sporadic SMARCA4-driven disease.
- **Relevance:** High — germline dimension and near-universal SMARCA4 involvement.
- **Source type:** Primary

#### Jelinic et al. 2014 — recurrent biallelic somatic SMARCA4
- **Citation:** Jelinic P, Mueller JJ, Olvera N, et al. (2014). "Recurrent SMARCA4 mutations in small cell carcinoma of the ovary." *Nature Genetics* 46:424–426. https://www.nature.com/articles/ng.2922
- **Key finding:** **Biallelic inactivating somatic SMARCA4 mutations in 100% (12/12)** of SCCOHT samples analyzed. Independent confirmation that SMARCA4 loss is the recurrent, defining event.
- **Relevance:** High — third simultaneous confirmation; biallelic (two-hit) inactivation.
- **Source type:** Primary

#### Karnezis et al. 2016 — dual SMARCA4/SMARCA2 loss is the diagnostic IHC signature
- **Citation:** Karnezis AN, Wang Y, Ramos P, et al. (2016; online 2015). "Dual loss of the SWI/SNF complex ATPases SMARCA4/BRG1 and SMARCA2/BRM is highly sensitive and specific for small cell carcinoma of the ovary, hypercalcaemic type." *Journal of Pathology* 238(3):389–400. doi:10.1002/path.4633 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC4832362/)
- **Key finding:** SMARCA2 (BRM) is **epigenetically silenced** in SCCOHT, so tumors lose BOTH mutually exclusive SWI/SNF ATPases. **All 45 SCCOHT cases** showed dual SMARCA4/SMARCA2 loss by IHC (43 SMARCA4-deficient + 2 SMARCB1-deficient). Among **2,324 ovarian tumors**, the dual-loss phenotype was seen ONLY in SCCOHT (clear cell carcinomas, 15/360, lose SMARCA4 but retain SMARCA2). SMARCA4-loss-alone sensitivity ~91% (46/50), specificity ~99%; the **dual-loss pattern is completely specific for SCCOHT among ovarian tumors**.
- **Relevance:** High — the practical diagnostic marker; explains why SMARCA2 (as well as SMARCA4) is absent in SCCOHT models.
- **Source type:** Primary

#### Karnezis et al. 2021 — COV434 and TOV-112D re-identification (the reclassification seed, CONFIRMED)
- **Citation:** Karnezis AN, Chen SY, Chow C, Yang W, Hendricks WPD, Ramos P, Briones N, Mes-Masson AM, Bosse T, Gilks CB, Trent JM, Weissman B, Huntsman DG, Wang Y (2021). "Re-assigning the histologic identities of COV434 and TOV-112D ovarian cancer cell lines." *Gynecologic Oncology* 160(2):568–578. doi:10.1016/j.ygyno.2020.12.004 (PubMed: https://pubmed.ncbi.nlm.nih.gov/33328126/)
- **Key finding:** **COV434 — historically labeled a granulosa cell tumor — is reclassified as SCCOHT**: complete loss of SMARCA4 protein, focal cytoplasmic SMARCA2, confirmed SMARCA4 mutations, and xenografts producing elevated serum calcium (hypercalcemia hallmark); SMARCA4 re-expression dramatically suppressed growth. **TOV-112D — historically a grade 3 endometrioid carcinoma — is reclassified as dedifferentiated ovarian carcinoma** (SMARCA4 + SMARCA2 loss in undifferentiated regions, transition from glandular to solid areas, TP53 missense with diffuse p53). Both share gene-expression profiles and EZH2-inhibitor sensitivity.
- **Relevance:** High — DIRECTLY confirms the seed; critical assay-aware caveat for the resource (see Breadth Flag). BIN67, SCCOHT1 and COV434 are the standard SCCOHT models (SMARCA4-mutant, SMARCA2-null).
- **Source type:** Primary

#### Genomic simplicity of SCCOHT — single driver, diploid, rhabdoid-like
- **Citation:** Pautier-adjacent comprehensive analysis: Lang JD, Hendricks WPD, et al. (2020). "Small Cell Carcinoma of the Ovary, Hypercalcemic Type (SCCOHT) beyond SMARCA4 Mutations: A Comprehensive Genomic Analysis." *Cells* 9(6):1496. doi:10.3390/cells9061496 (PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7349095/)
- **Key finding:** SCCOHT shows **remarkable genomic stability — diploid profiles and low mutation load (mean ~5.43 mutations/Mb)**; **SMARCA4 is essentially the only recurrently mutated gene**. Oncogenesis is driven by epigenetic dysregulation rather than genomic instability. Related work notes SCCOHT genomic/epigenomic signatures resemble atypical teratoid/rhabdoid tumor (ATRT) more than HGSC, supporting the view that SCCOHT is a **malignant rhabdoid tumor of the ovary**.
- **Relevance:** High — for a multi-omic descriptor, SCCOHT is the extreme "one-driver, quiet-genome" contrast to HGSC; useful framing for WES/CNV interpretation.
- **Source type:** Primary (genomic analysis)

---

### Q8 — Defining molecular features/markers of the other rare subtypes

#### OVARIAN CLEAR CELL — ARID1A discovery (Wiegand 2010, NEJM)
- **Citation:** Wiegand KC, Shah SP, Al-Agha OM, et al. (2010). "ARID1A Mutations in Endometriosis-Associated Ovarian Carcinomas." *New England Journal of Medicine* 363(16):1532–1543. doi:10.1056/NEJMoa1008433 https://www.nejm.org/doi/full/10.1056/NEJMoa1008433
- **Key finding:** Somatic **ARID1A** mutations (encoding BAF250a, a SWI/SNF subunit) in **46% (55/119) of ovarian clear cell**, **30% (10/33) of endometrioid**, and **0% (0/76) of high-grade serous** carcinomas. ARID1A mutation / BAF250a loss found in contiguous atypical endometriosis — an early event in malignant transformation of endometriosis.
- **Relevance:** High — canonical clear-cell/endometrioid driver; ties both subtypes to SWI/SNF (bridges to Q9).
- **Source type:** Primary

#### OVARIAN CLEAR CELL — ARID1A discovery (Jones 2010, Science, parallel report)
- **Citation:** Jones S, Wang TL, Shih IeM, et al. (2010). "Frequent Mutations of Chromatin Remodeling Gene ARID1A in Ovarian Clear Cell Carcinoma." *Science* 330(6001):228–231. doi:10.1126/science.1196333
- **Key finding:** Exome sequencing identified **ARID1A** as one of the most frequently mutated genes in ovarian clear cell carcinoma (secondary sources quote **~50%**; the companion oncogene **PPP2R1A** was also recurrently mutated). Independent, simultaneous discovery to Wiegand 2010. *(Exact fraction in the original paper not confirmed in this session's searches — see Gaps.)*
- **Relevance:** High — co-founding ARID1A paper.
- **Source type:** Primary

#### OVARIAN CLEAR CELL — PIK3CA and HNF1B / oxidative-stress biology
- **Citation (PIK3CA):** Kuo KT, Mao TL, Jones S, et al. (2009). "Frequent Activating Mutations of PIK3CA in Ovarian Clear Cell Carcinoma." *Am J Pathol* (via ResearchGate record: https://www.researchgate.net/publication/24260408). **Key finding:** frequent activating PIK3CA mutations in OCCC; later work shows **ARID1A + PIK3CA co-mutation cooperates** to drive tumorigenesis via sustained IL-6 (Chandler et al. 2015, *Nat Commun* 6:6118, https://www.nature.com/articles/ncomms7118).
- **Citation (HNF1B/oxidative stress):** Amano Y, et al. and Okamoto/Konishi group work: "HNF1β drives glutathione (GSH) synthesis underlying intrinsic carboplatin resistance of ovarian clear cell carcinoma." *Tumor Biology* (2015), doi via PubMed https://pubmed.ncbi.nlm.nih.gov/26520442/ ; and "Ovarian clear cell carcinoma meets metabolism; HNF-1β confers survival benefits through the Warburg effect and ROS reduction" (PMC4741562, https://pmc.ncbi.nlm.nih.gov/articles/PMC4741562/).
- **Key finding:** **HNF1B (HNF-1β) expression is a hallmark of OCCC** (a ~320-gene OCCC signature centers on HNF1β signaling and oxidative-stress genes). HNF1β drives **glutathione synthesis, anaerobic glycolysis (Warburg), ROS reduction**, underpinning intrinsic platinum resistance. OCCC arises from iron-rich endometriotic cysts → oxidative stress/DNA damage; low ARID1A correlates with high 8-hydroxyguanosine (oxidative-stress marker).
- **Relevance:** High — HNF1B, glutathione/oxidative-stress and PIK3CA/ARID1A are the expected clear-cell "validation" signals in the resource (n=7 clear cell lines).
- **Source type:** Primary

#### MUCINOUS — KRAS, HER2/ERBB2 amplification, CDX2/intestinal markers
- **Citation (KRAS/HER2):** Anglesio MS-adjacent / "The status of Her2 amplification and Kras mutations in mucinous ovarian carcinoma." *Human Genomics* (2016) 10:40. doi:10.1186/s40246-016-0096-9 https://humgenomics.biomedcentral.com/articles/10.1186/s40246-016-0096-9 . Also Mackenzie R et al. (2015), "Targeted deep sequencing of mucinous ovarian tumors...RAS-pathway activating mutations." *BMC Cancer* 15:415 (PMC4494777).
- **Key finding:** **KRAS mutation is the dominant alteration — 64.9% of mucinous carcinomas** (92.3% of mucinous borderline tumors); literature range 13–60% across cohorts. **HER2/ERBB2 amplification in ~19% of invasive mucinous** (6% of borderline; over-expression estimates 18–35%). KRAS mutation and HER2 amplification are near **mutually exclusive**.
- **Citation (CDX2/mucins):** Groisman/others: "Immunohistochemical expression of CDX2 in primary ovarian mucinous tumors..." *Modern Pathology* (https://www.nature.com/articles/3800698) and mucin panel study (PMC2386514).
- **Key finding:** **CDX2 (intestinal-differentiation TF) positive in ~79% of primary mucinous ovarian tumors**; intestinal-type immunophenotype MUC2+/CK20+/CDX2+ with MUC5AC (gastric pattern decreases adenoma→carcinoma as MUC2 rises).
- **Relevance:** High — KRAS/ERBB2 and CDX2/MUC markers are the expected mucinous signals (n=3 mucinous lines); note strong overlap with GI/colorectal phenotype (metastasis-mimic caveat).
- **Source type:** Primary

#### LOW-GRADE SEROUS — MAPK pathway (KRAS/BRAF/NRAS), distinct from HGSC
- **Citation:** Review/genomic sources: "Current concept of low-grade serous ovarian carcinoma" (PMC10894346) and "Molecular changes driving low-grade serous ovarian cancer..." *Int J Gynecol Cancer* (PMC11503204); "Genomic Analysis of Low-Grade Serous Ovarian Cancer" (Cureus).
- **Key finding:** **Mutually exclusive activating mutations in KRAS, BRAF, or NRAS in ~50–60% of LGSC**, arising early in serous borderline tumors. Reported ranges: KRAS 16–44%, BRAF 2–20%, NRAS up to 26%; one genomic series had **MAPK mutations in 60% (KRAS 33%, BRAF 11%, NRAS 11%)**. Hotspots KRAS G12R/C/D/V, BRAF V600E, NRAS Q61R/K. **Distinct from HGSC**, which is defined by near-universal TP53 mutation and high genomic instability; KRAS mutation and MAPK activation (IHC positive 63.6% LGSC vs 17.1% HGSC) are far more frequent in LGSC.
- **Relevance:** High — MAPK/RAS activation is the LGSC signature and the key discriminator from HGSC on cell lines.
- **Source type:** Review + primary genomic

#### ENDOMETRIOID — CTNNB1/β-catenin, ARID1A, PIK3CA/PTEN, MMR/POLE
- **Citation:** Hollis RL, Thomson JP, Stanley B, et al. (2020). "Molecular stratification of endometrioid ovarian carcinoma predicts clinical outcome." *Nature Communications* 11:4995. doi:10.1038/s41467-020-18819-5 https://www.nature.com/articles/s41467-020-18819-5
- **Key finding:** WES of 112 endometrioid ovarian carcinomas: **CTNNB1 43%, PIK3CA 43%, ARID1A 36%, PTEN 29%, KRAS 26%, TP53 26%, SOX8 19%**; **POLE 6% and mismatch-repair genes 18%** (co-occurring). TP53 and CTNNB1 mutations are **largely mutually exclusive** (CTNNB1 rate 56.6% in TP53-wildtype); ARID1A co-occurs with PI3K/AKT activation. Defines a TCGA-like four-group molecular stratification (POLE / MMRd / p53-abnormal / no-specific-molecular-profile).
- **Relevance:** High — CTNNB1/Wnt + ARID1A + PIK3CA are the expected endometrioid signals (n=2 EC lines); MMR/POLE relevant to hypermutation status in WES.
- **Source type:** Primary

---

### Q9 — Shared SWI/SNF (BAF) disruption across gynecologic malignancies; the EC–SCCOHT convergence question

#### McCluggage & Stewart 2021 — SWI/SNF-deficiency as a UNIFYING theme (the key synthesis)
- **Citation:** McCluggage WG, Stewart CJR (2021). "SWI/SNF-deficient malignancies of the female genital tract." *Seminars in Diagnostic Pathology* 38(3):199–211. doi:10.1053/j.semdp.2020.08.003 (PubMed: https://pubmed.ncbi.nlm.nih.gov/32978032/)
- **Key finding:** Frames diverse gynecologic entities as sharing **convergent SWI/SNF-complex disruption**: **SCCOHT** (SMARCA4 mutation ~98%, + SMARCA2 loss); **endometrioid & clear cell carcinomas** (ARID1A); **undifferentiated/dedifferentiated carcinomas of endometrium AND ovary** (SMARCA4, SMARCB1, ARID1A, ARID1B loss); SMARCB1-deficient vulvar epithelioid sarcoma/myoepithelial carcinoma; SMARCA4-deficient undifferentiated uterine sarcoma; anaplastic mural nodules in mucinous neoplasms. Establishes "SWI/SNF-deficient gynecologic malignancy" as an **established, named organizing concept**.
- **Relevance:** High — this is the literature backbone for the resource's SWI/SNF-convergence narrative; shows the theme is ESTABLISHED, not novel.
- **Source type:** Review

#### Karnezis et al. 2016 — SWI/SNF loss drives dedifferentiation in endometrial carcinoma
- **Citation:** Karnezis AN, Hoang LN, Coatham M, et al. (2016). "Loss of switch/sucrose non-fermenting complex protein expression is associated with dedifferentiation in endometrial carcinomas." *Modern Pathology* 29:302–314. doi:10.1038/modpathol.2015.155
- **Key finding:** In dedifferentiated endometrial carcinoma, loss of **SMARCA4 (or SMARCB1) together with concurrent SMARCA2 loss** is associated with histologic **dedifferentiation in ~half of cases**; the undifferentiated component carries frameshift/nonsense SMARCA4 mutations while the adjacent low-grade endometrioid component retains SMARCA4 — i.e., an ARID1A/SWI/SNF-mutant endometrioid tumor can **progress to a SMARCA4/SMARCA2-null (SCCOHT-like) state**.
- **Relevance:** High — the mechanistic bridge: endometrioid → SWI/SNF-null dedifferentiation, the same dual-ATPase-loss state that defines SCCOHT.
- **Source type:** Primary

#### Dedifferentiated/undifferentiated OVARIAN carcinoma — SWI/SNF inactivation
- **Citation:** (Modern Pathology 2023). "Dedifferentiated and Undifferentiated Ovarian Carcinoma: An Aggressive and Molecularly Distinct Ovarian Tumor Characterized by Frequent SWI/SNF Complex Inactivation." *Modern Pathology* — https://www.modernpathology.org/article/S0893-3952(23)00279-X/fulltext
- **Key finding:** Dedifferentiated/undifferentiated OVARIAN carcinomas are a molecularly distinct, aggressive group defined by **frequent SWI/SNF complex inactivation** (SMARCA4/SMARCA2/SMARCB1/ARID1A/ARID1B). Extends the endometrial dedifferentiation concept to the ovary — directly relevant to TOV-112D's reclassification (Q7).
- **Relevance:** High — the ovarian counterpart establishing that an "EC" ovarian line can be a SWI/SNF-deficient dedifferentiated tumor.
- **Source type:** Primary

#### Tessier-Cloutier et al. — prognostic weight of SWI/SNF deficiency
- **Citation:** Discussed in the SWI/SNF review context (McCluggage & Stewart 2021) and Tessier-Cloutier et al. Histopathology 2022 (https://onlinelibrary.wiley.com/doi/10.1111/his.14639); ARID1B IHC diagnostic work (PMC10486746).
- **Key finding:** Among stage III/IV dedifferentiated/undifferentiated endometrial carcinomas, detected SWI/SNF alteration associates with **median survival ~4 months vs ~36 months** when absent; ~half of SWI/SNF-deficient DDC/UDC show ARID1B loss; ARID1B IHC is an important diagnostic test.
- **Relevance:** Medium — supports clinical/biological coherence of the SWI/SNF-deficient dedifferentiated group.
- **Source type:** Primary

#### Shared therapeutic vulnerability — EZH2 synthetic lethality (mechanistic convergence)
- **Citation:** Kim KH & Roberts CWM and follow-on: "SWI/SNF catalytic subunits' switch drives resistance to EZH2 inhibitors in ARID1A-mutated cells." *Nature Communications* (2018) 9:4363 (PMC6175882, https://www.nature.com/articles/s41467-018-06656-6); Sasaki et al. (2020) "Synthetic lethal therapy based on targeting the vulnerability of SWI/SNF...-deficient cancers." *Cancer Science*, doi:10.1111/cas.14311.
- **Key finding:** SWI/SNF loss creates a **shared EZH2 dependency**: ARID1A inactivation is synthetically lethal with EZH2 inhibition; in SMARCA4-loss tumors, unopposed EZH2 deposits H3K27me3 to silence SMARCA2, and **EZH2 inhibitors reactivate SMARCA2 and restore differentiation in SMARCB1-deficient rhabdoid tumors and SMARCA4-deficient SCCOHT**. ARID1A-mutant (clear cell/endometrioid) and SMARCA4-deficient (SCCOHT) tumors thus converge on overlapping epigenetic vulnerabilities (EZH2, and SMARCA2/BRG1 re-expression; BET/ATR as orthogonal strategies).
- **Relevance:** High — functional (not just descriptive) evidence that ARID1A-mutant subtypes and SCCOHT converge; the COV434/TOV-112D pair share EZH2-inhibitor sensitivity (Karnezis 2021).
- **Source type:** Primary + review

#### Recent preprint — cell-state dynamics of SMARCA4-null dedifferentiated endometrial cancer
- **Citation:** (2025 preprint). "Cell State Chaos Underpins the Evolution of SMARCA4-Deficient Dedifferentiated Endometrial Cancer." *bioRxiv* — https://www.biorxiv.org/content/10.64898/2025.12.01.690470
- **Key finding:** SMARCA4-knockout cells move erratically between transcriptional states ("cell-state chaos") vs ordered trajectories in wild-type — proposes dedifferentiation arises from transient, disordered state acquisition after SWI/SNF loss.
- **Relevance:** Medium — frontier mechanistic framing; NOT peer-reviewed (preprint), flag accordingly.
- **Source type:** Primary (preprint, unreviewed)

---

## Search Log
- Query 1: "Ramos 2014 Nature Genetics SMARCA4 mutations small cell carcinoma ovary hypercalcemic type SCCOHT" — 9 examined, 5 relevant
- Query 2: "Witkowski Jelinic 2014 Nature Genetics SMARCA4 germline somatic mutations SCCOHT" — 9 examined, 4 relevant
- Query 3: "SCCOHT dual SMARCA4 SMARCA2 BRG1 BRM loss immunohistochemistry Karnezis diagnostic marker" — 8 examined, 4 relevant
- Query 4: "COV434 cell line reclassified granulosa cell tumor SCCOHT SMARCA4 BIN67" — 7 examined, 4 relevant
- Query 5: "ovarian clear cell carcinoma ARID1A PIK3CA mutation frequency HNF1B oxidative stress marker" — 7 examined, 5 relevant
- Query 6: "mucinous ovarian carcinoma KRAS mutation CDX2 HER2 ERBB2 amplification frequency" — 9 examined, 4 relevant
- Query 7: "low-grade serous ovarian carcinoma KRAS BRAF NRAS MAPK mutation frequency distinct high-grade" — 10 examined, 5 relevant
- Query 8: "ovarian endometrioid carcinoma CTNNB1 beta-catenin ARID1A PIK3CA mutation frequency mismatch repair POLE" — 10 examined, 3 relevant
- Query 9: "Wiegand 2010 NEJM ARID1A mutations ovarian clear cell endometrioid carcinoma" — 9 examined, 2 relevant
- Query 10: "Jones 2010 Science ARID1A frequent mutations ovarian clear cell carcinoma exome" — 8 examined, 2 relevant
- Query 11: "SCCOHT diploid genome single driver SMARCA4 low mutation burden genomically stable malignant rhabdoid" — 8 examined, 4 relevant
- Query 12: "SMARCA4 SMARCB1 deficient dedifferentiated undifferentiated endometrial carcinoma Kolin Karnezis ARID1A" — 8 examined, 5 relevant
- Query 13: "Tessier-Cloutier SWI/SNF deficient gynecologic carcinoma ARID1A SMARCA4 co-loss dedifferentiated" — 7 examined, 4 relevant
- Query 14: "mucinous ovarian carcinoma CDX2 intestinal differentiation MUC2 MUC5AC immunohistochemistry markers" — 8 examined, 3 relevant
- Query 15: "ovarian clear cell carcinoma MET amplification glutathione metabolism hepatocyte HNF1B gene expression signature" — 7 examined, 4 relevant (MET not found)
- Query 16: "SWI/SNF deficient cancer EZH2 inhibitor synthetic lethality shared vulnerability ARID1A SMARCA4 ovarian therapeutic" — 8 examined, 4 relevant
- Query 17: "Jones ARID1A ovarian clear cell Science 2010 57% / 24 of 42 exome PPP2R1A" — 9 examined, 1 relevant (citation confirmed; exact fraction not surfaced)
- Query 18: "endometrioid carcinoma dedifferentiation SMARCA4 loss SCCOHT-like convergent SWI/SNF ovarian shared cell state" — 8 examined, 4 relevant
- Query 19: "Karnezis 'Loss of SWI/SNF...dedifferentiation in endometrial carcinomas' Modern Pathology 2016 authors volume" — 8 examined, 3 relevant
- WebFetch A: PubMed 33328126 (Karnezis 2021 COV434/TOV-112D) — full citation + reclassification details confirmed
- WebFetch B: PubMed 32978032 (McCluggage & Stewart 2021 review) — full citation + convergence framing confirmed
- WebFetch C: PMC4832362 (Karnezis 2016 dual-loss J Pathol) — full citation + sensitivity/specificity numbers confirmed
- WebFetch D: nature.com/modpathol2015155 — paywall redirect (citation instead confirmed via Query 19)

## Gaps
1. **Jones et al. 2010 (Science) exact ARID1A frequency** — citation confirmed (Science 330:228–231) but the precise fraction (commonly cited as ~57%, 24/42) was NOT surfaced in this session's searches; multiple secondary sources give "~50%". Report as approximate pending primary-text verification.
2. **MET in ovarian clear cell carcinoma** — the seed mentioned MET, but searches did not surface primary MET-amplification frequency data for OCCC. HNF1B/glutathione/oxidative-stress biology is well-supported; MET-specific claim unverified — do not assert.
3. **Mucinous CDX2 as a mutation/CNV target** — CDX2 is supported as an IHC intestinal-differentiation marker (~79% positive), not as a genomic alteration. Frequency data are IHC-based.
4. **Direct "EC–SCCOHT convergence" as an explicitly named concept** — I found NO paper that names an "endometrioid-carcinoma-to-SCCOHT" convergence per se. What IS established: (a) shared SWI/SNF/BAF disruption across gynecologic cancers (McCluggage & Stewart 2021), (b) ARID1A-mutant endometrioid tumors dedifferentiating to SMARCA4/SMARCA2-null states (Karnezis 2016), (c) dedifferentiated OVARIAN carcinoma with SWI/SNF inactivation (Mod Pathol 2023), (d) shared EZH2/SMARCA2 vulnerability. Treat the resource's specific transcriptomic EC–SCCOHT convergence as an EXTENSION of an established theme, not a wholly novel concept.
5. Low-grade serous **seminal two-tier / KRAS-BRAF origin papers** (Singer/Kurman; Gershenson) were not directly retrieved; LGSC frequencies came from recent reviews/genomic series rather than the original discovery papers.

## Breadth Flag
CRITICAL assay-aware caveat for the Data Descriptor team, arising directly from Karnezis et al. 2021 (Gynecol Oncol 160:568–578):

- **The n=2 "EC" lines and the EC–SCCOHT convergence signal may be partly an identity artifact.** TOV-112D — one of the most widely used "endometrioid" ovarian lines — was RECLASSIFIED as **dedifferentiated ovarian carcinoma with SMARCA4/SMARCA2 loss**, and COV434 (historically "granulosa") is actually **SCCOHT**. If the prior analysis's EC–SCCOHT convergence is driven by TOV-112D, the convergence would substantially RECAPITULATE a known misannotation/biology (a SWI/SNF-null dedifferentiated line clustering with SWI/SNF-null SCCOHT lines) rather than reveal a novel EC-intrinsic link. **Recommend the team explicitly check the identity/SWI/SNF status (ARID1A, SMARCA4, SMARCA2 expression) of each EC and SCCOHT line before interpreting the convergence as biology.** This is the single most important cross-cutting finding for this cluster and overlaps the "cell-line panels/authentication" cluster — flag to that team.
# Q10 — ADC targets in ovarian cancer: clinical landscape & value of subtype-resolved target-expression atlases

Literature specialist raw findings. Every claim below is grounded in a WebSearch/WebFetch performed this session (July 2026). Papers reported individually; no synthesis. Calibrated confidence and gaps flagged at the end.

---

## Findings

### FOLR1 (folate receptor alpha) — mirvetuximab soravtansine (Elahere)

- **Citation:** Matulonis UA, et al. (2023). "Efficacy and Safety of Mirvetuximab Soravtansine in Patients With Platinum-Resistant Ovarian Cancer With High Folate Receptor Alpha Expression: Results From the SORAYA Study." *Journal of Clinical Oncology* 41(13). DOI: 10.1200/JCO.22.01900 (https://ascopubs.org/doi/10.1200/JCO.22.01900)
- **Key finding:** Single-arm phase II SORAYA (NCT04296890) in FRα-high platinum-resistant ovarian cancer (PROC), 1–3 prior lines including required bevacizumab. Investigator-assessed confirmed ORR 32.4% (95% CI 23.6–42.2; 5 CR + 29 PR among 105 efficacy-evaluable), median duration of response 6.9 months (95% CI 5.6–9.7). FRα-high defined as ≥75% of tumor cells with ≥2+ IHC intensity by the VENTANA FOLR1 (FOLR1-2.1) RxDx companion assay. Basis for accelerated approval.
- **Relevance:** High — primary registrational efficacy dataset; defines the FRα-high companion-diagnostic cutoff our resource should benchmark against.
- **Source type:** Primary trial

- **Citation:** [SORAYA final OS follow-up] Matulonis UA, et al. (2024). "Mirvetuximab soravtansine in folate receptor alpha (FRα)–high platinum-resistant ovarian cancer: final overall survival and post hoc sequence of therapy subgroup results from the SORAYA trial." PMC11347190 (https://pmc.ncbi.nlm.nih.gov/articles/PMC11347190/)
- **Key finding:** Final median overall survival in SORAYA was 15.0 months (95% CI 11.5–18.7). Confirms the ORR 32.4% / DOR 6.9-month primary readout.
- **Relevance:** Medium — mature survival data for the accelerated-approval cohort.
- **Source type:** Primary trial (follow-up)

- **Citation:** Moore KN, et al. (2023). "Mirvetuximab Soravtansine in FRα-Positive, Platinum-Resistant Ovarian Cancer." *New England Journal of Medicine* 389(23):2162–2174. PMID: 38055253. DOI: 10.1056/NEJMoa2309169 (https://www.nejm.org/doi/full/10.1056/NEJMoa2309169)
- **Key finding:** Phase III MIRASOL (GOG 3045/ENGOT-ov55), 453 FRα-high PROC pts randomized 1:1 to mirvetuximab vs investigator's-choice chemotherapy (ICC). ORR 42% vs 16% (P<0.0001); median PFS 5.62 vs 3.98 months (HR 0.65, P<0.0001); statistically significant OS benefit — the first novel agent to show an OS benefit in PROC in a phase III trial. Updated/full-approval dataset (median follow-up 30.5 months): median OS 16.85 vs 13.34 months (~32% reduction in risk of death). FRα-high selection by the same VENTANA FOLR1 RxDx assay.
- **Relevance:** High — confirmatory phase III; the pivotal efficacy + subtype/biomarker-selection evidence.
- **Source type:** Primary trial

- **Citation:** U.S. FDA / AbbVie (2022, 2024). Accelerated approval press coverage (CancerNetwork) and full-approval release: "FDA Grants Full Approval for ELAHERE (mirvetuximab soravtansine-gynx), Mar 22, 2024." (https://news.abbvie.com/2024-03-22-U-S-Food-and-Drug-Administration-FDA-Grants-Full-Approval-for-ELAHERE-R-mirvetuximab-soravtansine-gynx-for-Certain-Ovarian-Cancer-Patients) and FDA drug page (https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-mirvetuximab-soravtansine-gynx-fra-positive-platinum-resistant-epithelial-ovarian)
- **Key finding:** Accelerated approval 14 Nov 2022 (based on SORAYA); full/regular approval 22 Mar 2024 (based on MIRASOL) for adults with FRα-positive, platinum-resistant epithelial ovarian, fallopian tube, or primary peritoneal cancer with 1–3 prior systemic regimens. First ADC approved in ovarian cancer.
- **Relevance:** High — regulatory milestone; the "target X has an approved drug + companion Dx" anchor for the FOLR1 utility narrative.
- **Source type:** Regulatory

- **Citation:** FDA Approval Summary. Nguyen AM/Matulonis-context review (2023). "FDA Approval Summary: Mirvetuximab Soravtansine-Gynx for FRα-Positive, Platinum-Resistant Ovarian Cancer." *Clinical Cancer Research* 29(19):3835. (https://aacrjournals.org/clincancerres/article/29/19/3835/729087/) — plus subtype review: "Folate Receptor Alpha in Advanced Epithelial Ovarian Cancer: Diagnostic Role and Therapeutic Implications of a Clinically Validated Biomarker" (2025), PMC12154392 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12154392/)
- **Key finding — subtype association:** FRα/FOLR1 is expressed in ~90% of ovarian carcinomas and is most consistently over-expressed in the high-grade serous histotype (the predominant subtype studied). Companion Dx positivity = ≥75% of viable tumor cells with moderate/strong (≥2+) membrane staining.
- **Relevance:** High — direct histotype-expression association for FOLR1; frames why a subtype-resolved map is clinically meaningful.
- **Source type:** Regulatory / Review

- **Citation:** McHenry A, et al. (2026). "Folate receptor alpha (FRα/FOLR1) and HER2 immunohistochemical staining in high-grade endometrial carcinoma with aberrant p53 expression." *Histopathology*. (https://onlinelibrary.wiley.com/doi/10.1111/his.15533)
- **Key finding:** Co-examines FRα and HER2 IHC in high-grade (p53-aberrant) endometrial carcinoma — relevant read-across for endometrioid/serous-like gynecologic histology and for pairing FOLR1 + HER2 target maps.
- **Relevance:** Medium — adjacent-histology expression context.
- **Source type:** Primary (pathology cohort)

---

### NaPi2b / SLC34A2 — upifitamab rilsodotin (UpRi; XMT-1536)

- **Citation:** Richardson DL, et al. (2024). "UPLIFT (ENGOT-Ov67/GOG-3048): Results from the phase II trial of upifitamab rilsodotin (UpRi; XMT-1536), a NaPi2b-directed dolaflexin antibody-drug conjugate in platinum-resistant ovarian cancer." *Gynecologic Oncology* (SGO 2024 abstract). (https://www.gynecologiconcology-online.net/article/S0090-8258(24)00440-2/abstract)
- **Key finding:** Single-arm phase II, 268 HGSOC PROC pts (141 NaPi2b-high, 127 NaPi2b-low). Did NOT meet primary endpoint: ORR 15.6% (95% CI 10.0–22.7) in NaPi2b-high, 10.2% in NaPi2b-low, 13.1% overall; median DOR 7.4 months in NaPi2b-high. NaPi2b-high defined as tumor proportion score (TPS) ≥75. NaPi2b (SLC34A2), a sodium-dependent phosphate transporter, is highly expressed in HGSOC/fallopian-tube/primary-peritoneal cancer.
- **Relevance:** High — target validated at the expression level but ADC failed efficacy bar; illustrates that high target expression ≠ clinical benefit (a caution for the atlas narrative).
- **Source type:** Primary trial (abstract)

- **Citation:** Hamilton E, et al. (2023). "UP-NEXT (GOG-3049/ENGOT-Ov71-NSGO-CTU): A study of upifitamab rilsodotin (UpRi), a NaPi2b-directed ADC, in platinum-sensitive recurrent ovarian cancer." *Journal of Clinical Oncology* (TPS abstract). NCT05329545. (https://ascopubs.org/doi/10.1200/JCO.2023.41.16_suppl.TPS5614) and ClinicalTrials.gov (https://clinicaltrials.gov/study/NCT05329545)
- **Key finding:** Phase III, double-blind, placebo-controlled 2:1 maintenance study of UpRi monotherapy in NaPi2b-high (TPS ≥75) platinum-sensitive recurrent HGSOC; primary endpoint PFS. FDA placed a partial clinical hold on UP-NEXT and the phase I UPGRADE-A trial over bleeding events; the upifitamab rilsodotin program was subsequently discontinued in ovarian cancer.
- **Relevance:** High — documents that NaPi2b-directed ADC development in ovarian cancer has ended (safety + efficacy), important context for how our resource frames NaPi2b.
- **Source type:** Primary trial (protocol) / Regulatory (clinical hold, discontinuation — CancerNetwork/OncLive, https://www.onclive.com/view/journey-ends-for-upifitamab-rilsodotin-in-napi2b-platinum-resistant-ovarian-cancer)

- **Citation:** [Phase 1 expansion] (2022). "Updated Results from the Phase 1 Expansion Study of Upifitamab Rilsodotin (UpRi; XMT-1536), a NaPi2b-directed Dolaflexin ADC in Ovarian Cancer (076)." *Gynecologic Oncology* abstract. (https://www.sciencedirect.com/science/article/abs/pii/S009082582201294X)
- **Key finding:** Earlier phase I expansion that supported advancing UpRi in NaPi2b-expressing ovarian cancer (pre-UPLIFT signal).
- **Relevance:** Medium — early-phase target-expression rationale.
- **Source type:** Primary trial (abstract)

---

### TROP2 / TACSTD2 — sacituzumab govitecan; datopotamab deruxtecan (+ related TROP2 ADCs)

> Note: "sacituzumab govitecan" (SG / Trodelvy; SN-38 payload) and "sacituzumab tirumotecan" (sac-TMT / SKB264 / MK-2870; a distinct belotecan-derivative TROP2 ADC) are DIFFERENT agents — kept separate below.

- **Citation:** [SG retrospective] (2024). "Sacituzumab govitecan in heavily pretreated, platinum-resistant high grade serous ovarian cancer." *Gynecologic Oncology Reports*. PMID: 39108617. (https://pmc.ncbi.nlm.nih.gov/articles/PMC11300917/)
- **Key finding:** Retrospective series of sacituzumab govitecan (TROP2-SN38 ADC) used off-label in heavily pretreated platinum-resistant HGSOC — reports feasibility/activity signal; SG is NOT FDA-approved in ovarian cancer (approved in TNBC and urothelial).
- **Relevance:** Medium — real-world TROP2-ADC activity signal in HGSOC.
- **Source type:** Primary (retrospective)

- **Citation:** Yale investigators (2026). "A phase II evaluation of sacituzumab govitecan in platinum-resistant ovarian cancer (NCT06028932)." *Journal of Clinical Oncology* 44(16_suppl):5575 (ASCO 2026 abstract). (https://ascopubs.org/doi/10.1200/JCO.2026.44.16_suppl.5575)
- **Key finding:** Prospective single-arm phase II (n≈20, median age 67) of SG 10 mg/kg d1/d8 q21d in PROC; primary endpoint ORR (RECIST 1.1). Abstract notes >60% of HGSOC over-express (≥2+ IHC) TROP2.
- **Relevance:** Medium — ongoing prospective SG trial; TROP2 prevalence statement useful.
- **Source type:** Primary trial (abstract)

- **Citation:** Perrone E, et al. (Santin lab) (2020). "Preclinical Activity of Sacituzumab Govitecan, an ADC Targeting Trophoblast Cell-Surface Antigen 2 (Trop-2)... in Ovarian Cancer." PMC7028697. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7028697/)
- **Key finding:** Trop-2 assessed by IHC (FFPE tumors) and flow cytometry (cell lines): moderate-to-strong staining in 47% of ovarian tumors by IHC, while 89% of primary ovarian cancer cell lines over-expressed Trop-2 by flow cytometry — SG showed preclinical antitumor activity.
- **Relevance:** High — direct precedent for using an ovarian cell-line panel to profile a TROP2 ADC target; supports our reuse narrative.
- **Source type:** Preclinical primary

- **Citation:** Oaknin A, et al. (2024). "Datopotamab deruxtecan (Dato-DXd) in patients with endometrial (EC) or ovarian cancer (OC): Results from the phase II TROPION-PanTumor03 study." *Annals of Oncology* 35(suppl); ESMO 2024 abstract 714MO. (https://www.annalsofoncology.org/article/S0923-7534(24)02295-6/fulltext)
- **Key finding:** TROP2-directed ADC Dato-DXd 6 mg/kg q3w. Ovarian cohort (n=35): confirmed ORR 42.9% (95% CI 26.3–60.6; 2.9% CR + 40.0% PR) at 14.5-month median follow-up; endometrial ORR 27.5%. Ocular surface events in 40% (ovarian); grade 3 drug-related ILD in 1 pt/cohort.
- **Relevance:** High — strongest prospective TROP2-ADC efficacy signal in ovarian cancer.
- **Source type:** Primary trial (abstract)

- **Citation:** [Sacituzumab tirumotecan / sac-TMT] ESMO 2024 phase II (distinct agent). Reported via OncLive/ASCO coverage. (https://www.onclive.com/... ; secondary)
- **Key finding:** Sacituzumab tirumotecan monotherapy in advanced PROC: ORR 40%, median PFS 6.0 months, DCR 75%; ORR rose to 61.5% in high-TROP2 patients — a subtype/expression-response association for a TROP2 ADC.
- **Relevance:** Medium — separate TROP2 ADC; the high-TROP2 → higher-ORR association is directly relevant to expression-guided model selection. (Secondary source; primary abstract not fetched.)
- **Source type:** Primary trial (via secondary coverage)

- **Citation:** [TROP2 preclinical, endometrial/carcinosarcoma read-across] (2025). "Preclinical Activity of Datopotamab Deruxtecan (Dato-DXd)... in Poorly Differentiated Endometrial Carcinomas." *Cancer Research Communications* 5(9):1611. (https://aacrjournals.org/cancerrescommun/article/5/9/1611/765111/) and "Dato-DXd... against primary and metastatic uterine and ovarian TROP2 over-expressing carcinosarcoma." PMID: 40582040.
- **Key finding:** TROP2 over-expression and Dato-DXd activity extend to poorly differentiated endometrial carcinoma and uterine/ovarian carcinosarcoma models — TROP2 relevant beyond HGSOC.
- **Relevance:** Medium — broadens TROP2 subtype relevance to rarer/mixed histologies.
- **Source type:** Preclinical primary

---

### HER2 / ERBB2 — trastuzumab deruxtecan (T-DXd; Enhertu)

- **Citation:** Meric-Bernstam F, et al. (2024). "Efficacy and Safety of Trastuzumab Deruxtecan in Patients With HER2-Expressing Solid Tumors: Primary Results From the DESTINY-PanTumor02 Phase II Trial." *Journal of Clinical Oncology* 42(1):47–58. PMID: 37870536. DOI: 10.1200/JCO.23.02005 (https://ascopubs.org/doi/10.1200/JCO.23.02005)
- **Key finding:** Basket phase II, 267 pts across 7 HER2-expressing (IHC 3+/2+) solid-tumor cohorts. Ovarian cohort investigator-assessed ORR 45.0%; in HER2 IHC 3+ ovarian, ORR ~64%. Overall (all cohorts) ORR 37.1% (95% CI 31.3–43.2), and overall HER2 IHC 3+ ORR 61.3% (95% CI 49.4–72.4) — response strongly enriched in IHC 3+. Supported the tumor-agnostic HER2 IHC 3+ approval.
- **Relevance:** High — pivotal HER2-ADC efficacy in ovarian cancer with a clear IHC-expression → response gradient.
- **Source type:** Primary trial

- **Citation:** [DESTINY-PanTumor02 gynecologic analysis] (2026). "Trastuzumab deruxtecan in HER2-expressing gynecologic cancers from DESTINY-PanTumor02: antitumor activity, safety, and exploratory biomarker analyses." *J Gynecol Oncol* 37:e65. PMID: 42135959. (https://ejgo.org/DOIx.php?id=10.3802%2Fjgo.2026.37.e65)
- **Key finding:** Gynecologic-focused subgroup + exploratory biomarker analysis (endometrial/cervical/ovarian) with HER2 IHC status stratification.
- **Relevance:** High — gyn-specific HER2 activity and biomarker detail.
- **Source type:** Primary trial (subgroup)

- **Citation:** McAlpine JN/Köbel-context; Chao A, et al. (2018). "HER2 Is Frequently Over-expressed in Ovarian Clear Cell Adenocarcinoma..." PMC5926901 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5926901/)
- **Key finding — subtype association:** HER2 over-expressed in 42.9% of ovarian clear cell carcinoma vs 20.8% serous, 23.1% endometrioid, 30.0% mucinous (clear cell significantly higher, P=0.026). A companion series reports OCCC ~43% vs serous 21%, endometrioid 23%, mucinous 30%.
- **Relevance:** High — quantitative HER2-by-histotype map; directly supports "pick a clear-cell/mucinous model to study HER2-ADC."
- **Source type:** Primary (IHC cohort)

- **Citation:** McAlpine JN, et al. (2009). "HER2 overexpression and amplification is present in a subset of ovarian mucinous carcinomas and can be targeted with trastuzumab therapy." *BMC Cancer* 9:433. (https://bmccancer.biomedcentral.com/articles/10.1186/1471-2407-9-433) and Chapel DB, et al. (2023) "HER2/ERBB2 IHC Expression and Copy Number Status in Ovarian Mucinous Tumors." PMID: 37406458.
- **Key finding — subtype association:** HER2 amplification in ~18.2% (6/33) mucinous carcinomas and 18.8% of borderline mucinous tumors; up to ~25% of mucinous ovarian tumors show HER2 IHC over-expression with strong IHC–FISH concordance (gastric/uterine-serous scoring). Mucinous is the classic HER2-amplified ovarian histotype.
- **Relevance:** High — mucinous HER2 quantification; model-selection relevance.
- **Source type:** Primary (IHC/FISH cohorts)

- **Citation:** [Case report] (2024). "Exceptional Response to Trastuzumab Deruxtecan in a Patient With Recurrent Ovarian Clear Cell Carcinoma With HER2 Expression." PMC11371105. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11371105/)
- **Key finding:** Clinical proof-of-concept: T-DXd response in HER2-expressing OCCC, linking the clear-cell HER2-high phenotype to ADC benefit.
- **Relevance:** Medium — anecdotal but mechanistically aligned with the OCCC/HER2 association.
- **Source type:** Case report

---

### Mesothelin (MSLN) — anetumab ravtansine (and MSLN landscape)

- **Citation:** Alqaisi HA, et al. (2025). "Randomized Phase II Study of Bevacizumab with Weekly Anetumab Ravtansine or Weekly Paclitaxel in Platinum-Resistant/Refractory High-Grade Ovarian Cancer (NCI Trial #10150)." *Clinical Cancer Research* 31(6):993–1001. DOI: 10.1158/1078-0432.CCR-24-3128 (https://aacrjournals.org/clincancerres/article/31/6/993/753252/) ; PMC11911801
- **Key finding:** 57 pts randomized (28 anetumab ravtansine+bevacizumab [ARB] vs 29 paclitaxel+bevacizumab [PB]); 88% MSLN-positive. Primary endpoint PFS FAVORED the control arm — paclitaxel/bevacizumab was superior (ORR up to 66%, median PFS 12.7 months). NEGATIVE result for the mesothelin ADC in this setting.
- **Relevance:** High — the key randomized readout; anetumab ravtansine did not beat standard therapy despite high MSLN expression (another "expression ≠ benefit" caution).
- **Source type:** Primary trial

- **Citation:** Santin AD, et al. (2022). "Safety and activity of anti-mesothelin ADC anetumab ravtansine in combination with pegylated-liposomal doxorubicin in platinum-resistant ovarian cancer: multicenter, phase Ib dose escalation and expansion study." *International Journal of Gynecological Cancer* 33(4):562–570. DOI: 10.1136/ijgc-2022-003927 (https://pmc.ncbi.nlm.nih.gov/articles/PMC10086500/)
- **Key finding:** Anetumab ravtansine (MTD 6.5 mg/kg q3w) + PLD in 65 PROC pts; overall ORR 27.7% (95% CI 17.3–40.2; 17 PR + 1 CR). In an exploratory subgroup with high mesothelin expression and ≤3 prior lines (n=19), ORR 42.1%, median PFS 8.5 months, DOR 8.3 months — an MSLN-high enrichment signal.
- **Relevance:** High — best MSLN-directed ADC efficacy signal and explicit MSLN-high → higher-ORR association.
- **Source type:** Primary trial

- **Citation:** [First-in-human phase I] (2020). "First-in-Human, Multicenter, Phase I Dose-Escalation and Expansion Study of Anti-Mesothelin Antibody-Drug Conjugate Anetumab Ravtansine in Advanced or Metastatic Solid Tumors." PMID: 32213105. (https://pubmed.ncbi.nlm.nih.gov/32213105/)
- **Key finding:** 148 pts in expansion cohorts (mesothelioma; ovarian, pancreatic, NSCLC, breast): 1 CR, 11 PR, 66 SD — established tolerability and MSLN-positive target population. (First author not captured from snippet — see Gaps.)
- **Relevance:** Medium — dose/target foundation for anetumab ravtansine.
- **Source type:** Primary trial

- **Citation:** Dum D, et al. (2021). "Mesothelin Expression in Human Tumors: A Tissue Microarray Study on 12,679 Tumors." PMC8067734. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8067734/) and review: "The Role of Mesothelin Expression in Serous Ovarian Carcinoma" *Cancers* 14(9):2283 (2022) (https://www.mdpi.com/2072-6694/14/9/2283)
- **Key finding — subtype association:** Serous ovarian carcinoma shows among the highest MSLN positivity of any tumor type (~97% positive in the TMA study; 90% by 5B2 and 94% by MN-1 antibody in other series), with strong membranous/apical staining. One HGSOC series: 55.1% of primary tumors positive; MSLN not prognostic. MSLN binds MUC16 (CA-125).
- **Relevance:** High — MSLN is a near-universal serous-ovarian antigen; quantitative histotype anchor.
- **Source type:** Primary (TMA) / Review

- **Note (MSLN broader landscape):** Search surfaced ongoing MSLN CAR-T and other MSLN ADC efforts in ovarian cancer but no positive pivotal ovarian trial was identified this session — see Gaps.

---

### DPEP3 (dipeptidase 3) — tamrintamab pamozirine (SC-003)

- **Citation:** Hamilton E, et al. (2020). "Tamrintamab pamozirine (SC-003) in patients with platinum-resistant/refractory ovarian cancer: Findings of a phase 1 study." *Gynecologic Oncology* 158(3):640–645. PMID: 32513564. DOI: 10.1016/j.ygyno.2020.05.038 (https://www.sciencedirect.com/science/article/abs/pii/S0090825820311239)
- **Key finding:** SC-003 = anti-DPEP3 ADC with a pyrrolobenzodiazepine (PBD) payload. First-in-human phase 1a/1b (NCT02539719), 74 platinum-resistant/refractory ovarian pts. Overall ORR only 4%; post hoc, higher DPEP3 expression correlated with better response. Program DISCONTINUED — lacked the safety/activity profile to warrant further development.
- **Relevance:** High — the definitive (negative) DPEP3-ADC clinical readout; expression-response correlation still supports biomarker selection.
- **Source type:** Primary trial

- **Citation:** [Preclinical] Abstract NT-113 (Stemcentrx): "SC-003, an ADC targeting Dipeptidase 3, exhibits potent anti-tumor activity in PDX models of high-grade serous ovarian cancer." (https://www.researchgate.net/publication/352306991) ; and review PMC12673879.
- **Key finding — subtype association:** DPEP3 is frequently over-expressed in high-grade serous ovarian carcinoma (HGSC); SC-003 induced tumor regression (incl. platinum-resistant PDX) preclinically.
- **Relevance:** Medium — DPEP3 = HGSC-enriched target; preclinical PDX precedent.
- **Source type:** Preclinical (abstract) / Review

---

### Why a subtype-resolved ADC-target expression atlas across models is useful (preclinical model selection & companion-Dx work)

- **Citation:** (2026). "Spatial, temporal, and molecular heterogeneity of ADC targets in high-grade serous ovarian carcinoma." *British Journal of Cancer*. DOI: 10.1038/s41416-026-03482-2 (https://www.nature.com/articles/s41416-026-03482-2)
- **Key finding:** ADC targets in HGSOC vary spatially, temporally, and molecularly; patient selection for ADC therapy depends on tumor target expression, so characterizing this heterogeneity is essential.
- **Relevance:** High — direct literature justification that mapping ADC-target expression (as our resource does) has clear translational value.
- **Source type:** Primary/analysis

- **Citation:** (2023). "Elucidating Novel Targets for Ovarian Cancer ADC Development: Integrating In Silico Prediction and Surface Plasmon Resonance to Identify Targets with Enhanced Antibody Internalization Capacity." PMC10594448. (https://pmc.ncbi.nlm.nih.gov/articles/PMC10594448/)
- **Key finding:** Systematic target-discovery workflow for ovarian-cancer ADCs; expression pre-screening in cell lines is a practical way to obtain broad expression profiles while sparing scarce patient samples (a library of ~19 ovarian + 2 normal lines used for mRNA screening).
- **Relevance:** High — explicitly endorses cell-line expression panels for ADC target profiling — the exact reuse case for our resource.
- **Source type:** Primary (methods/discovery)

- **Citation:** (2016). "Characterization of ovarian cancer cell lines as in vivo models for preclinical studies." *Gynecologic Oncology*. (https://www.sciencedirect.com/science/article/abs/pii/S0090825816307570)
- **Key finding:** Cell-line selection for preclinical work is typically based on patient history, histological subtype, mutation patterns, and signaling pathways — reinforcing that subtype-annotated model panels are the standard unit for model choice.
- **Relevance:** High — supports "pick a model by subtype + target expression" framing.
- **Source type:** Primary (model characterization)

- **Citation:** (2025). "Therapeutic Applications and Target Strategies of Antibody-Drug Conjugates in Ovarian Cancer." PMC12673879. (https://pmc.ncbi.nlm.nih.gov/articles/PMC12673879/)
- **Key finding:** Comprehensive review of ovarian-cancer ADC targets (FOLR1, NaPi2b, TROP2, HER2, MSLN, DPEP3, etc.); ADCs deliver cytotoxic payloads based on cancer-cell surface biomarkers, so target-expression profiling underpins the whole class.
- **Relevance:** High — one-stop review to frame the ADC-target landscape section.
- **Source type:** Review

---

## Search Log
- Query 1: "mirvetuximab soravtansine SORAYA trial ovarian cancer folate receptor alpha ORR FDA approval" — 9 examined, 6 relevant
- Query 2: "mirvetuximab soravtansine MIRASOL phase 3 trial overall survival progression-free FRα-high ovarian" — 8 examined, 6 relevant
- Query 3: "upifitamab rilsodotin NaPi2b UpRi UPLIFT ovarian cancer trial results SLC34A2" — 9 examined, 7 relevant
- Query 4: "trastuzumab deruxtecan DESTINY-PanTumor02 HER2 ovarian cancer ORR results" — 9 examined, 7 relevant
- Query 5: "sacituzumab govitecan TROP2 ovarian cancer trial ORR TACSTD2 platinum-resistant" — 7 examined, 5 relevant (surfaced sac-TMT as distinct agent)
- Query 6: "datopotamab deruxtecan Dato-DXd ovarian endometrial cancer TROP2 TROPION trial" — 8 examined, 6 relevant
- Query 7: "anetumab ravtansine mesothelin ovarian cancer phase trial MSLN results" — 8 examined, 6 relevant
- Query 8: "tamrintamab pamozirine SC-003 DPEP3 ovarian cancer antibody-drug conjugate phase 1" — 8 examined, 6 relevant
- Query 9: "folate receptor alpha FOLR1 expression ovarian cancer histotype ... immunohistochemistry" — 9 examined, 6 relevant
- Query 10: "HER2 ERBB2 expression mucinous ovarian carcinoma clear cell carcinoma amplification" — 8 examined, 7 relevant
- Query 11: "upifitamab rilsodotin UP-NEXT maintenance trial platinum-sensitive ovarian cancer NaPi2b" — 9 examined, 6 relevant
- Query 12: "mirvetuximab soravtansine FDA full approval March 2024 MIRASOL platinum-resistant ovarian" — 9 examined, 7 relevant
- Query 13: "mesothelin MSLN expression epithelial ovarian cancer serous histotype immunohistochemistry prevalence" — 7 examined, 5 relevant
- Query 14: "ovarian cancer cell line panel ADC target expression preclinical model selection biomarker companion diagnostic" — 9 examined, 6 relevant
- Query 15: "anetumab ravtansine bevacizumab paclitaxel randomized phase II ... 2025" — 9 examined, 4 relevant
- Query 16: "anetumab ravtansine ... Clinical Cancer Research 2025 first author NCI 10150" — 8 examined, 3 relevant (first author Alqaisi confirmed)
- Query 17: "TROPION-PanTumor03 datopotamab ... ESMO 2024 first author" — 8 examined, 5 relevant (first author Oaknin confirmed; hypothesized "Raghavan" was wrong)
- WebFetch: SORAYA final-OS PMC (author/ORR/DOR/OS/cutoff confirmed); PubMed MIRASOL NEJM (Moore KN 2023, 389(23):2162-2174, PMID 38055253); PubMed DESTINY-PanTumor02 (Meric-Bernstam 2024 JCO 42(1):47-58); PMC anetumab+PLD (Santin 2022 IJGC 33(4):562-570); PubMed SC-003 (Hamilton 2020 Gynecol Oncol 158(3):640-645)
- WebFetch failures (403, publisher paywalls): JCO SORAYA direct, NEJM MIRASOL direct, Targeted Oncology MIRASOL-OS, AACR CCR anetumab direct — worked around via PMC/PubMed/secondary coverage.

## Gaps
- **MIRASOL primary OS medians:** confirmed a statistically significant OS benefit (first in PROC) and the updated/full-approval dataset medians (16.85 vs 13.34 mo, median follow-up 30.5 mo, ~32% risk reduction). The exact NEJM-primary OS medians/HR were not isolated from a primary-source snippet (NEJM/Targeted Oncology returned 403). Verify against Moore 2023 NEJM full text.
- **Anetumab first-in-human phase I first author** (PMID 32213105) not captured from snippet; likely a Hassan-led CCR 2020 paper but NOT confirmed this session — do not cite an author without checking.
- **Sacituzumab tirumotecan (sac-TMT) ovarian ORR 40% / high-TROP2 61.5%** came from secondary coverage (OncLive/ASCO), not the primary ESMO abstract text — confirm primary before quoting.
- **UPLIFT exact citation** (volume/pages) is an SGO 2024 *Gynecologic Oncology* supplement abstract; page-level cite not pinned. First author Richardson DL inferred from OncLive interview attribution — confirm.
- **SORAYA primary Matulonis 2023 volume/pages** given as JCO 41(13) via DOI 10.1200/JCO.22.01900; exact page range not separately confirmed.
- **MSLN CAR-T / other MSLN ADCs in ovarian:** flagged as existing but not individually retrieved — a dedicated search would be needed if the descriptor wants that breadth.
- **DPEP3 expression quantification** is qualitative ("frequently over-expressed in HGSC") from review/abstract; no histotype-stratified prevalence table found.

## Breadth Flag
All 6 requested targets covered with primary trial + subtype-expression sources, plus a 4-paper block substantiating the "expression-atlas-for-model-selection" reuse argument. Notable landscape items NOT deeply covered (candidates for a follow-up cluster if desired): (1) MSLN CAR-T and non-anetumab MSLN ADCs in ovarian; (2) CDH6-directed ADC raludotatug deruxtecan (R-DXd, DESTINY-related) — an emerging ovarian ADC target not on the assigned list; (3) B7-H4, TROP2 companion-diagnostic assay development; (4) the primary sac-TMT ovarian abstract. Two cross-cutting cautions worth surfacing to synthesis: NaPi2b (UPLIFT), MSLN (NCI#10150), and DPEP3 (SC-003) all show that high target expression did NOT translate into a positive pivotal trial — the atlas is useful for model selection, but expression alone is a necessary-not-sufficient predictor of ADC benefit.
# Lit Review 06 — Genomic HRD Assessment & Tumour-Only WES Methods

**Cluster:** Q11 (genomic HRD quantification vs expression proxies) + Q12 (tumour-only somatic calling / CNVkit reference best practices)
**Compiled:** 2026-07-23
**Rule:** every claim below is grounded in a WebSearch/WebFetch result (URLs given). Papers reported individually; no cross-paper synthesis.

---

## Findings

### ═══ Q11: GENOMIC HRD QUANTIFICATION ═══

### HRD-LOH component (Abkevich 2012) — the LOH scar
- **Citation:** Abkevich V, Timms KM, Hennessy BT, et al. (2012). "Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer." *British Journal of Cancer* 107(10):1776–1782. DOI: 10.1038/bjc.2012.451. https://www.nature.com/articles/bjc2012451 (PubMed: https://pubmed.ncbi.nlm.nih.gov/23047548/)
- **Key finding:** Defined the HRD-LOH score as the *number of LOH regions of intermediate size* (longer than a set minimum but not spanning a whole chromosome) in a tumour. LOH regions of intermediate size were significantly more frequent in tumours with defective BRCA1, BRCA2 or RAD51C. The score was validated against BRCA deficiency in two independent ovarian datasets and was also computed across 57 cancer cell lines, identifying breast and pancreatic lines with BRCA defects. Authors state the score detects HR defects "regardless of aetiology or mechanism."
- **Relevance:** High — foundational LOH component of the genomic HRD score; explicitly demonstrated in cell lines, which is directly germane to our ~31 ovarian line panel.
- **Source type:** Primary / Method

### HRD-TAI / NtAI component (Birkbak 2012) — the telomeric allelic-imbalance scar
- **Citation:** Birkbak NJ, Wang ZC, Kim J-Y, et al. (2012). "Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents." *Cancer Discovery* 2(4):366–375. DOI: 10.1158/2159-8290.CD-11-0206. https://aacrjournals.org/cancerdiscovery/article/2/4/366 (PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3806629/)
- **Key finding:** Defined NtAI = number of sub-chromosomal regions with allelic imbalance that extend to the telomere (but do not cross the centromere). NtAI predicted cisplatin sensitivity in vitro and pathologic response to pre-operative cisplatin in triple-negative breast cancer (TNBC); in serous ovarian cancer, higher NtAI forecast better initial platinum response.
- **Relevance:** High — the TAI component of the composite HRD score.
- **Source type:** Primary / Method

### HRD-LST component (Popova 2012) — the large-scale-state-transition scar
- **Citation:** Popova T, Manié E, Rieunier G, et al. (2012). "Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2 inactivation." *Cancer Research* 72(21):5454–5462. DOI: 10.1158/0008-5472.CAN-12-1470. https://aacrjournals.org/cancerres/article/72/21/5454 (PubMed: https://pubmed.ncbi.nlm.nih.gov/22933060/)
- **Key finding:** Defined LST as a chromosomal break between adjacent regions of at least 10 Mb (after smoothing/filtering short segments), with the threshold ploidy-adjusted. LST was a robust indicator of BRCA1 status, predicting BRCA1/2 inactivation in basal-like breast carcinomas with 100% sensitivity and 90% specificity (97% accuracy).
- **Relevance:** High — the LST component of the composite HRD score.
- **Source type:** Primary / Method

### Combined HRD score / Genomic Instability Score and the GIS ≥42 cutoff (Telli 2016)
- **Citation:** Telli ML, Timms KM, Reid J, et al. (2016). "Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy in Patients with Triple-Negative Breast Cancer." *Clinical Cancer Research* 22(15):3764–3773. DOI: 10.1158/1078-0432.CCR-15-2477. https://aacrjournals.org/clincancerres/article/22/15/3764
- **Key finding:** The HRD score is an **unweighted sum of the LOH + TAI + LST component scores**. HR deficiency was defined as **HRD score ≥42 OR a BRCA1/2 mutation** — this is the origin of the widely used "GIS ≥42" cutoff. In neoadjuvant platinum-based therapy for TNBC, HR-deficient tumours achieved 53% pathologic complete response (pCR) vs 18% in non-deficient (adjusted OR 4.64); among BRCA-negative patients, a high HRD score alone gave 50% pCR (≈4.5× odds vs non-deficient). Cutoff confirmed independently via ASCO Post coverage (https://ascopost.com/issues/january-25-2016/...).
- **Relevance:** High — defines the composite genomic HRD score and its canonical threshold; the reference point for "true genomic HRD."
- **Source type:** Primary / Method

### Clinical/regulatory genomic HRD test — Myriad myChoice CDx (GIS = genomic scars, NOT expression)
- **Citation:** Myriad Genetics. "MyChoice CDx" companion diagnostic (FDA-approved Oct 2019 as CDx for niraparib/Zejula in ovarian cancer). https://myriad.com/genetic-tests/mychoicecdx-tumor-test/ ; FDA-approval press release https://investor.myriad.com/news-releases/news-release-detail/27106/
- **Key finding:** The only FDA-approved HRD companion diagnostic determines HRD status from (a) BRCA1/2 sequencing including large rearrangements and (b) a tumour **Genomic Instability Score (GIS) composed of LOH + TAI + LST** — i.e., genomic scars, explicitly **not gene expression**. Vendor materials note the GIS cutoff has been clinically validated only in ovarian cancer.
- **Relevance:** High — establishes that the clinically/regulatory-accepted definition of "genomic HRD" is scar-based (LOH/TAI/LST) ± BRCA status, providing the benchmark against which an expression "HRD score" must NOT be equated.
- **Source type:** Guideline (regulatory / CDx)

### Mutational Signature 3 (Alexandrov 2013) — the HRD SBS signature (original)
- **Citation:** Alexandrov LB, Nik-Zainal S, Wedge DC, et al. (2013). "Signatures of mutational processes in human cancer." *Nature* 500(7463):415–421. DOI: 10.1038/nature12477. https://www.nature.com/articles/nature12477
- **Key finding:** Identified Signature 3, associated with failure of DNA double-strand-break repair by homologous recombination. Signature 3 is strongly associated with germline and somatic BRCA1/BRCA2 mutations in breast, pancreatic and ovarian cancers, and co-occurs with elevated numbers of large (>3 bp) insertions/deletions with overlapping microhomology at breakpoint junctions.
- **Relevance:** High — the mutational-signature readout of HRD; distinct from the copy-number scar score.
- **Source type:** Primary / Method

### SBS3 in the mutational-signature repertoire (Alexandrov 2020, PCAWG)
- **Citation:** Alexandrov LB, Kim J, Haradhvala NJ, et al. (2020). "The repertoire of mutational signatures in human cancer." *Nature* 578(7793):94–101. DOI: 10.1038/s41586-020-1943-3. https://www.nature.com/articles/s41586-020-1943-3 (PubMed: https://pubmed.ncbi.nlm.nih.gov/32025018/)
- **Key finding:** Reference re-derivation of signatures from 84,729,690 somatic mutations across 4,645 whole genomes and 19,184 exomes, yielding 49 SBS, 11 DBS, 4 clustered and 17 indel signatures. Signature 3 was refined and renamed **SBS3**, the canonical HRD single-base-substitution signature (COSMIC v3).
- **Relevance:** High — current nomenclature/definition of the HRD SBS signature (SBS3). NB: SBS3 is a "flat," featureless profile and is the hardest signature to call reliably from limited data such as exomes (see Gaps).
- **Source type:** Primary / Reference

### scarHRD (Sztupinszki 2018) — deriving the scar-based HRD score from WES/WGS
- **Citation:** Sztupinszki Z, Diossy M, Krzystanek M, et al. (2018). "Migrating the SNP array-based homologous recombination deficiency measures to next generation sequencing data of breast cancer." *npj Breast Cancer* 4:16. DOI: 10.1038/s41523-018-0066-6. https://www.nature.com/articles/s41523-018-0066-6 (tool: https://github.com/sztup/scarHRD)
- **Key finding (WebFetch-verified):** Introduces the **scarHRD R package** to compute the LOH, TAI(NtAI) and LST scores — and their sum — from NGS (WES/WGS) rather than SNP arrays. SNP-array vs WES correlations: **NtAI r=0.84, LST r=0.79, HRD-LOH r=0.73, HRD-sum r=0.87**. Reducing coverage to **30×** did not degrade the correlation; however the number of LOH events was somewhat lower in WES-based estimates (segmentation/quality effects). *(Caution: the r=0.82 correlation to Myriad myChoice that appears in some secondary summaries is NOT reported in this paper — it comes from other work; do not attribute it to Sztupinszki 2018.)*
- **Relevance:** High — this is the appropriate, WES-compatible tool for computing a genuine genomic HRD (scar) score from our tumour-only exomes; validated down to 30× WES.
- **Source type:** Primary / Method

### HRDetect (Davies 2017) — WGS mutational-signature predictor of BRCA1/2 deficiency
- **Citation:** Davies H, Glodzik D, Morganella S, et al. (2017). "HRDetect is a predictor of BRCA1 and BRCA2 deficiency based on mutational signatures." *Nature Medicine* 23(4):517–525. DOI: 10.1038/nm.4292. https://www.nature.com/articles/nm.4292 (PubMed: https://pubmed.ncbi.nlm.nih.gov/28288110/)
- **Key finding:** A lasso logistic-regression model integrating six whole-genome mutational features (including SBS3, rearrangement signatures, and the deletion-with-microhomology / HRD index). Detected BRCA1/2-deficient tumours with **98.7% sensitivity (AUC 0.98)** in 560 breast cancers, and revealed BRCA1/2 deficiency in up to ~22% of cases (vs ~1–5% by mutation alone).
- **Relevance:** Medium — gold-standard HRD signature classifier, but **requires whole-genome sequencing** (uses structural-rearrangement signatures); not directly applicable to WES-only data.
- **Source type:** Primary / Method

### CHORD (Nguyen 2020) — WGS random-forest HRD classifier
- **Citation:** Nguyen L, Martens JWM, Van Hoeck A, Cuppen E (2020). "Pan-cancer landscape of homologous recombination deficiency." *Nature Communications* 11:5584. DOI: 10.1038/s41467-020-19406-4. https://www.nature.com/articles/s41467-020-19406-4 (tool: https://github.com/UMCUGenetics/CHORD)
- **Key finding:** A random-forest **C**lassifier **o**f **H**omologous **R**ecombination **D**eficiency using relative counts of somatic mutation contexts — chiefly **deletions with flanking microhomology and 1–100 kb structural duplications**. Identifies BRCA1/BRCA2/RAD51C/PALB2 inactivation as the most frequent genetic cause of HRD pan-cancer (primary and metastatic).
- **Relevance:** Medium — like HRDetect, **relies on structural-variant contexts and therefore on WGS**; not designed for WES.
- **Source type:** Primary / Method

### HRD in cell-line libraries — retained genomic scars but do NOT track drug response (Takamatsu 2024, *Scientific Data*)
- **Citation:** Takamatsu S, Murakami K, Matsumura N (2024). "Homologous Recombination Deficiency Unrelated to Platinum and PARP Inhibitor Response in Cell Line Libraries." *Scientific Data* 11:171. DOI: 10.1038/s41597-024-03018-4. https://www.nature.com/articles/s41597-024-03018-4 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10847511/)
- **Key finding (WebFetch-verified):** Computed HRD three ways across **1,182 CCLE cell lines** (incl. **62 ovarian** and 54 breast) with COSMIC Cell Line Project validation: (a) **HRD score = TAI+LST+LOH via scarHRD** from SNP-array copy number; (b) **SBS3 via SigMA**; (c) BRCA1/2 mutation/methylation + HR-gene LOH. Genomic HRD markers **are retained and detectable in cell lines** — HRD scores and SBS3 were significantly higher in lines with HR-gene alterations. However, higher HRD scores and SBS3 were **significantly associated with resistance (not sensitivity)** to platinum/PARP inhibitors in vitro, and BRCA1/2 alteration did not correlate with sensitivity.
- **Relevance:** High — directly on-point and in our target journal. Confirms genomic HRD scars can be measured from ovarian cell-line genomic data, but is a strong cautionary tale against interpreting cell-line HRD as a functional/therapeutic-response biomarker.
- **Source type:** Primary / Method (Data Descriptor)

### Genomic HRD scars are a static historical record (retained despite HR restoration)
- **Citation:** Concept documented in HRD-signature literature; e.g., ovarian/uterine carcinosarcoma HRD-signature study, *Gynecologic Oncology* abstract (2022). https://www.sciencedirect.com/science/article/abs/pii/S0090825822012628 (see also Takamatsu 2024 above)
- **Key finding:** HRD generates specific genomic scars (SBS3/SBS8, indel signature ID6, structural-variant signatures, and the copy-number LOH/TAI/LST scars). These "scars remain part of the genome despite evolving drug resistance due to HR restoration" — i.e., the scar-based HRD score reflects a *past* HR-deficient state and can persist even after a cell line regains HR proficiency.
- **Relevance:** High — explains why a genomic HRD score is measurable in a long-passaged ovarian cell line, while also cautioning that it need not reflect *current* functional HR status (consistent with Takamatsu 2024).
- **Source type:** Primary (abstract) / concept

### Expression-based "HRD signature" (Peng 2014) — a FUNCTIONAL/transcriptional readout, not a genomic HRD scar
- **Citation:** Peng G, Chun-Jen Lin C, Mo W, et al. (2014). "Genome-wide transcriptome profiling of homologous recombination DNA repair." *Nature Communications* 5:3361. DOI: 10.1038/ncomms4361. https://www.nature.com/articles/ncomms4361 (PubMed: https://pubmed.ncbi.nlm.nih.gov/24553445/)
- **Key finding:** Derived a **230-gene HRD gene-expression signature** from a non-malignant mammary epithelial cell line in which HR was perturbed. The authors explicitly frame it as a way "to functionally assess HR repair status **without interrogating individual genetic alterations**" and show it predicts clinical outcome and PARP-inhibitor sensitivity across lineages (tested in NCI60 and 51 breast lines). This is almost certainly the "Peng HRD-related expression" signature referenced in the prior analysis.
- **Relevance:** High — this is a **transcriptomic/functional** signature, categorically different from the genomic-scar HRD score (LOH/TAI/LST) or SBS3. It measures an expression state, not accumulated genomic damage.
- **Source type:** Primary / Method

### Expression HRD signatures are dataset-dependent and not the accepted clinical HRD measure
- **Citation:** (illustrative) "Gene expression signature for predicting homologous recombination deficiency in triple-negative breast cancer." *npj Breast Cancer* (2024) 10:... DOI: 10.1038/s41523-024-00671-1. https://www.nature.com/articles/s41523-024-00671-1
- **Key finding:** Multiple published transcriptomic HRD signatures show **minimal gene overlap** with one another and with Peng's 230-gene set (a search summary reported as little as 1 gene shared between a 217-gene panel and Peng's 230 genes — see Gaps re: exact figure). No expression signature is an FDA-approved or guideline-endorsed HRD test; the accepted genomic HRD definition is scar-based (GIS = LOH+TAI+LST) ± BRCA status and/or SBS3 (see Myriad myChoice, Telli 2016, Alexandrov entries).
- **Relevance:** High — substantiates that expression signatures are **not** an accepted substitute for genomic HRD, and that labelling an expression score as "genomic HRD" is a category error.
- **Source type:** Primary (comparison) / supporting

---

### ═══ Q12: TUMOUR-ONLY WES / NO MATCHED NORMAL ═══

### GATK Mutect2 tumour-only caveats and best-practice recommendation for matched normals
- **Citation:** Broad Institute, Genome Analysis Toolkit (GATK) — "Mutect2" and "FAQ for Mutect2" documentation. https://gatk.broadinstitute.org/hc/en-us/articles/360037261691-Mutect2-BETA ; https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2
- **Key finding:** "Tumor-normal mode is much more reliable than tumor-only mode because tumor-only mode finds a lot of false positives." Any genome carries tens of thousands of germline variants, and roughly **~30,000 are too rare to appear in gnomAD**, so they are "hardly [distinguishable] from somatic variants" even with strict germline filtering. In tumour-normal mode Mutect2 uses population allele frequency (germline resource), the matched-normal reads, and tumour allele fraction to call germline; in **tumour-only mode the normal-read evidence is simply missing**. The panel-of-normals (PoN) filters recurrent technical artifact sites and common germline sites; gnomAD supplies the population-AF prior for germline filtering. Default germline-probability thresholds differ by mode (tumour-only 5e-8 vs tumour-normal 1e-6).
- **Relevance:** High — the authoritative statement that (i) matched normal is the recommended design, and (ii) even PoN + gnomAD filtering cannot remove rare/private germline variants — the exact mechanism behind spuriously high ATM/ATR/BRCA2 "mutation" rates in tumour-only cell-line calls.
- **Source type:** Guideline / Method (tool documentation)

### Quantified false-positive inflation in tumour-only calling + ancestry bias (Halperin 2017)
- **Citation:** Halperin RF, Carpten JD, Manojlovic Z, et al. (2017). "A method to reduce ancestry related germline false positives in tumor only somatic variant calling." *BMC Medical Genomics* 10:61. DOI: 10.1186/s12920-017-0296-8. https://pmc.ncbi.nlm.nih.gov/articles/PMC5649057/
- **Key finding (WebFetch-verified):** After strict filtering against dbSNP/ExAC/ESP6500/ARIC, an **average of 224 private germline variants per sample remained (range 126–319)**. Filtering alone gave a positive predictive value of only **35–62% for European-American** samples and **20–40% for non-European** samples — i.e., roughly a **~50% false-discovery rate for European ancestry and >70% FDR for African ancestry**. Their Bayesian caller **LumosVar** (models allelic copy number and clonal sample fractions using differences in germline vs somatic allele frequency in impure tumours, with 1000 Genomes + COSMIC priors) improved PPV to **67–91%**.
- **Relevance:** High — concrete FP/FDR magnitudes for tumour-only calling and the ancestry dependence (public germline databases under-represent non-European ancestry). Directly quantifies the risk in our tumour-only WES.
- **Source type:** Primary / Method

### UNMASC (Little 2021) — tumour-only calling with UNMATCHED (pooled/public) normal controls
- **Citation:** Little P, Jo H, Hoyle A, et al. (2021). "UNMASC: tumor-only variant calling with unmatched normal controls." *NAR Cancer* 3(4):zcab040. DOI: 10.1093/narcan/zcab040. https://pmc.ncbi.nlm.nih.gov/articles/PMC8494212/
- **Key finding (WebFetch-verified):** With ~10 unmatched normal controls, UNMASC achieved **94% sensitivity, 99% specificity, 76% PPV**, outperforming SomVarIUS/LumosVar 2.0 etc. Unmatched normals help because they capture **sequencing/alignment artifacts and population variants that stem from chemistry/platform rather than individual biology** (identifying hard-to-map regions via VAF clustering), which can then be applied to filter tumour variants. **Key limitation:** somatic vs germline becomes "indistinguishable in highly pure tumours where founder mutations present at high frequencies mimic germline heterozygosity."
- **Relevance:** High — directly analogous to our design (CNVkit/variant filtering against a pooled panel of public normals). The purity caveat is especially important: **cell lines are effectively ~100% tumour purity**, precisely the regime where high-VAF germline SNPs are hardest to separate from clonal somatic calls — reinforcing why tumour-only cell-line calls over-report germline as somatic.
- **Source type:** Primary / Method

### CNVkit primary paper + reference-construction best practices (Talevich 2016)
- **Citation:** Talevich E, Shain AH, Botton T, Bastian BC (2016). "CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing." *PLoS Computational Biology* 12(4):e1004873. DOI: 10.1371/journal.pcbi.1004873. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873 (docs: https://cnvkit.readthedocs.io/)
- **Key finding:** CNVkit builds a **pooled reference of per-bin copy-number estimates from several normal samples**, using both on-target and off-target reads. Best-practice guidance (paper + official docs): (i) a reference must be built **specifically for each target-capture panel** (from the panel's BED of baited regions) and from samples **sequenced on the same platform/capture kit**; (ii) ideally match sample type (FFPE vs fresh); (iii) a **flat reference** (neutral copy number) is possible when no normals are available but is inferior because it cannot correct capture/GC bias; (iv) if no normals exist, a reference may be built from tumours chosen for minimal, non-recurrent CNVs.
- **Relevance:** High — governs whether our "pooled panel of 5 public normal exomes" is methodologically sound. The **same-capture-kit / same-platform requirement is the critical caveat**: public normals captured with a different kit than our tumours will introduce systematic coverage bias that CNVkit can misread as copy-number change.
- **Source type:** Primary / Method + tool documentation

### Cell lines inherently lack matched normals — the CCLE precedent (Ghandi 2019)
- **Citation:** Ghandi M, Huang FW, Jané-Valbuena J, et al. (2019). "Next-generation characterization of the Cancer Cell Line Encyclopedia." *Nature* 569(7757):503–508. DOI: 10.1038/s41586-019-1186-3. https://www.nature.com/articles/s41586-019-1186-3
- **Key finding:** CCLE applied a **harmonized variant-calling pipeline** across WES/WGS/RNA-seq/targeted data. Because cell lines have no patient-matched normal, germline is handled by **filtering against population databases**; the study reports **high concordance of germline variant calls between CCLE and GDSC** for the same lines (with discordances flagging mislabelling/genetic drift). (The precise germline-filtering panel/thresholds were not fully specified in the retrieved text — see Gaps.)
- **Relevance:** High — the field-standard reference showing how large cell-line resources handle the unavoidable absence of matched normals: aggressive population-frequency filtering, accepting residual private-germline contamination.
- **Source type:** Primary

### Reference-sample choice materially affects tumour-only somatic calling (myeloid benchmark)
- **Citation:** (2025) "Evaluation of reference sample type for somatic variant calling in myeloid cancers." *Annals of Hematology*. https://link.springer.com/article/10.1007/s00277-025-06699-y
- **Key finding:** Benchmarked reference-sample types for somatic calling; **T-cells gave the highest sensitivity (0.91–1.00) and lowest false-positive rate**, while the **tumour-only pipeline produced a significantly higher average number of false positives** than any matched-normal reference. Reinforces that a true matched normal outperforms tumour-only regardless of downstream filtering.
- **Relevance:** Medium — independent, recent confirmation of tumour-only FP inflation vs matched-normal designs.
- **Source type:** Primary / benchmark

### Tumour-only filtering strategies for targeted/clinical panels (Sukhai / practical filtering)
- **Citation:** "Somatic Tumor Variant Filtration Strategies to Optimize Tumor-Only Molecular Profiling Using Targeted Next-Generation Sequencing Panels." *Journal of Molecular Diagnostics* (ScienceDirect). https://www.sciencedirect.com/science/article/pii/S1525157817305986
- **Key finding:** Describes in-silico germline-subtraction strategy for tumour-only panels: filter against population databases (dbSNP, ExAC/gnomAD, 1000 Genomes), retain/annotate against somatic databases (COSMIC), and apply VAF-based heuristics — while cautioning that dbSNP contains some pathogenic germline variants and COSMIC contains some germline/artifact variants (imperfect databases). 
- **Relevance:** Medium — describes the practical filtering cascade one applies when no matched normal exists (relevant to re-filtering our calls); notes database imperfection as a residual-error source.
- **Source type:** Primary / Method

---

## Search Log
- Query 1: "Abkevich 2012 loss of heterozygosity LOH score homologous recombination deficiency ovarian cancer BRCA" — 8 examined, 3 relevant (Abkevich BJC, PMID 23047548)
- Query 2: "Telli 2016 Clinical Cancer Research HRD score LOH TAI LST genomic instability score 42 cutoff" — 8 examined, 2 relevant (Telli CCR 2016; ASCO Post)
- Query 3: "Birkbak 2012 telomeric allelic imbalance TAI … Popova 2012 large-scale state transitions LST" — 16 examined, 4 relevant (Birkbak Cancer Discovery; Popova Cancer Research)
- Query 4: "Alexandrov 2013 Nature mutational signatures signature 3 … Alexandrov 2020 SBS3 repertoire" — 20 examined, 4 relevant (nature12477; s41586-020-1943-3)
- Query 5: "Davies 2017 HRDetect Nature Medicine …" — 10 examined, 2 relevant (nm.4292)
- Query 6: "scarHRD Sztupinszki 2018 WES … CHORD Nguyen 2020 Nature Communications" — 18 examined, 4 relevant (scarHRD npj BC; CHORD Nat Commun)
- Query 7: "HRD score ovarian cancer cell lines whole exome sequencing detectable BRCA" — 7 examined, 3 relevant (multiple ovarian HRD cutoff papers)
- Query 8: "Peng 2014 expression-based HRD signature …" — 6 examined, 2 relevant (Peng Nat Commun; recent overlap-comparison paper)
- Query 9: "Mutect2 tumor-only mode germline contamination false positives PoN GATK matched normal" — 9 examined, 5 relevant (GATK docs)
- Query 10: "CNVkit reference construction PoN vs flat reference same capture kit Talevich" — 7 examined, 3 relevant (CNVkit docs)
- Query 11: "tumor-only somatic variant calling false positive rate quantified concordance matched normal gnomAD" — 9 examined, 4 relevant (Halperin; UNMASC; myeloid benchmark)
- Query 12: "cell line WES no matched normal germline filtering common SNP CCLE somatic calling" — 10 examined, 3 relevant (JMD filtration; CCLE)
- Query 13: "Talevich 2016 CNVkit PLoS Comput Biol …" — 10 examined, 2 relevant (Talevich e1004873)
- Query 14: "Ghandi 2019 CCLE mutation calling without matched normal …" — 10 examined, 2 relevant (Ghandi Nature 2019)
- Query 15: "HRD genomic scars retained cell lines … SBS3 ovarian cell lines" — 9 examined, 4 relevant (Takamatsu Sci Data 2024; carcinosarcoma HRD abstract)
- Query 16: "Peng 2014 Nat Commun 230 gene HRD signature expression PARP" — 9 examined, 2 relevant (Peng ncomms4361, PMID 24553445)
- Query 17: "FDA approved HRD test GIS Myriad myChoice CDx genomic scar not expression ovarian" — 7 examined, 3 relevant (myChoice CDx materials)
- **WebFetch confirmations:** Halperin 2017 (FDR/PPV figures, LumosVar); scarHRD (correlations, 30× coverage); Telli 2016 GIS ≥42 (ASCO Post); Takamatsu 2024 (62 ovarian lines, scarHRD+SigMA, resistance association); Talevich CNVkit; UNMASC (94/99/76%, purity limitation). Two failed fetches: AACR Telli page (403; cutoff confirmed via ASCO Post + two search summaries) and a wrong-PMID guess (27140933 → unrelated melanoma paper; discarded).

## Gaps
1. **SBS3 detectability from WES specifically.** SBS3 is a flat, low-information signature that is well documented as difficult to call from exomes/low mutation counts; Takamatsu used **SigMA** (a tool built for panel/exome SBS3) but I did not independently pull the primary SigMA benchmark (Gulhan et al. 2019, Nat Genet) — recommend verifying SigMA's WES/panel sensitivity if SBS3 is to be computed from our exomes.
2. **Exact gene-overlap figure between expression signatures.** The "1 of 217 genes shared with Peng's 230-gene set" figure came from a WebSearch summary of a 2024 npj Breast Cancer paper; the specific number should be confirmed against that article before quoting.
3. **CCLE germline-filtering specifics.** Ghandi 2019's exact germline-filtering databases/thresholds were not fully captured in the retrieved text; the CCLE methods/supplement would confirm the precise cascade.
4. **Telli 2016 DOI** stated from the standard AACR pattern (10.1158/1078-0432.CCR-15-2477); the AACR landing page returned 403, but volume/issue/pages (22(15):3764–3773) and the ≥42 cutoff are triple-confirmed via URL path + ASCO Post + two search summaries.
5. **A single canonical "tumour-only false-positive rate."** No universal figure exists — it is strongly ancestry- and pipeline-dependent (Halperin: ~50–70% FDR by filtering alone). A secondary summary mentioned a "67% false positive rate" but I did not locate/verify its primary source, so it is not reported as a finding.
6. **Minimum number of normals for a CNVkit reference.** Docs emphasize same-capture-kit/platform over a strict count; whether 5 public normals suffices (and the bias from a mismatched capture kit) is a judgement the docs support qualitatively but do not give a hard n.

## Breadth Flag
Coverage is strong across both questions. All named seed papers were located and verified (Abkevich, Birkbak, Popova, Telli, Alexandrov 2013/2020, scarHRD, HRDetect, CHORD, Peng), plus a highly relevant target-journal cell-line HRD Data Descriptor (Takamatsu 2024, *Sci Data*) and a directly analogous unmatched-normal caller (UNMASC). The most decision-relevant methodological takeaways for our WES resource:
- **For genomic HRD from WES:** use **scarHRD** (LOH+TAI+LST; validated to 30× WES). HRDetect and CHORD are **WGS-only** and not applicable. SBS3 is computable from WES but only with SBS3-specialized tools and adequate mutation counts.
- **Expression signatures (Peng) are NOT genomic HRD** and are not clinically accepted as HRD measures — relabelling required.
- **Tumour-only calling** inflates false positives via unremovable private/rare germline variants (Halperin: 224/sample, ~50–70% FDR); the effect is **worst in ~100%-pure cell lines** (UNMASC purity caveat) — the mechanistic explanation for the implausible ATM/ATR/BRCA2 rates.
- **CNVkit** pooled reference is standard, but **must be same capture kit/platform** as the tumours — a mismatched public-normal kit is a real bias risk to flag.
