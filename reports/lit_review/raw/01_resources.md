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
