# 07 — Refine: Novelty Confirmation + MMMT Gap

Targeted literature refinement for the *Scientific Data* Data Descriptor of the OvCAN multi-omic
(bulk RNA-seq + TMT proteomics + WES) resource of ~40 human ovarian cancer cell lines.

Three novelty-confirmation questions. Every claim below is grounded in a WebSearch/WebFetch result
retrieved 2026-07-23; no claims from memory. Papers reported individually. Calibration is confined
to the Gaps and Breadth Flag sections.

---

## Findings

### Q1 — Is there any published integrated MULTI-OMIC resource specifically for an ovarian cancer cell-line PANEL (≥2 of transcriptomics / proteomics / genomics)?

**Direct answer:** No pan-subtype, ~40-line ovarian cell-line resource combining transcriptomics +
proteomics + genomics exists. But the novelty claim must be framed carefully: (a) one study
(Shrestha 2021) already combined all three omics on an ovarian cell-line panel, restricted to a
single subtype (LGSOC, n=14); (b) several ovarian-specific 2-omic panels exist (Coscia 2016
proteomics-primary; Beaufort 2014 transcriptomics + targeted genomics; Domcke 2013 genomics;
Ince 2015 genomics + pharmacology). OvCAN's novelty is BREADTH (all major + rare subtypes together),
SCALE (~40 lines), and the specific WES + RNA + TMT-proteomics combination on one uniform panel — not
"first multi-omic ovarian cell-line data." Closest prior art to cite explicitly: Shrestha 2021 and
Coscia 2016.

#### Shrestha et al. 2021 — LGSOC cell-line multi-omic resource (CLOSEST PRIOR ART)
- **Citation:** Shrestha R, Llaurado Fernandez M, Dawson A, et al. Multiomics Characterization of Low-Grade Serous Ovarian Carcinoma Identifies Potential Biomarkers of MEK Inhibitor Sensitivity and Therapeutic Vulnerability. *Cancer Res.* 2021;81(7):1681–1694. DOI: 10.1158/0008-5472.CAN-20-2222. PMID: 33441310. Preprint: bioRxiv 2020.06.18.135061. Data: EGA EGAS00001004724; PRIDE PXD019544.
- **Key finding:** Integrated **whole-exome sequencing + RNA-seq + mass-spectrometry proteomics** on **14 low-grade serous ovarian carcinoma (LGSOC) patient-derived cell lines**. This IS a genuine ovarian cell-line multi-omic resource combining all three omics layers — but restricted to a single rare subtype and 14 lines. Proteomics labelling strategy (label-free vs TMT) not stated on the abstract page (search-limitation).
- **Relevance:** High — the single most important pre-emption to address. Note personnel overlap: co-author Manuel Llaurado Fernandez is also on the 2025 gynecologic-carcinosarcoma-models paper (Q2), i.e. the BC Cancer / OVCARE group that supplies the VOA lines. OvCAN likely extends this group's own prior LGSOC work to a pan-subtype panel.
- **Source type:** Primary

#### Coscia et al. 2016 — proteomics of 26 ovarian cell lines (KNOWN)
- **Citation:** Coscia F, Watters KM, Curtis M, et al. Integrative proteomic profiling of ovarian cancer cell lines reveals precursor cell associated proteins and functional status. *Nat Commun.* 2016;7:12645. DOI: 10.1038/ncomms12645. PMID: 27561551. Data: PRIDE PXD003668.
- **Key finding:** Single-run label-free MS proteomics of **26 ovarian cancer cell lines** (+ HGSOC tumours, immortalized OSE and FTE), >10,000 proteins quantified; three proteomic groups (epithelial / clear-cell / mesenchymal); 67-protein signature. Proteomics-primary; integrates *public* genomic/transcriptomic data rather than generating matched WES/RNA on the same panel.
- **Relevance:** High — the other closest prior art. Subtype span is limited (common HGSC/epithelial lines + a clear-cell group; no carcinosarcoma, no SCCOHT, no mucinous/endometrioid breadth), label-free (not TMT), no matched WES panel.
- **Source type:** Primary

#### Beaufort et al. 2014 — Ovarian Cancer Cell Line Panel (OCCP)
- **Citation:** Beaufort CM, Helmijr JCA, Piskorz AM, et al. Ovarian Cancer Cell Line Panel (OCCP): Clinical Importance of In Vitro Morphological Subtypes. *PLoS One.* 2014;9(9):e103988. DOI: 10.1371/journal.pone.0103988. PMC4167545.
- **Key finding:** **39 ovarian cancer cell lines** profiled under uniform conditions for **mRNA + microRNA expression + exon sequencing + drug response** (multi-modal: transcriptomics + targeted genomics + pharmacology). 14 HGS, 4 serous-type, 1 LGS, 20 non-serous. **No proteomics.**
- **Relevance:** High — establishes that a multi-modal ovarian cell-line panel resource exists; OvCAN's differentiators are proteomics (TMT), full WES (vs exon panel), and rare-subtype breadth.
- **Source type:** Primary / Resource

#### Domcke et al. 2013 — genomic benchmarking of 47 ovarian lines (LANDMARK)
- **Citation:** Domcke S, Sinha R, Levine DA, Sander C, Schultz N. Evaluating cell lines as tumour models by comparison of genomic profiles. *Nat Commun.* 2013;4:2126. DOI: 10.1038/ncomms3126. PMID: 23839242.
- **Key finding:** Compared **47 ovarian cancer cell lines** to TCGA HGSOC by copy-number, mutation and mRNA profiles; showed the most-used lines are poor HGSOC models and nominated rarely used lines as better. Genomics-focused; reuses existing data (not a new multi-omic generation).
- **Relevance:** High — canonical citation for "cell-line models poorly represent ovarian tumours"; underpins the rationale for a better characterized panel.
- **Source type:** Primary

#### Ince et al. 2015 — 25 new ovarian cell lines spanning subtypes (OCI collection)
- **Citation:** Ince TA, Sousa AD, Jones MA, et al. Characterization of twenty-five ovarian tumour cell lines that phenocopy primary tumours. *Nat Commun.* 2015;6:7419. DOI: 10.1038/ncomms8419. PMID: 26080861. PMC4473807.
- **Key finding:** New medium enabling establishment of **25 ovarian cell lines across diverse subtypes** that retain genomic landscape, histopathology, molecular features and drug response of parental tumours. Genomics + pharmacology; not a proteomic resource.
- **Relevance:** Medium — another multi-subtype cell-line panel resource; comparator for breadth but not for proteomics/WES integration.
- **Source type:** Primary / Resource

#### Thu et al. 2017 — 18 HGSC lines + 3 PDX, multi-layer genomics
- **Citation:** Thu KL, Papari-Zareei M, Stastny V, et al. (Gazdar AF). A comprehensively characterized cell line panel highly representative of clinical ovarian high-grade serous carcinomas. *Oncotarget.* 2017 (article 9929); reported as 8(31):50489–50499. URL: https://www.oncotarget.com/article/9929/text/ (exact volume/pages to confirm — search-limitation).
- **Key finding:** **18 HGSC cell lines + 3 matched PDX** with WES + copy-number + DNA methylation + gene expression + spectral karyotyping, benchmarked to parental tumours. HGSC-only; **no proteomics.**
- **Relevance:** Medium — multi-omic (genomic + epigenomic + transcriptomic) but single-subtype and no proteomics.
- **Source type:** Primary

#### Nusinow et al. 2020 — CCLE proteomics (PAN-CANCER, includes ovarian)
- **Citation:** Nusinow DP, Szpyt J, Ghandi M, et al. Quantitative Proteomics of the Cancer Cell Line Encyclopedia. *Cell.* 2020;180(2):387–402.e16. DOI: 10.1016/j.cell.2019.12.023.
- **Key finding:** **TMT10-plex** proteomics of **375 CCLE cell lines across 22 lineages** (>12,000 proteins). Pan-cancer; ovarian lines included only as part of CCLE, not an ovarian-specific resource.
- **Relevance:** Medium — confirms the TMT-proteomics-of-cell-lines approach exists pan-cancer; OvCAN differs by being ovarian-focused, rare-subtype-inclusive, with matched WES + RNA on the same panel.
- **Source type:** Resource

#### Gonçalves et al. 2022 — ProCan pan-cancer proteomics (includes ovarian)
- **Citation:** Gonçalves E, Poulos RC, Cai Z, et al. Pan-cancer proteomic map of 949 human cell lines. *Cancer Cell.* 2022;40(8):835–849.e8. DOI: 10.1016/j.ccell.2022.06.010. PMC9387775.
- **Key finding:** SWATH/DIA-MS proteomes of **949 cancer cell lines across 28 tissue types** (8,498 proteins; ProCan-DepMapSanger). Ovarian composition follows the CCLE/GDSC panel (subtype breakdown not extracted; would require Table S1 — search-limitation). Pan-cancer, not ovarian-specific.
- **Relevance:** Medium — confirms pan-cancer proteomic resources include ovarian lines but are not ovarian-panel resources.
- **Source type:** Resource

> **No *Scientific Data* (or other) DATA DESCRIPTOR dedicated to an ovarian cancer cell-line multi-omic resource was found** across ≥6 varied queries, including preprint-oriented searches (bioRxiv, Research Square, Zenodo, Synapse). Labelled a literature-gap, with a residual search-limitation caveat for very recent/obscure preprints (see Gaps).

---

### Q2 — Are there established, molecularly-characterized ovarian CARCINOSARCOMA / MMMT 2D cell lines?

**Direct answer:** The gap is REAL but PARTIAL, and must not be overstated. Ovarian carcinosarcoma 2D
cell lines DO exist in the literature — but they are old (1990s–2009), one-off single-line
establishment reports, characterized only by dated methods (karyotype, IHC, a few oncogenes; NO
genomic/transcriptomic/proteomic omics), and essentially absent from major public repositories
(ATCC/DSMZ) and from CCLE/DepMap. A dedicated 2025 effort to build gynecologic-carcinosarcoma models
FAILED to establish any ovarian carcinosarcoma cell line (PDX only). So VOA5217/VOA5436 would be
genuinely novel as MODERN, OMICS-CHARACTERIZED, panel-integrated ovarian carcinosarcoma cell lines —
but the Descriptor should NOT claim "the first/only ovarian carcinosarcoma cell lines."

#### OV-MZ-22 — ovarian carcinosarcoma line (in vivo differentiation change)
- **Citation:** Möbus VJ, et al. Characterization of a human carcinosarcoma cell line of the ovary established after in vivo change of histologic differentiation. *Gynecol Oncol.* 2001;83(3):523–532. DOI: 10.1006/gyno.2001.6425. PMID: 11733966. Cellosaurus CVCL_E117.
- **Key finding:** Ovarian carcinosarcoma cell line; documented in vivo shift from papillary cystadenocarcinoma to carcinosarcoma (smooth-muscle differentiation). Characterized by morphology/ultrastructure, karyotype, IHC (intermediate filaments, actins), CA-125/CEA, MDR proteins, p53, topoisomerases, xenograft. No omics.
- **Relevance:** High — a bona fide ovarian carcinosarcoma line; the clearest counterexample to "none exist."
- **Source type:** Primary

#### NEYS — ovarian carcinosarcoma line
- **Citation:** Ide Y, Nakahara T, Nasu M, Tominaga N, Ohyama A, Tachibana T, Yasuda M. Establishment and characterization of the NEYS cell line derived from carcinosarcoma of human ovary with special reference to the susceptibility test of anticancer drugs. *Human Cell.* 2009;22(3):72–80.
- **Key finding:** Line from stage IIIc ovarian carcinosarcoma (56-yr-old); spindle/pleomorphic, no contact inhibition, not tumorigenic in nude mice; drug-susceptibility testing (resistant to CDDP, CPT-11, carboplatin, paclitaxel, docetaxel, 5-FU). No omics.
- **Relevance:** High — second clear ovarian carcinosarcoma line.
- **Source type:** Primary

#### LN1 — heterologous mixed Müllerian tumour of ovary
- **Citation:** Cytogenetic, morphologic and oncogene analysis of a cell line derived from a heterologous mixed mullerian tumor of the ovary. *In Vitro Cell Dev Biol Anim.* 1997. DOI: 10.1007/s11626-997-0001-x. PMID: 9196889.
- **Key finding:** Highly aneuploid heterologous MMMT line (undifferentiated mesodermal components); constitutive c-ras, c-erbB2, p53 mRNA expression. Cytogenetics/morphology/oncogene only.
- **Relevance:** High — third ovarian MMMT line.
- **Source type:** Primary

#### JoN (+ sublines LDF, HDF) — ovarian mixed Müllerian tumour / carcinosarcoma
- **Citation:** Development and characterization of a human cell line from an ovarian mixed mullerian tumor (carcinosarcoma). *In Vitro Cell Dev Biol.* DOI: 10.1007/BF02620867. (Authors/year not retrieved — Springer page paywalled; "BF" DOI prefix indicates a pre-1994 print. Search-limitation.)
- **Key finding:** Ovarian carcinosarcoma line JoN established in culture and nude mice; density-separated sublines LDF and HDF.
- **Relevance:** Medium — a fourth ovarian carcinosarcoma line; citation metadata incomplete.
- **Source type:** Primary

#### Ascites-derived in vitro ovarian carcinosarcoma model
- **Citation:** Presence of both Mesenchymal and Carcinomatous Features in an In-vitro Model of Ovarian Carcinosarcoma Derived from Patients' Ascitic Fluid. *Int J Hematol Oncol Stem Cell Res.* PMC4369227 (~2015).
- **Key finding:** In vitro ovarian carcinosarcoma model from ascites retaining mesenchymal + carcinomatous features.
- **Relevance:** Medium — an additional (non-immortalized-line) ovarian carcinosarcoma culture model.
- **Source type:** Primary

#### BUPH:OVSC — ovarian "sarcomatoid carcinoma" line (nomenclature-adjacent)
- **Citation:** Establishment and characterization of a human ovarian sarcomatoid carcinoma cell line BUPH:OVSC. *Int J Gynecol Cancer* (Elsevier S1048891X24032882, ~2024).
- **Key finding:** Recent ovarian **sarcomatoid carcinoma** cell line. Sarcomatoid carcinoma is related to but distinct from carcinosarcoma — flag nomenclature before treating as a carcinosarcoma line.
- **Relevance:** Medium — closest recent line, but histologic label differs; verify definition against the paper.
- **Source type:** Primary

#### Wong / Llaurado Fernandez et al. 2025 — gynecologic carcinosarcoma models (KEY: cell-line establishment FAILED for ovarian)
- **Citation:** Wong NKY, Llaurado Fernandez M, Kim H, et al. Establishment and characterization of preclinical models of human gynecologic tract carcinosarcomas demonstrates targetable FGFR1 alterations. *Transl Oncol.* 2025;63:102591. DOI: 10.1016/j.tranon.2025.102591. PMC12639385.
- **Key finding:** 6 PDX (5 uterine, **1 ovarian/tubo-ovarian**) + **1 cell line (uterine only)**. For the ovarian/tubo-ovarian case only a PDX was obtained; the authors state concurrent cell-line attempts for all cases were unsuccessful. WES/CNV/SNV, IHC, HRD scoring performed. All harboured TP53 mutations.
- **Relevance:** High — from the SAME BC Cancer group; directly documents that establishing an ovarian carcinosarcoma *2D cell line* remains unmet as of 2025. Strongly supports that VOA5217/VOA5436 fill a modern-model gap.
- **Source type:** Primary

---

### Q3 — Do CCLE / DepMap / ProCan contain ovarian RARE subtypes, or are they dominated by HGSC / unspecified adenocarcinoma?

**Direct answer:** The claim "CCLE/DepMap/ProCan lack rare ovarian subtypes" is TOO STRONG for the rare
EPITHELIAL subtypes. CCLE actually contains modest numbers of clear-cell (~10–12), mucinous (~5–7),
endometrioid (~5–6) and LGSOC (~7–8) lines — and, per Domcke, comparatively FEW good HGSC models. The
gap is well-supported specifically for **carcinosarcoma (absent)** and largely for **SCCOHT (only
COV434, historically misannotated)**. Also, CCLE/DepMap ovarian annotations contain documented errors.

#### Barnes et al. 2021 — transcriptional subtyping of CCLE ovarian lines
- **Citation:** Barnes BM, Nelson L, Tighe A, et al. Distinct transcriptional programs stratify ovarian cancer cell lines into the five major histological subtypes. *Genome Med.* 2021;13(1):140. DOI: 10.1186/s13073-021-00952-5. PMC8408985.
- **Key finding:** Classified **44 CCLE epithelial ovarian lines**: HGSOC 16, CCOC 10, LGSOC 8, MOC 5, ENOC 5. Explicitly flags ENOC as underrepresented (only 5/44). No carcinosarcoma or SCCOHT category.
- **Relevance:** High — quantifies rare-epithelial-subtype representation in CCLE; shows they ARE present in modest numbers.
- **Source type:** Primary

#### McCabe et al. 2023 — suitability of CCLE lines by subtype
- **Citation:** McCabe A, Zaheed O, McDade SS, Dean K. Investigating the suitability of in vitro cell lines as models for the major subtypes of epithelial ovarian cancer. *Front Cell Dev Biol.* 2023;11:1104514. DOI: 10.3389/fcell.2023.1104514. PMC9969113.
- **Key finding:** From 56 EOC lines: HGSOC 16, CCOC 12, MOC 7, ENOC 6, LGSOC 7 (unconfirmed), 14 poorly-classified. "Striking need to develop additional ENOC cell lines"; absence of LGSOC primary-tumour reference data. No carcinosarcoma/SCCOHT.
- **Relevance:** High — corroborates modest rare-epithelial representation and pinpoints ENOC (and LGSOC references) as the weakest epithelial categories.
- **Source type:** Primary

#### Karnezis et al. 2021 — CCLE/DepMap ovarian annotation ERRORS (incl. a CHUM line)
- **Citation:** Karnezis AN, et al. Re-assigning the histologic identities of COV434 and TOV-112D ovarian cancer cell lines. *Gynecol Oncol.* 2021;160(2):568–578. DOI: 10.1016/j.ygyno.2020.12.004. PMID: 33328126.
- **Key finding:** **COV434** (originally "granulosa cell tumor") → **SCCOHT** (SMARCA4 loss). **TOV-112D** (originally grade-3 endometrioid) → **dedifferentiated ovarian carcinoma**. Evidence: pathology review, SMARCA4/SMARCA2 IHC, WES/Sanger, SMARCA4 re-expression, hypercalcemic xenograft.
- **Relevance:** High — TOV-112D is a CHUM/Mes-Masson line in the OvCAN collection; its public annotation is outdated. Demonstrates OvCAN's value in providing corrected, uniformly re-characterized annotations. Also confirms SCCOHT enters CCLE only via a formerly-misannotated line.
- **Source type:** Primary

#### SCCOHT cell-line landscape (only ~3–4 lines worldwide)
- **Citation (BIN-67 / SCCOHT-1 model):** Otte A, et al. A tumor-derived population (SCCOHT-1) as cellular model for a small cell ovarian carcinoma of the hypercalcemic type. *Int J Oncol.* 2012. PMID: 22581215. Plus SCCOHT-CH-1: Establishment and characterization of a novel cell line (SCCOHT-CH-1) and PDX models... *(2023)* PMC10587334.
- **Key finding:** SCCOHT is modelled by ~3 established lines — **BIN-67, SCCOHT-1, COV434** (all SMARCA4-mutant, SMARCA2-null) — plus the newer SCCOHT-CH-1. Of these, essentially only COV434 sits in CCLE/ECACC (and was misannotated as granulosa; see Karnezis).
- **Relevance:** High — SCCOHT is barely represented in the standard pan-cancer resources.
- **Source type:** Primary / Review

#### Pan-cancer resources (CCLE proteomics, ProCan) inherit CCLE composition
- **Citation:** Nusinow 2020 (Cell 180:387–402; DOI 10.1016/j.cell.2019.12.023) and Gonçalves 2022 (Cancer Cell 40:835–849; DOI 10.1016/j.ccell.2022.06.010).
- **Key finding:** Both draw ovarian lines from the CCLE/GDSC panel, so they inherit the same rare-epithelial coverage and the same near-absence of carcinosarcoma/SCCOHT. (Per-subtype ovarian counts not extracted — would require their supplementary tables; search-limitation.)
- **Relevance:** Medium — extends the Q3 conclusion from transcriptomics/genomics to the proteomic resources.
- **Source type:** Resource

---

## Search Log

**Q1 (multi-omic ovarian cell-line resource):**
- Query 1: "multi-omic ovarian cancer cell line panel resource transcriptomics proteomics genomics integrated" — 9 examined, 2 relevant
- Query 2: "ovarian cancer cell lines proteogenomic characterization panel RNA-seq proteomics whole exome" — 8 examined, 5 relevant
- Query 3: "Coscia proteomics 26 ovarian cancer cell lines Nature Communications" — 9 examined, 1 relevant (confirmed Coscia 2016)
- Query 4: "ovarian cancer cell line resource data descriptor Scientific Data multi-omics" — 7 examined, 0 ovarian-cell-line descriptors (negative)
- Query 5: "comprehensively characterized cell line panel highly representative clinical ovarian high-grade serous carcinomas" — 5 examined, 2 relevant (Thu 2017; Papp/Ince-type panel)
- Query 6: "Beaufort Ovarian Cancer Cell Line Panel OCCP 39 cell lines morphological subtypes drug response" — 10 examined, 1 relevant (Beaufort 2014)
- Query 7: "ovarian cancer cell line multi-omic resource biorxiv preprint proteome genome transcriptome 2024 2025" — 10 examined, 0 pan-subtype ovarian cell-line resource (negative; surfaced unrelated preprints)
- Query 8: "Canadian ovarian cancer cell line collection TOV OV CHUM Mes-Masson molecular characterization" — 8 examined, 3 relevant (CHUM line-establishment papers)
- Query 9: "ProCan Goncalves 2022 pan-cancer proteomics cell lines DIA-MS ovarian lines included" — 10 examined, 1 relevant (Gonçalves 2022; ovarian breakdown not in results)
- Query 10: "ovarian cancer cell lines combined proteomic genomic transcriptomic profiling subtypes clear cell endometrioid mucinous" — 8 examined, 2 relevant
- Query 11: "Nusinow 2020 quantitative proteomics CCLE cell lines TMT pan-cancer Cell" — 7 examined, 1 relevant (Nusinow 2020)
- Query 12: "ovarian cancer cell line panel bulk RNA-seq TMT proteomics whole exome resource preprint Research Square Zenodo Synapse 2025 2026" — 8 examined, 1 highly relevant (surfaced Shrestha LGSOC multi-omics)
- Query 13: "multiomics low-grade serous ovarian carcinoma cell lines whole exome RNA-seq proteomics therapeutic vulnerabilities published journal" — 9 examined, 3 relevant (Shrestha 2021 Cancer Res + data repositories)
- Query 14: "Domcke 2013 evaluating cell lines tumour models comparison genomic profiles ovarian Nature Communications" — 10 examined, 2 relevant (Domcke 2013; Ince 2015)
- Query 15: "Ince characterization twenty-five ovarian tumour cell lines phenocopy primary tumours Nature Communications 2015" — 7 examined, 1 relevant (Ince 2015)

**Q2 (ovarian carcinosarcoma/MMMT cell lines):**
- Query 1: "ovarian carcinosarcoma cell line establishment characterization malignant mixed Mullerian tumor" — 8 examined, 4 relevant (JoN, LN1, OV-MZ-22, NEYS surfaced)
- Query 2: "uterine ovarian carcinosarcoma MMMT cell line model molecular characterization novel" — 6 examined, 2 relevant (2025 gyn-CS models; uterine CS TCGA)
- Query 3: "Cellosaurus ovarian carcinosarcoma cell line list" — 8 examined, 2 relevant (OV-MZ-22/CVCL_E117; CH1)
- Query 4: "NEYS cell line carcinosarcoma human ovary anticancer drug susceptibility 2009 established" — 10 examined, 3 relevant (NEYS; ascites model)
- Query 5: "ovarian carcinosarcoma cell line JoN LDF HDF OV-MZ-22 CH1 heterologous mixed mullerian tumor established" — 7 examined, 3 relevant
- Query 6: "novel ovarian carcinosarcoma cell line established 2021 2022 2023 2024 2025 molecular characterization" — 7 examined, 0 new carcinosarcoma lines (negative; only HGSC/SCCOHT/CC lines)
- Query 7: "\"JoN\" OR \"LN1\" ovarian carcinosarcoma cell line In Vitro Cellular Developmental Biology..." — 9 examined, 1 relevant (LN1/PMID 9196889)
- Query 8: "NEYS cell line Human Cell 2009 ovarian carcinosarcoma authors..." — 8 examined, 1 relevant (NEYS author list)

**Q3 (CCLE/DepMap/ProCan rare subtypes):**
- Query 1: "CCLE DepMap ovarian cancer cell lines list histological subtype high-grade serous clear cell annotation" — 9 examined, 3 relevant
- Query 2: "SCCOHT small cell carcinoma ovary hypercalcemic type cell lines BIN67 SCCOHT-1 COV434 SMARCA4" — 8 examined, 4 relevant
- Query 3: "ovarian clear cell carcinoma cell lines mucinous cell lines panel CCLE representation rare subtypes underrepresented" — 7 examined, 2 relevant
- Query 4: "DepMap ovarian cancer cell lines number OncoTree subtype serous carcinosarcoma small cell hypercalcemic" — 8 examined, 1 relevant (SCCOHT count; DepMap per-subtype not in results)
- Query 5: "COV434 reclassified SCCOHT SMARCA4 CCLE DepMap granulosa cell tumor misidentified" — 7 examined, 2 relevant (Karnezis 2021)
- WebFetch (Barnes 2021 PMC8408985) — subtype counts extracted
- WebFetch (McCabe 2023 PMC9969113) — subtype counts extracted

**Cross-cutting WebFetch confirmations:** Möbus/OV-MZ-22 (PMID 11733966), Karnezis 2021 (PMID 33328126), Shrestha 2021 (PMID 33441310), Wong 2025 gyn-CS models (PMC12639385).

---

## Gaps

- **Q1 — literature-gap (confirmed, with caveat):** No pan-subtype, ~40-line ovarian cell-line resource integrating transcriptomics + proteomics + genomics, and no dedicated *Scientific Data*/data-descriptor for an ovarian cell-line multi-omic resource, was found. **Caveat (search-limitation):** WebSearch is US-only and may miss very recent or non-indexed preprints on Research Square / Zenodo / Synapse; the closest matches were surfaced via keyword variation, so residual risk is low but non-zero. **Must-cite prior art:** Shrestha 2021 (all-three-omics but LGSOC-only, n=14) and Coscia 2016 (proteomics-primary, n=26). Frame novelty as breadth/scale/subtype-diversity + WES+RNA+TMT combination, not "first."
- **Q2 — partial literature-gap:** "Modern, molecularly/omics-characterized, publicly-deposited ovarian carcinosarcoma 2D cell lines" = genuine gap (supported by the failed 2025 establishment attempt). "Any ovarian carcinosarcoma cell line ever" = NOT a gap (OV-MZ-22, NEYS, LN1, JoN, ascites model, BUPH:OVSC exist). Descriptor language must reflect this distinction.
- **Q3 — largely NOT a gap for rare epithelial subtypes; IS a gap for carcinosarcoma/SCCOHT:** CCLE contains clear-cell, mucinous, endometrioid and LGSOC lines in modest numbers; carcinosarcoma is absent and SCCOHT enters only via a formerly-misannotated line (COV434).
- **Unresolved citation metadata (search-limitations):** JoN establishment paper authors/year (paywalled Springer, DOI 10.1007/BF02620867); exact Oncotarget volume/pages for Thu 2017; Shrestha 2021 proteomics labelling (label-free vs TMT); per-subtype ovarian counts within ProCan/Nusinow supplementary tables. None are load-bearing for the novelty conclusions.

---

## Breadth Flag

**Recalibrate the rare-subtype-breadth claim.** Do NOT assert broadly that "CCLE/DepMap/ProCan lack
rare ovarian subtypes." The evidence supports a sharper, defensible claim:

1. **Rare epithelial subtypes ARE present in CCLE** (clear cell ~10–12, mucinous ~5–7, endometrioid
   ~5–6, LGSOC ~7–8; Barnes 2021, McCabe 2023). OvCAN's edge here is UNIFORM, matched WES + RNA + TMT
   profiling of these subtypes together with CORRECTED annotations — not their mere inclusion.
2. **Carcinosarcoma is genuinely absent** from CCLE/DepMap/ProCan and lacks any modern
   omics-characterized 2D line (Q2). This is OvCAN's strongest rare-subtype novelty (VOA5217/VOA5436).
3. **SCCOHT is minimally represented** (only COV434 in CCLE, historically misannotated; ~3–4 lines
   exist worldwide). A well-annotated SCCOHT line in a multi-omic panel is a real, if smaller,
   contribution.
4. **Annotation quality is a distinct selling point:** documented CCLE/DepMap errors on ovarian lines,
   including the OvCAN-relevant TOV-112D (endometrioid → dedifferentiated carcinoma) and COV434
   (granulosa → SCCOHT) reassignments (Karnezis 2021).

Bottom line: novelty is defensible, but it is BREADTH-ACROSS-SUBTYPES + SCALE (~40) + the specific
WES/RNA/TMT combination + curated/corrected annotation — anchored most strongly by carcinosarcoma and
SCCOHT — rather than "first multi-omic ovarian cell-line resource." Cite Shrestha 2021, Coscia 2016,
Beaufort 2014, Domcke 2013 and Ince 2015 as the landscape the resource extends.
