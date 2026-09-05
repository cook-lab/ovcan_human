# Refinement 08 — Mucinous Ovarian Cell-Line Authenticity

**Scope:** Targeted refinement for the *Scientific Data* Data Descriptor. Our three mucinous lines are **TOV2414, VOA8762, VOA8771**. Question: are "mucinous ovarian" cell lines genuinely ovarian, or misidentified GI (colorectal/appendiceal/pancreatic) carcinomas — and are any of our three lines implicated?

**Bottom line up front:**
- **Cheasley et al. 2019 does NOT study cell lines.** It is a **primary-tumour** genomics study. Its relevance to our resource is that it *establishes MOC as a bona fide ovarian primary entity* (not extra-ovarian metastasis), and it defines the MOC genomic fingerprint (CDKN2A 76%, KRAS 64%, TP53 64%, ERBB2 26%) that can be used to authenticate a line as MOC. The task brief's framing that it addresses "cell-line authenticity/suitability" is incorrect — flagged below.
- **TOV2414 is well authenticated as genuine ovarian mucinous** (Sauriol et al. 2020; Cellosaurus CVCL_A1SR — no misidentification flag; STR profile on record). The discriminating evidence is **SATB2-negative + focal PAX8-positive** (+MUC5AC/MUC2, WT1-negative). CK20/CDX2 positivity is present but does NOT imply GI origin (intestinal differentiation is intrinsic to ovarian MC).
- **VOA8762 and VOA8771 could not be located** in Cellosaurus or in any accessible primary publication. The "VOA" prefix is the BC Cancer/OvCaRe (Huntsman lab, Vancouver) in-house convention (Anglesio et al. 2013). This is a real gap (see Gaps).

---

## Findings

### 1. Cheasley et al. 2019 — genomic taxonomy of MOC (PRIMARY TUMOURS, not cell lines)
- **Citation:** Cheasley D, Wakefield MJ, Ryland GL, et al. "The molecular origin and taxonomy of mucinous ovarian carcinoma." *Nature Communications* 2019;10(1):3935. DOI: 10.1038/s41467-019-11862-x. PMID: 31477716. PMC6718426.
  - **Title correction:** the brief's candidate title "The molecular *landscape* of mucinous ovarian carcinoma" is wrong; the actual title is "The molecular *origin and taxonomy* of mucinous ovarian carcinoma." DOI is correct.
- **Key finding (genomics, from PMC full text):**
  - **CDKN2A loss or mutation: 76%**; **KRAS mutation: 64%**; **TP53 mutation: 64%**; **ERBB2/HER2 amplification: 26%**; RNF43, BRAF, PIK3CA, ARID1A: ~8–12% each.
  - **9p13 amplicon** (minimal region ~33.785–35.159 Mb; genes linked to chromosome condensation/kinetochores) is a notable copy-number driver; **high copy-number-aberration burden associates with worse prognosis**.
  - MOC is "clearly genetically distinct from high-grade serous ovarian, endometrial, gastric and colorectal tumors, **including mucinous colorectal carcinomas and appendiceal neoplasms**." Distinction from pancreatic ductal adenocarcinoma: MOC has **ERBB2 amplification + RNF43 mutation** (absent in PDAC), whereas PDAC has frequent **SMAD4** loss (absent in MOC).
  - Abstract conclusion (verbatim): "Our data conclusively demonstrate that MOC arise from benign and borderline precursors at the ovary and **are not extra-ovarian metastases.**"
- **Cell lines:** **NONE.** The PMC full text contains no cell-line profiling, authentication, or misidentification analysis. It is exclusively primary/borderline/benign tumour tissue.
- **Relevance:** **High** — supplies the MOC genomic fingerprint used to authenticate a line as MOC, and the field's definitive evidence that MOC is a genuine ovarian primary. **Does not itself speak to cell-line authenticity.**
- **Source type:** Primary.

### 2. Meagher et al. 2025 — cellular origins + how to separate ovarian MOC from GI metastasis
- **Citation:** Meagher NS, Köbel M, Karnezis AN, Talhouk A, Anglesio MS, Berchuck A, Gayther SA, Pharoah PPD, Webb PM, Ramus SJ, Gorringe KL. "Cellular origins of mucinous ovarian carcinoma." *The Journal of Pathology* 2025. DOI: 10.1002/path.6407. PMC11985703. (Title/DOI confirmed; brief's candidate DOI correct.)
- **Key finding:**
  - **Best discriminating IHC markers = CK7 + SATB2 (95% accuracy vs lower-GI metastases):** primary ovarian mucinous = **CK7-positive / SATB2-negative**; colorectal/appendiceal metastasis = **CK7-negative / SATB2-positive**.
  - **Pancreatic** metastases are the hardest to exclude — CK7-positive like ovarian; "CK17 and SMAD4 have shown some promise" but individual-marker accuracy is limited.
  - PAX8 is generally low/absent in MOC (despite GWAS association); ER/PR generally low/absent.
  - **KRAS codon usage differs by site:** MOC G12D (39%) and G12V (33%); colorectal skews G13D; pancreatic G12R; lung G12C. **TP53+CDKN2A co-occurrence 35–54%** in MOC (mirrors PDAC) but MOC **lacks APC** (common colorectal/appendiceal) and **lacks SMAD4** (common PDAC). ERBB2 amplification 4–38% (more in expansile-pattern MOC; rare in PDAC).
  - Cell of origin = multifactorial: teratoma-associated (3–8%, germ-cell, now "LOMN"), Brenner/Walthard-nest-associated, and (most plausibly) mucinous/GI metaplasia within Müllerian inclusion cysts.
  - Clinical discriminators: primary MOC tends unilateral, large (mean ~20 cm); metastases often bilateral, smaller, nodular/surface, signet-ring (Krukenberg) → GI.
- **Cell lines:** none mentioned.
- **Relevance:** **High** — the authoritative, current marker panel for arguing ovarian vs GI origin; directly usable to defend line authenticity in the descriptor. (Note: senior authors Anglesio and Gorringe are the same groups behind the VOA lines and the MOC organoid biobank.)
- **Source type:** Review (expert, primary-literature-synthesizing).

### 3. Sauriol et al. 2020 — primary characterization of OUR line TOV2414
- **Citation:** Sauriol SA, Simeone K, Portelance L, et al. (senior author Mes-Masson AM). "Modeling the Diversity of Epithelial Ovarian Cancer through Ten Novel Well Characterized Cell Lines Covering Multiple Subtypes of the Disease." *Cancers* 2020;12(8):2222. DOI: 10.3390/cancers12082222. PMID: 32784519. PMC7465288.
- **Key finding (TOV2414):**
  - Patient 2414: **63-yo woman, stage IIIC mucinous carcinoma**, sampled 2005; chemoresistant; died 4 months post-diagnosis (shortest OS in the cohort).
  - **Authenticity statement (verbatim):** "Absence of SATB2 expression, focal expression of PAX-8, positive expression of MUC5AC and MUC2, and to a lesser extent, positive expression of CK20, strongly suggest that the tumor from patient 2414 is of the **mucinous subtype of ovarian origin, rather than a metastasis.**"
  - Markers: **TP53 wild-type**; **KRAS c.35G>C p.G12A**; **WT1 absent** (excludes HGSC); PAX8 focal in tissue (lost in WB/culture); CK7/8/18/19 positive; **MUC5AC+/MUC2+**; **CK20+** (lesser); **CDX2+** (noted present in ~38.3% of ovarian MC — i.e., intestinal markers do NOT indicate GI origin here); HER2 not detected.
  - **TOV2414 is the only mucinous line** among the 10 (8 HGSC, 1 clear cell TOV3392D).
  - Phenotype: fastest doubling (1.3 d), highest migration (103 µm/h), most carboplatin-resistant (IC50 11.2 µM); forms SC and IP xenografts.
  - Explicit patient-to-line STR matching is not detailed in the main text, but Cellosaurus lists an STR profile (Finding 4).
- **Relevance:** **High** — this is the citable authenticity dossier for TOV2414; the SATB2−/focal-PAX8+ argument is exactly the Meagher/Köbel framework applied to our line.
- **Source type:** Primary.

### 4. Cellosaurus record — TOV-2414 (CVCL_A1SR)
- **Citation:** Cellosaurus, TOV-2414 (CVCL_A1SR). https://www.cellosaurus.org/CVCL_A1SR (Bairoch A. "The Cellosaurus, a cell-line knowledge resource." *J Biomol Tech* 2018;29:25–38; DOI 10.7171/jbt.18-2902-002.)
- **Key finding:** Category cancer cell line; donor female, 63 y; disease **"Ovarian mucinous adenocarcinoma" (NCIt C5243)**; derived from ovary; **NO "Problematic cell line" / ICLAC "Misidentified" / "Contaminated" flag**; **STR profile present** (Amelogenin X; CSF1PO 10,12; D5S818 9,11; D7S820 9,10; D13S317 12; D16S539 13; D21S11 28,30; TH01 8,9.3; TPOX 8; vWA 14); contact Mes-Masson; reference = Sauriol et al. 2020.
- **Relevance:** **High** — independent registry confirmation that TOV2414 is annotated ovarian mucinous and is NOT on any misidentification register; the on-record STR profile supports authentication.
- **Source type:** Primary (database record).

### 5. Cellosaurus — VOA8762 and VOA8771 NOT FOUND
- **Citation:** Cellosaurus search (https://www.cellosaurus.org/), queries "VOA8762" and "VOA8771", July 2026.
- **Key finding:** **Zero hits for both.** No CVCL accession, no histology annotation, no misidentification flag — because there is no entry at all. The "VOA" numbering convention is the BC Cancer/OvCaRe (Huntsman lab) in-house designation (see Finding 7).
- **Relevance:** **High** — establishes that neither line is registered/authenticated in the standard registry; the descriptor should therefore carry the originating lab's own authentication data.
- **Source type:** Primary (database, negative result).

### 6. Barnes et al. 2021 — transcriptomic reclassification of "mucinous" ovarian lines
- **Citation:** Barnes BM, Nelson L, Tighe A, Burghel GJ, Lin I-H, Desai S, McGrail JC, Morgan RD, Taylor SS. "Distinct transcriptional programs stratify ovarian cancer cell lines into the five major histological subtypes." *Genome Medicine* 2021. DOI: 10.1186/s13073-021-00952-5. PMC8408985.
- **Key finding:** NMF on RNA-seq of 44 CCLE ovarian lines + mutation validation. Of 5 lines annotated MOC, **4 cluster as genuine MOC: MCAS, RMUG-S, COV644, JHOM-2B.** **OV-90** (originally undesignated) clusters with MOC and carries MOC mutations. **JHOM-1 reclassified to (tentative) LGSOC.** MOC lines show KRAS ~60%, TP53 ~60%, ERBB2 amp ~26%, plus BRAF/PIK3CA/ARID1A. **"Most reliable MOC models" = MCAS and OV-90** (both TFF3-mRNA-positive). Does **not** mention TOV2414, VOA8762, or VOA8771.
- **Relevance:** **High** — the template method for authenticating a mucinous line from expression + mutation data; also the current shortlist of trustworthy public MOC lines to benchmark ours against.
- **Source type:** Primary.

### 7. Anglesio et al. 2013 — origin of the "VOA" convention; broad misidentification caution
- **Citation:** Anglesio MS, Wiegand KC, Melnyk N, et al. (senior author Huntsman DG). "Type-Specific Cell Line Models for Type-Specific Ovarian Cancer Research." *PLoS One* 2013;8(9):e72162. DOI: 10.1371/journal.pone.0072162. PMC3762837.
- **Key finding:** BC Cancer/OvCaRe study that uses **in-house "VOA#" identifiers** (confirms VOA = this group's convention). Only one mucinous line detailed (MCAS): validated by STR genotyping, IHC (COSP algorithm), TP53 (homozygous 127-bp deletion) + KRAS G12V, and TFF3-mRNA positivity. Authors warn cell-line records "have recently come into question with a number contaminated and redundant cell lines acknowledged" and found 3 contaminated/mislabeled lines in their own stocks. **VOA8762 and VOA8771 are NOT named in this paper.**
- **Relevance:** **High** for provenance of the VOA naming (points to where VOA8762/VOA8771 likely originate) and for the authentication standard (STR + IHC + mutation + TFF3) expected of BC Cancer lines. **Medium** as direct evidence about our two VOA lines (they are not in it).
- **Source type:** Primary.

### 8. McCabe et al. 2023 — model-suitability audit for MOC
- **Citation:** McCabe A, Zaheed O, McDade SS, Dean K. "Investigating the suitability of in vitro cell lines as models for the major subtypes of epithelial ovarian cancer." *Frontiers in Cell and Developmental Biology* 2023;11:1104514. DOI: 10.3389/fcell.2023.1104514.
- **Key finding:** Ranked lines by expression correlation to primary tumours + mutation oncoprint (27 genes) + NMF. Best MOC-correlated lines: **COV644, MCAS, OVCA420, OVCA429, RMUG-S** (all MOC-annotated lines except JHOM-2B ranked in top-20 correlated to MOC tumours). **GTFR230** (from a stage IC MOC tumour) ranked only 48th and was reclassified to LGSOC. Does not mention TOV2414/VOA8762/VOA8771.
- **Relevance:** **High** — second independent -omic framework for judging MOC-model authenticity; reinforces MCAS/RMUG-S/COV644 as trustworthy benchmarks.
- **Source type:** Primary.

### 9. Korch et al. 2012 — documented misidentification among ovarian lines (STR imperative)
- **Citation:** Korch C, Spillman MA, Jackson TA, Jacobsen BM, Murphy SK, Lessey BA, et al. "DNA profiling analysis of endometrial and ovarian cell lines reveals misidentification, redundancy and contamination." *Gynecologic Oncology* 2012;127(1):241–248. DOI: 10.1016/j.ygyno.2012.06.017. PMID: 22710073.
- **Key finding:** STR + p53 SNP + MSI profiling of 51 ovarian lines: **10 redundant**, and **5 (A2008, OV2008, C13, SK-OV-4, SK-OV-6) are actually cervical cancer cells.** Establishes the scale of ovarian-line misidentification and the requirement to STR-authenticate/re-authenticate every lab's stock. (The famous documented reassignments here are to cervical, not GI; no specific classic *mucinous* line is named as GI-origin in this study.)
- **Relevance:** **Medium** — supports the authentication imperative and the descriptor's QC narrative; not directly about mucinous/GI.
- **Source type:** Primary.

### 10. Craig et al. 2026 (preprint) — MOC organoid biobank; authentication methodology (NOT the VOA source)
- **Citation:** Craig O, Salazar C, Abdirahman S, ... Simpson KJ, McNally OM, Gorringe KL. "Comprehensive drug efficacy data for mucinous ovarian carcinoma using a novel and extensive biobank of patient-derived organoid models." *bioRxiv* 2026.04.06.716848 (posted 2026-04-09). URL: https://www.biorxiv.org/content/10.64898/2026.04.06.716848v1 (full text 403-blocked to WebFetch; details from indexed abstract/search).
- **Key finding:** Largest MOC patient-derived organoid cohort (~10× prior; n=19 long-term lines; 70% success), "highly similar to the tumours of origin for genomic and immunohistochemical markers." **Institution = Peter MacCallum Cancer Centre / University of Melbourne (Gorringe lab)** — i.e., **NOT** BC Cancer/Vancouver, so this is **not** the source of the VOA lines. Illustrates current best practice for authenticating MOC models (match organoid to tumour by genomic + IHC markers).
- **Relevance:** **Medium** — methodological/benchmarking value and confirms MOC-model authentication norms; ruled out as VOA provenance.
- **Source type:** Primary (preprint — not peer-reviewed).

### 11. Ince et al. 2015 — validation-marker precedent (context)
- **Citation:** Ince TA, Sousa AD, Jones MA, Harrell JC, et al. "Characterization of twenty-five ovarian tumour cell lines that phenocopy primary tumours." *Nature Communications* 2015;6:7419. DOI: 10.1038/ncomms8419. PMC4473807.
- **Key finding:** Establishes 25 "OCI"-named ovarian lines; validated ovarian identity by PAX8, ER, WT1, p53, p16, HNF1β, ARID1A IHC (~86% positive for ER/PAX8/CK7) + mutation. Uses OCI nomenclature — **no VOA lines** (rules this out as the VOA source). Does not deeply resolve mucinous vs GI (no CK20/CDX2/SATB2 panel).
- **Relevance:** **Medium** — marker-panel precedent for ovarian identity; a source of alternative authenticated lines.
- **Source type:** Primary.

### 12. ICLAC Register of Misidentified Cell Lines (register context)
- **Citation:** International Cell Line Authentication Committee (ICLAC), Register of Misidentified Cell Lines, Version 14 (released 2026-02-15). https://iclac.org/databases/cross-contaminations/
- **Key finding:** 560 misidentified lines with no authentic stock; 163 contaminants (HeLa most common, 145 entries). I did not enumerate the register line-by-line for a mucinous-ovarian→GI entry (search-limitation); however, Cellosaurus (which mirrors ICLAC flags) shows **no misidentification flag on TOV2414**, and VOA8762/VOA8771 are absent from Cellosaurus entirely.
- **Relevance:** **Medium** — establishes the authoritative misidentification registry; used to confirm TOV2414 is unflagged.
- **Source type:** Primary (database).

---

## Synthesis: how to confirm a mucinous line is genuinely OVARIAN (vs GI) from -omic data

Converging recommendations from Findings 1, 2, 3, 6, 8:

1. **Identity/authentication (do first):** STR profile the line and match to the patient tumour of origin; register in Cellosaurus; cross-check ICLAC. (Korch 2012; Anglesio 2013.) TOV2414 has an on-record STR profile; VOA lines need this documented from the originating lab.
2. **IHC / protein (translatable to bulk expression):** ovarian MOC = **CK7+ / SATB2−**, **focal PAX8+**, **WT1−**, ER/PR low; MUC5AC/MUC2+. Colorectal/appendiceal = **CK7− / SATB2+ / CDX2+**. **Do NOT use CK20 or CDX2 alone** — both are frequently positive in genuine ovarian MC (intestinal differentiation; CDX2+ in ~38% of ovarian MC per Sauriol 2020). (Meagher 2025; Sauriol 2020.)
3. **Mutation profile:** MOC fingerprint = **KRAS (~60–64%; codon G12D/G12V), TP53 (~60–64%), CDKN2A loss (up to 76%), ERBB2 amp (~26%), RNF43**. **Presence of APC → colorectal**; **presence of SMAD4 loss → pancreatic**; both should be **absent** in genuine MOC. (Cheasley 2019; Meagher 2025.)
4. **Expression classification:** correlate the line's transcriptome to MOC primary-tumour references and/or NMF-cluster it against the five histotypes; genuine MOC co-clusters with MCAS/RMUG-S/COV644/OV-90 and expresses **TFF3 mRNA**. (Barnes 2021; McCabe 2023.)
5. **Copy number:** high CNA burden + the **9p13 amplicon** is characteristic of high-grade MOC. (Cheasley 2019.)

**Trustworthy public MOC benchmark lines** (for comparison in the descriptor): **MCAS, RMUG-S, COV644, OV-90** (MCAS and OV-90 most reliable per Barnes 2021; MCAS/RMUG-S/COV644/OVCA420/OVCA429 top-correlated per McCabe 2023). Lines to treat with caution: **JHOM-1** and **GTFR230** (reclassified to LGSOC); **EFO-27** (ambiguous endometrioid/mucinous).

---

## Search Log
- Query 1: "mucinous ovarian cancer cell lines misidentified gastrointestinal origin authentication" — 6 examined, 4 relevant (surfaced Meagher 2025, Korch 2012, Barnes 2021).
- Query 2: "Cheasley 2019 molecular landscape mucinous ovarian carcinoma cell line models" — 9 examined, 3 relevant (corrected title; Sauriol/TOV2414; 3AO/ES2).
- Query 3: "TOV2414 mucinous ovarian carcinoma cell line CHUM Mes-Masson" — 9 examined, 2 relevant (Sauriol 2020 = TOV2414 primary paper).
- Query 4: "VOA8762 VOA8771 ovarian cancer cell line mucinous BC Cancer" — 6 examined, 0 relevant (no hit on either line).
- Query 5: "Anglesio type-specific cell line models ovarian cancer VOA mucinous PLoS One 2013" — 10 examined, 2 relevant (Anglesio 2013 = VOA convention).
- Query 6: "Cheasley mucinous ovarian carcinoma PMC full text cell lines authentic contaminated colorectal" — 8 examined, 1 relevant (review context).
- Query 7: "ICLAC misidentified cell lines register mucinous ovarian gastrointestinal colon" — 10 examined, 2 relevant (ICLAC v14).
- Query 8: "Cellosaurus TOV2414 ovarian mucinous problematic cell line" — 5 examined, 1 relevant (led to CVCL_A1SR).
- Query 9: "\"mucinous\" ovarian cell line actually colorectal MCAS RMUG-S GI origin transcriptomic classification" — 9 examined, 4 relevant (Barnes 2021 MOC cluster; RMUG-S origin).
- Query 10: "\"Modeling the Diversity of Epithelial Ovarian Cancer\" ten novel cell lines TOV2414 mucinous PMC Fleury Mes-Masson" — 6 examined, 2 relevant (Sauriol 2020 PMC7465288).
- Query 11: "\"VOA8762\" OR \"VOA8771\" ovarian cancer cell line characterization mucinous" — 7 examined, 0 relevant.
- Query 12: "VOA8762 VOA8771 ovarian mucinous organoid patient-derived Huntsman OvCaRe 2024 2025" — 9 examined, 1 partly relevant (Craig 2026 preprint, later ruled out).
- Query 13: "Korch 2012 DNA profiling endometrial ovarian cell lines misidentification redundancy contamination Gynecologic Oncology PMID DOI" — 8 examined, 2 relevant (pinned PMID/DOI).
- Query 14: "\"Comprehensive drug efficacy data for mucinous ovarian carcinoma\" patient-derived organoid biobank authors institution" — 9 examined, 1 relevant (PeterMac provenance).
- WebFetch (verifications): Meagher 2025 (PMC11985703 ✓); Cheasley 2019 (PMC6718426 ✓ full text — no cell lines); Sauriol 2020 (PMC7465288 ✓ TOV2414); Barnes 2021 (PMC8408985 ✓); Anglesio 2013 (PMC3762837 ✓); Ince 2015 (PMC4473807 ✓); McCabe 2023 (Frontiers ✓); Cellosaurus CVCL_A1SR ✓; Cellosaurus VOA8762/VOA8771 ✓ (0 hits); PubMed 31477716 ✓. Failed fetches: nature.com (auth redirect), MDPI /htm (403), bioRxiv Craig preprint full text (403), ScienceDirect Korch (403) — all worked around via PMC / PubMed / search.

## Gaps
- **VOA8762 and VOA8771 provenance — PRIMARY GAP (mixed search-limitation + likely literature-gap).** Neither line appears in Cellosaurus (searched, 0 hits) nor in any accessible primary publication. The "VOA" prefix is the BC Cancer/OvCaRe (Huntsman lab) convention (Anglesio 2013), so a Vancouver origin is likely, but no characterization paper was found. They may be recently established and not yet published/registered, or described only in a source not indexed/accessible to web search. **Action for the descriptor: obtain STR profiles, histotype IHC (CK7/SATB2/PAX8/WT1), and mutation data for these two lines directly from the originating lab and report them in-house; do not rely on external authentication that does not exist.**
- **Direct GI-origin evidence for any specific classic "mucinous ovarian" line — literature-gap.** The documented ovarian-line problems I found are (a) misidentification to *cervical* lines (Korch 2012) and (b) *histotype* reclassification within ovarian (JHOM-1, GTFR230 → LGSOC; EFO-27 ambiguous). I did not find a specific named "mucinous ovarian" cell line firmly proven to be colorectal/appendiceal/pancreatic in origin. The GI-metastasis concern is strongest at the *primary-tumour/diagnostic* level (Meagher 2025), and Cheasley 2019 argues MOC as a class is genuinely ovarian.
- **ICLAC register not enumerated line-by-line** (search-limitation); relied on Cellosaurus mirroring of ICLAC flags to confirm TOV2414 is unflagged.
- **Craig 2026 is a non-peer-reviewed preprint** and its full text was inaccessible (403); MOC-organoid details are from the indexed abstract only.

## Breadth Flag
Coverage of the four commissioned questions is **strong for 3 of 4**: Cheasley genomics (Q1, with an important framing correction — it is not a cell-line paper), GI-vs-ovarian discrimination markers (Q2), and how to authenticate from -omics (Q4) are well sourced from primary literature. **Q3 is partially answered:** TOV2414 is resolved (authenticated, unflagged), but **VOA8762/VOA8771 are unresolved** — a genuine provenance gap that should be closed by the originating lab's own data rather than further web search. Recommend the team request BC Cancer/OvCaRe documentation for the two VOA lines.
