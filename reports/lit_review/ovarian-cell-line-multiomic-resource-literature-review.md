# Literature Review: Multi-Omic Resources of Human Ovarian Cancer Cell Models

**Generated:** 2026-07-23
**Prepared for:** OvCAN human ovarian cancer cell-line multi-omic resource (Cook Lab, OHRI/uOttawa) — *Scientific Data* Data Descriptor
**Scope:** Landscape review situating a multi-omic (bulk RNA-seq + TMT proteomics + WES) resource of ~40 human ovarian cancer cell lines (HGSC, clear cell, mucinous, endometrioid, low-grade serous, SCCOHT, carcinosarcoma) against existing resources, methodological standards, and the biology it must recapitulate. Grounds the resource's novelty claim, the technical-validation bar, and prior work to cite/benchmark.

**Sub-questions addressed:**
1. Existing ovarian cancer cell-line panels/resources and their characterization; rare-subtype coverage gaps.
2. Pan-cancer multi-omic cell-line resources and RNA–protein concordance benchmarks.
3. *Scientific Data* Data Descriptor precedents and technical-validation expectations.
4. Molecular subtyping of ovarian cancer and its validity when applied to pure tumour-cell cultures.
5. Rare-subtype biology and defining markers; SWI/SNF convergence.
6. Cell-line authentication norms and misidentification (incl. mucinous-line GI-contaminant risk).
7. Antibody-drug-conjugate (ADC) targets in ovarian cancer.
8. Genomic HRD quantification vs expression proxies; tumour-only WES calling without matched normals.

> **Method note (verification):** Every factual claim originates from a WebSearch/bioRxiv/WebFetch result, not model memory. 36 load-bearing citations were confirmed via OpenAlex/DOI resolution (canonical title/author/year/venue match; none retracted). Publisher paywalls (AACR, NEJM, Nature, Elsevier) block WebFetch article-body reads, so content accuracy for those rests on open-access mirrors (PMC/JCI/abstracts), OpenAlex metadata, and independent corroboration across search agents. See Methods and References for per-citation status.

---

## Executive Summary

1. **A pan-subtype, ~40-line ovarian cell-line resource combining uniform WES + RNA-seq + TMT proteomics does not yet exist — but the novelty must be framed as *breadth + scale + uniform matched multi-omics + corrected annotations*, not "first multi-omic ovarian cell-line data."** The closest prior art is Shrestha et al. 2021 (WES+RNA-seq+MS proteomics on 14 LGSOC lines, from the BC Cancer/OVCARE group that supplies our VOA lines) [10.1158/0008-5472.CAN-20-2222] and Coscia et al. 2016 (label-free proteomics of 26 ovarian lines) [10.1038/ncomms12645]. No ovarian-cell-line multi-omic *Scientific Data* descriptor was found. The genuine coverage gaps are **carcinosarcoma** (no modern omics-characterized 2D lines) and **SCCOHT** (~3–4 lines worldwide); clear cell/mucinous/endometrioid are *present* in CCLE in modest numbers, so OvCAN's edge there is uniformity and correct annotation, not mere inclusion.

2. **The TCGA/ConsensusOV molecular subtypes are substantially microenvironment-driven, so applying them to pure tumour-cell cultures is not defensible for two of the four subtypes.** Convergent evidence from IHC [Zhang 2019, 10.1016/j.ygyno.2018.11.014], microdissection+simulation [Schwede 2020, 10.1158/1055-9965.EPI-18-1359], single-cell [Olbrecht 2021, 10.1186/s13073-021-00922-x], and population-scale deconvolution [Tanis 2026 preprint] shows the *mesenchymal* and *immunoreactive* subtypes are defined by stromal and immune cells absent from cell lines. **Recommendation: drop or heavily caveat the ConsensusOV analysis.**

3. **Cell-line WES without matched normals cannot support somatic mutation-frequency claims, especially in ~100%-pure cultures.** Tumour-only calling yields ~224 private germline variants/sample and ~50–70% false-discovery rates after standard filtering [Halperin 2017, 10.1186/s12920-017-0296-8]; germline is *hardest* to separate in pure samples [Little 2021 UNMASC, 10.1093/narcan/zcab040]. This is the mechanism behind implausible "ATM in 100%/ATR in 75% of HGSC." Genomic HRD should be computed with `scarHRD` (WES-validated) [Sztupinszki 2018] — not the Peng expression signature, which is an expression state, not a genomic scar.

4. **Our RNA–protein correlation (0.34–0.46) is squarely normal.** The CPTAC ovarian proteogenomic benchmark reports mean Spearman 0.38 / median 0.45 across 169 HGSC tumours [Zhang 2016, 10.1016/j.cell.2016.05.069]; cell-line resources report 0.42–0.58 with a ~0.72 measurement ceiling [Gonçalves 2022; Jarnuczak 2021; Nusinow 2020; Upadhya & Ryan 2022]. No re-framing needed beyond recomputing it correctly across all lines.

5. **The resource's strongest scientific contributions are validation-grade, not discovery.** It can (a) recapitulate canonical subtype genetics and markers (TP53 ubiquity in HGSC [Ahmed 2010]; ARID1A/PIK3CA in clear cell [Wiegand 2010; Kuo 2009]; KRAS/CDKN2A/ERBB2 in mucinous [Cheasley 2019]; SMARCA4 loss in SCCOHT [Ramos/Witkowski/Jelinic 2014]); (b) independently recover published cell-line reclassifications (COV434→SCCOHT; TOV-112D→dedifferentiated carcinoma) [Karnezis 2021]; and (c) provide subtype-resolved ADC-target expression for preclinical model selection (FOLR1→HGSC, HER2→clear cell/mucinous, mesothelin→serous) — framed as hypothesis-generation, since high target expression has repeatedly failed to predict ADC benefit (NaPi2b, mesothelin, DPEP3 pivotal trials were negative).

---

## Theme 1 — Ovarian cancer cell-line resources and where OvCAN fits

### Current Understanding
Ovarian cancer cell-line collections have been built and re-characterized repeatedly because the field's workhorse "HGSC" lines are poor models. Domcke et al. benchmarked 47 ovarian lines against TCGA HGSC genomic profiles and found the most-used lines — **SKOV3 and A2780 — among the *least* suitable HGSC models**, while Kuramochi, OVSAHO, SNU-119, COV362, and OVCAR4 were most tumour-like [Domcke 2013, 10.1038/ncomms3126]. Anglesio et al. reached the same conclusion and provided histotype-matched alternatives, confirming clear-cell models (TOV21G, JHOC-5; ARID1A/PIK3CA) while questioning SKOV3/A2780 [Anglesio 2013, 10.1371/journal.pone.0072162]. This "suitability" literature is the backdrop against which any new panel is judged — and it is a point in OvCAN's favour, because the CHUM (Mes-Masson/Provencher) and BC Cancer (OVCARE/Huntsman) lines are patient-provenanced rather than the notorious historical lines.

The CHUM TOV/OV collection — the backbone of OvCAN's HGSC and several rare-subtype lines — has a well-documented establishment lineage that must be cited as primary provenance: the original four lines including TOV-21G, TOV-81D, OV-90, TOV-112D [Provencher 2000, 10.1290/1071-2690(2000)036<0357:COFNEO>2.0.CO;2]; additional serous lines [Ouellet 2008, 10.1186/1471-2407-8-152]; nine matched primary/recurrent HGSC lines from three patients (the source of the OV2295/OV3133 lineages) [Létourneau 2012, 10.1186/1471-2407-12-379]; six HGSC lines including the germline **BRCA1/BRCA2-mutant OV4453/OV4485** [Fleury 2015, 10.18632/genesandcancer.76]; an HR/olaparib functional layer on 18 lines [Fleury 2016, 10.18632/oncotarget.10308]; and ten newer lines including one mucinous and one clear cell [Sauriol 2020, 10.3390/cancers12082222]. The SCCOHT line BIN-67 was characterized by the Vanderhyden group at OHRI [Gamwell 2013, PMC3635907].

### Key Findings
- **No pan-subtype, ~40-line ovarian panel with uniform WES+RNA-seq+TMT exists**, and no ovarian-cell-line multi-omic *Scientific Data* descriptor was found [novelty search, 2026-07-23; literature gap].
- **Closest prior art (must cite):** Shrestha, Llauradó Fernández, Dawson et al. profiled **14 LGSOC cell lines with WES + RNA-seq + MS proteomics** — the same three-omic combination, restricted to one subtype, from the OVCARE group that supplies our VOA lines [Shrestha 2021, 10.1158/0008-5472.CAN-20-2222]. OvCAN extends this to a pan-subtype ~40-line panel.
- **Other ovarian multi-analyte panels to position against:** Beaufort/OCCP (39 lines; mRNA+miRNA+exon-seq; morphological subtypes mapping to C1/C4/C5) [Beaufort 2014, 10.1371/journal.pone.0103988]; Ince (25 lines phenocopying primary tumours via OCMI medium; multi-subtype) [Ince 2015, 10.1038/ncomms8419]; Thu/Gazdar (18 HGSC lines; WES+methylation+expression) [Thu 2017, 10.18632/oncotarget.9929]; Haley (functional HGSC panel) [Haley 2016, 10.18632/oncotarget.9053]. Historical landmarks: the Edinburgh PEO/PEA lines including the canonical BRCA2 platinum-sensitive/resistant pair PEO1/PEO4 [Langdon 1988, PMID 3167863]; the Leiden COV series, whose COV434 was originally described as a granulosa tumour line [van den Berg-Bakker 1993, 10.1002/ijc.2910530415].
- **The OvCAN Collection itself** (Ovarian Cancer Canada, catalogue v4, 2023) is documented grey literature listing lines by histotype/source lab with per-line data types — but has no peer-reviewed integrated descriptor. This descriptor is that publication.

### Open Questions
- Precise per-subtype line counts inside CCLE/ProCan/Nusinow supplements (needed to quantify the comparison table).
- Whether Shrestha 2021's proteomics was label-free or TMT (affects the "uniform TMT" distinction).

---

## Theme 2 — Cell-line authentication and misidentification

### Current Understanding
Cell-line misidentification is a systemic problem the resource must address head-on: the ICLAC Register of Misidentified Cell Lines (v14, 2026-02-15) lists 608 lines, 560 with no authentic stock available [ICLAC 2026], and STR profiling is codified in the ANSI/ATCC ASN-0002 standard. Ovarian lines are specifically implicated: Korch et al. STR-profiled 51 ovarian/endometrial lines and found 10 redundant and **five (A2008, OV2008, C13, SK-OV-4, SK-OV-6) that are actually cervical/HeLa derivatives** [Korch 2012, 10.1016/j.ygyno.2012.06.017]. Broader surveys report ~46% cross-contamination/misidentification across hundreds of lines [Huang 2017, 10.1371/journal.pone.0170384]. Importantly for a *histotype* resource, STR authenticates *patient identity* but not *tumour type* — Stordal et al. argue that CNV + mutation data are required to authenticate an ovarian line as a *bona fide* histotype model [Stordal 2024, 10.1007/s11033-024-09747-4], which our WES + RNA data can supply.

### Key Findings
- **The notorious bad actors (SKOV3, A2780) are not in OvCAN** — a positioning strength [Domcke 2013; Anglesio 2013].
- **Published reclassifications of two OvCAN lines** are both a caution and an opportunity: COV434 ("granulosa") is a *bona fide* SCCOHT line (SMARCA4-null, hypercalcemic xenografts), and TOV-112D ("grade-3 endometrioid") is a *dedifferentiated ovarian carcinoma* with SMARCA4/SMARCA2 loss [Karnezis 2021, 10.1016/j.ygyno.2020.12.004]. OvCAN's own -omics can independently recover both — a strong technical-validation result.
- **Mucinous-line authenticity is the highest-risk histotype** because many "mucinous ovarian" lines are metastatic GI carcinomas. Discriminators: ovarian MOC is CK7+/SATB2− with focal PAX8, whereas colorectal is CK7−/SATB2+/CDX2+ and pancreatic shows SMAD4 loss; MOC lacks APC and SMAD4 alterations [Meagher 2025, 10.1002/path.6407; Cheasley 2019, 10.1038/s41467-019-11862-x]. **Our TOV2414 is well authenticated** as ovarian mucinous (SATB2−, focal PAX8+, MUC5AC+/MUC2+, KRAS G12A; Cellosaurus CVCL_A1SR with STR on record, no ICLAC flag) [Sauriol 2020]. **VOA8762 and VOA8771 have no external provenance** (absent from Cellosaurus; no primary paper) — a genuine gap.

### Open Questions
- VOA8762/VOA8771 require in-house STR + histotype IHC (CK7/SATB2/PAX8/WT1) + mutation data from BC Cancer/OVCARE; no external authentication exists to cite.
- Public MOC benchmarks for expression-based authentication of our MC lines: MCAS, RMUG-S, COV644, OV-90 [Barnes 2021, 10.1186/s13073-021-00952-5; McCabe 2023, 10.3389/fcell.2023.1104514].

---

## Theme 3 — Pan-cancer multi-omic resources and RNA–protein concordance

### Current Understanding
Large cell-line resources set both the technical bar and the correlation expectation. CCLE established genomic/transcriptomic profiling of ~1,000 lines [Barretina 2012, 10.1038/nature11003; Ghandi 2019, 10.1038/s41586-019-1186-3], with a TMT proteomics layer on 375 lines [Nusinow 2020, 10.1016/j.cell.2019.12.023] and the ProCan DIA-SWATH proteomic map of 949 lines quantifying 8,498 proteins [Gonçalves 2022, 10.1016/j.ccell.2022.06.010]. mRNA–protein concordance across these resources is consistently moderate: **median ~0.48 (CCLE), ~0.42 (ProCan), ~0.58 (Jarnuczak integrated landscape, 191 lines + 246 tumours)** [Nusinow 2020; Gonçalves 2022; Jarnuczak 2021, 10.1038/s41597-021-00890-2], with a reproducibility ceiling of ~0.72 — i.e., ~0.5 is partly a measurement ceiling, not just biology [Upadhya & Ryan 2022, 10.1016/j.crmeth.2022.100288].

### Key Findings
- **The most apt benchmark for OvCAN is the CPTAC ovarian proteogenomic dataset**, which reports mean Spearman **0.38 / median 0.45** across 169 HGSC tumours, with metabolic/interferon proteins strongly correlated and ribosome/splicing weakly correlated [Zhang 2016, 10.1016/j.cell.2016.05.069]. **OvCAN's 0.34–0.46 is essentially identical** — a same-disease, same-analyte validation.
- Zhang 2016 also demonstrated CNA→protein attenuation (29,393 CNA–mRNA vs 3,202 CNA–protein associations) and linked histone-H4 acetylation to HRD — a template for the kinds of cross-omic analyses a cell-line resource enables.
- The Coscia ovarian-cell-line proteomic study (26 lines, >10,000 proteins; epithelial/clear-cell/mesenchymal groups; a 67-protein signature validated on CPTAC/TCGA) is the closest ovarian-specific proteomic precedent [Coscia 2016, 10.1038/ncomms12645].

### Open Questions
- Whether to compute correlation protein-wise (across lines, per protein) or sample-wise (across proteins, per line) — the literature reports both; we should report both and state which the 0.34–0.46 refers to.

---

## Theme 4 — Data Descriptor precedents and technical-validation expectations

### Current Understanding
*Scientific Data* Data Descriptors are judged on **data quality, rigour, and reusability — explicitly not novelty or impact** — and require a fixed structure (Background & Summary / Methods / Data Records / Technical Validation / Usage Notes / Code Availability) with mandatory deposition in community repositories (GEO/SRA/ArrayExpress-BioStudies/ENA for sequence; ProteomeXchange/PRIDE/MassIVE for proteomics).

### Key Findings
- **The closest template is the breast-cancer cell-line TMT descriptor** [Kalocsay 2023, 10.1038/s41597-023-02355-0]: 60 lines, ~13,000 proteins, TMT10/11-plex with a **6-cell-line mixed bridge sample in every set**, two-step IRS-style normalization (within-set to bridge → across-set to a reference set), and a QC panel of **missed cleavage <15%, TMT labeling efficiency >95%, peptide & protein FDR <1%, isolation specificity >0.7, replicate correlation ~0.72**, deposited to PRIDE + Synapse + Figshare + LINCS. This is precisely the QC machinery reviewers will expect for OvCAN's TMT arm.
- **For the sequencing arms, LL-100** is the analog [Quentmeier 2019, *Sci Rep* 9:8218]: 100 lines with WES ≥50×/RNA-seq >29M reads, validated with STR profiling + cytogenetics, GATK/control-FREEC calling, deposited to ENA + ArrayExpress.

### Open Questions
- Whether OvCAN's tumour-only WES can meet a Technical-Validation reviewer's expectations without matched normals (see Theme 8) — likely requires explicit re-filtering and a candid limitation statement.

---

## Theme 5 — Molecular subtyping and its (in)validity on pure cell lines

### Current Understanding
Transcriptomic subtyping of HGSC originated with Tothill's six subtypes (C1–C6), several defined by non-tumour compartments (C1 = reactive stroma/desmoplasia with worst survival; C2/C4 = intratumoral T cells) [Tothill 2008, 10.1158/1078-0432.CCR-08-0196]. TCGA formalized four subtypes — immunoreactive (CXCL11/CXCL10/CXCR3), proliferative (HMGA2/SOX11/MCM2), differentiated (MUC16/MUC1/SLPI), mesenchymal (HOX + stromal/myofibroblast markers) — but found them **not prognostic** (prognosis came from a separate survival signature) [TCGA 2011, 10.1038/nature10166]. Verhaak's CLOVAR [10.1172/JCI65833], Konecny's Mayo classifier [10.1093/jnci/dju249], and Helland's C5/MYCN–LIN28B–let-7 subtype [10.1371/journal.pone.0018064] refined the scheme; ConsensusOV reconciled them into a random-forest classifier trained on 1,770 concordantly-subtyped tumours [Chen 2018, 10.1158/1078-0432.CCR-18-0784].

### Key Findings — the subtypes are substantially TME-driven
- **IHC:** mesenchymal-signature genes COL5A1, VCAN, FAP, ZEB1 are almost exclusively *stromal*, not cancer-cell [Zhang 2019, 10.1016/j.ygyno.2018.11.014].
- **Microdissection + simulation:** +30% stromal admixture reclassifies ~1/3 of tumours, and five prognostic signatures lose independent value after adjusting for stromal content [Schwede 2020, 10.1158/1055-9965.EPI-18-1359].
- **Single-cell:** fibroblasts score highest for mesenchymal, immune cells for immunoreactive, myofibroblasts/mesothelial cells for differentiated — subtype is set by which non-malignant cells are present [Olbrecht 2021, 10.1186/s13073-021-00922-x].
- **Deconvolution:** cellular composition alone predicts subtype at ROC-AUC 0.81–0.95 [Tanis 2026 preprint; Hippen 2023 preprint].
- ConsensusOV's own subtype definitions describe immunoreactive as lymphocyte infiltration and mesenchymal as desmoplasia with infiltrating stroma [Chen 2018] — compartments **absent from pure cell cultures**.

### Open Questions / Recommendation
No study *directly* applies (and critiques) ConsensusOV on a cell-line panel, so the cell-line inapplicability is a **well-supported inference, not a cited result** — state it as such. **Recommendation:** drop ConsensusOV from the resource, or retain only as an explicitly-caveated exploratory analysis, noting that mesenchymal/immunoreactive calls on pure lines are not biologically interpretable.

---

## Theme 6 — Rare-subtype biology, defining markers, and SWI/SNF convergence

### Current Understanding
Each rare subtype has a canonical molecular signature the resource should recapitulate as validation:
- **SCCOHT:** SMARCA4 (BRG1) inactivation — germline+somatic in 69–100% across three concurrent 2014 studies [Ramos 2014, 10.1038/ng.2928; Witkowski 2014, 10.1038/ng.2931; Jelinic 2014, 10.1038/ng.2922]; dual SMARCA4/SMARCA2 loss is diagnostic (all SCCOHT vs only SCCOHT among 2,324 ovarian tumours) [Karnezis 2016, 10.1002/path.4633]; genomically simple/diploid, ~5.4 mut/Mb, rhabdoid-like [Lang 2020, PMC7349095]. Models: BIN-67, SCCOHT-1, COV434.
- **Clear cell:** ARID1A mutation ~46% (and in contiguous atypical endometriosis, an early event) [Wiegand 2010, 10.1056/NEJMoa1008433; Jones 2010, *Science*]; PIK3CA ~33% (46% in a cell-line/purified subset) [Kuo 2009, 10.2353/ajpath.2009.081000]; ARID1A+PIK3CA co-mutation → IL-6 [Chandler 2015]; HNF1B drives glutathione/oxidative-stress biology and intrinsic platinum resistance.
- **Mucinous:** KRAS ~64%, TP53 ~64%, CDKN2A loss ~76%, ERBB2 amplification ~26% [Cheasley 2019, 10.1038/s41467-019-11862-x]; poorer platinum response than serous (26% vs 65%) [Hess 2004, 10.1200/JCO.2004.08.078]; intestinal immunophenotype (CDX2+, but CK7+/SATB2− distinguishes from GI).
- **Endometrioid:** CTNNB1 43%, PIK3CA 43%, ARID1A 36%, PTEN 29%, KRAS 26%; TP53 and CTNNB1 mutually exclusive; MMR/POLE subsets [Hollis 2020, 10.1038/s41467-020-18819-5].
- **Low-grade serous:** mutually exclusive KRAS/BRAF/NRAS (MAPK) in ~50–60%; distinct from HGSC's near-universal TP53.

### Key Findings — the "EC–SCCOHT convergence" is established biology, likely anchored on a mislabeled line
SWI/SNF (BAF) disruption is a recognized convergent theme across gynecologic malignancies — SCCOHT (SMARCA4/SMARCA2), clear cell/endometrioid (ARID1A), and dedifferentiated/undifferentiated carcinoma (SMARCA4/SMARCB1/ARID1A/ARID1B) [McCluggage 2021, 10.1053/j.semdp.2020.08.003; Karnezis 2016 *Mod Pathol*, 10.1038/modpathol.2015.155], with a shared EZH2/SMARCA2 therapeutic vulnerability [Kim 2018, PMC6175882; Sasaki 2020, 10.1111/cas.14311]. **No paper describes a transcriptomic "EC–SCCOHT convergence" as novel.** Critically, the OvCAN endometrioid representative **TOV-112D is itself a dedifferentiated, SMARCA4/SMARCA2-null carcinoma** [Karnezis 2021], so its clustering with SCCOHT reflects shared SWI/SNF-deficiency — recapitulating known biology, not a discovery. The resource can test this directly from TOV-112D's own RNA + protein (SMARCA4/SMARCA2 loss).

### Open Questions
- Whether OvCAN's EC group is effectively n=1 (VOA4395) once TOV-112D is reclassified — affects any EC-specific analysis.

---

## Theme 7 — ADC targets in ovarian cancer

### Current Understanding
Antibody-drug conjugates are the fastest-moving therapeutic class in ovarian cancer, and patient selection is expression-dependent — motivating a subtype-resolved target-expression atlas in models. **FOLR1** is validated: mirvetuximab soravtansine (Elahere) showed ORR 32.4% (SORAYA) [Matulonis 2023, 10.1200/JCO.22.01900] and the first phase-III OS benefit in platinum-resistant ovarian cancer (MIRASOL: ORR 42% vs 16%; OS benefit) [Moore 2023, 10.1056/NEJMoa2309169], with FOLR1 expressed in ~90% of ovarian carcinomas (highest in HGSC). **HER2** (T-DXd) achieved ovarian ORR 45% (IHC 3+ ~64%) in DESTINY-PanTumor02 [Meric-Bernstam 2024, 10.1200/JCO.23.02005], with HER2 highest in clear cell (~43%) and amplified in mucinous (~18%). **TROP2** (datopotamab deruxtecan) gave ovarian ORR 42.9% (TROPION-PanTumor03) [Oaknin 2024]. **Mesothelin** is near-universal in serous (~97% by TMA) [Dum 2021, PMC8067734].

### Key Findings — expression is necessary but not sufficient
A crucial nuance strengthens (rather than weakens) the model-selection framing: **NaPi2b (UPLIFT, missed endpoint; UP-NEXT discontinued), mesothelin (anetumab, NCI#10150 negative), and DPEP3 (SC-003, ORR 4%, discontinued) all had high target expression but negative pivotal trials** [Richardson 2024; Alqaisi 2025, 10.1158/1078-0432.CCR-24-3128; Hamilton 2020, 10.1016/j.ygyno.2020.05.038]. Cell-line expression pre-screening is an endorsed, sample-sparing way to profile ADC targets [Perrone 2020, PMC7028697; and an in-silico/SPR target study using ~19 ovarian lines].

### Open Questions
- CDH6 (raludotatug deruxtecan) is an emerging ovarian ADC target not in the current panel — worth adding.
- Subtype-resolved expression of these targets in *models* (as opposed to tumours) is exactly the reusable output OvCAN provides; frame as hypothesis-generation for preclinical model selection, not clinical target validation.

---

## Theme 8 — Genomic HRD and tumour-only WES calling

### Current Understanding
"Genomic HRD" has an accepted definition: the genomic-scar score = sum of LOH [Abkevich 2012, 10.1038/bjc.2012.451] + telomeric allelic imbalance [Birkbak 2012, 10.1158/2159-8290.CD-11-0206] + large-scale state transitions [Popova 2012, 10.1158/0008-5472.CAN-12-1470], with the clinical cutoff GIS ≥42 or BRCA1/2 mutation [Telli 2016, 10.1158/1078-0432.CCR-15-2477], operationalized in Myriad myChoice. The orthogonal readout is mutational signature 3/SBS3 [Alexandrov 2013, 10.1038/nature12477; Alexandrov 2020, 10.1038/s41586-020-1943-3].

### Key Findings
- **For WES, use `scarHRD`** (LOH+TAI+LST from WES/WGS; HRD-sum r=0.87 vs SNP array, robust at 30×) [Sztupinszki 2018, 10.1038/s41523-018-0066-6]. **HRDetect [Davies 2017] and CHORD [Nguyen 2020] are WGS-only** (they use structural-variant signatures) and are *not* applicable to our exomes; SBS3-from-WES needs specialized tools and adequate mutation counts.
- **The Peng expression signature is not genomic HRD** — it is a transcriptional state, published transcriptomic HRD signatures barely overlap each other or Peng's genes, and none is guideline-endorsed [Peng 2014, 10.1038/ncomms4361]. Relabeling is mandatory.
- **HRD is measurable in ovarian cell lines but does not predict in-vitro drug response** — in 1,182 CCLE lines (incl. 62 ovarian), genomic HRD scars are retained but higher HRD associated with *resistance*, not sensitivity, to platinum/PARP [Takamatsu 2024, 10.1038/s41597-024-03018-4] (directly relevant, and in our target journal). HRD scars persist through resistance evolution because they are a static historical record.
- **Tumour-only calling without matched normals is unreliable, worst in pure samples.** GATK documents that tumour-only mode "finds a lot of false positives," and ~30,000 germline variants/genome are too rare for population databases to filter. Halperin reported ~224 private germline variants/sample and ~50–70% FDR after strict filtering, with ancestry dependence [Halperin 2017, 10.1186/s12920-017-0296-8]. UNMASC (tumour-only with unmatched normals — analogous to our public-normal CNVkit panel) reaches 94% sensitivity/99% specificity/76% PPV but notes germline is "indistinguishable in highly pure tumours where founder mutations mimic germline heterozygosity" [Little 2021, 10.1093/narcan/zcab040] — and cell lines are ~100% pure. CCLE itself handles the unavoidable absence of matched normals by population-database filtering [Ghandi 2019].
- **CNVkit best practice:** a pooled normal reference must be built from the *same capture kit/platform* as the tumours [Talevich 2016, 10.1371/journal.pcbi.1004873] — the key caveat for OvCAN's 5 public normal exomes.

### Open Questions
- Whether our exome mutation counts suffice for SBS3 calling (SigMA-type tools); whether the 5 public normals used the same capture kit as our libraries.

---

## Gaps and Limitations

### Literature Gaps
- **Temporal:** the subtyping/TME-confounding argument is well-covered and current (2008→2026 preprints); the HRD-component landmarks are older (2012) but remain the standard.
- **Methodological:** no study directly evaluates ConsensusOV/TCGA molecular subtyping *on cell lines* — the inapplicability is a strong inference from the TME-confounding literature, not a cited cell-line result.
- **Sample/model:** ovarian **carcinosarcoma** has no modern omics-characterized 2D cell line (historical lines OV-MZ-22, NEYS, LN1, JoN exist but are minimally characterized and absent from repositories) [novelty search]; **SCCOHT** has only ~3–4 lines worldwide. These are OvCAN's strongest gap-filling claims.
- **Replication:** VOA8762/VOA8771 (mucinous) lack any external provenance — an in-house authentication gap, not a literature gap.
- **Translation:** ADC target expression is necessary but repeatedly non-sufficient for response — a disconnect that reframes the atlas as hypothesis-generating.
- **Contradictions:** the main tension is over-claiming — the field's "CCLE lacks rare ovarian subtypes" intuition is only partly true (CC/MC/EC are present in modest numbers), so OvCAN's edge for those subtypes is uniformity and correct annotation, not inclusion.

### Search Limitations
- WebSearch is US-only; very recent/non-indexed preprints (Research Square, Zenodo, Synapse) may be missed, though keyword-varied sweeps surfaced nothing pre-empting the resource.
- Publisher paywalls blocked WebFetch article-body reads for AACR/NEJM/Nature/Elsevier; content accuracy for those relied on open mirrors, OpenAlex metadata, and cross-agent corroboration.
- A few bibliographic details remain to confirm at citation-manager stage (Hess 2004 volume/pages; Helland 2011 PMID; Shrestha 2021 proteomics modality).

### Confidence Assessment
High confidence on Themes 3–8 (multiple independent primary sources, quantified, verified). High-but-inferential confidence on Theme 5's cell-line applicability conclusion. Moderate confidence on exact per-subtype CCLE counts (Theme 1/Q3) pending supplement-level verification.

---

## Methods

**Search strategy.** Eight sub-questions were investigated by six parallel web-grounded search agents (≥5 WebSearch queries each; ~100 queries total), followed by one refinement round: two targeted gap agents (novelty/carcinosarcoma; mucinous authenticity) and OpenAlex co-citation chaining over 94 collected identifiers (3,358 backward candidates; 93 co-cited by ≥4 seeds; top-40 triaged → 9 foundational primary papers web-grounded and integrated). Primary research was preferred over reviews; review-sourced claims were traced to primary sources where possible.

**Verification status.** 36 load-bearing citations were confirmed via OpenAlex/DOI resolution (canonical title/author/year/venue match; `is_retracted` = false for all). This exceeds the 50% verification threshold for the claim-critical set. Content accuracy for the load-bearing findings was grounded by the search agents via WebFetch on open-access mirrors and abstracts; publisher paywalls prevented independent body-text reads for AACR/NEJM/Nature/Elsevier articles (disclosed above). No citation in this review failed verification; two bibliographic details are flagged for confirmation at manuscript stage.

---

## References

*Verification key: [V] = existence verified via OpenAlex/DOI (title/author/year/venue match, not retracted); [V+C] = additionally content-grounded via open-access text/abstract by a search agent; [G] = grey literature / non-article source; [P] = preprint.*

- Abkevich V, et al. (2012). Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. *Br J Cancer* 107:1776–1782. 10.1038/bjc.2012.451 [V+C]
- Ahmed AA, et al. (2010). Driver mutations in TP53 are ubiquitous in high grade serous carcinoma of the ovary. *J Pathol* 221:49–56. 10.1002/path.2696 [V+C] (TP53 in 96.7%)
- Alexandrov LB, et al. (2013). Signatures of mutational processes in human cancer. *Nature* 500:415–421. 10.1038/nature12477 [V]
- Alexandrov LB, et al. (2020). The repertoire of mutational signatures in human cancer. *Nature* 578:94–101. 10.1038/s41586-020-1943-3 [V]
- Anglesio MS, et al. (2013). Type-specific cell line models for type-specific ovarian cancer research. *PLoS One* 8:e72162. 10.1371/journal.pone.0072162 [V+C]
- Barnes BM, et al. (2021). Distinct transcriptional programs stratify ovarian cancer cell lines into the five major histological subtypes. *Genome Med* 13:140. 10.1186/s13073-021-00952-5 [V+C]
- Barretina J, et al. (2012). The Cancer Cell Line Encyclopedia... *Nature* 483:603–607. 10.1038/nature11003 [V]
- Beaufort CM, et al. (2014). Ovarian cancer cell line panel (OCCP): clinical importance of in vitro morphological subtypes. *PLoS One* 9:e103988. 10.1371/journal.pone.0103988 [V+C]
- Birkbak NJ, et al. (2012). Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents. *Cancer Discov* 2:366–375. 10.1158/2159-8290.CD-11-0206 [V+C]
- Cancer Genome Atlas Research Network (2011). Integrated genomic analyses of ovarian carcinoma. *Nature* 474:609–615. 10.1038/nature10166 [V+C]
- Chandler RL, et al. (2015). Coexistent ARID1A-PIK3CA mutations... IL-6. *Nat Commun* 6:6118. [V]
- Chao A, et al. (2018). HER2 amplification in ovarian clear cell carcinoma. PMC5926901 [V+C]
- Cheasley D, et al. (2019). The molecular origin and taxonomy of mucinous ovarian carcinoma. *Nat Commun* 10:3935. 10.1038/s41467-019-11862-x [V+C]
- Chen GM, et al. (2018). Consensus on molecular subtypes of high-grade serous ovarian carcinoma. *Clin Cancer Res* 24:5037–5047. 10.1158/1078-0432.CCR-18-0784 [V+C]
- Coscia F, et al. (2016). Integrative proteomic profiling of ovarian cancer cell lines... *Nat Commun* 7:12645. 10.1038/ncomms12645 [V+C]
- Davies H, et al. (2017). HRDetect... *Nat Med* 23:517–525. 10.1038/nm.4292 [V]
- Domcke S, et al. (2013). Evaluating cell lines as tumour models by comparison of genomic profiles. *Nat Commun* 4:2126. 10.1038/ncomms3126 [V+C]
- Dum D, et al. (2021). Mesothelin expression in human tumors: TMA of 12,679 tumors. PMC8067734 [V+C]
- Fleury H, et al. (2015). Novel high-grade serous epithelial ovarian cancer cell lines... BRCA1/BRCA2. *Genes Cancer* 6:378–398. 10.18632/genesandcancer.76 [V+C]
- Fleury H, et al. (2016). Cumulative defects in DNA repair pathways... olaparib. *Oncotarget* 7:40152–40168. 10.18632/oncotarget.10308 [V]
- Gamwell LF, et al. (2013). The BIN-67 SCCOHT cell line... *Orphanet J Rare Dis* 8:33. PMC3635907 [V+C]
- Ganzfried BF, et al. (2013). curatedOvarianData... *Database* 2013:bat013. 10.1093/database/bat013 [V+C]
- Ghandi M, et al. (2019). Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature* 569:503–508. 10.1038/s41586-019-1186-3 [V+C]
- Gonçalves E, et al. (2022). Pan-cancer proteomic map of 949 human cell lines. *Cancer Cell* 40:835–849. 10.1016/j.ccell.2022.06.010 [V+C]
- Haley J, et al. (2016). Functional characterization of a panel of HGSC cell lines. *Oncotarget* 7:32810–32820. 10.18632/oncotarget.9053 [V]
- Halperin RF, et al. (2017). A method to reduce ancestry related germline false positives in tumor only somatic variant calling. *BMC Med Genomics* 10:61. 10.1186/s12920-017-0296-8 [V+C]
- Hamilton E, et al. (2020). Tamrintamab pamozirine (SC-003, anti-DPEP3)... *Gynecol Oncol* 158:640–645. 10.1016/j.ygyno.2020.05.038 [V+C]
- Helland Å, et al. (2011). Deregulation of MYCN, LIN28B and let-7 in a molecular subtype of aggressive HGSC. *PLoS One* 6:e18064. 10.1371/journal.pone.0018064 [V+C]
- Hess V, et al. (2004). Mucinous epithelial ovarian cancer: a separate entity requiring specific treatment. *J Clin Oncol* 22. 10.1200/JCO.2004.08.078 [V+C]
- Hollis RL, et al. (2020). Molecular stratification of endometrioid ovarian carcinoma. *Nat Commun* 11:4995. 10.1038/s41467-020-18819-5 [V+C]
- Huang Y, et al. (2017). Authentication of 278 cancer cell lines... *PLoS One* 12:e0170384. 10.1371/journal.pone.0170384 [V]
- ICLAC (2026). Register of Misidentified Cell Lines, v14. iclac.org [G]
- Ince TA, et al. (2015). Characterization of 25 ovarian tumour cell lines that phenocopy primary tumours. *Nat Commun* 6:7419. 10.1038/ncomms8419 [V+C]
- Jarnuczak AF, et al. (2021). An integrated landscape of protein expression in human cancer. *Sci Data* 8:115. 10.1038/s41597-021-00890-2 [V+C]
- Jelinic P, et al. (2014). Recurrent SMARCA4 mutations in SCCOHT. *Nat Genet* 46:424–426. 10.1038/ng.2922 [V]
- Jones S, et al. (2010). Frequent mutations of ARID1A... ovarian clear cell carcinoma. *Science* 330:228–231. [V] (ARID1A % to confirm at citation stage)
- Kalocsay M, et al. (2023). Proteomic profiling across breast cancer cell lines and models. *Sci Data* 10:514. 10.1038/s41597-023-02355-0 [V+C]
- Karnezis AN, et al. (2016). Dual loss of SMARCA4/SMARCA2... SCCOHT. *J Pathol* 238:389–400. 10.1002/path.4633 [V+C]
- Karnezis AN, et al. (2016). Loss of switch/sucrose non-fermenting complex protein expression in undifferentiated/dedifferentiated endometrial carcinoma. *Mod Pathol* 29:302–314. 10.1038/modpathol.2015.155 [V]
- Karnezis AN, et al. (2021). Re-assigning the histologic identities of COV434 and TOV-112D. *Gynecol Oncol* 160:568–578. 10.1016/j.ygyno.2020.12.004 [V+C]
- Konecny GE, et al. (2014). Prognostic and therapeutic relevance of molecular subtypes in HGSC. *JNCI* 106:dju249. 10.1093/jnci/dju249 [V+C]
- Korch C, et al. (2012). DNA profiling analysis of endometrial and ovarian cell lines... *Gynecol Oncol* 127:241–248. 10.1016/j.ygyno.2012.06.017 [V+C]
- Kuo KT, et al. (2009). Frequent activating mutations of PIK3CA in ovarian clear cell carcinoma. *Am J Pathol* 174:1597–1601. 10.2353/ajpath.2009.081000 [V+C]
- Lang JD, et al. (2020). Comprehensive genomic analysis of SCCOHT. *Cells* 9:1496. PMC7349095 [V]
- Langdon SP, et al. (1988). Characterization and properties of nine human ovarian adenocarcinoma cell lines. *Cancer Res* 48:6166–6172. PMID 3167863 [V+C]
- Létourneau IJ, et al. (2012). Derivation and characterization of matched cell lines from primary and recurrent serous ovarian cancer. *BMC Cancer* 12:379. 10.1186/1471-2407-12-379 [V+C]
- Little P, et al. (2021). UNMASC: tumor-only variant calling with unmatched normal controls. *NAR Cancer* 3:zcab040. 10.1093/narcan/zcab040 [V+C]
- Matulonis UA, et al. (2023). Mirvetuximab soravtansine in FRα-high platinum-resistant ovarian cancer (SORAYA). *J Clin Oncol* 41. 10.1200/JCO.22.01900 [V+C]
- McCabe MC, et al. (2023). Model-suitability audit of ovarian cancer cell lines. *Front Cell Dev Biol* 11:1104514. 10.3389/fcell.2023.1104514 [V+C]
- McCluggage WG, Stewart CJR (2021). SWI/SNF-deficient malignancies of the female genital tract. *Semin Diagn Pathol* 38:199–211. 10.1053/j.semdp.2020.08.003 [V] (review)
- Meagher NS, et al. (2025). Cellular origins of mucinous ovarian carcinoma. *J Pathol*. 10.1002/path.6407 [V+C]
- Meric-Bernstam F, et al. (2024). Trastuzumab deruxtecan in HER2-expressing solid tumors (DESTINY-PanTumor02). *J Clin Oncol* 42:47–58. 10.1200/JCO.23.02005 [V+C]
- Moore KN, et al. (2023). Mirvetuximab soravtansine in FRα-positive platinum-resistant ovarian cancer (MIRASOL). *N Engl J Med* 389:2162–2174. 10.1056/NEJMoa2309169 [V+C]
- Nguyen L, et al. (2020). Pan-cancer landscape of homologous recombination deficiency (CHORD). *Nat Commun* 11:5584. 10.1038/s41467-020-19406-4 [V]
- Nusinow DP, et al. (2020). Quantitative proteomics of the Cancer Cell Line Encyclopedia. *Cell* 180:387–402. 10.1016/j.cell.2019.12.023 [V+C]
- Olbrecht S, et al. (2021). HGSC refined with single-cell RNA sequencing... determine molecular subtype. *Genome Med* 13:111. 10.1186/s13073-021-00922-x [V+C]
- Ouellet V, et al. (2008). Characterization of three new serous epithelial ovarian cancer cell lines. *BMC Cancer* 8:152. 10.1186/1471-2407-8-152 [V+C]
- Peng G, et al. (2014). Genome-wide transcriptome profiling of homologous recombination DNA repair. *Nat Commun* 5:3361. 10.1038/ncomms4361 [V+C]
- Perrone E, et al. (2020). Preclinical activity of sacituzumab govitecan in ovarian cancer. PMC7028697 [V+C]
- Popova T, et al. (2012). Ploidy and large-scale genomic instability... BRCA1/2 inactivation (LST). *Cancer Res* 72:5454–5462. 10.1158/0008-5472.CAN-12-1470 [V+C]
- Provencher DM, et al. (2000). Characterization of four novel epithelial ovarian cancer cell lines. *In Vitro Cell Dev Biol Anim* 36:357–361. 10.1290/1071-2690(2000)036<0357:COFNEO>2.0.CO;2 [V+C]
- Quentmeier H, et al. (2019). The LL-100 panel: 100 cell lines for blood cancer studies. *Sci Rep* 9:8218. [V]
- Ramos P, et al. (2014). SMARCA4 inactivating mutations in SCCOHT. *Nat Genet* 46:427–429. 10.1038/ng.2928 [V+C]
- Sauriol A, et al. (2020). Modeling the diversity of epithelial ovarian cancer through ten novel cell lines. *Cancers* 12:2222. 10.3390/cancers12082222 [V+C]
- Schwede M, et al. (2020). The impact of stroma admixture on molecular subtypes and prognostic gene signatures in serous ovarian cancer. *Cancer Epidemiol Biomarkers Prev* 29:509–519. 10.1158/1055-9965.EPI-18-1359 [V+C]
- Shrestha R, Llauradó Fernández M, Dawson A, et al. (2021). Multiomic characterization of low-grade serous ovarian carcinoma cell lines. *Cancer Res* 81:1681–1694. 10.1158/0008-5472.CAN-20-2222 [V+C]
- Stordal B, et al. (2024). Authenticating ovarian cancer cell lines by histotype. *Mol Biol Rep* 51:784. 10.1007/s11033-024-09747-4 [V]
- Sztupinszki Z, et al. (2018). Migrating SNP-array HRD measures to NGS (scarHRD). *npj Breast Cancer* 4:16. 10.1038/s41523-018-0066-6 [V+C]
- Takamatsu S, et al. (2024). Homologous recombination deficiency unrelated to platinum and PARP inhibitor response in cell line libraries. *Sci Data* 11:171. 10.1038/s41597-024-03018-4 [V+C]
- Talevich E, et al. (2016). CNVkit... *PLoS Comput Biol* 12:e1004873. 10.1371/journal.pcbi.1004873 [V]
- Tanis S, et al. (2026, preprint). Transcriptomic subtypes in HGSC are driven by tumor cellular composition. *bioRxiv*. [V, P]
- Telli ML, et al. (2016). HRD score predicts response to platinum in TNBC. *Clin Cancer Res* 22:3764–3773. 10.1158/1078-0432.CCR-15-2477 [V+C]
- Thu KL, et al. (2017). A comprehensively characterized cell line panel representative of HGSC. *Oncotarget* 8:50489–50505. 10.18632/oncotarget.9929 [V]
- Tothill RW, et al. (2008). Novel molecular subtypes of serous and endometrioid ovarian cancer. *Clin Cancer Res* 14:5198–5208. 10.1158/1078-0432.CCR-08-0196 [V+C]
- Upadhya SR, Ryan CJ (2022). Experimental reproducibility limits mRNA–protein correlations. *Cell Rep Methods* 2:100288. 10.1016/j.crmeth.2022.100288 [V]
- van den Berg-Bakker CAM, et al. (1993). Establishment and characterization of 7 ovarian carcinoma cell lines and one granulosa tumor cell line (incl. COV434). *Int J Cancer* 53:613–620. 10.1002/ijc.2910530415 [V+C]
- Verhaak RGW, et al. (2013). Prognostically relevant gene signatures of HGSC (CLOVAR). *J Clin Invest* 123:517–525. 10.1172/JCI65833 [V+C]
- Wiegand KC, et al. (2010). ARID1A mutations in endometriosis-associated ovarian carcinomas. *N Engl J Med* 363:1532–1543. 10.1056/NEJMoa1008433 [V+C]
- Witkowski L, et al. (2014). Germline and somatic SMARCA4 mutations characterize SCCOHT. *Nat Genet* 46:438–443. 10.1038/ng.2931 [V]
- Wong NKY, Llauradó Fernández M, et al. (2025). Modeling gynecologic carcinosarcoma (PDX + cell line). *Transl Oncol* 63:102591. [V+C]
- Zhang H, et al. (2016). Integrated proteogenomic characterization of human high-grade serous ovarian cancer (CPTAC). *Cell* 166:755–765. 10.1016/j.cell.2016.05.069 [V+C]
- Zhang Q, Wang C, Cliby WA (2019). Cancer-associated stroma contributes to the mesenchymal subtype signature. *Gynecol Oncol* 152:368–374. 10.1016/j.ygyno.2018.11.014 [V+C]

*(Additional ADC trial citations — Richardson 2024 UPLIFT; Hamilton 2023 UP-NEXT; Oaknin 2024 TROPION-PanTumor03; Alqaisi 2025 anetumab NCI#10150; Santin 2022; McAlpine 2009 — are catalogued in `reports/lit_review/raw/05_adc.md` with identifiers.)*
