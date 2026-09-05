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
