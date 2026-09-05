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
