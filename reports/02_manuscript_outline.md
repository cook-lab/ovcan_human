# Manuscript Outline — OvCAN Multi-Omic Ovarian Cancer Cell-Line Resource

**Target venue:** *Scientific Data* — Data Descriptor
**Prepared:** 2026-07-24
**Source material:** `reports/01_multiomic_characterization_results.md` (results backbone), `PROJECT_SPEC.md`, `output/supplement_per_line.csv`, `reports/assets/f_*.png`
**Status:** **SUPERSEDED — historical.** This outline preceded the first descriptor draft; the live manuscript is `docs/manuscript/v4/OvCAN_data_descriptor_v4.md`. It is retained for the structural mapping in §0, which is still useful. **Do not quote statistics from this file.** Several claims here were corrected by the methodological audit of 2026-07-27 (`reports/06_methodological_validity_and_submission_readiness_review.md`) and are annotated inline as *[CORRECTED 2026-07-27]*; the authoritative versions are in `reports/01_multiomic_characterization_results.md` and the v4 draft.

---

## 0. Framing note — how our analysis maps onto the Data Descriptor skeleton

A Data Descriptor is **not** a results paper. It has a fixed structure and, by editorial policy, **draws no biological conclusions and tests no hypotheses** — it describes a dataset and validates its technical quality so others can reuse it. Our current report is organized as a results narrative (F1–F5); it must be re-poured into the Descriptor skeleton. The key re-mapping:

| Our report content | Data Descriptor home | Why it moves there |
|---|---|---|
| F1 cohort/provenance/patient structure | **Methods** (Cell lines) + **Fig 1** overview | Describes *how the resource was assembled* |
| F2 QC (RNA/proteomics/WES) | **Technical Validation** §1 | QC *is* technical validation |
| F3 "recapitulates subtype biology" | **Technical Validation** §2 | Recovering known biology = evidence the data are trustworthy, **not** a discovery |
| F4 canonical drivers, CNV, TP53 | **Technical Validation** §3 | Recovering canonical genetics validates the WES |
| F4 DepMap concordance, SWI/SNF, STR, mucinous | **Technical Validation** §4 (identity) | Authentication = validation of *what the samples are* |
| F4 TOV21G MSI-high | **Usage Notes** (reuse feature) + brief TV note | A rare-model highlight for reusers; DepMap-corroborated part is validation |
| F5 ADC atlas, within-HGSC strata | **Usage Notes** | Worked reuse examples, not results |
| Decisions / limitations / reproducibility | **Usage Notes** (caveats) + **Code Availability** | Reuse guidance + code statement |

**Consequence for tone:** every "we found X" becomes "the data recover X, consistent with [ref], supporting their reliability" or "reusers can select X." No causal or novel-biology claims. Discovery biology (heterogeneity mechanisms, single-line stories) stays deferred to the companion paper, per `PROJECT_SPEC.md`.

---

## 1. Title (options)

1. *A uniform multi-omic resource of 42 human ovarian cancer cell-line models spanning common and rare histotypes* (recommended — leads with scope + the rare-subtype hook)
2. *Matched RNA-seq, quantitative proteomics, and whole-exome sequencing of the OvCAN ovarian cancer cell-line panel*
3. *Multi-omic characterization of ovarian cancer cell-line models across seven histotypes, including SCCOHT and carcinosarcoma*

Guidance: Descriptor titles describe the **dataset**, not a finding. Lead with "multi-omic … cell-line models" + the differentiator (breadth across rare subtypes). State the number.

---

## 2. Abstract (≤170 words, unreferenced)

One paragraph, ~5 sentences:
1. **What & why** — ovarian cancer is histologically diverse; existing omics cell-line resources (CCLE/DepMap) under-represent rare subtypes. 
2. **What we generated** — uniform bulk RNA-seq (n=31), TMT proteomics (n=31), and WES (n=23) for **42 models from 34 patients** across HGSC, clear cell, mucinous, endometrioid, low-grade serous, SCCOHT, and carcinosarcoma; 13 tri-omic.
3. **Provenance** — CHUM (TOV/OV), BC Cancer (VOA), OHRI (BIN67); data deposited in GEO/SRA, ProteomeXchange, and figshare.
4. **Validation** — data are technically sound, recapitulate canonical subtype biology and genetics (TP53 in 11/11 HGSC patients), and molecularly corroborate line identity against DepMap/Cellosaurus.
5. **Reuse** — a consolidated per-line table plus worked examples (ADC-target atlas; a candidate MSI-high clear-cell model) support model selection.

*(Draft to exact word count at writing time.)*

---

## 3. Background & Summary (~500–800 words; cites Fig 1)

- **Para 1 — clinical/biological motivation.** Epithelial ovarian cancer is not one disease; histotypes differ in genetics, origin, and therapy. Rare subtypes (clear cell, mucinous, LGSC, SCCOHT, carcinosarcoma) are clinically distinct and understudied.
- **Para 2 — the resource gap.** Existing cell-line omics (CCLE/DepMap; CCLE proteomics/Coscia 2016; Shrestha 2021) are HGSC/common-subtype-weighted; the Domcke 2013 "suitability" debate and the frequent misannotation of ovarian lines motivate a curated, uniformly processed, identity-checked panel. Rare-subtype coverage (esp. **SCCOHT** ~3–4 lines worldwide; **carcinosarcoma** with no modern 2D omics) is the strongest differentiator.
- **Para 3 — what this resource is.** The OvCAN / Ovarian Cancer Canada collection: 42 models, 34 patients, three matched layers, one processing pipeline, curated + corrected annotations (incl. two published reclassifications). We generated the *omic data*; the lines are pre-existing/published (cite originators). Calibrated novelty (per `PROJECT_SPEC.md`): breadth + scale + uniform matched WES/RNA/TMT + curation — **not** "first ovarian cell-line omics."
- **Para 4 — brief data-generation + validation overview**, pointing to **Fig 1** (workflow + coverage) and previewing the Technical Validation logic (quality → biology recapitulation → identity).
- **Landscape citations to lock via `literature-review`** (already scoped): CCLE 2012/2019; DepMap; Coscia 2016; Shrestha 2021; Domcke 2013; Beaufort 2014; TCGA-OV 2011; Verhaak/Tothill; ConsensusOV (Chen 2018); Karnezis 2021 (SWI/SNF); Mes-Masson/Provencher TOV/OV derivation series; Gamwell 2013 (BIN67); van den Berg-Bakker 1993 (COV434).

---

## 4. Methods (comprehensive; no length limit; subsectioned)

> In a Descriptor, Methods documents how the **data and the validation analyses** were produced — detailed enough to reproduce. Pull directly from scripts `00–20`.

**4.1 Cell lines, provenance, and annotation**
- Sources and original derivation (Provencher 2000; Ouellet 2008; Létourneau 2012; Fleury 2015/16; Sauriol 2020; Gamwell 2013; van den Berg-Bakker 1993).
- Empirical provenance verification (in-house run `NS.1676.003`; 5-plex TMT) vs excluded external Carey/LGSC data.
- **Patient-family structure** and the collapse rule (`metadata/line_family_map.csv`; script 15): 42 lines → 34 patients; families 1369/2295/3133/3291/3121; representative selection by omics score; variant-level confirmation.
- Subtype annotation and the consensus/curation process (`metadata/samples.csv`); histotype counts.

**4.2 RNA-seq generation and processing** (script 01)
- Library prep / sequencing (run NS.1676.003; read config).
- kallisto 0.46 pseudoalignment; **pinned Ensembl 105** tx2gene (cached; biomaRt lazy) — replaces the archived live-biomaRt call; `ignoreTxVersion`.
- tximport → gene level (39,568 genes); DESeq2 object; expressed-gene filter (≥10 counts in ≥2 lines → 22,544); VST.

**4.3 Proteomics generation and processing** (scripts 05/06; Morin lab GSC)
- TMT 5-plex design; sample-to-plex/channel layout.
- **Normalization = Pooled Internal Standard (PIS) common-reference hub** (channel 1 in every plex); ch10 bridge is a separate technical-replicate QC daisy-chain.
- Peptide/protein identification counts (8,430 proteins / 146,830 peptides); presence filter (≥50% → 7,734); explicit statement of **structural per-plex block-missingness** (`output/prot_block_missingness.csv`).

**4.4 Whole-exome sequencing and processing** (F2/F4; scripts 09–13)
- Sarek pipeline; **tumor-only** Mutect2 + CNVkit; capture kit.
- **Variant filtering (root-cause fix):** `FILTER==PASS` + gnomAD/1000G/ESP population-AF removal (document the 25,914→493 PASS example).
- **Somatic-confidence tiering** (`output/wes_driver_tiers.csv`) replacing the uninformative `germline_like_vaf` flag; rationale for tumor-only limits (LOH → VAF→1).
- CNV: autosome-restricted, per-sample median-centered baseline; **chrX excluded** (pooled-normal sex artifact); FGA definition; unmatched-normal + capture-kit + median-centering caveats.
- Build note: archived MAFs are **GRCh38** despite an erroneous `NCBI_Build=GRCh37` header (to correct on deposition).

**4.5 Computational analysis supporting validation** (scripts 02–08, 14–20)
- PCA / silhouette; **joint variance model + commonality decomposition** (`PC ~ subtype + site`); genome-wide variance decomposition (lme4 REML equivalent — note `variancePartition` unavailable in this R build).
- Marker scoring (top-2-of-6 rule) + effect sizes (Cohen's d, AUC; `output/rna_marker_effectsizes.csv`).
- RNA–protein concordance (per-gene / per-line Spearman; dynamic-range terciles).
- ConsensusOV classification (with the TME-driven caveat).
- Mutational signatures (MutationalPatterns; COSMIC SBS cosine; GRCh38 BSgenome).
- External validation: DepMap Public 24Q4 (Figshare 27993248) RNA/mutation concordance; Cellosaurus STR/records API.
- Software versions + seed (`scripts/00_setup.R`).

---

## 5. Data Records (describe every deposited record + repository/accession)

> Gate for submission: data must be in approved repositories with accessions **before** publication. Draft the table now with placeholders; fill accessions when deposited.

**5.1 Repositories (proposed)**
- **GEO** (SuperSeries): RNA-seq raw (FASTQ→SRA) + processed matrices (TPM, counts, VST); WES processed (VCF/MAF, CNV segments) — or WES under a linked SRA/BioProject.
- **ProteomeXchange / PRIDE (or MassIVE):** raw MS + search results — **requires raw MS from Morin lab (open item)**.
- **figshare / Zenodo:** processed multi-omic matrices, `metadata/samples.csv`, `metadata/line_family_map.csv`, `output/supplement_per_line.csv` (Table S1), derived tables.

**5.2 Data Records table** (main **Table 2**) — columns: *Record · Assay · Repository · Accession · Format · n samples · File(s)*.

**5.3 File organization & formats** — directory tree, file naming keyed to `cell_line`, units/columns for each processed matrix, and the sample sheet as the join key.

**⚠ Deposition-completeness items to resolve (affects what "Data Records" can claim):**
- WES **raw reads/BAMs archived?** Report notes tumor BAMs are *not* archived → confirm FASTQs exist for SRA deposit, else Data Records is processed-only (VCF/MAF/CNV) with a stated limitation.
- Raw MS availability for ProteomeXchange (Morin).
- MAF `NCBI_Build` header correction; strip bundled third-party normal exomes from CNVkit outputs.

---

## 6. Technical Validation (the core section; heavily figure-supported)

**§1 — Sequencing and proteome data quality** → **Fig 2**
- RNA: pseudoalignment median 91.1% (85.8–93.1), ~20,119 genes/line, library depth; modest site effect that does not affect detection.
- Proteomics: bridge technical-replicate reproducibility (Pearson 0.991–0.994; CV ~5.3%); protein/peptide yield; **quantified block-missingness** (81.3% complete; 18.7% plex-conditional).
- WES: filtering restores a credible landscape (waterfall panel); on-target/coverage metrics *if raw data permit* (else stated as not computable from archived outputs).

**§2 — The data recapitulate known subtype biology (⇒ trustworthy)** → **Fig 3**
- Subtype separation in RNA (PC1 20.7% / PC2 10.4%) and protein.
- **Biology, not batch:** joint model — site adds ≤0.2% of PC1 beyond subtype; clean cross-site control within clear cell (4–6% site).
- Marker recovery 16/22 (median AUC 0.69, |d| 0.72); GO/pathway recovery; the six "misses" explained as small-effect/shared-lineage.
- RNA–protein concordance per-gene median Spearman **0.40**, on the CPTAC/CCLE/ProCan benchmark (tercile monotonicity 0.30/0.40/0.53). *[CORRECTED 2026-07-27: the narrower protein spread is **consistent with** TMT ratio compression; these data cannot isolate a compression factor or show that it sets a concordance ceiling, because post-transcriptional buffering, turnover, pooled-standard ratioing and structured missingness contribute as well.]*

**§3 — Genomic fidelity** → **Fig 4** (panels A–B)
- **TP53 in 11/11 HGSC patients** (positive control survives patient collapse).
- Patient-level driver frequencies; pseudoreplication correction (CDK12 35%→18–27%); somatic-confidence tiers; defensible somatic BRCA1/2 = 0 (tumor-only caveat).
- Textbook HGSC CNV landscape robust to collapse (3q gain 82%, 17p loss 82%); autosome FGA ordering with honest n=1 caveats for rare subtypes.

**§4 — Cell-line identity / authentication** → **Fig 4** (panels C–D) + text
- **External:** 5 DepMap lines self-match at rank 1/67 (Spearman 0.74–0.88, reciprocal-best); driver cross-check corroborates (incl. TOV21G hypermutation independently seen in DepMap: 568 vs 7–18 damaging).
- **Reclassifications corroborated:** COV434 & BIN67 → SCCOHT (DepMap SMARCA4-damaging; Cellosaurus misclassification flag; Karnezis 2021); TOV112D → dedifferentiated (SWI/SNF-null, multi-omic).
- **STR/mycoplasma statement:** not done in-house; 30/42 have citable originator STR (Cellosaurus accessions in Table S1); 12 (incl. the identity-doubt VOA lines) lack public STR → in-house STR/IHC requested.
- **Histotype consistency + flagged discordances:** OV90 serous-identity caveat; VOA5436/VOA4841/VOA4395; mucinous VOA8762/VOA8771 read intestinal (Fig 5).

---

## 7. Usage Notes (reuse guidance + worked examples; figure-supported)

- **Model selection by subtype/genotype** — entry via the consolidated **Table S1** (`supplement_per_line.csv`): assay availability, QC, patient-family, tiered drivers, identity calls in one row per line.
- **ADC-target expression atlas** (**Fig 6A**; `output/adc_expression.csv`) — MSLN→HGS, HER2→CC/MC recover; FOLR1 bimodal within HGSC; **lead shortlists with RNA, protein confirmatory** (TMT 3–5× compression); expression is necessary-not-sufficient.
- **Rare-subtype reuse features** — SCCOHT (COV434/BIN67) and carcinosarcoma models; **candidate MSI-high / MMR-deficient clear-cell model TOV21G** (Fig 5) for immunotherapy/MMR biology (pending MMR IHC/MSI-PCR).
- **Within-HGSC heterogeneity** (**Fig 6B**) as a model-selection *example only* (n=15; no survival/discovery claims).
- **Reuse caveats (essential):** tumor-only WES → use tiered drivers, not a burden metric; **genomic HRD not computable** from archived data; TMT ratio compression + per-plex missingness; unmatched normals / capture-kit concordance before quantitative CNV reuse; obtain authenticated stocks (esp. VOA lines). *[CORRECTED 2026-07-27: no in-house STR or mycoplasma testing was performed in this project, so the panel must not be described as authenticated; add histotype–centre non-identifiability and the patient-level denominator rule to this list.]*

---

## 8. Code Availability (mandatory statement)

- Numbered scripts `00–20` (System R 4.5.2), runnable top-to-bottom; seed 1234; pinned Ensembl 105 (cached tx2gene); external data cached (`output/external/`); **`renv` lockfile to be finalized last** (environment still settling this round).
- Public repository (GitHub/Zenodo DOI) — decide host; note the `variancePartition` REML fallback.

---

## 9. Consolidated FIGURE PLAN

**Legend:** ✅ ready (asset exists) · 🔧 generatable now from existing outputs · 🎨 schematic/hand-drawn (BioRender/Illustrator) · 🔒 needs raw data we may not have.

### Main figures (recommend 6; `PROJECT_SPEC` targeted 5 — see note)

| Fig | Working title | Section | Panels (source) | Status |
|---|---|---|---|---|
| **1** | Resource overview & data generation | Background/Methods | (A) workflow schematic tissue→line→3 assays→processing→deposition 🎨; (B) **sample×assay coverage matrix**, 42 lines grouped by subtype with patient-family brackets 🔧; (C) subtype/patient count bars 🔧 | NEW |
| **2** | Data quality across three assays | Tech Val §1 | `f_rna_qc` ✅ · `f_prot_bridge` ✅ · `f_prot_compression` ✅ · **WES filtering waterfall** (25,914→493 PASS) 🔧 | Mostly ready |
| **3** | Data recapitulate subtype biology | Tech Val §2 | `f_rna_pca_subtype` ✅ · `f_rna_pca_site` ✅ · `f_variance_partition` ✅ · `f_rna_markers` ✅ · `f_concordance` ✅ | Ready |
| **4** | Genomic fidelity & line identity | Tech Val §3–4 | `f_wes_oncoplot` ✅ · `f_wes_cnv` ✅ · `f_external_concordance` ✅ · `f_auth_swisnf` ✅ | Ready |
| **5** | Rare-subtype & flagged models | Tech Val §4 / Usage | `f_wes_hypermutation` (TOV21G MSI) ✅ · `f_auth_mucinous` ✅ · **SCCOHT SWI/SNF mini-panel** (could reuse from Fig 4) | Ready |
| **6** | Reuse examples | Usage Notes | `f_adc` ✅ · `f_hgs_het` ✅ | Ready |

**5-vs-6 note:** Fig 1 (overview) is near-mandatory for a Descriptor and currently missing — recommend adding it. If a hard 5-figure cap is imposed, merge Fig 5 into Fig 4 (identity/rare-subtype) and keep Fig 1.

### Supplementary figures

| Fig | Title | Purpose | Status |
|---|---|---|---|
| **S1** | Proteomic PCA by subtype | HGS carries separation; plex/site secondary | ✅ `f_prot_pca` |
| **S2** | Passage sensitivity (within-assay + cross-assay) | Passage collinear w/ site; not an independent driver | ✅ `f_passage_check`, `f_passage_check_crossassay` |
| **S3** | RNA sample–sample correlation heatmap | Identity/no-swap QC; families cluster, subtype blocks | 🔧 (from `rna_vst.rds`) |
| **S4** | Per-plex protein overlap (UpSet) + block-missingness detail | Makes the 18.7% structural missingness explicit | 🔧 (from `prot_block_missingness.csv`) |
| **S5** | Genome-wide per-line CNV heatmap (lines × bins) | Shows actual CN profiles for reuse | 🔧 (from CNVkit segments) |
| **S6** | ConsensusOV calls vs intrinsic HGSC strata | Visualizes the TME-driven-label caveat | 🔧 (from `consensusov_calls.csv`) |
| **S7** | Full marker effect-size panel (all 22; d + AUC) | Backs the 16/22 headline with distributions | 🔧 (from `rna_marker_effectsizes.csv`) |
| **S8** | WES coverage / on-target per line | Standard seq-QC | 🔒 needs raw BAM/metrics — include only if recoverable |

### Alternative / supportive visualizations to consider (author's choice)
- **Fig 3 concordance** — add a Bland-Altman/MA-style RNA-vs-protein panel or the tercile-stratified violin as an alternative to the scatter (communicates the narrower protein spread more directly). *[CORRECTED 2026-07-27: "compression ceiling" retired — see §Technical Validation of the v4 draft.]*
- **Fig 1B** — a UMAP of all RNA lines (subtype + family colored) as an alternative/companion to the coverage matrix for a visual "the panel makes sense" cue.
- **Fig 4 CNV** — a per-arm frequency lollipop (per-patient) alongside the landscape, to foreground the collapse-robustness result.
- **Fig 2** — a single combined "depth/coverage per assay" small-multiple if we want one QC figure instead of assay-specific panels.

---

## 10. TABLE PLAN

| Table | Location | Content | Status |
|---|---|---|---|
| **Table 1** | Main | Cohort inventory (condensed): line, subtype, source, patient/family, assay availability, authentication call | 🔧 condense from `supplement_per_line.csv` |
| **Table 2** | Main | **Data Records** (repository/accession/format/n/files) | draft w/ placeholders |
| **Table S1** | Supp | Full per-line supplement (42 × 41) | ✅ `output/supplement_per_line.csv` |
| **Table S2** | Supp | WES driver somatic-confidence tiers | ✅ `output/wes_driver_tiers.csv` |
| **Table S3** | Supp | Marker effect sizes (d, AUC) | ✅ `output/rna_marker_effectsizes.csv` |
| **Table S4** | Supp | Cellosaurus/STR status per line | ✅ `output/cellosaurus_str_status.csv` |
| **Table S5** | Supp | ConsensusOV calls | ✅ `output/consensusov_calls.csv` |
| **Table S6** | Supp | Per-plex protein block-missingness | ✅ `output/prot_block_missingness.csv` |

---

## 11. Open items that gate submission (for PI)

1. **Deposition** — confirm what raw data can be deposited (WES FASTQs? raw MS from Morin?); mint GEO/PRIDE/figshare accessions; correct MAF build header; strip third-party normals.
2. **HRD decision** — recover WES BAMs → Sequenza → scarHRD, or state HRD as not computable (current default).
3. **Mucinous handling** — flag-and-keep VOA8762/VOA8771 as "ovarian origin unconfirmed," or hold out of ovarian-mucinous claims.
4. **STR / IHC** — request for the 12 lines without public STR (priority: the identity-doubt VOA lines; SMARCA4 IHC for COV434/BIN67).
5. **TOV21G** — MMR IHC / MSI-PCR to upgrade the candidate MSI-high call (or keep "candidate").
6. **Authorship / repository host** — confirm author list (`PROJECT_SPEC`) and code host (GitHub/Zenodo DOI).
7. **`renv` lockfile** — finalize last, then reference in Code Availability.
8. **Landscape citations** — run the scoped `literature-review` to lock Background & Summary references.

---

*Next step on your word: (a) generate the 🔧 figures (Fig 1 coverage matrix, WES waterfall, S3–S7) so the figure set is submission-complete, and/or (b) begin the first prose draft section-by-section from this skeleton.*
