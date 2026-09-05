# Synthesis & Recommendations — OvCAN Human Ovarian Cancer Cell-Line Multi-Omic Resource

**Date:** 2026-07-23
**Prepared by:** Cook Lab analyst (Claude)
**Inputs integrated:** 4 assessment reports (`reports/assessment/01–04`), the literature review (`reports/lit_review/…`), and the verified sample sheet (`metadata/samples.csv`).
**Purpose:** Convert the assessment + literature landscape into (a) a recommended reframe, (b) a proposed *Scientific Data* manuscript structure, (c) a phased analysis plan, and (d) the decisions that need your input.

---

## 1. Bottom line

The **data are a genuinely valuable, publishable resource**; the **current write-up is mis-framed** as hypothesis-driven discovery and contains several claims that are technical artifacts or over-reach. The re-scope is mostly **subtractive** (strip discovery), **corrective** (fix the WES over-claim and a few analysis bugs), and **additive** (the Data-Records / QC / authentication machinery a descriptor requires and the manuscript currently lacks).

**Recommended framing:** a lean *Scientific Data* **Data Descriptor** whose spine is *"uniform multi-omic profiling of a curated, correctly-annotated, rare-subtype-enriched ovarian cell-line panel that recapitulates known subtype biology — reuse it to pick and interrogate models."* Discovery claims (novel biomarkers, EC–SCCOHT convergence, complement biology, chemoresistance) move to a **future companion paper** or are demoted to caveated usage examples. This matches what the venue rewards (data quality + reusability, explicitly *not* novelty) and what the data can actually support.

**Novelty, stated defensibly** (the literature forces precision here): breadth across subtypes **+** scale (~40 lines) **+** uniform matched WES/RNA-seq/TMT **+** curated/corrected annotations, anchored most strongly by **carcinosarcoma (VOA5217/5436 — no modern omics-characterized 2D lines exist)** and **SCCOHT (~3–4 lines worldwide)**. *Not* "first multi-omic ovarian cell-line data" (Shrestha 2021 did WES+RNA+MS on 14 LGSOC lines; Coscia 2016 proteomed 26 lines), and *not* "CCLE lacks rare subtypes" (it has modest clear-cell/mucinous/endometrioid numbers) — for those subtypes our edge is uniformity + correct annotation.

---

## 2. What the resource actually is (verified)

| Assay | Coverage (analysis set, post-LGSC-exclusion) | Notes |
|---|---|---|
| Bulk RNA-seq (kallisto 0.46) | **31 cell lines**, 86–93% pseudoalignment | 8 Carey LGSOC excluded (FASTQ-corrupt + PI call) |
| TMT proteomics (5× 11-plex) | **31 lines**, 8,430 proteins | PIS + daisy-chain bridge channels; 2 sites (Mes-Masson/Huntsman); Carey LGS excluded |
| WES — CNVkit CN | **23 lines** | one shared pooled diploid-normal reference (5 public exomes) — sound |
| WES — Mutect2 SNV | **23 lines** | **tumor-only, no matched normal — artifact-prone** |
| **All three omics** | **13 lines** (all Mes-Masson/CHUM: 9 HGS incl. OV90/OV3331, 2 CC, 1 EC[TOV112D], 1 MC) | multi-omic integration realistically rests here |

Unique subtypes in the analysis set: HGS ~22, clear cell 8, mucinous 3, endometrioid 2 (→ likely 1 after TOV112D reclassification), MMMT 2, SCCOHT 2, LGS 1 (TOV81D, WES-only — annotated as the sole LGSC in the collection, per PI). **Report coverage per assay; do not claim "31 × 3 omics."** Source of truth: `metadata/samples.csv`.

---

## 3. Claims to fix, drop, or reframe (assessment × literature)

| Current claim | Verdict | Action |
|---|---|---|
| "RNA-seq + proteomics + WES on 31 models" | Over-claim | Report per-assay N (31/31/23; 13 tri-omic) |
| ATM 100% / ATR 75% / BRCA2 "majority" in HGSC | **Artifact** (tumor-only, no matched normal; Halperin ~50–70% FDR; worst in pure lines) | Re-filter (gnomAD/PoN); restrict to canonical drivers; state limitation |
| KRAS in clear-cell/endometrioid lines | Questionable (same cause) | Re-filter; don't narrate as biology |
| TP53 in 100% of HGSC | **Supported** (Ahmed 2010: 96.7%) | Keep as positive-control validation |
| ConsensusOV / TCGA subtypes on 15 HGSC lines | **Not defensible** (subtypes are TME-driven; mesenchymal/immunoreactive need stroma/immune absent from lines) | Drop from main (already cut from ms); supplement only, heavily caveated |
| "Novel EC–SCCOHT convergence (SWI/SNF)" | Established biology, anchored on a mislabeled line | Reframe: recapitulates known SWI/SNF-deficient group; TOV112D is dedifferentiated (Karnezis 2021) |
| "Concordant RNA-protein signature" (2–12 genes) | Over-sold | Relabel "cross-assay-validated markers"; show the genes honestly |
| Expression-based "HRD" (Peng signature) | **Category error** (expression state ≠ genomic scar) | Relabel; if HRD wanted, compute genomic HRD with scarHRD |
| RNA-protein r = 0.34–0.46 "consistent with literature" | **Supported** (CPTAC ovarian 0.38/0.45; ceiling ~0.72) | Recompute properly; benchmark to Zhang 2016 |
| Subtype separation (PCA/clustering/GO/markers) | **Supported for well-sampled subtypes** | Keep as core validation; annotate n; de-emphasize n=2 silhouettes |
| ADC targets (subtype-associated) | Useful reuse utility | Keep as usage example; frame as hypothesis-generation (expression ≠ response) |
| C4BPB / single-line stories | Unsupported (n=1) | Drop from descriptor |

Analysis bugs to fix in re-implementation (from `03_methods_notebooks.md`): CC/MC signature `Gene` overwrite bug; DESeq2 fed counts (not TPM) — verify; unadjusted p in volcanoes; hardcoded `/Volumes/Eevee/` paths + live biomaRt; proteomics `na.omit` listwise deletion; consensusOV hardcoded remap.

---

## 4. Proposed manuscript structure (*Scientific Data* Data Descriptor)

**Background & Summary** — the OvCAN collection; the model-quality problem (SKOV3/A2780 are poor HGSC models; misidentification is rife); the gap (no pan-subtype uniform multi-omic ovarian cell-line resource; carcinosarcoma/SCCOHT scarcely modeled); what we generated. Cite Shrestha 2021 as the LGSOC precedent we extend.

**Methods** — provenance per line (Provencher 2000, Ouellet 2008, Létourneau 2012, Fleury 2015/16, Sauriol 2020, Gamwell 2013); RNA-seq (kallisto); TMT proteomics (**needs Morin-lab pipeline + QC**); WES (Sarek/Mutect2/CNVkit; the pooled-normal reference; **explicit no-matched-normal limitation + re-filtering**); analysis (subtype-recap, gene-set scoring, authentication).

**Data Records** — accessions: GEO/SRA (RNA-seq, WES), PRIDE/MassIVE via ProteomeXchange (proteomics); file-level descriptions; `metadata/samples.csv` as the master table.

**Technical Validation** (the heart of the paper) — RNA depth/alignment; proteomics identifications/missingness/bridge-channel reproducibility (Kalocsay-2023-style); WES coverage/on-target; **STR authentication**; subtype separation + marker/GO recovery of known biology; RNA-protein concordance vs CPTAC (0.38/0.45); canonical genomics (TP53 ubiquity, HGSC CNV/3q); **independent recovery of the COV434→SCCOHT and TOV112D→dedifferentiated reclassifications** (a standout validation).

**Usage Notes** — model-selection guidance by subtype/genotype/ADC-target expression/HRD; caveats (tumor-only mutations; passage differs across assays; ConsensusOV inapplicable to pure lines; mucinous authentication).

**Code Availability** — the reproducible pipeline (renv, numbered scripts).

### Proposed figures (5 main)
1. **Cohort & design** — provenance/source labs, subtype composition, per-assay coverage grid (from `samples.csv`), workflow schematic.
2. **QC / technical validation** — RNA depth+alignment; proteomics protein/peptide counts, missingness, bridge-channel replicate correlation; WES coverage/on-target.
3. **Recapitulates known biology** — PCA/clustering/silhouette (RNA + protein, n annotated); subtype marker + GO recovery; RNA-protein concordance vs CPTAC benchmark.
4. **Genomics & authentication** — CNV landscape (canonical HGSC gains/losses); canonical drivers (TP53) with tumor-only caveat; STR + histotype authentication incl. recovery of the two reclassifications; genomic HRD (scarHRD, caveated).
5. **Usage example** — subtype-resolved ADC-target expression atlas (RNA+protein) + within-HGSC pathway heterogeneity for model selection.

*(Discovery-flavored panels — ConsensusOV, per-subtype PROGENy for MC/CC, EC–SCCOHT overlap, volcanoes, single-line stories — go to Supplement or the companion paper.)*

---

## 5. Analysis plan (phased re-analysis)

**Phase 0 — Setup** *(done/near-done):* project scaffold, `metadata/samples.csv`, `ANALYSIS_LOG.md`. Add `renv` + pinned annotation/gene-set versions; `00_setup.R`.

**Phase 1 — RNA-seq re-processing:** rebuild tximport→**counts**→DESeq2 (verify counts not TPM); pinned Ensembl annotation (no live biomaRt); recompute PCA/clustering/silhouette (annotate n; model site batch as covariate/diagnostic); recompute RNA-protein correlation properly across all lines; restrict formal DE to HGS (n=15) or present descriptively; fix CC/MC Gene bug.

**Phase 2 — Proteomics re-processing:** obtain/《reconstruct》 + document the Morin-lab TMT pipeline from `peptide_ratio.xlsx`; model plex + site batch; use PIS + bridge channels for normalization/QC; replace `na.omit` with principled missingness handling; recompute PCA/clustering/markers.

**Phase 3 — WES re-analysis:** re-filter Mutect2 tumor-only calls (gnomAD population AF + a panel-of-normals / the pooled normals); restrict oncoplots to well-supported canonical drivers; compute **genomic HRD (scarHRD)**; CNV QC + capture-kit concordance check for the public normals; resolve TOV3121D missing MAF + index/compression issues.

**Phase 4 — Authentication (new, high-value):** compile STR profiles (request where missing — VOA8762/8771, others); histotype-authenticate from CNV + mutation + expression (Stordal approach); **independently test COV434→SCCOHT and TOV112D→dedifferentiated (SMARCA4/SMARCA2 from RNA+protein)**; mucinous CK7/SATB2/PAX8/WT1 + cluster vs public MOC benchmarks.

**Phase 5 — Validation analyses:** subtype-recapitulation (separation + marker/GO recovery); RNA-protein concordance vs CPTAC; ADC-target atlas; within-HGSC heterogeneity as the usage example.

**Phase 6 — Deposition & writing:** deposit to GEO/SRA + PRIDE/MassIVE; assemble Data Records; write the descriptor; final figures to `docs/manuscript/figures/`.

**Environment:** system R for most; `mamba activate scverse` for any Python (scarHRD is R; SBS3 tooling if used). WES re-filtering may need HPC (`shellscripts/`).

---

## 6. Decisions — RESOLVED (PI, 2026-07-23)

**A. Framing — RESOLVED: lean Data Descriptor, with the ADC-target expression atlas retained as a featured *usage/validation* example.** Discovery biology stays out (supplement or a future companion paper); the ADC map is presented as resource-usage validation — a widely-wanted example — not clinical target discovery.

**B. WES mutations — RESOLVED: re-filter + canonical drivers only.** gnomAD + panel-of-normals filtering, canonical-driver reporting, a prominent no-matched-normal limitation statement, and genomic HRD via scarHRD. CNV retained as-is (sound).

**C. Deposition / cohort — RESOLVED (with correction): the LGSC/Carey models are EXTERNAL and removed from the resource entirely.** Their data were pulled from a published Carey/OVCARE source (not generated by this project), so they are dropped from the resource — not merely analysis-excluded. Deposit/document the full *generated* dataset (our CHUM / Huntsman / Vanderhyden lines, including single-omic ones, with QC flags). LGSC is therefore **not** a subtype of this resource; cite the Carey/OVCARE LGSOC work (e.g., Shrestha 2021) rather than including it.

**Secondary (proceeding on defaults; confirm when convenient):**
- **D. OV90 / OV3331:** HGS (default applied; Adenocarcinoma conflict recorded).
- **E. TOV81D:** a CHUM/Mes-Masson LGS line (Provencher 2000) with WES we generated → kept in the *generated* deposit, but LGS is not featured as a subtype (n=1, WES-only).
- **F. TOV112D:** verify the Karnezis-2021 dedifferentiated reassignment from our own SWI/SNF (SMARCA4/SMARCA2) RNA+protein; EC likely n=1 (VOA4395) if confirmed.
- **G. Collaborator asks:** Morin-lab TMT pipeline + per-sample passages + QC (GSC); BC Cancer/OVCARE STR + histotype IHC for VOA8762/VOA8771.
- **Open provenance check:** confirm no other samples are external/pulled-from-publication (PI was unsure if others besides LGSC/Carey were removed); the proteomics LGS channels sit in our TMT layout — confirm whether they were physically run in our TMT or merged from external data.

---

## 7. Risks & open items
- **Proteomics reproducibility** is the biggest descriptor risk: the normalization pipeline and QC are not in the archive and must come from the Morin lab, or the TMT arm can't pass Technical Validation.
- **Tumor-only WES** will draw reviewer scrutiny regardless of re-filtering; the honest limitation statement + canonical-driver focus is essential.
- **Passage differs across assays** for nearly every line — must be stated per line per assay (it's in `samples.csv`).
- **VOA8762/8771 provenance** unresolved pending BC Cancer data.
- **LGSC now absent** — Background must not pitch LGS coverage (a subtype the field, and Shrestha 2021, cares about).
