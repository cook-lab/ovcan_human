# Figures & Narrative Reconstruction — OvCAN Human Cell-Line Multi-Omic Resource

**Assessment date:** 2026-07-23
**Scope:** Reconstruct the current figure plan, the intended manuscript story arc, how it
evolved, and the fitness of each figure for a *Scientific Data* **Data Descriptor**.
**Sources read:** all 16 `judy_archive/figs/Figure N` folders (via individual PDFs +
`docs/Figures_aug15.pdf`), both table PDFs, `docs/Figure Captions.docx`,
`docs/Manuscript draft #1.docx` (Dec 4 2025), `judy_archive/docs/Final Thesis Feb 26 2026 JS.docx`
(TOC + figure/table lists), `TAC Meeting 1.pptx`, `TAC #2.pptx`,
`jumbo_data_dump_june2025.pptx`, and all 5 lab-meeting PDFs.

---

## 0. The single most important structural fact: there are TWO deliverables with TWO figure numberings

This project has produced **two parallel outputs that must not be conflated**:

| | **MSc Thesis** (Judy Sobh) | **Resource manuscript** |
|---|---|---|
| Latest file | `Final Thesis Feb 26 2026 JS.docx` | `Manuscript draft #1.docx` (Dec 4 2025) |
| Figures | **16** (1:1 with `judy_archive/figs/`) | **4 main** multi-panel figures |
| Framing | "Analysis of the heterogeneity of ovarian carcinoma" (3 aims) | "A multi-omic **resource** of human ovarian cancer cell models" → *Scientific Data* Data Descriptor |
| Audience | thesis committee | ovarian-cancer research community |

**The `judy_archive/figs/` 16 folders ARE the thesis figure set.** They map 1:1 to the thesis
"List of Figures." `docs/Figures_aug15.pdf` is an intermediate compiled superset of the same
content (with a visible work-in-progress note: *"Still need to fix subtype colours and rna vs
protein axes"*). The **manuscript** consolidates that 16-figure thesis set down to **4 dense
multi-panel figures** (Fig Captions.docx, Oct 16 2025; realized in Manuscript draft #1).

The team-lead brief's figure list uses the **thesis/`figs/` numbering** (Fig12 HGSC PROGENy,
Fig13 ConsensusOV, etc.). The manuscript uses **different numbers** (its Fig 1 = QC+WES; Fig 4 =
HGSC PROGENy). Below I anchor on the thesis/`figs/` numbering and cross-reference manuscript panels.

---

## 1. Figure-by-figure table (thesis / `figs/` numbering)

Verdict legend for a *Scientific Data* Data Descriptor: **ESSENTIAL** (data characterization /
technical validation) · **SUPPORTING** (fine as validation-of-known-biology or usage example) ·
**WEAK/DISCOVERY** (speculative, small-n, or belongs in supplement/another paper).

| # | Files (`figs/Figure N/`) | Analysis | Claim it supports | Manuscript panel | Resource verdict |
|---|---|---|---|---|---|
| **1** | `Subtypes_summary.pdf` | Schematic: H&E histology + clinical/molecular facts for 7 subtypes (HGSC, LGSC, MC, EC, CCC, SCCOHT, MMMT) | Background — what the subtypes are | Intro schematic | **SUPPORTING** — good background/overview graphic; not data |
| **2** | `2a_pca_plot_rna_lgs.pdf`, `2b_tsne_plot_rna_lgs_paper.pdf`, `2c_progeny_heatmap.png` | PCA/tSNE + PROGENy incl. **LGSOC P1–P11** columns | "LGSC has a distinct molecular profile" (thesis: *Initial exploratory analyses of LGSC*) | **DROPPED** | **WEAK/DISCOVERY** — LGSC was *excluded for batch effect*; keep only as a documented-exclusion note. See §5 gap. |
| **3** | `3a_readcount_plot.pdf`, `3a_pseudoalignedvalues_plot.pdf`, `3b_rna_vs_protein_OV3331/TOV21G.pdf`, `3c_cnvkit_heatmap.pdf`, `3d_HGSOC_oncoplot.pdf`, `3e_NONHGSOC_oncoplot.pdf` | RNA-seq QC (reads, %pseudoaligned), RNA-vs-protein scatter, CNVkit CNV heatmap (13 lines), Mutect2/GATK oncoplots (HGSC n=9; non-HGSC n=4) | Data is high quality; recapitulates canonical genetics (TP53 in 100% HGSC; KRAS in non-HGSC; HGSC genomic instability) | **Fig 1 b–g** | **ESSENTIAL** (QC + WES). But oncoplot mutation frequencies need matched-normal caveat; non-HGSC oncoplot is n=1–2/subtype (see §4) |
| **4** | `spearman_cor_heatmap.pdf` | Hierarchical clustering of pairwise Spearman correlation (RNA) | Samples cluster by subtype | **Fig 1h** | **ESSENTIAL** |
| **5** | `5a-c` RNA PCA/tSNE/silhouette, `5d-f` protein PCA/tSNE/silhouette | Dimensionality reduction + mean silhouette width per subtype (PC1–10) | Subtypes separate in RNA & protein; MC/CC most variable | **Fig 2 a–d, f–g** | **ESSENTIAL** |
| **6** | `euclidean_distances.pdf` (+ per-subtype panels in Aug15) | Pairwise Euclidean distance (PC1–10), clustered | Quantifies inter-sample separation | **Fig 2e** | **ESSENTIAL** |
| **7** | `RNAsig_heatmap_top100.pdf` | DESeq2 subtype-vs-Other, top-100 genes/subtype, z-scored | Models express their subtype's signature (with named exceptions) | **Fig 3a** | **ESSENTIAL** — core "recapitulates known biology" evidence |
| **8** | `cc/ec/hgs/mc/mmmt/sccoht_volcano.pdf` | Per-subtype DE volcano (6) | Which genes are DE per subtype | supplement/text | **SUPPORTING** — redundant with Fig 7; move to supplement |
| **9** | `CC/EC/HGS/MC/MMMT/SCCOHT_goterms.pdf` | GO over-representation (topGO elim) per subtype | GO terms match known biology (HGSC: DSB repair/nucleosome; MC: oligosaccharide; CC: glutathione) | text/supplement | **SUPPORTING** — good validation framing |
| **10** | `protein_merged_sig_heatmap.pdf`, `rna_merged_sig_heatmap.pdf` | Intersection of RNA (DESeq2) & protein (t-test) DE at p<0.05,|logFC|>1 | "Concordant multi-omic signature" reproduces known markers (GDA/GGT1 CC; RBP1 HGSC) | **Fig 3b** | **SUPPORTING but OVERSOLD** — intersection = **2–12 genes/subtype** (HGS 3, EC 2, MC 3, CC 5, MMMT 7, SCCOHT 12). Reframe as "cross-assay-validated markers," not a signature (see §4) |
| **11** | `adc_antigen_dotplot_protein.pdf`, `adc_antigen_dotplot_RNA.pdf` | Expression of clinical ADC targets (FOLR1, MSLN, TACSTD2/TROP2, ERBB2/HER2, TF, DPEP3, SLC34A2/NaPi2b) | Resource enables target/model selection (ERBB2 in CC+EC; TACSTD2 in HGSC) | **Fig 3c** | **ESSENTIAL/SUPPORTING** — strong "usage" figure for a resource |
| **12** | `HGSC_progeny_heatmap.pdf` | PROGENy + signaling-metagene + functional + metabolic gene-set activity (singscore), HGSC only | Within-HGSC heterogeneity → model selection | **Fig 4** | **SUPPORTING** — good "usage/heterogeneity" example (n=15 is defensible) |
| **13** | `HGSC_ConsensusOV_marker_heatmap.pdf` (+ **Table 2**) | ConsensusOV → TCGA molecular subtype (Differentiated/Immunoreactive/Mesenchymal/Proliferative) for 15 HGSC | HGSC lines span all 4 TCGA subtypes | **DROPPED** | **WEAK/DISCOVERY** — TCGA subtypes are stroma/immune-driven; applying to pure cell lines is questionable, and the marker heatmap does **not** separate cleanly by assigned label. Already cut from manuscript. Keep out (or supplement w/ heavy caveat) |
| **14** | `MC_progeny_heatmap.pdf` | Gene-set activity, mucinous (n=3) | Within-MC heterogeneity | supplement | **WEAK** — n=3; supplement only |
| **15** | `CCC_progeny_heatmap.pdf` | Gene-set activity, clear cell (n=7) | Within-CC heterogeneity | supplement | **SUPPORTING/WEAK** — n=7 borderline; supplement |
| **16** | `EC_SCCOHT_overlap.pdf` (3-panel: PCA, sig heatmap w/ highlighted boxes, Euclidean) | Signature/expression overlap between EC (n=2) and SCCOHT (n=2) | "Partial molecular convergence between EC and SCCOHT (SWI/SNF)" | 2 sentences in text | **WEAK/DISCOVERY** — n=2 vs n=2; interesting hypothesis, not resource-grade. Supplement or drop |

**Tables:** Table 1 = 31-line cohort summary (subtype, stage, chemo status, RNA/proteomics/WES
availability) — **ESSENTIAL**. Table 2 = ConsensusOV assignments — tied to the WEAK Fig 13; drop with it.

---

## 2. Reconstructed story arc

### Manuscript arc (current target — *Scientific Data* Data Descriptor, 4 figures)
A clean technical-validation narrative in the Data Descriptor mould:

1. **Overview & QC** *(Fig 1)* — "We generated RNA-seq, TMT proteomics, and WES for 31 OvCAN
   gold-standard ovarian-cancer cell models across 6 subtypes; the data are high quality
   (45–97M reads, >85% pseudoalignment, RNA–protein R≈0.34–0.46) and WES recapitulates canonical
   genetics (TP53 in 100% of HGSC, KRAS in non-HGSC, HGSC genomic instability incl. 3q gain)."
2. **Data recapitulates known subtype biology** *(Fig 2)* — Spearman clustering, PCA/tSNE,
   Euclidean distance and silhouette widths show subtypes separate in both RNA and protein; MC and
   clear cell are the most internally variable.
3. **Subtype-specific signatures, known + novel** *(Fig 3)* — DE top-100 signatures are consistent
   within subtype (exceptions: VOA6861, OV90, TOV2414); a concordant RNA-protein set reproduces
   known biomarkers (GGT1/GDA, RBP1); GO terms match known biology; ADC targets flagged for usage.
4. **Within-subtype heterogeneity → model selection** *(Fig 4)* — HGSC gene-set/PROGENy activity
   defines a growth-factor cluster, an inflammatory outlier (OV4453), a low-signaling group, and
   genome-transcriptome concordance (OV2295_R2 PI3K; OV3133-R SMARCA4/WNT).

The throughline: **"this dataset captures known subtype biology (so trust it), reveals some novel
patterns, and lets you pick the right model."** That is exactly the right spine for a Data Descriptor.

### Thesis arc (16 figures, 3 aims)
Same spine but expanded and more discovery-flavored: Ch3 collect/QC (incl. LGSC-exclusion figure) →
Ch4 subtype separation + signatures + GO + concordant signature + ADC → Ch5 within-subtype variation
(HGSC gene sets → **ConsensusOV molecular subtyping** → MC → CC) → Discussion (EC–SCCOHT convergence,
intra-subtype variability, limitations).

---

## 3. Narrative-evolution timeline

| Date | Artifact | Framing / what changed |
|---|---|---|
| **Nov 5 2024** | Lab meeting | Earliest. "Analyzing the heterogeneity of ovarian carcinoma." OvCAN table still carries **carboplatin response, olaparib response, tumours-in-mice (days), TP53/BRCA status, PMIDs** — rich functional metadata later dropped from Table 1 |
| **~May 2024–Jun 2025** | `jumbo_data_dump` (61 slides) | Kitchen-sink exploration: per-subtype GSEA (incl. LGSC), marker validation (AGR2/CDX2 for MC, PAX8/FOXM1 for HGSC), CNV split **by source lab** (Masson/Huntsman/LGSC/SCCOHT), volcano+GO for all subtypes. "LGSC excluded due to batch effect" already present |
| **Mar 4 2025** | Lab meeting | Objective = "evaluate and characterize a collection of 31 human cell models" — **heterogeneity/characterization** framing (not yet "resource") |
| **Jun 9 2025** | Lab meeting | Same heterogeneity framing; analyses maturing |
| **TAC Meeting 1** | Proposal | 3 aims: (1) collect/preprocess, (2) evaluate within-subtype similarity [hypothesis-driven: expects TP53/BRCA, MAPK, ARID1A; validate CA-125/HE4/CEA/CA19-9], (3) distinguish subtypes. LGSC already excluded. Not yet "resource" |
| **Aug 5 2025** | Lab meeting + `Figures_aug15.pdf` | Full compiled figure superset (WIP note on colors/axes). Content complete-ish |
| **TAC Meeting 2** | Update | **PIVOT:** objective becomes **"Develop a multi-omic *resource* for ovarian carcinoma models"**. Aims 2&3 swapped. Adds QC, silhouette, ADC, "Unexpected patterns in WES" slide, EC–SCCOHT overlap, EC DLK1. WES/CNV shown by source lab |
| **Oct 16 2025** | `Figure Captions.docx` | Consolidation to **4 manuscript figures**. (Fig 1 caption still lacks WES — oncoplots not yet "sorted") |
| **Nov 6 2025** | `Copy of Manuscript draft #1` | Earlier manuscript draft |
| **Dec 4 2025** | `Manuscript draft #1` | Full Data Descriptor draft. **WES oncoplots + CNV now folded into Fig 1** (late addition). ConsensusOV, LGSC, per-subtype PROGENy (MC/CC), volcanoes = **dropped** to supplement/thesis |
| **Dec 15 2025** | Lab meeting | Title reverts to **"A multi-omic *analysis* of the heterogeneity"** with the 3-aim structure = **thesis-defense track**. RNA–protein example lines updated (OV3331 R=0.38, TOV21G R=0.41); Spearman heatmap re-themed to colorblind-safe teal/navy |
| **Feb 26 2026** | `Final Thesis` | Comprehensive 16-figure, 2-table version. Re-includes everything the manuscript trimmed (LGSC, ConsensusOV, all PROGENy, volcanoes, EC–SCCOHT) |

**What got ADDED over time:** the word "resource"; QC panel; silhouette quantification; WES
(oncoplots+CNV, late); ADC-target usage figure; EC–SCCOHT convergence; ConsensusOV subtyping.
**What got DROPPED/reframed for the manuscript:** LGSC (batch effect); ConsensusOV/TCGA subtyping
(Fig 13/Table 2); per-subtype PROGENy for MC & CC; standalone volcanoes; and — from the earliest
tables — the **functional/clinical metadata** (drug response, xenograft growth).

---

## 4. Headline findings the figures push — and robustness

1. **"WES recapitulates canonical genetics: TP53 in 100% of HGSC, KRAS across non-HGSC, HGSC genomic
   instability."** — *Medium, with caveats.* TP53-in-all-HGSC and HGSC CNV burden are the expected
   hallmark validations (good). **But** the oncoplot also reports **ATM 100%, ATR 78%, BRCA2 67%**
   in HGSC and **KRAS 100%** in "non-HGSC" — these are suspicious. WES on cell lines **without
   matched normals** inflates germline/artifact calls; ATM/ATR at those rates are not credible as
   somatic drivers. The **non-HGSC oncoplot lumps 4 lines of 4 different subtypes** (2 CC, 1 MC, 1 EC),
   so "KRAS 100%" is n=1–2 per subtype, not a subtype claim. **Needs matched-normal handling / a
   documented germline filter before publishing frequencies.**
2. **"Subtypes separate in RNA and protein."** — *High for the well-sampled subtypes.* Spearman/PCA/
   tSNE/Euclidean/silhouette all agree. Honest reporting that MC and clear cell have silhouette ≈0
   (high intra-subtype variability) is a strength. **Caveat:** MMMT silhouette = 0.8 and SCCOHT = 0.46
   are inflated because **n=2** — silhouette on n=2 is near-meaningless. State n on the plot.
3. **"Models express their diagnosed subtype's signature (with documented exceptions)."** — *High.*
   The top-100 heatmap is convincing and the exceptions (VOA6861, OV90, TOV2414) are handled
   transparently. This is the strongest resource-validation result.
4. **"An integrated RNA-protein signature defines each subtype."** — *Low as stated.* The concordant
   set is **2–12 genes** (HGS = CTHRC1/RBP1/WFDC2; EC = COL3A1/ELOVL4). Real cross-assay-validated
   markers, but not a "signature." **Reframe.**
5. **"Novel patterns: EC–SCCOHT convergence; mucinous complement (C4BPB)."** — *Low / hypothesis.*
   n=2 vs n=2 (EC/SCCOHT) and n=3 (MC). Fine as flagged hypotheses in a thesis; not resource headlines.

---

## 5. Proposed figure triage for the resource paper

**KEEP (essential data characterization / technical validation)**
- Cohort **Table 1** + subtype-overview schematic (thesis Fig 1)
- **QC**: RNA-seq reads/alignment; RNA–protein correlation (thesis Fig 3 b/d/e)
- **WES**: CNV heatmap + HGSC oncoplot (thesis Fig 3 c/d) — *with matched-normal caveat*
- **Spearman clustering** (Fig 4); **PCA/tSNE/silhouette** RNA+protein (Fig 5); **Euclidean** (Fig 6)
- **Top-100 subtype signature heatmap** (Fig 7) — the core validation
- **ADC-target expression** (Fig 11) — the best "usage" figure
- **HGSC gene-set/PROGENy heterogeneity** (Fig 12) as the model-selection example

**REWORK**
- Concordant RNA-protein "signature" (Fig 10) → relabel "cross-assay-validated markers," show the
  handful of genes honestly; consider a genome-wide RNA-vs-protein logFC concordance scatter
  (already exists in Aug15 pp.13–14) as the quantitative version
- Non-HGSC oncoplot (Fig 3e) → either per-subtype (impossible at n=1–2) or fold into one WES panel
  labeled by line, not by "subtype frequency"
- Silhouette/summary plots → annotate n per subtype; de-emphasize n=2 values
- Volcanoes (Fig 8) + GO (Fig 9) → supplement; keep 2–3 GO validations in main text

**DROP from main (supplement or leave to thesis)**
- LGSC exploratory (Fig 2) — excluded for batch effect
- **ConsensusOV / TCGA molecular subtyping** (Fig 13 + Table 2) — already cut; assay-inappropriate for pure cell lines
- Per-subtype PROGENy for **MC (Fig 14)** and **CC (Fig 15)** — small n
- **EC–SCCOHT convergence (Fig 16)** — n=2 vs n=2 discovery claim

**NEEDS-NEW (missing for a credible multi-omic Data Descriptor)**
- **Proteomics QC** panel (protein/peptide counts, PSMs, mass accuracy, missing-value/coverage) — flagged "???" in the outline and "future directions" in TAC2; **currently absent**
- **WES QC** panel (mean target coverage, mapping/on-target rate) — flagged "???"; **absent**
- **Per-assay workflow/schematic** diagram (in the outline, not yet drawn)
- **Cell-line authentication (STR)** — expected in a cell-line resource; not present
- **Source/batch structure** figure — data come from ≥3 labs (CHUM/Mes-Masson, Huntsman/VOA, others); LGSC was dropped *specifically* for batch effect, so reviewers will ask whether source confounds the retained lines
- Consider **re-adding the dropped functional metadata** (carboplatin/olaparib response, xenograft
  growth, source, PMIDs) as a metadata table — high resource value, already partially in hand

---

## 6. Key inconsistencies / risks to flag before re-scope

1. **Two figure numberings** (thesis 16 vs manuscript 4) — keep strictly separate in any re-scope doc.
2. **LGSC is half-in, half-out.** Excluded for batch effect and absent from Table 1 (31 lines, no LGS),
   yet the **protein silhouette plot lists "Low-grade serous,"** and LGSOC P1–P11 appear in the
   PROGENy figure and in `output/geneset_score_matrix_LGS.csv` / `protein_abundance_with_subtypes_LGS.csv`.
   A resource built on the OvCAN "gold-standard collection" that **omits LGSC** is a real gap (LGSC is
   named as an understudied target in the Background). Decide: reintroduce with batch correction, or
   explicitly scope it out.
3. **Pearson vs Spearman mislabel:** the manuscript calls the 0.34–0.46 RNA–protein correlation
   "Pearson" (Results) and "Spearman" (Discussion); plots just say "R=". Pin this down.
4. **RNA–protein example lines drifted** across versions (Aug15 OV2295/OV3331 → Dec TOV21G/OV3331) —
   cosmetic, but reconcile.
5. **WES without matched normals** → inflated mutation frequencies (ATM 100%, ATR 78%); and the
   "non-HGSC KRAS 100%" claim rests on 4 mixed-subtype lines. Highest-priority validation issue.
6. **Small-n subtypes** (EC/SCCOHT/MMMT n=2, MC n=3) drive silhouette, "signature," and "convergence"
   claims — quantify honestly and move discovery claims to supplement/thesis.
7. **Multi-omic completeness varies:** proteomics missing for ~5 lines, WES for ~18/31. The resource
   should state per-line coverage plainly (Table 1 already encodes it — surface it in text).

---

*Prepared for the figures/narrative work-stream. Cross-references: cohort details in `metadata`/Table 1;
methods reconstruction lives with the `methods-notebooks` stream (kallisto/DESeq2/singscore/PROGENy/
CNVkit/Mutect2); manuscript prose with the `thesis-manuscript` stream.*
