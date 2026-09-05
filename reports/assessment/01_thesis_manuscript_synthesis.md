# Thesis & Manuscript Synthesis — OvCAN Human Ovarian Cancer Cell-Model Resource

**Assessment date:** 2026-07-23
**Author of assessment:** Cook Lab analyst (Claude)
**Documents synthesized:**
- MSc thesis: `judy_archive/docs/Final Thesis Feb 26 2026 JS.docx` (Judy Sobh, ~20k words, full read)
- Latest manuscript: `docs/Manuscript draft #1.docx` (Dec 4 2025, full read)
- Older manuscript: `docs/Copy of Manuscript draft #1.docx` (Nov 6 2025, skimmed for framing drift)
- `docs/Figure Captions.docx` (Oct 16 2025)
- Cross-check: `docs/Multiomic samples summary .xlsx` (sample roster)

**Target venue:** *Scientific Data*, "Data Descriptor" article type.
**Authorship (latest draft):** Judy Sobh (first) · Gian Negri · Euridice Carmona · Barbara C. Vanderhyden · Gregg Morin · David Huntsman · Anne-Marie Mes-Masson · David P. Cook (last/corresponding). Proteomics = Morin group (Michael Smith Genome Sciences Centre, BC Cancer); cell lines = Mes-Masson/Carmona (CRCHUM) + Huntsman (BC Cancer) + OHRI/Vanderhyden.

---

## 1. Executive Summary

This is a multi-omic characterization of **31 human ovarian cancer cell-line models** from the "OvCAN Collection" assembled by Ovarian Cancer Canada, drawn from three Canadian sources (OHRI, BC Cancer/BCCRC, CRCHUM). The panel spans six histological subtypes with **highly unequal representation**: HGSC (n=15), clear cell (n=7), mucinous (n=3), endometrioid (n=2), SCCOHT (n=2), MMMT (n=2). Three assays were generated: **bulk RNA-seq** (kallisto → tximport → TPM), **TMT multiplexed proteomics** (Morin lab, SPS-MS3), and **whole-exome sequencing** (Sarek/CNVkit/Mutect2/VEP). The selling point is breadth across *rare* subtypes (CC, MC, EC, SCCOHT, MMMT) that are underrepresented in resources like CCLE, combined with three molecular layers on matched models.

The work is being repackaged from an MSc thesis ("A multi-omic **analysis of the heterogeneity** of ovarian carcinoma") into a *Scientific Data* descriptor ("A multi-omic **resource** of human ovarian cancer cell models"). **Both documents are framed substantially as hypothesis-driven discovery**, not as a data descriptor: they argue for novel biomarkers, a novel EC–SCCOHT molecular convergence, tumor-cell-intrinsic complement/immunoregulatory mechanisms, chemoresistance biology, and ADC target discovery. The technical-validation scaffolding a descriptor requires (Data Records / accessions, per-assay QC, code/data availability, usage notes) is largely **absent or marked "???" / "To do"** in the current outline. Most biology-level claims rest on **n=2–3 per rare subtype**, **one-vs-all differential expression dominated by HGSC**, and **WES without matched normals**, so they are best treated as suggestive-to-unsupported. The strongest, defensible content is the "does the data recapitulate known biology" validation layer (subtype separation, known biomarkers, known GO terms, canonical HGSC genomics). A re-scope toward that validation framing — plus honest documentation of what data actually exists per assay — is the central editorial task.

**Two factual issues to resolve immediately (both verifiable in code/data):**
1. The topline claim "we systematically generated RNA-seq, TMT proteomics, **and WES** from **31 cell models**" is **not supported by the figures**: WES/CNV appears to cover only **~13 lines** (9 HGSC + 2 CC + 1 EC + 1 MC per Fig 3c caption). WES is far from complete and no SCCOHT/MMMT/most-CC lines are sequenced.
2. **ATM in 100% and ATR in 75% of HGSC lines**, and **BRCA2 in the "majority"** of HGSC, are almost certainly **artifacts of calling variants without a matched normal** (published somatic ATM/ATR rates in HGSC are ~2%). The thesis caveats this for ATM/ATR; the **manuscript states it with no caveat at all**.

---

## 2. What Exactly Is This Resource?

### 2.1 Cohort (as described)
- **31 cell models**, six subtypes: HGSC n=15, CCC n=7, MC n=3, EC n=2, SCCOHT n=2, MMMT n=2. (15+7+3+2+2+2 = 31; internally consistent.)
- Sources: OHRI, BC Cancer Research Institute (BCCRC), CRCHUM (Montreal). Subtype = **clinical diagnosis of the source patient** (not molecularly re-derived).
- HGSC lines are the CHUM/Mes-Masson OV*/TOV* collection; SCCOHT = **BIN67, COV434**; CC/MC/MMMT largely VOA* (BC Cancer) + TOV* lines.

### 2.2 Assays
| Assay | Pipeline (as described) | Apparent coverage |
|---|---|---|
| Bulk RNA-seq | Genome Quebec; NovaSeq S4 2×100 PE; kallisto→GRCh38→tximport→TPM | **31 analyzed** (+8 LGSC generated but excluded → 39 total; see §2.4) |
| TMT proteomics | Morin lab (GSC, BC Cancer); SP3, trypsin, TMT, high-pH fractionation, SPS-MS3, UniProt | N **not stated** in either doc; PCA/silhouette imply most subtypes covered |
| WES | Nimblegen SeqCap EZ Exome 3.0; HiSeq 2500 ~50M reads; **Sarek** (CNVkit, Mutect2, VEP); GSE85671 normals as CNV reference | **~13 lines only** (9 HGSC, 2 CC, 1 EC, 1 MC per Fig 3c); **no matched normals** |

### 2.3 What makes the panel notable (the intended pitch)
- **Rare-subtype breadth**: CC, MC, EC, SCCOHT, MMMT are poorly represented in CCLE/Pan-Cancer; this panel deliberately includes them.
- **Three matched molecular layers** on the same models (RNA + protein + exome), enabling concordance analysis.
- **Curated "gold-standard" provenance** (OvCAN / Ovarian Cancer Canada), positioned against the field problem that models are made piecemeal by many labs with uneven validation.

### 2.4 Reconciliation with the raw sample roster (`Multiomic samples summary .xlsx`) — IMPORTANT
The sample sheet contains **many more entries than 31**:
- The **core 31** match the manuscript roster.
- **LGSC / LGSOC_P1–P11** (patient-derived; 8 visible: P1,P2,P3,P4,P5,P7,P10,P11) are present but **excluded** (see §4, batch effect). **31 + 8 LGSC = 39**, which reconciles exactly with the "39 RNA-seq samples" noted in project recon.
- A **large additional block** of lines/replicates (e.g., VOA10841, VOA14202, VOA10816.1, VOA3993/.1, VOA14993, VOA7681, TOV1369.1, TOV3133G.1, OV2978, OV3133R2, OV3291, TOV81D, TOV1946, TOV2223G, TOV2295R, TOV2835EP, TOV2881EP, TOV2929D, TOV2978G, TOV3121D/EP, and `CL-1_1`…`CL-12_1`, `14CL`, `15CL`) appears in the sheet. These likely reflect **replicates, additional/failed samples, or assay-specific plex labels** — i.e., the *generated* data are broader than the *analyzed* 31.
- **Implication for a data descriptor:** *Scientific Data* expects the descriptor to document **all deposited data**, including excluded/replicate/failed samples with rationale. The current "31-model" scope is a curated analysis subset, not the full data record. This mismatch needs a deliberate decision (see §7).

---

## 3. Scientific Framing / Narrative

### 3.1 Stated purpose
Systematically characterize a controlled collection of human OvCa cell models across three omics layers, to (a) confirm the data captures known subtype biology, (b) surface subtype-specific molecular features/biomarkers/targets, and (c) quantify **intra-subtype** heterogeneity to guide model selection. Thesis aims (verbatim structure): collect/preprocess → identify distinguishing features → evaluate within-subtype variation. Thesis **hypothesis**: models of the same subtype are more similar to each other than to other subtypes.

### 3.2 Argued significance
- Rare subtypes are understudied → HGSC-derived treatment paradigms translate poorly → need subtype-specific models and molecular maps.
- Models are made inconsistently across labs; 15–20% of cancer lines are misidentified → need systematic, controlled characterization.
- Prior lab precedent cited: systematic characterization of **syngeneic mouse** HGSC models (immunoregulatory variation drives tumor immune composition) — the human panel is framed as the analogous human-model effort.

### 3.3 Claimed novelty / "what's new"
1. Multi-omic (RNA+protein+WES) resource **enriched for rare subtypes** absent from CCLE.
2. **Concordant RNA–protein subtype signatures** ("more robust" markers), including "hundreds of novel markers."
3. A **novel EC–SCCOHT molecular convergence** attributed to shared SWI/SNF disruption (ARID1A loss in EC; SMARCA4 loss in SCCOHT).
4. **Intra-subtype heterogeneity maps** (esp. HGSC pathway-activity clusters) framed as a model-selection tool.
5. **Subtype-associated ADC-target expression** (translational hook: TROP-2/TACSTD2, FOLR1, ERBB2/HER2, SLC34A2, MSLN, NaPi2b, etc.).

> Note the tension: the *summary* explicitly promises the data "**can also reveal novel molecular patterns**" and to "**highlight diversity among models**" — discovery language that a data descriptor generally should avoid.

---

## 4. A Pivotal Preprocessing Decision: LGSC Excluded for Batch Effect

The thesis (§3.2) removed **all LGSC** models before analysis:
- LGSC formed a tight cluster fully separated from everything else in PCA/t-SNE.
- LGSC showed **uniform high activity across *all* PROGENy pathways** — biologically implausible (LGSC should be MAPK-high, inflammation-low).
- **LGSC was the only set not RNA-sequenced with the rest of the cohort** → subtype is **perfectly confounded with batch**.
- Conclusion: cannot separate technical from biological; exclude.

**My assessment:** The reasoning is sound and the "all pathways uniformly high" pattern is a classic technical signature — exclusion is the defensible conservative call **for a discovery analysis**. But several things deserve scrutiny:
- Because batch and subtype are perfectly confounded, this is genuinely unrecoverable by standard correction; however that also means we cannot *prove* it is purely technical. State it as a confound, not a verdict.
- **For a resource/data descriptor, deleting data outright is the wrong default.** The LGSC RNA-seq (and any LGSC proteomics/WES, which may not share the batch problem) should still be **deposited and documented** with the caveat, not silently dropped. Removing LGSC also removes one of the rare subtypes the paper claims to champion (acknowledged as "unfortunate" in the thesis limitations).
- **Open question:** were LGSC proteomics and WES also excluded, or only RNA-seq? The batch argument is RNA-seq-specific.

---

## 5. Results Claims — Claims Ledger

Evidence tiers: **Supported** (recapitulates known biology or robust QC) · **Suggestive** (plausible, under-powered / partially confounded) · **Unsupported/Questionable** (likely artifact, n=1, or over-interpreted). "One-vs-all DE" throughout = each subtype vs. all others pooled, with the pooled "other" dominated by HGSC (n=15) — the authors themselves flag this as a bias.

| # | Claim (theme) | Evidence given | Skeptical assessment | Tier |
|---|---|---|---|---|
| **QC / Technical Validation** |
| 1 | RNA-seq high quality: 45–97M reads (median 64.3M), pseudoalign >85% (median 91.1%) | Fig 1b/c boxplots; ENCODE >30M threshold | Legitimate, standard QC. Depth & alignment are genuinely fine. Fair as validation. | **Supported** |
| 2 | RNA–protein concordant, Pearson R = 0.34–0.46, "consistent with literature" | Fig 1d/e — but caption shows only **2 lines** (TOV21G, OV3331) | R≈0.4 is typical mRNA–protein agreement. But (a) is the range across all 31 or just illustrated on 2? (b) within-sample-across-genes vs across-samples-per-gene not specified. Reasonable QC point, mildly over-framed as validating "suitability." Verify computation in code. | **Supported** (as QC) |
| 3 | "We systematically generated RNA-seq, TMT proteomics, **and WES** from **31 cell models**" | Table 1 / Fig 1a | **Over-claim.** WES/CNV covers **~13 lines** (Fig 3c: 9 HGSC+2CC+1EC+1MC); proteomics N unstated. Hedged elsewhere ("among the sequenced whole-exome lines") but topline says 31. Must be corrected to actual per-assay N. | **Unsupported (as written)** |
| **WES / Genomics** |
| 4 | All HGSC have TP53 mutations | Fig 1f/3d | Expected (~96% in HGSC); good positive control. | **Supported** |
| 5 | "Majority" of HGSC have BRCA2 (thesis: **5** of 9 sequenced); BRCA1 also found | Fig 3d | 5/9 ≫ population BRCA2 rate (~10–15%). Without matched normal, likely inflated by germline/unfiltered variants. Directionally consistent with HRD but the *rate* is not credible as somatic. | **Suggestive/Questionable** |
| 6 | **ATM mutated in 100% of HGSC, ATR in 75%** | Fig 3d | **Almost certainly artifact.** Somatic ATM/ATR ~2% in HGSC. Thesis caveats (no matched normal → germline). **Manuscript reports it with NO caveat.** Do not present as biology. | **Unsupported** |
| 7 | Non-HGSC lines (2 CC, 1 MC, 1 EC) "all exhibited KRAS mutations, linking to MAPK" | Fig 3e/1g | KRAS is canonical for **MC only**. KRAS in CC (ARID1A/PIK3CA-driven) and EC (CTNNB1/ARID1A/PIK3CA) is **unexpected**. Thesis flags this ("only known…in Mucinous"); **manuscript drops the caveat** and frames the anomaly as validation. Likely germline/unfiltered or genuinely atypical lines. | **Questionable** |
| 8 | Both CC lines have ARID1A mutations | Fig 3e | Canonical for CC (~46–58%). Good positive control. | **Supported** |
| 9 | EC line TOV112D has SMARCA4 mutation ("typically SCCOHT") → basis for EC–SCCOHT link | Fig 3e | n=1 line. Interesting but anecdotal; anchors an over-built convergence story (claim #17). | **Suggestive** |
| 10 | HGSC show widespread CNV / genomic instability; 3q amp in all HGSC except OV90; recurrent 11p deletions | Fig 3c | Global CNV burden & 3q gain are canonical HGSC (3q ~50% of cases) → supportive. "All except OV90" from ~9 lines. 11p-loss "recurrent/uncharacterized" is more novel → needs CN-calling scrutiny (segmentation, GSE85671 normal reference adequacy). | **Supported** (3q/instability) / **Suggestive** (11p) |
| **Subtype separation (core validation)** |
| 11 | Subtypes separate by Spearman clustering, PCA, t-SNE (RNA & protein) | Fig 1f, 2, 4–6 | Broadly holds, esp. HGSC. This is the legitimate "recapitulates biology" core. But driven by tiny rare-subtype n and no batch modeling beyond LGSC removal. | **Supported** (HGSC) / **Suggestive** (rare) |
| 12 | Silhouette: MMMT/SCCOHT/HGSC/EC positive (0.3–0.8); MC & CC ≈0 (high intra-subtype variability) | Fig 2c | Fair *descriptively*. But silhouette on n=2 groups (MMMT/SCCOHT/EC) is statistically unstable; high values partly reflect tiny, tight n. | **Suggestive** |
| 13 | Proteomic silhouette: EC & SCCOHT **negative** (group with other subtypes more than own) | Fig 2f (thesis) | **Undercuts** the clean-separation narrative for the very subtypes used in the convergence story (#17). Honestly reported in thesis. | **Supported** (as a caveat) |
| 14 | Outliers: OV90 (HGSC) resembles no subtype / clusters with mucinous; TOV2414 (MC) resembles HGSC/CC; VOA6861 (CC) lacks signature | Fig 2a/e, 3a, 4, 6 | Real and consistent across analyses. Interpretation "misannotation vs phenotypic divergence" is appropriate as hypotheses. **Caution:** OV90 is a widely STR-authenticated HGSC line — "misannotation" is a strong claim; check sample identity/STR & possible sample swap before publishing. | **Supported** (pattern) / **Suggestive** (cause) |
| **Signatures / DE / GO** |
| 15 | One-vs-all DE yields subtype signatures (padj<0.05, logFC>1, top 100); expression "consistent across models of each subtype" | Fig 3a, 7, 8 | Signatures exist and look block-structured. But top-100-by-p from one-vs-all with HGSC-dominated "other" and n=2 rare groups → many "signature" genes may be **line-specific, not subtype-specific**. Authors acknowledge the one-vs-all bias. | **Suggestive** |
| 16 | Concordant RNA–protein signatures reproduce known biomarkers (GGT1, GDA for CC; RBP1, WFDC2/HE4 for HGSC) **and hundreds of novel markers** | Fig 3b, 10 | Reproducing known biomarkers = good validation (**Supported**). "**Hundreds of novel markers**" and single-gene stories (e.g., **C4BPB → "novel immunoregulatory mechanism of secreted mucus"** in MC, n=3) are **over-reach** — unreplicated, possibly one-line-driven, no functional data. | **Supported** (known) / **Unsupported** (novel) |
| 17 | GO over-representation recapitulates known biology: MC glycosylation/oligosaccharide (B3GALT2, FUT3, ST8SIA4, TREH); HGSC proliferation (nucleosome assembly) + DSB repair (RAD1, XRCC2, PARP1, TP53, MCM); CC glutathione catabolism | Fig 9 | Good "recapitulates known biology" validation. Directionally solid. Layered speculation (hypoxia retention, complement immune-modulation) goes beyond the data. | **Supported** (as validation) |
| 18 | **EC–SCCOHT molecular convergence** via shared SWI/SNF disruption (ARID1A loss EC; SMARCA4 loss SCCOHT) | PCA proximity, reciprocal signature overlap, TOV112D SMARCA4 (Fig 16, thesis) | **n=2 vs n=2**, and proteomic silhouette says these very samples don't even separate cleanly (#13). Thesis correctly labels "hypothesis-generating"; manuscript states more boldly. Real chromatin-biology rationale, but evidence is minimal. | **Suggestive** (thesis) / **over-stated** (ms) |
| 19 | ADC targets show subtype-associated expression: TACSTD2/TROP-2 high in HGSC; FOLR1 high in MC; ERBB2/HER2 high in CC/EC (esp. TOV3392D); SLC34A2 high in VOA12539 | Fig 3d(ms)/11 | Useful **descriptive** resource utility ("pick a model"). But several are **n=1** line observations dressed as findings, and this is the main lever pulling the paper toward translational-discovery framing. | **Suggestive** (utility) |
| **Intra-subtype heterogeneity** |
| 20 | HGSC pathway-activity clusters (PROGENy/singscore): growth-factor/RTK cluster (OV1946, OV866, OV1369.R2, OV2295.R2, TOV3041G); OV4453 = inflammatory (NFκB/JAK-STAT/WNT); low-signaling cluster; chemoresistance high in OV1946/OV90 (low E2F/proliferation) | Fig 4/12 | Descriptively reasonable, but Z-scored gene-set scores + Ward clustering on n=15 are **unstable and gene-set-dependent**; no parameter-sensitivity or stability analysis. "Chemoresistance signature" is a transcriptional proxy, not measured drug response. | **Suggestive** |
| 21 | Genotype–phenotype concordance anecdotes: OV2295_R2 PIK3CA mut + high PI3K activity; "OV3313R" SMARCA4 frameshift + low WNT | Fig 3d + 12 | **n=1 cherry-picks** presented as validation; no across-line correlation test. **"OV3313R" is not in the 31-line roster** (roster has OV3331 and OV3133-R) → probable **line-labeling/typo error** to resolve. | **Unsupported** (as evidence) + **data-integrity flag** |
| 22 | ConsensusOV maps HGSC lines to TCGA subtypes (IMR/DIF/PRO/MES) | Table 2, Fig 13 (**thesis only**) | Applying a **bulk-tumor** classifier (heavily driven by stroma/immune) to **pure tumor-cell cultures** is of dubious validity; marker genes only weakly corroborate (CXCR3 absent, MES markers inconsistent). Thesis caveats this well. **Wisely dropped from the manuscript.** | **Suggestive/Low** |
| 23 | MC & CC per-line PROGENy narratives (TOV2414 HGSC-like; VOA8762 inflammatory; CC hormone-signaling patterns; TOV3392D Notch-high "stem-like") | Fig 14/15 (thesis) | Deep n=1 storytelling on 3 MC / 7 CC lines. Interesting for model selection but not generalizable. Largely relegated to thesis. | **Suggestive/Low** |
| **Overall** |
| 24 | Hypothesis "same-subtype models more similar than cross-subtype" is **partially** supported | All of the above | Honest, calibrated conclusion. Appropriate given confounds/n. | **Supported** (as stated, i.e., "partial") |

---

## 6. Methods (as described in prose — for later code comparison)

Capture these so the re-analysis can verify against the actual scripts:

**RNA-seq processing.** kallisto pseudoalignment → GRCh38 transcriptome; `tximport` **v1.36.1** → gene-level; PCA on **log-transformed TPM**; genes with all-zero expression removed; subtype-specific PCA uses only that subtype's samples; **Euclidean distances from first 10 PCs (~80% variance)**.

**Differential expression.** `DESeq2` **v1.48.2**; **one-vs-all** (subtype vs "Other"); Wald test; signature = padj<0.05 **&** log2FC>1, then **top 100 by padj**. (Note: heatmap caption says |log₂FC|>1 but signatures are described as up-regulated only — clarify one-sided vs two-sided.) *Concern: running DESeq2 on **TPM/kallisto-tximport** vs raw counts — verify the count matrix fed to DESeq2 is appropriate (DESeq2 expects counts, not TPM).*

**GO ORA.** `pcaExplorer::topGOtable` **v3.2.0** wrapping `topGO` **v2.60.1**, `org.Hs.eg.db`; **elim** algorithm; run per subtype on its signature.

**Gene-set scoring / signaling.** `singscore` **v1.28.1** on log-TPM, Z-scored across samples; `PROGENy` **v1.30.0** (14 pathways), **top 500 genes per model**; heatmaps via `ComplexHeatmap` **v2.24.1**, **Ward** clustering on columns.

**Proteomics.** Morin lab (GSC, BC Cancer): SP3 on paramagnetic beads → trypsin → **TMT** isobaric labels → offline high-pH RP fractionation → LC-MS/MS on Orbitrap with **SPS-MS3**; searched vs UniProt human; quant from TMT reporter ions. Protein-side DE = **unpaired t-test**, same thresholds (padj<0.05, logFC>1). *No proteomics-specific QC described (normalization, missingness, batch/plex structure, ratio compression) — a gap for a descriptor.*

**WES.** Nimblegen SeqCap EZ Exome 3.0; HiSeq 2500, ~50M reads, 100bp PE; **Sarek** Nextflow (CNVkit, Mutect2, VEP). **No matched normal**; **GSE85671 (Lai et al.) normals used as CNV reference.** *No matched normal is the root cause of the inflated mutation frequencies (#5–7).* Recon indicates "old" and "new" WES versions exist — determine which was used.

**Not described in Methods but present in Results:** ConsensusOV (thesis), the specific ADC-target list and its source, silhouette computation, t-SNE parameters, and the gene-set catalog used for Fig 4/12. Cell-line culture/derivation and RT-qPCR headers are present but empty ("Cell lines"/"RT-qPCR" are stubs).

---

## 7. Weaknesses, Over-claims & Gaps (skeptical pass)

**A. Framing mismatch with a Data Descriptor (biggest issue).** The narrative is discovery-driven. *Scientific Data* descriptors are judged on **data generation, quality, technical validation, and reuse potential**, and reviewers explicitly discourage extended interpretation/discovery. Claims of "novel biomarkers," "novel immunoregulatory mechanism," EC–SCCOHT convergence, chemoresistance biology, and ADC discovery belong in a research article. As written, this would likely be criticized as "a research paper mislabeled as a descriptor."

**B. Missing descriptor machinery.** Outline still shows **"Data Records → To do: data availability,"** **"Proteomics QC — ???,"** **"WES QC — ???,"** and no code-availability statement. These are **mandatory** for the venue. Package versions are commendably tracked (good start on reproducibility), but per-assay technical validation is thin (proteomics/WES essentially unvalidated in-text).

**C. Statistical power / pseudoreplication.** n=2 for EC, SCCOHT, MMMT; n=3 MC. Silhouette, DE, and "convergence" on n=2 are fragile. One-vs-all DE with an HGSC-dominated reference biases every non-HGSC signature. No donor/passage structure modeled (**passage number** — the pXX suffixes — is never used as a covariate, and some lines have R/R2 "recurrent"/re-derived versions, e.g., OV2295 vs OV2295-R2, that are treated as independent).

**D. WES artifacts stated as biology.** #6 (ATM 100%/ATR 75%), #5 (BRCA2 majority), #7 (KRAS in CC/EC) are the clearest cases where **technical limitation (no matched normal + permissive calling) is narrated as findings**, especially in the manuscript where caveats present in the thesis were dropped. This is the highest-risk credibility problem.

**E. Incomplete/uneven assay coverage not disclosed up front.** WES ≈13/31; proteomics N unstated; LGSC dropped; broader generated sample set undocumented. A descriptor must state exactly what exists for each sample × assay (the Table 1 availability grid is the right instinct — it just needs to be accurate and complete, including excluded/replicate samples).

**F. Single-line stories.** C4BPB/complement (MC), TOV3392D HER2, VOA12539 SLC34A2, TOV3392D Notch "stem-like," OV4453 inflammatory — repeatedly, **one line = a claim**. Fine as "this resource lets you find a model with property X"; not fine as subtype biology.

**G. Data-integrity flags to chase in code:** (i) "OV3313R" (not in roster) in manuscript §"Expression variability"; (ii) line-name inconsistencies (OV3133 vs OV3133-R; TOV1369 vs OV1369-R2/TOV1369); (iii) DESeq2 apparently fed log-TPM rather than counts — verify; (iv) whether Pearson 0.34–0.46 is computed across all samples or only the 2 shown; (v) TOV3392D is listed as **clear cell** in the roster but discussed with CC — confirm subtype label source.

**H. Causal/over-strong language.** "confirm validity," "these findings confirm," "supporting their suitability as ideal preclinical models" — calibrate to "consistent with," "recapitulates," "suggests."

**What is genuinely solid:** RNA-seq depth/alignment QC; subtype separation for HGSC; recovery of canonical biomarkers and canonical GO terms; canonical HGSC genomic instability & 3q gain; honest limitations section in the **thesis**; commendable version tracking.

---

## 8. Resource-Paper Fit (*Scientific Data* Data Descriptor)

**What the venue rewards:** clear description of *what was measured and how*; rigorous *technical validation of data quality*; complete *Data Records* with repository accessions and file-level detail; *reusability* (metadata, code, usage notes). Interpretation/biological discovery is explicitly de-emphasized.

**Where this stands:**
- **Fits well:** cohort table + per-assay workflow diagrams; RNA-seq QC; "recapitulates known biology" checks (subtype separation, known biomarkers/GO, canonical HGSC genomics) — these are ideal *technical validation* content.
- **Needs to be added/finished:** Data Records + accessions (GEO/SRA for RNA-seq & WES; PRIDE/MassIVE for proteomics); **proteomics QC** (identifications, missingness, normalization, plex/batch, ratio compression) and **WES QC** (coverage, on-target, contamination, and an honest statement on no-matched-normal + how variants were filtered); **code availability**; **usage notes** (recommended: keep them — they aid reuse and are where model-selection guidance belongs).
- **Needs to be removed or heavily demoted to "validation" framing:** novel-marker discovery, EC–SCCOHT convergence as a finding, chemoresistance/pathway biology, ConsensusOV, per-line PROGENy narratives, ADC "discovery." Retain ADC-target expression only as a **reuse example** ("the data can be queried for target expression"), not a claim.
- **Scope decision required:** document the **full generated dataset** (all lines × assays, incl. LGSC and replicates, with QC-based inclusion/exclusion flags), not just the 31-model analysis subset.

**Bottom line on fit:** The *data* are a legitimate, useful resource and a good match for the venue. The *manuscript* currently reads as a discovery thesis. The re-scope is primarily **subtractive** (strip discovery claims), **additive** (Data Records, per-assay validation, availability/code), and **corrective** (fix the WES-coverage over-claim and the no-matched-normal mutation artifacts).

---

## 9. Thesis vs. Manuscript — Framing & Scope Differences

**Title/frame.** Thesis: "A multi-omic **analysis of the heterogeneity** of ovarian carcinoma" (discovery thesis, full Intro→Aims→Hypothesis→Results→Discussion). Manuscript: "A multi-omic **resource** of human ovarian cancer cell models" (repackaged toward *Scientific Data*). The word changed; the content largely did not.

**Content dropped from thesis → manuscript (good editorial instincts):** ConsensusOV/TCGA mapping (Table 2, Fig 13); dedicated EC–SCCOHT figure (Fig 16); per-line MC/CC PROGENy deep-dives (Fig 14/15); volcano plots (Fig 8); most per-subtype GO narrative detail; explicit ATM/ATR **caveat** (dropped — this is a *bad* loss; see below).

**Content retained/promoted:** QC, RNA–protein concordance, PCA/t-SNE/silhouette, one-vs-all signatures + concordant RNA-protein signature, GO validation, HGSC gene-set heterogeneity, ADC targets. WES validation (TP53/BRCA/3q) is **more prominent** in the manuscript's Technical Validation.

**Caveat erosion (concerning).** The manuscript **drops** the thesis's honest hedges: ATM/ATR "may include germline" caveat is gone; the KRAS-in-CC/EC "only known in mucinous" caveat is gone. The thesis is consistently more calibrated than the manuscript.

**Older (Nov 6) → latest (Dec 4) manuscript drift:**
- LGSC removed from the subtype list in the Overview (Nov 6 still listed "Low-grade serous" as included — inconsistent with its exclusion); Dec 4 gives clean per-subtype n.
- Author affiliations complete in Nov 6 (5 affiliations incl. GSC, CRCHUM, BC Cancer Molecular Oncology); truncated to 2 in Dec 4 (likely formatting loss — restore).
- Dec 4 **adds**: explicit silhouette values; OV90/TOV2414/VOA6861 signature-discrepancy discussion; the WES validation paragraph (incl. uncaveated ATM/ATR 100%/75%); C4BPB/complement "novel mechanism"; EC–SCCOHT SWI/SNF paragraph; ADC ERBB2/TACSTD2 line; genotype–phenotype anecdotes (OV2295_R2, "OV3313R").
- Dec 4 **narrows** the heterogeneity section to HGSC (MC/CC → supplemental), per the outline.
- **Net drift:** even as it is relabeled a "resource," the Dec 4 draft leans *further* into discovery/novelty than Nov 6.

---

## 10. Top Decisions / Questions for the PI

1. **Descriptor vs. research article.** Commit to the *Scientific Data* descriptor frame and **strip discovery claims** to a "recapitulates known biology" validation layer? Or split into (a) a lean descriptor + (b) a separate discovery paper that can carry the EC–SCCOHT / heterogeneity / ADC stories with proper power?
2. **WES scope & the no-matched-normal problem.** WES covers only ~13 lines and lacks matched normals. Options: (a) reframe WES as partial/targeted validation and **remove mutation-frequency claims** (ATM/ATR/BRCA2/KRAS); (b) re-call variants with a panel-of-normals / population-frequency filtering; (c) sequence the remaining lines. Which?
3. **Cohort definition for the record.** Document the **full generated dataset** (all lines × assays, incl. LGSC + replicates + the extra roster entries, with QC-based include/exclude flags), or publish only the curated 31? The venue expects the former.
4. **LGSC.** Deposit + document the excluded LGSC (with the batch caveat) rather than delete it? Were LGSC proteomics/WES also excluded, or only RNA-seq (the batch argument is RNA-seq-specific)?
5. **OV90 identity.** Given OV90 (a canonical HGSC line) fails to resemble HGSC across every analysis — pursue **STR/identity verification** / possible sample swap before implying "misannotation"?
6. **Proteomics validation.** Who supplies the proteomics QC + normalization/plex details (Morin lab)? Needed for the Technical Validation section and PRIDE/MassIVE deposition.
7. **Data-integrity clean-up.** Resolve line-name errors ("OV3313R"; OV3133 vs OV3133-R; TOV1369 vs OV1369-R2) and confirm DESeq2 input is counts (not TPM) as part of the re-analysis.

---

*Report file:* `/Users/dpcook/Analysis/ovcan_human/reports/assessment/01_thesis_manuscript_synthesis.md`
