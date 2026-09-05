# Analysis Plan: OvCAN Human Ovarian Cancer Cell-Line Multi-Omic Resource

**Last Updated:** 2026-07-23
**Status:** Planning complete → beginning re-analysis
**Anchors:** `PROJECT_SPEC.md` · `reports/00_synthesis_and_recommendations.md` · `metadata/samples.csv`

## Objective
Produce a clean, reproducible re-analysis of the generated multi-omic data supporting a *Scientific Data* Data Descriptor: technical validation, cell-line authentication, recapitulation of known subtype biology, and a subtype-resolved ADC-target usage example. Correct the artifacts and analysis bugs identified in the archive assessment; drop external (LGSC/Carey) data.

## Input Data
- **RNA-seq:** kallisto abundances, `judy_archive/data/rna_seq/` (re-quantify from raw where needed). Generated lines only.
- **Proteomics:** TMT tables, `judy_archive/data/proteomics/` (+ pipeline/QC to be obtained from Morin lab).
- **WES:** Mutect2 MAFs + CNVkit segments, `judy_archive/data/wes*/` (23 lines).
- **Reference:** MSigDB Hallmark v7.4; PROGENy; ConsensusOV (exploratory only); gnomAD + PoN (WES filtering); scarHRD.
- **Sample sheet:** `metadata/samples.csv` (single source of truth; `provenance` flag separates generated vs external).

## Analysis Phases

### Phase 0 — Setup ✅ (near-complete)
Project scaffold, sample sheet, analysis log. **Remaining:** `renv` lockfile + pinned annotation/gene-set versions; `scripts/00_setup.R`.

### Phase 1 — RNA-seq re-processing
**Goal:** Trustworthy expression matrices + subtype-separation validation.
**Approach:** tximport → **counts** (verify, not TPM) → DESeq2; pinned Ensembl annotation (no live biomaRt); drop external LGSC; PCA/clustering/silhouette with n annotated and site modeled as covariate/diagnostic; recompute RNA–protein correlation across all lines (benchmark to CPTAC 0.38/0.45); restrict formal DE to HGS (n≈15) or present descriptively; **fix the CC/MC signature `Gene` bug**.
**Success:** subtypes separate as expected; markers/GO recover known biology; correlation in the 0.4 range; runs top-to-bottom.
**Output:** `output/` matrices; QC + separation figures.

### Phase 2 — Proteomics (document + batch-aware re-analysis; NOT full re-processing)
**Constraint (PI, 2026-07-23):** the TMT data were processed by the Morin lab; we have the processed relative-abundance + peptide tables and, separately, their **CV and peptide-coverage distributions** — but not the raw MS files or a reproducible upstream pipeline. So we *document* (not reconstruct) the processing.
**Goal:** Documented, batch-aware protein matrix + the QC we can actually assemble.
**Approach:** describe the TMT workflow from the methods provided; assemble Technical-Validation QC from (a) Morin's CV + peptide-coverage distributions and (b) our own computed metrics — protein/peptide counts, missingness pattern, **bridge-channel replicate correlation**, PCA/clustering separation; model plex + site batch; replace `na.omit` listwise deletion with principled missingness handling; recompute PCA/clustering/markers.
**Success:** proteomics QC assembled and batch structure addressed from available data.
**Output:** documented protein matrix; proteomics QC figure.
**Open (collaborator; PI handling):** raw MS files + fuller pipeline/QC from Morin for PRIDE/MassIVE deposition — pending; note as a limitation if unavailable.

### Phase 3 — WES re-analysis
**Goal:** Defensible genomics under the no-matched-normal constraint.
**Approach:** re-filter Mutect2 tumor-only calls (gnomAD population-AF + panel-of-normals); restrict oncoplots to well-supported canonical drivers; compute **genomic HRD (scarHRD)**; CNV QC + capture-kit concordance check for the 5 public normals; resolve TOV3121D MAF + index/compression issues.
**Success:** canonical drivers recovered (TP53 in HGSC); implausible germline-driven rates removed; HRD reported as genomic (not expression); limitation stated.
**Output:** filtered variant tables; CNV + oncoplot + HRD figures.

### Phase 4 — Authentication (new, high value)
**Goal:** Establish identity + histotype fidelity.
**Approach:** compile STR (request missing — VOA8762/8771 etc.); histotype-authenticate from CNV + mutation + expression (Stordal approach); independently test COV434→SCCOHT and TOV112D→dedifferentiated (SMARCA4/SMARCA2 from RNA+protein); mucinous CK7/SATB2/PAX8/WT1 + cluster vs public MOC benchmarks (MCAS/RMUG-S/COV644/OV-90).
**Success:** each line's identity/histotype documented; reclassifications independently recovered.
**Output:** authentication table + figure.

### Phase 5 — Validation & usage analyses
**Goal:** The descriptor's core evidence + usage example.
**Approach:** subtype-recapitulation (separation + marker/GO recovery); RNA–protein concordance vs CPTAC; **ADC-target expression atlas** (RNA + protein, by subtype — the featured usage example); within-HGSC pathway heterogeneity as a model-selection example.
**Success:** figures 3–5 (below) complete and quantified.
**Output:** `docs/manuscript/figures/`.

### Phase 6 — Deposition & writing
Deposit to GEO/SRA (RNA-seq, WES) + PRIDE/MassIVE (proteomics); assemble Data Records; write the descriptor; finalize figures. **Deposition hazard:** the archived CNVkit outputs bundle third-party public normal exomes (SRR4039087/88/89/96/97 + a 461 MB BAM, from PRJNA339046) — strip/replace before deposit, or deposit only our tumors' derived `.cnr/.cns`, so we don't redistribute others' data. Cite (don't deposit) external analysis references (MSigDB Hallmark, Han metagenes, Peng list, TCGA/ConsensusOV).

## Key Parameters & Decisions
| Parameter / Decision | Value | Rationale |
|---|---|---|
| Framing | Lean Data Descriptor + ADC usage example | Venue rewards data quality/reuse; PI decision 2026-07-23 |
| WES SNV | Re-filter (gnomAD+PoN), canonical drivers only | Tumor-only calling → germline FP (~50–70% FDR); PI decision |
| Genomic HRD | scarHRD (WES) | HRDetect/CHORD are WGS-only; Peng signature is expression, not HRD |
| LGSC/Carey | Removed as external | Not generated here (published Carey/OVCARE); PI decision |
| DE model | DESeq2 on counts; HGS-restricted or descriptive | n=2 rare-subtype groups can't support generalizable DE |
| Subtyping (ConsensusOV) | Supplement only, caveated | TCGA subtypes are TME-driven; invalid on pure lines |
| RNA–protein corr benchmark | CPTAC ovarian 0.38/0.45 | Same-disease/analyte benchmark |
| Batch covariates | site (Mes-Masson/Huntsman), TMT plex, passage | Tracked in `samples.csv`; modeled/diagnosed |

## Proposed manuscript figures (5 main)
1. Cohort & design (provenance, subtype composition, per-assay coverage grid, workflow).
2. QC / technical validation (RNA depth+alignment; proteomics counts/missingness/bridge reproducibility; WES coverage/on-target).
3. Recapitulates known biology (PCA/clustering/silhouette RNA+protein; marker/GO recovery; RNA–protein concordance vs CPTAC).
4. Genomics & authentication (CNV landscape; canonical drivers w/ caveat; STR + histotype + reclassification recovery; genomic HRD).
5. Usage example (subtype-resolved ADC-target atlas; within-HGSC heterogeneity).
*(Discovery panels → supplement/companion.)*

## Known Limitations
- WES is tumor-only (no matched normal); somatic SNV frequencies unreliable → canonical-driver focus + explicit caveat.
- Proteomics were externally processed (Morin lab); only CV + peptide-coverage QC and processed tables are available (no raw MS / reproducible pipeline). Technical Validation uses those + our computed downstream metrics; raw-MS deposition is pending Morin (PI handling).
- Multi-omic overlap is 13 lines; per-assay coverage varies — report honestly.
- Rare subtypes are small n (EC likely n=1 after TOV112D; SCCOHT/MMMT n=2) — descriptive, not powered.
- Passage differs across assays for most lines — stated per line per assay.
- VOA8762/8771 provenance unresolved pending BC Cancer data.
