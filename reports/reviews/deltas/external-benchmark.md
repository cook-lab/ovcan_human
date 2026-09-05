# Delta — External identity & benchmarking (CCLE/DepMap, ConsensusOV, Cellosaurus/STR)

**Workstream:** `external-benchmark` · **Script:** `scripts/18_external_benchmarking.R` (runnable top-to-bottom; sources `00_setup.R`)
**Addresses:** peer-review §4 essential #5 (external identity/technical validation) and §3.2.
**Date:** 2026-07-23

## TL;DR
- **5 of our 42 lines exist in CCLE/DepMap** (OV90, TOV21G, TOV112D, BIN67, COV434) — verified by DepMap `Model.csv` **and** Cellosaurus RRID. The 11 BC Cancer/OVCARE **VOA lines are absent from both DepMap and Cellosaurus**; only 3 Mes-Masson lines are in DepMap.
- **Every overlapping line self-matches** its DepMap namesake: Spearman **0.74–0.88** on 2000 HVGs, **rank 1 of 67** DepMap ovarian lines, **reciprocal-best in both directions**. Strong cross-lab, cross-platform (kallisto vs RSEM) identity confirmation.
- **Driver cross-check concordant:** DepMap independently recovers TOV21G hypermutation (568 damaging vs 7–18), TOV112D TP53-R175H + SMARCA4-damaging, all 5 TOV21G clear-cell drivers, OV90 SMAD4 — and **independently shows SMARCA4-damaging mutations in COV434 and BIN67** (lines we have no WES for), reinforcing the SCCOHT calls.
- **Cellosaurus independently flags COV434 as "Misclassified … SCCOHT"** (Karnezis 2021, PubMed 33328126) — external corroboration of our reclassification. TOV112D is **not** yet flagged (still listed endometrioid).
- **30/42 lines have a documented STR profile in Cellosaurus** (originator-deposited); the 11 VOA lines + TOV3121D have **no public STR record**. **No in-house STR or mycoplasma testing was performed** (STR_status=NA for all 42) — draft statements below.
- **ConsensusOV** calls surfaced for 15 HGS lines (DIF 6 / MES 4 / IMR 3 / PRO 2). Report **with the caveat** that IMR/MES are TME/stroma-driven and of dubious validity on pure tumor cultures.

**Figure:** `figs/f_external_concordance.pdf` + `reports/assets/f_external_concordance.png`
**Outputs:** `output/external_ccle_concordance.csv`, `output/consensusov_calls.csv`, `output/cellosaurus_str_status.csv`
**External data:** `output/external/` (DepMap Public 24Q4; Cellosaurus API JSON cache) — provenance at end.

---

## Task 1 — CCLE/DepMap expression concordance

### 1a. Overlap set (verified identifiers — none invented)
Determined by matching stripped cell-line names against DepMap `Model.csv`, then confirming each RRID against Cellosaurus. The overlap is small **and honest**: only 3 Mes-Masson lines (OV90, TOV21G, TOV112D) plus BIN67 and COV434 are in DepMap; **no VOA/BC Cancer line and no other Mes-Masson line is present.**

| Our line | Our subtype | DepMap ModelID | RRID (Cellosaurus) | DepMap Oncotree subtype |
|---|---|---|---|---|
| OV90 | HGS | ACH-000291 | CVCL_3768 | Serous Ovarian Cancer |
| TOV21G | CC | ACH-000885 | CVCL_3613 | Clear Cell Ovarian Cancer |
| TOV112D | EC | ACH-000048 | CVCL_3612 | Endometrioid Ovarian Cancer |
| BIN67 | SCCOHT | ACH-001278 | CVCL_S987 | **Small Cell Carcinoma of the Ovary** |
| COV434 | SCCOHT | ACH-000123 | CVCL_2010 | **Small Cell Carcinoma of the Ovary** |

> DepMap's own Oncotree labels COV434 and BIN67 as *Small Cell Carcinoma of the Ovary* (SCCOHT) — a third-party endorsement of the SCCOHT calls the report makes.

### 1b. Expression identity (the headline result)
Our RNA (linear TPM → collapsed to Entrez → log2(TPM+1)) vs DepMap `OmicsExpressionProteinCodingGenesTPMLogp1` (log2(TPM+1)); Spearman over the top-2000 most-variable shared genes (19,062 shared Entrez genes). Each overlap line was correlated against **all 67 DepMap ovarian lines**.

| Our line | Self Spearman (HVG) | Self Pearson | Rank of correct DepMap line (of 67) | Runner-up (subtype) | Margin self−next | Reciprocal best? |
|---|---|---|---|---|---|---|
| OV90 | **0.879** | 0.878 | **1** | OVCAR5 (0.445) | 0.434 | ✅ |
| TOV21G | **0.805** | 0.811 | **1** | OVTOKO (0.574) | 0.231 | ✅ |
| TOV112D | **0.839** | 0.844 | **1** | A2780 (0.635) | 0.204 | ✅ |
| COV434 | **0.821** | 0.841 | **1** | BIN67 (0.683) | 0.139 | ✅ |
| BIN67 | **0.739** | 0.721 | **1** | COV434 (0.703) | 0.036 | ✅ |

Interpretation: all 5 lines are their DepMap namesake's top match out of 67 ovarian lines **and** vice-versa. The cross-platform ceiling (kallisto/tximport here vs RSEM in DepMap, different labs/passages) caps absolute values well below 1; the **specificity** (self ≫ all others) is the identity signal. **Honest nuance:** BIN67 and COV434 are transcriptional near-twins (both SCCOHT/SMARCA4-null) — BIN67's self-margin over COV434 is only 0.036 — but self still wins for each, and the driver layer (below) separates them.

### 1c. Driver-mutation cross-check (DepMap somatic-mutation matrices vs our WES)
DepMap hotspot/damaging-mutation matrices (24Q4). Our WES exists for OV90, TOV21G, TOV112D only; COV434/BIN67 are DepMap-only checks.

| Gene | Concordance |
|---|---|
| **Hypermutation (TOV21G)** | DepMap **568 damaging** mutations vs 7–18 for the other four; our WES **1,416 PASS variants** vs 133–240. Both independently flag TOV21G (MSI-high) as the hypermutator. |
| **TP53** | OV90 (DepMap hot/dam=2/2; our p.S215R VAF 0.99) and TOV112D (DepMap 2/2; our **p.R175H** VAF 0.95) both TP53-mutant. TOV21G, BIN67, COV434 = TP53-**wildtype** in DepMap — correct for clear-cell and SCCOHT. |
| **SMARCA4** | TOV112D DepMap damaging=2 ↔ our truncating **p.L639X**. **COV434 damaging=2 and BIN67 damaging=1 in DepMap** (we have no WES for either) — independent genomic support for SWI/SNF loss underpinning the SCCOHT calls. |
| **SMARCA2** | 0 damaging in all 5 (not in hotspot matrix) — consistent with SMARCA2 loss being **epigenetic silencing**, not mutation (matches our mRNA-silencing finding). |
| **TOV21G clear-cell drivers** | 5/5 concordant: PIK3CA (DepMap hotspot ↔ our **H1047Y**), KRAS (↔ **G13C**), CTNNB1 (↔ **A5V**), ARID1A (DepMap damaging ↔ our 2 frameshifts), PTEN (↔ **L265X**). |
| **OV90 SMAD4** | DepMap hot/dam=2/2 ↔ our nonsense **p.R445\*** — corroborates the co-occurring SMAD4 lesion. |

**Minor, honestly-flagged discordances:** (i) our TOV112D **KRAS p.A59T** (VAF 0.36) is *not* called by DepMap — a rare, likely subclonal/low-confidence call in our tumor-only WES; (ii) our OV90 **BRAF in-frame deletion** is not in DepMap's hotspot/damaging matrices (those do not capture in-frame indels), so it is neither confirmed nor refuted there. Neither affects identity.

---

## Task 2 — ConsensusOV calls (surfaced, with critical caveat)

Extracted from `samples.csv notes` → `output/consensusov_calls.csv` (all 42 rows; call present for the 15 HGS RNA lines; the 9 WES-only HGS lines and all non-HGS lines are NA/not-applicable — ConsensusOV is an HGSC classifier).

| ConsensusOV | n (of 15 HGS RNA lines) | Lines |
|---|---|---|
| DIF (differentiated) | 6 | OV1369-R2, OV2295, OV3133-R, OV4485, TOV1369, TOV3133G |
| MES (mesenchymal) | 4 | OV1946, OV866-2, TOV3041G, TOV3291G |
| IMR (immunoreactive) | 3 | OV2085, OV2295-R2, OV90 |
| PRO (proliferative) | 2 | OV4453, OV3331 |

**CAVEAT (must accompany the table):** the TCGA/ConsensusOV HGSC subtypes are substantially **microenvironment-driven** — *immunoreactive (IMR)* reflects immune infiltrate and *mesenchymal (MES)* reflects stromal/CAF content (Verhaak 2013; Chen 2018 consensusOV). Our lines are **pure tumor-cell cultures with no TME**, so an IMR or MES label cannot carry its usual biological meaning. That 7/15 HGS lines are called MES or IMR — subtypes defined by compartments these cultures lack — is itself evidence the labels should be read as an **orthogonal external annotation, not ground truth**. Recommend presenting the calls transparently as a provenance/legacy annotation with this caveat, not as a validated subtype axis. (This complements the report's own within-HGSC pathway-activity strata in F5, which are tumor-cell-intrinsic.)

---

## Task 3 — Cellosaurus accessions + STR / mycoplasma documentation

Full table: `output/cellosaurus_str_status.csv` (per line: accession, identifier, in-Cellosaurus, STR documented, n STR markers, problematic flag, match confidence, RRID↔DepMap check). All 30 resolved accessions are **high-confidence exact name matches**; the 5 DepMap lines' accessions equal their DepMap RRID (cross-checked TRUE).

**Coverage:**
- **In Cellosaurus with a documented (originator-deposited) STR profile: 30/42.** All Mes-Masson OV/TOV lines resolved (e.g., TOV2414 = **CVCL_A1SR**, STR "from personal communication of Mes-Masson"), plus BIN67 (**CVCL_S987**, STR from originator Garson) and COV434 (**CVCL_2010**).
- **Absent from Cellosaurus (no accession, no STR, no flag): 12** — all 11 BC Cancer/OVCARE **VOA** lines (VOA10816, VOA12539, VOA14993, VOA295, VOA4841, VOA6861, VOA4395, VOA8762, VOA8771, VOA5217, VOA5436) and **TOV3121D** (its patient-matched sibling TOV3121EP = CVCL_A1SN *is* catalogued).

**Problematic/misidentification flags:** only **COV434** carries a Cellosaurus flag — *"Misclassified. Originally thought to be an ovarian granulosa cell tumor but seems to be a small cell carcinoma of the ovary, hypercalcemic type (SCCOHT) (PubMed=33328126)."* This is external, third-party corroboration of the report's COV434→SCCOHT reclassification. No misidentification/contamination flag exists for any other line — notably **TOV112D is still listed as endometrioid** (its dedifferentiated/SMARCA4-null reclassification is *not* yet reflected in Cellosaurus), so that call rests on our multi-omic SWI/SNF evidence + Karnezis 2021, not on an external flag.

**Direct relevance to the report's identity-doubt lines:** VOA8762, VOA8771 (GI-leaning "mucinous"), VOA5436 (clear-cell-like "MMMT"), and VOA4841 (atypical SMARCA4) have **no Cellosaurus entry and no public STR profile at all** — there is no external record to fall back on, which strengthens the report's recommendation to obtain STR + IHC from BC Cancer/OVCARE. SCCOHT lines BIN67/COV434 both have STR + DepMap SMARCA4 support; the reclassification burden is met externally for COV434 and partly for BIN67.

### Draft STR + mycoplasma statement (for Methods / authentication)
> **Cell-line authentication.** No short-tandem-repeat (STR) profiling or mycoplasma testing was performed in-house on the culture stocks used to generate the data in this resource. For 30 of the 42 models a documented STR profile is available in Cellosaurus, deposited by the originating laboratories (Mes-Masson/CHUM for the OV/TOV lines; K. Garson for BIN67), and each model's Cellosaurus accession is provided in Table S[X] (`cellosaurus_str_status.csv`). The 11 BC Cancer/OVCARE VOA lines and TOV3121D are not catalogued in Cellosaurus and have no publicly available STR profile. For the five models also present in CCLE/DepMap (OV90, TOV21G, TOV112D, BIN67, COV434), identity is further supported by high, line-specific concordance between this resource's RNA-seq and DepMap 24Q4 expression (Spearman 0.74–0.88; each line the top match among 67 DepMap ovarian lines) and by concordant driver-mutation calls (Fig. [F], `external_ccle_concordance.csv`). **We recommend that users requiring authenticated stocks obtain STR profiling and mycoplasma clearance, particularly for the VOA lines and for the lines flagged here for identity review (VOA8762, VOA8771, VOA5436, VOA4841).** Cellosaurus records one relevant flag: COV434 is annotated as reclassified to SCCOHT (Cellosaurus CVCL_2010; Karnezis et al. 2021), consistent with this study's reassignment.

*(If any in-house STR/mycoplasma data can still be recovered from the source labs before submission, replace the first sentence accordingly. As stated it is accurate to STR_status=NA in the sample sheet.)*

---

## Where to put this in the report

1. **F4 (Genomics & cell-line authentication)** — add a subsection **"External identity validation (CCLE/DepMap)"**: the overlap table (1a), the concordance result (1b) with Fig. F, and the driver cross-check (1c). Fold the SMARCA4/SCCOHT DepMap corroboration into the existing SWI/SNF paragraph (COV434 + BIN67 now have *genomic* support, not only protein). Add the Cellosaurus flag for COV434 to the reclassification bullet.
2. **F4 / Methods** — add the **STR + mycoplasma authentication** paragraph (draft above) and cite `cellosaurus_str_status.csv` as a supplementary table; update the "Outstanding items #3 (STR/IHC request)" to note the VOA lines are absent from Cellosaurus entirely.
3. **F2 (Technical validation)** — one sentence cross-referencing the CCLE concordance as external, cross-platform technical validation of the RNA pipeline.
4. **F3 or F1** — surface the ConsensusOV table **with the TME caveat** as an orthogonal legacy annotation (not a validated axis).
5. **Reproducibility/Data availability** — cite DepMap Public 24Q4 (Figshare article 27993248) and Cellosaurus (api.cellosaurus.org); accessions/DepMap IDs listed in the output CSVs.

## Confidence & caveats
- **Identifiers (High):** every DepMap ModelID and Cellosaurus accession verified against source metadata; 5/5 RRID↔DepMap cross-checks TRUE; 30/30 Cellosaurus matches are exact-name.
- **Expression concordance (High):** robust to Spearman/Pearson and to all-genes vs HVG; specificity holds against 67 ovarian lines both directions.
- **Driver cross-check (Medium-High):** canonical drivers concordant; two non-canonical variants in our tumor-only WES not corroborated (flagged, not identity-relevant).
- **ConsensusOV (report as external annotation only):** TME-driven subtypes on pure cultures — validity caveat is the point, not a limitation of our extraction.
- **STR (documentation, not in-house authentication):** "documented STR profile" = originator-deposited in Cellosaurus; the resource itself performed no STR/mycoplasma testing.

## Data provenance (for exact reproducibility)
- **DepMap Public 24Q4**, Figshare+ article **27993248**. Files (ndownloader IDs): `Model.csv` (51065297), `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (51065489), `OmicsSomaticMutationsMatrixHotspot.csv` (51065750), `OmicsSomaticMutationsMatrixDamaging.csv` (51065747). Cached in `output/external/`; overlap/ovarian subsets pre-extracted (`depmap_expr_overlap5.csv`, `depmap_expr_ovarian.csv`, `depmap_*_overlap5.csv`).
- **Cellosaurus** via `https://api.cellosaurus.org` (search + cell-line endpoints); 42 JSON responses cached in `output/external/cellosaurus/` (parsed offline — no network needed to re-run parsing).
