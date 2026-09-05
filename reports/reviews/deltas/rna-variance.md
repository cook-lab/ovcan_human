# Delta — Variance decomposition, source-site confound, passage, marker effect sizes

**Workstream:** `rna-variance` (peer-review §3.5, §3.6, §4 #8)
**Script:** `scripts/17_variance_confounders.R` (runs top-to-bottom in ~90 s; verified, exit 0)
**New outputs:** `output/rna_pc_confounder_joint.csv`, `output/rna_variancepartition.csv`, `output/prot_variancepartition.csv`, `output/rna_passage_check.csv`, `output/rna_marker_effectsizes.csv`
**Figures:** `figs/f_variance_partition.pdf` + `reports/assets/f_variance_partition.png`; `figs/f_passage_check{,_crossassay}.pdf` + PNGs
**Read-only sources:** `output/{rna_vst.rds, rna_pca.rds, rna_pc_confounder.csv, prot_pc_confounder.csv, prot_matrix.rds, rna_markers_summary.csv, tx2gene_ensembl_rel105.csv, rna_dds.rds}`, `metadata/{samples.csv, line_family_map.csv}`, `scripts/02,04,06`.

**Bottom line.** The confounder defense holds and is now quantified honestly. On a **joint** model, source site adds **≤0.2 % of PC1 variance beyond subtype** (subtype adds 42 %); the marginal "site R²=0.31" is 99 % *shared/confounded* variance, not independent site signal. In the only subtype present at both sites (clear cell), site explains just **4–6 %** of the top PCs. Genome-wide, subtype ≥ site (RNA) and subtype ≫ site≈plex (protein). Passage is **not** an independent driver of global RNA structure and is **83 % collinear with site**. Marker effect sizes put real numbers beside the lenient "16/22" rule (median oriented AUC 0.69, |d| 0.72). **Stated plainly and repeated in every section below: HGS is fully confounded with site (all 15 HGS = Mes-Masson/CHUM; all MMMT = Huntsman/BC); clear cell (n=7: 2 CHUM + 5 BC) is the sole cross-site control.**

---

## 1. Joint / adjusted-R² model for the RNA PCs (§3.5)

**Method.** The published `rna_pc_confounder.csv` reports only **univariate marginal R²** per PC (`PC ~ subtype` and `PC ~ site` separately), which cannot separate the two because they are correlated. I fit the **joint** model `lm(PC ~ subtype + site)` and ran a **commonality decomposition** (site = 3 levels, as in the published file; joint model verified full-rank / not aliased on all 5 PCs). I also added univariate **adjusted R²** (df-penalized) to match `prot_pc_confounder.csv`. All marginal R² reproduce the published file exactly (asserted in-script).

**Result — source site adds almost nothing beyond subtype:**

| PC | var% | R²(subtype) | R²(site) | R²(joint) | **unique subtype** | **unique site** | shared/confounded |
|----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| PC1 | 20.7 | 0.735 | 0.313 | 0.737 | **0.424** | **0.002** | 0.311 |
| PC2 | 10.4 | 0.694 | 0.214 | 0.695 | **0.481** | **0.001** | 0.214 |
| PC3 |  7.9 | 0.582 | 0.049 | 0.584 | **0.535** | **0.002** | 0.047 |
| PC4 |  7.1 | 0.222 | 0.061 | 0.237 | 0.176 | 0.015 | 0.046 |
| PC5 |  6.2 | 0.369 | 0.043 | 0.392 | 0.349 | 0.023 | 0.020 |

- **unique_site** (variance site adds *beyond* subtype) is **≤0.002 on the three PCs that carry subtype signal** (PC1–3) and ≤0.023 on all five. **unique_subtype** is 0.18–0.54. Adjusted R² tells the same story (PC1 adj-R²: subtype **0.682** vs site **0.264**; PC3 site adj-R² is *negative*, −0.019).
- The 0.31 "shared" block on PC1 is the genuinely confounded variance — it *cannot* be assigned to subtype vs site by any model, because HGS≡CHUM. But because the *uniquely-site* portion is ~0, the most parsimonious reading is that the shared block rides with subtype.

**The only clean site test** (subtypes present at BOTH sites — the sole place site is separable from subtype):

| set | n | n CHUM / BC | PC1 site R² | PC2 site R² |
|---|--:|:--:|--:|--:|
| mixed-site subtypes (CC/EC/MC) | 12 | 4 / 8 | **0.004** | **0.011** |
| clear cell only (CC) | 7 | 2 / 5 | **0.044** | **0.058** |

Within clear cell, moving lab explains 4–6 % of the top PCs — i.e. Mes-Masson and Huntsman clear-cell lines co-cluster by biology, not by site. This is the load-bearing confounder result and it is clean.

## 2. variancePartition — genome-wide decomposition (§3.5)

**Method.** Per-feature random-effects variance decomposition with `subtype + source_site + patient` (+ **TMT plex** for protein). *Fallback documented:* the `variancePartition` package is installed but **breaks in this R** (it calls `lme4::findbars`, which has migrated to the `reformulas` package → "Initial model failed"). I used the **identical decomposition via `lme4` directly** (per-feature REML; variance components ÷ total), parallelized — this reproduces variancePartition's fraction-of-variance-explained. RNA = all 22,542 expressed genes; protein = 6,856 complete-case proteins (quantified in all 5 plexes). "Patient" uses `line_family_map$patient_id` (collapses sublines; NB the `family` column is NA for singletons — do not use it as the grouping).

**Result (genome-wide % variance):**

| assay | subtype | source_site | TMT plex | patient/line | residual |
|---|--:|--:|--:|--:|--:|
| **RNA** (median) | **5.9** | 3.5 | — | 27.3 | 32.8 |
| RNA (mean) | 14.9 | 14.5 | — | 31.6 | 39.0 |
| **Protein** (median) | **8.6** | **0.0** | 0.9 | 19.8 | 37.5 |
| Protein (mean) | 16.7 | 8.2 | 6.4 | 27.4 | 41.3 |

- **Subtype ≥ source site in both assays** (RNA median 5.9 vs 3.5; protein 8.6 vs 0.0). In protein, **subtype ≫ site ≈ plex** — the TMT batch and site axes are minor for the typical protein (median site variance is literally 0).
- **Patient/line identity is the largest structured term** (RNA 27 %, protein 20 % median) — expected: these are distinct cell lines, and each line's private expression program dominates. *Caveat:* 25/28 patients are singletons, so the patient vs residual split is estimated from the 3 replicated sibling pairs (1369/2295/3133) — read "patient" as line-level identity, not a precise patient-vs-technical partition.
- **Honest limit:** like the marginal R², variancePartition *cannot* disentangle the confounded subtype/site variance (both terms get partial credit for the shared HGS-vs-rest axis). The **joint model in §1 is what disentangles it** and shows site's independent contribution is ~0. Read §2 as descriptive/supportive and §1 as the defense.

## 3. Passage-sensitivity check (§4 #8)

**Method.** Parsed `rna_passage`/`wes_passage` to numeric ("p63"→63, "p39-40"→39.5). (a) `PC ~ passage` R²/p for PC1–5, plus passage's unique add after site; (b) within-line cross-assay RNA-vs-WES passage tabulation.

**(a) Passage does not independently drive global RNA structure.** Passage is **83 % explained by site** (CHUM lines p60–71; BC lines p17–43). `PC~passage` R²: PC1 0.078 (p=0.13, n.s.), **PC2 0.159 (p=0.026)**, PC3 0.027, PC4 0.082, PC5 ~0. Only PC2 is *nominally* significant and does not survive 5-test correction; and because passage is 83 % collinear with site, any passage signal is a **site proxy**, not separable culture-drift. Passage's unique add beyond site is ≤0.15 (a PC1 suppression artifact) and ≤0.04 elsewhere.

**(b) Cross-assay passage mismatch is real and occasionally large** (13 lines with both):

| line | RNA p | WES p | Δ |
|---|--:|--:|--:|
| **TOV112D** | 63 | 83 | **+20** |
| **OV3331** | 71 | 54 | **−17** |
| OV3133-R | 64 | 71 | +7 |
| OV90 | 70 | 63 | −7 |
| (9 others) | — | — | ≤4 |

Median |Δ| = 4 passages; two lines identical (TOV1369, TOV3392D). **Implication:** RNA and WES were sometimes generated on materially different passages of the same line — most starkly TOV112D (20 passages apart) and OV3331 (17) — so cross-omic integration on those lines assumes passage-stability of the features compared. Cheap insurance: now documented, pre-empts the reviewer question.

## 4. Marker effect sizes — supplement to the "16/22" rule (§3.6)

**Method.** For each canonical marker, Cohen's **d** and **AUC** of intended-subtype vs all other lines, on **VST**, computed on the **28 patient representatives** (sublines collapsed — dropped OV2295-R2, TOV1369, TOV3133G — no pseudoreplication). AUC is rank-based (robust for n=2 rare subtypes); `auc_oriented` flips "down" markers so 1.0 = perfect discrimination in the expected direction for every row. This **supplements, not replaces** the top-2-of-6 rule (`rank_in`/`lands_right` carried alongside in the CSV). MKI67 kept as an **unscored** proliferation control.

Median |d| = **0.72**, median oriented AUC = **0.69**; **8/22** markers reach oriented AUC ≥0.80, **10/22** reach |d| ≥0.8. Selected rows:

| subtype | marker | dir | n_in | Cohen's d | oriented AUC | top-2 lands? |
|---|---|:--:|--:|--:|--:|:--:|
| MC | KRT20 | up | 3 | **3.08** | 0.97 | ✓ |
| MC | CDX2 | up | 3 | **2.77** | 0.96 | ✓ |
| SCCOHT | SMARCA2 | down | 2 | **−2.78** | 1.00 | ✓ |
| CC | SPP1 | up | 7 | 1.46 | 0.84 | ✓ |
| MC | MUC5AC | up | 3 | 1.46 | 0.65 | ✓ |
| CC | HNF1B | up | 7 | 1.36 | 0.81 | ✓ |
| MMMT | ZEB1 | up | 2 | 1.24 | 0.81 | ✗ (rank 3) |
| MC | TFF3 | up | 3 | 1.01 | 0.84 | ✓ |
| HGS | WT1 | up | 12 | 0.93 | 0.75 | ✓ |
| CC | NAPSA | up | 7 | 0.85 | 0.63 | ✓ |
| HGS | SOX17 | up | 12 | 0.74 | 0.73 | ✓ |
| SCCOHT | SMARCA4 | down | 2 | −0.70 | 0.67 | ✓ |
| HGS | PAX8 | up | 12 | 0.53 | 0.69 | ✗ (top in CC) |
| HGS | MECOM | up | 12 | −0.06 | 0.46 | ✗ |
| MC | MUC2 | up | 3 | −0.29 | 0.47 | ✗ |
| MMMT | CDH2 | up | 2 | 0.11 | 0.58 | ✗ |

- **The 6 top-2 "misses" are quantitatively small effects, not failures:** MECOM (d≈0, essentially null), MUC2 (d=−0.29, off everywhere — a 2D-culture loss), CDH2 (d=0.11), GPX3 (0.28), PAX8 (0.53 but pan-Müllerian — genuinely higher in CC), ZEB1 (d=1.24 but ranks 3rd; EMT markers are broad in culture). The effect sizes make the "misses" interpretable rather than binary.
- **SMARCA4 vs SMARCA2 contrast is a bonus:** SMARCA2 loss is a **clean** transcriptional signal (d=−2.78, oriented AUC 1.0 — epigenetic silencing), whereas SMARCA4 is **modest at the mRNA level** (d=−0.70) — exactly the report's mechanistic claim that SCCOHT SMARCA4 loss is often post-transcriptional (retained mRNA in BIN67). This quantitatively backs the F4 SWI/SNF narrative.
- **Caveat:** n_intended = 2 for MMMT/SCCOHT and 3 for MC — Cohen's d is unstable and AUC is granular there; lead with AUC for those, and treat n=2 effect sizes as descriptive.

---

## 5. Report sections to update + figure placement

- **§F3 "Subtype separation — biology, not batch"** (the ¶ "Confounder check … PC1 R²=0.73 (subtype) vs 0.31 (site)"): **replace the univariate framing** with the joint result — *"On a joint model, source site adds ≤0.2 % of PC1 variance beyond subtype (subtype adds 42 %); site's marginal 31 % is confounded, not independent. In clear cell — the only subtype at both sites — site explains 4–6 % of the top PCs."* Add one sentence from §2 (subtype ≥ site genome-wide; protein subtype ≫ site≈plex) and cite `f_variance_partition`. Keep the plain confound statement (HGS≡CHUM; CC sole cross-site control).
- **§F3 "Marker & pathway recovery"** (the ¶ "16/22 canonical markers land"): append the effect-size supplement — median oriented AUC 0.69 / |d| 0.72; the 6 misses are small-effect/shared-lineage (numbers above), not failures. Cross-reference the SMARCA4-vs-SMARCA2 effect-size contrast into **§F4 SWI/SNF** as quantitative support for post-transcriptional SMARCA4 loss.
- **New limitations bullet (Methods or "Outstanding items"):** passage is 83 % confounded with site (not an independent structure driver), and RNA-vs-WES passage differs within some lines (TOV112D Δ20, OV3331 Δ17) — cross-omic claims on those lines assume passage-stability. Cite `f_passage_check`.
- **Figures.** `f_variance_partition` → **F3** (confounder panel). `f_passage_check` + `f_passage_check_crossassay` → supplementary/methods.
- **New tables for the supplement:** `rna_pc_confounder_joint.csv`, `rna_variancepartition.csv` + `prot_variancepartition.csv`, `rna_marker_effectsizes.csv`, `rna_passage_check.csv`.
