# Data-Science Fidelity Audit — Report ↔ Code ↔ Output

**Report audited:** `reports/01_multiomic_characterization_results.md`
**Code:** `scripts/00_setup.R` … `scripts/14_hgs_heterogeneity.R`
**Ground truth:** `output/*.csv`, `metadata/samples.csv`, and (for two headline "before" numbers) the raw archive under `judy_archive/`.
**Auditor:** independent code-to-report reconciliation. Read-only on all project files.
**Date:** 2026-07-23

---

## 1. Overall fidelity assessment

**The report is an unusually faithful representation of the analysis the code actually performed.** I extracted ~60 checkable quantitative/methodological claims and traced each to a materialized output, a recomputation from the primary matrices, or the executing script. **Every substantive number matches** — cohort/coverage counts, QC medians, PCA/confounder R², silhouettes, marker recovery (16/22), RNA–protein concordance, the full WES driver landscape (including the pre-filter→post-filter collapse), CNV arm-level frequencies, authentication tallies, ADC shortlists, and the HGS-heterogeneity strata. Critically, the assay-aware framing is real, not cosmetic: the Mutect2 `FILTER==PASS` gate **is** applied (`07:151–155`), the tumor-only caveat is honored (drivers reported, no burden metric), genomic HRD is **not** fabricated (`09`), the archived CC/MC signature bug is documented and asserted-fixed (`03:25–31,119–124`), seeds are set before every stochastic step, and Ensembl 105 is pinned/cached. I found **no Critical or Major discrepancy** — no wrong number, inverted filter, swapped group level, stale value, or method described-but-not-implemented. The residual issues are Minor: one headline metric-labeling imprecision in the Summary, one loose gene-count juxtaposition, one slightly over-general "families co-cluster" phrase, and routine reproducibility gaps the report already flags (renv pending). **Confidence: HIGH.**

---

## 2. Discrepancy ledger (worst first)

| # | Report claim (location) | Code/output says | Match? | Severity |
|---|---|---|---|---|
| 1 | **Summary (l.13):** "RNA–protein concordance matches the CPTAC-ovarian benchmark (**median Spearman 0.41** vs 0.38–0.45)" | 0.41 is the **per-line** median (`integ_rnaprot_cor_summary.csv`: per_line_spearman median 0.4076). Script 12 explicitly states per-line is **NOT** benchmark-comparable and that the gene-wise median (0.3971≈**0.40**) is "the headline number" (`12:16–24`). The benchmark-comparable value (0.40) is also in range, so the conclusion holds; the Summary cites the less-appropriate metric. F3 (l.49) does report both. | Conclusion ✓; metric label imprecise | **Minor** |
| 2 | **F3 (l.49):** "30 lines / **8,212 shared genes (97.4% mapping)** … per-gene median 0.40" | 8,212 shared genes / 97.4% is correct (recomputed: 8212 / 8429 = 97.43%). But the per-gene median (0.40) is computed on the **7,893** genes with protein in ≥10 lines (`12:145–163`; summary file per_gene n=7893), not 8,212. Juxtaposition could mislead; neither number is wrong. | Both numbers ✓; adjacency loose | **Minor** |
| 3 | **F5 (l.84):** "same-patient families **co-cluster**" | 2295 (OV2295, OV2295-R2) and 3133 (OV3133-R, TOV3133G) co-cluster, but the **1369 family splits** (OV1369-R2 → cluster 1 Inflammatory; TOV1369 → cluster 3 Hypoxic-glycolytic) per `hgs_heterogeneity.csv`. Script 14's own console note is more careful: "families (e.g. 2295, 3133) **mostly** co-cluster" (`14:230`). Low stakes (explicitly descriptive, n=15). | Over-general | **Minor** |
| 4 | **F2 (l.30) & Repro:** "Ensembl **105**" | Executing code uses `useEnsembl(version = 105)` (`01:46`) and report matches. Only a **stale header comment** says "release 104" (`01:34`); inline comment (`01:36–38`) explains 104 is retired → 105 used. Cosmetic, code-only. | Report ✓ (comment stale) | **Minor (cosmetic)** |

No rows rated Critical or Major.

---

## 3. Code correctness & reproducibility concerns

**Correctness — clean.** Spot-checks of the high-risk areas all passed:
- **Mutect2 FILTER honored:** `raw[, is_pass := FILTER == "PASS"]` and `retained := is_pass & !is_common_germ` (`07:151–155`); the report's whole "root-cause fix" narrative is truthfully implemented. Pre-filter vs retained both verified against raw MAFs (§4 numbers).
- **DESeq2 counts path:** `DESeqDataSetFromTximport` with length offset; design `~ subtype` (`01:63`), one-vs-rest refits `~ grp` reusing the tximport normalizationFactors (`03:70–82`) — correct, not a re-created counts matrix.
- **PCA/silhouette on the right matrix:** PCA on top-2000-var VST genes (`02:26–27`); silhouette on Euclidean PC1–10 (`02:53–54`); confounder R² = `lm(PC ~ factor)$r.squared` for subtype and site (`02:62–66`) — exactly "variance per PC explained by subtype vs site."
- **Concordance ordering safe:** `R` and `P` matrices are intersected on shared genes/lines and `stopifnot(identical(dimnames(R), dimnames(P)))` (`12:123–125`); Spearman with pairwise-complete — no mismatched-order correlation.
- **CNV FGA:** `fga_0.2` = fraction of assessed genome with |median-centred log2| > 0.2 (`08:113–121`); report's "FGA 0.63" for HGS = mean/median of `fga_0.2` (both ≈0.627). Arm-event recovery recomputed exactly (§4).
- **No pseudoreplication in headline stats:** same-patient sublines (2295/3133/1369/3121/3291 families) are retained as distinct lines, but the report never treats them as independent replicates for a significance claim; per-family concordance is noted, not tested.

**Reproducibility — solid, with minor gaps (most already disclosed):**
- **Seeds:** `SEED <- 1234; set.seed(SEED)` in `00_setup.R:63–64`, re-set in every stochastic script — verified `set.seed` before `Rtsne` (`02:43–44`, `06`), and in `14` (progeny/hclust). Good.
- **Annotation pinned/cached:** `tx2gene_ensembl_rel105.csv` cached in `output/`; live biomaRt is bypassed when the cache exists (`01:39–53`). Good.
- **`library(biomaRt)` loaded unconditionally** at `01:15` even though the cached CSV means it is never queried, and biomaRt is **not** in the `.required` list (`00_setup.R:30–38`). Latent failure: on a machine without biomaRt, `01` errors at load despite not needing it. **Fix:** gate the load behind the cache-miss branch.
- **Machine-specific default path:** `PROJ <- Sys.getenv("OVCAN_PROJ", unset = "/Users/dpcook/Analysis/ovcan_human")` (`00_setup.R:17`) — overridable via env var, but the default is hardcoded to one machine. Acceptable; worth parameterizing for deposition.
- **Run order:** numbered scripts have a clean dependency DAG through `output/` (`02`←`01`; `03`,`04`←`01`; `06`←`05`; `12`←`01`+`05`; `14`←`04`); `01→14` in order will run top-to-bottom. `renv` lockfile is **not** present — the report states this openly ("`renv` lockfile still to be added").
- **Package versions** are checked for presence (`check_pkgs`) but not pinned; env is R 4.5.2 / BiocManager 3.21 (consistent with `wes_hrd_feasibility.md`).

---

## 4. Claims verified accurate (traced to output/recompute)

**Cohort & coverage** (`metadata/samples.csv`, recomputed in R): 42 generated + 13 external = 55; RNA 31, proteomics 31, WES 23 (CNVkit); 13 tri-omic, **all Mes-Masson** ✓. WES MAFs: 23 lines called but **22 usable** (TOV3121D flagged "no MAF"), matching report's "23 CNVkit; 22 Mutect2 MAFs" and asserted in `07:72,82`. RNA subtypes 15/7/2/3/2/2 (HGS/CC/EC/MC/MMMT/SCCOHT) ✓.

**RNA QC** (`rna_qc_metrics.csv`): 39,568 total genes (`rna_tpm.csv` rows) → 22,544 expressed ≥10 in ≥2 (`rna_counts.csv` rows) ✓; pseudoalign median 91.1%, range 85.8–93.1 ✓; lib-size median 56.5M ✓; detected median 20,119 ✓; site medians Mes-Masson 92.2% vs Huntsman 88.1% ✓.

**Separation** (`rna_pc_confounder.csv`, `rna_silhouette.csv`, `prot_*`): PC1 20.7%/PC2 10.4% ✓; confounder R² PC1 0.73/0.31, PC2 0.69/0.21, PC3 0.58/0.05 ✓; silhouettes HGS 0.16/CC 0.12/MC 0.15/MMMT 0.74/SCCOHT 0.82/EC −0.01 ✓; proteomic HGS 0.18, CC −0.003 (not cohesive) ✓; proteomic PC1 subtype 0.55 > site 0.36 > plex 0.13 (real-but-secondary) ✓.

**Proteomics** (`prot_qc.csv`, raw xlsx): 8,430 proteins ✓ and 146,830 peptides ✓ (both recomputed from the source xlsx); presence≥50% → 7,734, complete/na.omit 6,856, +878 = **+12.8%** ✓; CV median 5.31% ✓; bridge Pearson 0.991–0.994 (≈0.99) ✓ (`prot_bridge_cor.csv`).

**Markers** (`rna_markers_summary.csv`, `04:153–156`): 16/22 land (criterion = up-marker in top-2 of 6 subtype means AND mean log2TPM>1; down-marker in bottom-2) ✓; the 6 misses (PAX8/MECOM/GPX3/MUC2/ZEB1/CDH2) are coherent — PAX8 rank-3 because shared HGS+CC Müllerian, MUC2 rank-2 but mean 0.02<1 (off) ✓. Hallmark = MSigDB v7.4 singscore (`04:220–224`).

**GO recovery** (`rna_de_gsea_recovery.csv`): HGS DNA-repair recovered, CC glutathione recovered, SCCOHT cell-cycle recovered, MC glycan suggestive, MMMT EMT not recovered ✓; top HGS GO axis = mitochondrial/cytoplasmic translation + OXPHOS (`rna_de_gsea_go.csv`) — report's "OXPHOS/ribosome-biogenesis … confounded with site" is fair ✓.

**Concordance** (`integ_rnaprot_cor_summary.csv`, recomputed): 30 shared lines, 8,212 shared genes, 97.4% mapping ✓; per-line median 0.41 (IQR 0.36–0.44), per-gene median 0.40 ✓; all six benchmark values (CPTAC 0.45/0.38, CCLE 0.48, ProCan 0.42, Jarnuczak 0.58, ceiling 0.72) transcribed correctly from `12:66–74`.

**WES SNV** (`wes_mutations_filtered.csv`, `wes_driver_freq_by_subtype.csv`, raw MAFs): **TP53 17/17 HGSC (100%)**, 0 HGS without ✓; per-line retained coding median **206.5** (~206) ✓; artefact collapse **verified end-to-end** — pre-filter ATM 18/22=82%, ATR 17/22=77%, BRCA2 15/22=68%, POLE 6/22=27% (raw MAFs) → retained 9.1%/4.5%/9.1%/0% ✓; OV2295 **25,914 rows → 493 PASS** ✓; TP53 hotspots consistent within families (2295→I195T, 3133→Q192*, 1369→G244C) ✓; KRAS in all CC/EC/MC ✓; TOV112D SMARCA4 **p.L639X** + TP53 **R175H** ✓; OV90 **SMAD4 p.R445\*** nonsense + BRAF ✓; TOV2414 KRAS G12A + SMAD4 frameshift ✓.

**CNV** (`wes_cnv_fga.csv`, `wes_cnv_segments.csv` recomputed): arm-event recovery 3q26 89%, MYC/8q24 78%, CCNE1/19q12 72%, 20q 72%, TP53/17p-loss 78%, RB1/13q 56%, PTEN/10q 28% — **all exact** (n=18 HGS) ✓; FGA HGS 0.63 > CC 0.40 > MC 0.36 > EC 0.27 > LGS 0.07, TOV81D quietest ✓; OV90 lowest HGS FGA 0.308 ✓.

**Authentication** (`auth_perline_table.csv`, `auth_swisnf_panel.csv`): 26 expression-consistent, 20 genomics-consistent ✓; BIN67 "instructive case" — SMARCA4 mRNA retained (rank 23/31) but protein 2nd-lowest (rank 2), SMARCA2 mRNA lowest (rank 1) ✓; COV434 SMARCA4 low protein + SMARCA2 mRNA low ✓; VOA4841 lowest SMARCA4 both layers (rank 1 RNA & protein) ✓; VOA5436 CC program (HNF1B/NAPSA z≈2.0, discordant) ✓; VOA4395 sole EC ✓; OV90 lowest serous (z=−1.25) → "HGS-family, serous not confirmed" ✓; OV3331 clean HGS ✓.

**Mucinous** (`auth_mucinous.csv`): TOV2414 KRT7+/PAX8+/MUC5AC+, KRAS G12A, SMAD4 loss, externally authenticated ✓; VOA8762/VOA8771 CDX2-high/CK7-low/PAX8-low, GI-leaning, VOA8771 most GI-leaning (lowest ovarian_index −1.91) ✓; genuine ovarian mucinous n=1 ✓.

**HRD** (`wes_hrd_feasibility.md`, `09`): NOT FEASIBLE; 0/23 CNVkit runs used `--vcf`; 0 tumor BAMs archived (only public normal SRR4039087.bam); no score fabricated; Peng expression signature dropped ✓.

**ADC** (`adc_expression.csv`, `adc_subtype_summary.csv`): mesothelin→HGS ✓; HER2/ERBB2→CC/MC top TOV3392D (cross-assay) ✓; NaPi2b/SLC34A2→VOA12539 (cross-assay) ✓; CDH6→OV3331 (cross-assay) ✓; FOLR1 bimodal within HGSC (9.48→0.06), top line TOV3133G ✓; DPEP3 RNA-only (no protein row) ✓.

**HGS heterogeneity** (`hgs_heterogeneity.csv`): 3 strata with exact membership — Inflammatory/NF-κB-EMT (OV90, OV1369-R2, OV866-2, OV1946), Low-signaling (OV2085, OV2295, OV2295-R2, OV3133-R, TOV3133G, OV4485), Hypoxic-glycolytic (OV3331, OV4453, TOV1369, TOV3041G, TOV3291G); 4/6/5 = 15 ✓; Hallmark-singscore ward.D2 k=3 + PROGENy orthogonal corroboration (`14`) ✓.

**Method framing:** formal DE restricted to HGS (n≥10; `03:60`), rare subtypes descriptive ✓; CC/MC signature bug documented + `stopifnot` guards (`03:119–124`) ✓; LGSC/Carey excluded as external ✓.

**Figures:** all 35 `figs/*.pdf` exist; all 32 report-referenced panels are present — "every panel referenced here exists in `figs/`" holds ✓.

---

## 5. Could not verify (and why)

- **"kallisto 0.46"** (F2, l.30): the aligner version is provenance metadata in the archive's `run_info.json`, not carried in any `output/*.csv`. The scripts import pre-computed `abundance.tsv`; they do not run kallisto. The pseudoalignment percentages that depend on that run **are** in `metadata/samples.csv` and were verified.
- **Sequencing run "NS.1676.003" / "5-plex TMT"** (F1, l.21–22): the run ID is provenance not in outputs; the 5-plex TMT structure **is** consistent with the data (plex levels 1–5 in `prot_sample_qc.csv`, `05`).
- **Published-benchmark literature values** (CPTAC/CCLE/ProCan/Jarnuczak): verified to match the values hardcoded in `12:66–74`, i.e., report ↔ code is consistent; I did **not** re-verify those against the primary publications (out of code/output scope).
- **`.rds`-only internals** (e.g., `rna_pca.rds`, `prot_matrix.rds`): not opened directly; every derived quantity was checked via the corresponding CSV or recomputed from `rna_tpm.csv` / `prot_abundance_matrix.csv` / `wes_*` CSVs.
- **06 proteomic-marker panel figure** content not inspected; the numeric proteomic-separation claims (PC confounder, silhouette) were verified from `prot_pc_confounder.csv` / `prot_silhouette.csv`.
