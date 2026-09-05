# Peer Review — *Multi-Omic Characterization of the OvCAN Ovarian Cancer Cell-Line Panel*

**Venue:** *Scientific Data* (Nature Portfolio), Data Descriptor
**Reviewer role:** Cancer genomics / proteomics / ovarian-cancer biology / multi-omic analysis
**Material reviewed:** `reports/01_multiomic_characterization_results.md`; figures in `reports/assets/`; result tables in `output/`; `metadata/samples.csv`. All numbers below were recomputed or spot-checked against the underlying tables (see §6).

---

## 1. Summary assessment

This is a well-organized, unusually self-aware characterization of a genuinely valuable resource: uniform in-house bulk RNA-seq, TMT proteomics, and whole-exome sequencing across 42 ovarian-cancer cell-line models spanning seven histotypes, with 13 lines profiled on all three assays. For a *Scientific Data* Data Descriptor — where the bar is **technical soundness, reuse value, and honest metadata**, not biological novelty — the core of this work is close to what the venue wants. The cohort accounting is exact and reproducible (I reconstructed 42 generated / RNA 31 / proteomics 31 / WES 23 (22 with MAFs) / all-three 13 directly from `samples.csv`), the QC is quantified, and the technical-validation logic (subtype separation, marker recovery, RNA–protein concordance benchmarked against CPTAC/ProCan/CCLE, canonical WES drivers) is the right backbone. The report is also refreshingly candid about its own limits — tumor-only WES, non-computable genomic HRD, the OXPHOS/site confound, the chrX artifact, structural proteomic missingness — and it explicitly avoids fabricating an HRD score. That candor is the manuscript's single greatest asset and should be preserved.

That said, **it is not yet publishable as-is**, and the gap is mostly interpretive/statistical rather than a failure of the underlying data. The most serious problem is **pseudoreplication**: the panel contains several sets of sublines/isolates from the same patient (families 1369, 2295, 3133, 3291, 3121), and I confirmed at the variant level that these share *identical* mutations (all four 3133 lines carry TP53 p.Q192\* and CDK12 p.L123X; all three 2295 lines carry TP53 p.I195T; the 1369 pair shares p.G244C). Every genomic frequency in the report is computed over lines, not independent patients, which inflates denominators (HGSC n=17 lines → 11 patients) and manufactures apparent recurrence (CDK12 "35% of HGSC" is essentially one patient's frameshift counted four times). The oncoplot (`f_wes_oncoplot.png`) displays these four identical columns adjacently, so the redundancy is visually undeniable and a reviewer will seize on it. Second, the paper repeatedly uses the word **"authenticate,"** but no STR profiling or mycoplasma testing was performed — and *Scientific Data* requires cell-line authentication and mycoplasma statements. What was actually done is molecular *corroboration and discordance-flagging*, which is valuable but must be named accurately, especially since the analysis itself raises new identity doubts (VOA8762/VOA8771). Third, several **tumor-only WES driver calls are not defensible as somatic** (a BRCA2 truncation called in the single, near-diploid LGSC line; a BRCA2 A22E missense; CDK12 calls with germline-range VAF), and the `germline_like_vaf` flag does not do the work implied — in ~100%-pure lines it flags bona fide somatic TP53 hotspots (R175H, Q192\*) as "germline-like" because LOH drives VAF→1.

My overall assessment: **major revision.** The evidence level is **Medium-to-High for the descriptive/technical-validation claims** (subtype separation, marker/pathway recovery, concordance benchmarking, the SMARCA4/SMARCA2 reclassification triad), and **Low-to-Medium for the genomic-frequency, "authentication," and rarer-driver claims** until pseudoreplication and tumor-only artifacts are handled. None of the required fixes are fatal; most are re-tabulation, re-wording, one figure edit, and a few added diagnostics.

---

## 2. Major strengths

1. **Correct framing for the venue.** "Technically sound + recapitulates known biology" is exactly the Data-Descriptor standard; the report does not overreach into discovery claims (the HGSC heterogeneity strata are explicitly labelled "model-selection example, no survival/discovery claims"). This is the right instinct.

2. **Exceptional transparency about limitations.** Tumor-only WES is flagged and drivers are (mostly) presented descriptively rather than as a burden metric; genomic HRD is *declared non-computable rather than faked* (`output/wes_hrd_feasibility.md`); the OXPHOS/ribosome GO axis is disclosed as site-confounded; the chrX pooled-normal artifact and structural proteomic missingness are named. This is better than most published descriptors.

3. **The multi-omic SWI/SNF reclassification is genuinely convincing and well-visualized.** `f_auth_swisnf.png` and `auth_swisnf_panel.csv` support TOV112D (SMARCA4 RNA z=−2.14, truncating p.L639X, protein loss; SMARCA2 mRNA-low), COV434, and — instructively — BIN67 (SMARCA4 mRNA *retained*, rank 23/31, but protein 2nd-lowest; SMARCA2 mRNA lowest). The mechanistic distinction drawn (SMARCA4 = mutation/post-transcriptional; SMARCA2 = epigenetic silencing) is biologically correct and the "RNA-only would miss BIN67 → SMARCA4 IHC is confirmatory" point is a real, useful insight.

4. **Concordance is benchmarked honestly.** Per-line Spearman median 0.408 (verified in `integ_rnaprot_cor_summary.csv`) is placed against external benchmarks in `f_concordance.png`, and the figure makes clear OvCAN sits at the *low* end of the cell-line range (below ProCan/CCLE/Jarnuczak) rather than spinning it. Replacing a 2-line correlation with an n=30 distribution is a real improvement.

5. **Small-n honesty in the clustering.** Silhouette widths are reported with n and the n=2 inflation is called out (MMMT 0.74/SCCOHT 0.82 flagged "inflated"; EC ≈0 with the note that the two "EC" lines don't co-cluster). `rna_silhouette.csv` matches the text exactly.

6. **The re-filtering of the WES SNVs is a substantive, correct fix.** Recovering the `FILTER==PASS` field (OV2295: 25,914 → 493) and applying population-AF removal collapses implausible archived rates (ATM 82→9%, ATR 77→5%) to a credible landscape and yields the TP53 positive control. This is the kind of methods-level correction a descriptor should document.

---

## 3. Major concerns (ranked most → least serious)

### 3.1 Pseudoreplication: genomic frequencies are computed per line, not per patient *(most serious)*
`samples.csv` documents five same-patient families, and I confirmed in `wes_mutations_filtered.csv` that family members carry **identical** somatic variants:
- **3133 family** (OV3133-R, OV3133-R2, TOV3133D, TOV3133G): all four share TP53 p.Q192\* *and* CDK12 p.L123X (frame-shift) + p.CP1255-1256C.
- **2295 family** (OV2295, OV2295-R2, TOV2295-R): all three share TP53 p.I195T.
- **1369 pair** (OV1369-R2, TOV1369): both TP53 p.G244C.

Consequences I verified:
- **HGSC denominators are inflated:** 17 HGSC lines with MAFs represent **11 independent patients**. The report's "TP53 17/17 = 100%" survives collapse (11/11), so the positive control is fine — but this is the *only* frequency that is robust.
- **CDK12 "6/17 = 35.3%"** (`wes_driver_freq_by_subtype.csv`; oncoplot renders it 27% over 22 lines) collapses to **3/11 patients ≈ 27%**, and of those the 3133 event is one frameshift counted four times and the OV3291 event has germline-range VAF (`germline_like_vaf`≈0.985). The true independent CDK12 signal is ~1–2 events — versus the ~3% expected in TCGA HGSC. As written, the resource implies a striking CDK12 enrichment that is an artifact of subline counting.
- By contrast, I recomputed **CNV frequencies and they are *robust* to the collapse** (3q26 gain 83% lines → 82% patients; 17p loss 83% → 91%), because HGSC arm-level events are near-universal. This asymmetry should itself be stated: pseudoreplication distorts *mutation* frequencies badly and *CNV* frequencies little.

**Remedy (essential):** report all genomic frequencies at the **patient level** (collapse families to one representative, or count "≥1 line mutated per patient"), state n_patients alongside n_lines everywhere, and add a **patient-family annotation track** to the oncoplot, CNV heatmap, PCA, and concordance scatter so the redundancy is explicit. The same collapse should be at least acknowledged for the RNA PCA/silhouette (n=31 lines ≈ 27 patients) and the n=30 concordance distribution, even if you keep all lines for the atlas.

### 3.2 "Authentication" is claimed but STR / mycoplasma were not done *(serious — also an editorial requirement)*
The Summary and F4 say the resource "independently authenticate[s] cell-line identity." In cell-line practice (ICLAC/ANSI), *authentication* means **STR profiling**; none was performed — indeed STR is listed as an *outstanding request* (Outstanding items #3), and `auth_perline_table.csv` has `STR_status = NA` for all 42 lines. *Scientific Data* (and Nature Portfolio generally) **require** a cell-line authentication statement and a mycoplasma-contamination statement; neither is present. This is compounded by the fact that the analysis actively raises new identity doubts (VOA8762/VOA8771 look intestinal; VOA5436 looks clear-cell), which STR would resolve.

**Remedy (essential):** (i) replace "authenticate" with "molecularly corroborate identity / flag identity discordances" throughout; (ii) add the STR-authentication and mycoplasma statements the journal requires — at minimum cite the originators' STR records (Cellosaurus IDs; e.g. TOV2414 = CVCL_A1SR is already noted) and state which lines have documented STR vs. none; (iii) frame the discordance flags as *hypotheses requiring STR/IHC* (the report mostly does this — keep it consistent with the softened top-line wording).

### 3.3 Tumor-only WES: specific driver calls are not defensibly somatic, and the germline flag is uninformative in pure lines *(serious)*
- **The `germline_like_vaf` flag conflates germline with somatic-LOH.** Of 50 driver calls, 21 are flagged `germline_like_vaf=TRUE`, including canonical somatic TP53 hotspots (TOV112D R175H, the 3133 Q192\*, OV3291 R249W). In ~100%-pure lines every clonal somatic mutation sits at VAF≈0.5 or ≈1.0 (after LOH), so a VAF-threshold "germline" flag mostly flags real somatic events and cannot separate germline from somatic. The report should stop implying this flag adds somatic confidence.
- **Biologically incoherent calls persist.** `auth_perline_table.csv` lists **BRCA2 as the key driver for TOV81D**, the single LGSC line — which is also the *quietest* genome in the panel (FGA 0.073). A functional biallelic BRCA2 loss is incompatible with an HR-proficient, near-diploid genome; this is almost certainly a monoallelic/germline-carrier or passenger variant, and listing it as a "key driver" is misleading. The other BRCA2 call (TOV3133D p.A22E) is an N-terminal missense of no established pathogenicity. So the residual "BRCA2 9%" is not evidence of HR deficiency in any line.
- **BRCA1/2 germline problem specifically:** in HGSC patients germline BRCA1/2 variants are common and clinically meaningful; tumor-only calling with population-AF filtering *cannot* classify a rare pathogenic BRCA1/2 variant as somatic. Any BRCA1/2 statement must carry this caveat.

**Remedy (essential):** (i) exclude `germline_like_vaf`-flagged variants from driver *frequencies* (or present them as a clearly separated "cannot exclude germline" tier); (ii) add per-variant somatic-vs-germline caveats and drop/annotate the TOV81D BRCA2 and TOV3133D BRCA2 A22E calls; (iii) restrict headline driver claims to the events that are robust in tumor-only pure lines — TP53 (universal, LOH-supported), the KRAS hotspots (G12/G13, canonical), and SMARCA4 truncation in TOV112D.

### 3.4 CNV: reference choice and method affect *validity*, not just deposition; chrX artifact is still shown *(serious)*
- **Pooled unmatched public normals.** CNVkit used "5 unmatched PUBLIC exomes (PRJNA339046)" as the reference (figure footnote). If the capture kit differs from the tumor exomes, kit-specific coverage differences masquerade as copy-number changes genome-wide — this is a validity concern for the CNV layer, not merely a "strip before deposit" hygiene item. The report flags kit concordance as "unconfirmed."
- **Per-sample median-centering in high-FGA genomes** (figure footnote: "per-sample probe-weighted median-centred") is a known pitfall: when >50% of the genome is altered (several HGSC lines here have FGA>0.7), the median baseline sits on altered copy-number, systematically miscalling the neutral level and biasing FGA and gain/loss fractions.
- **The chrX artifact is disclosed in text but still displayed in `f_wes_cnv.png`,** where it produces the most eye-catching pattern in the figure (HGS lines blue/loss, non-HGS red/gain on chrX) — a pure pooled-normal sex-composition artifact. chrY was dropped; chrX should be too. As shown, the figure contradicts its own caveat.

**Remedy:** (i) drop chrX from `f_wes_cnv.png` (as chrY already is) or grey it with an explicit "artifact" annotation; (ii) confirm capture-kit concordance or, better, rebuild the reference from a kit-matched panel-of-normals; (iii) state that FGA is autosome-restricted (the ~2.9 Gb `total_assessed_mb` in `wes_cnv_fga.csv` implies it already is — say so) and caveat median-centering for high-FGA lines; (iv) report CNV frequencies at the patient level (§3.1).

### 3.5 Subtype/site confounding is only partially deconvolved on RNA *(moderate)*
The confounder defense rests on per-PC R² for subtype vs. site. But `rna_pc_confounder.csv` reports **only univariate marginal R²** (PC1: subtype 0.735 vs site 0.313), whereas the protein table (`prot_pc_confounder.csv`) correctly adds a **joint/adjusted R²** (PC1 adj: subtype 0.464 vs site 0.315 — much closer). Because all HGSC are Mes-Masson, subtype and site are correlated, so the RNA marginal comparison overstates subtype's independent contribution. The genuinely persuasive evidence is visual and narrow: in `f_rna_pca_subtype.png`/`f_rna_pca_site.png`, **clear cell is the only subtype present at both sites, and it co-clusters by subtype** (TOV3392D/TOV21G [Mes-Masson] sit among VOA10816/VOA295/VOA12539 [Huntsman]). That is a real internal control, but it is n=7 and it is the *only* place site can be separated from subtype at all.

**Remedy:** (i) add a joint model / `variancePartition` (subtype + site + patient-family as random effects) for the RNA PCs, matching the protein table; (ii) state plainly that HGSC subtype is fully confounded with site and that clear cell is the sole cross-site control; (iii) lean the "biology not batch" claim on marker recovery (lineage-specific, not a plausible batch artifact) and the CC cross-site control, as the report already half-does.

### 3.6 Validation criteria are lenient and under-specified; n≤2 subtypes are anecdotal *(moderate)*
- The "**16/22 markers land in the expected subtype**" headline uses a `lands_right` rule that (from `rna_markers_summary.csv`) is effectively "intended subtype is in the **top 2 of 6** subtype means" (e.g. WT1 lands "right" though its top subtype is SCCOHT and HGS is only rank 2; MUC2 is FALSE despite rank 2 because it is unexpressed). Being in the top-2-of-6 is a ~1/3 chance criterion, and several "subtypes" are n=1–2 lines so a single line sets the mean. The count is defensible directionally but the criterion must be stated, and ideally replaced/supplemented with a marker-vs-rest effect size (or AUC), plus the per-line heatmap (`04_rna_markers_heatmap`) as the honest display.
- **GO/pathway "recovery," FGA "ordering," and silhouettes for MMMT (n=2), SCCOHT (n=2), EC (n=2→1 after reclassification), MC (n=3, WES n=1), LGS (WES n=1)** are anecdotal. The FGA "ordering HGS>CC>MC>EC>LGS" rests on n=1 for MC/EC/LGS and n=2 for CC — and the two CC lines differ wildly (TOV3392D 0.688 vs TOV21G 0.121), so a CC "mean 0.40" is close to meaningless.

**Remedy:** state the marker criterion explicitly and add effect sizes; downgrade all n≤2 subtype "recovery"/"ordering" statements to per-line descriptions (show the points, not the subtype mean); avoid the word "ordered" for FGA across singleton subtypes.

### 3.7 TMT ratio compression is not carried into interpretation; bridge design is a chain *(moderate)*
The protein dynamic range is heavily compressed — FOLR1 spans 11.2–13.6 log2 at protein vs 0–9.5 at RNA (`adc_expression.csv`) — the classic TMT ratio-compression signature. This (a) caps achievable RNA–protein concordance (relevant to interpreting the 0.41), and (b) makes cross-line protein "shortlists" in the ADC atlas less discriminating than the RNA. Separately, `prot_bridge_cor.csv` shows only **4 bridge measurements forming a chain** (plex1↔2↔3↔4↔5, one of them the external VOA3993), not a common-reference design; Pearson≈0.99 per link is excellent, but cross-plex comparisons between distant plexes propagate through up to three bridges, so error can accumulate. Structural per-plex missingness (confirmed: `pct_missing` is identical within a plex in `prot_sample_qc.csv`) further means some proteins are block-missing across plexes.

**Remedy:** state the ratio-compression caveat where the 0.41 and the ADC shortlists are interpreted; describe the inter-plex normalization actually used (bridge-based?) and how block-missing proteins are handled in cross-line comparisons; note the bridge is a chain, not a hub.

### 3.8 An unflagged hypermutator (TOV21G, 1416 coding variants) *(moderate — could be a feature or an artifact)*
TOV21G carries **1,416 coding variants — >3× the next line** (413) — and 7 driver hits including multiple SWI/SNF genes (ARID1A×2, ARID1B, SMARCA2) plus KRAS/PIK3CA/PTEN/CTNNB1 (`wes_mutations_filtered.csv`; visible as the towering TMB bar in `f_wes_oncoplot.png`). This is either a real **MSI/mismatch-repair-deficient** clear-cell line (clinically meaningful — MSI occurs in clear-cell/endometrioid OC and would be a notable reuse feature) or a tumor-only artifact. Either way it should not pass unremarked.

**Remedy:** check MSI status / MMR-gene alterations and mutational context (even a caveated SBS/indel-signature look), and either report TOV21G as an MSI model or flag it as a QC outlier.

---

## 4. Additional analyses / figures / statistics to consider

**Essential (needed for a sound descriptor):**
1. **Patient-level re-tabulation of all genomic frequencies + family annotation tracks** (see §3.1). Payoff: removes the single biggest methodological objection; small effort (re-aggregation + one annotation row per figure).
2. **STR + mycoplasma statements** (see §3.2), even if only citing originators' records. Payoff: satisfies a hard journal requirement.
3. **Somatic-confidence tiering of WES drivers** and removal/annotation of germline-VAF and biologically-incoherent calls (§3.3). Payoff: makes the genomics defensible.
4. **Drop chrX from the CNV figure and confirm the CNV reference** (§3.4). Payoff: removes a self-contradiction and a validity risk.
5. **External identity/expression benchmarking against CCLE/DepMap and Cellosaurus.** Several lines are in public resources (e.g. TOV21G, TOV112D, COV434, OV90, OV-90; the CHUM TOV/OV lines and COV434 have CCLE/Cellosaurus entries). Correlating this RNA-seq (and mutation calls) against CCLE for the overlapping lines is the strongest possible **technical-validation and identity check** for a descriptor, and directly tests the discordance flags. Payoff: high; moderate effort. The `ConsensusOV` subtype calls already in `samples.csv` should also be surfaced as an orthogonal, external subtype validation.

**Nice-to-have (strengthen reuse / robustness):**
6. **`variancePartition` on RNA** with subtype + site + patient-family (§3.5). Payoff: turns the confounder argument from marginal-R² into a defensible decomposition.
7. **MSI / mutational-signature analysis** (TOV21G and panel-wide, caveated for tumor-only) (§3.8). Payoff: potential reuse feature; identifies artifacts.
8. **Passage-number sensitivity check.** `samples.csv` has RNA and WES passages (some differ substantially, e.g. TOV112D RNA p63 vs WES P83); a quick test that passage does not drive top PCs, and disclosure of within-line passage mismatch across assays, would pre-empt a reviewer question. Payoff: cheap insurance.
9. **Subtype-stratified and patient-collapsed concordance** (one line per patient; report per-subtype medians — HGS lines visibly sit lower in `f_concordance.png` panel B). Payoff: honest, and interesting.
10. **A proper genomic-HRD path or an explicit closed decision.** Either recover BAMs → Sequenza/ASCAT → allele-specific CN → HRD scar score (the report scopes this), or state clearly it is out of scope. LOH-segment counts from total CN are *not* a substitute and should not be offered as one.
11. **A consolidated per-line QC/metadata supplement** (pseudoalignment %, library size, genes detected, proteomic missingness/plex/channel, WES passage, variant count, purity assumption, ConsensusOV call, provenance, STR status). Most of this exists in `output/` already — surfacing it as one table is the highest-value low-effort reuse artifact.
12. **Deposition plan made concrete:** RNA FASTQ/counts (SRA/GEO), WES BAMs/VCFs (SRA/EGA, with the third-party normals stripped), proteomics raw + search results (PRIDE/MassIVE), and processed matrices; add the `renv` lockfile (noted as still-to-do). Payoff: this is what makes it a *Data* paper.

---

## 5. Minor points
- **Denominator inconsistency:** the oncoplot renders TP53 at "82%" (over ~22 all-subtype lines) while the text says "17/17 = 100%" (HGSC only). Both are correct; state the denominator in each place.
- **Gene-count mismatch in the concordance figure:** panel A says "n = 7893 genes," panel B says "8212 shared genes"; the text says 8,212 (97.4%). Reconcile (per-gene correlations require cross-line variance, hence fewer) in the caption.
- **"BRCA2 68→9%"** in the text vs. **5.9% (1/17)** in `wes_driver_freq_by_subtype.csv`. Reconcile the number and the denominator.
- **Weak "secondary SWI/SNF" calls:** OV2085 is called `swisnf_deficient=TRUE` on a lone SMARCA2 mRNA rank-3 z-score; OV2295 on a truncating SMARCA2 call that coexists with *high* SMARCA2 protein (z=+1.1) — internally inconsistent. Relabel these as "subunit variant of uncertain significance," clearly separated from the four reclassification-grade cases (TOV112D, COV434, BIN67, and the flagged VOA4841).
- **Protein-level cohesion:** the text says clear cell is "not cohesive" at protein level (silhouette −0.003, verified); note that SCCOHT is also non-cohesive at protein (0.028) — "HGS carries it" is fair but should mention SCCOHT too.
- **FOLR1 assay dependence:** the text says "highest FRα expressers panel-wide are HGSC (top TOV3133G)" — true for **RNA**; at **protein** the top line is a clear-cell line (VOA14993, 13.6). Specify the assay; the bimodality-within-HGSC claim (RNA) is otherwise well supported.
- **"Atlas" for n=1–2 subtypes** is grandiose; "panel"/"survey" is more accurate for MMMT/EC/MC/SCCOHT.
- **Define the audit criteria:** `expression_consistent` and `genomics_consistent` (and the marker `lands_right` rule) should have one-line definitions in the methods; the 26/20 counts are reproducible but the reader can't see the rule.
- **MKI67** appears in `rna_markers_summary.csv` unscored (correctly, as a proliferation control) — say so.
- **Silhouette space/metric** (which PCs, Euclidean?) should be stated.

---

## 6. Verification notes (what I recomputed and whether it matched)

Everything I checked reconciled with the report unless noted; the concerns above are interpretive, not arithmetic.

- **Cohort counts — MATCH.** Reconstructed from `samples.csv`: 42 generated; RNA 31 (HGS 15/CC 7/EC 2/MC 3/MMMT 2/SCCOHT 2), proteomics 31, WES 23 (22 with MAFs, TOV3121D "no MAF"), all-three 13 (all Mes-Masson). Confirmed RNA-set and protein-set both n=31 but *not identical* (VOA6861 RNA-only; VOA14993 protein-only).
- **Concordance — MATCH.** `integ_rnaprot_cor_summary.csv`: per-line Spearman median 0.4076 (report 0.41), per-gene median 0.3971 (report 0.40), n=30 lines.
- **PC confounder — MATCH (with caveat).** `rna_pc_confounder.csv` PC1 subtype 0.735 / site 0.313, PC2 0.694/0.214, PC3 0.582/0.049 — as reported. But RNA gives only marginal R² (no joint model), unlike `prot_pc_confounder.csv` which has adjusted R² (§3.5).
- **Silhouettes — MATCH.** `rna_silhouette.csv`: HGS 0.162, CC 0.121, EC −0.014, MC 0.154, MMMT 0.738, SCCOHT 0.818 — exactly as text.
- **FGA — MATCH numerically.** `wes_cnv_fga.csv` (fga_0.2): HGS mean ≈0.625 (report 0.63), CC 0.405 (0.40), MC 0.356 (n=1), EC 0.265 (n=1), LGS 0.073 (n=1). Note the non-HGS "subtype FGA" are singletons/pairs (§3.6).
- **Driver frequencies — MATCH the table, but pseudoreplicated.** `wes_driver_freq_by_subtype.csv`: TP53 HGS 17/17, CDK12 6/17 (35.3%). Confirmed via `wes_mutations_filtered.csv` that these counts double-count sublines with identical variants; patient-level HGS = 11, TP53 11/11 (100%), CDK12 3/11 (≈27%) (§3.1).
- **Coding variants/line — MATCH.** Median 206 (report "~206"); range 133–1416; TOV21G = 1416 outlier (§3.8). All rows in the filtered MAF are coding-consequence.
- **germline_like_vaf — flag is not doing the implied work.** 495/6036 variants (8.2%) and 21/50 drivers flagged TRUE, including canonical somatic TP53 hotspots (§3.3).
- **CNV arm frequencies — approximately MATCH; patient-robust.** My crude recompute (log2c>|0.3|): 3q26 gain 15/18 lines (83%) → 9/11 patients (82%); 17p loss 15/18 (83%) → 10/11 (91%). Report's 89%/78% are in the same range (threshold/region-definition dependent).
- **chrX — present in `wes_cnv_segments.csv` (194 segments) and displayed in `f_wes_cnv.png`** despite being called an artifact (§3.4). FGA's ~2.9 Gb `total_assessed_mb` implies FGA itself is autosome-restricted (good, but unstated).
- **Marker recovery — MATCH (16/22)** in `rna_markers_summary.csv`; criterion is lenient (top-2-of-6) (§3.6).
- **ADC associations — MATCH.** `adc_subtype_summary.csv`: MSLN→HGS and ERBB2→CC/MC recover in both assays (`top_mean_matches=TRUE`); FOLR1→HGS is `FALSE` by subtype mean (rank 3) but the max line is HGS — consistent with the "bimodal within HGSC" narrative (verified in `adc_expression.csv`: HGS FOLR1 RNA spans 0.06–9.5).
- **SWI/SNF panel — MATCH.** `auth_swisnf_panel.csv`: TOV112D SMARCA4 z=−2.14 + truncating; BIN67 SMARCA4 mRNA rank 23 / protein rank 2 (post-transcriptional); COV434 both low; VOA4841 lowest SMARCA4 (RNA & protein). Secondary HGSC calls (OV2085, OV2295) are weak/inconsistent (§4 minor).
- **Authentication audit counts — MATCH.** `auth_perline_table.csv`: 26 `expression_consistent="consistent"`, 20 `genomics_consistent="consistent"`; all `STR_status=NA` (§3.2). Mucinous flags (`auth_mucinous.csv`): TOV2414 ovarian_index +1.26 (externally authenticated), VOA8762 −0.76, VOA8771 −1.91 (GI-leaning) — as reported.
- **Proteomics QC — MATCH.** `prot_bridge_cor.csv` Pearson 0.991–0.994 (report ≈0.99); `prot_sample_qc.csv` missingness ~8–10%, *identical within plex* → confirms structural per-plex missingness; `prot_abundance_matrix.csv` 8,430 proteins × 31 lines, compressed dynamic range (§3.7).

---

*Bottom line:* a valuable resource and a careful, honest analysis whose descriptive/technical-validation core is largely sound, but which currently over-counts sublines as independent patients, over-labels molecular corroboration as "authentication," and carries a few indefensible tumor-only driver calls. Address §3.1–3.4 (mostly re-tabulation, re-wording, and one figure edit), add the CCLE/consensusOV external checks and the STR/mycoplasma statements, and this becomes a solid *Scientific Data* Data Descriptor.
