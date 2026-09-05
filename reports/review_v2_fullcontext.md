# Full-context accuracy audit — `docs/manuscript/v2/OvCAN_data_descriptor_v2.md`

**Auditor:** full-context reviewer (code, outputs, blind referee report, code audit, both fix reports, both figure reports, integration certification, and the curated writer package).
**Method:** every substantive quantitative claim traced to a deposited CSV or to the code that writes it; recomputed with `Rscript` where feasible; 8 of the 14 figures opened as images (`fig1`–`fig6`, `figs2`, `figs6`) and checked against their legends and source files. Read-only throughout; no manuscript or script was edited.

---

## 1. Verdict

The draft is substantially accurate: I checked roughly 200 individual numbers and the overwhelming majority reproduce exactly from the deposited outputs, including every composition count, every protein feature count, the full variant cascade, all ten arm frequencies, the commonality decomposition, the silhouettes, the marker effect sizes and the bridge-agreement statistics. The repairs are reported at their *after* values everywhere I checked, the reproducibility claim in Code Availability is scoped correctly, and the figure reports' "misreading, no fix needed" items are handled without inventing corrections.

But nine claims are wrong, and the two most serious are inferential rather than typographical: a **residual-germline estimate that is actually the median coding burden relabelled**, and a **stated impossibility about multiple-testing that the manuscript's own result refutes two sentences earlier**. Five of the nine trace to errors in the curated evidence package that the writer had no way to detect; four are the writer's own. Separately, the figure legends (7,882 words, 32% of the document) largely duplicate caption blocks already printed inside the panels.

**Recommendation:** fix the nine errors before the length reduction, because three of them will otherwise be preserved by a trim that cuts on length alone.

---

## 2. Critical errors

Each is a wrong number, wrong unit, wrong method description, or invalid inference.

### C1. "~206 residual germline coding variants per model" has no source and is the median coding burden relabelled

**Locations:** Methods §Whole-exome sequencing; Table 7 row 7; Usage Notes §Reuse limitations ("Somatic calls are tumour-only…").
**What it says:** "Residual germline contamination is estimated at approximately 206 coding variants per model, **the same order as the median coding burden itself**"; Table 7 row 7 cites `wes_pipeline_parameters.csv`.
**What is true:** No residual-germline estimate exists anywhere in `output/`. `wes_pipeline_parameters.csv` has 15 rows (CNVkit, capture, panel of normals, Mutect2 matched-normal status, build, proteomics) and contains no such quantity. The only "206" in the pipeline is `med_line = 206.5` at `scripts/07_wes_mutations.R:529` — the median coding burden, used as an assertion constant. `reports/review_code.md:251` verifies the previous draft's "~206 coding variants per line" against exactly that: "✅ median 206.5 (133–1,416)". The one germline-adjacent quantity that does exist, `germline_like_vaf`, flags 495 of 6,036 calls (≈22.5 per model) and was deliberately retired as a somatic filter.
**Why it is invalid, not just unsourced:** the sentence compares the estimate to the median burden and finds them "the same order" — but the number *is* the median burden. The comparison is tautological, and it presents an unmeasured quantity as measured, in a Methods section and in the limitations table a reuser is told to check.
**Verification:** confirmed. `grep` over all 102 CSVs; `wes_pipeline_parameters.csv` printed in full; `07_wes_mutations.R` inspected.
**Fix:** delete the number. Replace with: "No matched normal was sequenced, so the residual germline fraction is not estimable; the median coding burden of 206.5 candidates per model is an upper bound on it, and the tiering scheme is the mechanism for handling this." Correct Table 7 row 7's check-file to `wes_mutations_filtered.csv` (the empty `n_depth`/`n_ref_count`/`n_alt_count` columns) and drop the "~206" cell.

### C2. "No n = 2 marker can reach adjusted significance once more than about nine markers are tested" — refuted by the manuscript's own result

**Location:** Technical Validation §Histotype marker recovery, second paragraph.
**What it says:** `SMARCA2` "achieves exactly that floor (p = 0.00529 on n = 2 versus 26) … it cannot do better, and **no n = 2 marker can reach adjusted significance once more than about nine markers are tested**."
**What is true:** `SMARCA2` **does** reach adjusted significance. Recomputed: `p.adjust(BH)` over the 25 raw Wilcoxon p-values gives SMARCA2 BH = **0.04375**, and the deposited `rna_marker_effectsizes.csv$wilcox_p_bh` stores 0.0437. The manuscript states this itself two sentences earlier: "**4 of 25 survive Benjamini–Hochberg correction at 0.05** (`HNF1B`, `SPP1`, `KRT20`, **`SMARCA2`**)". The "~9" bound is the **Bonferroni** bound (0.05 / 0.00529 = 9.45), not the BH bound; under BH a marker at the floor survives when other markers are comparably small, which is what happened (KRT20 0.00242, SPP1 0.00385, SMARCA2 0.00529, HNF1B 0.00700 all adjust to 0.0437).
**Verification:** confirmed by recomputation; the four raw p-values and the stored BH values agree to 4 dp.
**Source of the error:** `EVIDENCE_DOSSIER.md:1897` — "no n = 2 marker can ever reach BH significance if more than ~9 markers are tested". Input error, faithfully reproduced.
**Fix:** "…it cannot do better. Under Bonferroni correction such a marker could not survive beyond nine tests; under the Benjamini–Hochberg procedure actually used it survives here only because three other markers have comparably small p-values, which is itself a reason to lead with the effect size and AUC."

### C3. BIN67's *SMARCA4* protein margin is not "to the next-lowest model"

**Locations:** Technical Validation §Model identity and annotation support; Table 5 SCCOHT row; Fig. 5 legend panel F.
**What it says:** "Because the margin **to the next-lowest model** is small (0.105 log2), the SCCOHT protein calls are consistent with, rather than proof of, complete loss."
**What is true:** From `auth_swisnf_long.csv`, the three lowest *SMARCA4* protein values are **VOA4841 13.30634 (rank 1)**, BIN67 13.31522 (rank 2), COV434 13.42020 (rank 3). BIN67's margin to the next-lowest model (VOA4841) is **0.0089 log2**, not 0.105. The 0.105 is the gap to the model **above** BIN67 (COV434). The code confirms the direction: `scripts/32_fig5_rare.R:353` computes `prot_sorted$value[d_rank + 1] - prot_sorted$value[d_rank]` — the gap upward — while line 356 labels it "margin to the next lowest line".
**Why it matters:** this is the case the manuscript calls "the strongest justification in this resource for depositing both expression layers". The corrected reading is arguably *stronger* — BIN67 and VOA4841 sit together at the panel floor, separated from the third-lowest model by 0.105 — but the sentence as written asserts a fact about the wrong pair, and it conceals that the panel minimum belongs to a clear-cell model (VOA4841), not to BIN67.
**Verification:** confirmed against `auth_swisnf_long.csv`, `auth_swisnf_panel.csv`, the rendered `fig5.png` panel F (which prints VOA4841 at rank 1 and BIN67 at rank 2), and the code.
**Source:** the code's own label, propagated through `FIGURE_SPECS_B.md:329`, `fig_report_B.md:41`, `fix_report_wes.md` item 10d and `integration_report.md:343`. The writer could not have caught it.
**Fix:** "BIN67's *SMARCA4* protein sits 0.009 log2 above the panel minimum (`VOA4841`, itself called *SMARCA4*-deficient), and the two lowest models are separated from the third-lowest by only 0.105 log2, so the SCCOHT protein calls are consistent with, rather than proof of, complete loss." Correct the label in `32_fig5_rare.R:356`.

### C4. "Protein interquartile ranges are 4–7× narrower … for every target" — false, and on the wrong denominator

**Location:** Usage Notes §Reuse leads, first paragraph.
**What it says:** "protein interquartile ranges are **4–7× narrower** than RNA interquartile ranges **for every target** (for example `TACSTD2` 8.27 → 1.16 log2; `ERBB2` 0.85 → 0.20)."
**What is true:** the folds are **3.52–9.44×** on the 31-models-per-layer basis and **3.84–9.64×** on the 30 dual-assay models. Two targets fall outside 4–7× on either basis (CD276 ≈ 3.5–3.8×; MSLN 9.44–9.64×), and FOLR1 (3.61×) falls outside on the 31-model basis. MSLN's 9.64× is stated by the manuscript itself in Technical Validation and printed in Fig. 6C. Separately, the two example values (8.27, 0.85) are the **31-model** IQRs; the Fig. 6A legend's IQR column gives the **30-model** values (TACSTD2 8.20 → 1.17; ERBB2 0.88 → 0.19), and Table 3 states that the 30 dual-assay models are the denominator for "**Every** transcript-versus-protein statistic".
**Verification:** confirmed by recomputation on both denominators from `adc_expression.csv` and `prot_dynamic_range.csv`.
**Source:** `EVIDENCE_DOSSIER.md:2449` states "The protein IQRs are 4–7× narrower … for every target" with those two examples, above a table from which the range 3.5–9.4× is directly computable. Input error.
**Fix:** "protein interquartile ranges are 4–10× narrower than RNA interquartile ranges for every target (`TACSTD2` 8.20 → 1.17 log2; `ERBB2` 0.88 → 0.19; range 3.8× for `CD276` to 9.6× for `MSLN`, on the 30 dual-assay models)."

### C5. TOV2414's SATB2 is not higher than *either* gastrointestinal-leaning model

**Location:** Technical Validation §Model identity and annotation support, mucinous paragraph.
**What it says:** TOV2414's "measured SATB2 mRNA z here is +0.17, rank 19 of 31: mid-panel, and **higher than either gastrointestinal-leaning model**."
**What is true:** `auth_mucinous.csv` / `auth_mucinous_marker_ranks.csv`: TOV2414 SATB2 z = +0.17 (rank 19), VOA8762 = −0.59 (rank 15), **VOA8771 = +0.80 (rank 23)**. TOV2414 is higher than VOA8762 only. The manuscript's own Fig. 5E legend gets this right: "mid-panel and *higher* than `VOA8762`'s."
**Verification:** confirmed; the two statements in the manuscript contradict each other.
**Source:** `EVIDENCE_DOSSIER.md:2247` — "higher than **both** other mucinous models (VOA8762 −0.59, VOA8771 +0.80)" — which contradicts its own parenthetical. `FIGURE_SPECS_B.md:279` is correct. Input error.
**Fix:** "…mid-panel, and higher than `VOA8762`'s though lower than `VOA8771`'s." The argument (SATB2-negative is a literature IHC value, not a measurement here) is unaffected.

### C6. The permutation p is not at the attainable floor, and the Fig. 3B legend contradicts itself

**Locations:** Technical Validation §Expression structure; Fig. 3 legend panel B.
**What it says:** "p < 0.002, **the minimum attainable at 1,000 draws**"; the legend reads "p = 0.002, the minimum attainable at 1,000 draws, so p < 0.002".
**What is true:** `scripts/17_variance_confounders.R:201` computes `(1 + sum(null >= obs)) / (N_PERM + 1)` with `N_PERM = 1000`, so the floor is **1/1001 = 0.000999**. The deposited value is 0.001998002 = **2/1001**, meaning exactly one of the 1,000 permutation draws met or exceeded the observed statistic. The observed value is therefore *not* at the floor, and 0.002 is not the floor. The legend's "p = 0.002 … so p < 0.002" is internally contradictory as written.
**Verification:** confirmed by reading the code and by 0.001998002 × 1001 = 2.000.
**Source:** `EVIDENCE_DOSSIER.md:1580` — "`0.00200` is `1/1001` rounded — the minimum attainable p at 1,000 permutations". Input error (1/1001 rounds to 0.001, not 0.002).
**Fix:** "permutation p = 0.002 (1 of 1,000 draws reached the observed value; floor 0.001 at this number of draws)". The substantive claim — observed 0.42–0.53 against a null 95th percentile of 0.24–0.39 — is unaffected and remains well supported.

### C7. Table 2's segment counts are autosome-restricted but attributed to a whole-genome file

**Location:** Table 2, whole-exome row.
**What it says:** "`wes_cnv_segments.csv` 5,428 × 10, row = copy-number segment, across 23 models (**median 185 segments per model, range 88–481**)."
**What is true:** the file has 5,428 rows over all chromosomes; per-model counts are **median 192, range 92–504**. The 185 / 88–481 figures are `wes_cnv_fga.csv$n_segments_auto` — autosomes only. Supplementary Fig. S5's legend uses 185 / 88–481 **correctly**, because that figure plots autosome profiles.
**Verification:** confirmed; both quantities recomputed.
**Fix:** either quote 192 (92–504) in Table 2, or write "median 185 autosomal segments per model (range 88–481); 192 (92–504) including sex chromosomes".

### C8. "The five targets with no expected histotype … three of them peak in mucinous"

**Location:** Usage Notes §Reuse leads, histotype-association paragraph.
**What is true:** **six** of the nine targets have no `expected` histotype in `adc_subtype_summary.csv` (CD276, CDH6, DPEP3, SLC34A2, TACSTD2, VTCN1). Of those six, **two** peak in mucinous (CD276 in both layers, TACSTD2 in RNA) — three (target, layer) rows, not three targets.
**Source:** `EVIDENCE_DOSSIER.md:2492-2494` says "the five targets with no `expected` histotype" and then counts "Three of them (CD276, TACSTD2 RNA, **FOLR1**)" — but FOLR1 *has* an expected histotype (HGS), so it cannot be one of them. Input error.
**Fix:** "For the six targets with no expected histotype recorded, only descriptive statements with the n attached are warranted; two of them (`CD276` in both layers, `TACSTD2` in RNA) peak in mucinous on n = 3."

### C9. "Our counts are 13–19× higher for the three shared models"

**Location:** Technical Validation §Model identity and annotation support, external driver cross-check paragraph.
**What is true:** `external_depmap_burden.csv`: OV90 240 / 18 = 13.3×; TOV112D 133 / 7 = 19.0×; **TOV21G 1,416 / 568 = 2.5×**. The range across the three shared models is 2.5–19×, not 13–19×. That TOV21G is the closest of the three is itself informative — it is the corroboration the same sentence goes on to claim.
**Source:** `EVIDENCE_DOSSIER.md:2400` states 13–19× immediately above the numbers that refute it. Input error.
**Fix:** "our counts are 2.5–19× higher across the three shared models (13× and 19× for the two ordinary-burden models, 2.5× for the hypermutator)".

---

## 3. Major issues — accurate but misleading, unsupported as stated, or missing a necessary qualifier

### Mj1. The floor-artefact result is under-claimed and mis-framed as a single-target curiosity

Technical Validation says compression "does **not** always shrink when the floor artefact is removed: for `MSLN` … 9.64× and 7.07×, whereas `FOLR1` shows **the ordinary pattern** (3.24× / 4.01× / 3.61×) … for at least one target it deflates it." Usage Notes repeats "For `MSLN` the floor-insensitive statistics show **more** compression".

Recomputed: the floor-insensitive interquartile ratio implies **more** compression than the range ratio for **all eight** ADC targets (including FOLR1: 3.24× → 4.01×) and for **65.6%** of all 7,896 paired genes; panel-wide the median IQR ratio (0.298) is below the median range ratio (0.342). Fig. 6C plots this — every rust circle sits left of its open circle — and its own in-panel text makes the general claim. FOLR1 is not a contrast case; it shows the same direction less dramatically, exactly as `EVIDENCE_DOSSIER.md:1010` states ("where compression also grows but far less dramatically"). The manuscript followed `INTERPRETIVE_NOTES.md:138`'s phrase "the ordinary pattern" and dropped the dossier's correction.

This is a directional under-claim in a place where the general statement is stronger, simpler and shorter. **Fix:** "The floor artefact does not inflate the compression estimate; removing it *increases* the estimated compression for all eight targets and for two-thirds of all paired genes, most extremely for `MSLN` (range ratio 5.30× against interquartile 9.64× and standard-deviation 7.07×)."

### Mj2. "Only 8 of 42 models have a positive lineage programme" conflicts with the manuscript's own Table 5

`auth_perline_table.csv$expression_basis`: **10** models rest on a positive lineage programme (7 CC + 3 MC); 8 of those are called `consistent` and 2 `partial`. The 8 / 16 / 2 = 26 decomposition is correct *for the 26 consistent models*, and the manuscript says so. But the headline — repeated in Technical Validation and in Table 7 row 13 — drops the qualifier, and Table 5's own rows say "All 7 assessed [CC] rest on a positive lineage programme" and "All 3 [MC]". A referee adds 7 + 3 and gets 10.
**Fix:** "Only 10 of 42 models have any positive lineage programme, and only 8 of those clear the consistency bar."

### Mj3. The exome raw-read claim asserts an existence that is nowhere confirmed

Table 2 lists "Raw sequence reads" for the exome record and commits it to NCBI SRA; Data Records says "The exome is a full layer with both raw reads and processed calls, parallel to the other two." Table 7 row 6 says "tumour BAM files not archived". The blind referee's M1(b) asked directly whether the WES FASTQs exist, and the previous draft's own note read "confirm that FASTQs exist for SRA deposit, otherwise the WES record is processed-only". Nothing in the manuscript states that this is unconfirmed. `WRITING_NOTES.md` §2 records that the SRA choice was the writer's own and that the directive required the full-layer description — so this is a directive-driven claim, not a writer error, but it is an unevidenced existence claim about a deposit that has not happened.
**Fix:** add one clause to Table 2 or the Data Records lead: "⟨TO CONFIRM — PI: that exome FASTQ files survive; tumour BAMs are not archived, and if the reads do not survive the exome record is processed-calls-only.⟩"

### Mj4. The "proliferation axis" agreement statistic compares two unrelated quantities, and the comparator is undisclosed

Technical Validation reports three agreement statistics between the HGS Hallmark strata and "an independent pathway-scoring method": "inflammatory axis … ρ = 0.818", "hypoxia axis … ρ = 0.632", "and **the proliferation axis not at all (ρ = −0.061, p = 0.83)**". The first two are `inflammatory_z` vs `progeny_NFkB` and `hypoxic_glycolytic_z` vs `progeny_Hypoxia` — matched pathways, computed in `37_supp_rnaprot.R:679-682`. The third is `proliferation_z` vs **`progeny_p53`** (`INTERPRETIVE_NOTES.md:329`, `EVIDENCE_DOSSIER.md:2009`); PROGENy has no proliferation pathway. It is not computed in any script or deposited file.

Reporting a p53-activity score as the comparator for a proliferation theme, under a sentence pattern that establishes matched-pathway comparisons, invites the reader to assume a matched comparator that does not exist. The manuscript separately and correctly notes that the proliferation theme was not used to name the strata — so a failed agreement on it is doubly uninformative.
**Fix:** delete the proliferation clause, or name the comparator and say why it is not a proliferation readout. Deleting it costs nothing: the two axes that matter are reported.

### Mj5. Referee M3(c) is only partially addressed — the largest centre effect in the analysis is never reported

The referee's complaint was selective reporting of the reassuring protein statistic. The manuscript addresses the isobaric-set half thoroughly. It does **not** report `prot_pc_confounder.csv`'s **site adjusted R² = 0.315 on protein PC1** (raw 0.361) — the largest centre R² anywhere in the analysis, larger than the RNA PC1 value of 0.313 that the manuscript *does* quote. The claim "Marginal R² for histotype exceeds marginal R² for centre on **every** component in both layers" is true (verified for all five components in both layers) but is illustrated only with the smaller number.

Also: "the isobaric set term … reaches **22–23%** adjusted R² on protein PC2 and PC5" — the values are **21.5%** and 22.9%. `EVIDENCE_DOSSIER.md:103` says 22–23%; input error, minor.
**Fix:** add six words — "(RNA PC1: 0.735 versus 0.313; protein PC1: 0.553 versus 0.361, adjusted 0.464 versus 0.315)" — and change 22–23% to 21–23%.

### Mj6. A computed, referee-requested null was silently dropped

Referee minor 8: "The top-2-of-6 marker rule has no stated null. With six subtype columns, a marker lands in the top 2 by chance one time in three, so 16/22 should be quoted against that null." The analysis computed it (`fix_report_rna.md` item 12: binomial against p = 1/3, **p = 0.0016** all graded, 0.0035 for up markers only; recorded in `rna_markers_summary.csv`). The manuscript quotes "Sixteen of 25 land in their intended histotype" with no null anywhere.
**Fix:** one clause — "sixteen of 25 land … against a one-in-three chance rate under the rank rule (binomial p = 0.0016)". This *strengthens* a result the manuscript currently reports as bare arithmetic.

### Mj7. `countsFromAbundance` is never stated

Referee minor 21 asked for it because it determines whether the DESeq2 counts are length-corrected. `scripts/01_rna_load_qc.R:123` sets `countsFromAbundance = "no"`. The Methods give `ignoreTxVersion = TRUE` and spend a full paragraph distinguishing `assigned_gene_counts` from `n_processed_fragments`, but omit the one tximport setting that governs what those counts are.
**Fix:** three words in the tximport sentence.

### Mj8. Tables S1–S4 have legends but are cited nowhere in the body

Body citation counts: Table S1 **0**, S2 **0**, S3 **0**, S4 **0**. All eight supplementary figures are cited (S1 ×3, S2 ×3, S3–S8 ×1). Referee minor 23 raised exactly this for Table S3, and the count has gone from one uncited supplementary table to four. `WRITING_NOTES.md` §7's citation check covered main Tables 1–7 only.
**Fix:** cite S1 in Data Records §Metadata, S2 with the arm frequencies, S3 with the marker recovery, S4 in Methods §Environment.

### Mj9. `EPCAM p.2_3fs` is quoted in violation of the manuscript's own rule

Methods: "**Only strings flagged `hgvsp_canonical` should be quoted.**" `wes_mutations_filtered.csv` records `hgvsp_canonical = FALSE` for `EPCAM p.2_3fs` — and Technical Validation quotes it twice. Every other HGVS string the manuscript quotes (`p.L639fs`, `p.Q192*`, `p.I195T`, `p.G244C`, `p.Y334F`, `p.L123fs`, `p.P1256del`) is flagged canonical; I checked all of them.
**Fix:** "`EPCAM` p.2_3fs (a reconstructed non-canonical string; the variant is a frameshift insertion at codon 2–3)".

### Mj10. Methods misattributes the "descriptive (n ≤ 3)" code-label defect

Methods: "(A label in **the figure-generation code** describes the descriptive class as 'n ≤ 3'…)". The label is in `scripts/03_rna_de_signatures.R` at lines 59–63, 160, 165–166 and 319 — the differential-expression analysis script, including its console messages. `grep` over the figure scripts `30`–`37` returns nothing. `WRITING_NOTES.md` §4.8 correctly identifies the DE script; the manuscript changed it. Beyond the misattribution, documenting an internal console-message defect in a published Methods section is questionable — see §7.

---

## 4. Minor issues and internal inconsistencies

| # | Location | Says | Should be |
|---|---|---|---|
| 1 | Data Records §The full analysis record; Code Availability; Table S4 legend | "13 text files (**12** of them per-script session records)"; "the **12** per-script session records"; "**12** `session_info_<script>.txt` files" | **11** per-script records exist (01, 02, 03, 04, 05, 06, 12, 13, 17, 19, 21). The 12th session-info file is `session_info.txt`, which is order-dependent — the very file the manuscript warns against. The 13th text file is `rna_variancepartition_dropped_genes.txt`. |
| 2 | Data Records §The full analysis record | "They group as: … (11); (24); (7); (4); (17); (13); (6); (6)" presented as a partition of the 102 CSVs | The counts sum to **88**, not 102. They are counts of *inventory rows* in `EVIDENCE_DOSSIER.md` §4.3, two of which are six-file wildcards; the RNA group has **25** rows (35 files), not 24. The grouping does not partition the deposit. |
| 3 | Technical Validation §Pipeline reproducibility | "with only **four** scripts exceeding one minute" | **Five**: `03` (520 s), `22` (571 s), `21` (253 s), `17` (243 s), `18` (**67 s**). Inherited from `integration_report.md` §7. Total 1,775 s verified exactly. |
| 4 | Technical Validation §Passage; Usage Notes; Table 7 row 4 | "differences reach **20 passages in both directions**" | −17 to +20. The maximum is +20 (TOV112D) and the minimum −17 (OV3331). Two of the three occurrences give the exact numbers adjacent; the Usage Notes occurrence does not. Fig. S2B's "spans −17 to +20" is precise. |
| 5 | Technical Validation §Site–histotype collinearity | "of the **four** multi-line families present, only family 3133 co-clusters" | Only **three** families have ≥2 RNA models (1369, 2295, 3133 — each 2). Family 3291 contributes a single RNA model (TOV3291G) and cannot co-cluster by construction. The Fig. S3 legend avoids the denominator and is correct. |
| 6 | Fig. 2E legend | "n = **836** features per decile" | 7 deciles of 836 and 3 of 835 (sum 8,357). Self-correcting because the denominator is given. |
| 7 | Methods §Environment | `mclust 6.1.2` cited; `package_versions.csv` described as "the correct single citation for the environment" | mclust is **not** among the 39 packages in `package_versions.csv`. Its version is recorded only in `session_info_13_adc_atlas.txt`. Either add it or soften the "single citation" claim. |
| 8 | Fig. 1 legend panel A | "Box fill distinguishes stage type only (**grey input/process**, rust RNA, teal protein, navy exome)" | The "Cell lines n = 42" box is rust-filled, identical to the RNA-seq box, not grey. |
| 9 | Fig. 1A vs Table 2 | Fig. 1A's deposition box reads `GEO · PRIDE · figshare` (three) | Table 2 commits to **four** repositories (adds NCBI SRA for exome raw reads). `WRITING_NOTES.md` §2 flags it; the manuscript does not. Add the fourth to the figure, or note the mismatch. |
| 10 | Throughout | Text says "models" and "isobaric sets"; figures say "lines" and "plexes" | The HGS/HGSC equivalence is stated once in Methods; these two are not. Fig. 1's "Unit throughout: cell line model" partly bridges the first. |
| 11 | Table 4 | HGS and MC rows give the nominal p in parentheses; the MMMT row does not | The MMMT grade "Suggestive at model level" rests on nominal p = 0.0457, which is not shown. |
| 12 | Methods §Structure | "an **equivalent** per-gene `lme4` REML decomposition is used instead" | The word "equivalent" is asserted without the agreement demonstration referee minor 20 asked for. The corrected description of *why* `variancePartition` fails ("installed but fails at run time … calls a function that has moved to another package") is right and is an improvement on the previous draft's wrong claim. |
| 13 | Fig. 5D legend | "on the **half-row below** `TOV21G`" | In the render the SBS6 interval sits on the TOV21G row. Cosmetic. |
| 14 | Data Records | No column-level data dictionary | Referee minor 24 asked for a data-dictionary supplementary table in the paper. Table 3 + Table S1's legend + §Metadata cover part of it; column definitions are still deferred. |
| 15 | Technical Validation §Model identity | Mycoplasma absence stated; no mention of the cheap partial substitute | Referee M14 suggested an unmapped-read mycoplasma screen, which the RNA-seq permits. Not done and not mentioned as an option. |

**Non-problems I checked and confirmed the manuscript handles correctly** (recorded so a trim does not "fix" them):
- The `figs6` intrinsic-stratum footnote is **already correct** in the rendered figure — it reads "ward.D2 hierarchical clustering, k = 3, on the 50 Hallmark set scores", with an explicit `[correction]` comment at `33_supp_genomics.R:336`. Only the stale `FIGURE_SPECS_B.md:422` still says "k-means". `WRITING_NOTES.md` §4.1's action item ("correct the in-panel text of `figs6`") is stale and **must not be executed**. The manuscript's S6 legend is right.
- Fig. 1's WES-SNV marks were never wrong (`fig_report_A` §1.7). The manuscript states the assertion rather than claiming a fix.
- MKI67 was never boxed (`fig_report_A` §1.6). The manuscript says it "was never treated as an expected marker" — correct, no invented correction.
- The Fig. 4C arithmetic the writer flagged as unreconciled (`WRITING_NOTES.md` §4.9) **closes**: 45 distinct gene × model cells, 5 carrying two calls (ARID1A/TOV21G, CDK12 × 3, TP53/TOV2881EP) = 50; the 5th is ARID1A/TOV21G with two frameshifts, one variant class, correctly not labelled multi-hit. Verified.
- The 6,689 / 6,688 distinction (compression vs correlation complete-case sets) is real and correctly assigned.
- The vendor-CV-versus-bridge-CV apparent discrepancy (`WRITING_NOTES.md` §4.10) is verified consistent: 100 × (2^SD − 1) with SD 0.430 → 34.7% and 0.106 → 7.6%.
- All repair *after* values are used: 8,427 (not 8,429/8,430), 70 (not 71), 7,733 (not 7,734), 22,542 (not 22,544 for the decomposition), 7,894 (not 8,212), 19 genomics-consistent (not 20), 6 SWI/SNF-deficient (not 7), DIF 5 / MES 5 / IMR 3 / PRO 2, 8 of 15 (not 7), `p.L639fs` (not `p.L639X`), SBS15 (not SBS20), 15.7 / 12.6 / 27.3% within clear cell (not 4–6%), 18–46 pp sweep (not 10–27).

---

## 5. Original-referee coverage

| Point | Status | Where / note |
|---|---|---|
| **M1** Abstract asserts deposition; unresolved placeholders | **Addressed** | "**No accessions have been issued**"; Table 2's ⟨ACCESSION PENDING⟩ ×4; future tense throughout; the seven bracketed notes replaced by ⟨TO OBTAIN⟩ markers each naming who closes it. |
| **M1(b)** Do the WES FASTQs exist? | **Not addressed** | See **Mj3**. Raw reads asserted; existence never confirmed; Table 7 row 6 says BAMs are not archived. |
| **M2** Proteomics not reproducible; missing versions and the arm threshold | **Addressed** | ⟨TO OBTAIN — Morin laboratory⟩ enumerates every missing parameter and states why acquisition mode matters. Arm rule stated (\|log2\| > 0.20 over > 50%, autosome median centring). COSMIC v3.2/GRCh38, MutationalPatterns 3.18.0, CNVkit 0.9.10, consensusOV 1.30.0, PROGENy 1.30.0, singscore 1.28.1, MSigDB h.all.v7.4 all given and verified against `package_versions.csv`. Sarek/GATK versions flagged as PI items. |
| **M3(a)** Median-only reporting; means tell the opposite story | **Addressed** | Mean, median and IQR reported for every term; the manuscript states plainly that the RNA means are "near-indistinguishable" and the IQRs "almost entirely overlap". |
| **M3(b)** Raw R² reported; no adjusted, no permutation null | **Addressed** | Raw *and* adjusted throughout; 1,000-draw nulls with means, 95th percentiles and p. (See **C6** for the p-floor mislabel.) |
| **M3(c)** Protein PC confounder table never referenced | **Partially** | Plex term reported (as "22–23%"; true values 21.5 / 22.9). Protein PC1 centre R² 0.361 raw / **0.315 adjusted** never reported. See **Mj5**. |
| **M3(d)** Fig. S3 uncited; CHUM/VOA split unaddressed | **Addressed** | Dedicated subsection "Site–histotype collinearity, seen directly"; the centre-wise top-level split is stated as the design, not a nuisance. |
| **M4** Pseudoreplication uncorrected for expression; patient term unidentifiable | **Addressed** | 28-representative sensitivity analysis for PCs, silhouettes, DE and enrichment; the patient term reported only with 28 levels / 31 observations / 3 replicated patients and the 0.76% restricted value; absence of a protein equivalent stated. |
| **M5** Stated null contradicted by the QC table | **Addressed** | +948 genes (+4.8%), Wilcoxon p = 0.0082, r = −0.63 and +0.50 all reported; depth described as "plausible but undemonstrated" with both non-significant p values. |
| **M6(a)** SWI/SNF claims cite a figure with no panel | **Addressed** | Fig. 5F added. |
| **M6(b)** BIN67 claim unquantified | **Addressed with a defect** | Quantified (rank 23/31 vs 2/31, z +0.62 vs −1.96) — but see **C3**. |
| **M6(c)** 5 × 5 submatrix cannot show rank 1 of 67 | **Addressed** | Fig. 4B rebuilt: 66 non-self ticks per row, named best non-self, Δρ and z per model. |
| **M7** CDK12 "matching the ~3% TCGA rate" | **Addressed by deletion, not over-corrected** | The TCGA comparison is gone; 6/22 models vs 3/16 patients reported with the full Tier composition and an explicit refusal to call it a somatic frequency. |
| **M8** Four uncited supplementary figures; passage confounder absent; bad Fig. S8 cross-ref | **Addressed** | All eight supplementary figures cited; dedicated passage subsection plus Table 7 row 4; the bogus cross-reference is gone. (Supplementary *tables* are still uncited — **Mj8**.) |
| **M9** Protein count wrong; 71 unquantified proteins; four feature counts unassigned | **Addressed and improved** | 8,430 − 3 = 8,427; 70 zero-plex characterised and retained; Table 3 assigns a denominator to each of the eleven counts. |
| **M10** Bridge "reproducibility" measures the wrong quantity; VOA3993 external | **Addressed** | Bland–Altman with bias, SD, LoA width and repeatability CV per link; CV by abundance decile; four different samples named; VOA3993's external status disclosed in the legend, the Methods and Table 7 row 20; Pearson r explicitly demoted. |
| **M11(a)** Fig. 6 protein panel cannot support the claim | **Addressed** | Row-scaled with an absolute IQR column and a stated clamp. |
| **M11(a2)** "FOLR1 is strongly bimodal" | **Addressed** | Tested and denied; "bimodal" appears only in the negative. |
| **M11(b)** Compression inflated by the RNA zero floor | **Addressed, then mis-framed** | Floor-insensitive statistics made primary and a floor-stratified check deposited — but see **Mj1** and **C4**. |
| **M11(c)** Fig. 6 claims carry no n or effect size | **Addressed** | Per-histotype n printed in-panel and in the legend; every "highest in X" statement carries its n. |
| **M12(a)** No figure or supplementary table legends | **Over-corrected** | 14 figure legends (7,882 words) and 11 table legends (1,163 words). See §7. |
| **M12(b)** n = 8,212 wrong for the median | **Addressed** | 7,894 used; 8,212 retained as a separately labelled quantity in Table 3 and the Fig. 3C legend. |
| **M12(c)** Per-gene n distribution undisclosed; 10.3% negative | **Addressed** | 84.7% at n = 30, 13 other values, complete-case sensitivity, 10.3% negative with the extreme value. |
| **M13** Mucinous de-authentication thin and circular | **Addressed except the external anchor** | Neutral column headers; z reference set stated and correct; 60-combination sweep; CDX2-alone for VOA8762; provenance decomposed; "expression can raise a flag but cannot call origin". No external TCGA/CPTAC anchor (`fix_report_wes` §4.5 declines it). See **C5** for the SATB2 comparison. |
| **M14** Culture conditions, ethics, mycoplasma, no EC markers | **Addressed** | Two ⟨TO OBTAIN⟩ markers; mycoplasma absence in the abstract and Table 7 row 15; EC marker set added and its failure reported at length with the culture-drift alternative left open. Mycoplasma read-screen not mentioned (minor 15). |
| **M15** "Corroborated orthogonally by PROGENy"; Fig. S8A circular; "mostly co-clustering" | **Addressed** | "consistency check with a different method on the same data, **not** an orthogonal assay"; partition-versus-labels distinction stated precisely; hypoxia label reported as not corroborated with three quantifications. Family co-clustering denominator is off by one (§4 item 5); the proliferation axis is problematic (**Mj4**). |

**Minor points:** 1–7, 9–19, 22 addressed. **8 not addressed** (binomial null dropped — **Mj6**). **20 partially** (§4 item 12). **21 not addressed** (`countsFromAbundance` — **Mj7**). **23 not addressed and worsened** (**Mj8**). **24 partially** (§4 item 14).

---

## 6. Figure legend accuracy

Checked against the rendered image and the source CSV. Every printed statistic in `fig1`–`fig6`, `figs2`, `figs6` reproduces from the deposited files unless noted.

| Figure | Verdict | Detail |
|---|---|---|
| **Fig. 1** | Accurate; two cosmetic defects | Counts, totals, family ramp order (1369 lightest → 3291 darkest), the `All 3` column and the TOV3121D assertion all verified. §4 items 8 and 9. |
| **Fig. 2** | **Accurate throughout** | Every one of ~35 printed values verified: medians 19,817 / 20,765 / 20,119; +948 (+4.8%) at −4.1 pp; p = 0.0082 (n = 30); r = 0.50 (p = 0.0038) and −0.63; depths 60.0 / 77.8 / 50.5 M with p = 0.37 and 54.4 / 68.3 M with p = 0.77; IQR medians 1.14 / 0.34; 98.0% of 7,896 and 99.2%; 70/264/360/366/512/6,855; vendor 11.0 → 3.0 and bridge 34.7 → 7.6; all four bridge facets' n / bias / SD / CV / LoA. The "836 per decile" rounding (§4 item 6) is the only imprecision. Panel C's normalisation-scale caveat is exactly right and load-bearing. |
| **Fig. 3** | **Accurate throughout** | PC1 20.7% / PC2 10.4%; commonality 42.4/39.3, 0.2/−2.5, 31.1/28.9, 73.7/65.7; n = 7,894 with the ≥10-model rule; 10.3%; all ten median/mean pairs in panel D; 22,542 / 6,855; the 0.76% restriction; 26 rows in seven blocks (5/5/3/6/4/2/1) with solid/dashed direction boxes; z reference set stated as the six histotype means. Only the panel-B permutation-p wording is wrong (**C6**). The explicit disclosure that the 0.38–0.48 band and the 0.72 line are external literature values not computed here is essential and must survive any cut. |
| **Fig. 4** | **Accurate throughout** | 25,914 → 493 → 170 and 557,392 → 15,692 → 6,036; ATM 18→2 (82%→9%), ATR 17→1 (77%→5%); all five Δρ and z; 45 cells / 50 calls with the five two-call cells correctly resolved; 7 dagger genes exactly matching the no-Tier-1 set; TP53 = 12 = 11 HGS patients + TOV112D; the six arm frequencies and the printed call rule; 44-row file with 5 acrocentric p-arms omitted → 39 plotted. |
| **Fig. 5** | Accurate except the SMARCA4 margin | 1,416 / 6.9× / z 21 / indel 0.279; median 206.5; 2,417 SNVs; cosines 0.877 / 0.651 / 0.315 and the top three 0.877 / 0.835 / 0.809; refit 0.733 vs 0.000–0.292, SBS6 0.428 [0.116, 0.454], 99.5%/200, cosine 0.977 (next 0.960); mucinous z values and the four RNA-only markers; panel F ranks, z limits ±3.19, tier text and the ‡ rule. **C3** is the one error. §4 item 13 is cosmetic. |
| **Fig. 6** | Accurate; contradicts the body text | All 17 IQR values verified exactly on the 30 dual-assay models; DPEP3 × glyph and the IQR = 0.00 caveat; clamp ±3 with observed max 4.94; per-histotype n; the 0.73 / 0.34 double denominator; FOLR1 BC 0.440, ΔBIC −1.06, range 0.06–9.48, 4 at mean 0.65 and 11 at 6.49; MSLN 0.189 / 0.104 / 0.141. Panel C's general claim about the floor-insensitive statistics is **more accurate than the Technical Validation text** — see **Mj1**. |
| **Supp. Fig. S1** | Accurate | PC1 12.7% / PC2 9.4% verified; all 14 silhouette values verified to the stored precision (protein CC −0.0035 → −0.004); n printed per histotype. |
| **Supp. Fig. S2** | **Accurate throughout** | 7.8 → 14.8 → 0.4 and 15.9 → 3.8 → 2.0 verified; the five `passage alone` p values (0.13 / 0.03 / 0.37 / 0.12 / 0.97); PC5 raw R² 0.00004; site median 3.53 → 0.00%; r = −0.33 [−0.75, 0.27] p = 0.26 recomputed exactly; −17 to +20 with median 4 and mean 5.8. |
| **Supp. Fig. S3** | Accurate | Off-diagonal range 0.757–0.969 median 0.834 verified against `fig_report_A` §1.1 and the untruncated scale; all 22,544 genes; the three-family co-clustering statement avoids the denominator error the body text makes (§4 item 5). |
| **Supp. Fig. S4** | Accurate | 32 patterns summing to 8,427 verified; the 1,572-includes-70 → 1,502 arithmetic verified. |
| **Supp. Fig. S5** | Accurate | 23 models; **185 (88–481) is correct here** because the panel is autosome-restricted — the same numbers are wrong in Table 2 (**C7**); TOV81D 0.021 as the panel minimum with TOV21G next lowest verified; TOV3121D as the sole ‡ verified; the pooled-normal caveat is load-bearing. |
| **Supp. Fig. S6** | **Accurate, and its method statement is already correct in the render** | 5/2/3/5; IMR+MES = 8/15; five margins < 0.10; median 0.160; the two input-set flips with correct directions; the one inherited-label disagreement. The in-panel footnote already says ward.D2, k = 3, Hallmark-only — see the non-problem note in §4. |
| **Supp. Fig. S7** | Accurate | Unit = 28 representatives with group sizes 12/7/2/3/2/2; three arrow-clipped upper bounds (9.5 / 19.0 / 20.6) verified; 4 of 25 BH-starred; SMARCA2 as the sole floor-limited marker; the signed-vs-oriented explanation. The legend does **not** repeat the body's false "~nine markers" claim (**C2**). |
| **Supp. Fig. S8** | Accurate except the proliferation comparator | Stratum sizes 4/5/6; ρ = 0.63 (p = 0.014) and 0.82 (p = 0.0003) recomputed exactly; 2 of 5 versus 1 of 4 above z = +1; inflammatory median +0.64 above hypoxic +0.56; max |z| 2.69. The partition-versus-labels distinction and "not an orthogonal assay" are both exactly right. **Mj4** applies to the body text, not this legend. |

---

## 7. Length-reduction guidance

Current: 24,399 words total — Abstract 183, Background 598, Methods 2,928, Data Records 1,620, Technical Validation 6,881, Usage Notes 2,635, Code Availability 491, **figure legends 7,882 (32.3%)**, table legends 1,163. Per-figure legends: Fig. 5 1,127, Fig. 3 1,025, Fig. 4 1,001, Fig. 2 927, Fig. 6 734, Fig. 1 341; supplementary 272–471 each.

### 7a. Load-bearing — must survive any cut

Removing any of these makes a surviving claim unsupportable or converts a calibrated statement into an over-claim.

- **Every unit and denominator.** 42 / 34 / 31 / 30 / 28 / 23 / 22 / 18 / 16 / 13 / 11 / 10, and **Table 3 in full**. Table 3 is 94 words of legend plus a table and it is the cheapest instrument in the document for keeping eleven feature counts and eight sample sets honest. Cutting it re-creates referee M9 in its entirety.
- **The double denominator for protein noise:** median IQR **0.344** log2 over 7,896 genes *and* **0.728** for the eight targets, with the derived bands **2.4–3.1×** and **1.1–1.4×**. Keeping only 0.344 understates the targets the paper reasons about; keeping only 0.728 understates the general compression. The manuscript says so explicitly — keep that sentence.
- **The within-clear-cell control in its full form:** 15.7 / 12.6 / 27.3%, E[R²] ≈ 16.7%, two of three below it, ANOVA and permutation p, and "no evidence … and equally no evidence against one". Trimming any element turns a correctly-framed underpowered null into either reassurance or alarm.
- **The structural-missingness paragraph:** identical detected-feature counts per set (7,738 / 7,758 / 7,562 / 7,694 / 7,653), 32 presence patterns, 8.796% overall, the 7.94–10.27% per-set spread, and "set-conditional presence must never be interpreted as biology". This is the single most consequential reuse caveat for the proteome.
- **Raw *and* adjusted commonality components**, the negative adjusted unique-centre values on all five components, the shared 0.311, and the Fig. S3 centre-wise split. These four together are the answer to referee M3.
- **Mean *and* median in the per-feature decomposition** (RNA 14.94 vs 14.50 means; overlapping IQRs). Deleting the mean is precisely the defect referee M3(a) identified.
- **The patient-term replication caveat:** 28 levels on 31 observations, 3 genuinely replicated patients, 0.76% restricted, and that dropping the term preserves the ordering.
- **Tier composition:** SMARCA2 / BRCA2 / NF1 Tier-3-only, CDK12 with no Tier 1, and the 7 → 6 SWI/SNF-deficient count with its reason.
- **Passage discordance** with its range, its median, and the concentration in two models one of which is the sole endometrioid model.
- **Arm frequencies with their threshold and the 18–46 pp / 18–28 pp sweep**, and the HGS-only 18-models / 11-patients unit.
- **Evidence class** (formal n ≥ 10, descriptive n < 10) and n per contrast.
- **Wilcoxon floors** (0.00529 / 0.000611 / 1.7 × 10⁻⁶) and **MUC5AC d = 1.92 [−0.44, 20.64]**.
- **Table 7 in its entirety.** Every row is a limitation, its number and its check-file. It is only 59 words of legend plus the table, and it is the densest honest content in the paper.
- **"No accessions have been issued"** and every future tense in Data Records.
- **Code Availability's offline-versus-certified distinction** and the figure-script scope limit. Both are precisely correct against `integration_report.md` and both are easy to lose in a trim.

### 7b. Redundant — same number or caveat in more than one place; keep the instance named

| Content | Where it repeats | Keep | Cut |
|---|---|---|---|
| **In-panel caption blocks** | Every printed value in `fig1`–`figs8` appears both inside the panel and in the legend | The legend (journals typeset legends; readers of the PDF get both) | The legend's *restatement* of panel text. `fig_report_A` §3.4 says the in-panel blocks are single `labs(caption=)` calls and "are the first thing to trim". Either way, choose one location. **Largest single saving available: 3,500–4,000 words.** Reduce Fig. 2–5 legends from ~1,000 to ≤250 words each by keeping only unit, n, test with its parameters, z/scale reference set, axis truncation or clamp disclosures, and one interpretive limit per panel. |
| Passage mismatch | Technical Validation §Passage; Usage Notes §Reuse limitations para 1; Table 7 row 4; Fig. S2B legend | Technical Validation paragraph + Table 7 row 4 | The Usage Notes paragraph ("Cross-assay analyses inherit a passage mismatch", ~55 w) |
| The five emphasised reuse limitations | Usage Notes §Reuse limitations prose (~330 w) restates Table 7 rows 1, 2, 4, 6, 7, 8 | Table 7 | The whole prose block; replace with one sentence pointing at Table 7 |
| "Single-protein, single-model differences are not interpretable" | Abstract, Background, Technical Validation, Usage Notes prose, Table 7 row 1 | Abstract + the Technical Validation derivation + Table 7 row 1 | Background and Usage Notes instances |
| 8,430 − 3 = 8,427 derivation | Methods §Proteomics, Table 2, Table 3, the Table 3 arithmetic line, Fig. 2D legend, Fig. S4 legend | Table 3 + Fig. 2D legend | Methods sentence, the arithmetic line, the Fig. S4 repeat |
| No STR / no mycoplasma | Abstract, Technical Validation, Table 7 row 15 | Abstract + Table 7 row 15 | The Technical Validation restatement (the surrounding Cellosaurus detail stays) |
| LGS has no expression data | Background, Table 1, Table 1 legend, Table 3, Table 5, Table 7 row 12, Fig. 3A legend | Background sentence + Table 1 | Table 1 legend's note, Table 7 row 12's LGS clause |
| The 23-vs-22 exome trap | Methods, Data Records §two traps, Table 3, Table 6, Table 1 legend, Fig. 1B legend | Data Records §two traps + Table 3 | Methods, Table 6 row, Table 1 legend clause |
| chrX pooled-normal artefact | Methods, Technical Validation, Table 6, Fig. S5 legend | Methods + Table 6 | The Technical Validation sentence (keep only "autosomes only is the metric to quote") |
| Bridge samples are four different samples, one external | Methods, Technical Validation, Table 7 row 20, Fig. 2F legend | Methods + Table 7 row 20 | Fig. 2F legend's four-sentence restatement (one clause suffices) |
| MSLN 9.64× | Technical Validation, Usage Notes, Fig. 6C legend | Technical Validation (in its corrected general form, **Mj1**) + Fig. 6C | The Usage Notes paragraph collapses to one sentence |

### 7c. Over-specified for what it establishes

- **Runtime** — "1,775 s (≈ 29.6 min), with only four scripts exceeding one minute". Cut entirely (and it is wrong; §4 item 3). Save ~20 w.
- **The canonical execution order plus the load-bearing-constraints list** in Code Availability (~130 w). Keep "numeric script order is not dependency order; the canonical order is recorded in the deposited integration report" and delete the chain and the constraint list. Save ~110 w.
- **The eight-way file-group breakdown** in Data Records (~90 w). Cut entirely — it is arithmetically broken (§4 item 2) and tells a reuser nothing actionable.
- **"Arithmetic that closes"** under Table 3 (~40 w). The table already closes; the Table 3 legend repeats it a third time.
- **The HGSC/HGS axis-label note** (Methods, ~25 w) and **the "descriptive (n ≤ 3)" code-label note** (Methods, ~35 w). These are production defects, not paper content. Harmonise the axis labels and fix the code label; delete both sentences. (Also removes **Mj10**.)
- **Three peptide totals** — Table 7 row 23's tail enumerating 146,830 / 146,629 / 158,885. Keep the row's protein counts; move the peptide enumeration to Table 3's existing note.
- **Fig. 1 legend's panel-A fill description and layout disclaimers** (~60 w; also inaccurate, §4 item 8).
- **Fig. S4B's justification for points rather than bars on a log axis** (~40 w) — a correct plotting argument, not a result.
- **Fig. S4A's full 32-value pattern list** (~65 w) — the counts are in `prot_feature_accounting.csv`; keep the first, the last and the sum.
- **Fig. 4C's column- and row-ordering rules** (~30 w) and **Fig. 1B's sort key** (~25 w).
- **Fig. 5F's deficiency rules restated in full** (~60 w) — they are already in Methods verbatim; the legend can point.
- **Molecular-subtype paragraph** — five margin values plus a median plus a minimum probability plus the class distribution plus the 8-of-15 count (~140 w). Keep "5 of 15 below 0.10, 2 of 15 unstable to the input set, 8 of 15 called to a microenvironment class" and drop the value list (it is in Fig. S6 and Table 7 row 14).

### 7d. Must NOT be cut although it looks like a candidate

- **All seven ⟨TO OBTAIN⟩ markers**, and in particular the proteomics one. Deleting the acquisition-mode marker while keeping the compression interpretation converts a bounded statement into an unsupported mechanistic claim — the manuscript's own sentence "*Acquisition mode bears directly on the magnitude of isobaric ratio compression, which this descriptor quantifies*" is what makes the compression discussion admissible.
- **The capture-kit marker's second half** — that it is unresolved whether the five pooled normals share the tumours' capture design. That is a **validity** caveat on the entire copy-number layer, not a documentation gap, and it is the reason no trinucleotide renormalisation was applied.
- **Every n = 2 / n = 3 annotation** on a silhouette, effect size, histotype mean or FGA value. These read as clutter and are the only thing preventing the rare-histotype numbers from being over-read.
- **"Screening, not attribution"** and the SBS1/5/6/15/44 mutual-similarity sentence.
- **"Not an orthogonal assay"** for the PROGENy check and **"the partition is not circular but the labels are"**.
- **The circularity disclosure for the two SCCOHT annotation calls** and the 8 / 16 / 2 decomposition (in its corrected form, **Mj2**).
- **Fig. 3C's disclosure that the 0.38–0.48 band and the 0.72 line are external literature values not computed here.** These two annotations still lack citations (`WRITING_NOTES.md` §5.9); the disclosure is the only thing keeping them from reading as this resource's own benchmarks. If the citations cannot be supplied, remove the annotations from the panel — but never remove the disclosure while keeping the annotations.
- **The two disclosed asymmetries in the DepMap comparison** and "specificity, not magnitude, is the signal".
- **The ad-hoc / uncalibrated framing of the mucinous ovarian index** and its verbatim definition, plus the hand-chosen thresholds.
- **"There is no dependency lockfile"** and the distinction between recording versions and pinning an environment.
- **The abstract's closing sentence** on STR and mycoplasma. It is the qualifier that stops "technical validation quantifies what each layer supports" from reading as a blanket quality claim.
- **Table 5's evidence-basis column.** Without it, "26 expression-consistent" is exactly the bare aggregate referee M-minor 15 objected to.

### 7e. Sequencing advice

Fix §2 and §3 **before** trimming. Three of the nine critical errors (**C1**, **C2**, **Mj1**) are in sentences a length-driven trim would plausibly keep because they are short and assertive, and two of them (**C1**, **C2**) currently read as the paper's most confident statistical statements. Cutting §7b's figure-legend redundancy first is safe and recovers 3,500–4,000 words — enough that no load-bearing body content need be touched at all.

---

## 8. What I could not verify

| Item | Why | What would settle it |
|---|---|---|
| Variant cascade stages 557,392 / 15,692 / 15,609 | Derived at render time from `judy_archive/data/wes - old/mutect2/`; not deposited as a table | Already independently recomputed by `review_code.md` §M-verification and re-checked in `integration_report.md` §6 (25/25 PASS). I accepted that evidence rather than re-deriving it. |
| Fig. 4A's pre-filter ATM (18) and ATR (17) model counts | Same — computed at render from the raw MAFs | Re-run `31_fig4_genomics.R` and read its console assertions. |
| Peptide total 146,830 | The peptide table is not in `output/`; only `prot_feature_accounting.csv` records the count | The archived `peptide_ratio.xlsx` row count. |
| Total fitted mutations 2,339.6 for TOV21G | Not in `wes_signature_refit_summary.csv`; presumably in `wes_signature_refit_exposures.csv` | I confirmed the conversion is internally exact (271.2/2,339.6 = 0.116; 1,061.7/2,339.6 = 0.454; 0.428 × 2,339.6 ≈ 1,001), which makes a wrong denominator very unlikely. |
| Whether the SBS6 point estimate lying near its own bootstrap upper bound (0.428 in [0.116, 0.454]) is expected | A property of `fit_to_signatures_bootstrapped` under strict selection | Inspect the 200 replicate exposures in `wes_signature_refit_bootstrap.csv`. The manuscript states the interval without comment; that is defensible, but a referee may ask. |
| Whether the exome FASTQ files exist | Not knowable from the archive; see **Mj3** | The PI. |
| Whether A6NIZ1 / A6NNZ2 / Q6ZSR9 are the three no-symbol accessions | Recorded in `fix_report_rna.md` and the integration logs, not in a deposited CSV | The archived `protein_relative_abundance.xlsx`. |
| The pre-correction fgsea p-values (0.0498 / 0.0500) | Not retained in any output. The manuscript correctly does **not** quote them and makes only the qualitative permutation-depth point | Nothing; the run was not archived. `WRITING_NOTES.md` §5.1 is right that the sentence can be dropped without consequence. |
| Fig. 5D's "half-row below" layout description | Judgement call from the render at the available resolution | Open `fig5.pdf` at high zoom. |
| Whether the ↓ glyph renders on all three loss markers in Fig. 3E | `35_fig3_biology.R:292` adds it to every `direction == "down"` marker and all three qualify, so it should; I could not resolve `SMARCA4`'s arrow in the PNG | Open `fig3.pdf`. |
