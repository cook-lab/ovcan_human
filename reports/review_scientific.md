# Referee report — "A uniform multi-omic resource of 42 human ovarian cancer cell-line models spanning common and rare histotypes"

**Manuscript type:** Data Descriptor / Resource
**Reviewer:** External referee (no prior knowledge of the authors or their work)
**Materials reviewed:** Manuscript text (converted from `docs/manuscript/OvCAN_data_descriptor_draft.docx`) and all 14 figures (`fig1`–`fig6`, `figs1`–`figs8`). I did not read the analysis code.

**Disclosure of table checks.** Where a number in the text looked internally inconsistent, I opened the corresponding summary table in `output/` to determine which value is correct. I did this for exactly six files: `prot_block_missingness.csv`, `wes_cnv_arm_freq_patient.csv`, `rna_qc_metrics.csv`, `rna_silhouette.csv` / `prot_silhouette.csv`, `rna_variancepartition.csv` / `prot_variancepartition.csv` / `rna_pc_confounder_joint.csv` / `prot_pc_confounder.csv`, and `integ_rnaprot_cor.csv` / `prot_bridge_cor.csv` / `supplement_per_line.csv`. Every such check is flagged in situ below. I did not otherwise inspect the analysis outputs.

---

## 1. Recommendation

**Major revision.**

The underlying resource is genuinely useful and the authors' instincts about its weaknesses are unusually good — patient collapse for genomic frequencies, tiered somatic confidence, the TMT compression ceiling, and explicit de-authentication flags are all handled better than in most cell-line resource papers I review. But the manuscript is not yet a submittable Data Descriptor, for three reasons that are independent of one another. First, the deliverable is absent: the abstract states that data "are deposited in public repositories", while every accession in Table 1 is a placeholder, the raw mass-spectrometry files are not yet in the authors' possession, and the WES raw reads may not exist at all. Second, the methods for the proteomics record are not reproducible even in principle — no search engine, database, FDR, or quantification scheme is given. Third, the central validation claim ("*Subtype separation is driven by biology, not source site*") is stated more strongly than the authors' own analysis supports, and in two places the manuscript reports the one summary statistic of several available that most favours the desired conclusion. None of these is fatal to the resource; all require substantive work rather than rewording.

---

## 2. Summary of the manuscript

The authors assemble 42 pre-existing ovarian cancer cell lines from three Canadian centres (CHUM/Mes-Masson, BC Cancer/OVCARE, OHRI) representing 34 patients and seven histotypes, and generate three uniformly processed omic layers on them: bulk RNA-seq in a single sequencing run (31 lines, kallisto/Ensembl 105), TMT 11-plex proteomics in five sets with a pooled internal standard in channel 1 and a technical bridge in channel 10 (31 lines), and tumour-only whole-exome sequencing via nf-core/Sarek with Mutect2 and CNVkit (23 lines CNV, 22 lines mutation calls); 13 lines have all three. They re-curate histotype annotations, map subline/patient family structure so that genomic frequencies are reported per patient rather than per line, and re-filter an archived MAF set that had ignored the Mutect2 `FILTER` field. Technical validation proceeds through data quality (alignment rates, bridge reproducibility, structural per-set protein missingness), recovery of expected histotype biology (PCA and silhouettes, a 22-marker panel scored by Cohen's *d* and AUC, GO/gene-set recovery, RNA–protein concordance benchmarked against CPTAC/CCLE/ProCan), and identity checks against DepMap (5 overlapping lines) and Cellosaurus STR records (30 of 42 lines). Usage Notes provide a consolidated per-line table, an ADC-target expression atlas, a candidate MMR-deficient clear-cell model (TOV21G), a descriptive within-HGSC stratification, and five reuse caveats.

---

## 3. Strengths

These are real and I want them on the record before the criticism.

- **Patient-level collapse is done, and the consequences are quantified.** The CDK12 example — one family's single frameshift inflating a per-line frequency fourfold — is the clearest demonstration of cell-line pseudoreplication I have seen in a resource paper, and the observation that arm-level CNV frequencies are near-invariant to collapse while point-mutation frequencies are not is a genuinely useful methodological point (Technical Validation §3, Fig. 4C–D).
- **The archived-MAF bug is disclosed and fixed rather than buried.** "Variant filtering corrects a root-cause error in the archived analysis, which ignored the Mutect2 `FILTER` field" (Methods) with the OV2295 waterfall in Fig. 4A is exactly the right way to handle this.
- **Retiring `germline_like_vaf` as a somatic filter** because "in near-pure lines it flags bona fide somatic hotspots (for example TP53 R175H and Q192\*) whose variant-allele fraction approaches 1 through loss of heterozygosity" is a subtle and correct call.
- **Negative and null results are reported.** "Defensible somatic BRCA1/2 is zero"; "Genomic HRD is not computable"; "no STR profiling or mycoplasma testing was performed in-house"; "the two endometrioid lines do not co-cluster". Most authors would have omitted these.
- **The ConsensusOV caveat is correctly reasoned.** Arguing that 7/15 lines being called mesenchymal or immunoreactive is itself evidence the labels are invalid for pure cultures is the right inference, and rare.
- **The marker effect-size analysis (Fig. S7) is done properly** and is internally consistent to the digit: I verified all six quoted *d* values, the 16/22 pass count, the 6 failures, and the 8 markers at AUC ≥ 0.80 directly against the figure.
- **Structural TMT missingness is characterised honestly**, including the warning that set-conditional presence/absence "would be an artifact".

---

## 4. Major concerns

### M1. The abstract asserts a deposition that has not happened, and the manuscript still contains unresolved drafting placeholders

Abstract: "**Raw and processed data are deposited in public repositories.**" This is contradicted by the body. Data Records: "All raw and processed data **are being deposited** ... accessions will be finalised before publication." Every one of the eight rows of Table 1 reads `[PLACEHOLDER]`. Methods: "raw mass-spectrometry files are required for ProteomeXchange/PRIDE deposition and **are being obtained from the Morin laboratory**". Data Records: "the analysis notes indicate **tumour BAMs are not archived; confirm that FASTQs exist** for SRA deposit, otherwise the WES record is processed-only."

For a Data Descriptor the deposition *is* the contribution, so this is not a wording problem. Three distinct issues:

(a) The abstract's tense must be corrected, and it cannot be corrected to something acceptable — a Data Descriptor with no accessions is not reviewable at the level of "can a reuser get the data".

(b) If the WES FASTQs do not exist, then one of the three advertised omic layers is deposited as derived calls only, from a pipeline whose capture kit is unknown (see M2), and whose CNV reference is "unmatched and their capture-kit concordance is unconfirmed". That is a materially weaker record than the abstract's "whole-exome sequencing (23 lines)" implies, and the abstract must say so.

(c) The following editorial placeholders remain in the submitted text and must all be resolved: `[CHECK: read configuration and library-prep kit to be confirmed for Methods]` (Methods, RNA); `[NOTE: raw mass-spectrometry files ...]` (Methods, proteomics); `[CHECK: capture kit identity to be stated in Methods]` (Methods, WES); `[NOTE: three deposition-completeness items gate what Data Records can claim and require PI resolution ...]` (Data Records); `[NOTE: on-target and coverage metrics can be reported only if ...]` (Technical Validation §1); `[NOTE: genomic HRD is not computable ...]` (§4); `[NOTE: the renv lockfile will be finalised last ...]` (Code Availability). Two of these ("read configuration and library-prep kit", "capture kit identity") are not editorial notes at all — they are missing mandatory methods.

**Fix:** obtain the accessions and the raw MS; determine definitively whether WES FASTQs exist; state the library-prep kit, read length/configuration, and exome capture kit; delete all bracketed notes and convert the substantive ones (HRD not computable; WES coverage metrics not computable) into plain declarative limitation sentences in Usage Notes.

### M2. The proteomics record is not reproducible as described

The entire description of how 8,430 protein quantities were produced is: "isobaric tandem-mass-tag (TMT) labelling across five 11-plex sets", channel 1 = PIS, channel 10 = bridge, and "The search results comprise 8,430 proteins across 146,830 peptides."

Missing, and all mandatory: instrument and acquisition mode (MS2 vs MS3/SPS-MS3 — this determines the magnitude of the ratio compression the paper's interpretation rests on); LC and fractionation scheme; search engine and version; sequence database and version; enzyme, missed cleavages, fixed/variable modifications; precursor and fragment tolerances; PSM/peptide/protein FDR thresholds and the method used; reporter-ion isolation-interference or co-isolation filter; whether quantification is at MS2 or MS3 and how peptide reporter intensities were rolled up to protein; whether a common-peptide or total-intensity normalisation preceded PIS ratioing; how "log2 abundance" in Fig. 6 is defined. The phrase "scripts `05`/`06` apply no further inter-set centring" describes the authors' downstream handling but not the upstream search.

Similar but smaller gaps elsewhere: nf-core/Sarek version and the GATK/Mutect2 version; CNVkit version and the bin/segmentation parameters; the log2 threshold used to call an arm gained or lost in Fig. 4D (never stated anywhere, yet six arm frequencies are quoted in §3); MutationalPatterns version and, critically, **which COSMIC SBS release** the cosine similarities are computed against; ConsensusOV, PROGENy, and MSigDB Hallmark versions; the gene-set collections used for the GO recovery. "R v4.5.2 with a fixed seed (1234)" and a promised `renv` lockfile do not substitute for these.

**Fix:** a complete proteomics methods paragraph (or an MIAPE-style supplementary table), plus versions and key parameters for every tool named. State the arm-level gain/loss threshold.

### M3. "Subtype separation is driven by biology, not source site" is over-claimed, and the supporting statistics are selectively reported

This is the load-bearing validation claim of §2 and it is stated causally ("driven by"). Three specific problems.

**(a) The genome-wide claim reverses if a different summary statistic is used.** Text: "Genome-wide variance decomposition is concordant: subtype ≥ site in RNA (median 5.9% vs 3.5%)". I checked `output/rna_variancepartition.csv`. The medians are as stated (5.95% vs 3.53%), but the same table gives **mean subtype 14.9% vs mean source_site 14.5%** — indistinguishable — and the interquartile ranges are 0–24.8% (subtype) versus 0–25.2% (site), i.e. the site upper quartile *exceeds* the subtype upper quartile. Fig. 3D plots the IQR bars, and they overlap almost completely; a reader looking only at the figure cannot conclude subtype > site. Reporting only the median, in a right-skewed per-gene distribution where the mean tells the opposite story, is not defensible. The protein comparison is more robust (subtype median 8.6%/mean 16.7% vs site median ~0%/mean 8.2%) and should be presented as the stronger of the two.

**(b) The commonality decomposition is presented in raw R² when the authors computed adjusted R².** Text: "source site adds ≤0.2% of PC1 variance beyond histotype (unique-subtype 42%, unique-site 0.2%, shared 31%)". These match `output/rna_pc_confounder_joint.csv` exactly. But the same file records `adjr2_joint = 0.657` against `r2_joint = 0.737` for PC1 — an 8-point shrinkage, because a 6-level subtype factor plus a 3-level site factor is ~8 parameters fitted to ~30 observations, and four of the six histotypes have n ≤ 3. Fitting a group with n = 2 costs one parameter and buys near-perfect fit; that is what inflates unique-subtype. Report the adjusted values, and add a permutation null (shuffle subtype labels within site, and vice versa) so the reader can see how much of the 42% is attainable by chance at this n.

**(c) The metric that shows the largest batch effects is the one not reported.** The manuscript reports the protein batch effect only as a per-protein median: "subtype ≫ site ≈ set in protein (8.6% vs 0.0% vs 0.9%)". `output/prot_pc_confounder.csv`, which is never referenced in the manuscript, gives for the leading protein PCs: **site adjusted R² = 0.315 on PC1**, and **TMT plex adjusted R² = 0.215 on PC2 and 0.229 on PC5**. So on the leading protein components, plex explains ~20% and site ~30%. Both statements can be true simultaneously (a batch effect concentrated in a few high-variance proteins has a low per-protein median and a high PC loading), but reporting only the reassuring one is selective. Fig. S1B is the figure that would show this — and it is never cited (see M8).

**(d) An uncited figure argues against the claim.** Fig. S3 (never cited) is a line × line RNA correlation heatmap with a dendrogram whose **top-level split is essentially CHUM (TOV/OV) versus BC Cancer (VOA)**. That is the single most direct visual assessment of site structure in the dataset and it is neither shown to the reader nor discussed. It should be cited and confronted.

**Fix:** retitle the paragraph to something the analysis supports — e.g. "Source site adds little to the leading components beyond histotype, but the two are strongly collinear and cannot be fully separated in this design." Report adjusted R², the mean as well as the median in Fig. 3D, the protein PC confounder table, and a permutation null. Cite Fig. S3 and address the CHUM/VOA split.

### M4. Pseudoreplication is corrected for the genomics and left uncorrected for the expression and protein analyses

The manuscript's own standard is explicit: "Genomic frequencies throughout are reported per patient" and "Patient-level tabulation also corrects pseudoreplication that inflates per-line frequencies." That standard is applied to the marker effect sizes ("across **28 patient-level RNA representatives**") but not to the other expression analyses, and the manuscript does not say so.

I checked `output/rna_silhouette.csv`: the silhouettes quoted in §2 are computed on **31 lines** (HGS n = 15, CC 7, EC 2, MC 3, MMMT 2, SCCOHT 2), not 28 patients. The 15 HGSC "samples" come from 12 patients (families 1369, 2295, and 3133 each contribute two RNA lines). Same-patient sublines are by construction each other's nearest neighbours — Fig. S3 shows TOV3133G/OV3133-R as the darkest off-diagonal pair in the matrix — so pooling them inflates the HGSC silhouette. The PCA in Fig. 3A/B and the variance decomposition in Fig. 3D are likewise per line. Two units of analysis are used within one paragraph without flagging the switch.

**More seriously, the Fig. 3D model is not identifiable as specified.** It fits `patient` as a term on 31 RNA samples drawn from 28 patients, i.e. **exactly three within-patient replicate pairs**. The patient variance component (27.3%) and the residual (32.8%) are therefore estimated from 3 degrees of freedom, and the statement "with line/patient identity the largest structured term, as expected for distinct cell lines" is close to tautological — with 28 levels on 31 observations, `patient` absorbs essentially all reproducible variance regardless of the biology. The protein model is worse in the same way (30 samples, 27 patients). This matters because Fig. 3D is one of the two pillars of the M3 claim.

**Fix:** recompute the PCA, silhouettes and variance decomposition on the 28 patient representatives and report those numbers; or, if the per-line version is retained, report both and state which is which. Drop `patient` from the variance model or replace it with `family` restricted to the multi-line families, and say plainly that within-patient replication is n = 3 pairs.

### M5. A stated null result is contradicted by the authors' own QC table

Technical Validation §1: "A modest source-site difference in pseudoalignment rate (Mes-Masson 92.2% vs Huntsman 88.1%) **does not translate into a difference in the number of genes detected, so it does not affect downstream detection.**" No test, no effect size, and no uncertainty is given for this null.

Fig. 2A shows the opposite by eye: nearly every Huntsman triangle sits above nearly every Mes-Masson circle on the *y* axis. I checked `output/rna_qc_metrics.csv`:

| site | n | median pseudoalignment | median genes detected | median library size |
|---|---|---|---|---|
| Mes-Masson | 19 | 92.2% | **19,817** | 54.4 M |
| Huntsman | 11 | 88.1% | **20,765** | 68.3 M |
| Huntsman/Vanderhyden | 1 | 87.9% | 20,119 | 43.8 M |

Mann–Whitney on genes detected, Mes-Masson vs Huntsman: **p = 0.0098** (my calculation on the authors' table). Across all 31 lines, pseudoalignment rate and genes detected are **negatively** correlated (Pearson r = −0.63, p < 0.001). And library size differs between sites by ~14 M reads, which is a second site-confounded technical variable that appears nowhere in the manuscript and is the likely mechanism.

The direction is benign for the resource (the lower-alignment site detects *more* genes, consistent with deeper libraries), so this is fixable by rewriting rather than by new data. But as written the sentence asserts a null that the authors' data reject, and it removes from the reader a real site-confounded covariate.

**Fix:** replace with the effect sizes and the test; add median library size per site; note that pseudoalignment rate is inversely related to gene detection here and that depth, not alignment rate, is the operative difference. Add library depth to the confounder models in Fig. 3D and to the Reuse caveats.

### M6. Three identity claims cite figures that do not contain the stated evidence

(a) §4, "Independent recovery of two published reclassifications": "COV434 ... and BIN67 are consistent with SCCOHT: both show **SMARCA4/SMARCA2 loss** ... **(Fig. 5**; Karnezis et al. 2021)" and "TOV112D ... SMARCA4 truncation (p.L639X), **low protein, silenced SMARCA2 mRNA**". **Fig. 5 contains no SWI/SNF panel.** Its four panels are the TOV21G variant burden, the TOV21G SBS-96 spectrum, TOV21G cosine similarities, and the three mucinous lines. There is no per-line SMARCA4/SMARCA2 RNA or protein anywhere in the 14 figures — Fig. 3E shows subtype *means* only. This is the claim the Background & Summary advertises ("including independent recovery of two published reclassifications") and it currently has no visual support at all.

(b) The most interesting single observation in the paper — "BIN67 is the instructive case: **SMARCA4 mRNA is retained but protein is second-lowest** (post-transcriptional loss), so SMARCA4 IHC/protein is confirmatory and RNA-only authentication would miss it" — has no figure, no quantification beyond "second-lowest", and no uncertainty. For a resource paper this is a directly actionable authentication lesson and it deserves a panel.

(c) §4: "Each self-matches its DepMap namesake at Spearman 0.74–0.88, **ranking first of 67 DepMap ovarian lines and reciprocal-best in both directions**; ... **specificity is the signal (Fig. 4)**." Fig. 4B is a 5 × 5 submatrix. A 5 × 5 panel cannot demonstrate a rank of 1 out of 67, and it cannot show the margin between the self-match and the best non-self match — which is the entire point when the authors themselves say specificity rather than magnitude is the signal. In Fig. 4B the off-diagonal TOV112D column is visibly comparable to the BIN67 diagonal (0.74), which is biologically sensible but visually undercuts the panel.

**Fix:** add a SWI/SNF panel (per-line SMARCA4/SMARCA2 RNA and protein, with the WES truncation annotated) and cite it for the reclassification claims. Replace or supplement Fig. 4B with a self-match-versus-best-other rank plot across all 67 DepMap ovarian lines.

### M7. The CDK12 sentence contains a numerical contradiction

§3: "the apparent CDK12 enrichment (6/17 lines, 35%) collapses to **2–3 of 11 patients (18–27%), matching the ~3% TCGA-HGSC rate**".

18–27% does not match 3%; it is six- to nine-fold higher. The rhetorical point — that pseudoreplication inflated the estimate — survives, but the conclusion drawn from it does not. With Tier 3 excluded from headline frequencies as stated in Methods, the figure is 2/11 = 18%, still ~6× TCGA. Note also that all CDK12 calls in Fig. 4C are Tier 2 (salmon) or Tier 3 (light), i.e. "plausible somatic" or "cannot exclude germline"; a frequency estimate built entirely from non-Tier-1 calls should probably not be quoted as a frequency at all.

**Fix:** state the corrected number honestly ("collapses to 2 of 11 patients (18%), still above the ~3% TCGA-HGSC rate and based only on Tier 2 calls, so it should not be treated as a somatic frequency estimate"), or drop the TCGA comparison.

### M8. Four of eight supplementary figures are never cited, and one of them contains a confounder that is absent from the entire manuscript

**Fig. S1, S2, S3 and S5 are cited nowhere in the text.** (S4, S6, S7, S8 are cited.) Table S3 is also never cited, though S1, S2, S4, S5 and S6 are. Separately, the cross-reference at §1 — "on-target and coverage metrics can be reported only if the raw reads or alignment metrics are recoverable ... (**Fig. S8** is contingent on this)" — points at the within-HGSC heterogeneity figure, which has nothing to do with WES coverage. Either a WES QC supplementary figure is missing or the numbering is wrong.

The most consequential of these is **Fig. S2**, which reveals two facts that appear nowhere in the manuscript:

- **Passage number explains ~7.8% of RNA PC1 raw and ~14.7% conditional on site** (Fig. S2A). That is larger than the genome-wide subtype variance component reported in Fig. 3D (5.9%). Passage is not in the Fig. 3D model, is not in the PC confounder model, is not discussed in §2, and is not in the Reuse caveats.
- **Passage is not matched between assays for the same line** (Fig. S2B): TOV112D is ~63 for RNA and ~83 for WES; OV3331 is ~71 for RNA and ~54 for WES. So the "matched multi-omic" layers are, for some lines, ~20 passages apart. Given the size of the passage effect in panel A, this is a genuine limitation for every cross-assay analysis in the paper — the RNA–protein concordance, the RNA-versus-WES identity checks, and the ADC atlas.

**Fix:** cite all eight supplementary figures or remove the uncited ones. Add passage to the confounder models and to the Reuse caveats, and state the RNA-versus-WES passage discordance explicitly with its range. Fix the Fig. S8 cross-reference.

### M9. The protein count is wrong, and 71 proteins in the matrix are quantified nowhere

Methods: "The search results comprise **8,430 proteins** ... Of 8,430 proteins, 6,856 (81.3%) are quantified across all five sets, and 1,573 (18.7%) are absent from at least one whole set."

6,856 + 1,573 = **8,429**, not 8,430. I checked: `output/prot_block_missingness.csv` has 8,429 rows, and its `present_n_plex` column tallies 5:6856, 4:512, 3:366, 2:360, 1:264, **0:71**. Fig. 2C reproduces these exactly and sums to 8,429. So the correct total is 8,429 and "8,430" is an error propagated through Methods and §1. The percentages (81.3/18.7) are computed against the wrong denominator but round to the same values.

Separately: **71 proteins are present in 0 of 5 plexes** — they are in the matrix with no quantification anywhere. This is not mentioned in the text or in Fig. S4 (which shows only the 31 non-empty patterns), yet these 71 are silently inside the "1,573 absent from at least one whole set" count. A reuser needs to know what these rows are and whether they are in the deposited matrix.

**A related clarity problem:** the manuscript uses at least four different protein feature counts without saying which analysis uses which — 8,430/8,429 (search output), 7,734 ("A presence filter requiring quantification in at least 50% of lines retains 7,734 proteins for analysis"), 7,893 (the actual per-gene concordance set, see M12), and 6,856 (Fig. 3D panel title and the variance decomposition). Each analysis should name its feature set.

### M10. Bridge "reproducibility" measures the wrong quantity, and one bridge sample is not part of the resource

§1: "Bridge technical-replicate reproducibility is high: per-link Pearson correlation is 0.991–0.994 across the four daisy-chain links, with a median coefficient of variation of 5.3%."

(a) A Pearson correlation computed on log2 abundances spanning ~12 log2 units (Fig. 2D *x* axis runs 5 to ~18) is dominated by dynamic range, not by reproducibility; r ≈ 0.99 is close to unavoidable for any two measurements of the same proteome and cannot discriminate good from mediocre. This is the classic correlation-versus-agreement confusion. The median CV of 5.3% is the informative statistic and should lead. Add Bland–Altman (MA) plots and the full CV distribution stratified by abundance decile, since compression and noise are both abundance-dependent — which is the paper's own thesis in §2 and Fig. 6.

(b) Fig. 2D labels the third link "plex 3→4 · **VOA3993\***". VOA3993 is not one of the 42 lines in Fig. 1B, and the asterisk is never explained — there is no footnote in the figure and no caption anywhere. I checked `output/prot_bridge_cor.csv`: the four links use four *different* lines (VOA10816, TOV1369, VOA3993, TOV3133G) and VOA3993 carries `external = TRUE`. So (i) the Methods description "channel 10 carries **a** separate technical-replicate 'bridge' sample" is singular and misleading — it is a different sample per link, which makes this a chain of four independent one-off technical replicates rather than a single reference threaded through all five sets, and (ii) one of the four links is validated using material outside the deposited resource. Both must be stated in Methods. The same table gives n = 7,276–7,436 proteins per link, which should be reported with the correlations.

### M11. Fig. 6's protein panel cannot support the claim it is cited for, and the compression factor is inflated by a floor artefact

Usage Notes: "Known associations recover **in both RNA and protein**: mesothelin is highest in HGSC, and HER2/ERBB2 is higher in clear cell and mucinous (Fig. 6)."

The protein heatmap in Fig. 6 is plotted on an absolute log2-abundance scale spanning **9 to 16**, while the paper's own Fig. 2B establishes that the median cross-line protein IQR is **0.34 log2 units**. The consequence is visible: the protein panel is almost uniformly mid-orange, and the MSLN protein row in particular shows no discernible variation. A reader cannot verify "mesothelin is highest in HGSC" in protein from this panel; the claim is unsupported by the cited figure. Row-scaling (z per protein) or a per-protein colour scale is required. The same applies to B7-H3 and B7-H4.

Two further problems in the same paragraph:

(a) "**FOLR1 is strongly bimodal** within HGSC (RNA range 0.06–9.5)". A range is not evidence of bimodality. Show the distribution (a per-line dot plot or histogram) and, if the word "bimodal" is retained, a dip test or a two-component mixture fit. As drawn, the two mucinous lines (VOA8762, VOA8771) look comparable in FOLR1 to TOV3133G, which makes "the highest FRα expressers panel-wide are HGSC lines" hard to verify from the figure.

(b) "TMT ratio compression narrows protein target ranges 3–5× relative to RNA (for example **MSLN protein range 1.8 vs RNA 9.3**)". These are not comparable quantities. The RNA scale is log2(TPM+1) with a hard floor at 0 for undetected transcripts, so a gene with any zero value automatically has a large "range"; the protein scale has an effective detection floor near 9 and no zeros. Much of the apparent 5× is the RNA zero-floor, not TMT compression. Use a floor-insensitive statistic on both (IQR or SD ratio, both of which exist in `output/prot_dynamic_range.csv`) and quote that instead. The same criticism applies to the Fig. 2B RNA-versus-protein IQR comparison, which additionally compares a VST/log-TPM scale against a log ratio-to-reference scale — some of the difference is normalisation, not compression. This weakens a mechanistic claim that is otherwise well argued.

(c) None of the Fig. 6 claims carries an effect size, a test, or an *n*. "Highest in HGSC" and "higher in clear cell and mucinous" are group statements about groups of n = 2–6. Either quantify them (with the caveat that n = 2 groups cannot support them) or state them as descriptive observations only.

### M12. No figure legends, no supplementary table legends, and a reported *n* that does not match the analysis

**(a) There are no figure captions anywhere.** I checked the source `.docx` directly to be sure this was not a conversion artefact: it is not. Not one of the 14 figures has a legend. This has substantive consequences throughout — I could not determine, for example, whether the *z*-scores in Fig. 5D are computed across the 3 mucinous lines or all 30 (which decides whether that panel is circular), what "98% of genes compressed" in Fig. 2B means operationally, what "artefact collapse ATM 82%→9%, ATR 77%→5%" in Fig. 4A is a percentage *of*, what the asterisk in Fig. 2D denotes, which modality the Fig. 3E *z*-scores use, whether Fig. 3A is 28 or 31 samples, or what the two indistinguishable grey swatches in the Fig. 4C "Variant class" legend encode. Supplementary table legends are also absent.

**(b)** §2: "On 30 lines and **8,212 shared genes**, per-gene median Spearman correlation is 0.40". I checked `output/integ_rnaprot_cor.csv`: there are **7,893** per-gene correlations (and 30 per-line values). If 8,212 is the count of shared symbols before some downstream filter, say so; as written the *n* attached to the reported median is wrong by 319 genes.

**(c)** The same file shows the per-gene *n* is **not** 30 for all genes — it ranges from 10 to 30, with 622 genes at n < 20 and 353 at n < 15, because of the structural missingness. A Spearman coefficient on n = 10 is nearly uninformative, and pooling those with n = 30 into a single median is inappropriate without disclosure. To the authors' credit the conclusion is robust: restricting to the 6,685 complete-case genes gives a median of 0.408 versus 0.397 overall (my calculation). Report the complete-case value as the headline, disclose the *n* distribution, and note that **10.3% of genes have negative RNA–protein correlation** — a fact worth a sentence.

### M13. The mucinous de-authentication is consequential for reusers but rests on a very thin panel that is presented circularly

§4: "VOA8762 and VOA8771 show an intestinal/GI-leaning profile (CK7-low, PAX8-low, CDX2-high) ... **their ovarian origin is not supported by expression** (Fig. 5)" → "**Genuine ovarian mucinous coverage may therefore be n=1** (TOV2414) pending STR and IHC."

This is the kind of statement that will determine whether anyone uses these two lines again, so the evidence bar should be high. Currently: six genes, one sample per line, no replicates, no uncertainty, no external comparison. Three specific problems.

(a) **Fig. 5D's column headers pre-assert the conclusion** — "TOV2414 **ovarian**", "VOA8762 **GI-leaning**", "VOA8771 **GI-leaning**". Labelling the columns with the interpretation under test is circular presentation. Use neutral line names and put the interpretation in the caption.

(b) **The z-scoring reference set is unstated.** If *z* is computed across only these three columns, then one line being high and two low is guaranteed by construction and the panel is uninformative. This must be specified (and if it is the three-column version, recomputed against all 30 lines).

(c) **The intestinal call is weaker than described for VOA8762.** The text cites CDX2-high plus SATB2 status, but in Fig. 5D VOA8762's SATB2 is essentially at the reference mean (near-white); only VOA8771 shows any SATB2 elevation, and it is modest (~+1.5). So for VOA8762 the "intestinal" assignment rests on CDX2 alone. Note also that Fig. 3E shows MUC2 highest in the **HGSC** column, and the text attributes this to 2D silencing — the same argument (culture-induced drift in mucinous differentiation markers) is a live alternative explanation for the VOA lines' profile that is not considered.

The authors do hedge appropriately ("pending STR and IHC"), which is why this is M13 rather than higher. But the panel should be strengthened: more markers (e.g. the Cheasley/Meagher panels the authors already cite), an external anchor (TCGA/CPTAC mucinous ovarian versus colorectal/gastric), and per-line values against the full panel distribution rather than three columns.

### M14. Missing resource metadata that reusers will need, and one histotype with no marker validation at all

For a cell-line resource these are not optional:

- **Culture conditions are described nowhere.** No medium, serum concentration, supplements, or incubation conditions for any of the 42 lines. Since expression phenotypes in these lines depend on medium (and since the lines come from three labs with different house protocols — a further site confounder), this is required metadata and should be a column in `samples.csv`.
- **No ethics or consent statement.** The lines are human patient-derived. Even for pre-existing published lines, a statement of the governing approvals under which they were derived and are used is required.
- **Mycoplasma status is absent by admission** ("No STR profiling or mycoplasma testing was performed in-house"). This is honest, but the abstract's "the data are of high quality" should be qualified accordingly, and a mycoplasma-contamination expression check (e.g. unmapped-read screening, which the RNA-seq permits) is cheap and would partially close the gap.
- **There is no endometrioid marker set.** Fig. S7 covers HGS, CC, MC, MMMT and SCCOHT — five intended subtypes, 22 markers, no EC. So the one histotype whose membership the paper actively revises is the one with no marker-level validation. This matters because the paper concludes "VOA4395 becomes the sole endometrioid line once TOV112D is reclassified" — a single line, whose endometrioid identity is nowhere tested. Add EC markers (ARID1A, ESR1, PGR, vimentin per Hollis et al. 2020, which is already cited) or state explicitly that endometrioid identity is unvalidated in this resource.

### M15. "Corroborated orthogonally by PROGENy" is not orthogonal, and Fig. S8 colours points by their own axes

Usage Notes: "the 15 HGSC lines resolve into three descriptive pathway-activity strata (inflammatory/NF-κB-EMT; low-signaling; hypoxic-glycolytic), **corroborated orthogonally by PROGENy**".

PROGENy scores are computed from the same RNA matrix, and the pathways shown in Fig. S8B include TNFa, NFkB, MAPK and Hypoxia — i.e. the same biology as the Hallmark axes that *define* the strata. This is not orthogonal corroboration; it is the same signal via a second gene-set resource. Furthermore Fig. S8A plots the strata against the two scores used to construct them and colours the points by the resulting label, so the panel cannot fail.

Looking at Fig. S8B, the corroboration is also weaker than claimed: of the five lines in the "Hypoxic-glycolytic" stratum, only OV4453 and OV3331 show elevated PROGENy Hypoxia; TOV1369, TOV3041G and TOV3291G are near or below zero, while OV1946 (Inflammatory stratum) is elevated.

Relatedly, "with same-patient families **mostly co-clustering**" rests on **three** assessable families in the RNA set: family 2295 and family 3133 co-cluster; family 1369 does not (OV1369-R2 = Inflammatory, TOV1369 = Hypoxic-glycolytic in Fig. S6). "2 of 3" should be stated rather than "mostly". Fig. S3 makes the same point independently: OV1369-R2 and TOV1369 are not adjacent in the correlation dendrogram.

**Fix:** replace "corroborated orthogonally" with "consistent with an independent pathway-scoring method applied to the same data"; drop or heavily caveat the corroboration claim; use a genuinely independent axis (protein, or copy number) if orthogonality is wanted; state 2 of 3 families.

---

## 5. Minor concerns

1. **Novelty claim rests on an uncitable source.** "no ovarian-cell-line multi-omic *Scientific Data* Data Descriptor was found (**literature review, this project**)". An internal unpublished search cannot support a negative-novelty claim in a published paper. Either describe the search (databases, dates, query) in Methods so it is reproducible, or soften to "we are not aware of".
2. **The abstract over-states the identity corroboration.** "molecularly corroborate line identity against DepMap and Cellosaurus" — DepMap covers 5 of 42 lines (12%), and Cellosaurus provides STR records deposited by *other* laboratories for 30 of 42, with no in-house verification. The body is candid about this ("no BC Cancer VOA line exists publicly"; "The remaining 12 models ... have no public STR record"); the abstract should be too.
3. **The abstract's "high quality" is not established for the WES arm.** By the authors' own note, on-target and coverage metrics may not be computable at all. The abstract should not make a blanket data-quality claim covering an assay whose QC metrics are unavailable.
4. **WES *n* is quoted as 23 in the abstract** but is 23 for CNV and 22 for mutations (Table 1 and Fig. 4C both show 22). State both.
5. **Site nomenclature drifts.** Methods name "three centres: CHUM ..., BC Cancer / OVCARE ..., and OHRI (Vanderhyden BIN67)". Fig. 2A uses "Mes-Masson / Huntsman / Huntsman/Vanderhyden"; Fig. 3B collapses to two ("Mes-Masson / Huntsman"), silently folding BIN67 into Huntsman. I checked `supplement_per_line.csv`, which uses `Mes-Masson` (29), `Huntsman` (12), `Huntsman/Vanderhyden` (1) — so the metadata label for the OHRI line is "Huntsman/Vanderhyden", not OHRI. Pick one scheme (institutional labels are preferable to investigator surnames for a public resource), use it in text, figures and `samples.csv`, and explain what "Huntsman/Vanderhyden" means. Also note that Fig. 3B's two-level site variable is not the three-level variable used for the R² in `rna_pc_confounder_joint.csv` (which stores both `r2_site` and `r2_site_2lvl`, 0.313 vs 0.227 on PC1); the text should say which was used for the quoted numbers.
6. **Histotype abbreviations drift between text and figures.** Text uses HGSC, MMMT/carcinosarcoma, low-grade serous; figures use HGS, MMMT, LGS. For a resource whose users will join on these strings, define the abbreviations once, note that the figures and `samples.csv` use the short forms, and be consistent.
7. **Fig. 3D's panel title says 22,544 genes; the variance decomposition ran on 22,542** (`rna_variancepartition.csv`, `n_features`). Trivial, but reconcile.
8. **The top-2-of-6 marker rule has no stated null.** With six subtype columns, a marker lands in the top 2 by chance one time in three, so 16/22 should be quoted against that null (a binomial test gives p ≈ 1×10⁻⁴). Also state why the denominator is six columns when only five subtypes have intended markers, and why LGS is absent (no RNA for TOV81D — worth one clause).
9. **No confidence intervals on any effect size.** Cohen's *d* and AUC in Fig. S7, the silhouettes in §2, all six arm-level frequencies in §3, and the Spearman correlations in §4 are point estimates. Bootstrap CIs are cheap and would let the reader see immediately that SMARCA2 *d* = −2.8 (n = 2 vs 26) and SCCOHT silhouette 0.82 (n = 2) are not comparable to KRT20 *d* = 3.1. The authors already caveat the n = 2 silhouettes verbally; do the same quantitatively and extend it to the effect sizes.
10. **Protein silhouettes are computed but never reported.** I checked `output/prot_silhouette.csv`: HGS 0.18, CC **−0.003**, EC 0.13, MC 0.14, MMMT 0.41, SCCOHT **0.028**. The text's "and, more weakly, in protein" substantially understates this — clear cell is *negative* and SCCOHT is ~zero in protein, and the quoted silhouettes (HGSC 0.16, CC 0.12, MC 0.15, MMMT 0.74, SCCOHT 0.82) are the RNA values, which the text never says. Report both sets and label the modality.
11. **The EC silhouette is available and should be quoted.** "the two endometrioid lines do not co-cluster" is supported by RNA silhouette = −0.014; give the number.
12. **Mutational signature analysis lacks uncertainty and discrimination.** "SBS6 cosine 0.88, plus SBS44/SBS15/SBS20 ... and no POLE signal". Cosine similarity to a single reference is not a fit; SBS1/5/6/15/44 are mutually similar (Fig. 5C shows the clock group at 0.65, not far below 0.88), so 0.88 does not discriminate strongly. Use a signature refit (e.g. NNLS or `fit_to_signatures_bootstrapped`) with bootstrap intervals, report the full cosine profile, and state the COSMIC version. Also, "no POLE signal" inferred from a cosine of 0.32 is a null asserted from a low similarity; frame it as "no evidence of" rather than "no". Note also that Fig. 5C lists SBS6/44/**15** while the text lists SBS44/SBS15/**SBS20** — reconcile. Finally, since ~206 residual coding variants per line are germline-derived, state how many of TOV21G's 1,416 variants were used for the spectrum and whether germline contamination could shape it.
13. **"Enriched for indels" is unquantified.** Fig. 5A encodes indel fraction as colour but the text gives no value for TOV21G or the panel. Give the fraction and the comparison.
14. **Fig. 4A's "artefact collapse ATM 82%→9%, ATR 77%→5%"** is not explained in the text — percentages of what? These look like informative numbers; define them.
15. **"26 expression-consistent, 20 genomics-consistent"** (§4) is given without denominators. Out of how many assessable lines each? And what were the consistency criteria?
16. **The clear-cell cross-site control is under-specified and under-powered.** "within clear cell (2 CHUM + 5 BC Cancer), source lab explains only 4–6% of the leading components". n = 2 versus n = 5 gives essentially no power to detect a site effect, so "explains only 4–6%" is not reassurance; report a permutation p or a CI so the reader can see the power. Also state that n = 7 is the RNA subset of the 8 clear-cell lines (VOA14993 has no RNA), since Fig. 1C shows 8.
17. **"Its strongest gap-filling coverage is for ... SCCOHT ... and carcinosarcoma"** — the resource contributes 2 lines of each. That is a real contribution for SCCOHT, but "two carcinosarcoma models (VOA5217, VOA5436)" should be read alongside the paper's own flag that "VOA5436 (annotated MMMT) expresses a strong clear-cell program". Effective validated carcinosarcoma coverage may be n = 1, exactly as the authors concede for mucinous. Say so.
18. **Fig. 4C's `TP53` row shows 12 patients** whereas §3 says 11 HGSC patients. I verified this is consistent (11 HGSC + TOV112D), but it will read as a contradiction without a caption stating that the right-hand bars count all patients, not HGSC patients.
19. **Fig. 1B appears to mark WES-SNV for both TOV3121EP and TOV3121D**, but Fig. 4C and Fig. 5A show 22 MAF lines including only TOV3121EP, while Fig. S5 shows 23 CNV lines including both. Verify the Fig. 1B SNV column against `samples.csv`.
20. **`variancePartition` substitution.** "an equivalent per-gene `lme4` REML variance decomposition is used" — "equivalent" needs support. `variancePartition` uses a specific parameterisation and reports fractions of total variance; state the exact model formula and how percentages were computed, and ideally show agreement with `variancePartition` on a subset in a compatible R session.
21. **`ignoreTxVersion = TRUE` and 39,568 → 22,544 genes**: state the tximport `countsFromAbundance` setting, since it determines whether the DESeq2 counts are length-corrected.
22. **Table 1 lists "figshare / PRIDE" and "GEO/SRA or figshare"** as alternative repositories for the same record. Commit to one per record; a reuser cannot follow "or".
23. **Table S3 is never cited.** Either cite it or renumber.
24. **The Data Records section never gives column definitions.** "per-matrix column and unit definitions are provided in the repository README" — for a Data Descriptor, a data-dictionary supplementary table is expected in the paper itself, not deferred to a README.

---

## 6. Figure-by-figure assessment

### Fig. 1 — Resource overview
**Readable: yes**, at the size provided; panel B's row labels are small but legible.
- Panel B: the **`Sub.` colour strip has no legend**. It is decodable only by cross-referencing the right-hand text labels (HGS/CC/…) or panel C. Add a subtype legend, or move the subtype key into panel B.
- Panel B: "Assay present" orange is the **same hue as the HGS subtype swatch and the HGS bar in panel C**. Three different variables share one colour. Use a neutral fill (grey/black) for presence.
- Panel B: family colours here (1369 purple, 2295 teal, 3121 orange, 3133 blue, 3291 magenta) match Fig. 4C but **not Fig. S3** (see below).
- Panel B: consider marking the **13 lines with all three assays**, since that subset is quoted in the abstract and is the most valuable slice of the resource.
- Panel C: encoding lines as bar length and patients as an overplotted diamond on the same axis works, but the printed numerals are line counts only. Print both, or label them.
- Panel C: verify the numeral placement — for HGS the "24" sits at the bar end while the diamond is at 16; for the other histotypes the two coincide, which could read as though the numeral labels the diamond.
- I verified panel B (24+8+3+2+2+2+1 = 42 rows) and panel C (patients summing to 34) against the Methods composition. Both are consistent.

### Fig. 2 — Sequencing and proteome quality
**Readable: yes.**
- **Panel A contradicts the text** (see M5): Huntsman triangles sit systematically above Mes-Masson circles on genes detected. Add the per-site medians and the test to the panel or the caption, and add a library-size panel, since depth is the likely driver.
- Panel A: the *y* axis spans ~18.7–21.7k, a narrow window that visually amplifies a ~5% difference. That is acceptable for a scatter, but the caption should give the relative magnitude.
- Panel A: "median 91" annotates the vertical reference line; add an equivalent reference for genes detected, and state *n* = 31 on the panel.
- Panel B: **"98% of genes compressed" is undefined**. Define the criterion.
- Panel B: comparing RNA VST/log-TPM IQR against protein log ratio-to-PIS IQR conflates dynamic range with normalisation scale (see M11b). Note the caveat.
- **Panel C sums to 8,429, not the 8,430 stated in the text** (M9), and the "0 plexes present" bar (71 proteins) is unexplained.
- Panel D: **"VOA3993\*" is not one of the 42 lines and the asterisk is never defined** (M10b). Each of the four links uses a different line, contradicting the Methods' singular "a separate technical-replicate 'bridge' sample".
- Panel D: report *n* proteins per panel (7,276–7,436) and replace or supplement Pearson r with MA plots and CV distributions (M10a).

### Fig. 3 — Recovery of canonical subtype biology
**Readable: yes**, though panel B's annotation text overlaps a data point.
- Panels A/B: **state *n***. I could not tell from the figure whether the PCA is on 31 lines or the 28 patient representatives; the silhouette table implies 31, while the marker analysis in the same section uses 28 (M4).
- Panel B: the annotation "Site adds 0.2% of PC1 variance beyond subtype" **overlaps a data point** at the top of the panel — reposition.
- Panel B: visually, Mes-Masson occupies the lower-left and Huntsman the upper and far-right regions. The 0.2% figure is the *unique* site contribution; the 31% shared component is what a reader is seeing. Annotate all three components (unique-subtype / unique-site / shared) on the panel so the collinearity is visible, not just the reassuring number.
- Panel B collapses three sites to two (M-minor 5).
- Panel A: I could not confirm "the two endometrioid lines do not co-cluster" from this panel — both purple points sit at PC1 ≈ 100–115 with PC2 ≈ −35 to −45, which reads as reasonably close; one of them overlaps the SCCOHT cluster. Annotate the two EC lines by name (the TOV112D-with-SCCOHT placement is itself part of the reclassification story) and quote the EC silhouette (−0.014).
- Panel D: **IQR bars for subtype and source site overlap almost entirely in the RNA panel**, which undercuts the "subtype ≥ site" claim the text draws from it (M3a). Add the means, and state that per-gene percentages need not sum to 100 (they sum to 69.5% RNA / 66.8% protein as plotted).
- Panel D: the `Patient` term is not identifiable at n = 3 within-patient pairs (M4). Either remove it or annotate the replication.
- Panel D: `TMT plex` appears only in the protein panel; note that RNA had no batch term because all libraries were one run.
- Panel E: 23 rows are shown but the text says 22 markers. MKI67 (right-labelled "Prolif", boxed under SCCOHT) is evidently not one of the 22 — Fig. S7 confirms 22 rows without MKI67. State this, and clarify why MKI67's box is in the SCCOHT column when its row group is "Prolif".
- Panel E: **the boxes do not encode expected direction.** SMARCA4/SMARCA2 are boxed in SCCOHT and appear dark blue (low), while every other boxed cell is orange (high). A reader who does not know these are loss markers will read the SCCOHT block as a failure. Add a direction glyph or split the panel.
- Panel E: **no endometrioid marker row group** despite an EC column (M14).
- Panel E: shows subtype **means** only, so it cannot reveal whether a subtype-level signal is driven by one line. Fig. S7 gives effect sizes but also not distributions. Add per-line values (a dot-per-line overlay, or move Fig. S7 to per-line dot plots).
- Panel C: good panel. Add *n* (which should be 7,893, not 8,212 — M12b) and note that 10.3% of genes are negatively correlated.
- Colour: the diverging blue–orange in panel E and the categorical palette in A/B/D are colour-vision-safe. No rainbow, no red–green. Good.

### Fig. 4 — Genomic fidelity
**Readable: mostly.** Panels A, B and D are clear; panel C's top bar is not interpretable (see below).
- **Panel C's "coding cand." bar track has no axis and no scale** — only the TOV21G value (1416) is labelled. A reader cannot read any other line's burden. Add an axis.
- **Panel C's "Variant class" legend shows two grey swatches ("truncating", "missense") that are visually identical to me**, and I could not find any element in the plot that uses them. Either the legend is orphaned or the encoding is invisible at this rendering. Fix or remove.
- Panel C: the right-hand "patients" bars have no axis. Add one, and state in the caption that they count all patients (hence TP53 = 12, not 11 — M-minor 18).
- Panel C: the Tier 2 salmon, the HGS subtype orange and the LGS subtype salmon are close in hue. Recolour the tier scale to a hue not used for subtype.
- Panel C: CDK12 is 6 lines / 3 patients as plotted, consistent with the text — but all calls are Tier 2/3 (M7).
- Panel A: define "artefact collapse ATM 82%→9%, ATR 77%→5%" (M-minor 14). Also note in the caption that this waterfall is for one line (OV2295) and give the panel-wide equivalent.
- Panel B: a 5 × 5 submatrix cannot show "first of 67" (M6c). The colour scale saturates at 0.8, so 0.82 and 0.88 are indistinguishable — extend it, and add the self-minus-best-other margin.
- Panel D: **good panel** — two-sided gain/loss, *n* = 11 stated on the axis. I verified all six frequencies quoted in §3 against `output/wes_cnv_arm_freq_patient.csv` (3q gain 82%, 20q gain 91%, 17p loss 82%, 8q gain 73%, 13q loss 64%, 19q gain 55%): **all six are exactly correct.**
- Panel D: with n = 11, each patient is 9.1%; say so, and state the log2 threshold used to call an arm gained or lost (M2).

### Fig. 5 — TOV21G hypermutation and mucinous authenticity
**Readable: yes**, except that panel B's sample identity is not stated.
- **Panel B is not labelled with the line it belongs to.** From context it must be TOV21G, but nothing in the figure says so and there is no caption. Add the line name and the variant count used.
- Panel A: the "indel fraction" colour legend shows only 0.1 and 0.2 with no endpoints; give the full range and the TOV21G value (M-minor 13).
- Panel C: bars are labelled by signature *group* but it is not stated whether the value is the maximum or mean cosine within the group. State it. Add bootstrap intervals and a proper refit (M-minor 12). Note also that panel C lists SBS6/44/15 whereas the text lists SBS44/SBS15/SBS20.
- Panel C: the clock group at 0.65 versus MMR-d at 0.88 is a modest separation in cosine space given how similar these signatures are; the panel should not be read as decisive, and the text's "no POLE signal" from 0.32 should be softened.
- **Panel D's column headers assert the conclusion** ("ovarian" / "GI-leaning") — circular (M13a). Use neutral labels.
- Panel D: the *z*-score reference set is unstated and decides whether the panel is informative (M13b). The legend runs 0→4 with unlabelled negative extent; give symmetric, labelled limits.
- Panel D: VOA8762's SATB2 is at the reference mean, so its intestinal call rests on CDX2 alone (M13c).
- **Fig. 5 is cited for the SWI/SNF reclassifications but contains no SWI/SNF panel** (M6a). Either add the panel here or correct the citation.

### Fig. 6 — ADC target atlas
**Readable: the RNA panel yes; the protein panel effectively no.**
- **The protein heatmap is nearly uniform** because it is plotted on an absolute 9–16 log2 scale while cross-line spread is ~0.34 log2 (Fig. 2B). It cannot support "Known associations recover in both RNA and protein" (M11). Row-scale it.
- The **DPEP3 protein row is a distinct light blue/grey (not detected) with no legend entry** for missing data. Add one.
- The heatmap fill (Reds) shares hue with the HGS subtype annotation bar; the subtype strip is hard to separate from the data. Use a different hue family for the annotation.
- "protein log2 abund." should specify what the abundance is relative to (PIS ratio?).
- Add *n* per subtype block to the panel (HGS 15, CC 6, EC 2, MC 3, MMMT 2, SCCOHT 2 as drawn), so a reader sees immediately that the non-HGSC group statements rest on 2–6 lines.
- I could verify the three "cross-assay-consistent" examples the text names (HER2→TOV3392D, NaPi2b→VOA12539, CDH6→OV3331): **all three read correctly in both panels.** I could not verify "the highest FRα expressers panel-wide are HGSC lines" — VOA8762/VOA8771 look comparable to TOV3133G.

### Fig. S1 — Protein PCA by subtype and TMT plex
**Readable: yes.**
- **Never cited in the text.**
- **The two panels reuse the identical palette for two different variables**: HGS orange = plex 1 orange, CC teal = plex 3 teal, MC gold = plex 4 gold, EC purple = plex 5 purple, SCCOHT navy = plex 2 navy. In adjacent panels this is a direct invitation to misread. Recolour panel B.
- Panel A shows essentially no subtype structure. That is an honest and useful negative result — and it is the panel that would let a reader calibrate "more weakly, in protein". Cite it, and report the protein silhouettes alongside (CC = −0.003, SCCOHT = 0.028 — M-minor 10).
- Panel B is a genuinely useful batch control and supports the low per-protein plex variance. But it should be shown together with the PC-level plex R² (0.215 adjusted on PC2), not instead of it (M3c).
- Add *n* to both panels. Note the single CC line at PC2 ≈ 20 that appears to drive PC2.

### Fig. S2 — Passage confounding
**Readable: yes.**
- **Never cited — and it should be one of the most prominent figures in the paper** (M8). Panel A shows passage explaining ~15% of PC1 conditional on site, exceeding the reported genome-wide subtype component. Panel B shows RNA and WES passages differing by up to ~20 for the same line.
- Panel A: add uncertainty (bootstrap or permutation) and explain why PC5 has a conditional value (2.4%) but no raw bar.
- Panel A: label the *y* axis units fully ("% of PC score variance explained by passage, univariate vs conditional on site").
- Panel B: give *n* (13 labelled lines), state why only these lines appear, and add a correlation coefficient with CI. Label axes "passage number".
- Both panels must be reflected in the text: add passage to the confounder models, and add "omic layers are not passage-matched for all lines (range of RNA-vs-WES difference: X to Y passages)" to the Reuse caveats.

### Fig. S3 — Line × line RNA correlation with dendrogram
**Readable: yes**; row/column labels are legible.
- **Never cited — and it bears directly on the central site-confounding claim** (M3d): the top-level dendrogram split separates CHUM (TOV/OV) from BC Cancer (VOA) lines almost perfectly. This must be shown and addressed, not omitted.
- **The colour scale is truncated at 0.7**, which exaggerates the block structure. State the actual range (all pairs are ~0.75–0.90) so the reader can calibrate; consider a scale starting at the observed minimum.
- **The family palette differs from Fig. 1B and Fig. 4C** (here 1369 = orange-red, 2295 = navy, 3133 = teal, 3291 = gold; there 1369 = purple, 2295 = teal, 3121 = orange, 3133 = blue, 3291 = magenta). Same variable, three figures, two palettes. Harmonise.
- Add a site annotation strip — it is the most informative annotation for this panel and it is missing.
- State the feature set (all 22,544 genes? HVGs? VST or TPM?) and the *n* (31 lines, which I verified by counting the labels).
- Substantive observation the authors should address: **families 1369 and 2295 do not co-cluster** here (OV1369-R2 at position 13 versus TOV1369 at position 8; OV2295 at 5 versus OV2295-R2 at 9), while family 3133 does. For sublines of one patient, that is worth a sentence — it is either biology (divergent primary/recurrent isolates, which is the point of those lines) or an identity concern, and the paper should say which it believes.

### Fig. S4 — Per-set protein presence patterns
**Readable: yes.**
- Cited. Good, complete panel: all 31 non-empty presence patterns are shown, and the counts sum correctly within each plex-count group (512 / 366 / 360 / 264 / 6856), which I verified against Fig. 2C and `prot_block_missingness.csv`.
- **The 71 proteins present in zero plexes are omitted without comment** (M9), yet they are inside the 1,573 "absent from at least one set" figure.
- **Panel B uses bars on a log10 axis.** Bar length is then not proportional to value, which is a known distortion; use points/lollipops, or a linear axis with a break.
- Panel A has no row labels or counts; alignment with panel B is implicit. Add the *n* to each row or a connecting rule.

### Fig. S5 — Genome-wide copy-number heatmap
**Readable: yes for broad patterns; individual bins are speckled at this rendering.**
- **Never cited** (M8) — despite being the natural support for all of §3's CNV claims and for the FGA ordering. Cite it.
- I verified the row count: 23 lines, matching "Twenty-three lines have CNVkit copy-number profiles", including TOV3121D (which is absent from Fig. 4C's 22 MAF lines). This is the figure that makes the 23-vs-22 distinction legible; use it.
- Colour: symmetric diverging blue–white–red centred at 0, −2 to +2. Appropriate and colour-vision-safe.
- The driver-locus labels are colour-coded (3q26/8q24/19q12/20q in orange, 10q/13q/17p in black), presumably gain versus loss, but this is never stated. Explain.
- State the bin size / segment resolution, and add a genomic-coordinate reference within chromosomes.
- The `family` strip annotates OV3291 alone as "family 3291" because its sibling has no WES. Note this, or mark such rows as "family member, sibling not profiled".
- The near-blank TOV81D row is a nice visual anchor for FGA = 0.02; point it out in the caption.

### Fig. S6 — ConsensusOV calls versus intrinsic strata
**Readable: yes.**
- Cited. Fully internally consistent: I verified DIF 6, PRO 2, IMR 3, MES 4 (total 15) and IMR+MES = 7 against the text's "differentiated 6, mesenchymal 4, immunoreactive 3, proliferative 2" and "7 of 15".
- **The stratum colours are incompatible with Fig. S8** — here Inflammatory = teal, Hypoxic-glycolytic = purple, Low-signaling = gold; in Fig. S8 Inflammatory = orange, Hypoxic-glycolytic = **teal**, Low-signaling = grey. The same word maps to two colours and teal maps to two different strata across two adjacent supplementary figures. I confirmed the *assignments* are identical between the two figures; only the palette differs. Harmonise — this is the single most confusing colour problem in the set.
- The figure conveys only categorical labels. Since the paper's argument is that the ConsensusOV labels are unreliable for pure cultures, **show the ConsensusOV per-class probabilities or margins** — a line called MES with a 0.34/0.33 margin makes the point far better than the label alone.
- "Intrinsic stratum" is used without definition in the figure; define it (and note that it is derived in Fig. S8).
- IMR (light orange) and MES (dark orange) are the two microenvironment-driven classes and are the closest pair in the palette; the in-bar text labels rescue this, but consider separating the hues.

### Fig. S7 — Marker effect sizes
**Readable: yes. This is the best-executed figure in the set.**
- Cited. I verified every number the text draws from it: 22 markers; 16 circles / 6 triangles; KRT20 *d* ≈ 3.05, CDX2 ≈ 2.75, SPP1 ≈ 1.45, HNF1B ≈ 1.35, PAX8 ≈ 0.5, MUC2 ≈ −0.28, MECOM ≈ 0; and exactly 8 markers at AUC ≥ 0.80 (SMARCA2, ZEB1, KRT20, CDX2, TFF3, TFF1, SPP1, HNF1B). All correct.
- **No confidence intervals** (M-minor 9). SMARCA2 *d* = −2.8 and AUC = 1.00 rest on n = 2 positives; KRT20 *d* = 3.1 on n = 3. Without CIs these are not comparable, and the largest effects come from the smallest groups.
- **Panel A plots signed *d* while panel B plots oriented AUC**, so SMARCA4/SMARCA2 point left in A and right in B. The text quotes median |*d*|. Harmonise, or annotate the loss-marker rows.
- Add *n* per intended subtype to the panel.
- **No EC markers** (M14). State that endometrioid identity is not marker-validated.
- Panel B's axis lower limit is unlabelled (first tick 0.6) though the reference line at 0.5 is clear.
- "Intended subtype" here versus "expected histotype" in the text — pick one term.

### Fig. S8 — Within-HGSC heterogeneity
**Readable: yes.**
- Cited (three times correctly, once incorrectly — the WES-coverage note in §1 should not point here; M8).
- **Panel A is circular**: points are plotted against the two scores used to define the strata and coloured by the resulting label (M15). Present the strata as a partition of this 2D space, not as an independent result — or derive them by an independent method (clustering on all Hallmark scores, protein, or CNV) and show the agreement.
- **Colour mapping conflicts with Fig. S6** (see above).
- Panel A's size encoding ("Prolif. z", legend −1/0/1) cannot encode sign in area, has no visual zero reference, and at least one point (OV2085) is smaller than the legend minimum. Use a separate small panel or a colour-plus-size double encoding.
- Panel A: several labels sit ambiguously between two points (OV2295-R2 / OV866-2; OV3133-R / TOV3133G). Add leader lines throughout.
- Panel B: **the corroboration is weaker than claimed** — only 2 of 5 "Hypoxic-glycolytic" lines show elevated PROGENy Hypoxia, while OV1946 (Inflammatory) does (M15). Quantify the agreement (e.g. a stratum-wise test or silhouette on PROGENy space) rather than asserting it.
- Panel B: colour scale is symmetric, diverging, centred at 0, labelled. Good.
- Both panels: add *n* = 15 and state that no clustering was performed on the PROGENy matrix (columns are ordered by the assigned stratum).

---

## 7. Internal inconsistencies and factual errors

| # | Location | Stated | Should be / conflict |
|---|---|---|---|
| 1 | Abstract | "Raw and processed data **are deposited** in public repositories" | Nothing is deposited; all Table 1 accessions are `[PLACEHOLDER]`; raw MS not yet obtained; WES FASTQs may not exist |
| 2 | Methods (proteomics), §1 | "8,430 proteins" | **8,429** (6,856 + 1,573 = 8,429; `prot_block_missingness.csv` has 8,429 rows; Fig. 2C and Fig. S4 both sum to 8,429) |
| 3 | §3 | CDK12 "2–3 of 11 patients (18–27%), **matching the ~3% TCGA-HGSC rate**" | 18–27% is 6–9× the TCGA rate; it does not match |
| 4 | §1 | pseudoalignment site difference "**does not translate into a difference in the number of genes detected**" | Mes-Masson 19,817 vs Huntsman 20,765 genes; Mann-Whitney p = 0.0098 (`rna_qc_metrics.csv`) |
| 5 | §2 | "On 30 lines and **8,212 shared genes**" | 7,893 per-gene correlations in `integ_rnaprot_cor.csv` |
| 6 | §2 | Silhouettes "HGSC 0.16, clear cell 0.12, mucinous 0.15 ... MMMT (0.74) and SCCOHT (0.82)" | These are the **RNA** values; modality is never stated. Protein values differ sharply (CC −0.003, SCCOHT 0.028) and are never reported |
| 7 | §2 / Fig. 3D | "subtype ≥ site in RNA (median 5.9% vs 3.5%)" | True for the median only; means are 14.9% vs 14.5% and IQRs overlap (`rna_variancepartition.csv`) |
| 8 | §2 | "site adds ≤0.2% of PC1 variance ... unique-subtype 42%" | Raw R² values; the adjusted joint R² (0.657 vs raw 0.737) is in the same file and not reported |
| 9 | §4 heading + text | Reclassification claims cited to "**Fig. 5**" | Fig. 5 has no SWI/SNF panel; SMARCA4/SMARCA2 per-line RNA/protein appear in no figure |
| 10 | §4 | "ranking **first of 67** DepMap ovarian lines ... (Fig. 4)" | Fig. 4B is a 5 × 5 matrix and cannot show this |
| 11 | §1 (WES note) | "(**Fig. S8** is contingent on this)" | Fig. S8 is the HGSC heterogeneity figure; either a WES QC figure is missing or the number is wrong |
| 12 | Text vs figures | Figs. S1, S2, S3, S5 exist | None is cited anywhere in the text; Table S3 is also uncited |
| 13 | Fig. 3D panel title vs analysis | "RNA (22,544 genes)"; Methods "22,544 genes" | `rna_variancepartition.csv` `n_features` = 22,542 |
| 14 | Fig. 3E vs §2 and Fig. S7 | 23 marker rows in Fig. 3E | 22 markers in the text and in Fig. S7; MKI67 ("Prolif") is the extra row and is not one of the 22 |
| 15 | Methods (proteomics) vs Fig. 2D | "channel 10 carries **a** separate technical-replicate 'bridge' sample" | Four *different* lines across the four links (VOA10816, TOV1369, VOA3993, TOV3133G), and VOA3993 is external to the resource (`prot_bridge_cor.csv`, `external = TRUE`) |
| 16 | Fig. 2D | "VOA3993\*" | Not among the 42 lines in Fig. 1B; asterisk undefined |
| 17 | Methods vs figures/metadata | "three centres: CHUM ..., BC Cancer / OVCARE ..., and **OHRI** (Vanderhyden BIN67)" | Fig. 2A and `supplement_per_line.csv` label the third site "**Huntsman/Vanderhyden**"; Fig. 3B collapses to two sites |
| 18 | Abstract / Table 1 / Fig. 1A | "whole-exome sequencing (**23 lines**)" | 23 lines CNV, **22** lines mutation calls (Table 1, Fig. 4C, Fig. 5A) |
| 19 | Fig. 5C vs §"Rare-subtype" | Fig. 5C: "SBS6/44/**15**" | Text: "SBS6 cosine 0.88, plus SBS44/SBS15/**SBS20**" |
| 20 | Fig. 4C vs §3 | TP53 patient bar = **12** | §3 says "all **11** HGSC patients"; consistent only if the bar counts non-HGSC patients too (TOV112D), which no caption states |
| 21 | Fig. S6 vs Fig. S8 | Stratum colours: Inflammatory = teal / Hypoxic = purple / Low-sig = gold | Fig. S8: Inflammatory = orange / Hypoxic = **teal** / Low-sig = grey (assignments agree; palettes do not) |
| 22 | Fig. 1B vs Fig. 1B/S3/Fig. 4C | Family palette in Fig. 1B and Fig. 4C | Fig. S3 uses a different family palette for the same five families |
| 23 | Methods §"Cell lines" vs Fig. 1C | clear cell "n=7: 2 CHUM + 5 BC Cancer" (Methods, confounder section) | Fig. 1C and the histotype composition give clear cell = 8 lines / 8 patients; n = 7 is the RNA subset and this is never stated |
| 24 | Fig. 1B vs Fig. 4C / Fig. S5 | Fig. 1B appears to mark WES-SNV for both TOV3121EP and TOV3121D | Fig. 4C (22 MAF lines) includes only TOV3121EP; Fig. S5 (23 CNV lines) includes both. Verify Fig. 1B |
| 25 | Throughout | HGSC / MMMT-carcinosarcoma / low-grade serous | Figures use HGS / MMMT / LGS; abbreviations never defined |
| 26 | Whole manuscript | — | **No figure legends and no supplementary table legends exist** (verified in the source `.docx`) |

---

## 8. Questions for the authors

These are the questions whose answers would change my assessment.

1. **Are the accessions obtainable before publication, and do the WES FASTQs exist?** If the WES record is processed-calls-only from an unspecified capture kit against unmatched public normals, what exactly can a reuser do with it, and is the abstract's "whole-exome sequencing (23 lines)" the right way to describe it?
2. **What is the full proteomics search and quantification workflow** (instrument, MS2 vs MS3, search engine and version, database, FDR, interference filter, peptide-to-protein rollup)? The paper's central proteomics interpretation — that the RNA–protein correlation of 0.40 is a TMT compression ceiling rather than data quality — depends on the acquisition mode, since MS2-based TMT compresses substantially more than SPS-MS3.
3. **What do the analyses look like on the 28 patient representatives** rather than 31 lines? Specifically: does the PC1 commonality decomposition still give unique-site ≈ 0.2%, and do the HGSC and clear-cell silhouettes hold up? And what does a permutation null (subtype labels shuffled within site) give for unique-subtype on PC1?
4. **How much of the leading-PC structure is passage rather than biology?** Fig. S2A puts passage at ~15% of PC1 conditional on site. What are the subtype and site components when passage is added to the joint model, and what is the range of RNA-vs-WES passage discordance across the 13 lines in Fig. S2B?
5. **What is the SMARCA4/SMARCA2 evidence, per line, in RNA and protein?** Specifically the BIN67 claim that mRNA is retained while protein is second-lowest: what are the actual values, what is the rank margin, and what is the uncertainty given the ~0.3 log2 cross-line protein spread the paper itself reports?
6. **What is the self-match margin against all 67 DepMap ovarian lines** — self-match Spearman minus best non-self Spearman, per line? "Reciprocal-best" with a margin of 0.01 and with a margin of 0.15 are very different claims about identity.
7. **What are the 71 proteins present in zero TMT plexes**, and are they in the deposited matrix?
8. **In Fig. 5D, across what set are the mucinous z-scores computed** — the three mucinous lines or all 30? And can the VOA8762/VOA8771 GI-leaning call be reproduced against an external anchor (TCGA/CPTAC mucinous ovarian versus colorectal/gastric), given that the paper's own MUC2 result shows mucinous differentiation markers drift in 2D culture?
9. **What log2 threshold defines an arm-level gain or loss** in Fig. 4D, and how sensitive are the six quoted frequencies to it — particularly given the flagged median-centring problem in TOV3121D and TOV2929D and the unconfirmed capture-kit concordance of the pooled normals?
10. **What are the culture conditions for each line**, and under what ethics approvals are these patient-derived lines held and distributed?
11. **Why is there no endometrioid marker set?** With TOV112D reclassified, VOA4395 is the sole EC line and its histotype is nowhere validated — is there evidence for it, or should the EC coverage be described as n = 1, unvalidated?
12. **Is a signature refit with bootstrap intervals consistent with the SBS6 call for TOV21G**, and how many variants entered the spectrum after PASS and population-AF filtering? Given ~206 residual germline coding variants per line, what fraction of TOV21G's 1,416 could be germline?
13. **Given that VOA5436 is flagged as expressing "a strong clear-cell program"**, is the carcinosarcoma coverage effectively n = 1 (VOA5217)? If so, the gap-filling claim in Background & Summary should be revised in the same way the mucinous claim was.
