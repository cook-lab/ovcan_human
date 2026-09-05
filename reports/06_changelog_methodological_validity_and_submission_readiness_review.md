# Audit response changelog — methodological validity and submission-readiness review

**Companion to:** `reports/06_methodological_validity_and_submission_readiness_review.md`
**Response date:** 27 July 2026 (round 1), same date (round 2)
**Responder:** Cook Lab (Claude)
**Scope:** every finding in both review rounds, with the action taken, the evidence for it, and the
reasoning where no action was taken
**Verification:** all edited analysis/figure scripts rerun to exit 0; report rebuilt (38 pp); manuscript
`.docx` rebuilt. **A full clean-room re-execution has not been rerun** — see §11.

> **Round-2 note.** A second review round accepted the round-1 response and raised six refinements: the
> protein "reliability" claim was still overclaimed, the ADC comparison was not truly per-class, the
> reference analysis was mislabelled a sensitivity analysis, the Spearman inference and "proxy"
> terminology needed narrowing, one disagreement was rhetorical rather than substantive, and a
> protein-key join defect needed repair. **All six were accepted and implemented.** Sections below are
> updated in place; **§12 records the round-2 changes and their effect on the numbers.** Where a
> round-1 statement was itself imprecise, it is corrected in place and marked.

---

## 1. Executive summary of the response

The review made 13 numbered methodological findings, an 11-row claim-validity matrix, 8 reporting
principles, and 12 submission blockers. Of these:

| Disposition | Count | Character |
|---|---:|---|
| **Accepted and implemented as code changes** | **5** | Rerun scripts, new deposited tables, regenerated figure |
| **Accepted and implemented as text changes** | **9** | Claim language, unit discipline, caveat placement |
| **Accepted — factual error corrected** | **2** | One wrong citation scope, one false statement of fact |
| **Already satisfied before the review** | **5** | Wording tightened; no substantive revision needed |
| **Accepted in substance, implemented differently** | **3** | Reasoning given in §5 |
| **Accepted, evidential boundary clarified** | **1** | §6 — recategorised from "not accepted" in round 2 |
| **Ran a diagnostic the review suggested; it flagged a risk rather than clearing one** | **1** | §4 |
| **Round-2 refinements, all accepted** | **6** | §12 |

**Nothing in either round was rejected.** Round 1 recorded one item as "not accepted"; round 2 showed
that was a mischaracterisation of the review rather than a disagreement with it, and it is now §6.

**Two findings were the serious ones.** The proteomic repeatability conversion was a genuine
statistical error whose correction *reverses the direction* of a headline claim, and the v4 draft
contained a sentence asserting that authentication work had been performed when it had not. Both are
detailed below, and for the second I traced and closed the upstream cause rather than only the
sentence.

**One observation on scope.** Several review items were aimed at text that
`reports/01_multiomic_characterization_results.md` had already qualified correctly but that
`docs/manuscript/v4/OvCAN_data_descriptor_v4.md` had not — the report and the manuscript had drifted.
Those items are recorded below as "already satisfied (report) / implemented (manuscript)", because
crediting them as fully satisfied would misrepresent the manuscript's state.

---

## 2. Accepted and implemented as code changes

### 2.1 Proteomic bridge repeatability lacked the replicate-difference adjustment

**Review finding 5.** Upheld in full. This was a real error, not a wording problem.

Each bridge link re-runs one sample in the adjacent isobaric set, so `primary − bridge` is a
*difference of two measurements*. Where the two carry equal independent error variance,
`Var(difference) = 2 × Var(single)`, so a per-measurement SD requires `SD_single = SD_diff / √2`. The
deposited quantity was `100 × (2^SD_diff − 1)`, which applies neither that adjustment nor the standard
lognormal conversion, and it was labelled a repeatability CV.

| Link | Plexes | SD of paired differences | Old "CV" | **Per-measurement SD** | **Per-measurement CV** |
|---|---|---:|---:|---:|---:|
| `VOA10816` | 1→2 | 0.2683 | 20.44% | **0.1897** | **13.21%** |
| `TOV1369` | 2→3 | 0.2106 | 15.72% | **0.1489** | **10.35%** |
| `VOA3993` *(external)* | 3→4 | 0.2536 | 19.22% | **0.1793** | **12.48%** |
| `TOV3133G` | 4→5 | 0.2304 | 17.32% | **0.1629** | **11.33%** |
| **Range** | | 0.211–0.268 | *15.7–20.4%* | **0.149–0.190** | **10.4–13.2%** |

**Implemented in** `scripts/19_proteomics_dynamic_range.R` §4b. The corrected columns `sd_single` and
`cv_pct_per_measurement` are deposited; the old value is retained as `sd_diff_cv_pct_legacy` and
flagged deprecated in a `units_note` column, so the superseded number in earlier drafts stays
traceable rather than vanishing. A `replicate_caveat` column records that the four links are four
non-identical samples, not replicated aliquots of one reference.

The abundance-decile series carried the same error and is corrected the same way
(`prot_cv_by_abundance.csv` gains `bridge_sd_single` and `bridge_cv_pct_per_measurement`): the bridge
series across deciles 1→10 is now **21.3% → 5.2%**, not 34.7% → 7.6%, and is therefore on the same
per-measurement scale as the vendor CV plotted beside it (11.0% → 3.0%).

### 2.2 "Technical noise exceeds biological spread by 2.4–3.1×" was an artefact of the comparison

**Review finding 5, second half.** Upheld in full, and the consequence is larger than the review
states: **the corrected comparison reverses the claim's direction.**

The retired ratio divided a 95% limits-of-agreement span — which is `3.92 × SD` on the
paired-difference scale — by a median cross-model IQR, which is `1.35 × SD` on the single-measurement
scale. That inflates the ratio by roughly 2.9× from the choice of dispersion measure alone, which
accounts for the entire "2.4–3.1×" effect.

Like-for-like, SD against SD, on the 7,896 paired genes:

| Quantity | Value |
|---|---:|
| Median cross-model protein SD | **0.287 log2** |
| Per-measurement technical SD | 0.149–0.190 log2 |
| **Observed cross-model SD ÷ technical SD** | **1.51–1.93×** |
| **Technical share of observed cross-model variance** | **27.0–43.8%** |
| Implied biological SD (by subtraction) | 0.215–0.245 log2 |
| Approximate variance ratio (1 − technical share) | 0.56–0.73 |

**Implemented in** `scripts/19_proteomics_dynamic_range.R` §4b-ii, deposited as
`output/prot_noise_vs_biology.csv`.

**⚠ Round-1 phrasing corrected in round 2.** This changelog originally called the last row a
"single-measurement reliability" and the column was named `reliability_icc`. **It is not an intraclass
correlation** and is now named `approx_reliability_ratio`, with the four approximations that separate it
from one stated in the script and in both documents: homoscedastic error across the abundance range
(false — see below), a technical SD pooled *across proteins* compared against an observed SD computed
per gene and then medianed *across genes*, four non-identical bridge samples rather than replicated
aliquots, and independent additive error so that the biological component is reached by subtraction.
Defensible phrasing: *under a homoscedastic independent-error approximation, bridge variability
corresponds to roughly 27–44% of the variance of a typical protein's observed cross-model measurements.*

**The homoscedasticity assumption was the weak part, and dropping it changes the answer.** Matching the
technical estimate to the observed spread *within* mean-abundance deciles — both sides decile-local —
replaces one number with a nine-fold gradient (`output/prot_cv_by_abundance.csv`):

| Abundance decile | Technical SD | Observed SD | **Technical variance share** | Approx. variance ratio |
|---|---:|---:|---:|---:|
| 1 (lowest) | 0.303 | 0.369 | **67.5%** | 0.33 |
| 3 | 0.211 | 0.317 | 44.4% | 0.56 |
| 5 | 0.165 | 0.273 | 36.5% | 0.63 |
| 7 | 0.119 | 0.257 | 21.3% | 0.79 |
| 10 (highest) | 0.075 | 0.272 | **7.5%** | 0.92 |

The global 27–44% figure sits in the middle of that gradient, understating the problem at low abundance
and overstating it at high abundance. **The reuse conclusion is now explicitly abundance-conditional:
there is no dataset-level answer to "can I interpret one protein in one model".**

**Retained and quantified rather than dropped.** The review's recommendation was to remove the ratio
entirely. I kept a corrected version, because the review's own recommended replacement (report bias,
SD and LoA and stop) leaves a reuser without an answer to the question the claim existed to answer.
The abundance-matched gradient answers it, conditioned on the one thing a reuser always knows about
their protein — roughly how abundant it is.

**⚠ The ADC comparison was also corrected in round 2.** Round 1 swapped in the ADC-target median
*observed* SD but kept the **global** technical SD, so it was not a per-class statement — a valid
objection. The bridge differences restricted to the ADC target rows now give a direct estimate:
technical SD **0.134 log2** from **32 differences** against a median cross-model SD of **0.600**, i.e. a
technical share near **5%** rather than the 6–10% reported in round 1. Because it comes from the targets
themselves this is genuinely per-class, and it is quoted with its n because 32 differences is few. The
abundance-matched deciles are the cross-check and agree, these eight being abundant proteins.

### 2.3 Marker recovery needed a correlation-preserving null

**Review finding 3.** Upheld. The exact binomial treated 25 markers as independent Bernoulli trials at
p = 1/3; markers intended for the same histotype are co-regulated (`PAX8`, `WT1` and `MUC16` all track
one serous programme), so a single relabelling moves several together and the binomial p is
anti-conservative by an unknown amount.

**Implemented in** `scripts/04_rna_markers_genesets.R` §2a — 20,000 joint permutations of the histotype
labels at the project seed 1234, preserving marker–marker correlation, group sizes and the
absolute-expression floor in the up rule. Deposited as
`output/rna_marker_recovery_permutation.csv`.

| Analysis unit | Observed | Null mean | Null 95th pct | **Permutation p** | Binomial p |
|---|---:|---:|---:|---:|---:|
| 31 line models | 16 / 25 | 6.70 | 11 | **0.0013** | 0.0016 |
| **28 patient representatives** | **15 / 25** | 6.82 | 12 | **0.0046** | 0.0056 |

The qualitative conclusion strengthens: observed counts exceed the null 95th percentile on both units.
One detail worth recording — the null mean of 6.7–6.8 of 25 is **~27%, not the binomial's 33%**,
because the expression floor makes landing strictly harder than a rank test alone. The binomial was
therefore mis-specified in two directions at once, and the permutation absorbs both.

The review ran this test independently and obtained p = 0.0019 and 0.0052. The difference from the
deposited 0.0013 and 0.0046 is the seed alone; the project seed 1234 is now canonical, and the review's
values should be treated as the independent replication they were.

### 2.4 Negative transcript–protein correlations were reported as established

**Review finding 8, second half.** Upheld. The claim that for 813 genes "a transcript measurement is
not a proxy for protein in any direction" asserts a relationship the estimates do not support: with 10
to 30 paired models the sampling SE of a Spearman ρ is roughly `1/√(n−3)` — 0.38 at n = 10, 0.19 at
n = 30 — so most negative point estimates are indistinguishable from zero.

**Implemented in** `scripts/12_rna_protein_concordance.R` §5c. Every per-gene estimate now carries a
two-sided p and a BH-adjusted q, both deposited in `integ_rnaprot_cor.csv`, with tier counts in
`integ_rnaprot_negative_genes.csv`.

**⚠ Two round-2 corrections.** Round 1 used a **Fisher z transform scaled by 1/√(n−3)** — a
Pearson-derived approximation applied to a rank correlation. The conventional choice is the asymptotic
Spearman **t approximation with n − 2 degrees of freedom**, which is what `cor.test(method = "spearman")`
falls back to once ties are present, and ties are ubiquitous here. Switched. The reviewer's independent
check predicted the effect exactly:

| Tier | Round 1 (Fisher z) | **Round 2 (conventional)** | % of per-gene set |
|---|---:|---:|---:|
| Inverse point estimate (ρ < 0) | 813 | **813** | 10.30% |
| Inverse and nominally significant (p < 0.05) | 51 | **50** | 0.63% |
| **FDR-supported inverse association (q < 0.05)** | 30 | **29** | **0.37%** |
| FDR-supported positive association (q < 0.05) | 3,626 | **3,599** | 45.59% |

Nothing substantive turns on the change, but the conventional test is the one a reader can reproduce.

**And the terminology was still wrong in round 1.** Describing the FDR-supported subset as genes for
which "RNA fails as a proxy" conflates two different things: **an inverse association is not an absence
of proxy value.** A protein that reliably falls as its transcript rises is predictable — simply in the
other direction. Proxy value additionally depends on predictive error and external validation, neither
assessed here. The tiers are now named for the association they establish (`inverse point estimate`,
`FDR-supported inverse association`), and the deposited `interpretation` field says so explicitly. The
813 remain candidates for discordant regulation; the 29 are inverse associations.

### 2.5 Same-donor sublines in the concordance estimates

**Review finding 8, third part** — the review asked for a patient-representative sensitivity analysis
and did not predict its direction. **Implemented in** `scripts/12_rna_protein_concordance.R` §5d,
deposited as `integ_rnaprot_patientrep_sensitivity.csv`.

| Statistic | All 30 dual-layer models | 27 patient representatives | Shift |
|---|---:|---:|---:|
| Per-model Spearman (median) | 0.4080 | 0.4227 | **+0.0147** |
| Per-gene Spearman (median) | 0.3971 | 0.4168 | **+0.0197** |
| Fraction of per-gene estimates negative | 0.1030 | 0.0970 | −0.0060 |

Both medians move **up** under collapse, so repeated donors were very slightly *deflating* the
concordance rather than inflating it — the opposite of the concern. This is a clean null on the
review's point, and it is now reported as such rather than left as an assumption.

---

## 3. Accepted and implemented as text changes

| # | Review finding | Action | Where |
|---|---|---|---|
| 1 | Centre–histotype confounding limits causal interpretation | Language made associational throughout. New §4.2 subsection states precisely what commonality analysis establishes (centre adds nothing beyond histotype) and what it cannot (apportion the 0.311 shared component, since no design cell has the two varying independently). "Biology, not batch" and "Histotype, not the contributing centre, structures the expression data" both retired. New **limitation 5b** added to the register. | Report callout, §4 intro, §4.2, limitation 5b; manuscript Technical Validation §2, Fig. 3B legend |
| 2 | The patient, not the derived line, is the inferential unit | Split made explicit and load-bearing: **28 patient representatives for inference** (marker effect sizes, the recovery permutation, DE/enrichment sensitivity), **31 line models for description** (PCA displays, silhouettes, heatmaps, selection tables). Stated in the denominator inventory, in Methods, and as two separate rows of the reuse-defaults table. Rationale for collapse over a donor-aware mixed model recorded (three replicated RNA families will not estimate donor variance precisely). | Report §1.2, §5.3, §10.5; manuscript Methods (new paragraph), Technical Validation §2 |
| 4 | Gene-set "recovery" is uneven | The report already graded two of five programmes as robust. Added the missing disclosure: programme definitions are regular expressions and the reported term is the *best-scoring matching term*, so BH covers the GO result set and not the choice of term — the programme is prespecified, the term is not. Manuscript gained the grading it lacked entirely. | Report §5.2; manuscript Technical Validation §2 (new paragraph), Methods |
| 6 | Narrower protein spread is multifactorial | "Isobaric ratio compression of 3.34× sets the concordance ceiling" retired as a section heading and as a claim. Reframed: the matrix has a 3.3-fold narrower cross-model spread, **consistent with** known TMT compression and post-transcriptional regulation, with the four competing mechanisms named. The weaker true statement is retained — 0.40 is where this assay pair lands and matches published benchmarks. | Report §6.2; manuscript Technical Validation §2, Methods |
| 8 | Concordance is descriptive, not an accuracy certificate | The per-model/per-gene agreement is no longer called an "internal consistency check" — the estimands differ, so their numerical closeness is a coincidence of scale. The per-gene value is now named as the relevant reuse statistic. | Report §6.1; manuscript Technical Validation §2 |
| 9 | Tumour-only exomes are intrinsically qualified | Somatic-confidence tier now **defined as a prioritisation heuristic**, not a validated somatic classifier, with the reason (no matched normal, so even Tier 1 is a high-confidence candidate). *BRCA1/2* phrasing corrected to "no defensible *somatic* call can be made from this dataset", with the explicit consequence that absence of a call carries no information about the models' status. | Report §7.2; manuscript Methods, Technical Validation §3 |
| 10 | Copy-number conclusions depend on undocumented capture compatibility | 30-fold FGA contrast reframed as a **panel range and positive control** (18 models against 1) rather than a histotype effect estimate. Relative-total-copy and pooled-normal caveats attached at the point of claim rather than only in Methods. New **limitation 9b**. | Report §7.3 (heading and body), limitation 9b; manuscript Technical Validation §3 |
| 11 | Identity evidence is uneven | SWI/SNF framing changed to candidate subunit-loss flags with uncalibrated, project-specific thresholds, with the gene- and layer-specific evidence named as the primary record and the aggregate count as summary. `TOV21G` explicitly a **candidate** MMR-deficient/MSI-high model with the four reasons attribution is not possible. Mucinous verdicts explicitly provenance flags, not reclassifications. | Report §9.3; manuscript Technical Validation §3, §4 |
| — | *Scientific Data* genre and reporting fit | Added the missing required sections: **Data Availability, Ethics, Author Contributions, Competing Interests, Funding, Use of Generative AI, References**. Abstract cut **176 → 164 words** (title unchanged, 93 chars). Data Records gained the three properties that govern reuse (per-matrix denominators, structural missingness, donor structure). Usage Notes expanded from four to six analysis-limiting properties and made technical. Gap register grown **15 → 22 items**. Four new supplementary tables specified (**S14–S17**) for the permutation null, bridge agreement, per-gene concordance with q values, and the reference-mismatch sensitivity. | Manuscript throughout |

---

## 4. The reference-mismatch check — a design-alignment diagnostic, not a sensitivity analysis

**Review finding — claim-validity matrix row 2, and submission blocker 10.** The review recommended a
matched-index sensitivity analysis, on the expectation that it would *strengthen* the RNA quality claim.

**⚠ Round-1 framing corrected.** This section originally opened *"I implemented it and it did the
opposite."* **That is inaccurate, and the round-2 review was right to flag it.** What I implemented is
not a matched-index sensitivity analysis. It is a **design-alignment diagnostic**: it asks whether the
silently dropped transcripts are evenly distributed with respect to the terms the structure analyses use.

`scripts/01_rna_load_qc.R` §3c, deposited as `output/rna_reference_sensitivity.csv`. The
Ensembl-104-era index and release-105 transcript map disagree on 3,529 of 185,299 targets, and tximport
drops them without warning. Regressing each model's dropped-TPM fraction on the design terms:

| Term | Levels | R² | Kruskal–Wallis p |
|---|---:|---:|---:|
| Contributing centre | 3 | **0.330** | **0.0033** |
| Histotype | 6 | **0.388** | **0.0100** |

The two terms are collinear here, so these are not independent findings.

**What this establishes, stated precisely.** The omitted transcripts are not evenly distributed with
respect to centre and histotype, so the mismatch sits on the same axis §4 of the report interprets.
That is a reason for caution. It is **not** evidence that the mismatch manufactured any observed
structure, and it cannot be: **biologically distinct histotypes would be expected to express different
amounts of any transcript set, including this one**, so the association is exactly what a purely
biological model predicts as well. The diagnostic flags a risk; it does not measure a harm, and it
cannot be read in either direction as evidence about §4's conclusions. Round 1 over-read its own result
as a negative finding.

Three things bound the risk: the absolute quantity is small and tightly bounded (1.60–3.41% of TPM, a
1.81-percentage-point range across all 31 models), the reference pair is identical for every model so
none is advantaged by a different annotation, and the affected transcripts are absent from the
gene-level matrix entirely rather than mis-assigned to the wrong gene.

**Action taken.** Matched-index re-quantification is retained as the actual sensitivity analysis and
carried as a "would change the science" open item (report §12 item 5; manuscript gap 16), **and it has
not been run.** The diagnostic is labelled as such in the report (§2.1, limitation 6), in the
manuscript Methods and in the deposited file, so it cannot be mistaken for the sensitivity analysis it
is not.

---

## 5. Accepted in substance, implemented differently

### 5.1 `consensusOV` and the within-HGS strata

**Review recommendation:** remove from the main descriptor, or deposit as exploratory annotation with a
strong warning.

**What I did:** kept both in the main text, with the warnings strengthened, and added a gap item to
relocate Fig. 6B and its paragraph to the Supplementary Information at submission.

**Reasoning.** In the current draft `consensusOV` does not function as a validation claim — it appears
under "Two checks bound what the expression layers establish about identity", i.e. as a *limitation
disclosure*. Deleting it would remove a caveat, not a claim, and would leave a reuser who finds the
deposited `consensusov_calls.csv` with no guidance on why the calls are unreliable on pure cultures. I
therefore added the classifier's provenance (developed for bulk tumours containing stromal and immune
compartments, so ambiguity is expected), the statement that the calls do not validate identity, and
the note that they are deposited as an exploratory annotation with probabilities and margins rather
than bare labels. The HGS strata paragraph is now labelled **Exploratory** in both the text and the
figure legend, with all five reasons named (k not prespecified, no stability analysis, 15 models from
12 patients, labels drawn from the clustering space, partial PROGENy agreement).

**The relocation itself is deferred, and flagged rather than done.** Moving Fig. 6B means re-cutting a
validated composite figure, and doing that unprompted risks breaking a figure pipeline that has been
pixel-verified against a pre-edit backup. Gap 22 records the decision for the PI.

### 5.2 The limitations register

**Review principle 1 and layer-readiness table** imply reorganising limitations by layer. The register
is ordered by how much each item constrains reuse, and §11 already carries a number on every row. I
added the two missing items (**5b** centre–histotype non-identifiability, **9b** relative copy number
and unconfirmed capture compatibility) with letter suffixes and updated the header count 23 → 25,
rather than renumbering, because the existing numbers are cited across the report, the manuscript and
`ANALYSIS_LOG.md`. Renumbering would have silently invalidated those cross-references.

### 5.3 Retaining a corrected noise-to-signal statistic

Covered in §2.2 above. The review recommended deleting the ratio; I replaced it with a reliability
figure and an ADC-panel-specific contrast, on the grounds that the underlying reuse question is
legitimate and the original statistic's defect was the comparison, not the question.

---

## 6. Accepted, with the evidential boundary clarified

### 6.1 The commonality analysis — what it establishes, and what it does not

**⚠ Round-1 categorisation corrected.** This section was originally headed "Not accepted" and
characterised the review as saying the commonality analysis "adds nothing". **It said no such thing** —
the review explicitly called the analysis *"a useful improvement over a simple PCA colour plot"*. The
round-2 review was right that this was a mischaracterisation, and the disposition is recategorised as
**accepted, with the evidential boundary clarified**. It is not a disagreement.

The substantive position, which both rounds agree on: the partition establishes that **centre carries
no explanatory signal beyond histotype** — adjusted unique centre negative on all five leading
components, permutation p of 0.49–0.86 — and it **cannot apportion the 0.311 shared component**, because
no design cell has histotype and centre varying independently.

**One over-read of my own has also been fixed.** Round 1 wrote that "two independent alternatives are
excluded: centre as a source of variance beyond histotype and duplicate subline genomes." The second
half invites a wrong inference, and round 2 caught it: **donor collapse rules out duplicate donor
families as the source of the structure; it does not rule out a centre-aligned batch effect**, since
the 28 representatives still come from the same centres in the same proportions. Likewise the
commonality partition excludes centre as a source of *additional* variance, not a centre-aligned batch
effect, which by construction lives in the shared component. Report §4.2 now says exactly this, and
the summary callout no longer implies otherwise.

The net effect is that **neither of the two exclusions addresses the confounding**, and the report says
so rather than letting two narrow negative results accumulate into an implied positive one.

---

## 7. Already satisfied before the review

These items required wording tightening only. Recording them is not a defence — it is needed so a
future reader does not conclude that a substantive gap existed where one did not.

| Review finding | Pre-existing state |
|---|---|
| **4** Gene-set recovery uneven | Report §5.2 was already titled "Two of five expected pathway programmes recover robustly", already graded the other three "Suggestive"/"not recovered", and already warned that the grade is stable while the best term is not. The *manuscript* had no grading at all, which is where the fix landed. |
| **7** Proteomic structure retains plex/site structure | Report §4.1 already gave protein PC1 adjusted R² 0.464 histotype vs 0.315 site and flagged plex surviving adjustment on PC2 (0.215) and PC5 (0.229); §4.8 already gave the negative clear-cell (−0.004) and near-zero SCCOHT (0.028) protein silhouettes with the n = 2 caveat. The review's point that high bridge correlations cannot rule out plex effects was already §3.3's own opening argument ("Lead with agreement statistics, not correlation… r near 0.99 discriminates nothing"). |
| **9** Tumour-only exomes | Already "coding candidates" throughout; already stated the median burden of 206.5 is an **upper bound** on residual germline, not an estimate; already reported tier composition with every frequency. The tier-as-heuristic definition and the *BRCA1/2* phrasing were the genuine additions. |
| **10** Capture compatibility | Already disclosed in §2.3, in the Fig. S5 legend and in limitations 3–4, including that the missing BED is the stated reason no exome-to-genome renormalisation was applied. |
| **12** Variance partitioning and passage are descriptive | Already stated that the patient term is a design artefact fitted with 28 levels on 31 observations and "must not be reported as a finding", with the restricted-term collapse (27.27% → 0.758%) shown. |

---

## 8. Two factual errors corrected, and one root cause closed

### 8.1 `BIN67` was not reclassified by Karnezis 2021

**Review finding — claim-validity matrix row 11.** Upheld. Karnezis et al. 2021 (PMID 33328126)
reassigned **`COV434`** to SCCOHT and **`TOV112D`** to dedifferentiated ovarian carcinoma. `BIN67` was
already an established SCCOHT model and was not reclassified in that study, so counting it made the
tally two when it is one.

Corrected in the report callout, §1.4, §9.3 and §12; in the manuscript Technical Validation §4; and in
`docs/manuscript/v4/LEDGER.md`, `TRIAGE.md`, `CHECKLIST.md` and `VENUES_society_data.md`, each of
which had propagated it.

The verb was also wrong. Molecular profiles cannot *perform* a histopathological reclassification.
What the three layers establish is that `TOV112D` carries a SWI/SNF-null rather than an endometrioid
pattern — independent molecular evidence **consistent with** the published reassignment, arrived at
without reading the pathology. The multi-layer argument the claim existed to support is untouched:
`BIN67`'s *SMARCA4* mRNA is mid-panel (rank 23/31) while its protein is second-lowest (rank 2/31), so
an RNA-only authentication would return a false negative.

### 8.2 The v4 draft asserted authentication work that was never performed

**Review finding 11.** Upheld, and this was the most serious item in the review. The draft opened its
identity section: *"Short-tandem-repeat profiles and mycoplasma clearance were obtained for the panel
and are given in Supplementary Table S4."* **Neither was performed.** The sentence directly
contradicted the project record, `reports/01` §9.2, and the gap marker sitting immediately after it.

Removed and replaced with the limitation: no in-house STR or mycoplasma testing; 30 of 42 models carry
a third-party Cellosaurus *reference* profile and 12 have none; a reference profile documents what a
model should match, not that these stocks do; 25 of 26 named sources are a single personal
communication. The section now states explicitly that **the panel must not be described as
authenticated**, and Usage Notes tells a reuser to obtain STR and mycoplasma status before relying on
the annotations.

**Root cause, closed.** The sentence was not a slip. `docs/manuscript/v4/LEDGER.md` — the content
authority the draft was written against — contained, as a planning convenience:

> *"Assumption throughout: all four Bin-A gap classes are closed (STR + mycoplasma; Morin TMT
> parameters + raw MS; archive/PI recoverables; centre metadata). No `⟨TO OBTAIN⟩` markers exist."*

and, under §4, an **H**-tagged (headline) bullet reading *"STR profiling and mycoplasma testing were
performed on all models, with results in Table S4. Lead the section with this."* Writing the paper as
it would read once the requests landed made every dependent bullet indistinguishable from an
established fact, and one was drafted verbatim.

Both are now retired in place, with the old text quoted so the error cannot be silently reintroduced,
and a rule added to the ledger's tag legend: **a bullet anticipating a checklist item closing must be
tagged conditional, never as an H-tag headline.** The corresponding `CHECKLIST.md` and `TRIAGE.md`
items were re-tensed from accomplished to conditional.

---

## 9. Findings the review did not raise

Three issues surfaced while implementing the above.

1. **A two-number discrepancy between the report and the manuscript.** The report said 98.0% of paired
   genes have a smaller protein than transcript IQR; the manuscript said 97.8%. Both traced to
   `scripts/34_fig2_qc.R`, where `mean(dr$iqr_ratio < 1, na.rm = TRUE)` silently dropped 14 genes with
   a non-finite ratio, shrinking the denominator. Changed to `mean(dr$prot_iqr < dr$rna_iqr)` over all
   7,896 genes, and both documents now read **97.8%**.
2. **`docs/manuscript/v4/build_docx.sh` was fragile.** Its front-matter strip matched a fixed pattern
   ending at the first italic block and asserted that the result began with the Authors line, so adding
   a second provenance note broke the build. Re-anchored on the Authors line itself, with two guard
   assertions (the stripped text must start with the title; it must contain no section heading).
3. **The recorded output-file count is stale.** `output/` now holds 107 CSV + 6 R objects + 13 text + 2
   markdown = **128 files**, against the 122 recorded from the last clean-room run. Noted in report
   §2.4 and in `ANALYSIS_LOG.md` rather than quietly updated, because the 122 describes a verified
   execution and this does not.

---

## 10. Inventory of changes

**Scripts (6)**

| File | Change |
|---|---|
| `scripts/01_rna_load_qc.R` | New §3c: reference-mismatch design-alignment test → `rna_reference_sensitivity.csv` |
| `scripts/04_rna_markers_genesets.R` | New §2a: 20,000-draw joint permutation null, both units → `rna_marker_recovery_permutation.csv` |
| `scripts/12_rna_protein_concordance.R` | New §5c (per-gene p/q via the conventional Spearman t approximation, inverse-association tiers) and §5d (patient-representative sensitivity); `integ_rnaprot_cor.csv` gains `p_value`, `q_value` |
| `scripts/19_proteomics_dynamic_range.R` | §4b rewritten (√2 adjustment, lognormal CV, legacy column); new §4b-i (protein key built once and asserted unique); new §4b-ii like-for-like comparison with an ADC stratum on its own technical estimate; §4c decile join repaired and abundance-matched columns added |
| `scripts/34_fig2_qc.R` | Panel E series switched to per-measurement CV; panel F reference band switched to the cross-model *difference* span with the ratio recomputed; `frac_compressed` denominator fixed; console summary updated |
| `docs/manuscript/v4/build_docx.sh` | Front-matter strip re-anchored with guard assertions |

**New deposited tables (5)**

`rna_reference_sensitivity.csv` (2 rows) · `rna_marker_recovery_permutation.csv` (2) ·
`prot_noise_vs_biology.csv` (3) · `integ_rnaprot_negative_genes.csv` (4) ·
`integ_rnaprot_patientrep_sensitivity.csv` (3)

**Modified output schemas (4)**

`prot_bridge_agreement.csv` (+`sd_single`, `cv_pct_per_measurement`, `sd_diff_cv_pct_legacy`,
`loa_span`, `units_note`, `replicate_caveat`) · `prot_cv_by_abundance.csv` (+`bridge_sd_single`,
`bridge_cv_pct_per_measurement`, `n_bridge_diffs`, `prot_sd_median_decile`,
`observed_to_technical_sd_ratio`, `technical_variance_share_pct`, `approx_reliability_ratio`) ·
`integ_rnaprot_cor.csv` (+`p_value`, `q_value`) · `integ_rnaprot_negative_genes.csv` (+`test`, rewritten
`interpretation`, tiers renamed)

**Regenerated figure**

`fig2` in plain mode → `reports/assets2/fig2.{pdf,png,svg}`; in manuscript mode →
`docs/manuscript/figures/fig2.{pdf,png}`. Panel F's limits-of-agreement lines now render *inside* the
reference band, which is the visual form of the corrected conclusion.

**Documents (9)**

`reports/01_multiomic_characterization_results.{md,html,pdf}` (36 → 38 pp; ~20 sections edited, new
revision callout, limitations 23 → 25, open items 10 → 11) ·
`docs/manuscript/v4/OvCAN_data_descriptor_v4.{md,docx}` (abstract 176 → 164 words, 7 new sections, gap
register 15 → 22) · `docs/manuscript/v4/{LEDGER,CHECKLIST,TRIAGE,VENUES_society_data}.md` ·
`reports/02_manuscript_outline.md` and `reports/05_scientific_data_descriptor_draft.md` (both marked
**SUPERSEDED** with inline `[CORRECTED 2026-07-27]` annotations rather than edited wholesale, since v4
is the live draft) · `ANALYSIS_LOG.md`

---

## 11. What remains open

**Not done, and deliberately so.**

- **A full clean-room re-execution.** Six scripts changed across two rounds and five outputs were
  added. The 122-file, 25-numeric-check verification predates these edits and should be rerun before
  submission. Flagged in report §2.4 and `ANALYSIS_LOG.md`; not claimed as done. **This is now the
  first thing to do**, since round 2 confirmed the analytical content is settled — the round-2 review's
  own closing recommendation was to make the language changes, fix the protein key, and then rerun.
- **Relocating Fig. 6B to the Supplementary Information** (gap 22) — requires re-cutting a
  pixel-verified composite figure; flagged for the PI rather than done unprompted.

**Unchanged by this response — the review's Critical blockers are external requests.** None of
blockers 1–6 can be closed by analysis: repository deposition and reviewer access; TMT acquisition and
search parameters; STR and mycoplasma results; ethics approvals and the Human Data Checklist; culture
metadata; exome capture kit and target BED. All are represented in the v4 gap register, which now also
carries the four new administrative gaps (author contributions, competing interests, funding, AI
disclosure) and the two new scientific ones (RNA re-quantification, data dictionaries and checksums).

**The one Major-before-acceptance item that analysis could not close** is blocker 9 — `TOV21G` remains
a candidate. MSI-PCR and MMR immunohistochemistry are the only route to attribution, and the manuscript
now says so in place of implying the signature refit settles it.

---

## 12. Round-2 record

Six refinements, all accepted. Two changed deposited numbers; one was a defect.

| # | Round-2 finding | Disposition | Effect on numbers |
|---|---|---|---|
| 1 | `reliability_icc` is not an ICC — it mixes a pooled-across-proteins bridge error with a median-across-genes observed SD, from four different samples | **Accepted.** Renamed `approx_reliability_ratio`; four approximations enumerated in the script and both documents; phrasing changed to "under a homoscedastic independent-error approximation, bridge variability corresponds to roughly 27–44% of the variance of a typical protein's observed cross-model measurements" | Same global values, no longer called a reliability. **New abundance-matched table** replaces the single figure: technical share **67.5% → 7.5%** across deciles, variance ratio **0.33 → 0.92** |
| 2 | The ADC statement was not per-class — it swapped the observed SD but kept the global technical SD | **Accepted.** Bridge differences restricted to the ADC target rows give a direct estimate | Technical SD **0.134** from **32 differences**, observed 0.600 → technical share **≈5%**, revised from the round-1 "6–10%" |
| 3 | The reference analysis is a risk diagnostic, not a matched-index sensitivity analysis; the association is equally consistent with genuine biology | **Accepted.** Relabelled "design-alignment diagnostic" everywhere; the round-1 sentence "I implemented it and it did the opposite" is retracted in §4; matched-index re-quantification retained as the actual, undone sensitivity analysis | No numbers change; framing corrected in report §2.1, limitation 6, §12 item 5, manuscript Methods, and §4 of this changelog |
| 4 | Fisher z on Spearman is a non-standard approximation; and "RNA fails as a proxy" is the wrong description of an inverse association | **Accepted on both counts.** Switched to the conventional asymptotic Spearman t approximation (df = n − 2); tiers renamed for the association they establish | Nominal **51 → 50**, FDR **30 → 29**, FDR-positive **3,626 → 3,599**. Terminology: "FDR-supported inverse association", never "fails as a proxy" |
| 5 | The commonality disagreement is rhetorical — the review called the analysis useful; and donor collapse does not rule out a centre-aligned batch effect | **Accepted.** §6 recategorised from "not accepted" to "accepted, boundary clarified"; report §4.2 and the summary callout now state that neither exclusion addresses the confounding | No numbers change |
| 6 | The abundance-decile bridge join keys decorated matrix row ids against raw gene symbols — 13 symbols / 44 rows affected | **Accepted; this was a defect.** `(Symbol, Uniprot)` is unique on both sides and is now the join key, resolved to the canonical row id before the decile join, with the join asserted one-to-one and the `relationship = "many-to-many"` escape hatch removed | Decile SD of differences **0.430 → 0.429** (decile 1) and **0.106 → 0.105** (decile 10); Fig. 2E endpoints unchanged at **21.3% → 5.2%**. As predicted, no conclusion moves — but it was wrong |

**Where round 2 independently reproduced the work**, its checks are recorded here because they are
evidence and not commentary: the conventional Spearman test gives 50/29 against the Fisher-z 51/30, and
recomputing the donor-collapse comparison on the identical 7,832-gene set gives per-gene median
0.3976 → 0.4168 with the inverse fraction 0.1023 → 0.0970 — confirming that the reassuring
patient-representative result in §2.5 is not an artefact of the slightly different gene sets.

**Round-2 code changes.** `scripts/19_proteomics_dynamic_range.R` — new §4b-i (protein key built once,
asserted unique), §4b-ii renamed and re-caveated with an ADC stratum on its own technical estimate, §4c
join repaired and abundance-matched columns added. `scripts/12_rna_protein_concordance.R` — §5c test and
tier names changed. `prot_noise_vs_biology.csv` gains `stratum`, `technical_sd_basis` and
`approximation_note` and loses `reliability_icc`; `prot_cv_by_abundance.csv` gains `n_bridge_diffs`,
`prot_sd_median_decile`, `observed_to_technical_sd_ratio`, `technical_variance_share_pct` and
`approx_reliability_ratio`; `integ_rnaprot_negative_genes.csv` gains `test` and a rewritten
`interpretation`. Fig. 2 regenerated in both modes; report and `.docx` rebuilt.

---

## 13. Bottom line

Round 1's two substantive statistical findings were correct, and one inverted a headline claim: the
protein layer's technical error is a **minority share** of observed cross-model variance, not a quantity
exceeding biological spread. The marker-recovery conclusion survived a stricter null and strengthened.
The two factual errors were real, and one had a systemic cause in the manuscript's own planning
documents that has now been closed.

Round 2 was mostly about the difference between a defensible statistic and a defensible *claim about* a
statistic, and it was right on every count. Three of its findings narrowed language that was still
running ahead of the evidence — an approximate variance ratio called a reliability, a risk diagnostic
called a sensitivity analysis, an inverse association called a failure of proxy value. One found a real
defect in a join. One made a comparison genuinely per-class rather than nominally so. And one caught me
mischaracterising the review's own position in order to record a disagreement that did not exist.

The most useful outcome of round 2 is not a correction but a replacement: dropping the homoscedasticity
assumption turned one global noise figure into an **abundance-conditional gradient** (technical share
67.5% → 7.5% across deciles). That is a better answer to the question a reuser actually asks, and it
came from taking the objection seriously rather than hedging the sentence.

The resource's substance is unchanged across both rounds. What changed is the set of claims it licenses.
The remaining barriers to submission are the clean-room rerun and the external acquisition,
authentication, ethics, culture, capture and repository gaps — none of which is an analytical problem.
