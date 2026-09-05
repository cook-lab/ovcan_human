# The OvCAN multi-omic cell-line resource *— consolidated analysis report*

**Lede:** Matched transcriptome, proteome and exome data for 42 human ovarian cancer cell models from 34 patients across seven histotypes, generated at three Canadian centres and processed through one pipeline. This report is the internal reference for everything the re-analysis has examined: what the resource contains, how each layer was made, what each layer can and cannot support, which biology it recovers, what it says about model identity, and how to reuse it. Every number is traceable to a file in `output/`.

**Prepared** 2026-07-26 · **Revised** 2026-07-27 in response to the methodological validity review (`reports/06_methodological_validity_and_submission_readiness_review.md`) · **Analyst** Cook Lab (Claude) · **Scripts** `scripts/00–22` (analysis) + `30–37` (figures) · **Source of truth** `metadata/samples.csv`, `metadata/line_family_map.csv`, `output/*` · **Numeric authority** `docs/manuscript/v2/EVIDENCE_DOSSIER.md` · **Full trail** `ANALYSIS_LOG.md`

{{callout: What the 2026-07-27 revision changed}}
Five statistical corrections, each traceable to a rerun script rather than a wording change. **(1)** The protein bridge "repeatability CV" of 15.7–20.4% lacked the √2 replicate-difference adjustment; the per-measurement figure is **10.4–13.2%**, and the claim that technical noise exceeded biological spread by 2.4–3.1× was an artefact of comparing a 3.92-SD interval width with a 1.35-SD IQR — **the corrected like-for-like comparison reverses its direction** (§3.3). **(2)** Marker recovery is now tested against a **correlation-preserving joint label permutation** rather than an exact binomial that assumed the 25 markers independent, and the patient-representative unit is primary: 15 of 25, permutation *p* = 0.0046 (§5.3). **(3)** The 3.34× protein-versus-RNA spread ratio is reported as **consistent with** isobaric compression rather than as an isolated compression factor setting a concordance ceiling (§6.2). **(4)** The 813 inverse transcript–protein correlations are **candidates**; 29 carry an FDR-supported inverse association, which is not the same as transcript failing as a proxy (§6.1). **(5)** Histotype-versus-centre language is now associational throughout, because the design has no cell in which the two vary independently (§4.2, limitation 5b). Two factual corrections: `BIN67` was not reclassified by Karnezis 2021 and is no longer counted as such (§9.3), and the reference-mismatch loss is not evenly distributed with respect to the design (§2.1, new open item 5).

A second review round then tightened four of these. The protein variance ratio is **not an intraclass correlation** and is now named and caveated as an approximate diagnostic, with an **abundance-matched** version replacing the single global figure (technical share 67.5% in the lowest abundance decile to 7.5% in the highest) and the ADC-target comparison recomputed from **those targets' own** bridge differences (§3.3). The reference-mismatch analysis is a **design-alignment diagnostic, not a sensitivity analysis** — the association it finds is equally consistent with genuine biology, and the matched-index comparison remains undone (§2.1, §12 item 5). Spearman inference moved to the conventional t approximation (50 nominal / 29 FDR, from 51 / 30). And a protein-key join defect in the abundance-decile stratification was repaired.
{{endcallout}}

{{stats}}
42 / 34 || cell models / independent patients — RNA 31 · protein 31 · exome 23
100% || *TP53* in HGSC — 11 / 11 patients, every call Tier 1–2
0.40 || transcript–protein Spearman — per-gene median over 7,894 genes
5 / 5 || DepMap self-matches at rank 1 of 67, all reciprocal-best
{{endstats}}

{{callout: What this resource establishes, and where it stops}}
- **The three layers meet or exceed the quality benchmarks conventional for their assay**, and the numbers that bound each layer are measured rather than asserted: RNA pseudoalignment 91.1% median, per-measurement protein repeatability SD 0.149–0.190 log2 (10.4–13.2% CV), exome filtering from 557,392 raw records to 6,036 coding candidates.
- **Histotype labels are associated with the leading transcriptomic axes, and centre adds no independent explanatory signal.** On a joint model of RNA PC1 the unique-histotype component is 0.424 raw / 0.393 adjusted against 0.002 / −0.025 for unique centre, and the adjusted unique-centre component is negative on all five leading components. **What this cannot do is separate histotype biology from a site-, culture- or processing-linked effect aligned with histotype**, because all 15 high-grade serous RNA models come from one centre and the shared variance component is 0.311 (§4.2). Separately, collapsing sublines to one model per patient *strengthens* the structure, which excludes duplicate donor genomes as its cause — but not a centre-aligned batch effect, since the representatives come from the same centres in the same proportions.
- **The exome recovers the genetics each histotype requires.** *TP53* is mutated in 17 of 17 high-grade serous models and 11 of 11 high-grade serous patients; fraction of genome altered runs 0.635 in HGSC to 0.021 in the single low-grade serous model, a 30-fold panel range in the predicted direction on n = 1 outside HGSC.
- **Every model with an external reference matches it,** and the three layers are molecularly consistent with one published pathology reclassification, `TOV112D` → SWI/SNF-null dedifferentiated carcinoma (Karnezis 2021), and with the established SCCOHT identity of `COV434` and `BIN67`. `BIN67` is the case that justifies measuring more than one layer: its *SMARCA4* transcript sits mid-panel while its *SMARCA4* protein is second-lowest.
- **One model is a genuine reuse feature.** `TOV21G` carries 1,416 coding candidates (6.9× the panel median), an indel-rich spectrum, and the only bootstrap-supported mismatch-repair signature exposure in the panel. Candidate MMR-deficient clear-cell models are rare.
- **Four limits change what analyses are possible, not merely how they are worded:** bridge variability corresponds to roughly 27–44% of a typical protein's observed cross-model variance and, once matched on abundance, from 68% in the lowest decile to 8% in the highest — so whether one protein in one model can be interpreted depends on that protein's abundance; copy number is total rather than allele-specific, so no allele-specific quantity or genomic HRD score is derivable; variant calls are tumour-only; and the layers are not passage-matched. Each is quantified in §11.
{{endcallout}}

---

## 1 · Composition of the resource

### 1.1 Models, patients, histotypes and per-assay coverage

The panel is the 42 rows of `metadata/samples.csv` carrying `analysis_include = Y` **and** `provenance = generated`. The file has 55 rows; the other 13 are published low-grade serous data generated by another group (`source_site = Carey`) that this project did not produce. A naive `nrow()` returns 55 and is wrong.

- **42 cell models from 34 independent patients**, spanning **seven histotypes** and **three contributing centres**.
- **Thirteen models are sublines of five patients**: families 3133 (4 models), 2295 (3), 1369 (2), 3121 (2), 3291 (2). All five families are high-grade serous and all come from one centre, so subline duplication is confounded with both histotype and centre.
- **The subline structure is confirmed at the variant level, not assumed**: all four members of family 3133 carry *TP53* p.Q192\*, all three of family 2295 carry p.I195T, and both members of 1369 carry p.G244C. That is direct internal evidence of common donor origin.

| Histotype | Models | Patients | RNA | Protein | CNV | Mutations | All three |
|---|---|---|---|---|---|---|---|
| HGS (high-grade serous) | 24 | 16 | 15 | 15 | 18 | 17 | 9 |
| CC (clear cell) | 8 | 8 | 7 | 7 | 2 | 2 | 2 |
| MC (mucinous) | 3 | 3 | 3 | 3 | 1 | 1 | 1 |
| EC (endometrioid) | 2 | 2 | 2 | 2 | 1 | 1 | 1 |
| MMMT (carcinosarcoma) | 2 | 2 | 2 | 2 | 0 | 0 | 0 |
| SCCOHT (small cell, hypercalcaemic) | 2 | 2 | 2 | 2 | 0 | 0 | 0 |
| LGS (low-grade serous) | 1 | 1 | **0** | **0** | 1 | 1 | 0 |
| **Total, models** | **42** | **34** | **31** | **31** | **23** | **22** | **13** |
| **Total, patients** | — | — | **28** | **28** | **16** | **16** | **10** |

Two consequences are easy to get wrong. First, **the single low-grade serous model (`TOV81D`) has exome data only**, so every expression analysis in this report has **six** histotype groups, not seven; the level vector used in code and figures is `c("HGS","CC","EC","MC","MMMT","SCCOHT")`. Second, **the fully multi-omic subset is smaller than any single layer, and it is not a random subset of the panel**. All 13 models come from CHUM and 9 of 13 are HGSC, and because 6 of the 13 belong to only 3 patients the subset represents **10 patients**, not 13. Anything computed on it describes a single centre and predominantly one histotype.

{{figure: fig1 | **Resource overview.** (**A**) Workflow schematic — box fill encodes stage type only; sizes and arrow lengths are layout, not quantities. (**B**) Assay coverage, one row per model, sorted by histotype then patient family. Left strips give histotype and family; the four assay columns use a two-level fill (ink = present, pale = absent); the `All 3` column marks the 13 models with RNA, protein and exome; column totals are printed beneath. `TOV3121D` is the only model differing between the two exome columns — it has copy-number calls but no Mutect2 MAF. (**C**) Bar length = models per histotype; diamond on the same axis = distinct patients of origin.}}

### 1.2 Ten denominators coexist, and different analyses use different ones

Ten legitimate denominators coexist. A statistic quoted without its unit will be wrong somewhere in this report, because the correct unit differs between adjacent paragraphs.

| Unit | n | Where it is the denominator |
|---|---|---|
| Models in the resource | 42 | Composition, coverage, per-model tables, STR status |
| Independent patients | 34 | Any frequency that must not double-count a donor |
| RNA models · protein models | 31 · 31 | PCA, silhouettes, differential expression, gene-set scores, variance components |
| Models with **both** expression layers | **30** | Every transcript-versus-protein statistic |
| Patient representatives **with RNA** | **28** | **The inferential unit**: marker effect sizes, bootstrap intervals, the marker-recovery permutation, and the DE/pathway sensitivity analyses |
| Exome copy-number models · patients | 23 · 16 | Fraction of genome altered, segment profiles |
| Exome mutation models · patients | 22 · 16 | Driver calls, burden, signatures |
| HGSC exome models · patients | 18 · 11 | **Arm-level copy-number frequencies** (not the whole exome panel) |
| Fully multi-omic models · patients | 13 · 10 | Cross-layer passage comparison |

`scripts/15_patient_family_map.R` flags one model per patient by ordering deterministically on `omics_score = has_rna + has_prot + has_wes_maf + 0.1 × has_wes_cnv`, giving **34 representatives**, of which **28 carry RNA**.

**The division of labour between 31 and 28 is deliberate and follows the kind of claim being made.** Treating all 31 line models as independent observations overstates the precision of every one-versus-rest comparison, and it does so even in groups that contain no sublines, because duplicate donors also enter the comparator set — §4.7 measures that directly. So **anything inferential runs on the 28**: marker effect sizes and bootstrap intervals, the marker-recovery permutation test, and the differential-expression and enrichment sensitivity analyses. **Anything descriptive runs on the 31**: PCA displays, silhouettes, marker heatmaps, the per-model browser and the model-selection tables, where the line rather than the donor is the object a reuser will order. Both are deposited for every analysis where both are computable, and where they disagree the 28-model result governs the claim. With only three replicated RNA donor families a donor-aware mixed model is the other defensible option, but it would not estimate donor variance precisely, which is why collapse rather than modelling is used here.

### 1.3 The exome layer is single-centre; the expression layers span three centres

| Centre | Models | Patients | RNA | Protein | CNV | Mutations | Tri-omic |
|---|---|---|---|---|---|---|---|
| CHUM (`Mes-Masson`) | 29 | 21 | 19 | 19 | **23** | **22** | **13** |
| BC Cancer / OVCARE (`Huntsman`) | 12 | 12 | 11 | 11 | 0 | 0 | 0 |
| OHRI (`Huntsman/Vanderhyden`) | 1 | 1 | 1 | 1 | 0 | 0 | 0 |

**The entire exome layer is single-centre.** All 23 copy-number and all 22 mutation-call models come from CHUM, so no centre comparison exists in the genomic layer, and no genomic result in this resource can be separated from CHUM-specific culture or processing history. In the expression layers a comparison is possible but weak: the only histotype present at two centres in usable numbers is clear cell, at 2 models versus 5 (§4.5).

### 1.4 Three models carry unresolved histotype annotation conflicts

Three of 42 models have a `subtype_status` of CONFLICT rather than `consensus`:

- **`OV3331` and `OV90`** are labelled "Adenocarcinoma" in some source annotations and HGS in others. Both were assigned HGS, and the genomic layer supports that: *TP53* mutation plus a high fraction of genome altered.
- **`TOV112D`** is labelled endometrioid but was reassigned to dedifferentiated carcinoma with SWI/SNF loss (Karnezis 2021). The resource's own three layers are molecularly consistent with that reassignment on evidence generated independently of it (§9.3). It matters twice over: `TOV112D` is 1 of only 2 endometrioid models and 1 of the 13 tri-omic models, so the reassignment leaves endometrioid at effectively **n = 1** for marker validation.

---

## 2 · Data generation and processing

### 2.1 Transcriptome

All RNA libraries were sequenced in a **single run** (`NS.1676.003`), so the RNA layer carries no sequencing-batch term. That removes one candidate explanation for structure before any modelling.

- **kallisto 0.46.0** against an index of **185,299 targets** (`index_version 10`, built pre-2021-10-20 against an Ensembl-104-era annotation), summarised to genes with **tximport 1.36.1** (`ignoreTxVersion = TRUE`, `countsFromAbundance = "no"`), so gene-level counts are summed transcript estimates and are **not** length-corrected.
- The transcript-to-gene map is pinned to the **Ensembl 105 archive** (266,615 transcripts; md5 `119dfe0ab4f856e6b81efc6ce78f4ba7`, asserted at load) and lives at `data/reference/` as a protected input rather than a regenerable output.
- **That pin has a measured cost.** The index predates the map, so **3,529 of 185,299 index targets (1.90%)** are absent from the release-105 map and are dropped by tximport without warning, carrying a median **2.22%** of TPM per model (1.60–3.41%) and **1.53%** of estimated counts. A reuser who rebuilds the index against Ensembl 105 will therefore *not* reproduce these matrices, which is why the map must be deposited with the data.
- **The loss is not evenly distributed across the panel, and that is worth knowing before the structure analyses are read.** Regressing each model's dropped-TPM fraction on the design terms (`output/rna_reference_sensitivity.csv`) gives R² = **0.330** for contributing centre (Kruskal–Wallis p = 0.0033) and **0.388** for histotype (p = 0.010) — the two terms are collinear here, so these are not independent findings.

  **What this is: a design-alignment diagnostic. What it is not: a sensitivity analysis.** It establishes that the omitted transcripts are not distributed evenly with respect to the terms §4 interprets, which is a reason to be cautious. It does **not** show that the mismatch manufactured any observed structure, and it cannot: **biologically distinct histotypes would be expected to express different amounts of any transcript set, including this one**, so the association is exactly what a purely biological model predicts as well. The diagnostic flags a risk; it does not measure a harm, and it cannot be read in either direction as evidence about §4's conclusions.

  Three things bound the risk it flags: the absolute quantity is small and tightly bounded (1.60–3.41% of TPM, a 1.81-percentage-point range across all 31 models), the reference pair is identical for every model so none is advantaged by a different annotation, and the affected transcripts are absent from the gene-level matrix entirely rather than mis-assigned to the wrong gene. **The actual sensitivity analysis remains undone**: only re-quantifying against an index matched to the deposited transcript-to-gene map can show whether the QC, PCA, marker and pathway results are invariant to the mismatch. That is item 5 in §12.
- **39,568 genes** returned by tximport → **22,544** retained (≥10 estimated counts in ≥2 models) → variance-stabilised with **DESeq2 1.48.2**. **22,542** genes enter the per-gene variance decomposition; two failed to converge and their identifiers are deposited. The 2-gene difference is real, not a typo.

Two depth quantities are named separately in the deposited files because they are routinely confused. `assigned_gene_counts` is the column sum of the post-filter estimated-counts matrix (median 56,497,657); `n_processed_fragments` is kallisto's sequenced depth (median **64.3 M paired-end fragments**). Neither is a read count, and fragments are not reads.

### 2.2 Proteome

Protein abundance was measured by isobaric tandem mass tag labelling in **five 11-plex sets** (55 channels). Channel 1 of every set carries a **pooled internal standard** to which all reported abundances are ratioed, and channel 9 is empty. That leaves 45 channels carrying data, with the 31 proteomics models distributed **7 / 6 / 5 / 6 / 7** across sets 1–5.

- The inter-set links form a **chain** (1↔2↔3↔4↔5), and the four links carry **four different samples** (`VOA10816`, `TOV1369`, `VOA3993`, `TOV3133G`), not four replicates of one sample. One of them, `VOA3993`, is a Carey low-grade serous proteomics sample **external to the 42-model panel**, and its identity to any RNA sample is recorded as unconfirmed. This must be disclosed wherever bridge reproducibility is reported.
- The search returned **8,430 rows**; three carry a UniProt accession but no gene symbol and cannot enter a symbol-keyed analysis, leaving an **8,427-row analysis matrix**. Identifications are filtered at q ≤ 0.01 (median q = 0), with a median of 12 total and 11 unique peptides per protein and **146,830 peptides** quantified.
- **The proteomics processing workflow is not recoverable.** The vendor readme accompanying the archived data is a 1-byte empty file, so instrument, acquisition mode, search engine, database, FDR thresholds, rollup and the definition of the vendor-reported CV are all unknown. Acquisition mode (MS2 versus SPS-MS3) bears directly on the magnitude of isobaric ratio compression, which this report quantifies in §6.2, so that interpretation currently rests on an unknown.

### 2.3 Exome

- **nf-core/Sarek**, variants called with **GATK Mutect2 in tumour-only mode**. No matched normal was sequenced: `n_depth`, `n_ref_count` and `n_alt_count` are empty in every MAF row.
- The build is **GRCh38/hg38**, established from the Mutect2 panel of normals, VCF contig lengths (chr1 = 248,956,422) and hotspot coordinates. The MAF `NCBI_Build` field reads `GRCh37`, which is a spurious vcf2maf default; the deposited derived table carries no such column, and the string survives only in the raw vendor MAFs under `judy_archive/`, which were deliberately not modified.
- **Filtering is three-stage**: `FILTER == PASS`; population allele frequency ≤ 0.001 (max of `gnomADe_AF`, `AF`, `AA_AF`, `EA_AF`, NA treated as 0); then membership of the maftools non-synonymous class set.
- `HGVSp_Short` was empty in every row, so **all 6,036 protein-change strings were reconstructed** from position and alleles: 5,682 flagged canonical, 103 non-canonical, 251 not derivable. Only strings flagged `hgvsp_canonical` should be quoted. Frameshifts are written `fs` and nonsense `*`; an earlier convention wrote `X` for frameshifts, which matches nothing in ClinVar or COSMIC.
- Copy number: **CNVkit 0.9.10** against the GRCh38 GATK bundle with `--drop-low-coverage --diagram --scatter`, referenced to a panel of **five public healthy exomes (BioProject PRJNA339046) unmatched to these tumours**. Arm-level calls use |log2 copy ratio| > **0.20** over > **50%** of arm length after per-sample autosome probe-weighted median centring.
- **The exome capture kit is not recoverable.** The target file `intervals_sorted.bed` is referenced by every CNVkit invocation but is not archived. This is a validity question rather than only a documentation gap, because it leaves unresolved whether the five pooled normals, which come from an unrelated study, share the tumours' capture design. It is also the stated reason no exome-to-genome trinucleotide renormalisation was applied to the mutational signatures.

### 2.4 Reproducibility

- **R 4.5.2** / **Bioconductor 3.21** on aarch64-apple-darwin20, seed **1234**, re-set inside stochastic functions. The environment is pinned in **`renv.lock` (288 packages)** and restored with `renv::restore()`; `output/package_versions.csv` remains the citable per-package record, and should be preferred to `session_info.txt`, which records only the last script sourced.
- **`scripts/run_all.sh`** encodes the real dependency order in six stages, because numeric script order is *not* dependency order. Two places where the numbering misleads are called out in its header: 18 must precede 10, and 15 must precede 30.
- **A clean-room re-execution from the archived inputs regenerated all 122 output files** across 23 script invocations with zero non-zero exits: 109 byte-identical, 12 `session_info*` files differing only in an embedded timestamp, 1 (a bootstrap cache) legitimately larger because the reference had been written by a superseded code path. **All 25 headline numeric checks re-verified.** Total sequential runtime 1,775 s (~30 min), so a full clean run is a practical pre-submission check.
- Two classes of latent defect were found and closed by that exercise. Three dependencies were satisfied only because a file happened to exist on disk: the filtered-variant table, whose write had been commented out; the patient/family map, read bare by several scripts; and the transcript-to-gene cache. Separately, a bootstrap cache was reusable while incomplete. All are now generated or guarded by assertions. These are the defects that let a pipeline look reproducible until someone runs it clean.
- `output/` holds 107 CSV files, 6 R objects, 13 text files and 2 markdown files. The clean-room figure of 122 above predates the 2026-07-27 audit revision, which added five deposited tables (`rna_reference_sensitivity.csv`, `rna_marker_recovery_permutation.csv`, `prot_noise_vs_biology.csv`, `integ_rnaprot_negative_genes.csv`, `integ_rnaprot_patientrep_sensitivity.csv`); a fresh clean-room run should be executed before submission and its file count re-recorded. `output/external/` separately holds the 652 MB external reference set (55 files: DepMap 24Q4 and Cellosaurus), which is **fetched rather than generated** and is not regenerable offline.

---

## 3 · Data quality by layer

Each layer is reported with both the quantity that flatters it and the quantity that bounds it. Exome quality is inseparable from the filtering cascade and is therefore reported in §7.1.

{{figure: fig2 | **Data quality across three assays.** (**A**) Genes detected against pseudoalignment rate, n = 31 models; **both axes are truncated**. Fill and shape both encode contributing centre; hollow marks are per-centre medians. (**B**) Genes detected against processed fragments; band is an OLS fit ± 95% CI and is descriptive only. (**C**) Per-gene cross-model interquartile range for RNA and protein over the 7,896 paired genes; the two layers sit on different normalisation scales — RNA log2(TPM + 1) with a hard zero floor, protein log2 abundance relative to the plex pooled internal standard — so this compares each layer's spread on its own scale. (**D**) Protein features by number of isobaric sets in which they are quantified, summing to the 8,427-row matrix; the bar at 0 is the 70 zero-set features. (**E**) Technical coefficient of variation against mean-abundance decile (~836 features per decile, denominator 8,357); ink = vendor-reported replicate CV with q25–q75 ribbon, rust = inter-set bridge repeatability **per measurement**, i.e. after the √2 replicate-difference adjustment and the lognormal conversion of §3.3, so both series are on the same scale. (**F**) Bland–Altman agreement for the four inter-set bridge links — **four different samples, not four replicates of one sample**; solid line = bias, dashed = 95% limits of agreement. The shaded band is the **like-for-like** reference: the 95% span a difference between two randomly chosen models would occupy for a typical protein (1.59 log2 = 3.92 × √2 × the median cross-model protein SD of 0.287), so both quantities are 95% intervals of a difference and the technical span (0.83–1.05) is **0.52–0.66×** it. The band in earlier versions was a cross-model IQR drawn on a paired-difference axis, which is the dispersion-measure mismatch §3.3 retires. y clipped at ±1.6 log2, hiding 0.19% of points.}}

### 3.1 Transcriptome quality is uniform across the 31 models, with no outlier

| Metric | Median | Range | Unit |
|---|---|---|---|
| Pseudoalignment rate | **91.1** | 85.8–93.1 | % of fragments |
| Genes detected | **20,119** | 18,775–21,650 | genes per model |
| Sequenced depth | **64.3 M** | 45.1–97.0 M | paired-end fragments |
| Assigned gene counts | 56,497,657 | 41.0–86.0 M | estimated counts |

Named extremes (Fig. 2A, 2B): lowest pseudoalignment `VOA4841` (85.8%), highest `TOV1369` (93.1%); fewest genes `TOV3392D` (18,775), most `VOA8762` (21,650); shallowest `TOV3291G` (45.1 M fragments), deepest `OV2295-R2` (97.0 M).

**No model is an outlier on any of the three metrics**, so there is no candidate for exclusion on technical grounds and the analysis set is the full 31.

### 3.2 BC Cancer models detect 4.8% more genes than CHUM models; the mechanism is unestablished

- BC Cancer models detect **+948 more genes (+4.8%)** (Fig. 2A) than CHUM models (medians 20,765 versus 19,817; exact Wilcoxon, n = 11 versus 19, **p = 0.0082**) at **4.1 percentage points lower** pseudoalignment (88.1% versus 92.2%, p = 1.1 × 10⁻⁵).
- Depth is the obvious candidate mechanism, and both depth differences point the way a depth explanation requires (+17.7 M fragments, +13.9 M assigned counts), but **neither reaches significance** (p = 0.372 and 0.767).
- Panel-wide, pseudoalignment and gene detection are **negatively** correlated (Pearson r = −0.635, p = 1.3 × 10⁻⁴) while depth and detection are positively correlated (r = +0.504, p = 0.0038).

**The detection difference is established and its explanation is not**, so it would be wrong to write that deeper sequencing accounts for it. What matters for reuse is that the difference is confounded with centre: any depth-sensitive reanalysis, including detection thresholds, low-expression genes and presence/absence calls, must model centre. The centre ranges for pseudoalignment also barely overlap (CHUM minimum 89.3%, BC Cancer maximum 89.9%), so the alignment difference is close to a clean separation and cannot be dismissed as noise even though its mechanism is unknown.

### 3.3 Technical error accounts for 27–44% of the observed cross-model protein variance

The bridge channels measure the protein layer's technical noise directly, on the same material through the same workflow. **Two scales are involved and mixing them was the error in an earlier version of this section**, so the conversion is stated before the numbers.

Each bridge link is one sample re-run in the adjacent isobaric set, so `primary − bridge` is a **difference of two measurements**. If the two carry equal independent error variance, then Var(difference) = 2 × Var(single measurement), and a per-measurement SD requires the replicate-difference adjustment

`SD_single = SD_difference / √2`.

An earlier version of this report quoted `100 × (2^SD_difference − 1)` = 15.7–20.4% and called it a repeatability CV. That quantity applies no √2 adjustment and is not the standard lognormal conversion either, so it overstated per-measurement imprecision by about 1.5-fold. The corrected quantity uses both steps, `CV = 100 × √(exp((ln2 × SD_single)²) − 1)`.

| Quantity | Value |
|---|---|
| Bias across the four bridge links | −0.028 to +0.009 log2; negligible, with no systematic inter-set offset |
| SD of the paired differences | 0.211–0.268 log2 |
| **95% limits-of-agreement span** | **0.826–1.052 log2** (a width on the paired-difference scale) |
| **Per-measurement SD** | **0.149–0.190 log2** (median 0.171) |
| **Per-measurement CV** | **10.4–13.2%** |
| Vendor-reported replicate CV, median | **5.31%** (IQR 3.36–8.79%, n = 7,940 features) |
| Pearson r per link | 0.991–0.994 |
| *Retired:* `100 × (2^SD_difference − 1)` | *15.7–20.4% — deposited as `sd_diff_cv_pct_legacy` for traceability only* |

**Lead with agreement statistics, not correlation.** A Pearson r near 0.99 across a ~15-log2 dynamic range is close to unavoidable for any two measurements of the same proteome and discriminates nothing. The vendor CV remains two- to two-and-a-half-fold smaller than the bridge estimate, and **the bridge figure is still the honest one for cross-model comparisons**, because the vendor value's definition is among the unrecovered method parameters.

**Comparing the noise band with biological spread requires the same dispersion measure on both sides.** A 95% limits-of-agreement span is 3.92 × SD on the paired-difference scale; a cross-model IQR is 1.35 × SD on the single-measurement scale. Dividing the first by the second inflates the ratio about 2.9-fold for free, and that is exactly how the retired claim "technical noise exceeds a typical protein's biological spread by 2.4–3.1×" arose. **The claim was an artefact of the comparison, and the corrected comparison reverses its direction.**

| Like-for-like comparison | Value |
|---|---|
| Median cross-model protein SD, 7,896 paired genes | **0.287 log2** |
| Per-measurement technical SD | 0.149–0.190 log2 |
| **Observed cross-model SD ÷ technical SD** | **1.5–1.9×** |
| **Technical share of observed cross-model variance** | **27–44%** |
| Implied biological SD, by subtraction | 0.215–0.245 log2 |
| Approximate variance ratio (1 − technical share) | 0.56–0.73 |

**These are approximate variance-ratio diagnostics, not intraclass correlations,** and the last row is deliberately not called a reliability. Four approximations separate it from one. It assumes **homoscedastic error** across the abundance range, which is false and is corrected immediately below. It compares a technical SD **pooled across proteins** within a link against an observed SD computed **per gene and then medianed across genes**, so the two are not the same protein-level quantity. The four links are **four different samples**, one re-run each, not replicated aliquots of a common reference. And it assumes **independent additive error**, so the biological SD is reached by subtraction rather than estimated — the observed cross-model SD is not a pure biological quantity. The defensible reading of the table is: *under a homoscedastic independent-error approximation, bridge variability corresponds to roughly 27–44% of the variance of a typical protein's observed cross-model measurements.* It is not: *a single protein in a single model has a reliability of 0.56–0.73.*

**Dropping the homoscedasticity assumption changes the answer substantially, and this is the version to use.** Matching the technical estimate to the observed spread *within* abundance deciles — both sides decile-local — replaces one number with a nine-fold gradient (`output/prot_cv_by_abundance.csv`):

| Abundance decile | Technical SD (log2) | Observed cross-model SD (log2) | **Technical share of variance** | Approx. variance ratio |
|---|---:|---:|---:|---:|
| 1 (lowest) | 0.303 | 0.369 | **67.5%** | 0.33 |
| 2 | 0.280 | 0.356 | 62.1% | 0.38 |
| 3 | 0.211 | 0.317 | 44.4% | 0.56 |
| 5 | 0.165 | 0.273 | 36.5% | 0.63 |
| 7 | 0.119 | 0.257 | 21.3% | 0.79 |
| 10 (highest) | 0.075 | 0.272 | **7.5%** | 0.92 |

**There is no single answer to "can I interpret one protein in one model" — the answer is set by that protein's abundance.** In the lowest decile technical error accounts for about two-thirds of the observed cross-model variance and a single measurement carries little information; in the highest decile it accounts for under a tenth and a single measurement is informative. The global 27–44% figure sits in the middle of that gradient and understates the problem at low abundance while overstating it at high abundance. **Any protein-level claim must therefore carry its abundance context**, which is also visible in the two precision series of Fig. 2E: vendor CV falls monotonically from **11.0% to 3.0%** across the deciles and the per-measurement bridge CV in parallel from **21.3% to 5.2%** (SD of differences 0.429 → 0.105 log2).

**The ADC targets are a favourable stratum, and this now rests on their own technical estimate.** Restricting the bridge differences to the eight ADC target rows gives a technical SD of **0.134 log2** from **32 differences** against a median cross-model SD of **0.600 log2**, so technical error accounts for about **5%** of their cross-model variance. Because it is estimated from the targets themselves rather than from the global pool, this is a genuine per-class statement — but on few differences, so it should be quoted with its n, and the abundance-matched deciles above are the cross-check (these eight are abundant proteins, so the two agree). The typical protein is the harder case, not the ADC panel.

The corrected quantities are deposited in `output/prot_bridge_agreement.csv` (per link), `output/prot_noise_vs_biology.csv` (global and ADC strata, each with its approximation note) and `output/prot_cv_by_abundance.csv` (the abundance-matched comparison).

### 3.4 Protein missingness is a property of the isobaric set, not of the sample

{{figure: figs4 | **Protein presence patterns across isobaric sets.** Unit: protein feature. (**A**) 32 distinct presence patterns across the 5 sets; ink tile = quantified in that set, row label = feature count. The counts sum to 8,427; the bottom row is the 70 features quantified in no set. The 1,572 features absent from at least one whole set include those 70, leaving 1,502 partially observed. (**B**) The same rows on a log10 axis as labelled points rather than bars.}}

- **Every model within an isobaric set has exactly the same number of detected features**: 7,738 / 7,758 / 7,562 / 7,694 / 7,653 for sets 1–5. Presence is therefore a property of the **set**, not of the sample.
- Overall missingness is **8.796%**, and it is itself confounded with set membership (7.94% in set 2 to 10.27% in set 3).
- **1,572 features (18.7%) are absent from at least one whole set**; 6,855 (81.3%) are complete; **70 were identified but never quantified in any model** and are retained in the deposited matrix so the row count stays 8,427, flagged `zero_plex = TRUE`.

So **a protein "absent" in a model may simply be absent from that model's set**, and because set composition tracks histotype, a set-conditional presence/absence call can look like biology when it is design. The recommended reuse filter retains the **7,733** features quantified in at least half the models.

Six protein feature counts and three peptide totals legitimately coexist, which is the largest single source of arithmetic confusion in the project. The arithmetic closes: 8,430 search rows − 3 with no gene symbol = 8,427; 6,855 complete + 1,572 partial = 8,427; 8,427 − 70 zero-set = 8,357 (the abundance-decile denominator); 8,396 distinct symbols + 31 non-representative duplicate rows = 8,427. `output/README_protein_matrix.md` is the single pointer, and `prot_feature_accounting.csv` the reconciliation.

---

## 4 · Sources of variation in the expression data

All 15 high-grade serous RNA models come from one centre, so histotype and contributing centre are collinear by construction. If the apparent histotype separation were a centre effect, every downstream result would be affected, which makes this the section the rest of the report depends on. **The conclusion this section can and cannot reach is bounded by the design**: it can show that centre adds nothing beyond histotype, and it cannot show that the shared variance is biology rather than a site-aligned effect. §4.2 states the boundary explicitly, and every downstream biological claim in this report should be read as histotype-*associated* rather than histotype-*caused*.

### 4.1 Principal components and marginal R² for histotype and centre

| | RNA (22,544 genes) | Protein (6,855 complete-case features) |
|---|---|---|
| PC1 variance | **20.7%** | **12.7%** |
| Top 5 PCs | 52.3% | 44.1% |
| PC1 R², histotype vs centre (Fig. 3A, 3B) | **0.735** vs 0.313 | 0.553 vs 0.361 (adjusted 0.464 vs 0.315) |

**Histotype R² exceeds centre R² on every component in both layers.** Protein PC1 carries the largest centre R² anywhere in the analysis, which is why it is the value quoted.

The isobaric-set (plex) term should not be summarised by a single number. On protein PC1 its raw R² of 0.134 collapses to **0.0004** under adjustment, almost all of it having been the degrees of freedom a 5-level term buys, yet it survives adjustment on PC2 (0.215) and PC5 (0.229). A batch effect concentrated in a few high-variance directions behaves this way: low per-feature median, high loading on a couple of components. Both readings are true and neither alone is honest.

### 4.2 Commonality analysis: centre adds no independent signal, but the design cannot license a causal reading

Comparing two marginal R² values cannot answer a collinearity question, so the joint model was partitioned into unique and shared components instead.

| RNA PC | Unique histotype | Unique centre | Shared | Adjusted unique histotype | **Adjusted unique centre** |
|---|---|---|---|---|---|
| PC1 | 0.4238 | 0.00195 | 0.3112 | 0.3928 | **−0.0251** |
| PC2 | 0.4806 | 0.00051 | 0.2136 | 0.4438 | **−0.0312** |
| PC3 | 0.5348 | 0.00236 | 0.0468 | 0.4761 | **−0.0406** |

- **The adjusted unique-centre component is negative on all five leading components** (−0.025 to −0.062). A negative adjusted unique component means the centre term adds nothing beyond histotype once its degrees of freedom are charged for.
- Against **1,000-draw permutation nulls**, unique-histotype variance on PC1–PC3 (0.42–0.53) sits far above its null (means 0.11–0.21; 95th percentiles 0.24–0.39), at the minimum attainable p (1/1001 is the floor, so the correct form is **p < 0.002**). Unique-centre variance sits **below its null mean on every component** with permutation p of 0.49–0.86, i.e. indistinguishable from chance and, if anything, smaller than chance.
- `joint_aliased = FALSE` on all five components, so the two designs are collinear but not aliased, and the partition is estimable.

**Interpretation, and its boundary.** The marginal centre R² of 0.31 on PC1 was almost entirely variance that centre *shares* with histotype rather than independent signal, which is the failure mode that comparing marginal R² values produces. What the partition establishes is therefore precise and limited: **centre carries no explanatory signal that histotype does not already carry.**

It does **not** establish that the separation in Fig. 3A is biology rather than a site-linked effect, and this report should not be read as claiming that it does. The reason is the design, not the model. A **shared component of 0.311 on PC1 is large, and no analysis can apportion it**, because there is no design cell in which histotype and centre vary independently: all 15 high-grade serous RNA models come from CHUM, the rare histotypes are concentrated at BC Cancer, and the only within-histotype cross-centre comparison available is clear cell at 2 versus 5 models (§4.5, uninformative in both directions). A regression cannot estimate what a high-grade serous model grown at BC Cancer would look like when no such model exists. Tissue handling, growth medium, passage history, extraction and library preparation are all aligned with centre here, and §3.2 shows a real, unexplained centre difference in gene detection while §4.3 shows the unsupervised top-level split following centre.

**The defensible statement is that histotype labels are associated with the leading transcriptomic axes and retain unique explanatory signal under adjustment and under patient collapse, while source- and culture-associated effects aligned with histotype cannot be excluded.**

Two narrow alternatives *are* excluded, and it is worth being exact about what each one covers, because they are easy to over-read. This section excludes **centre as a source of variance beyond histotype** — not a centre-aligned batch effect, which by construction contributes to the shared component and is untouched by the partition. §4.7 excludes **duplicate subline genomes** as the cause of the structure — collapsing to one model per patient removes non-independent donors, which says nothing about a centre-aligned batch effect either, since the representatives still come from the same centres in the same proportions. **Neither result addresses the confounding, and neither should be cited as if it did.**

{{figure: fig3 | **The data recapitulate subtype biology.** (**A**) RNA principal components 1 and 2 by histotype, n = 31 models (one point per model, not per patient representative); `prcomp` on the top 2,000 most variable variance-stabilised genes, centred and unscaled. Six histotype levels are present, no low-grade serous model having RNA; the endometrioid models are labelled because `TOV112D` sits with the SCCOHT pair. (**B**) The same points by contributing centre. (**C**) Per-gene transcript–protein Spearman correlation across the 30 models with both layers — n = 7,894 genes, a gene receiving a correlation only if at least 10 models carry both measurements. The band at 0.38–0.48 and the line at 0.72 are external literature values, neither computed here; x clipped to [−0.5, 1]. (**D**) Per-feature variance decomposition by layer; filled circle = median, open diamond = mean, bar = interquartile range. Percentages need not sum to 100, and the RNA facet has no batch term, so the two facets are not the same model. `Patient` is asterisked as a design artefact — 28 levels on 31 observations with 3 replicated patients. (**E**) Canonical marker recovery, 25 graded markers plus an ungraded proliferation control, row-scaled: **for each marker, z is computed across the six histotype means of that marker**, not across models. Expected cells are boxed, solid = high and dashed = low. (**F**) The same rows with all 31 RNA models as points on an absolute axis; the dashed line at log2 TPM = 1 is the rank rule's expression floor.}}

### 4.3 Unsupervised clustering of all 22,544 genes splits by contributing centre

The reassuring 0.2% unique-centre figure of §4.2 should be read alongside the design it comes from. Clustering the 31 × 31 matrix of between-model Spearman correlations on **all 22,544** variance-stabilised genes, with no feature selection, splits at the top level **by contributing centre**.

{{figure: figs3 | **Model-to-model RNA correlation and the centre-wise split.** Clustered heatmap of the 31 × 31 between-model Spearman correlation matrix on all 22,544 variance-stabilised genes (not the top-2,000 set used for the principal components); Ward.D2 linkage. The colour scale is **not** truncated, spanning the observed off-diagonal range 0.757–0.969 (median 0.834). The top-level split is essentially by contributing centre — the collinearity of Fig. 3B seen without modelling. Of the three multi-line families with two or more RNA models, only 3133 co-clusters as an adjacent pair.}}

Two things temper this. The off-diagonal correlations span only **0.757–0.969** with a median of 0.834, so all 31 models are highly correlated and the centre-wise block structure is a modest modulation of a high baseline, not two separate populations. And of the three multi-line families with ≥2 RNA models, **only 3133 co-clusters as an adjacent pair**; 1369 and 2295 separate across sub-branches. Same-patient sublines are therefore not interchangeable at the transcriptome level, which matters when choosing between them.

### 4.4 Per-feature variance components and their sensitivity to model specification

Fitted per gene / per protein across 31 models by `lme4` REML decomposition. (`variancePartition` 1.38.1 is installed but fails at run time in this R because it calls a function that has moved to another package; the two approaches are not demonstrated to agree, and the substitution is recorded in the deposited files.)

| Term | RNA median / mean % | Protein median / mean % |
|---|---|---|
| Histotype | **5.95** / 14.94 (IQR 0–24.79) | **8.53** / 16.65 |
| Contributing centre | **3.53** / 14.50 (IQR 0–25.19) | **~0** (6 × 10⁻¹³) / 8.15 |
| Isobaric set | — | 0.85 / 6.42 |
| Patient\* | 27.27 / 31.57 | 19.80 / 27.44 |
| Residual | 32.84 / 38.99 | 37.53 / 41.34 |

**Report the median, the mean and the IQR, because the distributions are strongly skewed and the summaries disagree.** In RNA the medians order histotype above centre by 1.7-fold, but the *means* are nearly indistinguishable (14.94 versus 14.50) with almost entirely overlapping interquartile ranges. The histotype-over-centre claim in this layer therefore rests on the medians and must be reported with both. In protein the margin is far larger: histotype 8.53% against a centre median of essentially zero.

**The `patient` term is a design artefact and must not be reported as a finding.** It is the largest structured component in both layers, but it was fitted with **28 levels on 31 observations** and only **three patients are genuinely replicated**. Restricting the term to those three families collapses the median from 27.27% to **0.758%**; dropping it entirely preserves the ordering (RNA 7.09% versus 3.95%; protein 9.69% versus ~0%). The main conclusion therefore does not depend on the artefactual term, which is the point of testing it.

Adding passage as a scaled fixed covariate moves the RNA centre median from 3.53% to **exactly 0.00%** while passage itself takes 4.62%. That is a confounding disclosure, not a result about passage: **passage and centre are largely redundant in this design.**

### 4.5 The within-histotype cross-centre comparison is underpowered at 2 versus 5 models

Clear cell is the only histotype with models from two centres in usable numbers: 2 from CHUM, 5 from BC Cancer. Refitting the PCA *within* clear cell gives centre R² of **15.7 / 12.6 / 27.3%** on PC1–PC3.

Those look substantial and they are not informative. At a 2-versus-5 split the null expectation for the corresponding degrees of freedom is **E[R²] ≈ 16.67%**, two of the three observed values fall *below* it, and no test approaches significance (ANOVA p = 0.380 / 0.435 / 0.228; permutation p = 0.490 / 0.541 / 0.246). **There is no evidence of a within-histotype centre effect at this sample size, and equally no evidence against one.** This belongs in the report as a design limitation, not as evidence either way; §4.2 is the analysis that supports the histotype conclusion.

For the record, an earlier version of this analysis regressed *global* PCs on centre within clear cell and obtained 4.4–5.8%. That read as reassuring and used the wrong axes. Both versions are retained in `output/rna_within_cc_site.csv` so the difference is visible.

### 4.6 Passage number, and cross-assay passage mismatch

{{figure: figs2 | **Passage sensitivity and cross-assay passage discordance.** (**A**) RNA principal-component variance attributed to passage under three nested models, n = 31: passage alone, passage given centre, and passage given histotype + centre. The printed p belongs to `passage alone` only; the conditional quantities are decompositions with no p value in the deposited files. (**B**) Exome against RNA passage for the 13 models with both records; diagonal is y = x.}}

- Passage alone explains 7.8% of PC1 (p = 0.128) and reaches nominal significance on **PC2 only** (R² 0.159, p = 0.026).
- Conditioned on centre its partial R² on PC1 rises to 14.8%, but **that intermediate figure is not the headline**. With both histotype and centre in the model, the unique passage contribution to PC1 is **0.41%**.
- **The layers are not passage-matched, and this is a real reuse limitation independent of the variance decomposition.** Among the 13 models with both records the median absolute difference is 4 passages, but the extremes span **−17** (`OV3331`, RNA p71 versus exome P54) to **+20** (`TOV112D`, p63 versus P83). The rank correlation between the two records is **negative** (Spearman −0.251; Pearson r = −0.33, 95% CI −0.75 to 0.27, p = 0.26). Eleven of 13 are within 7 passages, so the problem is concentrated in two models. One of them, `TOV112D`, is the sole remaining endometrioid model after reassignment, so the worst-matched model is also one of the most load-bearing.

### 4.7 Collapsing sublines to one model per patient strengthens structure and reduces DE power

Both directions must be reported; showing only one would be selective.

| Quantity | 31 models | 28 patient representatives |
|---|---|---|
| PC1 R², histotype | 0.735 | **0.792** |
| PC1 R², centre | 0.313 | **0.272** |
| PC1 unique histotype | 0.424 | **0.523** |
| Overall silhouette | 0.220 | 0.208 |
| DE genes, HGS one-vs-rest | 7,913 | 6,372 |
| DE genes, CC one-vs-rest (n unchanged at 7) | 1,256 | **821** |

**Structure strengthens under collapse.** Duplicate genomes inflating an apparent separation would produce the opposite, so this is positive evidence that the histotype separation is real. Holding the 31-model gene set fixed rather than reselecting highly variable genes gives 0.781 and 0.510, so the improvement is not gene reselection.

**Differential expression weakens, and it weakens in contrasts that lose no members.** Only HGS shrinks (15 → 12), yet four of the five groups that lose nobody still lose differentially expressed genes, because the comparator "rest" group shrinks. Clear cell is the clearest case: an unchanged group of 7 loses **435 of 1,256** DE genes (−34.6%). **Power was being drawn from duplicate genomes in every contrast, not only the one containing the sublines.** The caution generalises to any cell-line panel with sublines.

### 4.8 Histotype separation differs between RNA and protein, and reverses for two groups

{{figure: figs1 | **Proteomic principal components and histotype separation.** (**A**) Protein PC1–PC2 by histotype, n = 31 models; top 2,000 most variable complete-case proteins of 6,855, centred and unscaled. There is essentially no histotype structure in these two components, which panel C quantifies. (**B**) The same points by isobaric set, on a ramp used for no other variable in this report. (**C**) Mean silhouette width per histotype for both layers: Euclidean distance on PC1–10 with groups set to the **annotated** histotype, so a negative value means that group's models sit closer to other histotypes than to each other. Unit: model.}}

| Histotype (n) | RNA silhouette | Protein silhouette | Difference |
|---|---|---|---|
| SCCOHT (2) | **0.819** | **0.028** | +0.790 |
| MMMT (2) | 0.738 | 0.410 | +0.328 |
| HGS (15) | 0.162 | 0.178 | −0.016 |
| MC (3) | 0.154 | 0.141 | +0.014 |
| CC (7) | 0.121 | **−0.004** | +0.124 |
| EC (2) | **−0.014** | 0.128 | −0.142 |
| **All (31)** | **0.220** | **0.135** | +0.085 |

The two layers should be reported separately rather than averaged. **SCCOHT separates almost perfectly in RNA and not at all in protein**, a 29-fold difference. The transcriptional identity of those two models is unmistakable; their proteomes are not distinguishable from the rest of the panel by this metric. Clear cell has a *negative* protein silhouette on 7 models, meaning the average clear-cell model sits marginally closer to another histotype's centroid than to its own in protein space. Endometrioid reverses sign between layers.

Two caveats attach to the table. Silhouettes here use the **annotated** histotype as the grouping, so they measure how well the annotation matches the geometry, not the quality of a clustering. And every n = 2 value is dominated by a single pairwise distance. The only groups with n ≥ 7 are HGS and CC, and for those two the RNA/protein difference is small (−0.016) and moderate (+0.124). The dramatic contrasts rest on n = 2 and should always be quoted with the n visible.

The practical reading is that **the two expression layers are not interchangeable for identity work**. A rare-histotype question is better asked of the transcriptome; a post-transcriptional-loss question can only be asked of the proteome (§9.3).

---

## 5 · Recovery of canonical histotype biology

Recovering canonical, histotype-specific biology from independently generated layers is evidence that the data are sound. Every result in this section is intended as validation, not discovery.

### 5.1 Differential expression, graded formal or descriptive by group size

Only the HGS contrast is `formal`; the other five are `descriptive`. The rule is in code (`03_rna_de_signatures.R:59-63`: `formal_ok <- names(n_by)[n_by >= 10]`), which selects HGS alone at n = 15, and the grading is carried in the deposited files rather than left to the reader.

| Contrast | n | Class | padj < 0.05 | Up | Down | Up **and** log2FC > 1 | Signature symbols |
|---|---|---|---|---|---|---|---|
| **HGS** | 15 | **formal** | **7,913** | 3,573 | 4,340 | 2,088 | 1,858 |
| CC | 7 | descriptive | 1,256 | 421 | 835 | 400 | 365 |
| SCCOHT | 2 | descriptive | 2,232 | 292 | 1,940 | 290 | 278 |
| EC | 2 | descriptive | 1,880 | 277 | 1,603 | 277 | 270 |
| MMMT | 2 | descriptive | 1,387 | 200 | 1,187 | 199 | 191 |
| MC | 3 | descriptive | 565 | 168 | 397 | 168 | 153 |

Down-regulated genes outnumber up-regulated in **every** contrast, overwhelmingly so in the n = 2 groups (SCCOHT 1,940 down versus 292 up). That is the expected signature of a small group against a large heterogeneous remainder, where most of the signal is genes the "rest" expresses and the small group does not, so **the up-regulated counts are the interpretable ones**. DE files hold one row per gene identifier while signature files hold one row per distinct gene symbol, so 2,088 and 1,858 are not the same quantity.

### 5.2 Two of five expected pathway programmes recover robustly

fgsea 1.34.2 on GO biological process with `nPermSimple = 50000`, because the default 1,000 permutations did not resolve p-values near the 0.05 boundary at this n and two of five programme grades changed status across seeds. At 50,000 permutations all five grades are stable across three seeds.

| Programme | Histotype (n) | Best GO-BP term | NES | BH p | Grade |
|---|---|---|---|---|---|
| Cell cycle / mitosis | SCCOHT (2) | replication fork processing | 2.33 | **8.9 × 10⁻⁵** | **Recovered robustly** |
| Oxidative / glutathione / detox | CC (7) | sulfur amino acid metabolic process | 2.14 | **0.0072** | **Recovered robustly** |
| DNA repair (HR / DSB) | HGS (15) | regulation of DSB repair via HR | 1.64 | 0.075 | Suggestive (nominal 0.003) |
| Glycan / oligosaccharide | MC (3) | oligosaccharide catabolic process | 1.77 | 0.12 | Suggestive (nominal 0.010) |
| EMT / migration / ECM | MMMT (2) | negative regulation of VSMC migration | 1.57 | 0.34 | Suggestive at model level; **not recovered** under patient collapse (0.46) |

**Recovery does not track group size.** The two programmes that recover robustly are in groups of n = 2 and n = 7, while the n = 15 group's programme does not reach adjusted significance. What it tracks is how *distinctive* the programme is: the SCCOHT cell-cycle signal has 26 positively significant terms because those two models are transcriptionally extreme (§4.8). **The grade is stable while the single best term is not.** At a different seed the SCCOHT programme's top term becomes "nuclear chromosome segregation", so the programme-level grade should be quoted rather than a named term.

**One selection effect is not corrected and should be disclosed rather than absorbed.** The programme definitions are regular expressions, and the term reported for each programme is the **best-scoring GO-BP term matching that regex** within a much larger result set. The Benjamini–Hochberg adjustment is applied across the GO result set, not across the choice of which matching term to report, so the two robust grades are confirmatory only in the weaker sense that a prespecified *programme* — not a prespecified *term* — recovered. **Clear cell and SCCOHT are the robust pathway checks; high-grade serous and mucinous are exploratory and MMMT is not recovered**, and the five should not be presented as an equivalently recovered set.

### 5.3 Marker recovery: 15 of 25 land on the inferential unit, against a correlation-preserving null

{{figure: figs7 | **Histotype marker effect sizes with bootstrap intervals.** Unit, different from every other figure here: symbol-summed log2(TPM + 1) on **28 patient representatives**, not 31 models; 25 markers graded, intended-group sizes printed. (**A**) Cohen's *d*, **signed**: point = *d*, bar = 2,000-resample bootstrap 95% interval, axis clipped at ±6 with three upper bounds printed as arrows. Glyphs: ↓ = loss marker; \* = survives Benjamini–Hochberg correction at 0.05 (only 4 of 25 do); † = `SMARCA2`, whose Wilcoxon p equals the smallest attainable for its group size. (**B**) AUC oriented to the expected direction, axis fixed to [0, 1]. Because *d* is signed in A and oriented in B, a loss marker succeeds at negative *d* in A but above 0.5 in B.}}

- **16 of 25 markers land in their intended histotype on all 31 models, and 15 of 25 on the 28 patient representatives**, under the rank rule of Fig. 3E and 3F: intended group among the top 2 of 6 histotype means for an up marker, bottom 2 for a loss marker, with a 1-log2-TPM expression floor.
- **Only 4 of 25 survive Benjamini–Hochberg correction**, all at adjusted p = 0.0437: `KRT20` (MC), `SPP1` (CC), `SMARCA2` (SCCOHT), `HNF1B` (CC). Eight reach nominal p < 0.05.
- **15 of 25 effect-size intervals cross zero.**

**The null model matters here, and the binomial version is not the right one.** An exact binomial test at p = 1/3 treats the 25 markers as independent Bernoulli trials. They are not: markers intended for the same histotype are co-regulated — `PAX8`, `WT1` and `MUC16` all track one serous programme — so a single relabelling of the models moves several markers together, and the binomial p is anti-conservative by an unknown amount. The appropriate null permutes the **histotype labels** and recomputes the entire landing count, which preserves the marker-marker correlation structure, the group sizes and the absolute-expression floor. 20,000 joint permutations at seed 1234 (`scripts/04_rna_markers_genesets.R`, deposited to `output/rna_marker_recovery_permutation.csv`):

| Analysis unit | Observed | Null mean | Null 95th pct | **Permutation p** | Binomial p |
|---|---:|---:|---:|---:|---:|
| 31 line models | 16 / 25 | 6.70 | 11 | **0.0013** | 0.0016 |
| **28 patient representatives** | **15 / 25** | 6.82 | 12 | **0.0046** | 0.0056 |

**The conclusion strengthens rather than weakens.** Recovery is far outside what jointly permuted labels produce on either unit — the observed counts exceed the null 95th percentile in both cases — and the permutation p is the quotable one. Note that the null mean is 6.7–6.8 of 25, i.e. about 27% rather than the 33% the binomial assumes, because the expression floor in the up rule makes landing strictly harder than a rank test alone. **The 28 patient representatives are the inferential unit** (§4.7 and §1.2), so **15 of 25 at permutation p = 0.0046** is the figure to quote; the 31-model count is the resource-level description.

The governing caveat is the Wilcoxon p-value floor at small n: the minimum attainable two-sided p is **0.00529** at 2 versus 26 and **0.000611** at 3 versus 25. `SMARCA2` in SCCOHT achieves the floor exactly, with *d* = −2.08 and a **perfect oriented AUC of 1.000**; it cannot do better. **Effect size and AUC, not the p-value, are the quantities to read for small groups.** It survives adjustment only because three other markers have comparably small p-values.

Equally, **small-group point estimates must not be ranked against one another.** `MUC5AC` on n = 3 gives *d* = 1.92 with a 95% interval of **[−0.44, 20.64]**, a large point estimate whose interval spans no effect to an effect of 20. The widest intervals are all in the mucinous panel (MUC5AC 21.1 units wide, KRT20 17.2, CDX2 8.6).

**One marker set is a clean negative, and it has two readings.** No endometrioid marker lands: `ESR1` ranks 6th of 6 histotype means (highest in HGS) and `PGR` 5th of 6, with mean values of 0.026 and 0.002 log2 TPM in the intended group and `VOA4395` reading **exactly 0.000** for both. Either the annotation is unsupported by expression, **or** these receptors do not survive two-dimensional culture. Panel-wide context keeps the second reading live: `PGR` has a median of exactly 0.000 across all 31 models with **22 of 31 at zero**, and `ESR1` a median of 0.286 with 3 at zero. The defensible conclusion is that **endometrioid identity is not marker-validated in this resource**, not that the annotation is wrong. The same culture-drift alternative applies wherever a lineage marker reads low, including the mucinous assessment in §9.5.

---

## 6 · Transcript–protein concordance

### 6.1 Median transcript–protein correlation is 0.40, within the published range

On the **30 models with both layers** (`VOA6861` is RNA-only, `VOA14993` protein-only; Fig. 3C):

| Metric | n | Median | Range / IQR |
|---|---|---|---|
| Per-model Spearman | 30 models | **0.408** | 0.330–0.479 |
| Per-gene Spearman | 7,894 genes | **0.397** | q25–q75 0.189–0.578 |
| Per-gene Spearman, complete case | 6,688 genes | 0.408 | — |
| Per-gene Pearson | 7,894 genes | 0.436 | — |

- **The two medians answer different questions and should not be read as replicating each other.** The per-model value correlates thousands of genes within one model and is dominated by stable gene-to-gene abundance differences; the per-gene value asks whether model-to-model differences agree between layers, and it is the relevant reuse statistic. Their numerical closeness (0.408 versus 0.397) is reassuring but is **not** an internal-consistency check, because the estimands differ — two different quantities landing at the same value is a coincidence of scale, not agreement.
- **The per-model range is narrow (0.330–0.479), so no model is an outlier.** The concordance limit is a property of the assay pair, not of particular samples, which is what to establish before anyone reads a low correlation as a sample-quality problem.
- **The estimate is insensitive to the inclusion rule.** Raising the minimum from 10 to 30 models moves the median only from 0.397 to 0.408 while shrinking the gene set by 15%.
- **It is also insensitive to same-donor sublines.** Three of the 30 dual-layer models are sublines, and collapsing to the 27 patient representatives moves the per-model median from 0.408 to 0.423 and the per-gene median from 0.397 to 0.417 (`integ_rnaprot_patientrep_sensitivity.csv`). Both move **up**, so repeated donors were very slightly deflating the concordance rather than inflating it — the opposite of the concern.
- These values sit within the published range for the same comparison: CPTAC ovarian 0.38 / 0.45, ProCan 0.42, CCLE 0.48, Jarnuczak 0.58, with a reproducibility ceiling around 0.72.
- **A tenth of the measured proteome has an inverse point estimate, and almost none of it is statistically established.** 813 genes (10.3%) have a negative Spearman correlation, the most extreme reaching −0.714 — but at 10 to 30 paired models most are indistinguishable from zero. Under the conventional asymptotic Spearman test (t approximation, df = n − 2), **50 are nominally significant and 29 (0.37% of the set) survive Benjamini–Hochberg correction**, against 3,599 genes (45.6%) with an FDR-supported positive association. **Two wording corrections travel with these counts.** An earlier version of this report said that for the 813 genes "a transcript measurement is not a proxy for protein in any direction" — that asserts a relationship the estimates do not support, and it also conflates two different things. What the 29 establish is an **FDR-supported inverse association**, which is *not* an absence of proxy value: a protein that reliably falls as its transcript rises is predictable, simply in the opposite direction. Proxy value additionally depends on predictive error and external validation, neither of which is assessed here. **So: treat the 813 as candidates for discordant regulation, quote the 29 as inverse associations, and do not describe either set as genes for which transcript fails as a proxy.** Per-gene p and q values are deposited in `integ_rnaprot_cor.csv` and the tier counts in `integ_rnaprot_negative_genes.csv`.

### 6.2 The protein matrix has a 3.3-fold narrower cross-model spread than RNA, consistent with isobaric ratio compression

Values below are on the **complete-case** subset (6,689 genes measured in all 30 dual-layer models), so the three statistics share one denominator; the ratio column is the median of the per-gene ratio, not the ratio of the medians.

| Statistic | RNA median | Protein median | Median per-gene ratio | Implied fold narrowing |
|---|---|---|---|---|
| **IQR** (primary, floor-insensitive) | 1.113 | 0.336 | **0.299** | **3.34×** |
| **SD** (primary, floor-insensitive) | 0.848 | 0.274 | 0.317 | 3.15× |
| Range (secondary, floor-sensitive) | 3.386 | 1.160 | 0.337 | 2.97× |

**97.8% of the 7,896 paired genes have a smaller protein than transcript interquartile range**, and 99.2% by SD (Fig. 2C). The narrowing varies by target from 3.84× for `CD276` to **9.64×** for `MSLN` (Fig. 6C).

**What this measures, and what it does not.** The ratio is a robust property of these two processed matrices. It is **not** an isolated estimate of a TMT compression factor, and it should not be quoted as one. The two layers sit on different scales with different normalisation — RNA as log2(TPM + 1) with a hard zero floor, protein as log2 reporter intensity relative to a plex pooled internal standard — and at least four mechanisms other than reporter-ion co-isolation narrow the protein axis relative to the transcript axis: post-transcriptional buffering of protein abundance against mRNA abundance, protein turnover and homeostasis, the pooled-internal-standard ratioing itself, and structured missingness that preferentially retains widely quantified (and typically abundant, low-variance) proteins. **The defensible statement is that the protein matrix has a narrower cross-model spread than the RNA matrix, consistent with known TMT ratio compression and with post-transcriptional regulation.** No design here can partition the 3.34× between them.

The obvious objection to a compression claim is that the RNA scale is inflated by zeros, and the artefact runs the other way. **Removing the floor increases the estimated compression rather than deflating it**: for all eight dual-layer ADC targets and for 65.6% of paired genes, the IQR ratio is below the range ratio. Stratifying by whether a gene's RNA scale reaches zero confirms it. The 531 floor-touched genes have nearly double the RNA range (6.28 versus 3.38 log2) yet a slightly *smaller* IQR ratio (0.275 versus 0.299), because floor-touched genes are also the more variable ones.

Finally, **genes with more protein spread show better cross-assay agreement**, not less (Spearman between protein IQR and per-gene ρ = +0.332). That is the pattern expected if a narrow protein axis limits concordance for low-spread genes, rather than concordance limiting the apparent narrowing. And the external benchmarks, which share the isobaric limitation, land at the same value.

**On "ceiling".** It is reasonable to say that **0.40 is where this assay pair lands in this configuration, and that the published cell-line and tumour benchmarks agree**, so a median of 0.40 is not a verdict on these data. It is not defensible to say that a 3.34-fold compression *sets* a concordance ceiling: that would require isolating compression from the other mechanisms above and demonstrating the causal link, and neither is possible here. The complementary constraint on how finely the protein layer can be read — technical error accounting for 27–44% of observed cross-model protein variance — is in §3.3, and the two together bound the resolution without either being a causal attribution.

---

## 7 · Genomic fidelity of the exome layer

{{figure: fig4 | **Genomic fidelity and model identity.** (**A**) Variant-filtering waterfall for one exemplar model and all 22 exome models, y axis log10; stages are all Mutect2 records with FILTER ignored, `FILTER == PASS`, then population-AF ≤ 0.001 plus coding non-synonymous. The `TP53 retained in 17/17` control is per **model**; the patient-level equivalent is 11 of 11. (**B**) Identity specificity against DepMap: 5 models with a namesake, each against 67 DepMap ovarian models, x = Spearman ρ on 2,000 highly variable genes. Grey tick = one non-self model; open circle = best non-self match, named; filled circle = self-match; segment = margin. Printed per row: Δρ = self ρ − best non-self ρ, and z = (self ρ − median of the 66 non-self ρ) ÷ their SD. (**C**) Oncoprint of retained candidate driver calls; unit 22 exome models / 16 patients, 14 genes, 50 calls. Cell fill = somatic-confidence tier; centre bar = variant class; a cell's tier is the most confident among its calls. Top track = coding candidate counts, truncated at 460 with a break glyph on the one model at 1,416. Right-hand bars count patients with at least one call out of all 16, so the *TP53* row reads 12 rather than 11 — that is the panel-wide denominator, not the HGSC one. **A dagger means no Tier-1 call exists in any model for that gene.** (**D**) Per-arm copy-number frequency, unit **% of 11 HGSC patients**, so one patient is 9.1%. Call rule |log2 copy ratio| > 0.20 over > 50% of arm length after per-sample autosome probe-weighted median centring; acrocentric p-arms omitted.}}

### 7.1 Filtering reduces 557,392 raw calls to 6,036 coding candidates and removes the tumour-only artefacts

| Stage | Variants | Filter |
|---|---|---|
| Raw Mutect2 records | **557,392** | — |
| `FILTER == PASS` | **15,692** | Mutect2 PASS only (2.82% of raw) |
| Rare | 15,609 | population AF ≤ 0.001 |
| **Coding non-synonymous** | **6,036** | maftools non-synonymous set (1.08% of raw) |

The diagnostic value of this cascade (Fig. 4A) is what it does to specific genes. **Models carrying a coding non-synonymous *ATM* call fall from 18 of 22 (82%) to 2 (9%), and *ATR* from 17 (77%) to 1 (5%), while *TP53* is retained in all 17 high-grade serous models.** The recurrent artefacts of tumour-only calling are removed and the events the histotype requires are preserved, which is the correct shape for a filter. It is also why the pre-fix driver landscape, with *ATM*, *ATR* and *BRCA2* apparently mutated in the majority of models, was uninterpretable.

**Two things the cascade does not show.** It is not a measure of true somatic content: tumour-only calling cannot be fully cleaned, and stage-3 calls remain candidates. And the median coding burden of **206.5 per model** is an **upper bound** on residual germline, not an estimate of it. No residual-germline estimate exists in any deposited file, and none should be quoted.

### 7.2 *TP53* is mutated in 11 of 11 high-grade serous patients; driver frequencies require tier composition

Every retained call in a canonical driver gene was assigned a somatic-confidence tier: **Tier 1** highest, **Tier 2** intermediate, **Tier 3** recorded but germline not excludable. Fifty calls across 14 genes and 22 models: **27 Tier 1 / 11 Tier 2 / 12 Tier 3**.

**The tier is a prioritisation heuristic, not a validated somatic classifier**, and the name should not be read as one. The rules are literature- and variant-class-based triage — hotspot or truncating variants in selected canonical genes receive higher confidence, rare *BRCA1/2* calls are deliberately held at Tier 3 — and they are project-specific and externally uncalibrated. They do not read each model's histotype label, which avoids the more serious circularity, but without a matched normal no rule can establish somatic origin, so a Tier 1 call is a *high-confidence candidate*, not a somatic mutation.

- ***TP53* is mutated in 17 of 17 high-grade serous models and 11 of 11 high-grade serous patients, 100%, every call Tier 1 or 2** (Fig. 4C). This is the strongest single genomic validation in the resource. Panel-wide the figure is 12 of 16 patients (75%), or 18 of 22 models (81.8%); the two units must not be swapped.
- **Frequencies are reported per patient**, because subline duplication inflates them: *CDK12* is called in 6 of 22 models (27.3%) but only 3 of 16 patients (18.8%), and only 2 patients at Tier 1–2 (12.5%).
- **Tier composition must travel with any frequency.** *SMARCA2* (3 calls), *BRCA2* (2) and *NF1* (1) rest **entirely on Tier 3**, and *CDK12*'s nine calls are five Tier 2 and four Tier 3 with **no Tier 1 call in any model**. A frequency computed entirely from non-Tier-1 calls in tumour-only data is a *candidate* frequency and is not a somatic frequency.

The retired `germline_like_vaf` flag was worse than nothing. In cultures of essentially 100% tumour purity, a bona-fide somatic hotspot under loss of heterozygosity has a variant allele fraction approaching 1, so the flag fired on *TP53* R175H and Q192\*, the calls it should have protected. The tiering scheme replaces it, and applying that scheme consistently is what removes both *BRCA2* calls from the headline. **The correct statement is that no defensible *somatic* *BRCA1/2* call can be made from this dataset** — not that no *BRCA1/2* alteration exists in these models. Tumour-only calling with population-AF filtering cannot classify a rare germline *BRCA1/2* variant as somatic in the first place, so absence of a somatic call here carries no information about the models' *BRCA1/2* status, and a reuser needing that status must obtain it another way.

### 7.3 Fraction of genome altered spans 30-fold across the panel, from 0.635 in HGSC to 0.021 in low-grade serous

{{figure: figs5 | **Genome-wide copy-number profiles.** 23 models × 300 columns, rows split by histotype. Each column is one consecutive 10 Mb autosome bin, p→q within each chromosome, valued as the segment-length-weighted mean CNVkit log2 copy ratio after per-sample autosome probe-weighted median centring (median 185 autosomal segments per model, range 88–481). Fill is symmetric about 0 and **clamped at ±1.5**, so more extreme values render identically; grey = no exome coverage. The right-hand bar gives fraction of the assessed autosome with |log2| > 0.20, which makes the near-blank row quantitative: `TOV81D`, the single low-grade serous model, at 0.021. `‡` marks the model with copy-number data but no mutation file. The five public healthy exomes in the panel of normals are not confirmed to share the capture design, so the relative pattern is more trustworthy than the absolute log2 level.}}

Median fraction of genome altered (autosomes, |log2| > 0.2) is **0.6223** across 23 models (Fig. S5), and the ordering by histotype is exactly what the histotypes predict:

| Histotype | n | Median FGA | Range |
|---|---|---|---|
| HGS | 18 | **0.6345** | 0.269–0.844 |
| CC | 2 | 0.371 | 0.071–0.671 |
| MC | 1 | 0.321 | — |
| EC | 1 | 0.226 | — |
| **LGS** | 1 | **0.0209** | — |

**The spread between high-grade and low-grade serous is 30-fold, in the predicted direction, on data that never saw a histotype label.** That is a positive control on the copy-number layer, and it should be reported as a **panel range rather than a histotype effect estimate**: it compares a well-sampled group of 18 HGSC models with **one** low-grade serous model, so there is no variance estimate on the low side and no test to run. It rests on n = 1 or 2 outside HGSC and should always be quoted that way; the two clear-cell models are themselves 0.071 and 0.671, which is most of the panel's range inside one histotype. Loss exceeds gain panel-wide (0.326 versus 0.261 of the autosome).

Two further constraints bound how precisely these values can be read, both stated in §2.3 and repeated here because they attach to every FGA number. The five pooled normals come from an unrelated study and are **not confirmed to share the tumours' capture design**, so target-set mismatch can create systematic log-ratio structure and segmentation artefacts. And per-sample autosomal median centring removes an unknown baseline in highly aneuploid models, so these are **relative total-copy profiles, not absolute copy states**. Both support coarse within-sample gain/loss comparison under a declared threshold and neither supports a claim of precision comparable to a matched, capture-compatible copy-number study.

**Chromosome X is excluded, and the exclusion is the artefact-free choice.** chrX behaves as a pooled-normal sex-composition artefact here: median chrX fraction altered is 0.964 and the per-model median log2 copy ratio swings from −1.23 to +0.83 with the sign flipping between models. Excluding it barely moves high-FGA genomes (all 23 models 0.6057 → 0.6223; HGSC 0.6271 → 0.6345) but corrects quiet ones sharply: **`TOV81D` moves from 0.073 to 0.021**. The X-inclusive column is retained only as a legacy field. Applying the corrected metric consistently is also what moves one model's genomic-consistency call from *consistent* to *partial* (§9.6), which is the stated rules being applied rather than a loss.

### 7.4 Arm-level frequencies recover the canonical HGSC events but move 18–46 pp with the calling threshold

Unit: **HGSC only, 18 models / 11 patients** (Fig. 4D). Anyone quoting these as panel-wide frequencies is wrong by construction; `arm_freq_at()` filters to HGS before tabulating, and every row of the deposited file carries `n_lines = 18, n_patients = 11`.

| Arm | Direction | % of 11 patients | Sweep range | Total span |
|---|---|---|---|---|
| 20q | gain | **91** | 73–100 | 27 pp |
| 3q | gain | **82** | 64–82 | **18 pp** |
| 17p | loss | **82** | 64–91 | 27 pp |
| 8q | gain | **73** | 55–82 | 27 pp |
| 13q | loss | **64** | 45–82 | 37 pp |
| 19q | gain | **55** | 27–73 | **46 pp** |

The canonical HGSC events all recover. But across a ten-combination threshold sweep (|log2| ∈ {0.10 … 0.40} at majority 0.5, plus majority ∈ {0.3 … 0.7} at |log2| 0.20) **the six headline frequencies move by 18–46 percentage points**. `3q` gain is the most stable arm and never exceeds its headline; `19q` gain is the least stable, moving from 27% to 73%. **A frequency quoted without its threshold is not reproducible**, so 19q's headline 55% should be quoted with its range or not emphasised.

The two genomic layers respond differently to patient collapse: **collapse distorts point-mutation frequencies badly and arm-level frequencies barely** (*CDK12* 27.3% → 18.8%, versus 3q gain 83% of models → 82% of patients). Arm-level events are trunk events shared across sublines, so they survive collapse; a subline-private point mutation does not.

### 7.5 Genomic HRD is not computable from total copy number

CNVkit produced **total** copy number only. No B-allele frequency is available, no VCF was passed to CNVkit, and the tumour BAM files are not archived. Consequently **loss of heterozygosity, allelic imbalance, biallelic inactivation and genomic homologous-recombination-deficiency scores cannot be derived from this record**, and none was estimated; `scarHRD` is correctly absent from the environment. A genuine score would require recovering the exome BAMs and running Sequenza → scarHRD. The expression-based "HRD" signature that appeared in earlier drafts was dropped: it measures an expression state, not a genomic scar, and labelling it HRD was the error.

---

## 8 · Mutational signatures and the TOV21G mismatch-repair candidate

{{figure: fig5 | **Rare-histotype and flagged models.** (**A**) Ranked coding-candidate count per exome model, n = 22; fill = indel fraction, dashed line = panel median 206.5. Hypermutator flag: robust z > 5 **and** > 3× median, robust z being (n − median) ÷ median absolute deviation. (**B**) The 96-context substitution spectrum of `TOV21G` on 2,417 exome-wide PASS single-nucleotide variants, GRCh38 — substrate is exome-wide PASS variants at population AF ≤ 0.001, **not** the coding candidate set, called tumour-only, with no exome-to-genome renormalisation applied. (**C**) Cosine screen against three COSMIC v3.2 signature groups; each bar is the **maximum** cosine within that group. This is a screen, not an attribution, and SBS1, SBS5, SBS6, SBS15 and SBS44 are mutually similar. (**D**) Relative exposure of the mismatch-repair signature group from a strict refit, n = 22, sorted ascending; 200 bootstrap replicates, unit = fraction of that model's total fitted mutations. The open circle with a segment is **SBS6's own** exposure and bootstrap interval, not the group sum. (**E**) Müllerian-versus-gastrointestinal markers across the three mucinous models in both layers; z is computed per gene across all 31 RNA and all 31 protein models, **not** across the three columns shown, since within-column z would guarantee a one-high/two-low pattern by construction. Fill clamped at z = ±3; pale grey = not measured in that layer. (**F**) SWI/SNF subunit standing, 7 models × 3 genes × 2 layers, with the exome call at right. **The number in each cell is that model's rank among the 31 panel models, 1 = lowest**; fill = z across those same 31 models. Rows are the 6 models called deficient plus `OV2295`, whose truncating call is Tier 3 only and is not admitted — which is why the count is 6 rather than 7.}}

### 8.1 TOV21G carries 6.9× the panel median coding burden

- Median coding burden is **206.5** per model (range 133–1,416; Fig. 5A). **`TOV21G` carries 1,416 coding candidates, 6.86× the panel median, robust z = +21.2.** The next-highest model has 413 (robust z = 3.6) and is not flagged.
- `TOV21G` is extreme on both composition metrics as well: indel fraction **0.279** against a panel median of 0.098 (next-highest 0.187) and Ts/Tv **2.72** against 0.77 (next-highest 1.21).
- All three quantities point the same way, and one model out of 22 is flagged. The classification therefore rests on a single very large outlier. A 6.9-fold excess with an indel-rich, transition-shifted spectrum is nonetheless not a marginal call.

### 8.2 Cosine screening is inconclusive; the bootstrap refit supports a mismatch-repair exposure in TOV21G alone

- **Cosine screening** gives only `TOV21G` a positive mismatch-repair-versus-clock margin (+0.225 against a panel median of −0.162). But its top three mismatch-repair matches sit close together (SBS6 0.877, SBS44 0.835, SBS15 0.809), and SBS1/5/6/15/44 are mutually similar, so **a moderate cosine margin does not discriminate strongly on its own.** Across the panel, the closest reference signature is SBS5 in 17 of 22 models.
- **The bootstrap refit separates the model far more cleanly.** `TOV21G` reaches a mismatch-repair-class relative exposure of **0.733** against 0.000–0.292 for every other model; its top signature SBS6 has a relative exposure of 0.428, was selected in **199 of 200 bootstrap replicates**, and has a bootstrap 95% interval that **excludes zero**. Its reconstruction cosine of **0.9765 is rank 1 of 22** (every other model falls in 0.892–0.960). **The four next-highest mismatch-repair exposures all have bootstrap lower bounds of 0.0**, so `TOV21G` is the only model whose exposure is bootstrap-supported.
- **Two caveats travel with every signature number.** No exome-to-genome trinucleotide renormalisation was applied, because the capture kit is unknown, so relative exposures are approximate and only within-model ranking and contrasts between identically processed models carry. And the calls are tumour-only, so germline cannot be excluded from the spectra.

In the other direction, SBS3, the homologous-recombination-associated signature, is the top refit signature in three HGSC models. Two of them are sublines of one patient, so that is **2 patients, not 3**, at bootstrap selection frequencies of 0.41–0.76, which is weak to moderate. It is not evidence of HRD, and §7.5 explains why no HRD score exists here.

### 8.3 No causal coding lesion is present, consistent with epigenetic *MLH1* silencing

- **No mismatch-repair enzyme (*MLH1*, *MSH2*, *MSH6*, *PMS2*) and no *POLE* coding mutation is present in any model.**
- Two mismatch-repair-panel coding hits exist and neither closes the case. `TOV21G` carries **`EPCAM` p.2_3fs**, which *is* a truncating variant in a panel gene. But *EPCAM* sits on mismatch-repair panels because deletions at its **3′** end silence the neighbouring *MSH2* in *cis*, and this is a frameshift insertion at codon 2–3 of a 314-residue protein, the extreme **5′** terminus. It truncates *EPCAM* itself and cannot act through 3′ read-through. Four further caveats apply: it is novel in dbSNP; it is an insertion, so there is no variant-allele-fraction support; the calls are tumour-only, so germline cannot be excluded, and germline *EPCAM* deletions exist; and it sits in the one hypermutated model, where passenger burden is highest. The second hit, `MSH3` p.Y334F in `TOV2881EP`, is a missense of uncertain significance and is not an established loss mechanism on its own.
- The absence of a convincing coding lesion is **consistent** with the mechanism this phenotype usually has: **epigenetic *MLH1* promoter methylation**, which is typical of MSI-high clear-cell carcinoma and which exome sequencing cannot detect.

**Interpretation.** `TOV21G` is best described as a **candidate mismatch-repair-deficient / MSI-high clear-cell model**, on converging evidence from burden, composition, spectrum, bootstrap-supported signature exposure, and independent corroboration of its burden-outlier status in DepMap (568 damaging variants against 7–18 for the other two models with both sources). The other clear-cell model with exome data, `TOV3392D`, is not hypermutated, so the signal is model-specific rather than a histotype property or a processing artefact. This is a **feature, not a defect**, since MSI-high clear-cell models are rare and useful for immunotherapy and mismatch-repair biology. It needs MMR immunohistochemistry, MSI-PCR and germline testing before any mechanistic claim.

---

## 9 · Identity and authentication

### 9.1 All five models with a DepMap counterpart self-match at rank 1 of 67

Five models have a namesake in DepMap Public 24Q4 (`OV90`, `TOV21G`, `TOV112D`, `COV434`, `BIN67`) and were compared against 67 DepMap ovarian expression models (Fig. 4B) on 2,000 highly variable genes drawn from 19,062 shared Entrez identifiers.

| Model | Self ρ | Rank | Best non-self | Margin Δρ | z vs non-self | Reciprocal best |
|---|---|---|---|---|---|---|
| `OV90` | 0.879 | **1 / 67** | OVCAR5 | **0.434** | 11.5 | yes |
| `TOV21G` | 0.805 | **1 / 67** | OVTOKO | 0.231 | 5.5 | yes |
| `TOV112D` | 0.839 | **1 / 67** | A2780 | 0.204 | 4.2 | yes |
| `COV434` | 0.821 | **1 / 67** | BIN67 | 0.139 | 4.8 | yes |
| `BIN67` | 0.739 | **1 / 67** | COV434 | **0.036** | 4.0 | yes |

**All five self-match at rank 1 of 67 and all five are reciprocal-best matches.** On all shared genes rather than highly variable genes the self-match correlations are higher still (0.916–0.947) and rank 1 in every case.

Report all five margins rather than a range, and explain the smallest rather than hiding it. **`BIN67`'s margin of 0.036 reflects histotype, not doubtful identity**: its best non-self match is `COV434`, the other model of the same rare entity, reciprocally so, and both carry the same SWI/SNF lesion. Two models of the same rare disease with the same pathway lesion are expected to be transcriptionally close, and `BIN67`'s self-match is still **4.0 standard deviations** above the non-self distribution.

**Specificity, not magnitude, is the interpretable quantity here**, which is why the two disclosed asymmetries in the comparison do not undermine it: our side sums linear TPM per Entrez identifier then takes log2 whereas DepMap averages log2(TPM+1), and the highly variable genes were selected on the combined 31 + 67 matrix. Both inflate absolute ρ for self **and** non-self pairs alike, and neither affects the rank or the margin.

A driver cross-check on the three models with both in-house exome data and a DepMap record found **10 of 14 gene–model pairs agreeing (71%)** where either source makes a call. The four disagreements form a coherent pattern rather than a contradiction: two of the three calls made only here are in-frame deletions and the third is Tier 3, which are the call classes least likely to be reproduced by a different pipeline. Burdens are **not** comparable between sources, our counts being 2.5–19× higher because the filters differ, but the **ordering** is preserved. That `TOV21G` is the burden outlier on both sides is itself the corroboration.

### 9.2 No in-house STR profiling was performed; 30 of 42 models have a public record

**No short-tandem-repeat profiling and no mycoplasma testing were performed in this project.** This is the resource's principal identity limitation, and the public record is partial:

- **30 of 42 models** are found in Cellosaurus by exact name match, and all 30 carry a documented STR profile; **12 are not found at all** (no fuzzy matches accepted).
- Those profiles were **deposited by other laboratories**: of the 26 with a named source, **25 are a single personal communication** (Mes-Masson) and one is Garson.
- The 12 models with no public record include the identity-doubt models below, which is why in-house STR and immunohistochemistry are requested specifically for them.
- `output/auth_perline_table.csv$STR_status` is now populated for all 42 rows (30 carry `reference profile in Cellosaurus <ACC> (n markers); no in-house profile`, 12 carry `no Cellosaurus record and no in-house profile`), derived by joining `cellosaurus_str_status.csv` with a `stopifnot` that no row is NA, and it picks up `output/str_inhouse.csv` automatically when in-house profiles arrive.
- One model carries a `problematic_flag`, verbatim: *"Misclassified. Originally thought to be an ovarian granulosa cell tumor but seems to be a small cell carcinoma of the ovary, hypercalcemic type (SCCOHT) (PubMed=33328126)."* That is `COV434`, and this resource's own evidence is consistent with it.

### 9.3 Six models are flagged SWI/SNF-deficient, and BIN67's loss is detectable only in protein

Per-gene loss calls across the panel: *SMARCA4* 4 models, *SMARCA2* 4, *ARID1A* 1, *SMARCB1* 0. **Six models are called SWI/SNF-deficient**: `BIN67`, `COV434`, `TOV112D`, `VOA4841`, `OV2085` and `TOV21G`, with per-model evidence strings deposited. The calling rules are layer-appropriate by design: *SMARCA4* loss on a Tier 1–2 truncating variant **or** protein z ≤ −1 (because the loss is often post-transcriptional, so RNA evidence is not required); *SMARCA2* loss on a truncating variant **or** RNA rank ≤ 4 of 31 **or** RNA z ≤ −1.5 (for epigenetic silencing).

**`BIN67` is the instructive case** (Fig. 5F). Its *SMARCA4* **mRNA sits mid-panel** (rank 23 of 31, z +0.62) while its *SMARCA4* **protein is second-lowest** (rank 2 of 31, z −1.96). The protein layer detects a loss the transcript layer misses entirely, so **an RNA-only authentication of this panel would return a false negative on a model whose defining lesion is unambiguous.** That is the strongest justification in the resource for depositing both expression layers, and it generalises to any gene whose loss is post-transcriptional.

Two features limit how strongly the `BIN67` result can be read. `BIN67` sits only **0.009 log2** above the panel minimum, which belongs to `VOA4841`, a clear-cell model itself called *SMARCA4*-deficient, and those two lowest models are separated from the third-lowest by only 0.105 log2. So the SCCOHT protein calls are **consistent with, rather than proof of, complete loss**, and the 0.009-log2 margin is equally the limit of the approach.

`OV2085` is the mirror case and a caution against reading the label as the complex: its *SMARCA4* protein is the panel **maximum** (z +2.48, rank 31 of 31) while its *SMARCA2* protein is the panel **minimum** (z −2.01, rank 1). "SWI/SNF-deficient" on this model means *SMARCA2* specifically.

**One model leaves the deficient set when the tiering rule is applied.** `OV2295` carries a truncating *SMARCA2* call that is **Tier 3 only**, so it is not admitted and the count moves from 7 to 6. A resource whose stated confidence rules are applied consistently, including against itself, is more trustworthy than one whose counts are maximal.

**Two corrections to how this was previously stated, both of which narrow the claim.**

First, the count. Karnezis et al. 2021 (PMID 33328126) re-assigned **two** lines: `COV434` to SCCOHT and `TOV112D` to dedifferentiated ovarian carcinoma. **`BIN67` is not among them** — it is a long-established SCCOHT model and was not reclassified in that study, so counting it as a third recovered reclassification, as earlier versions of this report did, is wrong. The resource's evidence on `BIN67` corroborates a known identity; it does not recover a reassignment.

Second, the verb. Molecular profiles cannot *perform* a pathology-based reclassification, and "independently recovers" overstates what expression, protein and exome data can do. What the three layers establish is that **`TOV112D` carries a SWI/SNF-null rather than an endometrioid molecular pattern** — a *SMARCA4* truncation with low protein, silenced *SMARCA2* mRNA and *TP53* R175H — which is independent molecular evidence **consistent with** the published reassignment, arrived at without reading the pathology. Cellosaurus still lists `TOV112D` as endometrioid, so on this point the resource is ahead of the databases while remaining behind the histopathology that actually made the call. `COV434` and `BIN67` are likewise molecularly consistent with SCCOHT and externally corroborated by DepMap *SMARCA4*-damaging calls and the Cellosaurus problematic flag.

**The "SWI/SNF-deficient" label is a prioritisation flag, and its thresholds are conventions.** The six-model set combines *SMARCA4* protein z ≤ −1, *SMARCA2* RNA rank ≤ 4 of 31 or z ≤ −1.5, and selected Tier 1–2 truncations — biologically motivated, but not externally calibrated, and the umbrella term merges loss of different subunits with different disease implications (compare `OV2085`, above, whose *SMARCA4* protein is the panel maximum). **Report the gene- and layer-specific evidence first, and describe newly flagged models as candidate [subunit]-loss models pending immunoblot or IHC**; the per-model evidence strings in `auth_swisnf_panel.csv` are the primary record and the aggregate count the summary.

### 9.4 Recomputed molecular subtype calls are unreliable on pure cultures

{{figure: figs6 | **Molecular subtype calls with posterior probabilities.** n = 15 HGSC models with RNA. Calls were **computed here** with consensusOV 1.30.0, seed 1234, on Entrez-collapsed log2(TPM + 1) expression rather than parsed from a metadata field. Left: stacked posterior probability, the four classes summing to 1, with the arg-max call and top probability printed. Middle: top minus second probability, dashed line at the median margin 0.160, band from 0 to 0.10 marking the thin-margin region as a display convention rather than an inferential cut-off; `→ CLASS` marks the 2 of 15 calls that change when the input model set changes. Right: each model's within-HGSC stratum from Ward.D2 clustering with k = 3 on the 50-set Hallmark z matrix, which involves no microenvironment signal.}}

Recomputing the calls with the actual classifier, rather than transcribing labels from a free-text field with no version recorded, is a correction that **strengthens** the argument rather than weakening it:

- **5 of 15 calls have a top-versus-second margin below 0.10** (0.022, 0.044, 0.072, 0.078, 0.088); the minimum top probability is 0.318 and the median margin 0.160.
- **2 of 15 calls change when the input model set changes** (`OV2085` IMR → PRO, `TOV3041G` MES → DIF).
- One call disagrees with the previously inherited label (`OV1369-R2`, inherited DIF, computed MES).
- **The immunoreactive and mesenchymal classes are called in 8 of 15 models**, and those two classes are defined by immune and stromal compartments that a pure culture does not contain.

**Interpretation.** These calls describe transcriptional similarity to a tumour-derived centroid; they are not tumour subtype assignments for a cell line, and a third of them sit within 0.10 of a coin flip. Report probabilities and margins, never bare labels. The tumour-cell-intrinsic strata of §10.3 are the appropriate resolution for within-HGSC structure in this panel.

### 9.5 Two of three mucinous models read gastrointestinal-leaning and have no external provenance

Marker standing across the three mucinous models in both layers (Fig. 5E), z computed per gene across all 31 models per layer:

| Model | KRT7 z | PAX8 z | SATB2 z (rank) | CDX2 z (rank) | MUC5AC z | Ovarian index | Verdict |
|---|---|---|---|---|---|---|---|
| **`TOV2414`** | **+2.58** | +0.28 | +0.17 (19/31) | +0.21 (25/31) | **+4.25** | **+1.26** | Ovarian-compatible |
| **`VOA8762`** | −1.55 | −1.15 | −0.59 (15/31) | **+2.52** (30/31) | −0.43 | −0.76 | GI-leaning |
| **`VOA8771`** | −0.98 | −1.23 | +0.80 (23/31) | **+3.64** (31/31) | −0.14 | −1.91 | GI / colorectal-leaning |

- **All three verdicts are stable across a 60-combination threshold sweep** in the sense that matters: every model stays on the same side of the ovarian/GI divide in 60 of 60 combinations. `VOA8771`'s *sub*-classification within GI moves in 12 of 60.
- **The composite index is ad hoc and its own deposited file says so**, verbatim: `mean(KRT7_z, PAX8_z) − SATB2_rna_z; ad-hoc unweighted composite mixing protein z (KRT7) with RNA z (PAX8, SATB2); uncalibrated`. The verdict thresholds are hand-chosen conventions. It should be described as ad hoc, never as calibrated.
- **External provenance differs sharply between the three, and that determines what each verdict can claim.** `TOV2414` has a Cellosaurus record (CVCL_A1SR) and a primary paper (Sauriol 2020) reporting KRAS G12A, so its ovarian call does **not** rest on expression alone, and KRAS G12A is confirmed here. `VOA8762` and `VOA8771` are **absent from Cellosaurus with no primary paper**, so for them expression is the only in-house check, and it can raise a flag but cannot call origin.

**Interpretation.** Genuine ovarian mucinous coverage in this panel may be **n = 1** (`TOV2414`) pending STR and histotype immunohistochemistry on the other two. `VOA8762`'s intestinal signal rests on *CDX2* alone.

One literature value must not be quoted as an in-house measurement. Sauriol 2020 reports `TOV2414` as SATB2-negative **by immunohistochemistry**, but its **measured** SATB2 mRNA z here is +0.17, rank 19 of 31: mid-panel, higher than `VOA8762` and lower than `VOA8771`. The literature value must not be quoted as an in-house measurement. `TOV2414` also carries a somatic *SMAD4* frameshift with low *SMAD4* mRNA; that co-occurs in a subset of mucinous ovarian carcinoma and, given the external STR record, is not evidence of pancreatic origin.

### 9.6 Only 8 of 42 models have a positive lineage programme supporting their label

The aggregate count must be decomposed rather than quoted bare, because the three kinds of evidence behind it are not equivalent.

**Expression consistency (all 42 models):** 26 consistent · 2 partial · 1 discordant · 13 not assessed.

Of the 26 consistent calls:

- **8 rest on a positive lineage programme**, with specific markers expressed;
- **16 rest only on the *absence* of a competing programme**, which is not a positive call for the model's own histotype;
- **2 rest on the same SWI/SNF evidence used to reclassify them**, which is circular.

By basis rather than by verdict, 10 models have a positive-lineage basis (8 consistent + 2 partial) and 17 an absence-only basis (16 consistent + 1 discordant), with 13 not assessed (11 have no RNA; 2 have no positive marker panel available).

**Genomic consistency (the 23 models with copy number):** 19 consistent · 3 partial (`OV90`, `TOV3121D`, `TOV3392D`) · 1 discordant (`TOV112D`) · 19 not assessed. Every one of the 19 consistent models is either *TP53*-mutant HGSC, or *TP53*-wild-type with a low-FGA genome matching its histotype (`TOV81D` LGS 0.021; `TOV21G` CC 0.071), or an intermediate-FGA histotype-appropriate genome (`TOV2414` MC 0.321).

**Interpretation.** **Expression-supported coverage of a histotype is smaller than its raw model count, for every histotype in this panel.** Only 8 of 42 models have a positive lineage programme supporting their label. That is not a criticism of the panel. It follows from two-dimensional culture, from small histotype groups, and from the absence of a positive marker panel that survives culture for several histotypes. It is nonetheless the number a reuser needs when deciding how much weight a label carries.

Three models carry specific flags:

- **`OV90`** is best read as "HGS-family carcinoma, serous identity not confirmed by expression": PAX8, WT1, SOX17 and KRT7 all sit near zero, the lowest serous programme of the HGSC models, alongside a *SMAD4* nonsense variant and the lowest HGSC FGA. `OV3331`, the other Adenocarcinoma-vs-HGS conflict, is the clean HGS.
- **`VOA5436`**, labelled MMMT, expresses a strong clear-cell programme (HNF1B / NAPSA, z = +2.04) and is the one discordant expression call.
- **`VOA4841`**, labelled clear cell, has the lowest *SMARCA4* of all 42 models in both layers and is called SWI/SNF-deficient, which is atypical for its histotype and worth verifying.

---

## 10 · Resources for reuse

### 10.1 Master per-model table (42 × 55)

`output/supplement_per_line.csv` is **42 rows × 55 columns**, one row per model, carrying metadata, assay coverage, per-layer quality metrics, patient-family membership, tiered driver calls, autosome-restricted fraction of genome altered, authentication verdicts **and their basis**, the recomputed subtype call with its margin, and STR record status. It is the recommended first stop for model selection.

### 10.2 Gene-level browser

`app/ovcan_viewer_standalone.html` is a single self-contained file (3.47 MB, no server or installation) that returns transcript and protein abundance across the panel for any gene symbol, coloured by the same histotype palette used in these figures, with honest not-detected and no-RNA states and flags for label conflicts. It covers 28,901 RNA symbols and 8,427 protein features across the 32-model union, and a view can be linked directly with `?gene=SYMBOL`. It is deployed at `cooklab.ca/ovcan_viewer`, so **the full transcript and protein matrices are already public**. That needs reconciling with the deposition and embargo plan.

Gene-identifier caveats are surfaced in the interface itself: 1,835 multi-mapping symbols are resolved by maximum mean TPM, and roughly 17% of Ensembl identifiers carry no symbol and are dropped.

### 10.3 Within-HGSC pathway strata

Ward.D2 clustering with k = 3 on the **full 50-set Hallmark z matrix** across the 15 HGSC models with RNA:

| Stratum | n | Models |
|---|---|---|
| Inflammatory / NF-κB-EMT | 4 | OV1369-R2, OV1946, OV866-2, OV90 |
| Low-signalling | 6 | OV2085, OV2295, OV2295-R2, OV3133-R, OV4485, TOV3133G |
| Hypoxic-glycolytic | 5 | OV3331, OV4453, TOV1369, TOV3041G, TOV3291G |

{{figure: figs8 | **Within-HGSC pathway strata.** (**A**) The two Hallmark theme axes used to *name* the clusters, with assigned stratum by colour, n = 15 HGSC models. (**B**) The proliferation theme, for which no independent pathway comparator exists. (**C**) PROGENy pathway activities for the same 15 models, grouped by stratum. The partition comes from clustering the full 50-set Hallmark z matrix, so a model's position on the naming axes cannot test its own label — the partition is not circular but the labels are.}}

**The clustering is independent of the axes the strata are plotted against, but the stratum names are not.** The strata come from clustering the full 50-set matrix, while those axes are the theme means used to *name* the resulting clusters, so a model's position on them cannot test its own label.

Checking the labels against PROGENy on the same expression matrix (a different method, **not an orthogonal assay**) gives a split verdict:

- **Inflammatory is corroborated**: `inflammatory_z` versus PROGENy NF-κB, Spearman ρ = **0.818**, p = 3.0 × 10⁻⁴ (and ρ = 0.825 against TNFα).
- **Hypoxia is not corroborated.** The correlation is only moderate (ρ = 0.632, p = 0.014), and **the inflammatory stratum's median PROGENy hypoxia score (+0.639) exceeds the hypoxic stratum's (+0.561)**. Only 2 of 5 hypoxic-stratum models exceed z > +1, and two are negative.
- **Proliferation has no comparator**, because PROGENy has no proliferation pathway.

**Interpretation.** Treat the labels as descriptive shorthand for cluster membership, and do not describe PROGENy as corroborating the strata as a set. The membership assignment is deposited in `hgs_heterogeneity.csv` and is usable for model selection on pathway state; the *names* should not be quoted as pathway claims. Two of the three subline pairs present fall in the same stratum and one does not, which gives a direct reading of how much within-patient variation the panel contains.

### 10.4 Antibody–drug-conjugate target expression atlas

{{figure: fig6 | **Antibody–drug-conjugate target atlas.** (**A**) Nine clinically pursued targets across the 30 models measured in both layers, RNA above and protein below, with models per histotype block printed above each block — so any non-HGSC statement rests on 2–6 models. **Row-scaled: each cell is z per feature across the 30 models, separately within each layer**, which is necessary because the median cross-model protein IQR is 0.34 log2 and an absolute-scale protein heatmap would be nearly uniform. Absolute information is retained in the right-hand `IQR log2` column, which gives how many log2 units one z unit is worth. `DPEP3` RNA has an interquartile range of 0.00, so a single model accounts for its entire row signal; a pale tile with × means not detected. (**B**) Every HGSC model's `FOLR1` value, RNA above and protein below, vertical jitter only, n = 15, separate x scales because the units differ. (**C**) Protein-to-RNA spread ratio for the eight targets measured in both layers under three statistics, n = 30; open circle = range ratio, filled circle = IQR ratio, triangle = SD ratio.}}

Nine targets under clinical development, tabulated in both layers on the 30 dual-layer models (Fig. 6A). Three properties determine what the atlas can be used for:

- **Coverage is not uniform.** `DPEP3` is essentially unexpressed: median exactly 0 log2 TPM, **25 of 31 models at zero**, maximum 4.59, and no protein measurement at all. `ERBB2` varies across less than a 1-log2 protein window (12.92–13.90). **The panel therefore contains neither a DPEP3-high nor an ERBB2-negative model.**
- **Protein interquartile ranges are 4–10× narrower than transcript interquartile ranges for every target measured in both layers** (3.84× for CD276 to 9.64× for MSLN). This is the compression of §6.2 seen target by target, and it has a concrete consequence: **protein-level *ranking* of targets is possible; protein-level *thresholding* is not.** Lead shortlists with RNA and use protein as rank-consistency confirmation.
- **Expected histotype associations hold target by target, not as a set.** `MSLN` → HGSC holds in both layers (top mean and top single model). `ERBB2` → clear cell / mucinous holds in both, with the top histotype swapping between layers on nearly identical protein means. **`FOLR1` → HGSC does not hold**: the highest histotype mean is mucinous on n = 3 in *both* layers, with HGSC third, and even the weaker version of the claim is assay-dependent (the top single model is HGSC in RNA but clear cell in protein).

The `FOLR1` distribution within HGSC (Fig. 6B) is easy to describe wrongly. It spans nearly the full expression range (0.057–9.479 log2 TPM, range 9.42, SD 3.09) and a forced two-component fit does split a small low group (4 of 15) from a larger high group. **But the formal tests do not support two components**: the bimodality coefficient is 0.4399 against a 0.555 threshold, and both mclust parameterisations prefer one component (ΔBIC −1.06 and −0.95). Hartigan's dip test was not run, the `diptest` package being absent from this R installation. **Describe it as a wide range with a subset of near-zero models, show the distribution, and do not use the word bimodal**, because bimodality would imply a natural cut point for patient-selection reasoning and the data do not establish one. Three of the four low models sit below 0.35 and the fourth at 1.88 sits in the gap, which is why a forced fit finds something while the formal tests prefer one component.

Several agents with high target expression had negative pivotal trials, so this atlas is a **model-selection aid, not clinical target validation**. Expression is necessary, not sufficient.

### 10.5 Recommended defaults for reusers

| Choice | Default | Why, with the number |
|---|---|---|
| Protein feature filter | `present_n_lines ≥ 16` → **7,733** features | The deposited accounting file calls this the recommended reuse filter |
| Protein features to exclude outright | the **70** `zero_plex = TRUE` rows | Identified but never quantified; retained only to preserve the 8,427 row count |
| FGA metric | **`fga_auto_0.2`** (autosome-restricted) | chrX is a pooled-normal artefact: median chrX fraction altered 0.964, per-model log2 −1.23 to +0.83 |
| Spread-ratio statistic | **IQR ratio and SD ratio** | Floor-insensitive; median RNA range nearly doubles (3.38 → 6.28) when the zero floor is engaged |
| Protein precision statistic | **`sd_single` / `cv_pct_per_measurement`** | Per-measurement, after the √2 replicate-difference adjustment; `sd_diff_cv_pct_legacy` is deprecated |
| Cross-assay gene set | **n ≥ 10 models → 7,894** genes | Median ρ moves only 0.397 → 0.408 across thresholds 10 → 30 |
| Cross-assay discordance filter | **`q_value < 0.05`**, not the sign of ρ | 813 genes have a negative point estimate; **30** are BH-significant |
| Mutation denominator | **`line_family_map.csv$has_wes_maf`** (22) | `samples.csv$wes_mut` has been corrected to 22, but the derived flag remains authoritative |
| Unit for **inference** (markers, pathways, DE validation) | **28 patient representatives** | 13 subline genomes come from 5 patients; label exchangeability holds only on the collapsed set |
| Unit for **description** (browsing, model selection, PCA display) | **31 models**, with the 28-representative sensitivity alongside | Both are deposited; the structure conclusion strengthens under collapse |
| Protein-change strings | only those flagged **`hgvsp_canonical`** (5,682 of 6,036) | All strings are reconstructed; 103 non-canonical, 251 not derivable |
| Environment citation | **`renv.lock`**, with `output/package_versions.csv` as the per-package record | `session_info.txt` records only the last script sourced |
| Execution order | **`scripts/run_all.sh`** | Numeric script order is not dependency order |

---

## 11 · Limitations

Twenty-five limitations, ordered by how much they constrain reuse. The first six (1–5b) change which analyses are possible; the rest change how results should be worded. Two entries carry letter suffixes because they were added at the methodological audit of 2026-07-27 and the established numbering is cited elsewhere.

| # | Limitation | The number |
|---|---|---|
| 1 | **Whether a single protein in a single model can be interpreted depends on that protein's abundance.** There is no single dataset-level answer. | Per-measurement bridge SD **0.149–0.190 log2** (10.4–13.2% CV) against a median cross-model protein SD of **0.287 log2**: observed spread **1.5–1.9×** technical, so bridge variability corresponds to **27–44%** of a typical protein's observed cross-model variance under a homoscedastic independent-error approximation. **Abundance-matched, that runs 67.5% (decile 1) to 7.5% (decile 10)** — approximate variance ratio 0.33 → 0.92. ADC targets, on their own technical estimate (32 bridge differences), **≈5%**. These are variance-ratio diagnostics, **not** intraclass correlations. Supersedes the retired "noise is 2.4–3.1× biological spread", which compared a 3.92-SD interval width with a 1.35-SD IQR |
| 2 | **Protein missingness is structural, not random.** Set-conditional presence must never be read as biology. | Exactly one distinct feature count per isobaric set (7,738 / 7,758 / 7,562 / 7,694 / 7,653); 8.796% NA; **1,572** features absent from ≥1 whole set; **70** never quantified |
| 3 | **Copy number is total, not allele-specific.** LOH, allelic imbalance, biallelic status and genomic HRD cannot be derived. | No BAF, no `--vcf` passed to CNVkit, tumour BAMs not archived; `scarHRD` correctly absent |
| 4 | **Variant calls are tumour-only against unmatched pooled normals of unconfirmed capture concordance.** Use tiered per-patient driver calls, not a burden metric. | `n_depth` / `n_ref_count` / `n_alt_count` empty in every MAF row; median coding burden 206.5 is an **upper bound** on residual germline, not an estimate; 12 of 50 driver calls are Tier 3 |
| 5 | **The layers are not passage-matched.** | Differences up to **20 passages in both directions** (max +20 `TOV112D`, min −17 `OV3331`); median \|Δ\| 4; Spearman between the two records **−0.251** |
| 5b | **Histotype and contributing centre are not separately identifiable in the expression layers.** No analysis in this resource can establish that histotype-associated structure is biological rather than a source-, culture- or processing-linked effect aligned with histotype. | All 15 HGS RNA models from CHUM and the rare histotypes concentrated at BC Cancer, so no design cell has the two varying independently; PC1 shared variance component **0.311**; the only within-histotype cross-centre comparison is clear cell at **2 versus 5** models and is uninformative in both directions (permutation p = 0.49 / 0.54 / 0.25) |
| 6 | **The RNA matrices cannot be reproduced from Ensembl 105 alone, and the omitted transcripts are not evenly distributed with respect to the design.** A design-alignment diagnostic, not a demonstration of harm — the association is equally consistent with genuine histotype differences in expression of the omitted set. | **3,529 of 185,299** index targets (1.90%) dropped silently, carrying median 2.22% of TPM (range 1.60–3.41%) and 1.53% of counts; dropped-TPM fraction R² = **0.330** on centre (p = 0.0033) and **0.388** on histotype (p = 0.010). **No matched-index sensitivity analysis has been run** — see §12 item 5 |
| 7 | **Sublines cost power in every contrast, not only the one that contains them.** | Collapsing to 28 representatives: five of six groups lose **no** members, yet CC (unchanged at n = 7) loses **435 of 1,256** DE genes (−34.6%) |
| 8 | **The exome layer is entirely single-centre and predominantly HGSC.** No centre comparison is possible in the genomic layer. | 23/23 CNV and 22/22 mutation models from CHUM; 18 of 23 are HGSC; all 13 tri-omic models from CHUM, 9 of 13 HGSC |
| 9 | **Arm-level frequencies are HGSC-only and threshold-dependent.** | Unit **18 models / 11 patients**, not 23/16; the six headline frequencies span **18–46 pp** across a ten-combination sweep (19q worst: 27–73%) |
| 9b | **Copy number is relative, not absolute, and pooled-normal capture compatibility is unconfirmed.** Coarse gain/loss and FGA under a declared threshold only. | Five public normals from an unrelated BioProject; `intervals_sorted.bed` not archived; per-sample autosomal median centring removes an unknown baseline; the 30-fold HGS-vs-LGS FGA contrast is **18 models against 1** and is a panel range, not an effect estimate |
| 10 | **Signature exposures are approximate.** No exome-to-genome renormalisation, because the capture kit is unknown; only within-model ranking carries. | `exome_renormalisation` = "none applied — capture kit not recoverable"; 303–2,417 SNVs per model |
| 11 | **A tenth of the proteome has an inverse transcript–protein point estimate, but almost none of it is established.** Use the q-value, not the sign — and read the survivors as inverse *associations*, not as failures of proxy value. | 813 of 7,894 genes (10.3%) inverse; **50** nominally significant, **29** (0.37%) BH-significant, against 3,599 with an FDR-supported positive association; most negative ρ = −0.714; median per-gene ρ 0.397. Conventional asymptotic Spearman test, df = n − 2 |
| 12 | **Two histotypes are transcriptomically absent or effectively n = 1.** | LGS: **0 RNA, 0 protein** models (exome only, n = 1). EC: 2 models, one reclassified, leaving 1 whose hormone-receptor markers read exactly 0.000 |
| 13 | **Only 8 of 42 models have a positive lineage programme supporting their histotype label.** | 16 rest on the absence of a competing programme; 2 are circular; 13 not assessed |
| 14 | **Molecular subtype labels are unreliable on pure cultures.** | 5 of 15 margins < 0.10; 2 of 15 change with the input set; minimum top probability 0.318; IMR + MES called in 8 of 15 |
| 15 | **No in-house STR profiling and no mycoplasma testing.** STR records are third-party depositions. | 30 of 42 with a documented profile; 25 of 26 named sources are one personal communication; 12 of 42 absent from Cellosaurus |
| 16 | **Two mucinous models have no external provenance at all**, so expression can flag but not call origin. | `VOA8762`, `VOA8771`: `external_provenance` = NONE; both read GI-leaning |
| 17 | **Formal differential expression is defensible for one histotype only.** | `evidence = formal` for HGS (n = 15); `descriptive` for the other five contrasts (n = 2–7) |
| 18 | **Small-group effect sizes are not comparable to large-group ones.** | Wilcoxon p floor 0.00529 at n = 2 and 0.000611 at n = 3; 15 of 25 marker CIs cross zero; widest MUC5AC *d* = 1.92, CI [−0.44, 20.64] |
| 19 | **Silhouette separation is modality-dependent and sometimes reverses.** | SCCOHT 0.819 (RNA) vs 0.028 (protein); CC +0.121 vs −0.004; EC −0.014 vs +0.128 |
| 20 | **One of the four isobaric bridge samples is external to the resource.** | `VOA3993`: Carey LGS, proteomics-only, identity to any RNA sample unconfirmed; carries the set 3↔4 link |
| 21 | **The centre difference in gene detection is real; its mechanism is not established.** | +948 genes (+4.8%), p = 0.0082; the two depth differences point the right way but are not significant (p = 0.37, 0.77) |
| 22 | **Six protein feature counts and three peptide totals legitimately coexist.** | 8,430 / 8,427 / 7,733 / 6,855 / 8,212 / 7,894, plus 8,357 and 6,688; peptides 146,830 / 146,629 / 158,885 |
| 23 | **The 652 MB external reference set is not regenerable offline.** | `output/external/`: 55 files, 652 MB (DepMap 24Q4 + Cellosaurus); the fetch script is verify-only unless `--download` is passed |

---

## 12 · Open items

Each item below closes a sentence that currently limits a section, and two would change what the resource can claim. They are grouped by consequence rather than by effort.

**Would change the science**

1. **STR profiles and mycoplasma clearance**, for all 42 models or at minimum the 12 with no Cellosaurus record. This is the highest-value outstanding item. It would invert the identity section from the resource's principal limitation into its strongest, stacking in-house STR on 5/5 DepMap self-matches at rank 1 of 67 and on molecular profiles consistent with the published `TOV112D` reassignment and the established SCCOHT identities. Priority models: `VOA8762`, `VOA8771` (no external provenance at all), plus `VOA5436` (MMMT versus clear cell), `VOA4841` (atypical *SMARCA4* loss), `VOA4395` (endometrioid markers) and `OV90` (serous identity). Histotype immunohistochemistry on the two mucinous models would shorten §9.5 to a single sentence.
2. **The exome capture kit and `intervals_sorted.bed`.** Load-bearing beyond documentation: it resolves the pooled-normal capture-concordance question for the entire copy-number layer, and it enables exome-to-genome renormalisation of the mutational signatures, currently the stated reason the `TOV21G` result can only be ranked rather than attributed.
3. **MMR immunohistochemistry, MSI-PCR and germline testing on `TOV21G`**, to upgrade a well-supported candidate call to a confirmed one.
4. **A decision on genomic HRD**: recover the exome BAMs and run Sequenza → scarHRD for a real score, or state that HRD is out of scope for this record. There is no third option that involves a number.
5. **RNA re-quantification against a transcript index matched to the deposited transcript-to-gene map.** **This is the sensitivity analysis, and it has not been done.** What §2.1 provides is a design-alignment diagnostic: the dropped-TPM fraction is associated with contributing centre (R² = 0.330) and histotype (R² = 0.388), which flags that the mismatch sits on the same axis §4 interprets. That association is equally consistent with genuine biological differences in expression of the omitted transcripts, so it neither demonstrates harm nor rules it out. The absolute quantity is small and bounded (1.60–3.41% of TPM) and every model shares the same reference pair, so the expected effect is modest — but only re-running kallisto against a matched index, or rebuilding the map from the index's own annotation release, can settle whether the QC, PCA, marker and pathway results are invariant. Requires the RNA FASTQ files, whose survival should be confirmed alongside the exome reads.

**Would complete the record**

6. **TMT method parameters and raw mass spectra from the Morin laboratory:** instrument, acquisition mode (MS2 versus SPS-MS3, which bears directly on the compression interpretation), LC and fractionation, search engine and database, FDR thresholds, rollup, any normalisation preceding the pooled-internal-standard ratio, the definition of the vendor CV, and the complete channel map for all five sets. Raw spectra are also required for PRIDE deposition.
7. **Culture and ethics metadata from the contributing centres:** base medium, serum, supplements and incubator conditions per model per centre, plus REB protocol numbers and consent framework. Expression phenotype depends on culture conditions, so this is required metadata for a cell-line resource rather than an optional annotation.
8. **Sarek / Mutect2 versions and the RNA library-prep kit and read configuration**, recoverable from a surviving Nextflow `pipeline_info/` directory and from the sequencing core respectively. **Confirmation that exome FASTQ files survive.** Tumour BAMs are known not to be archived, and if the reads do not survive, that record is processed-calls-only.

**Deposition hygiene**

9. Strip the bundled **third-party normal exomes** (SRR4039087/88/89/96/97 and a 461 MB BAM from PRJNA339046) from the CNVkit outputs before deposit, or deposit only our tumours' derived `.cnr` / `.cns`. Another group's data should not be redistributed here.
10. Decide where the **browser** is hosted so it can be cited as a named record (Zenodo DOI and/or GitHub Pages), and **reconcile the already-public `cooklab.ca/ovcan_viewer` deployment with the deposition and embargo plan**, since the full transcript and protein matrices are currently reachable without an accession.
11. If the **raw vendor MAFs** are ever deposited, correct the `NCBI_Build` header on a copy rather than in the archive.

---

## Appendix A · Figure inventory

| Figure | Content | Script |
|---|---|---|
| Fig 1 | Resource overview: workflow, per-model coverage matrix, histotype/patient counts | `30_fig1_overview.R` |
| Fig 2 | Quality across three assays: RNA QC, compression, protein accounting, CV by abundance, bridge agreement | `34_fig2_qc.R` |
| Fig 3 | Expression structure and biology: PCA by histotype and centre, concordance, variance components, marker recovery | `35_fig3_biology.R` |
| Fig 4 | Genomic fidelity and identity: filtering waterfall, DepMap specificity, driver oncoprint, per-arm CNV | `31_fig4_genomics.R` |
| Fig 5 | Rare-histotype and flagged models: burden, SBS-96 spectrum, cosine screen, refit exposures, mucinous markers, SWI/SNF | `32_fig5_rare.R` |
| Fig 6 | Target-expression atlas, FOLR1 within HGSC, per-target compression | `36_fig6_adc.R` |
| Fig S1 | Proteomic PCA and per-histotype silhouettes in both layers | `37_supp_rnaprot.R` |
| Fig S2 | Passage sensitivity and cross-assay passage discordance | `37_supp_rnaprot.R` |
| Fig S3 | Model-to-model RNA correlation and the centre-wise split | `37_supp_rnaprot.R` |
| Fig S4 | Protein presence patterns across isobaric sets | `37_supp_rnaprot.R` |
| Fig S5 | Genome-wide copy-number profiles with FGA | `33_supp_genomics.R` |
| Fig S6 | Recomputed molecular subtype calls with posterior probabilities | `33_supp_genomics.R` |
| Fig S7 | Marker effect sizes with bootstrap intervals | `37_supp_rnaprot.R` |
| Fig S8 | Within-HGSC strata and PROGENy pathway activities | `37_supp_rnaprot.R` |

Figures are embedded as **vector SVG** (converted from the `cairo_pdf` originals), so every panel is resolution-independent in both the HTML and the PDF; `reports/assets2/` also holds the source PDFs and 400 dpi PNGs, and `reports/assets2/panels/` the single-panel exports. Figure text is drawn as paths, so it renders identically without Arial installed but is not selectable in the PDF. These are the manuscript figures with the in-panel methodological footnotes suppressed (`OVCAN_FIG_PLAIN=1`, handled in `scripts/00b_figure_theme.R`); the geometry, palettes and statistics are identical, and the narrative here carries what those footnotes said.

## Appendix B · Where the numbers live

| Question | File |
|---|---|
| Which models, which patients, which assays | `metadata/samples.csv`, `metadata/line_family_map.csv` |
| Everything per model, in one row | `output/supplement_per_line.csv` (42 × 55) |
| RNA QC and the centre comparison | `output/rna_qc_metrics.csv`, `rna_qc_site_comparison.csv` |
| Reference-mismatch cost, and whether it is design-aligned | `output/rna_reference_reconciliation.csv`, `rna_reference_sensitivity.csv` |
| Protein feature accounting and per-feature annotation | `output/prot_feature_accounting.csv`, `prot_qc.csv`, `README_protein_matrix.md` |
| Protein precision, and technical noise against observed spread | `output/prot_bridge_agreement.csv`, `prot_noise_vs_biology.csv`, `prot_cv_by_abundance.csv` |
| Protein-versus-RNA spread ratio | `output/prot_dynamic_range.csv`, `prot_compression_floor_check.csv` |
| Structure, variance and confounders | `output/rna_pc_confounder_joint.csv`, `rna_pc_confounder_permutation.csv`, `rna_within_cc_site.csv`, `rna_variancepartition.csv`, `prot_variancepartition.csv`, `variancepartition_sensitivity.csv`, `silhouette_by_modality.csv` |
| Subline-collapse sensitivity | `output/sensitivity_patient_reps*.csv` |
| Differential expression and enrichment | `output/rna_de_*.csv`, `rna_de_gsea_recovery*.csv` |
| Marker effect sizes and the recovery null | `output/rna_marker_effectsizes.csv`, `rna_markers_summary.csv`, `rna_marker_recovery_permutation.csv`, `rna_ec_markers.csv` |
| Cross-assay concordance, with uncertainty and donor sensitivity | `output/integ_rnaprot_cor*.csv`, `integ_rnaprot_n_thresholds.csv`, `integ_rnaprot_negative_genes.csv`, `integ_rnaprot_patientrep_sensitivity.csv` |
| Variants, tiers and burden | `output/wes_mutations_filtered.csv`, `wes_driver_tiers.csv`, `wes_driver_freq_patient.csv`, `wes_mutation_load.csv` |
| Copy number | `output/wes_cnv_segments.csv`, `wes_cnv_fga.csv`, `wes_cnv_arm_freq_patient.csv`, `wes_cnv_arm_freq_sensitivity.csv` |
| Signatures and MSI | `output/wes_sbs_context.csv`, `wes_sbs_cosine.csv`, `wes_signature_refit_*.csv`, `wes_msi_mmr.csv` |
| Identity and authentication | `output/external_selfmatch_margin.csv`, `external_depmap_driver_crosscheck.csv`, `cellosaurus_str_status.csv`, `auth_swisnf_panel.csv`, `auth_perline_table.csv`, `consensusov_calls.csv`, `auth_mucinous.csv` |
| Reuse layer | `output/adc_expression.csv`, `adc_subtype_summary.csv`, `adc_folr1_bimodality.csv`, `hgs_heterogeneity.csv` |
| Environment and reproducibility | `renv.lock`, `output/package_versions.csv`, `scripts/run_all.sh`, `reports/integration_report.md` |
