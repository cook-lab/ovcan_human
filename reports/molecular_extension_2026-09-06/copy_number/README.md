# Expanded relative gene copy-number screen

6 September 2026. Scripts 46, 47 and 47b extend the CCNE1/ERBB2 review to 52 cancer, repair and chromatin-associated loci. This is a targeted research screen across 23 WES models from 16 patients, with RNA/protein measurements joined by exact model. It does not assign absolute amplification, homozygous deletion, biallelic inactivation or HRD. All 42 resource models remain in the combined table.

## Findings worth following up

| Locus / model | Evidence | Interpretation and next step |
| --- | --- | --- |
| **AKT2 / OV3331** | Gene segment log2 ratio 2.964; six positive gene-overlapping bins, median 3.219 and range 2.575–4.038; bin median is 3.029 above the arm median. RNA 1,177.011 TPM, highest of 31; protein 12.852, second-highest of 31, without an isoDoping flag. | Strong concordant same-model DNA/RNA/protein lead. Confirm absolute copy number and amplicon boundaries before calling clinical amplification. Assays are not assumed to use the same passage or aliquot. |
| **CCND1 / patient 3133** | OV3133-R, OV3133-R2 and TOV3133G pass the gain screen, with three measured bins each. Segment log2 ratios 1.655–1.702; bin medians 2.152–2.270. | Three related models, one patient family. TOV3133D instead has a CCNE1 gain; inspect family-specific structures and passages. RNA/protein support differs among measured sublines. |
| **NF1 / TOV2835EP** | Gene segment −7.105; 66/87 gene-span bins have zero depth. Of 65 bins overlapping the pinned canonical transcript's exons, 51 have zero depth; the remaining positive-bin median is −7.843. All five reference normals have positive coverage at these canonical-exon bins. | Strong candidate for deep locus depletion, supported beyond a single segment or intronic targets. Determine affected exons, residual reads, mapping specificity and allelic/structural status from alignments. No RNA/protein assay is available. |
| **NF1 / TOV3121D and TOV3121EP** | Canonical-exon bins have zero depth in 38/65 and 28/65, respectively; positive-bin medians −6.403 and −6.534. Some bins retain substantial coverage. | Related models with complex/partial-locus depletion. TOV3121D's segment coverage spans only 75.5% of the annotated gene, so it fails the strict complete-span event screen. This is a review trigger, not a negative result or proof of whole-gene deletion. |
| **CDKN2A / OV1369-R2 and TOV1369** | All four gene-span bins have depth <0.9; respectively three and two are zero. The pinned p16 transcript overlaps three bins, of which two and one are zero. All five normals cover these bins, although normal capture is modest. RNA is low (0.264 and 0.038 TPM); supplied protein remains quantified. | Candidate locus depletion in two related models. Resolve p16 and p14ARF transcript structure and residual sequence; do not infer complete protein loss from the DNA/RNA pattern. |
| **CDKN2A / TOV2414** | Three of four gene-span bins have zero depth, while the remaining bin has depth 86.327. Two of three pinned p16-exon bins are zero. The broad segment is near baseline. | Partial-locus dropout warrants local review. It is not a complete-gene deletion call. A segment or positive-bin median alone conceals this pattern. |

The complete screen has **41 model–gene gain rows and 15 model–gene loss rows** meeting the declared thresholds. These are neither independent-patient counts nor a genome-wide prevalence estimate. Broad 3q or other arm events may affect several panel genes together. Gene/histotype denominators and numbers of unique patients with any positive model are provided in `patient_event_counts.csv`; repeated sublines must not be counted as independent patients. Missing WES and a failed screen are distinct from evidence of a normal locus.

The count table separates WES availability (`models_with_WES`, `patients_with_WES`) from strict-screen eligibility (`models_screen_eligible`, `patients_any_model_screen_eligible`). Eligibility requires complete segment-span coverage and at least three positive-depth, positive-weight targets. Of 1,196 WES model–gene rows, 1,074 meet those requirements; the remaining 122 are not eligible negatives. In particular, very deep depletion can leave too few positive bins to enter the strict screen. Among HGSC models, NF1 has 18 WES models from 11 patients but 16 eligible models from 10 patients; CDKN2A has the same denominator split. Patients count once if any of their models meets the stated condition. These counts support transparent review, not disease-event prevalence estimates.

The previous CCNE1 and ERBB2 results are reproduced exactly for all 46 model–locus combinations. In particular, TOV3392D remains an ERBB2 expression lead without a comparable gene-bin gain; the broad segment must not be used to relabel it as focal ERBB2 amplification. See the [earlier assessment](../../clinical_classification_2026-09-06/README.md).

## Methods

**Panel and coordinates.** `gene_loci.csv` lists 52 selected genes, coordinate conventions, stable IDs and sources. GRCh38 gene spans were pinned from the public [Ensembl lookup API](https://rest.ensembl.org/documentation/info/lookup_post), release 116; raw responses and hashes are saved. CCNE1 and ERBB2 retain the earlier pinned NCBI spans for continuity. These boundaries are used only to intersect DNA loci. RNA remains quantified with the matched Ensembl release-93 map. The ovarian [TCGA study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3163504/) provides biological context for recurrent copy-number loci; it does not supply this screen's thresholds.

The remainder outside the prespecified oncogene-gain grouping uses the neutral label `other_locus_screen`, including CTNNB1. Group labels do not assign every remaining gene a tumour-suppressor role and do not change the gain/loss calculations.

**Relative DNA summaries.** The current script-29 target-only segmentation, followed by script 08's autosomal centring, is used. Original CNR files are checked against saved SHA-256 hashes; antitarget bins are excluded. The same autosomal offset is applied to target log2 ratios. Gene segment values are overlap-length-weighted means. Segment min/max and covered fraction expose breakpoints or incomplete span coverage. Gene-bin medians are unweighted medians of bins with positive depth and positive weight; raw-depth and zero-depth summaries retain all gene-overlapping target bins. Arm medians are segment-length-weighted, and gene-minus-arm contrasts are descriptive local-versus-broad evidence, not focality calls.

**Membership matters.** A target bin intersects a locus if its 0-based half-open interval overlaps the annotated gene span. This can include intronic or nested-gene targets, and includes a whole bin's measured depth when only part overlaps the span. `gene_overlap_bp` records the actual overlap; the mean-depth field explicitly weights full target-bin widths. These are not exon-specific gene-copy estimates. Therefore NF1, CDKN2A and AKT2 received a separate canonical-exon intersection using pinned GRCh38 transcripts: NF1 ENST00000358273.9, CDKN2A ENST00000304494.10 and AKT2 ENST00000392038.7. Exons include UTRs; this check is not a CDS-only analysis or an exhaustive all-isoform review. In particular, the CDKN2A check refers to the pinned p16 transcript and does not collapse p14ARF into it.

**Screens.** A relative gain/loss row requires complete segment coverage of the gene span, at least three positive-depth/positive-weight bins, and both the gene segment and bin median ≥1 for gain or ≤−1 for loss. These analyst-defined ±1 log2 cutoffs prioritize review and are not absolute-copy or clinical cutoffs. Separate flags identify a segment/bin difference ≥0.5, an extreme bin median without a concordant complete-span screen, and at least half of three or more bins with zero depth. Zero-depth bins are never replaced with a finite pseudo-log ratio or allowed to disappear from the review.

`relative_screen_eligible` records the coverage/bin-count requirement separately from either positive-event flag; it remains missing where WES is unavailable. Complete segment-span coverage means the annotated gene is spanned by the interpolated CNS intervals, not that every base or exon has been directly sequenced. A false event flag on an ineligible row must not be interpreted as absence of the event.

**Normal coverage.** Script 47 verifies all five recovered normal target-CNN hashes and joins exact bin coordinates for seven selected loci. Script 47b intersects three of those loci with canonical exons. Positive normal coverage supports the assertion that depleted model regions can be captured; model/normal processing equivalence has not been established sufficiently to interpret raw-depth ratios as calibrated absolute copy number. Low normal coverage is retained and reported, particularly at CDKN2A and some AKT2 targets. All NF1 canonical-exon normal/bin measurements are positive, with the minimum 18.340×; the corresponding CDKN2A minimum is 1.227×.

**Expression.** Raw RNA TPM and supplied protein log2 normalised abundance are joined by exact model name using existing gene/protein mappings. Descending assay ranks use minimum ranks for ties among measured models. They are descriptive ranks across this collection, not population reference ranges or histotype-adjusted expression calls. Patient, recorded histotype, assay passage, TMT plex and isoDoping fields remain visible. Missing measurements are not imputed. Global patient/family and histotype sensitivity for other selected targets is handled in the [expression extension](../expression/README.md).

## Files and reproduction

From the repository root, using the project's Python numpy/pandas and R plotting environment:

```bash
python3 scripts/46_gene_copy_number_screen.py
python3 scripts/47_locus_depth_review.py
python3 scripts/47b_locus_exon_review.py
Rscript scripts/48_molecular_cnv_figure.R
```

Normal runs use pinned reference JSON without network. The explicit `--fetch-reference` and `--fetch-exons` options refresh those references and should not be used for an exact rerun. `OVCAN_DATA` can redirect original CNR inputs; `OVCAN_COVERAGE_BUNDLE` can select the extracted September coverage bundle. No original input is modified.

| Files | Contents |
| --- | --- |
| `gene_model_evidence.csv` | 2,184 model–gene rows (42×52), including 1,196 WES rows, raw expression, screens and denominators |
| `gene_target_bins.csv` | Individual gene-span target intersections with raw depth, weights and centred ratios |
| `supported_relative_events.csv`, `locus_review_queue.csv` | Strict research-screen positives and broader review flags |
| `patient_event_counts.csv`, `within_patient_contrasts.csv` | Explicit model/patient denominators and related-model contrasts |
| `selected_loci_model_normal_bins.csv`, `selected_loci_depth_summary.csv` | Seven-locus depth review against five reference normals |
| `canonical_exon_bin_evidence.csv`, `canonical_exon_depth_summary.csv` | Three-locus canonical-exon sensitivity and retained/excluded bins |
| `*validation.json`, `reference_provenance.json`, `*lookup.json`, `canonical_transcript_exons.json` | Input/output checks, reference provenance and pinned annotations |
| `independent_review.md`, `independent_review.json` | Source-level audit and verification that denominator/group-label corrections preserved all scientific values and 41/15 screen-positive counts |

## Figure legend

**Expanded copy-number and AKT2 evidence** (`output/pdf/molecular_copy_number_extension.pdf`). **A:** Twenty selected loci across all 23 WES models. Rows identify model and patient; several rows belong to the same patient. Colour represents the median centred log2 ratio of positive-depth, positive-weight gene-overlapping target bins, saturated below −3 and above +3 for display. It does not encode the strict event screen or absolute copy number. Crosses identify loci with at least three bins and at least half at zero depth; those zeros are excluded from the finite colour median and retained explicitly in companion depth tables. Grey denotes an unavailable finite median. Gene-span intersections can include noncanonical or intronic targets; canonical-exon follow-up is provided for NF1, CDKN2A and AKT2. **B:** AKT2 overlap-weighted gene-segment log2 ratio versus RNA TPM in 13 exact models with both measurements. **C:** AKT2 RNA versus supplied protein log2 normalised abundance in 30 exact models with both assays. RNA is plotted as log2(TPM+1), with ticks labelled in TPM. OV3331 is highlighted in rust. Values are descriptive; related models and passage differences remain, and no correlation p-value or treatment-response inference is made.
