# WES completion and manuscript v8

The cluster handoff substantially improves the paper: WES acquisition, computational processing, capture definitions and direct sequencing QC are now described and supported by checked records. Manuscript v8 and Figures 2/4 incorporate the findings. **The paper still needs specific laboratory and author records beyond data deposition.** No additional discovery analysis is needed to address those gaps.

## What changed

| Area | Verified result and manuscript change |
| --- | --- |
| Acquisition | Provider worksheets describe HiSeq 4000 paired-end 100-bp reads, TruSeq LT adapters and Roche NimbleGen SeqCap EZ Exome v3. Eighteen records are from run 4316 (2017-08-11), five from run 4706 (2018-05-25). Twenty-two rows match model/passage names; `TPV81D_23-pool` remains an explicitly flagged candidate mapping to TOV81D. The 2014 acquisition was not substituted. |
| Workflow | The completed `angry_allen` run on 2025-10-23 records Sarek 3.5.1, Nextflow 24.10.2 and 38,541 successful tasks. The retained annotated VCF headers link to this run. The recorded 32-character revision field is not a verified Git commit. |
| Read processing | Recovered commands support BAM-to-read conversion, fastp 0.23.4 with adapter trimming disabled but filtering retained, BWA-MEM 0.7.18, GATK 4.5.0.0 duplicate marking/BQSR and tumour-only Mutect2 filtering, and VEP 113/cache 113. The variant panel of normals remains distinct from the five CNV-reference exomes. |
| Capture | Vendor 242,232 targets reconcile with 242,215 mapped source intervals plus 17 unmapped intervals; split mappings yield 242,421 GRCh38 intervals. Their sum is 63,709,951 bp and union is 63,514,049 bp. The 290,475 CNV target bins cover that same union. The signature analysis used the standard-chromosome subset, 63,465,385 bp. |
| Variant filtering | 582,474 annotated records include 19,816 VCF PASS records. MAF conversion flags 3,735 of these as `common_variant`, leaving 16,081 MAF PASS records. The additional population filter leaves 15,995; coding/consequence selection leaves the same **6,194 candidates**. Figure 4 now says **MAF PASS**, and Methods describe all population-filter stages. |
| WES technical validation | Independently parsed all 23 models' read counts, primary alignment/proper-pair fractions, duplicate fractions, target depth/coverage and GATK contamination estimates. Figure 2 adds depth/coverage panels G/H. No aggregate-QC rule justified dropping a model. |
| CNV interpretation | All 69 recovered segment files match the retained originals. Reference/bin-quality filtering accounts for 85,769 excluded target bins; 204,706 target bins remain in each CNR, of which 203,913–204,064 have positive depth. This is distinct from the full capture-BED denominator used for sequencing QC. In every model, 42,705 of 42,709 antitarget bins have zero depth; segment probe totals match positive-depth bins, consistent with their exclusion. Profiles are effectively target-supported. Segments interpolate between targets, so FGA is a segment-span summary and does not measure every intergenic base. Existing relative profiles and patient summaries are preserved with this qualification. |

## WES QC added to the paper

All summaries below use 23 models, not independent patients. Fractions are converted to percentages where labelled. Full precision, model-level values, denominators and source links are in the companion tables.

| Metric | Median | Range |
| --- | ---: | ---: |
| Read pairs after fastp | 61.6 million | 54.2–68.4 million |
| Primary mapped reads | 99.95% | 99.83–99.98% |
| Properly paired primary reads | 99.09% | 98.57–99.25% |
| Picard duplication | 19.8% | 14.9–24.0% |
| Mean target depth | 78.8× | 69.7–88.7× |
| Target bases at least 10× | 95% | 91–97% |
| Target bases at least 20× | 88% | 79–92% |
| Target bases at least 30× | 78% | 66–84% |
| GATK estimated contamination | 0.27% | 0–0.49% |

Depth uses the genome-wide duplicate-marked stage and mosdepth's default flag mask 1796: duplicate reads are excluded and overlapping mates are corrected. Its BED denominator is the **sum of interval lengths**, including overlapping contributions. Coverage fractions are rounded to 0.01 in the source, so exact threshold-covered base counts cannot be reconstructed. TOV2929D has the lowest depth/coverage and remains included; negative calls require locus-specific coverage checks. High mapping rates and low estimated contamination do not prove unique mapping, stock identity or uniform variant sensitivity.

Recalibrated alignments are interval restricted and cannot supply genome-wide on-target/mapping denominators. Their summaries are retained separately with stage-specific denominators. Numerical agreement with the contemporaneous MultiQC report is not proof that later rewritten alignments are byte-identical to the calling inputs.

## Corrections to the supplied handoff

- The handoff generator missed all 23 Picard files because they are nested one sample directory deeper than its expected path. The new parser recovers and independently checks those duplicate metrics.
- The handoff description of duplicate-inclusive mosdepth coverage was incorrect for the executed options.
- TOV3121D's purported uncompressed VCF hash equals the local **compressed file** hash. The other 22 reported hashes match uncompressed streams. The corrected audit retains both hash scopes and does not claim to have independently read absent cluster VCF bytes.
- The five normal depth values remain **reported-only**: no primary normal coverage CNN files were included. They are not presented as independently validated manuscript results.
- Source BEDs, interval-weighted denominators and standard-chromosome footprints are now distinguished explicitly.
- The original converter version for 22 MAFs remains unknown. Its frequency-flag behaviour is reproduced across all 582,474 rows by the vendored vcf2maf default; the version is established only for the local TOV3121D reconstruction.

## What remains, apart from deposition

1. **WES laboratory/provider:** DNA extraction/quality and complete library preparation; confirm TOV81D's provider alias and what was pooled. TruSeq LT adapter metadata is not a full protocol.
2. **WES cluster:** small five-normal coverage files and existing manual CNV execution/input evidence. Recover original converter commands, workflow/reference identities and liftover details where retained. The [targeted cluster request](../../docs/cluster/recovery/2026-09-05/FOLLOWUP.md) replaces the previous broad search list.
3. **Other laboratory methods:** actual culture/harvest conditions and aliquot relationships, RNA extraction/integrity, and proteomics acquisition/search/normalisation/CV definition.
4. **Stock records and authors:** stock-linked STR/mycoplasma records, ethics/consent determination, prior use of these exact datasets, authorship/contributions/interests/funding and AI-use verification.
5. **Publication presentation:** resolve remaining TODOs, finalise the author reference-manager/numbered-citation conversion, and verify any claim about the hosted browser. The public code URL is filled in; its final licence/version/archive still needs an author decision.

These requirements are detailed in [the readiness review](manuscript_readiness.md) and [v8 author confirmations](../../docs/manuscript/v8/author_confirmation.md). Optional HRD, allele-specific CNV, MSI or additional signature analyses are not necessary for the present Data Descriptor and are not substitutes for these records. If an unrecoverable historical detail is appropriately reported as a limitation, it need not become an endless search task.

## Reviewable outputs and reproduction

- [Manuscript v8 source](../../docs/manuscript/v8/OvCAN_Scientific_Data_draft_v8.md) and [Word draft](../../docs/manuscript/v8/OvCAN_Scientific_Data_draft_v8.docx) · [review PDF](../../docs/manuscript/v8/OvCAN_Scientific_Data_draft_v8.pdf); v7 is preserved.
- [Figure 2](../../docs/manuscript/figures/fig2.pdf), [Figure 4](../../docs/manuscript/figures/fig4.pdf), [all figures and legends](../../docs/manuscript/figures/README.md).
- [QC model table](../../output/wes_qc_model_summary.csv), [metric summaries](../../output/wes_qc_metric_summary.csv), [coverage curves](../../output/wes_qc_coverage_profile.csv), [QC source hashes](../../output/wes_qc_sources.csv); definitions and checks in [QC validation](qc_validation.md) and [JSON](qc_validation.json).
- [Acquisition records](../../output/wes_acquisition_records.csv), [cell-level worksheet provenance](../../output/wes_acquisition_sources.csv), [summary](../../output/wes_acquisition_summary.json).
- [Run/filter/interval audit](provenance_audit.md), [machine-readable provenance](../../output/wes_recovered_provenance.json), [model reconciliation](../../output/wes_recovered_provenance_models.csv), [CNV support](../../output/wes_recovered_provenance_cnv_support.csv).
- Builders: `scripts/23_wes_recovered_provenance.py`, `scripts/24_wes_recovered_qc.py`, `scripts/25_wes_acquisition_records.py`. They inspect existing records without executing recovered code or launching analyses on the cluster.

The archive SHA-256 is `42badc44c0195ef0919cf99650e8c4e9ed67d4bb9960fa496e75256490441f58`; it contains 969 regular files (499,092,422 uncompressed bytes). It was extracted only after rejecting unsafe paths and non-regular file types. The original tar and extracted provider/run records remain ignored by Git. No cohort, variant candidate, RNA/protein matrix or existing release checksum was changed. New QC companions are ready for review and later deposition packaging; the existing 49-file release was deliberately not reorganised.

## Final validation

All 19 Word/PDF pages were visually inspected, including the four embedded main figures. A blank page after the references was removed by attaching page breaks to headings. Figures 2/4 were rendered separately and checked; their Arial fonts are embedded, and all 14 individual figure PDFs remain single-page exports without PDF annotation objects. The combined PDF and figure ZIP were regenerated. The full variant script passed its pinned content check, all prior WES result CSVs and metadata were unchanged, and all 49 existing processed-release checksums passed. See [validation manifest](final_validation.json) and [script execution checks](script07_validation.md).
