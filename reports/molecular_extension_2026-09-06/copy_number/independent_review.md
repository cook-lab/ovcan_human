# Independent copy-number extension review

6 September 2026. Reviewed scripts 46, 47 and 47b and current source/output records without fetching references or rerunning alignments. After the read-only audit, made the two requested focused corrections to script 46 and reran it using the cached references. Scripts 47, 47b and 48, cluster documents, canonical scientific matrices and manuscript files were not edited.

## Corrections and validation

- Replaced `suppressor_repair_locus_screen` with neutral `other_locus_screen`; the previous label included CTNNB1 by list position. Grouping does not alter quantitative gain/loss calculations.
- Renamed event-count availability denominators to `models_with_WES` and `patients_with_WES`. Added `models_screen_eligible` and `patients_any_model_screen_eligible`, requiring complete segment-span coverage and at least three positive-depth/positive-weight targets. Added the row-level `relative_screen_eligible` field to preserve that distinction throughout the outputs.
- All 1,196 pre-existing WES model–gene measurements, their 988 non-WES companion rows, expression values, review flags and event flags are unchanged. The pre/post hash of every pre-existing `gene_model_evidence.csv` column except the intentionally renamed grouping is identical. There remain 41 gain-screen and 15 loss-screen rows.
- The corrected counts identify 1,074 eligible WES model–gene rows and 122 ineligible rows. Positive-event model and patient counts never exceed eligible counts, and eligible counts never exceed WES availability. Missing WES retains a missing eligibility field.
- `gene_target_bins.csv`, the source of scripts 47/47b, remains byte-identical. All their existing depth/exon tables remain byte-identical, so those scripts did not require rerunning.

## Independent numerical and source checks

The 2,184 model–gene rows contain exactly 52 genes x 42 distinct models, with the correct 34-patient mapping. There are 23 WES models from 16 patients. Current CNS intervals do not overlap within model/chromosome, and gene segment covered fractions do not exceed one. The 23,069 target-evidence rows have no duplicate model/gene/coordinate keys and contain no antitarget rows.

Directly reread original CNR files for OV3331/AKT2, both 1369 models/NF1 and CDKN2A, and TOV2835EP/TOV3121D/TOV3121EP/NF1: all 449 target rows match coordinates, raw depth and raw log2 values. Independently recomputed canonical-exon membership using chromosome-aware 0-based half-open overlap against the pinned exon records; all eight reviewed model–gene memberships agree. Canonical exons and bins are on their expected chromosomes (NF1 chr17, CDKN2A chr9, AKT2 chr19). This checks exon-overlapping target bins, not the number of affected exons or base-level coverage.

Reread each of the five primary normal targetcoverage CNN files and checked all 119 selected gene/bin records per normal, 595 records total. Every coordinate join succeeds, and depth values match exactly, with maximum absolute difference zero. These are source bin depths, not calibrated model/normal copy ratios.

| Finding | Independently supported result |
| --- | --- |
| AKT2 / OV3331 | Six positive targets, median centred log2 3.219126; gene-segment log2 2.964191. Five target bins overlap canonical exons. RNA 1,177.011454 TPM is rank 1/31; protein 12.851898 is rank 2/31, verified directly from canonical matrices. |
| NF1 / TOV2835EP | 66/87 gene-span targets at zero depth; 51/65 canonical-exon-overlapping targets at zero depth. |
| NF1 / TOV3121D | 54/87 gene-span targets at zero depth; 38/65 canonical-exon-overlapping targets at zero depth. |
| NF1 / TOV3121EP | 38/87 gene-span targets at zero depth; 28/65 canonical-exon-overlapping targets at zero depth. |
| NF1 / OV1369-R2 and TOV1369 | No zero-depth target among 87 gene-span or 65 canonical-exon-overlapping bins. Segment ratios -0.232240 and -0.516676, respectively, support the contrast with their low RNA abundance rather than a deep locus-depletion explanation. |
| CDKN2A / OV1369-R2 | Three of four gene-span targets at zero depth; two of three canonical-exon-overlapping targets at zero depth. Only one positive-depth bin, so ineligible for the strict three-positive-bin screen. |
| CDKN2A / TOV1369 | Two of four gene-span targets at zero depth; one of three canonical-exon-overlapping targets at zero depth. Only two positive-depth bins, so ineligible for the strict screen. |

## Interpretation boundaries

The AKT2 concordant multimodal lead and NF1/CDKN2A depletion-review findings are supported. Raw-depth denominators must stay explicit: the gene-target mean uses full overlapping target-bin widths, and canonical-exon summaries include any-overlap bins and UTRs of a single pinned transcript. CDKN2A's p16 transcript check does not establish the p14ARF transcript state. CNS span coverage is interpolated segment coverage, not base-level capture completeness.

The eligibility distinction is especially important for depleted loci. A positive-bin screen can deliberately exclude the most depleted regions, so its false flags cannot stand for validated absence of loss. Zero-depth review flags and the canonical-exon tables preserve the relevant evidence. The normal comparison supports capture feasibility; it does not remove preprocessing, mapping, structural or allelic uncertainty. No absolute copy number, homozygous deletion, protein loss or functional-deficiency call is established by this review.
