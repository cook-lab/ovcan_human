# Exact-variant evidence extension

**Update:** the author-approved [NCBI/ClinVar query](clinvar/README.md) is complete: 23 exact alleles / 29 model observations. Use that snapshot for current ClinVar fields; the manual tables retain pre-API review history and primary-paper context.

6 September 2026. This directory adds source read evidence and a separate interpretation layer to retained candidates. It does not change the canonical 6,194 variants, their research tiers, or manuscript v9.

**All 65 selected candidates match their original VCF alleles exactly.** They represent 53 distinct genomic alleles across all 23 WES models (16 patients). The original 36-gene molecular panel contributes 57 candidates; the approved extension adds eight candidates in CTNNB1, CDKN2A, NF1, ARID1B and SMAD4. The scan checked the hashes of all 23 VCF/MAF pairs and reconciled all 582,474 source record pairs. All selected records are PASS in both the VCF and MAF. Caller AF has been recovered for 24 indels whose earlier MAF-derived VAF field was unavailable.

## Files and reproduction

| File | Purpose |
| --- | --- |
| `variant_read_evidence.csv` | All 65 model/candidate records, with exact source alleles, AD/AF, depth, strand/orientation, local phase fields, original consequences/tiers/filters and file hashes. |
| `variant_transcript_consequences.csv` | Original VEP fields for each selected gene/candidate, including MANE identifier and coding/protein positions; no new transcript annotation is inferred. |
| `source_manifest.csv` | All 23 source VCF/MAF paths, SHA-256 hashes, record counts and selected-candidate counts. |
| `curated_assertions.csv` | Manually reviewed evidence for 30 distinct alleles. Includes exact public-record assertions where available, literature support, and explicit unresolved/read-review findings. This is not a table of 30 clinically classified alleles. |
| `curated_model_evidence.csv` | Exact genomic-allele join of those reviews to the current model observations. Original fields remain separate from review fields. |
| `public_evidence_sources.csv` | Public record IDs/versions and primary-source URLs used in this bounded review. |
| `compound_haplotype_review.csv` | Two local phase groups that should not be counted as independent hits. |
| `bam_review_candidates.csv` / `bam_review_regions.bed` | 25 flagged records and approximately 100-bp flanks for local/cluster BAM inspection. BED coordinates are zero-based, half-open; model labels select the corresponding BAM. |
| `external_annotation_request.tsv` | The approved, completed 53-allele query input; results and provenance are under `clinvar/`. Only build, chromosome, one-based VCF position, reference and alternate are present. |
| `validation.json` / `QA.md` | Reconciliation counts, input/script hashes, invariants and limitations. |

Run from the project root with standard-library Python:

```bash
python3 scripts/45_variant_read_evidence.py \
  --additional-genes CTNNB1 CDKN2A NF1 ARID1B SMAD4
```

`OVCAN_PROJ` selects the project; `OVCAN_DATA` can restore the archived data layout elsewhere. This is an offline extraction. It performs no web queries or variant calling. The manual curation CSV is an input to the final join and is hashed in the validation record. The script does not regenerate that manual assessment.

## What changed the interpretation

**OV2295 ATM needs joint haplotype review.** The eight-base deletion at VCF chr11:108345741 AGGAGGTGC>A and adjacent eight-base insertion at chr11:108345749 C>CAAAAAAAA share PID, PS and PGT. They have AD 181,7 and 174,8, caller AF 0.035/0.039, and alternate strand counts 0,7 and 0,8. If these events occupy the same haplotype, their combined net length change is zero. The deletion crosses a splice acceptor, so a splice effect remains possible, but the isolated insertion's frameshift annotation cannot be treated as an independent second ATM loss event. Joint realignment/haplotype and splice-junction review is the first priority; no ATM-deficient label is supported.

**TOV2881EP TP53 also has a composite representation.** Its three-base insertion and nearby R156P substitution share a local phase group and AD 0,24. Both have strong caller support. Review the combined transcript/protein consequence rather than count two independent TP53 hits. Their combined net length change is +3 bp. This affects representation and interpretation, not the existing model-level TP53-positive observation.

The remaining highest-value read checks include PALB2 C824S in OV3133-R (two alternate reads), BRCA2 A22E in TOV3133D (five alternate reads, AF 0.038, TLOD 5.10), CDK12 S325I in TOV2929D (three alternate reads, TLOD 4.72), and SMARCA2 P874fs in OV2295 (four alternate reads, all on one strand). These are candidates to inspect, not demonstrated false calls. Several receptor/MYC missense candidates also have only two to four alternate reads and should not drive model classifications.

## Supported biological leads

| Lead | Source evidence and interpretation |
| --- | --- |
| **TOV81D BRCA2** | AD 75,45; raw AD fraction 0.375 and Mutect AF 0.392. Exact allele matches ClinVar c.8537_8538del, p.Glu2846GlyfsTer22, expert-panel Pathogenic. The original p.R2845fs string is a repeat-position reconstruction and is retained as source provenance. Historical TOV81D work reported a retained, expressed wild-type allele; neither BRCA2-null nor present functional HRD is established. |
| **OV90 BRAF** | N486_P490del has AD 55,50, AF 0.481 and balanced alternate strands. The exact 15-bp coding deletion agrees with published OV90 sequencing and later deletion-peptide work. This is a useful non-V600 BRAF model lead, not a V600E-equivalent response category. |
| **KRAS / PIK3CA** | TOV21G KRAS G13C has AD 42,45; TOV2414 KRAS G12A has AD 0,84; TOV3392D KRAS G12S has AD 0,80. Exact public condition/oncogenicity assertions are recorded with their differing evidence scopes. TOV21G PIK3CA H1047Y has AD 43,29, an exact ClinVar match and primary gain-of-function evidence. TOV112D KRAS A59T has strong local support but no newly assigned clinical classification. |
| **ARID1A / SMARCA4** | Both TOV21G ARID1A frameshifts have two-strand support and roughly balanced allelic reads. Published TOV21G studies provide protein-position/model context; current trans phasing is unknown. TOV112D SMARCA4 L639fs has AD 0,57 and a source MANE early-frameshift consequence, supporting follow-up of the established SWI/SNF phenotype. |
| **CDK12** | The early L123 frameshift has AD 0,36–42 across four isolates from patient 3133. OV3291 has a strongly supported canonical splice-donor candidate (AD 0,72). Allelic predominance is clear, while absolute biallelic inactivation is unmeasured. CDK12-associated genomic instability is not interchangeable with BRCA-associated HRD. The P1256 in-frame deletion is not a second truncating hit. |
| **ATR / SMAD4 / ARID1B** | TOV2414 ATR R1653* has AD 28,26 and a source stop-gain consequence, the archived ClinVar string was uncertain significance, while the approved current lookup is Pathogenic/Likely pathogenic. Neither establishes ATR-null/HRD status. SMAD4 stop/frameshift candidates in OV90 and TOV2414 have strong support. ARID1B truncations are useful validation leads, with limited depth in TOV21G. |

The pre-API manual review was deliberately narrow: five alleles had public ClinVar nucleotide-level matches in that subset. The subsequent approved query found 23 exact alleles, documented in `clinvar/README.md`. The following paragraph describes the earlier selected-source evidence, not the current aggregate snapshot. Their assertions differ: the BRCA2 record is expert-panel pathogenic; KRAS G12S has a somatic oncogenicity assertion; the cited KRAS G12A/G13C records are condition-specific; PIK3CA H1047Y has a germline-disease aggregate plus separate primary functional evidence. These fields must not be flattened into one generic “clinically pathogenic” column. [BRCA2](https://www.ncbi.nlm.nih.gov/clinvar/variation/9328/), [KRAS G12S](https://www.ncbi.nlm.nih.gov/clinvar/variation/12584/), [KRAS G12A](https://ncbi.nlm.nih.gov/clinvar/RCV004549454.2/), [KRAS G13C](https://www.ncbi.nlm.nih.gov/clinvar/RCV000443868.1/), [PIK3CA H1047Y](https://www.ncbi.nlm.nih.gov/clinvar/variation/39705/).

OV90's source VEP coding interval is 1457–1471 on NM_004333.6; the reverse complement of its genomic deletion is ATGTGACAGCACCTA, matching the published coding deletion. Other frameshifts retain rough source labels until transcript-aware reference normalization is performed. [OV90 direct sequencing](https://doi.org/10.1371/journal.pone.0001279); [OV90 endogenous peptide/functional evidence](https://doi.org/10.1126/sciadv.ade7486). The historical BRCA2 and CDK12 distinctions are supported by [the TOV81D transcript study](https://doi.org/10.1158/1940-6207.CAPR-11-0547) and [the primary CDK12 genomic study](https://pubmed.ncbi.nlm.nih.gov/26787835/).

## Field definitions and review priorities

AD is the caller's reference/alternate allelic read support. `AD_alt_fraction_all` divides the selected alternate count by all AD entries; `caller_AF` is Mutect's estimated fraction and is not forced to equal AD fraction. All 65 selected records are biallelic VCF records, so the all-allele and reference-plus-selected fractions coincide here. FORMAT DP equals summed AD for these selected records. FAD reports fragments; F1R2/F2R1 report pair orientation; alternate forward/reverse counts come from the allele-specific strand table. They are different summaries and need not have identical denominators.

The VCF genotype, local PGT/PID/PS, allele fractions and zero reference reads do not prove tumour-normal origin, whole-gene copy number, distant trans phasing or biallelic loss. Source ClinVar strings date to VEP's archived 202404 resource; only the explicitly cited public records received an additional bounded review.

Heuristic flags are reproducible triage rules: alternate reads <5, AF <0.10, TLOD <10, all alternate reads on one strand when alternate support is at least five, indel alternate reads <10, a supplied STR flag, another selected call within 10 bp, or shared local phase identifiers. They are not new call filters and were not calibrated against truth. Four records belong to the two complex phase groups; 12 have other read-support/nearby concerns; nine are flagged solely for repeat context; 40 have none of these flags. A well-supported repeat indel, including TOV81D BRCA2, is not downgraded merely because it lies in a repeat. The four-read SMARCA2 event is flagged by low support even though its single-strand count does not meet the separate five-read strand-rule threshold.

Matching uses the exact archived MAF row joined to the corresponding VCF record, checked independently by common-prefix/suffix-trimmed genomic edits. This is not reference-based left normalization. The selected events all match uniquely without a reference lookup. Full HGVS normalization and combined transcript effects are separate work. No absence in this retained panel is a callable wild-type result.

## Prepared external lookup and cluster handoff

The author subsequently approved the coordinate-only query, and all 53 lookups completed. See [the current results](clinvar/README.md), which supersede pre-API and archived ClinVar classification strings. The original automatic approval block was resolved by this explicit authorization; no model names, read data or source paths were transmitted. Saved raw responses, URLs, dates and hashes now support offline replay.

The completed lookup does not replace local cluster work: use the BED and candidate CSV to select the corresponding BAMs, inspect the two compound haplotypes jointly, confirm low-support calls with strand/orientation/mapping and local-repeat context, and return compact per-candidate evidence. Normalize alleles against the exact GRCh38 reference and annotate combined events on versioned transcripts. Preserve the original PASS calls and report any revised candidate status in a separate layer. RNA junction evidence is particularly valuable for ATM and the OV3291 CDK12 donor candidate. This extraction did not produce a new whole-exome callset or assign biological loss categories.
