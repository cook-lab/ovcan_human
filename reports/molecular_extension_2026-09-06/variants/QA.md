# Variant extension QA

Final offline extraction completed successfully on 6 September 2026 with bundled Python 3.12 and the standard library. See `validation.json` for the exact script and input SHA-256 hashes.

| Check | Result |
| --- | --- |
| Source VCF and MAF hashes | All 46 files match the pre-existing input manifest; gzip file bytes and uncompressed file bytes are labeled separately. |
| Archived record linkage | All 582,474 VCF/MAF record pairs reconcile in count and chromosome order. |
| Selected allele identity | 65/65 uniquely match the corresponding minimal genomic edit after common-prefix/suffix trimming; 53 unique VCF alleles. |
| Scope | 57 candidates in the core 36 genes, plus eight in the five approved existing driver genes; 23 models from 16 patients. |
| Source filters | All 65 selected records retain VCF PASS and MAF PASS. No source filter or research tier was changed. |
| Read-field completeness | AD, AF, DP and TLOD available for all selected records; caller AF recovered for 24 indels. All selected VCF records have one alternate allele. |
| Internal read counts | Summed AD equals FORMAT DP for all 65 records; AD-derived fraction and caller AF are separate. |
| Earlier HRR audit | All 13 independently extracted HRR records reproduce exact AD and numerically identical AF/TLOD. |
| Curation linkage | 30 unique manually reviewed allele rows join by exact chromosome/VCF position/ref/alt to 35 model observations; gene identity also matches. |
| Compound groups | Two groups: OV2295 ATM (net length zero if the calls share a haplotype) and TOV2881EP TP53 (net +3 bp). Local phase fields support review, not final reconstructed haplotypes. |
| Review triage | Four complex-group records, 12 other read-support/nearby concerns, nine repeat-only records, 40 without these heuristic flags. |
| External queue | Exactly 53 deduplicated rows and five permitted fields; no model/patient identifiers, source paths, allele fractions or annotations. The author subsequently approved the query; all 53 lookups are complete, with results and response hashes under `clinvar/`. |
| Preservation | Script and all saved input hashes were independently rechecked after the final run. Canonical mutations, tiers and metadata remain unchanged. |

No raw BAMs or matched normals were available to this workstream, so no new caller sensitivity, truth-set performance, strand-bias significance test, somatic origin, biallelic copy state, clinical validity or treatment response is claimed. Simple triage cutoffs do not establish a false-positive call. Reference-based left normalization and transcript-aware annotation of joint haplotypes remain cluster follow-up work.

Public curation is deliberately bounded: five selected alleles have exact nucleotide-level public ClinVar matches; BRAF has exact published coding-deletion evidence; other reviewed alleles are clearly labeled as source consequences, published model context or unresolved findings. Record versions and sources are retained in the curation/source tables. Complete, independently hashed online record snapshots were not downloaded. The attempted broader coordinate query was rejected by automatic approval review; no retry or indirect coordinate lookup followed, and the prepared queue requires explicit authorization before external submission.
