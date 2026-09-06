# Approved NCBI/ClinVar coordinate annotation

Completed 6 September 2026 after the author explicitly approved querying the prepared 53-allele coordinate-only queue. All queries completed successfully. **23 of 53 distinct alleles matched directly classified ClinVar records**, representing **29 of the 65 model–variant observations**. The remaining 30 alleles had no exact indexed Canonical SPDI match; this is not evidence of benignity or absence from every external resource.

## Findings that add to the earlier manual review

| Model / allele | Current exact ClinVar evidence | Interpretation |
| --- | --- | --- |
| **TOV2414 / ATR R1653\*** | [VCV000224559.6](https://www.ncbi.nlm.nih.gov/clinvar/variation/224559/): germline aggregate **Pathogenic/Likely pathogenic**, multiple submitters without conflict. MANE NM_001184.4:c.4957C>T. | Updates the earlier archived uncertain-significance string. The local call has AD 28,26; this does not establish biallelic ATR loss, ATR-null status or HRD. |
| **TOV21G / PTEN frameshift** | [VCV000092828.45](https://www.ncbi.nlm.nih.gov/clinvar/variation/92828/): germline aggregate **Pathogenic**, multiple submitters without conflict. MANE NM_000314.8:c.800del, preferred p.Lys267fs. | Exact genomic allele confirmed despite the rough source protein label p.L265fs. Strengthens variant-level annotation; retained reference reads and unknown absolute/allelic CN preclude a PTEN-null assertion. |
| **OV90 / BRAF N486_P490del** | [VCV000666267.3](https://www.ncbi.nlm.nih.gov/clinvar/variation/666267/): **Oncogenic** aggregate. The separate germline-classification XML field contains the legacy category **drug response**. | Adds exact database oncogenicity support to the existing primary-paper evidence. The legacy field concerns dabrafenib response in pancreatic adenocarcinoma; it must not be relabelled germline pathogenicity or an ovarian treatment-response category. V600E equivalence is not inferred. |
| **TOV81D / BRCA2 deletion** | [VCV000009328.54](https://www.ncbi.nlm.nih.gov/clinvar/variation/9328/): expert-panel germline **Pathogenic**; exact VCF and normalized repeat-aware SPDI matches. | Confirms the existing c.8537_8538del annotation. Historical retained wild-type transcript evidence and present biallelic/functional uncertainty remain. |
| **TOV2929D / ATM L2390del** | [VCV002796934.1](https://www.ncbi.nlm.nih.gov/clinvar/variation/2796934/): **Uncertain significance**, single submitter. | No ATM-deficient label follows from this in-frame deletion. This is distinct from the adjacent low-fraction ATM calls in OV2295. |
| **TOV1369 / NF1 H1348R** | [VCV000184987.11](https://www.ncbi.nlm.nih.gov/clinvar/variation/184987/): **Conflicting classifications of pathogenicity**, specifically uncertain significance versus benign. | There is no pathogenic submission in this snapshot. This missense allele is a different observation from the NF1 depletion seen in TOV2835EP and the 3121 family. |
| **TP53 R249W and G244C** | [R249W / VCV000141881](https://www.ncbi.nlm.nih.gov/clinvar/variation/141881/) in OV3291 and [G244C / VCV000376599](https://www.ncbi.nlm.nih.gov/clinvar/variation/376599/) in the 1369 family have conflicting germline aggregates. G244C separately has **Likely oncogenic** support. | Discordance across classification frameworks is retained; a pathogenicity-conflict label is not proof of a biologically neutral TP53 allele. |

The full 23-record snapshot contains 17 germline P/LP-category aggregates, three conflicts, two uncertain-significance aggregates and one legacy drug-response aggregate. Separately, seven alleles have **Oncogenic** and three **Likely oncogenic** assertions. Seven have a somatic clinical-impact aggregate, all labelled Tier I—Strong in the returned records. These counts overlap and are counts of distinct alleles, not independent patients or validated model phenotypes. The exact raw review status, record version, evaluation date and contributing contexts accompany every classification.

TOV2881EP's TP53 R156P SNV also has an exact ClinVar record, but the local VCF phases it with a nearby three-base insertion. The matched classification describes the isolated SNV; it is not transferred to the combined allele. The model table contains an explicit local-haplotype caveat, and the existing joint read-review request remains open.

Clinical-impact tiers particularly require context. For example, the retrieved KRAS G12A Tier I assertion concerns diagnostic support in pilocytic astrocytoma/high-grade astrocytoma with piloid features; G12S Tier I concerns diagnostic support in rhabdomyosarcoma. Other submitted contexts do not all contribute to the highest-tier aggregate. These must not become ovarian drug-response labels. The table separates contributing from other conditions and retains SCV/RCV assertion attributes. [ClinVar classification definitions](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/)

## Query and identity methods

Script 50 reads only the approved five-column queue when constructing project-specific network requests. It verifies GRCh38 chromosome accession versions against NCBI's GRCh38.p14 assembly report, converts VCF positions to zero-based SPDI and requests a normalized canonical representation through NCBI Variation Services. This accounts for repeat-equivalent indels instead of treating a nearby coordinate or protein label as an exact match. The scientific basis is described in the [NCBI SPDI publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC7523648/).

It then searches ClinVar using a **quoted complete Canonical SPDI** and fetches the corresponding VCV XML. Search results are accepted only when the returned archive and directly classified `SimpleAllele` have the expected VariationID and identical Canonical SPDI. Haplotype components and `IncludedRecord` classifications are not transferred to individual alleles. As an independent identity check, **all 23 matches also agree exactly with the direct allele's GRCh38 `positionVCF`, `referenceAlleleVCF` and `alternateAlleleVCF`**, including the original anchored BRCA2 repeat deletion. [ClinVar search documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/help/), [programmatic access](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)

The parser uses separate direct paths for germline classification, somatic clinical impact and oncogenicity in the current VCV format. It retains each framework's review status, evaluation date and counts, VCV accession/version, VariationID, AlleleID, MANE coding HGVS, assembly coordinates and source response hash. Raw XML is retained for audit. A record's update date is distinct from the clinical evaluation date. Missing classification fields remain unavailable, never benign or uncertain significance by default. The production schema and review semantics are documented by [ClinVar VCV2.6](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/ClinVar_VCV_2.6.xsd) and [review-status definitions](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).

Local model names and read evidence are joined **after** querying. No model/sample names, source paths, allele fractions or patient identifiers were transmitted. The approved request is an annotation query, not a submission of new clinical assertions to ClinVar. It changes neither the original callset nor manuscript/release data.

## Files and reproduction

| File | Contents |
| --- | --- |
| `allele_lookup.csv` | All 53 input alleles, input/canonical SPDI, exact-match/no-match status, query strings and response hashes |
| `matched_allele_annotations.csv` | The 23 exact records with separate aggregate frameworks and contributing/noncontributing condition sets |
| `model_annotations.csv` | All 65 model observations joined to query status and, when present, current ClinVar fields; original caller evidence preserved |
| `submission_assertions.csv` | 226 SCV framework assertions, including review status, aggregate contribution and clinical-impact attributes |
| `condition_assertions.csv` | 152 RCV condition records, retaining framework-specific XML and assertion attributes |
| `requests.json`, `responses/` | 130 production requests/responses: assembly report, 53 normalizations, 53 searches and 23 VCV fetches; URL, time, status, bytes and SHA-256 |
| `validation.json`, `independent_review.md` | Query/annotation denominators, identity checks, authorization and independent review |

The earlier manual `curated_assertions.csv`/`curated_model_evidence.csv` retain primary-paper evidence and the pre-API review history. **For current ClinVar classifications, use this directory's snapshot**, which supersedes archived/pre-API classification strings. No-match records can still have useful literature or local sequence evidence, and none should be labelled wild-type from this lookup.

Some submitted conditions are represented by database identifiers without preferred names: 151 SCV assertion rows therefore have a blank condition-name display field. Their identifiers remain in `condition_xml`, with condition-level context additionally retained in the RCV table; a blank display name does not mean no submitted condition.

```bash
# Offline replay from saved response snapshots; no network.
python3 scripts/50_clinvar_coordinate_annotation.py

# Only when intentionally retrieving missing responses:
python3 scripts/50_clinvar_coordinate_annotation.py --fetch
```

Existing response hashes are checked and the script does not refresh them silently. Network requests are sequential and throttled. This snapshot covers the specific 53-allele queue; it is not whole-exome annotation, full reference-based HGVS curation of unmatched calls, a reclassification of the research tiers, or an assay for somatic origin, second hits or drug sensitivity. The cluster's read/haplotype, exon-depth, allele-specific CN and MSI tasks remain useful; repeating this completed query is unnecessary.
