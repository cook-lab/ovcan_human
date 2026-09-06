# Independent review of the ClinVar annotation snapshot

Reviewed 6 September 2026. This review used the saved NCBI responses and local tables only; no external requests were made. The production ClinVar VCV 2.6 schema and official classification/review-status documentation had been checked separately before the lookup. The final script SHA256 is `a21ee951626afb3a16351ac1bff2f9968dc4705daaff01594320f5672700f8eb`, matching the offline replay recorded at `2026-09-06T16:47:06.397436+00:00`.

**Outcome: no material identity, parsing, joining or interpretation defect remains in the reviewed outputs.** The annotations are suitable as a dated research annotation layer. They do not establish model-specific somatic origin, biallelic status, functional deficiency or ovarian treatment sensitivity.

## Independent checks

- Recomputed the URL-key hashes, response SHA256 values and byte counts for all 130 saved responses. All matched the request manifest. The approved queue contains exactly 53 unique coordinate-only alleles, and its SHA256 matches `validation.json`.
- Rechecked all 53 input SPDI strings against the saved GRCh38.p14 assembly accession report, all canonical-normalization responses, and all exact-field search results. The 23 matched alleles each resolve to one current, directly classified `SimpleAllele` with matching VariationID and CanonicalSPDI. Independently, **all 23 also match the exact GRCh38 chromosome accession, VCF position, reference and alternate tuple** in the XML. No haplotype or included-record classification was assigned to a component allele.
- The remaining 30 searches return no exact indexed hit, with the expected missing-phrase/no-items messages and no unknown-field error. This is a search outcome, not evidence of benignity, novelty or absence from every ClinVar representation.
- Checked all 69 framework-specific aggregate blocks against XML, including descriptions, review status, evaluation dates, submission counts, and contributing versus other condition sets. All match. Germline classifications are present for all 23 alleles; somatic clinical impact is absent for 16 and oncogenicity is absent for 13. Missing values remain blank.
- Checked all **226 classified SCV assertions** and **152 RCV records** against their XML source. The 232 total SCV records include six without a classification in these three frameworks; they remain available in the raw XML and are appropriately absent from the classified-assertion table. Flagged and noncontributing classified submissions are retained with their status.
- Checked the join to all **65 model-candidate rows**: 29 rows from 14 patients receive exact annotations, and every original local source field is unchanged. Related models sharing an allele are not additional independent clinical assertions. The explicit local-haplotype caveat is present for the TOV2881EP TP53 compound-event components.

## Biological conclusions checked

| Model / allele | Verified current annotation | Required interpretation |
|---|---|---|
| TOV2414 ATR p.Arg1653Ter | [VCV000224559.6](https://www.ncbi.nlm.nih.gov/clinvar/variation/224559/): germline Pathogenic/Likely pathogenic, multiple submitters without conflicts. The two contributing assertions were evaluated in 2024 and 2025; the older VUS assertion is flagged and noncontributing. | This updates the archived annotation. It does not establish ATR-null function, biallelic inactivation or HRD. Preserve the inherited-disease condition context. |
| TOV21G PTEN frameshift | [VCV000092828.45](https://www.ncbi.nlm.nih.gov/clinvar/variation/92828/): germline Pathogenic, multiple submitters without conflicts. The curated title is `NM_000314.8:c.800del (p.Lys267fs)`. | Preserve the source shorthand `p.L265fs` separately; the exact genomic match resolves the naming discrepancy. AF 0.329 does not establish biallelic loss. |
| OV90 BRAF p.Asn486_Pro490del | [VCV000666267.3](https://www.ncbi.nlm.nih.gov/clinvar/variation/666267/): Oncogenic, criteria provided by one submitter. The separate legacy germline field is `drug response`, without assertion criteria. | The legacy SCV concerns dabrafenib response in pancreatic adenocarcinoma. It is not an ovarian response label or an inherited pathogenicity assertion. |
| TOV1369 NF1 p.His1348Arg | [VCV000184987.11](https://www.ncbi.nlm.nih.gov/clinvar/variation/184987/): conflicting germline classifications, specifically VUS (two) versus Benign (one). | There is no pathogenic submission in this conflict. Do not promote the candidate to a pathogenic NF1 model classification. |
| TOV1369 / OV1369-R2 TP53 p.Gly244Cys; OV3291 TP53 p.Arg249Trp | [VCV000376599.19](https://www.ncbi.nlm.nih.gov/clinvar/variation/376599/) and [VCV000141881.15](https://www.ncbi.nlm.nih.gov/clinvar/variation/141881/) retain germline conflicts. Gly244Cys independently has a Likely oncogenic somatic assertion. | Preserve conflicts and framework separation rather than selecting only pathogenic submissions. |
| TOV2881EP TP53 p.Arg156Pro | [VCV000634768.14](https://www.ncbi.nlm.nih.gov/clinvar/variation/634768/): germline Pathogenic/Likely pathogenic for the isolated SNV. | The local call is phased with a nearby three-base insertion. Its combined consequence still requires joint read/haplotype review; the isolated-SNV classification is not transferred to that combined allele. |

The seven aggregate somatic clinical-impact descriptions are all Tier I, but the associated tumour and diagnostic/therapeutic contexts differ. The contributing-condition split and retained SCV/RCV attributes prevent a highest-tier statement from being applied to every listed condition. The README correctly gives the KRAS examples and retains this limitation.

## Reuse notes

The SCV `conditions` text field is blank for 151 assertions because their submitted XML does not include a preferred condition name in that location. This is not missing condition evidence: `condition_xml` retains identifiers, and the RCV table and aggregate condition sets provide additional context. Downstream summaries should use those sources rather than treating a blank display name as an unspecified disease.

Research tiers remain untouched and must not be renamed as clinical pathogenicity tiers. In particular, the expert-panel Pathogenic BRCA2 allele in TOV81D remains compatible with uncertain tumour/germline origin and unestablished biallelic status. Public annotation also cannot repair low read support, determine the combined ATM event, or resolve the pending cluster read-review and allele-specific copy-number tasks.

Reference semantics: [ClinVar classification definitions](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/), [review-status rules](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/), and [production VCV 2.6 schema](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/ClinVar_VCV_2.6.xsd). Accession versions above refer to the cached snapshot, not a claim that unversioned live pages will remain unchanged.
