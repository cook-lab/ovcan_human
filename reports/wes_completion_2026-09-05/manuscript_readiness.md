# Manuscript readiness after WES recovery

Review date: 5 September 2026. Baseline: manuscript v7, its author-confirmation document, `docs/PROJECT_STATUS.md`, the September analysis audit, and the cluster WES handoff and recovery report. This review does not alter v7. The parent workstream is preparing v8 and independently reconciling the recovered records. New WES facts below distinguish the handoff narrative from checks reported by the acquisition, QC and provenance reviewers.

The acquisition record is materialised in [model-level provider records](../../output/wes_acquisition_records.csv), [source-workbook provenance](../../output/wes_acquisition_sources.csv) and the [acquisition summary](../../output/wes_acquisition_summary.json). These retain the source sheet/row, source-file checksum and the unresolved TOV81D mapping status. The legacy [project status](../../docs/PROJECT_STATUS.md) describes pre-recovery requests and should be updated alongside v8 rather than used to re-open facts recovered here.

## Readiness judgment

**The paper is not yet complete apart from data deposition.** The WES recovery substantially closes the missing computational methods and adds direct sequencing-quality evidence. It does not supply the actual culture/harvest, RNA-extraction or proteomics protocols, stock-linked authentication and mycoplasma records, ethics determination, or agreed authorship/declarations. These are substantive completion items. WES DNA preparation and one provider pool/alias also require clarification. More exploratory biology is not the remedy for these gaps.

The existing Data Descriptor structure and restrained interpretation are appropriate. Keep the four main figures centred on resource composition, assay quality, expected expression patterns and genomic validation; the planned addition of WES depth/coverage panels addresses a real weakness. Extensive subtype discovery, pathway interpretation, signature fitting and therapeutic speculation should remain outside the main descriptor. The journal assesses technical quality and completeness, and does not require new biological hypotheses or extensive follow-up experiments. [Scientific Data guidance for referees](https://www.nature.com/sdata/policies/for-referees)

## Remaining non-deposition completion items

| Item | What is still needed | Completion criterion |
|---|---|---|
| **M01: culture, provenance and harvest** | Actual conditions used for each centre/model, assay-specific passage/stock, relation between RNA and protein aliquots, VOA derivation sources, OV3291 status and the specific TOV3121D derivation. Published oxygen and serum conditions differ. | Laboratory-approved methods and model metadata, with unavailable historical details identified explicitly. A related paper is a source to check, not evidence that its protocol was followed here. |
| **M02: RNA extraction and quality** | Extraction protocol; integrity metric, values or available summary, acceptance criterion and pre-library checks. Library preparation/sequencing details already supplied by the author should remain. | An actual extraction/QC account. The original kallisto index/build record is useful if retained, but the exact release-93 transcript identifiers and lengths have already been reconciled; a new quantification is not required solely because an index binary is unavailable. |
| **M03: proteomics acquisition and processing** | Generating laboratory's sample preparation, LC–MS acquisition, search/database, identification filters, reporter/channel processing, protein summarisation, pooled-standard composition and normalisation formula; definition of the supplied CV. | A verified protocol and processing account matching these files. The Orlando et al. protocol remains a proposal for laboratory review. Its two-unique-peptide rule cannot be copied onto a table containing lower-support features. |
| **M04: residual WES acquisition** | DNA extraction/quality, library construction protocol and relation of the provider's TOV81D pool/alias to the analysed stock. | Provider/laboratory confirmation; do not equate a TruSeq LT adapter annotation with a complete library-preparation protocol. Instrument, paired-end read length, run dates and capture design have now been recovered. |
| **M05: residual CNV provenance** | Recover the manual CNV execution log and five per-normal coverage files, or document unavailable evidence; the pooled reference is already local. The profiles are now established to be supported almost entirely by target bins. | Retain the target-support/interpolation qualification and finish the source-record audit. This is not a renewed question about whether the authors intended to use the same capture kit. That compatibility is author-confirmed. |
| **M07: ethics and consent** | Source-tumour approvals/consent, current institutional determination and permitted genomic sharing. | A record-based statement and Human Data Checklist. Do not substitute approval identifiers from derivation papers without checking their applicability. |
| **M09: stock testing** | STR/authentication and mycoplasma method, date, result and relationship to the profiled material. | Report what was actually done. If historical records cannot be recovered, state that limitation rather than declaring the tests were never performed or claiming authentication from a database match. |
| **A02: prior use of these exact data** | Which RNA, protein and WES samples/data files appeared in earlier publications, and whether those exact data already had published accessions. | An explicit previous-use/first-release account in Background & Summary, with corresponding primary citations. This is a suitability question, not merely reference formatting. |
| **A01, A03–A06 and M08: authors and declarations** | Agreed author list/order/affiliations, contributions, interests, funding, acknowledgements if applicable, and the final AI-assistance disclosure and author review. | Author-approved statements. These are not inferable from repository history. |
| **D02/D03: claims about public tools** | Verify the hosted browser serves the corrected data, or revise the claim to the available local/archive artifact; identify the public working code and intended final revision. | The text must describe an actually available version. The working repository exists; its URL should no longer be an unknown. Licence/immutable archive decisions remain separate. |

M06 is a lower-priority historical-documentation item: cached Cellosaurus entries, their entry versions and checksums, and the documented 23 July 2026 snapshot already support reproducible use. If the database release/original HTTP records were not retained, record that fact and close the request without inventing a release. DepMap 24Q4 Public version 1 is already resolved and cited.

Scientific Data requires sufficiently complete acquisition and processing methods and disclosure of earlier publications using these data. Its current prior-publication policy makes substantial overlap with previously disclosed datasets especially important to establish. The extent of overlap is presently unknown, so A02 cannot be treated as a cosmetic TODO. [Submission guidelines](https://www.nature.com/sdata/submission-guidelines), [related-publication policy](https://www.nature.com/sdata/policies/editorial-and-publishing-policies)

The Human Data Checklist applies to human-derived datasets. The institutional assessment must determine the appropriate sharing conditions; human origin alone is not a sufficient reason to assert that controlled access is required. [Scientific Data human-data policy](https://www.nature.com/sdata/policies/data-policies)

## WES revisions for v8

### Acquisition and computational methods

Replace the broad M04/M05 retrieval placeholders with verified facts and a short, precise residual request. The acquisition reviewer has parsed provider spreadsheets for 18 current samples from run 4316 dated 11 August 2017 and five from run 4706 dated 25 May 2018. These identify Illumina HiSeq 4000 paired-end 100-base sequencing, TruSeq LT adapters and the Roche NimbleGen SeqCap EZ Exome v3 design. Keep the TOV81D `TPV81D_23-pool`/`TOV81D_P23+_-2` relationship unresolved until its pool composition and model/passage interpretation are established. Do not extend the cohort to the additional historical exomes in the handoff.

The revised method should describe the actual sequence of operations:

1. Paired reads were recovered from provider GRCh37-aligned BAMs and realigned to the GRCh38 reference used by nf-core/Sarek 3.5.1, running through Nextflow 24.10.2. This distinction explains why an original alignment assembly and the current analysis assembly differ.
2. fastp 0.23.4 performed read preprocessing with adapter trimming disabled and a minimum read length of 15 bases. Quality, ambiguous-base and length filtering still operated. **Do not copy the handoff's “no trimming” as a claim of no read filtering.** Any other preprocessing behaviour should be taken from the executed command/log, not inferred from a current pipeline default.
3. BWA-MEM 0.7.18 alignment was followed by GATK 4.5.0.0 duplicate marking and base-quality recalibration. Marking duplicates is distinct from removing their records. State the reference/resource identifiers and link the archived versions, commands and checksums through the processing records; do not crowd the article with cluster paths or job identifiers.
4. Tumour-only Mutect2/FilterMutectCalls 4.5.0.0 used a population allele-frequency resource, a variant panel of normals, contamination estimates and orientation-bias modelling. The variant panel of normals is distinct from the five healthy exomes used as the CNV reference. No matched normal samples were used.
5. Annotation used Ensembl VEP 113.0/GRCh38.p14. The population-frequency flagging in the MAFs is reproduced by vcf2maf 1.6.22, which was used with `--inhibit-vep` to reconstruct TOV3121D. The original converter command/version for the other 22 MAFs was not retained in the handoff; do not imply that its historical version was independently established. Preserve the exact database versions verified from VCF headers and add the converter-associated population-frequency filtering to the method.
6. Manual CNVkit 0.9.10 processing supplied the retained copy-number profiles and superseded the flat-reference Sarek CNV output. Preserve the established relative-copy-number and chromosome-arm method with the target-support qualification below. Do not imply that a default Sarek CNV run produced the published segments.

Explain that the original hg19 vendor target design was lifted to GRCh38 and filtered for the analysis. Distinguish vendor intervals, calling/coverage intervals and derived CNVkit target/antitarget bins. The coverage interval sum is 63,709,951 bp; its merged union is 63,514,049 bp. The full CNVkit target-bin union is the same 63,514,049 bp; its standard-chromosome subset used for exploratory signature sensitivity is 63,465,385 bp. These are different footprints, not interchangeable versions of a single “exome size”. The provenance table should carry the interval counts, exclusions and checksums.

Later CRAM re-executions and the historical successful VCF-generating run must be distinguished in the file catalogue. Configuration agreement does not establish byte-identical historical alignments. Likewise, no unverified full reference-FASTA checksum, container digest, liftOver-chain version or manual-CNV command option should be supplied from inference.

### Variant-filtering description: material correction

The present v7 claim that 16,081 records passed “Mutect2 filters” is inaccurate. The provenance reviewer traced **19,816 VCF PASS records to 16,081 MAF PASS records after vcf2maf added `common_variant` to 3,735 records**. This is a filtering change during conversion, not simply a distinction between VCF and MAF counting units. Subsequent explicit population-frequency filtering retains 15,995 records and the coding/splice filter retains 6,194. The retained 6,194 candidates are unchanged by this clarification.

Suggested main-text wording:

> The converted tumour-only variant tables contained 582,474 records. Of these, 16,081 retained PASS status after caller filtering and the population-frequency filter applied during VCF-to-MAF conversion. Subsequent population-frequency and consequence filters retained 6,194 coding or splice candidates across 23 models. The individual filter states are supplied with the candidate records.

Use “MAF PASS” or an equally explicit label in Figure 4a and its legend. Include the 19,816-to-16,081 reconciliation in the detailed source table/methods; do not silently call the entire difference caller-quality filtering.

### Technical validation

The QC reviewer independently parsed all 23 primary samtools/Picard/mosdepth/contamination records and the fastp records embedded in MultiQC. These support direct WES validation in the main text and the new Figure 2 panels:

| Quantity across 23 models | Minimum | Median | Maximum | Interpretation |
|---|---:|---:|---:|---|
| Read pairs after fastp | 54,191,475 | 61,594,469 | 68,390,777 | The handoff's doubled numbers count R1 and R2 reads separately. |
| Mean target depth | 69.72× | 78.76× | 88.75× | mosdepth 0.3.8, duplicate-excluding; denominator is the sum of supplied BED interval lengths. |
| Target fraction at ≥30× | 0.66 | 0.78 | 0.84 | Archived fractions are rounded to two decimals. |
| Target fraction at ≥20× | 0.79 | 0.88 | 0.92 | Do not reconstruct exact covered-base counts from rounded fractions. |
| Target fraction at ≥10× | 0.91 | 0.95 | 0.97 | These are analysis-interval coverage fractions, not whole-genome fractions. |
| Picard duplicate fraction | 0.1491 | 0.1977 | 0.2400 | Report as 14.9–24.0% if expressed as percentages. |
| Estimated contamination fraction | 0 | 0.002668 | 0.004885 | Maximum 0.4885%; distinguish fractions from percentages. |

mosdepth used the standard FLAG 1796 exclusion, MAPQ ≥0 and mate-overlap correction in the verified command configuration; contrary to the handoff narrative, duplicate reads were excluded. Overlapping BED intervals contribute repeatedly to the interval-weighted depth summary. Provide the exact computation/denominator in the QC dictionary and a concise explanation in Methods or the legend. Genome-wide alignment percentages should come from the duplicate-marked alignment records, not the interval-restricted recalibrated CRAM statistics.

Suggested concise validation text:

> Across 23 models, mean depth over the analysis intervals was 78.8× at the median (range, 69.7–88.7×), and the median fraction of interval bases covered at least 30× was 0.78 (range, 0.66–0.84). Duplicate fractions ranged from 14.9% to 24.0%, and estimated contamination was below 0.5%. These summaries describe read-level quality and interval coverage; they do not establish sensitivity for every locus or variant class.

The pooled `reference.cnn` is already available locally in 35 byte-identical copies (SHA-256 `848c052274264ac544897977648860455ac742d89bc1de5edc7ece1ca185eeda`); it should not be requested again. The five per-normal coverage CNN files are absent from the transferred bundle. Their handoff-reported depth range of 75.9–95.7× has not been independently validated and should not be presented as a checked result. Retrieve the small source files if possible; otherwise state the provenance limit. The provenance review found 42,705 of 42,709 antitarget bins at zero depth in every one of the 23 models. Only four antitarget bins had positive depth (two autosomal and two on X). In each model, segment probe totals matched positive-depth bins, consistent with `--drop-low-coverage`; the zero-depth antitargets were not segmented as biological deletions. Thus the retained profiles are effectively supported by captured targets, and their segment spans interpolate between targets. Retain relative CN/arm summaries and the segment-span FGA with that explicit qualification; do not claim direct intergenic copy-number validation. Exact manual logs/input provenance and reference-normal QC remain incomplete.

## Necessary work versus optional analyses

**Necessary before claiming the manuscript is complete:** finish the exact run/input/filter reconciliation; make the new QC tables and figures agree with their primary records and denominators; retain the verified target-supported CNV interpretation and resolve any additional material input discrepancy identified by the final provenance audit; obtain or explicitly account for the laboratory/ethics/stock records above; verify prior use and declarations; and replace every manuscript TODO with a defensible statement. Recheck figure legends, references and release counts after those changes. The 23-model cohort and 6,194-candidate record do not need expansion to accomplish this.

**Useful if existing inputs make them straightforward:** targeted coverage/callability summaries for genes displayed in the oncoprint, independent sample-identity comparison across modalities, and missing reference-normal QC exports. These would strengthen reuse and any claim of absence. They are not a reason to label a model wild type when only the absence of a retained candidate is established. If identity checks uncover a discrepancy, resolving that discrepancy becomes necessary.

**Optional for this descriptor:** new allele-specific copy-number fitting, HRD scores, MSI testing, additional signature fitting, discovery of RNA subtypes, differential-expression/pathway studies, proteomic re-search for new biological questions, and drug/functional experiments. Strong claims about those endpoints would require suitable validation; the current descriptor can instead omit those claims. Existing limitations cannot be solved by adding more unsupervised plots.

Maintain the prior audit's interpretation limits: tumour-only tiers prioritise candidates rather than prove somatic origin; high VAF does not establish LOH; relative CN does not measure absolute/allele-specific CN or HRD; histotype and centre remain partly associated; related models do not add independent patients; some assays represent different passages; sparse histotypes cannot support broad rankings; protein missingness is plex- and abundance-dependent; bridge precision is technical repeatability, not absolute accuracy; and the RNA/protein spread comparison cannot identify ratio compression as its sole cause. TOV21G signature fitting is sensitive to reference/opportunity correction and cannot establish MSI. Keep the historical TOV112D/COV434 annotation qualifications.

## Concise author/laboratory request

The following can be sent as one coordinated request, with source records attached where available:

1. Confirm culture/harvest conditions and assay-stock relationships, including the listed derivation exceptions and VOA sources.
2. Supply the RNA extraction and integrity/QC records; identify any unrecoverable historical details.
3. Supply the actual proteomics acquisition/search/processing protocol, channel map, pooled-standard formula and CV definition.
4. Confirm WES DNA extraction/library construction and the provider's TOV81D pool/alias and passage composition. The recovered instrument/read-length/capture facts do not need to be requested again.
5. Supply stock-linked STR and mycoplasma records from each contributing laboratory, or confirm what cannot be documented.
6. Confirm the applicable ethics/consent and genomic-sharing determination and complete the Human Data Checklist.
7. Identify every publication/deposit using these exact datasets, then agree authorship, contributions, interests, funding and the AI-assistance statement.

Remaining manual CNV provenance, existing normal QC and manifest reconciliation are analysis/retrieval tasks for the computational team, not vague questions for the authors to reconstruct from memory.

## Bibliographic and submission cleanup

Use the exact historical software versions in Methods, even though newer releases exist. The official Sarek 3.5.1 page recommends the 2024 paper and retains the original 2020 paper. A pipeline reference does not replace documentation of the executed options. [Sarek 3.5.1 citation guidance](https://nf-co.re/sarek/3.5.1/)

| Action | Primary reference |
|---|---|
| Add the current Sarek methods paper; retain Garcia et al. 2020 as appropriate | Hanssen, F. et al. Scalable and efficient DNA sequencing analysis on different compute infrastructures aiding variant discovery. *NAR Genomics and Bioinformatics* **6**, lqae031 (2024). [DOI](https://doi.org/10.1093/nargab/lqae031) |
| Add for BWA-MEM 0.7.18; do not describe this as BWA-MEM2 | Li, H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv* 1303.3997 (2013). [Primary preprint](https://arxiv.org/abs/1303.3997) |
| Add for fastp 0.23.4 | Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* **34**, i884–i890 (2018). [DOI](https://doi.org/10.1093/bioinformatics/bty560) |
| Add for mosdepth 0.3.8 and the new coverage panels | Pedersen, B. S. & Quinlan, A. R. Mosdepth: quick coverage calculation for genomes and exomes. *Bioinformatics* **34**, 867–868 (2018). [Primary article](https://academic.oup.com/bioinformatics/article/34/5/867/4583630) |
| Identify vcf2maf 1.6.22 explicitly; cite its software record or URL without inventing a paper/DOI | [Official v1.6.22 release](https://github.com/mskcc/vcf2maf/releases/tag/v1.6.22) |
| Add Nextflow if named in the expanded Methods | Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. *Nature Biotechnology* **35**, 316–319 (2017). [Primary article](https://www.nature.com/articles/nbt.3820) |
| Add SAMtools if the BAM-to-read conversion/alignment-QC operations are described | Danecek, P. et al. Twelve years of SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021). [Primary article](https://academic.oup.com/gigascience/article/10/2/giab008/6137722) |
| Retain existing method references | Benjamin et al. 2019 Mutect2; Talevich et al. 2016 CNVkit; McLaren et al. 2016 VEP. Match the cited method to the implemented step. |

If BWA-MEM2 is ultimately specified for the five reference-normal alignments, document its actual version separately and add its own primary citation; the tumour BWA-MEM citation cannot establish which normal-alignment software ran. Add the relevant derivation/proteomics primary papers only after confirming their applicability to the present data. Keep the DepMap version-1 dataset citation already present.

Preserve the author's current author–year/DOI reference-manager workflow during substantive revision. Before submission, convert citations to the journal's sequential numbered style, remove placeholders and tracked changes, and check all newly added references. Place Code Availability immediately before References; handle contributions/interests/acknowledgements according to the submission system, and retain the required funding statement. Acknowledgements can be omitted if inapplicable. Prepare a readable manuscript PDF and separate supplementary information. [Scientific Data submission guidelines](https://www.nature.com/sdata/submission-guidelines)

The current author-confirmation companion still carries a version-6 heading despite accompanying v7. The v8 companion should have its heading and resolved WES requests updated. This is a routine documentation correction, not an additional scientific blocker.

## Read-only check of the developing v8 text

The developing v8 main manuscript was read after the QC and CNV findings were incorporated. The new numeric summaries and the target-supported CNV/segment-span qualification agree with the reviewers' findings. The remaining editorial suggestions sent to the parent were to restore the author-requested author–year plus DOI/URL callouts consistently for new WES references, move detailed interval split counts and optional historical archive gaps into the provenance records/checklist, and qualify “Descriptive plots retain all models” because some displayed histotype means use patient representatives. The first-release claim remains conditional on A02. No v7 manuscript content was changed by this review.
