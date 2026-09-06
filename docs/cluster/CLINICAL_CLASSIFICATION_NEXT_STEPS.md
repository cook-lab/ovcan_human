# HRD and molecular classification: focused next steps

**Updated operational handoff:** the author's broader follow-up is now specified in [MEX01–MEX11](molecular_extension_2026-09-06/AGENT_TASK.md), with exact model/locus tables, recipes and return templates. Use that file for the next cluster-agent task. This document retains the initial HRD/CCNE1 rationale and laboratory follow-up context.

This is the follow-up to the author's 6 September 2026 request to explore HRD, CCNE1 and other clinically relevant model annotations. Read [the findings](../../reports/clinical_classification_2026-09-06/README.md) and [HRD feasibility](../../reports/clinical_classification_2026-09-06/HRD_FEASIBILITY.md). The previous [WES provenance follow-up](recovery/2026-09-06/FOLLOWUP.md) remains useful; do not repeat retrieval of the already verified coverage CNNs, reference or target BED.

## First return: inspect existing resources

Keep original runs read-only. Search the operator-supplied project/run locations, not the whole filesystem. Return paths, sample aliases, sizes, versions, checksums of small files and a status table. Do not transfer large BAMs or submit new analysis jobs during initial discovery.

1. **Final alignments for all 23 exact models**: identify BAM/CRAM plus index, the complete reference FASTA/dictionary identity, input filtering/interval restriction and the corresponding run/sample alias. The prior coverage return suggests almost no independent model off-target evidence; establish whether alignments were interval-restricted before proposing an allele-count workflow.
2. **Existing allele-specific outputs**: inventory PureCN/ASCAT/Sequenza/FACETS or other total/minor CN fits, allelic pileups, heterozygous SNP counts, ploidy/purity alternatives, fit plots and logs. Confirm model identity, genome build and assay. Do not call historical total-CN-only `.call.cns` allele-specific data.
3. **SNP input eligibility**: find any unfiltered somatic/germline-aware Mutect2 VCF generated with germline and PoN sites retained, common-SNP genotypes/AD counts, or suitable pileups. The retained VCFs explicitly used `--genotype-germline-sites false`, `--genotype-pon-sites false`, and `--interval-padding 0`; PASS-only files retain very few common informative SNPs. Their existing variants do not establish an unbiased allelic denominator.
4. **Appropriate normal and reference resources**: determine whether matched donor normals exist at all (none are documented in current calling); otherwise identify capture-matched normal alignments/interval coverage sufficient for a tumour-only workflow, mapping-bias estimates, population allele frequencies and exact target definitions. The five public CNV-reference normals and the Mutect2 variant PoN are distinct resources; neither is a matched donor normal.
5. **Existing MSI or copy-number orthogonal results**: collect completed MSI locus analyses/MMR tests, CCNE1 or ERBB2 FISH/ddPCR/array measurements, karyotypes/DNA-content measurements and original assay definitions if available. Report missing evidence without interpreting absence as a negative result.

## Proposed analysis after input review

The initial response should propose commands and resource estimates; execution can be authorised by the cluster operator after review. A sensible pilot is a tumour-only allele-specific workflow such as PureCN with suitably prepared SNP/coverage inputs. Do not pair a tumour with an unrelated reference normal and describe it as a matched-normal Sequenza/FACETS analysis.

Start with a few contrasting models: TOV81D (near-flat relative profile and a curated pathogenic BRCA2 allele), OV2085 and TOV2929D (strong CCNE1 gains), and OV1369-R2 (known centring sensitivity). These are diagnostic cases, not known HRD-positive/negative controls. Include established positive/negative controls if compatible data are available. Retain alternative ploidy solutions, high-purity boundary behaviour, SNP/BAF plots, coverage residuals, minor-CN uncertainty and model-fit QC; do not force a diploid or 100%-purity solution merely because these are cell lines.

For credible allele-specific fits, return per-segment total/major/minor CN, fit estimates and alternatives, gene-level CCNE1/ERBB2 absolute CN and ploidy-normalised ratios, BRCA/HRR allele status and uncertainty. Only then evaluate LOH/TAI/LST scar components on the appropriate assembly with a pinned implementation. Report continuous scar scores first; clinical thresholds are assay- and validation-context-specific. Existing FGA or SBS3 refits must not be renamed HRD scores. Standard WGS-based HRDetect/CHORD workflows are not drop-in WES analyses.

Prioritise exact read-level review of TOV81D BRCA2 chr13:32371000 AAG>A (VCF representation; MANE c.8537_8538del), CDK12 candidates in the 3133 family and OV3291, and the small CCNE1 gene-bin peak in TOV2835EP. Return mapping/strand/read-position evidence, allele counts, local CN and transcript annotation. Germline versus somatic origin, pathogenicity, biallelic loss and HR function are separate columns.

## Laboratory records / experiments, separate from cluster computation

- CCNE1: absolute-copy FISH/ddPCR or an equivalent validated assay, plus Cyclin E1 protein assessment in the same stocks. TOV2929D and OV2085 are the strongest DNA leads; OV3291, TOV3133D and TOV2881EP are additional relative-gain candidates. TOV3291G is RNA-high but is a distinct isolate from OV3291.
- HR function: existing or new controlled RAD51-foci/HR repair assays, with cell-cycle and DNA-damage controls; genomic scars do not directly measure present repair function or PARP-inhibitor sensitivity.
- TOV81D: historical literature reports retention/expression of a wild-type BRCA2 allele. Its pathogenic variant does not itself establish current BRCA2-null status.
- MMR/MSI: TOV21G is a priority for existing MSI/MMR validation. No retained MMR coding variant does not exclude methylation or other MMR defects.
- HER2/FRalpha: validate membrane protein with appropriate assays in prioritised stocks; RNA/TMT measurements are not clinical IHC positivity.
- Proteomics: obtain the actual isoDoping peptide list, channel corrections, rollup and normalisation definitions. The current CCNE1 biological-channel imports are verified; the flag alone does not establish quantitative failure.

Return a compact handoff with received/unavailable/not-searched distinctions, commands and software identities, QC tables and plots, model/assay/passages, provenance and unresolved ambiguity. Keep raw sequence/alignment data and unreviewed provider logs out of the public repository. This request does not organise deposition or authorise new clinical labels before validation.
