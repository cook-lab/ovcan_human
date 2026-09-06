# WES coverage return and CNV correction — 6 September 2026

The second archive closes the request for coverage tables from all five CNV-reference exomes. Independent checks also identified a technical antitarget artifact in the archived CNV segments, prompting target-only resegmentation of all 23 models. Manuscript v9 incorporates the coverage evidence and correction. Original sequencing, coverage and ratio inputs remain unchanged.

## What arrived and what was verified

`ovcan_human_wes_cnv_coverage_2026-09-06.tar.gz` contains 109 files: 46 model coverage tables, ten normal coverage tables, the pooled reference, derived target BED, four command/wrapper files, 46 historical plots and a checksum manifest. There is **no handoff document or README**. Its SHA-256 is `c5d79958ae05115c975d63e3858fae7cd372b37c593295e0e2472bbf38c6706a`.

All 108 payload files match the manifest. The manifest's own self-entry is incorrect; this packaging defect is recorded separately and does not invalidate the payload checks. All returned model coverage files and the pooled reference match the retained historical copies. The command files are unchanged from the first handoff and do not supply new execution logs. See [archive validation](archive_validation.json), [file inventory](archive_inventory.csv) and [provenance assessment](provenance_assessment.md).

## Reference coverage now supported by primary evidence

| Reference exome | Mean target depth |
| --- | ---: |
| SRR4039087 | 95.47× |
| SRR4039088 | 83.75× |
| SRR4039089 | 75.85× |
| SRR4039096 | 92.40× |
| SRR4039097 | 95.72× |

These length-weighted means use 290,475 CNVkit target bins covering 63,514,049 bp; their median is **92.40×**. They verify the earlier rounded recovery summary. Before reference masking, 10.9–16.5% of each normal's target bins had zero mean depth. Zeros were distributed across chromosomes; the files do not establish their technical cause.

The existing reference mask excludes all 23,406 targets with zero mean depth in all five normals. Its **204,706 retained target bins cover 44,729,346 bp**: 195,873 have positive depth in all five normals, 5,478 in four and 3,355 in three. Thus each retained target is supported by at least three normals. Retaining a target in the reference mask does not imply sufficient depth in every model; model-specific low-coverage filtering follows during segmentation.

Normal mean antitarget depths are 0.452–0.633×. In contrast, every model has zero depth in 42,705 of 42,709 antitarget bins, and its four positive bins overlap capture targets. This does not supply independent off-target validation. An interval-restricted model alignment is a plausible explanation to check against the missing manual-workflow input records, not a verified cause.

The CNN target means and the model mosdepth metrics in Figure 2 measure different alignment stages/interval accounting. They must not be combined as directly comparable estimates. A mean depth per bin cannot establish the fraction of individual bases with ≥10×, ≥20× or ≥30× coverage. [Coverage definitions and validation](coverage_validation.md) · [Sample summary](../../output/wes_cnv_coverage_sample_summary.csv) · [Cross-normal support](../../output/wes_cnv_coverage_reference_support.csv).

## CNV artifact and corrective analysis

The four overlapping antitargets carry normalized log2 ratios of +12.48 to +19.75 with appreciable weights. A reconstruction from exact contributing-bin membership reproduces all 5,670 archived segment means to their printed precision. It identifies two 1.89-Mb chromosome 1 spans, in OV2295 and TOV2835EP, whose means were inflated by approximately 0.31 log2. Removing their antitarget contributors at the original centering baseline changes their gain classification. This diagnostic justified resegmentation; manually replacing two means would leave possible breakpoint effects unresolved.

Script 29 removes all antitarget rows from the original CNR files and runs native CNVkit 0.9.10 CBS on the unchanged target ratios and weights. It preserves the original reference mask. Parameters are alpha 0.0001, low-coverage filtering, the default rolling outlier threshold of 10, no additional CBS smoothing, and native chromosome-arm partitioning. DNAcopy uses the CNVkit seed 0xA5EED. New segments are then centred and summarised by script 08 using the established autosomal, patient and arm rules. No variant-calling, alignment or allele-specific analysis is part of this correction.

Original-CNR control replays reproduce OV2295 and TOV2835EP (202 segments each) and OV1369-R2 (248 segments) exactly in coordinates, probe counts and log2 ratios. The third control was added after the FGA sensitivity finding, excluding a mismatch in reproduction settings as its explanation. The corrected 23-model run retains all 4,692,292 positive target-bin observations (203,913–204,064 per model), producing 6,073 segments across all contigs, of which 5,808 are on chromosomes 1–22 and X. [Run manifest](../../output/wes_cnv_target_only/manifest.csv) · [Baseline replay](../../output/wes_cnv_target_only/baseline_replay.csv) · [Runtime verification](runtime_reproducibility.md).

The main numerical changes, measured against public pre-correction revision `854def3`, are:

| Measure | Archived profiles | Target-only profiles |
| --- | ---: | ---: |
| HGSC patient median FGA | 0.622317 | 0.624126 |
| 20q gain, HGSC patients | 10/11 | 10/11 |
| 3q gain, HGSC patients | 9/11 | 9/11 |
| 17p loss, HGSC patients | 9/11 | 9/11 |
| 19q gain, HGSC patients | 6/11 | 5/11 |
| OV1369-R2 FGA | 0.641469 | 0.859173 |
| TOV81D FGA | 0.021 | 0.023 |

There are 23 changed model/arm calls and 15 arm-frequency rows with changed model or patient counts. Across the declared parameter sweep, 13 patient-frequency rows change. Two model/point-locus classifications change; these exploratory point checks are not focal amplification or deletion validation. No authentication/consistency classification changes. The exact comparisons are in [correction_impact.json](correction_impact.json) and its linked-name CSV companions, including [arm calls](correction_arm_call_changes.csv), [arm frequencies](correction_arm_frequency_changes.csv) and [FGA by model](correction_fga_comparison.csv).

**OV1369-R2 demonstrates substantial baseline sensitivity.** Its autosomal probe-weighted median shifts from −0.18947 to −0.11787. Holding the corrected segments fixed, the old center gives FGA 0.67155 and the new center gives 0.85917. Bin-based alternatives also give high FGA (0.9114 with the unweighted target-bin median; 0.9383 with the CNR-weighted target-bin median), confirming baseline ambiguity rather than identifying a uniquely neutral state. These threshold effects do not represent newly acquired biological alterations or an absolute copy-number result. Removing antitarget rows also changes native gap-based partitioning and the segment-span denominator, so the full correction has effects beyond deleting the two conspicuous chromosome 1 means. The manuscript retains the declared centring convention and adds this example explicitly. [Independent sensitivity evidence](ov1369_centering_diagnostic.json) · [Reproduction script](../../scripts/29c_cnv_centering_sensitivity.py) · [Sensitivity-curve caption and data](ov1369_centering_sensitivity_caption.md).

Figure 4D, Supplementary Figure S5, model metadata, the local browser payload and the existing 49-file processed release were regenerated from the corrected summaries. The S5 FGA axis now accommodates the largest value without clipping. Both changed PDFs use embedded Arial and were rendered and inspected; legends remain external. The 6,194 variant candidates, cohort counts, RNA/protein matrices and model mosdepth QC are unchanged.

## Remaining work and current files

- [Manuscript v9](../../docs/manuscript/v9/OvCAN_Scientific_Data_draft_v9.md) and [author confirmations](../../docs/manuscript/v9/author_confirmation.md).
- [Responses to 13 comments in the author's v7 Word file](../../docs/manuscript/v9/user_comment_responses.md), carried into v9. They resolve the new-data and AI-wording questions, add highlighted collaborator-review protocols and explain the PC1 decomposition. The Cellosaurus provenance item is now an explicitly documented analysis-team result.
- [Current targeted cluster follow-up](../../docs/cluster/recovery/2026-09-06/FOLLOWUP.md). Coverage CNNs, pooled reference and target designs are received; do not request them again.
- Seek existing successful manual CNVkit execution records, exact input alignment identities and normal preprocessing/QC records; retain the narrower converter, workflow, reference and liftover provenance requests where recoverable.
- DNA extraction/preparation, the TOV81D pool alias, culture/harvest, RNA extraction/integrity, proteomics acquisition/scaling, stock tests and author/institutional declarations still require records or decisions. This archive does not resolve those separate paper requirements.

The raw archive, provider/run records and derived CNR scratch files remain outside Git. Curated evidence and analysis outputs are reviewable in the public repository. Organizing deposition, assigning accessions and redeploying the hosted browser remain separate work.

## Reproduction

With the two exact handoffs and original CNR inputs restored, scripts 26–28 regenerate the independent coverage, provenance and archive checks. Script 29 additionally needs the pinned CNVkit environment and R/DNAcopy. Script 29b compares regenerated summaries with public revision `854def3653918034055ac0e8e94907c852da2484` using local Git history. See [reproduction instructions](../../docs/REPRODUCIBILITY.md). These scripts never execute the recovered Slurm wrapper or submit jobs.
