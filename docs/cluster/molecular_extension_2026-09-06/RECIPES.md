# Corrected pilot recipes and command review

Use with [EXECUTION_PLAN](EXECUTION_PLAN.md). These are reviewed command specifications for the cluster agent to integrate into its scheduler wrappers, not claims of an executed or cluster-tested workflow. The incoming `proposed_commands/*.sbatch` are immutable historical proposals and must not be submitted unchanged. Preserve their source hashes; write revised wrappers in a new working directory and return their exact text.

## Defects to fix in the incoming wrappers

| Script | Required revision |
| --- | --- |
| All five | Read explicit manifests, validate exactly one model per array index, recheck paths/indexes and fail on missing inputs. Record actual commands, versions, failures, exit status and scheduler resources. Use fresh result directories; exclude checksum files from their own manifests. Stage substantial computation through Slurm. |
| `mex02_03_locus_review.sbatch` | Use the new 30-row variant list. Preserve complete alignments locally for sequence/quality/haplotype review; the nine-column SAM excerpt is insufficient. Keep candidate IDs when merging overlapping windows. Explicitly resolve the reference for every CRAM reader, not only `samtools view`. Distinguish permissive/quality-filtered depth and intervals/exons/CNR bins. Avoid unlabelled comparisons of duplicate-excluded models with unmarked normals. |
| `mex04_msisensor2_pilot.sbatch` | Keep complete md CRAM, staged hg38 models and `-c 20`. Check container/model hashes, BAM/index and flag/filter behavior. Return site accounting and depth sensitivity; a score alone is not enough. Fix manifest self-inclusion. |
| `mex05_mutect2_germline_sites.sbatch` | Use complete md CRAM from `execution_models.tsv`; remove wildcard-log path extraction. Precede the five-model run with the fixed-SNP md/recal comparison. Keep genotype flags true and padding 75. Produce full filtered VCF plus companions and retain non-PASS germline/PoN records. |
| `mex05b_normals_mutect2_for_mappingbias.sbatch` | Conditional sensitivity: mark duplicates in derived normal copies first, record aligner/BQSR differences, use explicit normal paths. Validate biallelic SNP normalization and sample AD/missingness in the multisample input. Call it a small unmatched, processing-mismatched panel. Substantial panel preparation belongs in a scheduled job. |
| `mex06_purecn_pilot.sbatch` | Use `--tumor CNR`; remove `--log-ratio-file CNR`, Mutect2 `--stats-file`, and calling-BED `--intervals`. Feed complete filtered VCF. Add `--model-homozygous`, use `--model beta`, align SNP padding to 75, and set minimum base quality 20. Hold ploidy/max-copy settings fixed when changing purity. Add BAF-aware segmentation before scar interpretation. Discover the wrapper via R, not hard-coded module internals. |

The original shared `IntervalFile.R` step used 242,421 calling intervals against 204,706 CNR rows and could race between array tasks. It is unnecessary for the CNVkit import and is removed. If interval annotation is later added, derive it from the exact imported-bin coordinates and validate a one-to-one relationship before fitting. Existing CNR `gene` fields include coordinate strings; annotate accepted segments/bins with a versioned gene/exon map before producing gene-level conclusions, rather than treating those labels as gene symbols.

## Alignment and focused-depth checks

Resolve `MODEL_CRAM`, `REF` and output paths from the returned manifests. `quickcheck` and index/header checks do not establish full read integrity, lineage or interval completeness. Compare existing task evidence and read-group/library identifiers as specified in the plan. For tools without an explicit reference switch, convert the bounded region to BAM using `samtools view -T "$REF"` and analyze that BAM; never let CRAM decoding depend on an unrecorded remote reference lookup.

For coverage, state what is counted. `bedcov` reports a sum of base depths (and optionally read counts); divide by interval length only when reporting the corresponding mean. Exclude CIGAR D/N operations for aligned-base coverage. Keep zero rows and count overlapping mates consistently. The permissive and quality-filtered tracks are diagnostic views, not new canonical CN ratios. The [samtools depth](https://www.htslib.org/doc/1.22/samtools-depth.html) and [bedcov](https://www.htslib.org/doc/1.22/samtools-bedcov.html) manuals define differing flag/quality/overlap semantics; use the installed 1.22.1 help and record actual settings.

For example, on a locally converted region BAM:

```bash
# Diagnostic high-quality track: MAPQ >=20, BQ >=20, remove overlapping mates,
# default excluded flags plus supplementary reads; explicit zero-depth output.
samtools depth -aa -s -Q 20 -q 20 -G 2048 -r "$REGION" "$REGION_BAM" > "$DEPTH_TSV"
```

Apply this as a declared sensitivity alongside the permissive track. It does not reproduce CNVkit coverage exactly and cannot by itself phase a complex indel. Keep full read bases/qualities/CIGAR and inspect fragment-level haplotypes separately. Validate the resolved version's cumulative excluded-flag behavior before use.

## MSI

The [official MSIsensor2 interface](https://github.com/niu-lab/msisensor2) uses BAM plus an adjacent index and reference-specific models. With verified paths:

```bash
msisensor2 msi -M "$MSI_MODEL_DIR" -t "$MODEL_BAM" \
  -o "$RESULT_PREFIX" -c 20 -b 4
```

Keep the score, `_dis`, `_somatic` and log files. Reconcile reported total/unstable sites with per-site outputs; record how many modeled loci overlap the capture and how many have useful data. Distinguish software-test controls from biological MSI controls. No matched positive/negative assay controls were found in the WES roots. Threshold/depth sensitivity and site resampling characterize robustness without making the cohort clinically calibrated.

## Mutect2 and PureCN input preparation

Retain GATK 4.5.0.0 and the hash-identified FASTA, PoN and AF-only gnomAD resource. MEX05 changes the complete alignment choice, both `--genotype-germline-sites true --genotype-pon-sites true`, and `--interval-padding 75`; it is a distinct SNP-input analysis. Keep orientation, contamination and FilterMutectCalls companions. The PureCN input is the **full filtered VCF**, not a PASS extraction. Do not convert a missing locus into a homozygous-reference call.

Verify autosomal subsets of target-only CNR, corrected CNS and VCF for the primary pilot. Filter by exact contig names `chr1`–`chr22`, preserving coordinates and values; save input/output hashes and excluded counts. A subset is a new derived input with its own hash, not the original 204,706-bin file. The VCF must retain its complete records/FILTER values within that scope. Export the matching CNS with installed CNVkit 0.9.10 `export seg`; inspect the actual output's chromosome encoding, 1-based starts, end coordinates, probe counts and ratio values before import. Never combine raw CNR with script 08's independently centred segments.

The wrapper below is verified against [PureCN 2.14.1 source](https://bioconductor.org/packages/3.21/bioc/src/contrib/PureCN_2.14.1.tar.gz) (archive SHA-256 `e6f08ff5a5a586ccc96f038d4919f921be18ebf1709222c6bf0de6d9c8d353d6`). `readCoverageFile()` handles CNR; `readLogRatioFile()` expects a different layout. `filterVcfMuTect()`'s stats-file argument is for MuTect1. `filterVcfMuTect2()` recognizes FILTER annotations without selecting only PASS. The package's `Quick.Rmd` specifies the GATK base-quality setting; `--model-homozygous` retains high-purity LOH candidates. Capture installed wrapper/source identity and help because an environment label alone is insufficient.

```bash
PCN=$(Rscript -e 'cat(system.file("extdata", package="PureCN"))')
test -s "$PCN/PureCN.R"
# CNR/SEG/FULL_FILTERED_VCF are the matching, validated autosomal inputs.
# FIT_PREFIX is a fresh path; do not overwrite an earlier fit.
Rscript "$PCN/PureCN.R" \
  --out "$FIT_PREFIX" --sampleid "$SAMPLE_ID" \
  --tumor "$CNR" --seg-file "$SEG" \
  --vcf "$FULL_FILTERED_VCF" --genome hg38 --sex F \
  --fun-segmentation none \
  --model beta --model-homozygous \
  --min-base-quality 20 --interval-padding 75 \
  --min-purity 0.30 --max-purity 0.99 \
  --min-ploidy 1.5 --max-ploidy 8 --max-copy-number 16 \
  --post-optimize --out-vcf --seed 20260906 --cores 4
```

Repeat with minimum purity 0.80 and all other parameters fixed. A 12-versus-16 maximum-copy sensitivity is a separate comparison. That setting controls allele fitting, not a hard upper bound on reported total copy number. Preserve alternative optima and ceiling/boundary warnings.

For a BAF-aware comparison, verify PSCBS is installed and use `--fun-segmentation PSCBS` on the same inputs. Keep the SEG argument for the external-input pathway; this method recomputes boundaries from bins/SNPs rather than treating the input CNS boundaries as fixed. `Hclust` only clusters existing segments and cannot alone establish new copy-neutral LOH boundaries. Record internal centring, native CNR use and any segment-derived synthetic-ratio fallback. Failed or unstable BAF segmentation does not yield a valid zero scar score.

## Optional mapping-bias panel

Use duplicate-marked derived normal alignments and explicit allele QC. PureCN's mapping-bias VCF reader uses AD/ALT and does not apply FILTER itself. For an AD-preserving multisample VCF: normalize against the same reference, retain biallelic SNPs, ensure distinct sample names, merge without combining alternate alleles (`bcftools merge -m none`) and never use `--missing-to-ref`. Check source AD values survive exactly, including missing observations. An appropriately populated GenomicsDB is another supported input. A sites-only final variant PoN is not a mapping-bias input.

Record normal support/missingness and allele-balance behavior per SNP before `NormalDB.R`. Five normals cannot support the default minimum panel size of ten normals for position-specific dispersion fitting; using `betabin` without such support can substitute imputed dispersion. The primary `beta` diagnostic and any panel/betabin sensitivities must be labelled separately. See PureCN's `calculateMappingBiasVcf()` and `.imputeBetaBin()` in the pinned source above. No panel makes the five public normals donor-matched or removes their alignment-processing differences.

## Scar computation and result handling

Use a pinned [scarHRD](https://github.com/sztup/scarHRD) version only after the plan's fit/BAF gates. Preserve total/minor CN, assembly, centromeres, exclusions, ploidy and competing fits. Return LOH/TAI/LST components and sensitivity without a clinical cutoff. Complete exome alignments do not add the structural-variant measurements needed by standard WGS HRDetect/CHORD.

For result checksums, write the temporary manifest **outside** the tree being hashed, exclude existing `SHA256SUMS`, then move the completed file into place. Preserve raw and failure logs on managed storage; return reviewed small evidence and diagnostic PDFs. Never include read slices, whole VCFs or alignments in the public task overlay.
