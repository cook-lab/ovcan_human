# Bounded pilot recipes and decision gates

These templates are not executed by preparing the handoff. Resolve inputs and installed versions first, then obtain the cluster operator's execution authorization. All outputs belong in a new analysis directory, with source runs read-only. Save the actual fully resolved command and environment for each accepted result.

## Alignment inventory

Read existing logs and manifests before touching large files. For a located BAM, inexpensive headers/index statistics help establish identity and completeness:

```bash
samtools quickcheck -v "$MODEL_BAM"
samtools view -H "$MODEL_BAM" > "$RESULT_DIR/alignment_header.sam"
samtools idxstats "$MODEL_BAM" > "$RESULT_DIR/alignment_idxstats.tsv"
```

Define `MODEL_BAM` and `RESULT_DIR` explicitly from verified paths. Quickcheck detects some file/truncation problems but does not validate every record; headers/idxstats do not establish the original interval filtering. For CRAM, establish the exact reference and do not convert or regenerate large files merely to perform discovery. Sanitise paths/read-group fields before placing header excerpts in the public repo.

## MSI pilot

The official [MSIsensor2 documentation](https://github.com/niu-lab/msisensor2) provides a tumour-only module using a BAM, adjacent index and reference-specific pretrained models. Its model choice must match actual reference/contig conventions. Record executable and model-directory versions/hashes and the exact number of callable microsatellite sites.

```bash
# Illustrative GRCh38-compatible, WES-depth pilot after input verification.
msisensor2 msi -M "$MSI_MODEL_DIR" -t "$MODEL_BAM" \
  -o "$RESULT_DIR/msisensor2" -c 20 -b 4
```

Retain the main score file, `_dis` and `_somatic` files. Inspect distribution and site depth before scaling out. Report any source threshold with provenance separately from the continuous score; model performance on this capture/passage context is not demonstrated by obtaining a number. Include suitable controls if available; otherwise say that calibration is absent. Prioritize TOV21G, then the contrasts listed in AGENT_TASK. RNA/TMT expression cannot substitute for this locus analysis or for MMR IHC.

## Allele-specific CN and genomic scars

Follow installed-version documentation, with [PureCN best practices](https://bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/Quick.html) and [CNVkit import guidance](https://bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/Quick.html#third-party-segmentations) as primary references. Core input requirements are retained germline-SNP alleles and compatible coverage/mapping-bias information. Mutect2 input preparation should retain germline/PoN sites (`--genotype-germline-sites true --genotype-pon-sites true`), including eligible non-PASS germline records. Proposed 50–75 bp caller padding is useful only if those flanking reads exist in the BAM. Establish this before rerunning anything.

For CNVkit import, use native corrected target-only CNR and matching CNS. Export a SEG file with the installed CNVkit's `export seg` command and confirm chromosome encoding agrees with the VCF and PureCN. The pilot should save the exact command assembled from the installed help, including input/SEG/VCF/mapping-bias files, reference, seed and filtering parameters. New coverage generation is unnecessary if existing compatible third-party inputs pass QC. The previous four-antitarget problem must not be reintroduced to generate nominally broader coverage.

Compare plausible high-purity priors and alternative ploidy solutions; document evidence for choices. The cell-line [DepMap PureCN workflow](https://depmap-omics.readthedocs.io/en/latest/workflows/PureCN.html) is a useful precedent for high-purity sensitivity, not a validated fit for these samples. Evaluate high-copy settings for AKT2/CCNE1 peaks rather than allowing a default to dictate biology. A flat log-ratio profile does not determine absolute ploidy. Return fitting alternatives, BAF distributions and uncertainty as well as preferred outputs.

Only acceptable allele-specific fits proceed to a pinned [scarHRD implementation](https://github.com/sztup/scarHRD) or equivalently documented LOH/TAI/LST definitions. Validate build, chromosome names, total/minor CN, centromeres, exclusions and ploidy. Report component sensitivity to alternative fits, sparse exome coverage and segment gaps. Do not manufacture minor CN, fill missing chromosomes with normal values, or apply a clinical positive/negative threshold without validation.

## Local sequence review

Use the model-specific variant BED and `locus_review_requests.tsv` to restrict inspection. Variant BED uses zero-based half-open coordinates; source VCF positions are one-based with an anchor base for indels. Preserve both representations. Inspect flanking sequence, read alignment and combined local haplotypes before normalizing consequences. Zero reference reads do not resolve allele-specific copy number in an aneuploid cell line. Zero-depth bins must stay visible, including bins excluded from segmentation.

For CDKN2A, distinguish the p16 and p14ARF transcripts. The supplied canonical-exon check is a starting point and includes UTRs; it is not a full isoform/CDS audit. For NF1, gene-span targets can overlap intronic/nested features, so exon-specific mapping and residual coverage are essential. For apparent gains, return gene/amplicon versus arm-level evidence, absolute total/minor CN only when fitted credibly, and passage-matched orthogonal evidence if it exists.

Do not generate TMB from retained candidate counts and a kit-size denominator. Recover callability and credible somatic filtering first. Fusion/splice analysis likewise requires RNA reads or existing alignments/junctions; the current gene-level RNA matrix contains neither.
