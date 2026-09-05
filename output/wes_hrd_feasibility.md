# Genomic HRD — feasibility assessment (WES)

**Date:** 2026-07-23  |  **Script:** `scripts/09_wes_hrd.R`  |  **Project:** OvCAN human ovarian cancer cell-line resource

## Verdict: **NOT FEASIBLE FROM ARCHIVED DATA**

Genuine genomic HRD (the genomic-scar score = HRD-LOH + telomeric allelic imbalance + large-scale state transitions) **cannot be computed from the archived WES data**, because the copy-number calls are **total copy number only — there is no allele-specific / B-allele-frequency information**. We therefore report **no HRD score** rather than a fabricated one.

## What genomic HRD requires
- The composite HRD-scar score combines LOH, TAI, and LST; LOH and TAI require allele-specific information, so an LST-like count from total-copy profiles cannot substitute for the composite score.
- `scarHRD` (Sztupinszki 2018, npj Breast Cancer; HRD-sum r=0.87 vs SNP array, robust at 30x) consumes **per-segment allele-specific copy number** (major `A_cn` + minor `B_cn`), normally produced by **Sequenza** or **ASCAT** from a BAF-bearing SNP VCF.
- **Total copy number / log2 (plain CNVkit) is insufficient:** LOH is invisible in total CN. HRDetect (Davies 2017) and CHORD (Nguyen 2020) are **WGS-only** (SV-signature based) and do not apply to exomes.

## Evidence from the archived files (this script)
- Allele-specific columns (`baf`/`cn1`/`cn2`/`A_cn`/`B_cn`) in CNVkit `.cnr`/`.cns`/`.call.cns`: **NONE**.
- `.call.cns` columns are total-CN only: `chromosome`, `start`, `end`, `gene`, `log2`, `cn`, `depth`, `p_ttest`, `probes`, `weight` (has `cn`/`log2`, no allele split).
- CNVkit `batch` calls passing a SNP `--vcf` (the BAF input that would enable allele-specific segmentation): **0 of 23** — none. CNVkit was run with `--diagram --scatter` but **no `--vcf`**.
- Tumor `recal.bam`s archived (needed to re-run Sequenza/ASCAT): **0**. The only archived BAM is a single *public normal* (SRR4039087.bam; PRJNA339046). The tumors' recal BAMs live on HPC scratch (`/scratch/asmab/...`), not in the archive.

## scarHRD tooling note (secondary)
- `scarHRD` installed in this environment: **FALSE**.
- Its dependency **`copynumber`** is installed locally (Bioconductor 3.21): **FALSE**.
- Install attempt outcome: `not attempted: software installation cannot resolve the missing allele-specific input`.
- This is a *secondary* obstacle only; even with scarHRD installed, the **data blocker above is dispositive**.

## Recommended path to a genuine genomic HRD score
1. Recover the tumors' recalibrated BAMs (nf-core/Sarek `recal.bam`, HPC scratch) — the CNV/SNV inputs.
2. Generate **allele-specific CN** per line with **Sequenza** (or ASCAT/FACETS): call heterozygous SNP BAFs (a population SNP panel suffices for tumor-only; a matched normal is better) + depth ratio, fit cellularity/ploidy, emit segments with `A_cn`/`B_cn`.
3. Run **`scarHRD`** on the Sequenza segments to obtain HRD-LOH + TAI + LST and the HRD-sum.
4. Interpretation caveat (Takamatsu 2024, *Sci Data*): in cell lines, HRD scars **persist** and do **not** predict in-vitro platinum/PARP sensitivity — report HRD as a genomic-scar descriptor, not a drug-response predictor.

## Note on the dropped expression 'HRD'
The Peng et al. 2014 expression signature previously labelled 'HRD' is a **transcriptional state, not a genomic scar** (published transcriptomic HRD signatures barely overlap and none is guideline-endorsed). It is **dropped**; 'HRD' in this resource must mean the genomic-scar score above.

