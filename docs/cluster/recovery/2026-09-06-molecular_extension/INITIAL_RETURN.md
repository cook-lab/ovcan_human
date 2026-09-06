# Molecular-extension handoff — first return (inventory and proposed pilots), 2026-09-06

Response to `docs/cluster/molecular_extension_2026-09-06/AGENT_TASK.md` (MEX01–MEX11), prepared
on Nibi (login node, read-only discovery; no job submitted). Companion files here:
`input_inventory.tsv` (204 verified rows), `task_status.tsv`, `target_only_cnr_check.tsv`, and
`proposed_commands/*.sbatch` (five job scripts, **not submitted**). Working area on the cluster:
`/project/6090753/active/ovcan_human/wes/molecular_extension_2026-09-06/` (inputs, scripts,
results). Nothing in the source runs or the repository's canonical outputs was modified.

## What is established (MEX01)

| Question | Answer | Evidence |
|---|---|---|
| Final alignments, all 23 | Three generations exist per model, all GRCh38 (`Homo_sapiens_assembly38`): **(a)** duplicate-marked genome-wide CRAM (`md.cram`, 3.1–4.0 GB, all reads incl. off-target, not BQSR'd); **(b)** published recalibrated CRAM (`recal.cram`, interval-restricted, 10 from Oct-24 and 13 from Oct-31 re-executions); **(c)** the Oct-24 recalibrated CRAM in `work/` that was converted to the BAM fed to the manual CNVkit (all 23 still present). Provider GRCh37 BAMs (23, ~154 GB) remain the Sarek input. | `input_inventory.tsv`; recovery 2026-09-05/06 |
| Reference identity | FASTA staged and hash-verified (SHA-256 `93157a16…5564`, MD5 `7ff13495…b5e`), with `.fai/.dict`; contigs `chr`-prefixed, 3,366 sequences incl. alt/decoy. | `wes/reference/WholeGenomeFasta/` |
| Interval restriction | Recalibrated CRAMs contain only reads overlapping the 405 Sarek shards (ApplyBQSR `-L`, padding 0). The `md.cram` is complete. Off-target/flanking evidence therefore exists **only** in `md.cram`. | `task_evidence/oct23_GATK4_APPLYBQSR_*`, samtools stats md vs recal |
| Target/bait definitions | `intervals_sorted.bed` (242,421; calling), CNVkit target bins (290,475) and antitarget bins (42,709); SeqCap EZ Exome v3 vendor targets/probes. | hashes in `input_inventory.tsv` |
| Read groups / mapping | `ID:<sample>_<lane> PU:<lane> SM:<sample>_<sample> LB:<sample>`; lanes 5–8 per model; 23 models ↔ 16 patients reconciled from `models.tsv`. | recovery 2026-09-05 `model_status.tsv` |
| Existing germline-aware VCF / allele counts | **None.** The executed Mutect2 used `--genotype-germline-sites false --genotype-pon-sites false --interval-padding 0` (VCF header). No CollectAllelicCounts, pileup-at-SNP or ASCAT/FACETS/Sequenza/PureCN output exists in the Sarek results, `work/`, the manual CNVkit tree or the scratch tar. | grep of tar index and results tree |
| Existing MSI results | **None.** The 2022 McGill PCGR run had `msi.run = false` ("Not determined"); Sarek never ran msisensorpro. | PCGR JSON `metadata/config/msi` |
| Matched donor normals | None found anywhere. The five public exomes are capture-matched CNV references (bwa-mem2, no duplicate marking, no BQSR); the GATK 1000G PoN is the variant panel. | `cnvkit/normal_samples/` scripts |
| CNV inputs for allele-specific pilot | All 23 cluster CNR files hash-match the repository's `source_cnr_sha256`; removing `gene == Antitarget` rows reproduces the repository's `target_only_cnr_sha256` for **23/23**. Target-only CNRs now exist on the cluster; corrected CNS are in the clone (`output/wes_cnv_target_only/`). | `target_only_cnr_check.tsv` |
| Why antitargets are empty | Processing artefact: the manual CNVkit inputs were interval-restricted recal BAMs, so off-target reads were absent; the normals kept 0.45–0.63× off-target depth. | recovery 2026-09-06 |

## Software resolved on Nibi

| Need | Resolution |
|---|---|
| MSIsensor2 | no module; biocontainer `msisensor2:0.1--h077b44d_4` pulled (`wes/containers/`, sha256 `b8eca5ed…1eb5`); pretrained **hg38 models (2,829 files, `chr`-prefixed)** staged from the tar (`wes/reference/MSIsensor2/models_hg38`) |
| GATK 4.5.0.0 (same as original run) | biocontainer `gatk4:4.5.0.0--py36hdfd78af_0` pulled; module `gatk/4.6.1.0` exists but requires an interactive licence prompt |
| Mutect2 resources | `1000g_pon.hg38.vcf.gz` and `af-only-gnomad.hg38.vcf.gz` staged from the tar, hash-verified |
| PureCN | **PureCN 2.14.1** in `r-bundle-bioconductor/3.21` with `r/4.5.0` (DNAcopy 1.82.0); `optparse 1.8.2` installed to `~/R` so the `PureCN.R`/`IntervalFile.R`/`NormalDB.R` wrappers run |
| CNVkit 0.9.10 | biocontainer pulled (wheelhouse CNVkit 0.9.6 lacks `pomegranate`) |
| samtools / bcftools | modules 1.22.1 / 1.22 |
| scarHRD | not available; would need a pinned GitHub install into `~/R` (MEX07, only after credible fits) |
| Slurm | account `def-dcook_cpu`; accounting retention is days only, so job records must be saved by the scripts themselves |

## Proposed pilots (each requires your authorisation before `sbatch`)

| Script | Task | What it does | Inputs | Resources |
|---|---|---|---|---|
| `mex02_03_locus_review.sbatch` | MEX02, MEX03, MEX08 | For the 25 review regions: region BAM slices, `bcftools mpileup` with AD/ADF/ADR/DP, per-base depth, raw alignment fields (for strand/MAPQ/read-position parsing). For the 25 model-locus requests: per-target `samtools bedcov` for the model **and the five normals**, per-base depth over each gene span. | `md.cram` ×13/12 models, normal BAMs, FASTA, calling BED | 1 job, 2 CPU, 8 GB, < 2 h |
| `mex04_msisensor2_pilot.sbatch` | MEX04 | CRAM→BAM in `$SLURM_TMPDIR`, `msisensor2 msi -M models_hg38 -c 20 -b 4`, keep score, `_dis`, `_somatic`, flagstat and provenance. | `md.cram` of TOV21G, TOV2414, TOV3392D, OV3331 | array 4 × (4 CPU, 16 GB, ≤ 2.75 h) |
| `mex05_mutect2_germline_sites.sbatch` | MEX05 | Mutect2 4.5.0.0 tumour-only with `--genotype-germline-sites true --genotype-pon-sites true --interval-padding 75`, same PoN/gnomAD/intervals; plus orientation model, pileups, contamination and FilterMutectCalls companions. | the **Oct-24 recal CRAM** (exact CNVkit input) of TOV81D, OV2085, TOV2929D, OV3331, OV1369-R2 | array 5 × (4 CPU, 16 GB, ≤ 7 h; expected 3–6 h unsharded) |
| `mex05b_normals_mutect2_for_mappingbias.sbatch` | MEX05/06 optional | Same mode on the five public normals → `NormalDB.R` mapping-bias RDS (process-matched approximation, clearly labelled). | normal sorted BAMs | array 5 × (4 CPU, 16 GB, ≤ 7 h) |
| `mex06_purecn_pilot.sbatch` | MEX06 | `IntervalFile.R` from the calling BED + FASTA; SEG from the corrected target-only CNS; `PureCN.R --log-ratio-file <target-only CNR> --seg-file <SEG> --fun-segmentation none --vcf <MEX05> --stats-file … --genome hg38 --model betabin --post-optimize --out-vcf --seed 20260906`; **two fits per model**: broad grid (purity 0.30–0.99, ploidy 1.5–6, max CN 12) and high-purity/high-CN sensitivity (purity 0.80–0.99, ploidy 1.5–8, max CN 16). | MEX05 VCFs + verified CNR/CNS | array 5 × (4 CPU, 16 GB, ≤ 3 h), after MEX05 |

Total if everything is approved: about 20 array tasks, none above 4 cores or 16 GB, longest 7 h;
well inside the account's normal QOS. All outputs go to
`wes/molecular_extension_2026-09-06/results/<task>/<model>/` with per-directory `SHA256SUMS` and
a `provenance.txt` (job id, node, input paths/sizes/mtimes, versions, exact command), because
Slurm accounting on Nibi will not retain the job records.

## Decisions I need from you

1. **Alignment choice.** I propose `md.cram` (complete reads) for MSI and read review, and the
   Oct-24 recal CRAM (exact CNVkit input) for MEX05/MEX06 so SNP allele counts and log-ratios
   come from the same reads. Alternative: `md.cram` for MEX05 too (adds flanking SNPs beyond
   ~100 bp of targets, at the cost of un-recalibrated base qualities and a different read set
   from the CNR). Say which.
2. **Which scripts to submit**, and whether the optional normals run (mapping bias) is wanted.
3. **MEX07 (scarHRD)**: install from a pinned GitHub commit into `~/R` only if MEX06 fits pass QC.
   Nothing else on Nibi provides LOH/TAI/LST.
4. **MEX10/MEX11** need locations outside the WES roots (RNA FASTQs or alignments; proteomics
   PSM/reporter files). I did not search for them.

## Limits and caveats carried into the pilots

- No MSI-positive or MSI-negative control alignment exists in the searched roots; the
  MSIsensor2 pilot returns continuous scores and site yields without calibration.
- No matched donor normals: PureCN runs tumour-only; the five public normals can supply mapping
  bias only, and were processed differently (bwa-mem2, no MarkDuplicates/BQSR).
- Recalibrated CRAMs cannot supply off-target or distant flanking SNPs; padded calling on them
  yields reads only within ~read length of targets.
- The TOV81D provider alias (`TPV81D_23-pool`) is unresolved on the cluster; its WES sample id
  `TOV81D_P23` is used as-is.
- Earlier 2026-09-06 findings (workflow-revision verification, reference hashes, antitarget
  cause, missing manual-CNVkit logs) are in `docs/cluster/recovery/2026-09-06/FOLLOWUP_INVENTORY.md`
  in this clone and answer several items still listed as outstanding upstream; they have not
  been pushed.
