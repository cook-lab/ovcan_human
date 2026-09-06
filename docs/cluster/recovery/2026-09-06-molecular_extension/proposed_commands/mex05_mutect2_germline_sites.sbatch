#!/bin/bash
#SBATCH --account=def-dcook_cpu
#SBATCH --job-name=mex05_m2gl
#SBATCH --array=1-5
#SBATCH --time=07:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%A_%a.out
# MEX05: regenerate a germline-SNP-retaining Mutect2 VCF for the allele-specific CN pilot.
# Same GATK 4.5.0.0 as the original run (biocontainer), same reference, PoN and germline
# resource, same calling intervals, but with --genotype-germline-sites/--genotype-pon-sites
# true and 75 bp interval padding. Input: the Oct-24 recalibrated CRAM in work/ that was the
# exact CNVkit input (so SNP allele counts and CNVkit log-ratios derive from the same reads).
# Pilot models: TOV81D_P23 OV2085_P64 TOV2929D_P57 OV3331_P54 OV1369_R2_P66.
# Note: the recal CRAMs hold only reads overlapping the calling intervals; flanking SNP
# coverage tapers within ~100 bp of targets. NOT SUBMITTED until authorised.
set -euo pipefail
module load StdEnv/2023 apptainer/1.4.5 samtools/1.22.1
BASE=/project/6090753/active/ovcan_human/wes
REF=$BASE/reference/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
PON=$BASE/reference/GATKBundle/1000g_pon.hg38.vcf.gz
GNOMAD=$BASE/reference/GATKBundle/af-only-gnomad.hg38.vcf.gz
TARGETS=$BASE/scratch_archive_extract/intervals_sorted.bed
IMG=$BASE/containers/gatk4_4.5.0.0--py36hdfd78af_0.sif
LIST=$BASE/molecular_extension_2026-09-06/inputs/ascn_pilot_models.txt
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
# exact CNVkit input alignment: resolved from the Oct-24 CNVKIT_BATCH task directory listing
CRAM=$(grep -h "$S.recal.cram ->" $BASE/sarek_3.5.1_run/task_evidence/oct24_CNVKIT_BATCH_$S.*/_dir_listing.txt | sed -E 's/.*-> //')
[ -s "$CRAM" ] || { echo "input CRAM not resolved for $S" >&2; exit 2; }
OUT=$BASE/molecular_extension_2026-09-06/results/mex05_mutect2_germline_sites/$S
mkdir -p "$OUT"; cd "$SLURM_TMPDIR"
# 3-column BED for GATK -L (drop name/score columns)
cut -f1-3 "$TARGETS" > targets3.bed

apptainer exec -B "$SLURM_TMPDIR","$BASE",/project/6090753/active/ovcan_gq_analysis "$IMG" gatk --java-options "-Xmx12g -XX:-UsePerfData" Mutect2 \
  -R "$REF" -I "$CRAM" -L targets3.bed --interval-padding 75 \
  --panel-of-normals "$PON" --germline-resource "$GNOMAD" \
  --genotype-germline-sites true --genotype-pon-sites true \
  --native-pair-hmm-threads "$SLURM_CPUS_PER_TASK" \
  --f1r2-tar-gz "$OUT/$S.f1r2.tar.gz" \
  -O "$OUT/$S.mutect2.germline_sites.vcf.gz" 2>&1 | tee "$OUT/$S.mutect2.log"

# Standard filtering companions (kept separate; PureCN can read the unfiltered VCF)
apptainer exec -B "$SLURM_TMPDIR","$BASE",/project/6090753/active/ovcan_gq_analysis "$IMG" gatk LearnReadOrientationModel -I "$OUT/$S.f1r2.tar.gz" -O "$OUT/$S.artifactprior.tar.gz"
apptainer exec -B "$SLURM_TMPDIR","$BASE",/project/6090753/active/ovcan_gq_analysis "$IMG" gatk GetPileupSummaries -I "$CRAM" -R "$REF" -V "$GNOMAD" -L targets3.bed --interval-padding 75 -O "$OUT/$S.pileups.table"
apptainer exec -B "$SLURM_TMPDIR","$BASE",/project/6090753/active/ovcan_gq_analysis "$IMG" gatk CalculateContamination -I "$OUT/$S.pileups.table" -tumor-segmentation "$OUT/$S.segmentation.table" -O "$OUT/$S.contamination.table"
apptainer exec -B "$SLURM_TMPDIR","$BASE",/project/6090753/active/ovcan_gq_analysis "$IMG" gatk FilterMutectCalls -R "$REF" -V "$OUT/$S.mutect2.germline_sites.vcf.gz" \
  --contamination-table "$OUT/$S.contamination.table" --tumor-segmentation "$OUT/$S.segmentation.table" \
  --orientation-bias-artifact-priors "$OUT/$S.artifactprior.tar.gz" -O "$OUT/$S.mutect2.germline_sites.filtered.vcf.gz"

{ echo "date_utc=$(date -u +%FT%TZ) job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID node=$(hostname) model=$S";
  echo "input_cram=$CRAM bytes=$(stat -c %s "$CRAM") mtime=$(stat -c %y "$CRAM")"; echo "reference=$REF"; echo "pon=$PON"; echo "germline_resource=$GNOMAD"; echo "intervals=$TARGETS padding=75";
  apptainer exec "$IMG" gatk --version 2>&1 | grep -E 'GATK v|Picard'; } > "$OUT/$S.provenance.txt"
(cd "$OUT" && sha256sum -- * > SHA256SUMS)
