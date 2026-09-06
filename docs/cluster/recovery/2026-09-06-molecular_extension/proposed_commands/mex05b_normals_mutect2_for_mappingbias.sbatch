#!/bin/bash
#SBATCH --account=def-dcook_cpu
#SBATCH --job-name=mex05b_normals
#SBATCH --array=1-5
#SBATCH --time=07:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%A_%a.out
# OPTIONAL companion to MEX05/MEX06: run Mutect2 in the same germline-site-retaining mode on the
# five public CNV-reference exomes so PureCN's NormalDB.R can build a mapping-bias file
# (calculateMappingBiasVcf). These are capture-matched public normals aligned with bwa-mem2,
# NOT matched donor normals and NOT the variant panel of normals; the resulting mapping bias
# is a process-matched approximation and must be labelled as such. NOT SUBMITTED until authorised.
set -euo pipefail
module load StdEnv/2023 apptainer/1.4.5
BASE=/project/6090753/active/ovcan_human/wes
REF=$BASE/reference/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
PON=$BASE/reference/GATKBundle/1000g_pon.hg38.vcf.gz
GNOMAD=$BASE/reference/GATKBundle/af-only-gnomad.hg38.vcf.gz
TARGETS=$BASE/scratch_archive_extract/intervals_sorted.bed
IMG=$BASE/containers/gatk4_4.5.0.0--py36hdfd78af_0.sif
N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<< $'SRR4039087\nSRR4039088\nSRR4039089\nSRR4039096\nSRR4039097')
BAM=/project/6090753/active/ovcan_gq_analysis/cnvkit/normal_samples/$N/$N.sorted.bam
OUT=$BASE/molecular_extension_2026-09-06/results/mex05b_normals_mutect2/$N; mkdir -p "$OUT"; cd "$SLURM_TMPDIR"
cut -f1-3 "$TARGETS" > targets3.bed
apptainer exec -B "$SLURM_TMPDIR","$BASE",/project/6090753/active/ovcan_gq_analysis "$IMG" gatk --java-options "-Xmx12g -XX:-UsePerfData" Mutect2 \
  -R "$REF" -I "$BAM" -L targets3.bed --interval-padding 75 \
  --panel-of-normals "$PON" --germline-resource "$GNOMAD" \
  --genotype-germline-sites true --genotype-pon-sites true --native-pair-hmm-threads "$SLURM_CPUS_PER_TASK" \
  -O "$OUT/$N.mutect2.germline_sites.vcf.gz" 2>&1 | tee "$OUT/$N.mutect2.log"
{ echo "date_utc=$(date -u +%FT%TZ) job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID normal=$N"; echo "input_bam=$BAM bytes=$(stat -c %s "$BAM")";
  echo "note=bwa-mem2 alignment without duplicate marking or BQSR (see cnvkit/normal_samples scripts)"; } > "$OUT/provenance.txt"
(cd "$OUT" && sha256sum -- * > SHA256SUMS)
# Afterwards (single short step, login node acceptable):
#   bcftools merge -Oz -o normals.mutect2.germline_sites.vcf.gz results/mex05b_normals_mutect2/*/*.vcf.gz && bcftools index -t ...
#   Rscript $PCN/NormalDB.R --normal-panel normals.mutect2.germline_sites.vcf.gz --genome hg38 --assay SeqCapEZv3 --out-dir results/mex06_purecn/normaldb
#   then add --mapping-bias-file results/mex06_purecn/normaldb/mapping_bias_SeqCapEZv3_hg38.rds to PureCN.R.
