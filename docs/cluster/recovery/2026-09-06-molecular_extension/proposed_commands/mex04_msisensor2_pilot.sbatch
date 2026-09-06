#!/bin/bash
#SBATCH --account=def-dcook_cpu
#SBATCH --job-name=mex04_msi
#SBATCH --array=1-4
#SBATCH --time=02:45:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%A_%a.out
# MEX04 tumour-only MSI pilot (MSIsensor2 0.1, pretrained hg38 models) on four contrasting
# models: TOV21G_P59, TOV2414_P65, TOV3392D_P69, OV3331_P54 (inputs/msi_pilot_models.txt).
# Input: Sarek duplicate-marked, genome-wide CRAM (not the interval-restricted recal CRAM).
# NOT SUBMITTED until the operator authorises execution.
set -euo pipefail
module load StdEnv/2023 apptainer/1.4.5 samtools/1.22.1
BASE=/project/6090753/active/ovcan_human/wes
REF=$BASE/reference/WholeGenomeFasta/Homo_sapiens_assembly38.fasta          # sha256 93157a16...5564
MODELS=$BASE/reference/MSIsensor2/models_hg38                                 # 2,829 pretrained xgboost site models (GATK hg38 iGenomes bundle)
IMG=$BASE/containers/msisensor2_0.1--h077b44d_4.sif                           # sha256 b8eca5ed...1eb5
LIST=$BASE/molecular_extension_2026-09-06/inputs/msi_pilot_models.txt
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
CRAM=/project/6090753/active/ovcan_gq_analysis/wes_sarek_results_20251022/preprocessing/markduplicates/$S/$S.md.cram
OUT=$BASE/molecular_extension_2026-09-06/results/mex04_msisensor2/$S
mkdir -p "$OUT"
cd "$SLURM_TMPDIR"

# 1. CRAM -> BAM on node-local disk (MSIsensor2 reads BAM); no read filtering applied here.
samtools view -@ "$SLURM_CPUS_PER_TASK" -T "$REF" -b -o "$S.md.bam" "$CRAM"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$S.md.bam"
samtools flagstat -@ "$SLURM_CPUS_PER_TASK" "$S.md.bam" > "$OUT/$S.md.flagstat.txt"

# 2. MSIsensor2 tumour-only: -c 20 (WXS default coverage threshold), -b threads.
#    Outputs: <prefix> (score), <prefix>_dis (per-site distributions), <prefix>_somatic.
apptainer exec -B "$SLURM_TMPDIR","$BASE" "$IMG" \
  msisensor2 msi -M "$MODELS" -t "$S.md.bam" -o "$OUT/$S.msisensor2" -c 20 -b "$SLURM_CPUS_PER_TASK" \
  2>&1 | tee "$OUT/$S.msisensor2.log"

# 3. Provenance record.
{
  echo "date_utc=$(date -u +%FT%TZ)"
  echo "slurm_job=$SLURM_JOB_ID array_task=$SLURM_ARRAY_TASK_ID node=$(hostname)"
  echo "model=$S"
  echo "input_cram=$CRAM bytes=$(stat -c %s "$CRAM") mtime=$(stat -c %y "$CRAM")"
  echo "reference=$REF"
  echo "models_dir=$MODELS n_model_files=$(ls "$MODELS" | wc -l)"
  echo "image=$IMG"
  echo "command=msisensor2 msi -M models_hg38 -t $S.md.bam -o $S.msisensor2 -c 20 -b $SLURM_CPUS_PER_TASK"
  samtools --version | head -1
  apptainer --version
} > "$OUT/$S.provenance.txt"
(cd "$OUT" && sha256sum -- * > SHA256SUMS)
