#!/bin/bash
#SBATCH --account=def-dcook_cpu
#SBATCH --job-name=mex02_03_locus
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=%x-%j.out
# MEX02 (25 variant review regions) and MEX03/MEX08 (25 model-locus requests):
# read-level and per-target depth evidence from the duplicate-marked genome-wide CRAMs,
# plus per-target depth of the five public CNV-reference normals over the same loci.
# Produces small region BAM slices, mpileup VCFs with AD/ADF/ADR/DP, per-base depth and
# per-target bedcov tables. Read-only on all sources. NOT SUBMITTED until authorised.
set -euo pipefail
module load StdEnv/2023 samtools/1.22.1 bcftools/1.22
BASE=/project/6090753/active/ovcan_human/wes
REF=$BASE/reference/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
TARGETS=$BASE/scratch_archive_extract/intervals_sorted.bed
IN=$BASE/molecular_extension_2026-09-06/inputs
OUT=$BASE/molecular_extension_2026-09-06/results/mex02_03_locus_review
MD=/project/6090753/active/ovcan_gq_analysis/wes_sarek_results_20251022/preprocessing/markduplicates
NORMALS=/project/6090753/active/ovcan_gq_analysis/cnvkit/normal_samples/all_bams
mkdir -p "$OUT"/slices "$OUT"/pileup "$OUT"/depth "$OUT"/bedcov

# cell_line -> cnv_sample_id map
declare -A SID
while IFS=$'\t' read -r sid cl cram; do [ "$sid" = "cnv_sample_id" ] && continue; SID["$cl"]="$sid"; done < "$IN/models_23.tsv"

# ---- MEX02: variant review regions (name = <cell_line>_<gene>)
while IFS=$'\t' read -r chr s e name prio; do
  cl=${name%_*}; sid=${SID[$cl]:-}
  [ -z "$sid" ] && { echo "no sample id for $cl" >&2; continue; }
  cram=$MD/$sid/$sid.md.cram; reg="$chr:$((s+1))-$e"; tag="${sid}_${name#*_}_${chr}_$((s+1))"
  samtools view -T "$REF" -b -o "$OUT/slices/$tag.bam" "$cram" "$reg"; samtools index "$OUT/slices/$tag.bam"
  bcftools mpileup -f "$REF" -r "$reg" -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,INFO/AD -Q 0 -q 0 -d 100000 -Oz -o "$OUT/pileup/$tag.mpileup.vcf.gz" "$cram"
  samtools depth -a -Q 0 -q 0 -r "$reg" "$cram" > "$OUT/depth/$tag.depth.tsv"
  # MAPQ / strand / read-position context: raw alignments for later pysam parsing
  samtools view "$cram" "$reg" | cut -f1-9 > "$OUT/slices/$tag.reads.tsv"
done < "$IN/bam_review_regions.bed"

# ---- MEX03/MEX08: gene loci, per-target depth for the model and the five normals
while IFS=$'\t' read -r chr s e name task; do
  cl=${name%_*}; gene=${name#*_}; sid=${SID[$cl]:-}
  [ -z "$sid" ] && { echo "no sample id for $cl" >&2; continue; }
  awk -v c="$chr" -v s="$s" -v e="$e" '$1==c && $3>s && $2<e' "$TARGETS" | cut -f1-4 > "$OUT/bedcov/${sid}_${gene}.targets.bed"
  samtools bedcov -c "$OUT/bedcov/${sid}_${gene}.targets.bed" "$MD/$sid/$sid.md.cram" --reference "$REF" > "$OUT/bedcov/${sid}_${gene}.model_md.bedcov.tsv"
  for n in SRR4039087 SRR4039088 SRR4039089 SRR4039096 SRR4039097; do
    samtools bedcov -c "$OUT/bedcov/${sid}_${gene}.targets.bed" "$NORMALS/$n.sorted.bam" > "$OUT/bedcov/${sid}_${gene}.normal_$n.bedcov.tsv"
  done
  samtools depth -a -Q 0 -q 0 -r "$chr:$((s+1))-$e" "$MD/$sid/$sid.md.cram" | gzip > "$OUT/depth/${sid}_${gene}.perbase.depth.tsv.gz"
done < "$IN/locus_review_regions.bed"

{ echo "date_utc=$(date -u +%FT%TZ) job=$SLURM_JOB_ID node=$(hostname)"; samtools --version | head -1; bcftools --version | head -1;
  echo "reference=$REF"; echo "targets=$TARGETS"; echo "model_cram_dir=$MD (duplicate-marked, genome-wide)"; echo "normals=$NORMALS (bwa-mem2 sorted BAMs)"; } > "$OUT/provenance.txt"
(cd "$OUT" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
