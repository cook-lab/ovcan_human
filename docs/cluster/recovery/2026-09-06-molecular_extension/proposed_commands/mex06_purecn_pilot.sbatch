#!/bin/bash
#SBATCH --account=def-dcook_cpu
#SBATCH --job-name=mex06_purecn
#SBATCH --array=1-5
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%A_%a.out
# MEX06 tumour-only allele-specific CN pilot with PureCN 2.14.1 (Bioconductor 3.21, R 4.5.0)
# importing the verified target-only CNVkit log-ratios (inputs/target_only_cnr/<S>.targetonly.cnr,
# hashes == repo models.tsv target_only_cnr_sha256) and the repository's corrected target-only
# CNS (output/wes_cnv_target_only/<cell_line>.cns) as a SEG file, with the MEX05 germline-site
# Mutect2 VCF. Segmentation is NOT redone (--fun-segmentation none). Two fits per model:
#   fit A: broad purity/ploidy grid (default ploidy 1.5-6, purity 0.30-0.99)
#   fit B: high-purity / high-copy sensitivity (purity 0.80-0.99, max copy number 16)
# Pilot models: TOV81D_P23 OV2085_P64 TOV2929D_P57 OV3331_P54 OV1369_R2_P66.
# Requires MEX05 outputs. NOT SUBMITTED until authorised.
set -euo pipefail
module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.5     # optparse 1.8.2 installed 2026-09-06
BASE=/project/6090753/active/ovcan_human/wes
REPO=$HOME/ovcan_human
PCN=/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcc12/r-bundle-bioconductor/3.21/PureCN/extdata
REF=$BASE/reference/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
TARGETS=$BASE/scratch_archive_extract/intervals_sorted.bed
IN=$BASE/molecular_extension_2026-09-06/inputs
RES=$BASE/molecular_extension_2026-09-06/results
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IN/ascn_pilot_models.txt")
CL=$(awk -F'\t' -v s="$S" '$1==s{print $2}' "$IN/models_23.tsv")
VCF=$RES/mex05_mutect2_germline_sites/$S/$S.mutect2.germline_sites.vcf.gz
STATS=$RES/mex05_mutect2_germline_sites/$S/$S.mutect2.germline_sites.vcf.gz.stats
CNR=$IN/target_only_cnr/$S.targetonly.cnr
CNS=$REPO/output/wes_cnv_target_only/$CL.cns
OUT=$RES/mex06_purecn/$S; mkdir -p "$OUT"
for f in "$VCF" "$CNR" "$CNS"; do [ -s "$f" ] || { echo "missing input $f" >&2; exit 2; }; done

# 0. Interval file (once; shared): PureCN-normalised targets from the same calling BED, GC from the same FASTA.
INT=$RES/mex06_purecn/purecn_intervals_hg38.txt
if [ ! -s "$INT" ]; then
  cut -f1-3 "$TARGETS" > "$OUT/targets3.bed"
  Rscript "$PCN/IntervalFile.R" --in-file "$OUT/targets3.bed" --fasta "$REF" --genome hg38 --out-file "$INT" --force 2>&1 | tee "$OUT/intervalfile.log"
fi

# 1. SEG from the corrected target-only CNS (ID chrom loc.start loc.end num.mark seg.mean); CNS is 0-based start.
awk -F'\t' -v id="$S" 'BEGIN{OFS="\t"; print "ID","chrom","loc.start","loc.end","num.mark","seg.mean"} NR>1{print id,$1,$2+1,$3,$7,$5}' "$CNS" > "$OUT/$S.targetonly.seg"
# (CNS columns: chromosome start end gene log2 depth probes weight; verify header at run time)
head -1 "$CNS" > "$OUT/cns_header_seen.txt"

# 2. PureCN fits
run_fit () {
  tag=$1; shift
  Rscript "$PCN/PureCN.R" --out "$OUT/$tag" --sampleid "$S" \
    --log-ratio-file "$CNR" --seg-file "$OUT/$S.targetonly.seg" --fun-segmentation none \
    --vcf "$VCF" --stats-file "$STATS" --genome hg38 --intervals "$INT" --sex F \
    --model betabin --post-optimize --out-vcf --seed 20260906 --cores "$SLURM_CPUS_PER_TASK" --force "$@" \
    2>&1 | tee "$OUT/$tag.purecn.log"
}
mkdir -p "$OUT/fitA_broad" "$OUT/fitB_highpurity"
run_fit fitA_broad        --min-purity 0.30 --max-purity 0.99 --min-ploidy 1.5 --max-ploidy 6 --max-copy-number 12
run_fit fitB_highpurity   --min-purity 0.80 --max-purity 0.99 --min-ploidy 1.5 --max-ploidy 8 --max-copy-number 16

{ echo "date_utc=$(date -u +%FT%TZ) job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID node=$(hostname) model=$S cell_line=$CL";
  Rscript -e 'cat("PureCN", as.character(packageVersion("PureCN")), R.version.string, "\n")';
  echo "cnr=$CNR sha256=$(sha256sum "$CNR" | cut -c1-64)"; echo "cns=$CNS sha256=$(sha256sum "$CNS" | cut -c1-64)"; echo "vcf=$VCF"; echo "intervals=$INT"; } > "$OUT/provenance.txt"
(cd "$OUT" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
