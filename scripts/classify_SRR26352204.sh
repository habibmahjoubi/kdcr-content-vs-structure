#!/usr/bin/env bash
# Classify SRR26352204 (RSV 10^5), which had been classified only with the alternative-pairing
# databases, against the fixed-content, merge-corrected and taxonomy-filtered (Design 1) databases,
# then report the cumulative counts at the RSV (11250) and MRV (351073) nodes.
set -uo pipefail
K=/home/hm/kdcr
KR=/home/hm/miniforge3/envs/kdcr/bin/kraken2  # v2.17.1, as in the study
OUT=$K/work/SRR26352204
R1=$K/data/fastq/SRR26352204_1.fastq
R2=$K/data/fastq/SRR26352204_2.fastq
TSV="$(dirname "$(readlink -f "$0")")/../results/SRR26352204_all_designs.tsv"
mkdir -p "$OUT"
echo -e "design\tversion\treads_RSV_11250\treads_MRV_351073" > "$TSV"
for design in fixed fixed_corrected design1; do
  case $design in
    fixed) base=$K/kraken2_dbs_fixed ;;
    fixed_corrected) base=$K/kraken2_dbs_fixed_corrected ;;
    design1) base=$K/kraken2_dbs ;;
  esac
  mkdir -p "$OUT/$design"
  for n in $(seq 29 41); do
    m=msl$n; rep="$OUT/$design/${m}_report.txt"
    if [ ! -s "$rep" ]; then
      $KR --db "$base/$m" --memory-mapping --paired --threads 4 --output /dev/null --report "$rep" "$R1" "$R2" 2>/dev/null
    fi
    rsv=$(awk -F'\t' '$5==11250{print $2}' "$rep"); mrv=$(awk -F'\t' '$5==351073{print $2}' "$rep")
    echo -e "$design\t$m\t${rsv:-0}\t${mrv:-0}" | tee -a "$TSV"
  done
done
