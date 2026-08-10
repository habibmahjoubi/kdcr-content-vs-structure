#!/usr/bin/env bash
# Confidence-threshold sensitivity sweep: reclassify all 9 case-study runs against all 13
# MSL-version Kraken2 databases at a non-default --confidence threshold, mirroring
# classify_all_versions.sh exactly except for the added --confidence flag.
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kdcr
cd ~/kdcr

CONF="$1"           # e.g. 0.1 or 0.5
OUTTAG="conf_${CONF}"

RUNS=(SRR26352217 SRR26352207 SRR26352216 SRR26352206 SRR26352195 SRR26352212 SRR26352205 SRR26352208 SRR26352202)
DBS=(msl29 msl30 msl31 msl32 msl33 msl34 msl35 msl36 msl37 msl38 msl39 msl40 msl41)

for db in "${DBS[@]}"; do
  mkdir -p "results/${OUTTAG}/${db}"
  for run in "${RUNS[@]}"; do
    if [ -f "results/${OUTTAG}/${db}/${run}_report.txt" ]; then
      continue
    fi
    kraken2 --db "kraken2_dbs/${db}" --paired --threads 4 --confidence "${CONF}" \
      --report "results/${OUTTAG}/${db}/${run}_report.txt" \
      --output "results/${OUTTAG}/${db}/${run}_output.txt" \
      "data/fastq/${run}_1.fastq" "data/fastq/${run}_2.fastq" >/dev/null 2>&1
  done
  echo "[done] db=${db} conf=${CONF}"
done
echo "ALL DONE conf=${CONF}"
