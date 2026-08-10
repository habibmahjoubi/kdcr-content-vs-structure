#!/usr/bin/env bash
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kdcr
cd ~/kdcr

RUNS=(SRR26352217 SRR26352207 SRR26352216 SRR26352206 SRR26352195 SRR26352212 SRR26352205 SRR26352208 SRR26352202)
DBS=(msl29 msl30 msl31 msl32 msl33 msl34 msl35 msl36 msl37 msl38 msl39 msl40 msl41)

for db in "${DBS[@]}"; do
  mkdir -p "results_fixed/${db}"
  for run in "${RUNS[@]}"; do
    if [ -f "results_fixed/${db}/${run}_report.txt" ]; then
      echo "[skip] ${db} / ${run} already classified"
      continue
    fi
    echo "=== ${db} / ${run} ==="
    kraken2 --db "kraken2_dbs_fixed/${db}" --paired --threads 4 \
      --report "results_fixed/${db}/${run}_report.txt" \
      --output "results_fixed/${db}/${run}_output.txt" \
      "data/fastq/${run}_1.fastq" "data/fastq/${run}_2.fastq" 2>&1 | tail -2
  done
done

echo "All fixed-library classifications complete."
