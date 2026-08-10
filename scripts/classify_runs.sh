#!/usr/bin/env bash
# Classify all 9 subsampled runs against the current Kraken2 viral database, then summarize
# RSV/REO read counts per run for the KDCR case study (dilution-series signal check).
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kdcr
cd ~/kdcr

RUNS=(
  SRR26352217
  SRR26352207
  SRR26352216
  SRR26352206
  SRR26352195
  SRR26352212
  SRR26352205
  SRR26352208
  SRR26352202
)

mkdir -p results/current_db

for run in "${RUNS[@]}"; do
  echo "=== $run ==="
  kraken2 --db kraken2_dbs/current --paired --threads 4 \
    --report "results/current_db/${run}_report.txt" \
    --output "results/current_db/${run}_output.txt" \
    "data/fastq/${run}_1.fastq" "data/fastq/${run}_2.fastq" 2>&1 | tail -3
done

echo "=== Summary: RSV / REO reads per run ==="
for run in "${RUNS[@]}"; do
  rsv=$(grep -i -E "syncytial" "results/current_db/${run}_report.txt" | awk -F'\t' '{sum+=$3} END{print sum+0}')
  reo=$(grep -i -E "orthoreovirus|reoviridae" "results/current_db/${run}_report.txt" | awk -F'\t' '{sum+=$3} END{print sum+0}')
  echo "$run: RSV_reads=$rsv REO_reads=$reo"
done
