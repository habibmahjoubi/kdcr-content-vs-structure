#!/usr/bin/env bash
# A2: classify the 4 large (10M-pair) runs against all 13 fixed-library (confound-controlled)
# Kraken2 databases, exactly as done for the original 2M-pair runs.
set -euo pipefail
cd ~/kdcr
mkdir -p results_fixed_large

RUNS=(SRR26352206 SRR26352205 SRR26352202 SRR26352216)
VERSIONS=(msl29 msl30 msl31 msl32 msl33 msl34 msl35 msl36 msl37 msl38 msl39 msl40 msl41)

for db in "${VERSIONS[@]}"; do
  mkdir -p "results_fixed_large/${db}"
  for run in "${RUNS[@]}"; do
    r1="data/fastq_large/${run}_1.fastq"
    r2="data/fastq_large/${run}_2.fastq"
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
      echo "  [skip] $run reads not found yet"
      continue
    fi
    out="results_fixed_large/${db}/${run}_output.txt"
    rep="results_fixed_large/${db}/${run}_report.txt"
    if [ -f "$rep" ]; then
      echo "  [skip] $db/$run already classified"
      continue
    fi
    echo "=== $db / $run ==="
    kraken2 --db "kraken2_dbs_fixed/${db}" --paired --threads 6 \
      --report "$rep" --output "$out" "$r1" "$r2" 2>&1 | tail -3
  done
done
echo "All large-run classification complete."
