#!/usr/bin/env bash
set -euo pipefail
BASE=~/kdcr
FIXED_FASTA="$BASE/kraken2_dbs_fixed/library_fixed.fna"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kdcr

for msl in msl29 msl30 msl31 msl32 msl33 msl34 msl35 msl36 msl37 msl38 msl39 msl40 msl41; do
  dbdir="$BASE/kraken2_dbs_fixed/$msl"
  echo "=== $msl ==="
  if [ -f "$dbdir/hash.k2d" ]; then
    echo "[skip] already built"
    continue
  fi
  mkdir -p "$dbdir/taxonomy"
  cp "$BASE/kraken2_dbs/$msl/taxonomy/nodes.dmp" "$dbdir/taxonomy/nodes.dmp"
  cp "$BASE/kraken2_dbs/$msl/taxonomy/names.dmp" "$dbdir/taxonomy/names.dmp"
  kraken2-build --add-to-library "$FIXED_FASTA" --db "$dbdir" > "$BASE/logs/fixed_add_${msl}.log" 2>&1
  kraken2-build --build --db "$dbdir" --threads 4 > "$BASE/logs/fixed_build_${msl}.log" 2>&1
  tail -5 "$BASE/logs/fixed_build_${msl}.log"
done
echo 'All fixed-library databases built.'
