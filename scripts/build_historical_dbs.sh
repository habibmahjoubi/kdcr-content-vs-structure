#!/usr/bin/env bash
# Builds one Kraken2 viral database per ICTV MSL version (37-41), each using the SAME set of
# tagged viral genome sequences (viral.1.1.genomic.taxid.fna, current-day taxid assignments)
# but a DIFFERENT historical NCBI taxonomy tree (nodes.dmp/names.dmp) matching that MSL's
# release year. Sequences whose taxid did not yet exist in a given year's taxonomy tree are
# excluded from that year's library -- this is the mechanism that lets classification results
# differ across "MSL versions" using literally the same reads and the same genome sequences.
set -euo pipefail

BASE=~/kdcr
TAGGED_FASTA="$BASE/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna"

declare -A MSL_DATE=(
  [msl29]=2014-12-01
  [msl30]=2015-12-01
  [msl31]=2016-12-01
  [msl32]=2017-12-01
  [msl33]=2018-06-01
  [msl34]=2018-12-01
  [msl35]=2019-12-01
  [msl36]=2020-12-01
  [msl37]=2021-12-01
  [msl38]=2022-12-01
  [msl39]=2023-12-01
  [msl40]=2024-12-01
  [msl41]=2025-12-01
)

source ~/miniforge3/etc/profile.d/conda.sh
conda activate kdcr

for msl in msl29 msl30 msl31 msl32 msl33 msl34 msl35 msl36 msl37 msl38 msl39 msl40 msl41; do
  date="${MSL_DATE[$msl]}"
  dbdir="$BASE/kraken2_dbs/$msl"
  echo "=== $msl (taxdmp_${date}) ==="

  if [ -f "$dbdir/hash.k2d" ]; then
    echo "[skip] $msl database already built"
    continue
  fi

  mkdir -p "$dbdir/taxonomy"
  if [ ! -f "$dbdir/taxonomy/nodes.dmp" ]; then
    unzip -o -q "$BASE/taxdump_archive/taxdmp_${date}.zip" nodes.dmp names.dmp -d "$dbdir/taxonomy"
  fi

  # Filter the tagged FASTA to only sequences whose taxid exists in THIS snapshot's nodes.dmp
  filtered_fasta="$dbdir/viral.filtered.fna"
  python3 "$BASE/scripts/filter_fasta_by_taxid_existence.py" \
    "$TAGGED_FASTA" "$dbdir/taxonomy/nodes.dmp" "$filtered_fasta" \
    > "$BASE/logs/filter_${msl}.log" 2>&1
  tail -5 "$BASE/logs/filter_${msl}.log"

  kraken2-build --add-to-library "$filtered_fasta" --db "$dbdir" > "$BASE/logs/add_${msl}.log" 2>&1
  kraken2-build --build --db "$dbdir" --threads 4 > "$BASE/logs/build_${msl}.log" 2>&1
  tail -10 "$BASE/logs/build_${msl}.log"
done

echo "All historical databases built."
