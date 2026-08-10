#!/usr/bin/env bash
# A2: retrieve larger subsamples (10,000,000 read pairs, 5x the original 2,000,000) for the
# 4 runs that carry any non-zero signal in the original study (RSV 1e6/1e3/1e1, REO 1e6),
# to increase the dynamic range of the classification-count outcome variable.
set -euo pipefail
cd ~/kdcr/data
mkdir -p fastq_large raw_chunks_large

N_READS=10000000
N_LINES=$((N_READS * 4))
RANGE_BYTES=800000000   # 800MB compressed chunk, ~5.3x safety margin over the 150MB/2M-read ratio

RUNS=(
  SRR26352206
  SRR26352205
  SRR26352202
  SRR26352216
)

get_ena_url () {
  local run="$1" mate="$2"
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp" \
    | tail -n1 | cut -f2 | tr ';' '\n' | grep "_${mate}.fastq.gz"
}

for run in "${RUNS[@]}"; do
  if [ -f "fastq_large/${run}_1.fastq" ] && [ -f "fastq_large/${run}_2.fastq" ]; then
    n1=$(( $(wc -l < "fastq_large/${run}_1.fastq") / 4 ))
    echo "[skip] $run already extracted ($n1 read pairs)"
    continue
  fi
  echo "=== $run ==="
  for mate in 1 2; do
    url_path=$(get_ena_url "$run" "$mate")
    url="https://${url_path}"
    echo "  mate $mate: $url"
    curl -s -r 0-${RANGE_BYTES} -o "raw_chunks_large/${run}_${mate}.gz.part" "$url"
    { zcat -f "raw_chunks_large/${run}_${mate}.gz.part" 2>/dev/null || true; } | head -n "$N_LINES" > "fastq_large/${run}_${mate}.fastq"
    got=$(( $(wc -l < "fastq_large/${run}_${mate}.fastq") / 4 ))
    echo "  mate $mate: extracted $got read records"
    rm -f "raw_chunks_large/${run}_${mate}.gz.part"
  done
done

echo "All large runs processed. Sizes:"
ls -la fastq_large/
