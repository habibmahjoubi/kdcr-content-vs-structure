#!/usr/bin/env bash
# Usage: download_one_run_large.sh <RUN_ACCESSION>
set -uo pipefail
cd ~/kdcr/data
mkdir -p fastq_large raw_chunks_large

RUN="$1"
N_READS=10000000
N_LINES=$((N_READS * 4))
RANGE_BYTES=800000000

get_ena_url () {
  local run="$1" mate="$2"
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp" \
    | tail -n1 | cut -f2 | tr ';' '\n' | grep "_${mate}.fastq.gz"
}

echo "=== $RUN ==="
for mate in 1 2; do
  if [ -f "fastq_large/${RUN}_${mate}.fastq" ]; then
    n=$(( $(wc -l < "fastq_large/${RUN}_${mate}.fastq") / 4 ))
    if [ "$n" -ge 1000000 ]; then
      echo "  mate $mate: already done ($n read records), skipping"
      continue
    fi
  fi
  url_path=$(get_ena_url "$RUN" "$mate")
  url="https://${url_path}"
  echo "  mate $mate: $url"
  curl -s -r 0-${RANGE_BYTES} -o "raw_chunks_large/${RUN}_${mate}.gz.part" "$url"
  echo "  mate $mate: downloaded $(stat -c%s "raw_chunks_large/${RUN}_${mate}.gz.part" 2>/dev/null || echo 0) bytes"
  zcat -f "raw_chunks_large/${RUN}_${mate}.gz.part" 2>/dev/null | head -n "$N_LINES" > "fastq_large/${RUN}_${mate}.fastq"
  got=$(( $(wc -l < "fastq_large/${RUN}_${mate}.fastq") / 4 ))
  echo "  mate $mate: extracted $got read records"
  rm -f "raw_chunks_large/${RUN}_${mate}.gz.part"
done
echo "$RUN done."
