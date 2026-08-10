#!/usr/bin/env bash
# Downloads a subsample (~200,000 read pairs) of each selected public SRA run for the KDCR
# case study, via ENA's plain fastq.gz mirror + an HTTP byte-range partial fetch (avoids
# pulling each full ~8 GB mate file when only a small subsample is needed).
# BioProject PRJNA1026487 (NIIMBL/FDA CBER adventitious-virus NGS study, public data) --
# see data/sra_runs_manifest.csv for the exact runs and their role in this study.
# Run inside WSL, kdcr conda env, from ~/kdcr/data.
set -euo pipefail
cd ~/kdcr/data
mkdir -p fastq raw_chunks

N_READS=2000000
N_LINES=$((N_READS * 4))
RANGE_BYTES=150000000   # 150 MB compressed chunk; ~50MB empirically decompresses to ~1M reads at ~200bp

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

get_ena_url () {
  local run="$1" mate="$2"
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp" \
    | tail -n1 | cut -f2 | tr ';' '\n' | grep "_${mate}.fastq.gz"
}

for run in "${RUNS[@]}"; do
  if [ -f "fastq/${run}_1.fastq" ] && [ -f "fastq/${run}_2.fastq" ]; then
    echo "[skip] $run already extracted"
    continue
  fi
  echo "=== $run ==="
  for mate in 1 2; do
    url_path=$(get_ena_url "$run" "$mate")
    url="https://${url_path}"
    echo "  mate $mate: $url"
    curl -s -r 0-${RANGE_BYTES} -o "raw_chunks/${run}_${mate}.gz.part" "$url"
    { zcat -f "raw_chunks/${run}_${mate}.gz.part" 2>/dev/null || true; } | head -n "$N_LINES" > "fastq/${run}_${mate}.fastq"
    got=$(( $(wc -l < "fastq/${run}_${mate}.fastq") / 4 ))
    echo "  mate $mate: extracted $got read records"
    rm -f "raw_chunks/${run}_${mate}.gz.part"
  done
done

echo "All runs processed. Sizes:"
ls -la fastq/
