#!/usr/bin/env bash
# Positive-control test for KDCR: Canine distemper virus (NC_001921.1, taxid 3139435) is a
# real viral genome whose taxid was newly added to NCBI taxonomy between the MSL39 (end 2023)
# and MSL40 (end 2024) snapshots -- part of the documented Mononegavirales/Paramyxoviridae
# binomial-nomenclature reorganization. Unlike RSV/REO (stable across all 5 versions), this
# taxon's classification SHOULD change sharply across that specific transition, since older
# reference databases literally do not contain it (excluded by our taxid-existence filter).
set -euo pipefail
BASE=~/kdcr
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kdcr
cd "$BASE"

mkdir -p positive_control

# 1. Extract the target genome sequence from the bulk tagged FASTA
python3 - <<'PYEOF'
target = "NC_001921.1"
in_fasta = "kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna"
out_fasta = "positive_control/canine_distemper.fna"
write = False
with open(in_fasta, encoding="utf-8", errors="replace") as fin, open(out_fasta, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            write = line[1:].split("|")[0] == target
            if write:
                fout.write(f">{target}\n")
            continue
        if write:
            fout.write(line)
PYEOF
echo "Extracted genome:"
grep -c '^>' positive_control/canine_distemper.fna
wc -l positive_control/canine_distemper.fna

# 2. Simulate paired-end reads (5000 pairs, 150bp, realistic error rate) with wgsim
wgsim -N 5000 -1 150 -2 150 -e 0.01 -r 0 \
  positive_control/canine_distemper.fna \
  positive_control/cdv_sim_1.fastq positive_control/cdv_sim_2.fastq \
  > positive_control/wgsim.log 2>&1
echo "Simulated reads:"
echo "$(( $(wc -l < positive_control/cdv_sim_1.fastq) / 4 )) read pairs"

# 3. Classify the simulated reads against all 5 MSL databases
for db in msl37 msl38 msl39 msl40 msl41; do
  echo "=== $db ==="
  kraken2 --db "kraken2_dbs/$db" --paired --threads 4 \
    --report "positive_control/${db}_report.txt" \
    --output "positive_control/${db}_output.txt" \
    positive_control/cdv_sim_1.fastq positive_control/cdv_sim_2.fastq 2>&1 | tail -3
done

echo "Positive-control test complete."
