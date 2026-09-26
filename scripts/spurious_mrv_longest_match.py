#!/usr/bin/env python3
"""Longest exact match to the MRV type 3 reference segments for every read assigned to the MRV clade
(351073) by the current database, except the one genuine MRV read (SRR26352216.1152145).
A genuine minimizer match with the default Kraken2 spaced seed needs >= 17 contiguous identical nt."""
import os, subprocess, csv
from collections import defaultdict

K = "/home/hm/kdcr"
W = f"{K}/work/allfp"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
RUNS = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195",
        "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]
GENUINE = {"SRR26352216.1152145"}
os.makedirs(W, exist_ok=True)

par = {}
for l in open(f"{K}/kraken2_dbs/current/taxonomy/nodes.dmp"):
    p = l.split("\t|\t"); par[int(p[0])] = int(p[1])
def in_mrv(t):
    while t in par:
        if t == 351073: return True
        if par[t] == t: return False
        t = par[t]
    return False

want = defaultdict(set)
for r in RUNS:
    for l in open(f"{K}/results/current_db/{r}_output.txt", errors="replace"):
        if l[0] != "C": continue
        x = l.split("\t", 3)
        if x[1] not in GENUINE and in_mrv(int(x[2])): want[r].add(x[1])
q = f"{W}/q.fa"
with open(q, "w") as fo:
    for r, ids in want.items():
        for m in (1, 2):
            with open(f"{K}/data/fastq/{r}_{m}.fastq") as f:
                while True:
                    h = f.readline()
                    if not h: break
                    s = f.readline().strip(); f.readline(); f.readline()
                    rid = h[1:].split()[0].split("/")[0]
                    if rid in ids: fo.write(f">{r}|{rid}|{m}\n{s}\n")
# MRV type 3 references (taxids 538123 and 10886) from the current library
ref = f"{W}/mrv.fa"
with open(f"{K}/kraken2_dbs/current/library/added/BsVO5XdUTo.fna") as fi, open(ref, "w") as fo:
    keep = False
    for l in fi:
        if l[0] == ">":
            tid = l[1:].split()[0].split("|")[-1]
            keep = tid in ("538123", "10886")
        if keep: fo.write(l)
res = subprocess.run(f"blastn -task blastn-short -word_size 7 -query {q} -subject {ref} "
                     f"-outfmt '6 qseqid length pident' -evalue 1000 2>/dev/null",
                     shell=True, capture_output=True, text=True).stdout
best = defaultdict(int)
for l in res.splitlines():
    qid, ln, pid = l.split("\t")
    if float(pid) == 100.0: best[qid] = max(best[qid], int(ln))
rows = []
for r, ids in want.items():
    for rid in sorted(ids):
        m = max(best.get(f"{r}|{rid}|1", 0), best.get(f"{r}|{rid}|2", 0))
        rows.append([r, rid, m])
with open(os.path.join(OUT, "spurious_mrv_longest_match.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["run", "read_id", "longest_exact_match_to_MRV3_nt"]); w.writerows(rows)
L = [x[2] for x in rows]
print("reads", len(rows), "max longest match", max(L) if L else None, ">=17:", sum(1 for x in L if x >= 17))
