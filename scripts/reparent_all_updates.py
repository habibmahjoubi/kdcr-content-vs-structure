#!/usr/bin/env python3
"""Reparenting at fixed content (Design 2), all twelve updates. For every taxon of the fixed library
whose parent node changed between consecutive snapshots, 1,000 read pairs are simulated per accession
(wgsim, 150 bp, 1% error, -r 0 -R 0, insert 250 +/- 25, seed derived from the accession) and
classified against the 13 fixed-content databases (pipeline stage s7). The summary (summarize_designs_panel.py)
decomposes each change read by read: same taxid vs different taxid, and clade membership."""
import os, sys, subprocess, zlib

K = "/home/hm/kdcr"
R2 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
W = f"{K}/work_designs/reparent"
MSL = [f"msl{n}" for n in range(29, 42)]


def parents(m):
    p = {}
    with open(f"{K}/kraken2_dbs_fixed/{m}/taxonomy/nodes.dmp") as f:
        for l in f:
            x = l.split("\t|\t", 2); p[int(x[0])] = int(x[1])
    return p


def simulate():
    os.makedirs(W, exist_ok=True)
    acc_tax = {}
    for l in open(f"{K}/kraken2_dbs_fixed/msl29/seqid2taxid.map"):
        s, t = l.split("\t"); acc_tax[s.split("|")[0]] = int(t)
    taxa = set(acc_tax.values())
    prev = parents(MSL[0]); per_update = {}
    for a, b in zip(MSL, MSL[1:]):
        cur = parents(b)
        rep = {t for t in taxa if t in prev and t in cur and prev[t] != cur[t] and t not in (11250, 351073)}
        per_update[f"{a}->{b}"] = sorted(acc for acc, t in acc_tax.items() if t in rep)
        prev = cur
    union = sorted(set(x for v in per_update.values() for x in v))
    print("accessions to simulate:", len(union), flush=True)
    seqs, cur_id = {}, None
    for l in open(f"{K}/kraken2_dbs_fixed/library_fixed.fna"):
        if l[0] == ">":
            cur_id = l[1:].split("|")[0]; seqs[cur_id] = [] if cur_id in union else None
        elif seqs.get(cur_id) is not None:
            seqs[cur_id].append(l.strip())
    o1, o2 = open(f"{W}/rep_1.fq", "w"), open(f"{W}/rep_2.fq", "w")
    skipped = []
    for acc in union:
        s = "".join(seqs[acc])
        if len(s) < 400: skipped.append(acc); continue
        with open(f"{W}/tmp.fa", "w") as f: f.write(f">{acc}\n{s}\n")
        seed = zlib.crc32(acc.encode()) % 2147483647
        subprocess.run(f"wgsim -N 1000 -1 150 -2 150 -e 0.01 -r 0 -R 0 -d 250 -s 25 -S {seed} "
                       f"{W}/tmp.fa {W}/t1.fq {W}/t2.fq > /dev/null 2>&1", shell=True, check=True)
        for src, dst in ((f"{W}/t1.fq", o1), (f"{W}/t2.fq", o2)):
            with open(src) as f:
                for j, l in enumerate(f):
                    if j % 4 == 0: l = f"@{acc}__{l[1:]}"
                    dst.write(l)
    o1.close(); o2.close()
    with open(f"{R2}/reparent_manifest.tsv", "w") as f:
        f.write("update\tn_reparented_taxa_accessions\taccessions\n")
        for u, v in per_update.items(): f.write(f"{u}\t{len(v)}\t{';'.join(v)}\n")
    with open(f"{R2}/reparent_skipped_short.tsv", "w") as f:
        f.write("\n".join(skipped) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "simulate": simulate()
