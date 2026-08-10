#!/usr/bin/env python3
"""Characterize what drives the MSL30->31 and MSL31->32 read-reassignment spikes: which
(before_taxid, after_taxid) pairs account for the bulk of discordant reads."""
import os
from collections import Counter

BASE = "/home/hm/kdcr"
RESULTS_DIR = os.path.join(BASE, "results_fixed")
RUNS = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195",
        "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]

def load_output(path):
    d = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line[0] != "C":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            read_id, assigned = parts[1], parts[2]
            try:
                d[read_id] = int(assigned)
            except ValueError:
                d[read_id] = assigned
    return d

def load_names(names_dmp_path, wanted):
    names = {}
    with open(names_dmp_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            if len(parts) < 4:
                continue
            taxid_s = parts[0].strip()
            if not taxid_s.isdigit():
                continue
            taxid = int(taxid_s)
            if taxid in wanted and "scientific name" in parts[3]:
                names[taxid] = parts[1].strip()
    return names

def characterize(a, b):
    pair_counts = Counter()
    for run in RUNS:
        da = load_output(os.path.join(RESULTS_DIR, a, f"{run}_output.txt"))
        db = load_output(os.path.join(RESULTS_DIR, b, f"{run}_output.txt"))
        for rid, ta in da.items():
            tb = db.get(rid)
            if tb is not None and ta != tb:
                pair_counts[(ta, tb)] += 1
    print(f"\n=== {a} -> {b}: top 15 (before_taxid, after_taxid) reassignment pairs ===")
    wanted = set()
    for (ta, tb), n in pair_counts.most_common(15):
        wanted.add(ta); wanted.add(tb)
    names = load_names(os.path.join(BASE, "kraken2_dbs_fixed", b, "taxonomy", "names.dmp"), wanted)
    total = sum(pair_counts.values())
    for (ta, tb), n in pair_counts.most_common(15):
        print(f"  {n:>6} ({100*n/total:5.1f}%)  {ta} [{names.get(ta,'?')}]  ->  {tb} [{names.get(tb,'?')}]")
    print(f"  TOTAL discordant pairs: {len(pair_counts)}, total discordant reads: {total}")

if __name__ == "__main__":
    characterize("msl30", "msl31")
    characterize("msl31", "msl32")
