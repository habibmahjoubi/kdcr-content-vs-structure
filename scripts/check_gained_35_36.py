#!/usr/bin/env python3
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

gained_counts = Counter()
for run in RUNS:
    da = load_output(os.path.join(RESULTS_DIR, "msl35", f"{run}_output.txt"))
    db = load_output(os.path.join(RESULTS_DIR, "msl36", f"{run}_output.txt"))
    for rid, tb in db.items():
        if rid not in da:
            gained_counts[tb] += 1
wanted = set(t for t, n in gained_counts.most_common(10))
names = load_names(os.path.join(BASE, "kraken2_dbs_fixed", "msl36", "taxonomy", "names.dmp"), wanted)
print("Top newly-classified (gained) taxids at MSL35->MSL36:")
for t, n in gained_counts.most_common(10):
    print(f"  {n:>5}  {t} [{names.get(t,'?')}]")
