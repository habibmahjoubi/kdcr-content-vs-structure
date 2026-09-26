#!/usr/bin/env python3
"""Design 3: fixed, complete reference content (all 19,582 RefSeq viral sequences), each sequence
labelled in each version with the taxid valid at that date:
  1. its current taxid, if present in the snapshot;
  2. otherwise a predecessor merged into it later (current merged.dmp, resolved transitively);
  3. otherwise the nearest ancestor of its current taxid present in the snapshot
     (e.g. the species above a virus-level node created later, or the genus above a new species).
Writes one FASTA per version and a label table. Parent links only, one snapshot at a time.
"""
import os, sys, csv
from collections import defaultdict, Counter

K = "/home/hm/kdcr"
W = f"{K}/work_designs/design3"
SRC = f"{K}/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna"   # unmasked source library
CUR_TAX = f"{K}/kraken2_dbs/current/taxonomy"
MSL = [f"msl{n}" for n in range(29, 42)]
os.makedirs(W, exist_ok=True)


def parents(path):
    p = {}
    with open(path) as f:
        for l in f:
            x = l.split("\t|\t", 2); p[int(x[0])] = int(x[1])
    return p


cur_par = parents(f"{CUR_TAX}/nodes.dmp")
merged = {}
with open(f"{CUR_TAX}/merged.dmp") as f:
    for l in f:
        a, b = l.replace("\t|\n", "").split("\t|\t")[:2]; merged[int(a)] = int(b)
def resolve(t, seen=None):
    seen = seen or set()
    while t in merged and t not in seen:
        seen.add(t); t = merged[t]
    return t
pred = defaultdict(set)          # current taxid -> all old taxids merged into it
for old in merged:
    pred[resolve(old)].add(old)

# accession -> current taxid, from the source headers
acc_tax = {}
with open(SRC) as f:
    for l in f:
        if l[0] == ">":
            h = l[1:].split()[0]; parts = h.split("|")
            acc_tax[parts[0]] = int(parts[2])
print("sequences", len(acc_tax), flush=True)

rows = []
for m in MSL:
    nodes = set(parents(f"{K}/kraken2_dbs/{m}/taxonomy/nodes.dmp"))
    lab, cat = {}, Counter()
    for acc, c in acc_tax.items():
        if c in nodes:
            lab[acc] = c; cat["current"] += 1; continue
        ps = sorted(p for p in pred.get(c, ()) if p in nodes)
        if ps:
            lab[acc] = ps[0]; cat["predecessor"] += 1; continue
        t = c
        while t in cur_par and t not in nodes and cur_par[t] != t:
            t = cur_par[t]
        lab[acc] = t if t in nodes else 1; cat["ancestor"] += 1
    out = f"{W}/{m}.fna"
    with open(SRC) as fi, open(out, "w") as fo:
        for l in fi:
            if l[0] == ">":
                acc = l[1:].split()[0].split("|")[0]
                rest = l.split(" ", 1)[1] if " " in l else "\n"
                fo.write(f">{acc}|kraken:taxid|{lab[acc]} {rest}")
            else:
                fo.write(l)
    with open(f"{W}/{m}_labels.tsv", "w") as f:
        for acc in sorted(lab): f.write(f"{acc}\t{acc_tax[acc]}\t{lab[acc]}\n")
    rows.append([m, len(lab), cat["current"], cat["predecessor"], cat["ancestor"]])
    print(rows[-1], flush=True)
    del nodes

OUT = sys.argv[1] if len(sys.argv) > 1 else W
with open(os.path.join(OUT, "design3_labelling_summary.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["version", "sequences", "current_taxid", "merge_predecessor", "nearest_existing_ancestor"])
    w.writerows(rows)
