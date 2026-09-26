#!/usr/bin/env python3
"""Characterise (1) the read-level jump between the official indexes of 17 May 2021 and 7 June 2022,
(2) the non-phiX reassignment on MSL32->MSL33 (Design 2): top taxid pairs, names from the current taxonomy,
and whether the reads are phiX174 (BLAST ids from stage s9 are available for the Design 2 runs only)."""
import os, csv
from collections import Counter
K = "/home/hm/kdcr"; W = f"{K}/work_designs"; R2 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
RUNS14 = ["SRR26352207", "SRR26352206", "SRR26352204", "SRR26352215", "SRR26352205", "SRR26352203", "SRR26352202",
          "SRR26352217", "SRR26352216", "SRR26352195", "SRR26352184", "SRR26352212", "SRR26352209", "SRR26352208"]
RUNS9 = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195", "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]


def cls(path):
    return {l.split("\t", 3)[1]: int(l.split("\t", 3)[2]) for l in open(path) if l[0] == "C"}


def names(ids):
    n = {}
    for l in open(f"{K}/kraken2_dbs/current/taxonomy/names.dmp", errors="replace"):
        x = l.split("\t|\t")
        t = int(x[0])
        if t in ids and "scientific name" in x[3]: n[t] = x[1]
    return n


rows = []
# (1) official 20210517 -> 20220607
pairs, lost, gained, n_a, n_b = Counter(), Counter(), Counter(), 0, 0
for r in RUNS14:
    a = cls(f"{W}/official/cls/20210517/{r}_output.txt"); b = cls(f"{W}/official/cls/20220607/{r}_output.txt")
    n_a += len(a); n_b += len(b)
    for k, t in a.items():
        if k in b and b[k] != t: pairs[(t, b[k])] += 1
        if k not in b: lost[t] += 1
    for k, t in b.items():
        if k not in a: gained[t] += 1
ids = set(x for p in pairs for x in p) | set(lost) | set(gained)
nm = names(ids)
rows.append(["official 20210517->20220607", "classified reads", n_a, n_b, ""])
for (t1, t2), c in pairs.most_common(8):
    rows.append(["official 20210517->20220607", "reassigned", c, f"{t1} {nm.get(t1, '?')}", f"{t2} {nm.get(t2, '?')}"])
for t, c in lost.most_common(5):
    rows.append(["official 20210517->20220607", "classified only in 2021", c, f"{t} {nm.get(t, '?')}", ""])
for t, c in gained.most_common(5):
    rows.append(["official 20210517->20220607", "classified only in 2022", c, f"{t} {nm.get(t, '?')}", ""])
# (2) Design 2 MSL32 -> MSL33, non-phiX
pairs2 = Counter()
for r in RUNS9:
    phx = set(l.strip() for l in open(f"{W}/phix/{r}.phix_ids"))
    a = cls(f"{K}/results_fixed/msl32/{r}_output.txt"); b = cls(f"{K}/results_fixed/msl33/{r}_output.txt")
    for k, t in a.items():
        if k in b and b[k] != t and k not in phx: pairs2[(t, b[k])] += 1
nm2 = names(set(x for p in pairs2 for x in p))
for (t1, t2), c in pairs2.most_common(8):
    rows.append(["Design 2 msl32->msl33 (non-phiX)", "reassigned", c, f"{t1} {nm2.get(t1, '?')}", f"{t2} {nm2.get(t2, '?')}"])
with open(os.path.join(R2, "characterise_jumps.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["comparison", "category", "reads", "taxon_before", "taxon_after"]); w.writerows(rows)
for r in rows: print("\t".join(map(str, r)))
