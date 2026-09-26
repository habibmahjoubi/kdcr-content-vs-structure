#!/usr/bin/env python3
"""Read-level reassignment in the merge-corrected fixed-content library (5,518 taxids), which
contains the phiX174 genome NC_001422.1, for comparison with Design 2."""
import os, csv
K = "/home/hm/kdcr"; B = f"{K}/results_fixed_corrected"
RUNS9 = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195",
         "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]
MSL = [f"msl{n}" for n in range(29, 42)]
def load(p):
    d = {}
    for l in open(p, errors="replace"):
        if l[0] == "C":
            x = l.split("\t", 3); d[x[1]] = int(x[2])
    return d
rows = []
for i in range(12):
    a, b = MSL[i], MSL[i + 1]; both = disc = 0
    for r in RUNS9:
        pa, pb = f"{B}/{a}/{r}_output.txt", f"{B}/{b}/{r}_output.txt"
        if not (os.path.exists(pa) and os.path.exists(pb)): continue
        da, db = load(pa), load(pb)
        for k, t in da.items():
            u = db.get(k)
            if u is None: continue
            both += 1; disc += (t != u)
    rows.append([f"{a}->{b}", both, disc, round(100 * disc / both, 3) if both else ""]); print(rows[-1])
with open(os.path.join(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"), "reassignment_merge_corrected.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["update", "classified_in_both", "discordant", "pct"]); w.writerows(rows)
