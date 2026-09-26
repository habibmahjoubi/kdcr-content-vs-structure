#!/usr/bin/env python3
"""Genuine RSV reads of the RSV 1e6 run lost at confidence 0.1 in the prebuilt databases: where do they go?"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exposure_host_prebuilt import parents_report, clade_fn, W2, W3, OFFICIAL, RSV, tsv

rows = []
for o in OFFICIAL:
    par0, nm0 = parents_report(f"{W2}/official/cls/{o}/SRR26352206_report.txt"); inc0 = clade_fn(par0)
    par1, nm1 = parents_report(f"{W3}/oconf/{o}/SRR26352206_report.txt")
    c0 = {}
    for l in open(f"{W2}/official/cls/{o}/SRR26352206_output.txt"):
        if l[0] == "C":
            c = l.rstrip("\n").split("\t")
            if inc0(int(c[2]), RSV) or int(c[2]) == 12814: c0[c[1]] = (int(c[2]), c[4])
    c1 = {}
    for l in open(f"{W3}/oconf/{o}/SRR26352206_output.txt"):
        c = l.rstrip("\n").split("\t")
        if c[1] in c0: c1[c[1]] = int(c[2])
    for r, (t, hits) in c0.items():
        t1 = c1.get(r, 0)
        rows.append([o, r, t, nm0.get(t, ""), t1, nm1.get(t1, "unclassified" if t1 == 0 else ""), hits[:200]])
tsv("rsv1e6_reads_conf0_vs_conf0.1_prebuilt.tsv", ["index", "read", "taxid_conf0", "name_conf0", "taxid_conf0.1", "name_conf0.1", "hits_conf0"], rows)
