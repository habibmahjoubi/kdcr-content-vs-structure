#!/usr/bin/env python3
"""B7: analyze the full --confidence 0.5 sweep (13 dbs x 9 runs) -- check whether the
RSV/REO rank-matched stability pattern holds, not just a current-db spot-check."""
import os

BASE = "/home/hm/kdcr"
MSL_ORDER = [f"msl{n}" for n in range(29, 42)]
RSV_TAXID = 11250
REO_TAXID = 351073
RUNS_RSV = {"SRR26352206": "1e6", "SRR26352205": "1e3", "SRR26352202": "1e1", "SRR26352207": "neg"}
RUNS_REO = {"SRR26352216": "1e6", "SRR26352195": "1e5", "SRR26352212": "1e3", "SRR26352208": "1e1", "SRR26352217": "neg"}

def clade_count(report_path, taxid):
    if not os.path.exists(report_path):
        return None
    with open(report_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                this_taxid = int(parts[4])
            except ValueError:
                continue
            if this_taxid == taxid:
                return int(parts[1])
    return 0

print("=== RSV (taxid 11250) at --confidence 0.5, all 13 dbs ===")
for run, dil in RUNS_RSV.items():
    series = [clade_count(os.path.join(BASE, "results", "conf_0.5", msl, f"{run}_report.txt"), RSV_TAXID) for msl in MSL_ORDER]
    print(f"{run} ({dil}): {series}")

print("\n=== REO (taxid 351073, rank-matched) at --confidence 0.5, all 13 dbs ===")
for run, dil in RUNS_REO.items():
    series = [clade_count(os.path.join(BASE, "results", "conf_0.5", msl, f"{run}_report.txt"), REO_TAXID) for msl in MSL_ORDER]
    print(f"{run} ({dil}): {series}")
