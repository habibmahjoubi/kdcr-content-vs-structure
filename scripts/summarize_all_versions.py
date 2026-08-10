#!/usr/bin/env python3
"""Summarize RSV/REO read counts for all 9 KDCR runs across all 5 per-MSL-version Kraken2
databases, then compute the change in read counts across each real MSL transition and
compare its magnitude to the previously computed delta_K (P1 test data).

IMPORTANT (2026-08-07): counts are matched by STABLE TAXID (cumulative clade count, Kraken2
report column 2), not by taxon name substring. An earlier name-substring version ("syncytial")
produced a spurious jump between MSL39 and MSL40 that turned out to be a pure nomenclature
change -- NCBI renamed taxid 11250 from "Human orthopneumovirus" to "human respiratory
syncytial virus" without any change in the underlying read classification. Matching by taxid
is robust to exactly this kind of rename, which is itself a real and worth-reporting knowledge-
drift phenomenon (see guide_reproductibilite.md).
"""
import csv
import os

RESULTS_DIR = "/home/hm/kdcr/results"
DELTA_K_CSV = "/home/hm/kdcr/results/delta_K_ICTV_MSL_transitions.csv"
OUT_WIDE = "/home/hm/kdcr/results/classification_by_msl_version.csv"
OUT_DELTA = "/home/hm/kdcr/results/classification_change_vs_delta_K.csv"

MANIFEST = {
    "SRR26352217": ("REO", "neg"),
    "SRR26352207": ("RSV", "neg"),
    "SRR26352216": ("REO", "1e6"),
    "SRR26352206": ("RSV", "1e6"),
    "SRR26352195": ("REO", "1e5"),
    "SRR26352212": ("REO", "1e3"),
    "SRR26352205": ("RSV", "1e3"),
    "SRR26352208": ("REO", "1e1"),
    "SRR26352202": ("RSV", "1e1"),
}
DBS = ["msl37", "msl38", "msl39", "msl40", "msl41"]

# Stable anchor taxids, identified once from a positive sample's report and confirmed present
# (existing, even if renamed) across all 5 snapshots -- see guide_reproductibilite.md.
RSV_TAXID = "11250"   # Human orthopneumovirus / human respiratory syncytial virus (genus-level species node)
REO_TAXID = "10882"   # Orthoreovirus (genus)

def cumulative_count(report_path, target_taxid):
    if not os.path.exists(report_path):
        return None
    with open(report_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            if cols[4].strip() == target_taxid:
                return int(cols[1])  # column 2: cumulative reads in this taxon's clade
    return 0

def main():
    rows = []
    for run, (virus, dilution) in MANIFEST.items():
        row = {"run": run, "virus": virus, "dilution": dilution}
        target = RSV_TAXID if virus == "RSV" else REO_TAXID
        for db in DBS:
            path = os.path.join(RESULTS_DIR, db, f"{run}_report.txt")
            row[db] = cumulative_count(path, target)
        rows.append(row)

    fieldnames = ["run", "virus", "dilution"] + DBS
    with open(OUT_WIDE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"Saved wide table: {OUT_WIDE}\n")
    print(f"{'Run':<14}{'Virus':<6}{'Dilution':<10}" + "".join(f"{db:<8}" for db in DBS))
    for r in rows:
        print(f"{r['run']:<14}{r['virus']:<6}{r['dilution']:<10}" + "".join(f"{r[db]:<8}" for db in DBS))

    delta_k = {}
    with open(DELTA_K_CSV, "r", encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            delta_k[rec["transition"]] = float(rec["delta_K_aggregate"])

    transitions = [("msl37", "msl38", "MSL37->MSL38"),
                   ("msl38", "msl39", "MSL38->MSL39"),
                   ("msl39", "msl40", "MSL39->MSL40"),
                   ("msl40", "msl41", "MSL40->MSL41")]

    delta_rows = []
    for db_a, db_b, trans_key in transitions:
        dk = delta_k.get(trans_key)
        total_abs_change_rsv = sum(abs(r[db_b] - r[db_a]) for r in rows if r["virus"] == "RSV")
        total_abs_change_reo = sum(abs(r[db_b] - r[db_a]) for r in rows if r["virus"] == "REO")
        delta_rows.append({
            "transition": trans_key,
            "delta_K_aggregate": dk,
            "total_abs_change_RSV_reads": total_abs_change_rsv,
            "total_abs_change_REO_reads": total_abs_change_reo,
        })

    with open(OUT_DELTA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(delta_rows[0].keys()))
        w.writeheader()
        w.writerows(delta_rows)
    print(f"\nSaved delta comparison: {OUT_DELTA}")
    print(f"{'Transition':<15}{'delta_K':<10}{'|Δ RSV reads|':<16}{'|Δ REO reads|':<16}")
    for r in delta_rows:
        print(f"{r['transition']:<15}{r['delta_K_aggregate']:<10}{r['total_abs_change_RSV_reads']:<16}{r['total_abs_change_REO_reads']:<16}")

if __name__ == "__main__":
    main()
