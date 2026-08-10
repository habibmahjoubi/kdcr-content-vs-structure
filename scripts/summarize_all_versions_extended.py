#!/usr/bin/env python3
"""Extended version of summarize_all_versions.py covering all 12 MSL databases (msl30-msl41,
11 real transitions) instead of the original 5 (4 transitions). Correlates classification
change against delta_K_species (available uniformly back to MSL30) rather than delta_K_aggregate
(only available for MSL37-41, since ICTV's "Taxon Counts" sheet doesn't exist before MSL38)."""
import csv
import os

RESULTS_DIR = "/home/hm/kdcr/results"
DELTA_K_EXT_CSV = "/home/hm/kdcr/results/delta_K_species_extended_MSL30_41.csv"
OUT_WIDE = "/home/hm/kdcr/results/classification_by_msl_version_extended.csv"
OUT_DELTA = "/home/hm/kdcr/results/classification_change_vs_delta_K_extended.csv"

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
DBS = [f"msl{n}" for n in range(29, 42)]
RSV_TAXID = "11250"
REO_TAXID = "10882"


def cumulative_count(report_path, target_taxid):
    if not os.path.exists(report_path):
        return None
    with open(report_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            if cols[4].strip() == target_taxid:
                return int(cols[1])
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
    print(f"{'Run':<14}{'Virus':<6}{'Dilution':<10}" + "".join(f"{db:<7}" for db in DBS))
    for r in rows:
        print(f"{r['run']:<14}{r['virus']:<6}{r['dilution']:<10}" + "".join(f"{r[db]:<7}" for db in DBS))

    delta_k = {}
    with open(DELTA_K_EXT_CSV, "r", encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            delta_k[rec["transition"]] = float(rec["delta_K_species"])

    transitions = [(f"msl{a}", f"msl{a+1}", f"MSL{a}->MSL{a+1}") for a in range(29, 41)]

    delta_rows = []
    for db_a, db_b, trans_key in transitions:
        dk = delta_k.get(trans_key)
        total_abs_change_rsv = sum(abs(r[db_b] - r[db_a]) for r in rows if r["virus"] == "RSV")
        total_abs_change_reo = sum(abs(r[db_b] - r[db_a]) for r in rows if r["virus"] == "REO")
        delta_rows.append({
            "transition": trans_key,
            "delta_K_species": dk,
            "total_abs_change_RSV_reads": total_abs_change_rsv,
            "total_abs_change_REO_reads": total_abs_change_reo,
        })

    with open(OUT_DELTA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(delta_rows[0].keys()))
        w.writeheader()
        w.writerows(delta_rows)
    print(f"\nSaved delta comparison: {OUT_DELTA}")
    print(f"{'Transition':<16}{'delta_K_sp':<12}{'|Δ RSV|':<10}{'|Δ REO|':<10}")
    for r in delta_rows:
        print(f"{r['transition']:<16}{r['delta_K_species']:<12}{r['total_abs_change_RSV_reads']:<10}{r['total_abs_change_REO_reads']:<10}")

    try:
        from scipy import stats
        dk_vals = [r["delta_K_species"] for r in delta_rows]
        rsv_vals = [r["total_abs_change_RSV_reads"] for r in delta_rows]
        reo_vals = [r["total_abs_change_REO_reads"] for r in delta_rows]
        rho_rsv, p_rsv = stats.spearmanr(dk_vals, rsv_vals)
        rho_reo, p_reo = stats.spearmanr(dk_vals, reo_vals)
        rho_c, p_c = stats.spearmanr(dk_vals * 2, rsv_vals + reo_vals)
        print(f"\nSpearman (n={len(dk_vals)} transitions):")
        print(f"  RSV: rho={rho_rsv:.3f} p={p_rsv:.3f}")
        print(f"  REO: rho={rho_reo:.3f} p={p_reo:.3f}")
        print(f"  combined (n={len(dk_vals)*2}): rho={rho_c:.3f} p={p_c:.3f}")
    except ImportError:
        print("scipy not available for correlation stats")


if __name__ == "__main__":
    main()
