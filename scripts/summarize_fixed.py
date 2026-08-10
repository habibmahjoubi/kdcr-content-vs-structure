#!/usr/bin/env python3
import csv
import os

RESULTS_DIR = "/home/hm/kdcr/results_fixed"
DELTA_K_CSV = "/home/hm/kdcr/results/delta_K_species_extended_MSL30_41.csv"
OUT_CSV = "/home/hm/kdcr/results_fixed/classification_rank_matched_fixed_library.csv"

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
REO_TAXID_NEW = "351073"


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
        target = RSV_TAXID if virus == "RSV" else REO_TAXID_NEW
        for db in DBS:
            path = os.path.join(RESULTS_DIR, db, f"{run}_report.txt")
            row[db] = cumulative_count(path, target)
        rows.append(row)
        print(row)

    delta_k = {}
    with open(DELTA_K_CSV, "r", encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            delta_k[rec["transition"]] = float(rec["delta_K_species"])

    transitions = [(f"msl{a}", f"msl{a+1}", f"MSL{a}->MSL{a+1}") for a in range(29, 41)]
    delta_rows = []
    for db_a, db_b, trans_key in transitions:
        dk = delta_k.get(trans_key)
        d_rsv = sum(abs(r[db_b] - r[db_a]) for r in rows if r["virus"] == "RSV")
        d_reo = sum(abs(r[db_b] - r[db_a]) for r in rows if r["virus"] == "REO")
        delta_rows.append({"transition": trans_key, "delta_K_species": dk,
                            "abs_change_RSV": d_rsv, "abs_change_REO_rank_matched": d_reo})

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(delta_rows[0].keys()))
        w.writeheader()
        w.writerows(delta_rows)

    print(f"\n{'Transition':<16}{'delta_K':<10}{'|dRSV|':<8}{'|dREO_rankmatched|':<10}")
    for r in delta_rows:
        print(f"{r['transition']:<16}{r['delta_K_species']:<10}{r['abs_change_RSV']:<8}{r['abs_change_REO_rank_matched']:<10}")

    from scipy import stats
    dk_vals = [r["delta_K_species"] for r in delta_rows]
    rsv_vals = [r["abs_change_RSV"] for r in delta_rows]
    reo_vals = [r["abs_change_REO_rank_matched"] for r in delta_rows]
    rho_rsv, p_rsv = stats.spearmanr(dk_vals, rsv_vals)
    rho_reo, p_reo = stats.spearmanr(dk_vals, reo_vals)
    print(f"\nFIXED-LIBRARY RESULTS:")
    print(f"RSV: rho={rho_rsv:.3f} p={p_rsv:.3f}")
    print(f"REO (rank-matched): rho={rho_reo:.3f} p={p_reo:.3f}")
    print(f"RSV total abs change across 12 transitions: {sum(rsv_vals)}")
    print(f"REO total abs change across 12 transitions: {sum(reo_vals)}")
    print(f"\nSaved: {OUT_CSV}")


if __name__ == "__main__":
    main()
