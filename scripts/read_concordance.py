#!/usr/bin/env python3
"""A3: read-to-read concordance analysis between consecutive MSL versions.

The cumulative clade-count proxy used throughout Sections 5.1-5.4 is insensitive to
within-clade reassignments and to any change in classification correctness that does not
move the target's own cumulative count. This script instead computes, for EVERY read
classified by Kraken2 in a run (not only reads assigned to RSV's or REO's own taxid), the
proportion whose assigned taxid changes between consecutive taxonomy versions -- a much
higher-power observable (~32,500 classified reads/run vs. 1-18 target-specific reads/run).

Run on the fixed-library (confound-controlled) results, since that is the valid comparison
once library-composition growth is controlled for (see delta_K_ncbi.py / results_fixed).
"""
import csv
import os

BASE = "/home/hm/kdcr"
RESULTS_DIR = os.path.join(BASE, "results_fixed")
MSL_ORDER = [f"msl{n}" for n in range(29, 42)]
RUNS = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195",
        "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]

def load_output(path):
    """Return dict read_id -> assigned_taxid, CLASSIFIED reads only (status C).
    Unclassified reads (~98% of each file) are simply absent from the dict -- this keeps
    memory bounded to the ~32K classified reads/run instead of all 2M reads/run."""
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

def main():
    print("Loading per-read outputs for all runs x all 13 versions (fixed library)...")
    data = {}  # (run, msl) -> {read_id: taxid_or_None}
    for msl in MSL_ORDER:
        for run in RUNS:
            path = os.path.join(RESULTS_DIR, msl, f"{run}_output.txt")
            data[(run, msl)] = load_output(path)
        print(f"  loaded {msl}")

    rows = []
    for i in range(len(MSL_ORDER) - 1):
        a, b = MSL_ORDER[i], MSL_ORDER[i + 1]
        trans = f"MSL{a[3:]}->MSL{b[3:]}"
        n_classified_both = n_concordant = n_discordant = n_gained = n_lost = 0
        per_run_discordant = {}
        for run in RUNS:
            da, db = data[(run, a)], data[(run, b)]
            run_disc = 0
            all_ids = set(da.keys()) | set(db.keys())
            for rid in all_ids:
                ta, tb = da.get(rid), db.get(rid)
                ca, cb = ta is not None, tb is not None
                if ca and cb:
                    n_classified_both += 1
                    if ta == tb:
                        n_concordant += 1
                    else:
                        n_discordant += 1
                        run_disc += 1
                elif cb and not ca:
                    n_gained += 1
                elif ca and not cb:
                    n_lost += 1
            per_run_discordant[run] = run_disc
        rate = n_discordant / n_classified_both if n_classified_both else float("nan")
        rows.append({
            "transition": trans, "n_classified_both": n_classified_both,
            "n_concordant": n_concordant, "n_discordant": n_discordant,
            "discordance_rate": round(rate, 6),
            "n_gained_classification": n_gained, "n_lost_classification": n_lost,
        })
        print(f"{trans}: classified_both={n_classified_both} discordant={n_discordant} "
              f"rate={rate:.5f} gained={n_gained} lost={n_lost}")

    out_csv = os.path.join(BASE, "results_fixed", "read_concordance_by_transition.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved: {out_csv}")

    # Correlate discordance_rate against delta_K (ICTV species-level and NCBI-based)
    from scipy import stats
    ictv = {}
    with open(os.path.join(BASE, "results", "delta_K_species_extended_MSL30_41.csv"), encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            ictv[rec["transition"]] = float(rec["delta_K_species"])
    # MSL29->MSL30 not in the extended species file (starts at MSL30); pull from original 12-pt file if present
    orig29 = os.path.join(BASE, "results", "delta_K_ICTV_MSL_transitions.csv")

    ncbi = {}
    with open(os.path.join(BASE, "results", "delta_K_ncbi_transitions.csv"), encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            ncbi[rec["transition"]] = (float(rec["delta_K_ncbi_full"]), float(rec["delta_K_ncbi_viral"]))

    merged = []
    for r in rows:
        t = r["transition"]
        if t in ictv and t in ncbi:
            merged.append({
                "transition": t,
                "delta_K_ictv": ictv[t],
                "delta_K_ncbi_full": ncbi[t][0],
                "delta_K_ncbi_viral": ncbi[t][1],
                "discordance_rate": r["discordance_rate"],
            })

    print(f"\nMerged n={len(merged)} transitions (some may be missing from ICTV species file)")
    if len(merged) >= 4:
        x_ictv = [m["delta_K_ictv"] for m in merged]
        x_full = [m["delta_K_ncbi_full"] for m in merged]
        x_viral = [m["delta_K_ncbi_viral"] for m in merged]
        y = [m["discordance_rate"] for m in merged]
        for label, x in [("delta_K_ictv", x_ictv), ("delta_K_ncbi_full", x_full), ("delta_K_ncbi_viral", x_viral)]:
            rho, p = stats.spearmanr(x, y)
            print(f"discordance_rate vs {label}: rho={rho:.4f} p={p:.4f} n={len(x)}")

    out_merged = os.path.join(BASE, "results_fixed", "read_concordance_vs_delta_K.csv")
    with open(out_merged, "w", newline="", encoding="utf-8") as f:
        if merged:
            w = csv.DictWriter(f, fieldnames=list(merged[0].keys()))
            w.writeheader()
            w.writerows(merged)
    print(f"Saved: {out_merged}")

if __name__ == "__main__":
    main()
