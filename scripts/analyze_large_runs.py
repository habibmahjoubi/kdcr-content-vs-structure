#!/usr/bin/env python3
"""A2: analyze the 10M-pair (5x) reclassification of the four informative runs against
all 13 fixed-library databases. Extract cumulative clade counts at the rank-matched
target taxids (RSV=11250, REO=351073) and recompute the delta_K correlation with more
statistical power than the original 2M-pair design."""
import csv
import os

BASE = "/home/hm/kdcr"
RESULTS_DIR = os.path.join(BASE, "results_fixed_large")
MSL_ORDER = [f"msl{n}" for n in range(29, 42)]
RUNS = {
    "SRR26352206": ("RSV", "1e6", 11250),
    "SRR26352205": ("RSV", "1e3", 11250),
    "SRR26352202": ("RSV", "1e1", 11250),
    "SRR26352216": ("REO", "1e6", 351073),
}

def clade_count(report_path, taxid):
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
                return int(parts[1])  # cumulative clade-count column
    return 0

rows = []
for msl in MSL_ORDER:
    row = {"version": msl}
    for run, (virus, dil, taxid) in RUNS.items():
        rep = os.path.join(RESULTS_DIR, msl, f"{run}_report.txt")
        row[f"{virus}_{dil}"] = clade_count(rep, taxid)
    rows.append(row)

out_csv = os.path.join(BASE, "results_fixed_large", "large_run_counts_by_version.csv")
fieldnames = ["version"] + [f"{v}_{d}" for v, d in [("RSV","1e6"),("RSV","1e3"),("RSV","1e1"),("REO","1e6")]]
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print(f"{'version':<8}" + "".join(f"{k:>10}" for k in fieldnames[1:]))
for row in rows:
    print(f"{row['version']:<8}" + "".join(f"{row[k]:>10}" for k in fieldnames[1:]))
print(f"\nSaved: {out_csv}")

# --- Statistics ---
from scipy import stats
import statistics as st

rsv_1e6 = [row["RSV_1e6"] for row in rows]
rsv_1e3 = [row["RSV_1e3"] for row in rows]
rsv_1e1 = [row["RSV_1e1"] for row in rows]
reo_1e6 = [row["REO_1e6"] for row in rows]

print("\n--- Summary stats (n=13 versions each) ---")
for label, series in [("RSV_1e6", rsv_1e6), ("RSV_1e3", rsv_1e3), ("RSV_1e1", rsv_1e1), ("REO_1e6", reo_1e6)]:
    print(f"{label}: min={min(series)} max={max(series)} mean={st.mean(series):.2f} "
          f"sd={st.stdev(series):.3f} range={max(series)-min(series)}")

# per-transition |delta| summed across the 4 runs (mirrors original Table 4 summary)
delta_sum = []
for i in range(len(MSL_ORDER) - 1):
    d = 0
    for series in (rsv_1e6, rsv_1e3, rsv_1e1, reo_1e6):
        d += abs(series[i+1] - series[i])
    delta_sum.append(d)
transitions = [f"MSL{MSL_ORDER[i][3:]}->MSL{MSL_ORDER[i+1][3:]}" for i in range(len(MSL_ORDER)-1)]

print("\n--- |delta| per transition (summed across 4 runs, 5x-depth data) ---")
for t, d in zip(transitions, delta_sum):
    print(f"  {t}: {d}")

# correlate against ICTV species-level delta_K and NCBI-based delta_K
ictv = {}
with open(os.path.join(BASE, "results", "delta_K_species_extended_MSL30_41.csv"), encoding="utf-8") as f:
    for rec in csv.DictReader(f):
        ictv[rec["transition"]] = float(rec["delta_K_species"])
# MSL29->30 from original 4/12-pt file if present
orig = os.path.join(BASE, "results", "delta_K_ICTV_MSL_transitions.csv")
if os.path.exists(orig):
    with open(orig, encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            key = rec.get("transition")
            if key and key not in ictv and "delta_K_species" in rec:
                try:
                    ictv[key] = float(rec["delta_K_species"])
                except (ValueError, KeyError):
                    pass

ncbi = {}
with open(os.path.join(BASE, "results", "delta_K_ncbi_transitions.csv"), encoding="utf-8") as f:
    for rec in csv.DictReader(f):
        ncbi[rec["transition"]] = float(rec["delta_K_ncbi_viral"])

x_ictv, x_ncbi, y = [], [], []
for t, d in zip(transitions, delta_sum):
    if t in ictv and t in ncbi:
        x_ictv.append(ictv[t]); x_ncbi.append(ncbi[t]); y.append(d)

print(f"\nMerged n={len(y)} transitions")
if len(y) >= 4:
    rho1, p1 = stats.spearmanr(x_ictv, y)
    rho2, p2 = stats.spearmanr(x_ncbi, y)
    print(f"delta_K_ictv_species vs |delta|(5x depth): rho={rho1:.4f} p={p1:.4f} n={len(y)}")
    print(f"delta_K_ncbi_viral vs |delta|(5x depth): rho={rho2:.4f} p={p2:.4f} n={len(y)}")

# Separately: RSV 1e6 alone vs REO 1e6 alone at 5x depth, matching original per-virus split
rsv_delta = [abs(rsv_1e6[i+1]-rsv_1e6[i]) + abs(rsv_1e3[i+1]-rsv_1e3[i]) + abs(rsv_1e1[i+1]-rsv_1e1[i]) for i in range(len(MSL_ORDER)-1)]
reo_delta = [abs(reo_1e6[i+1]-reo_1e6[i]) for i in range(len(MSL_ORDER)-1)]
xr, yr, xo, yo = [], [], [], []
for t, dr, do in zip(transitions, rsv_delta, reo_delta):
    if t in ictv:
        xr.append(ictv[t]); yr.append(dr)
        xo.append(ictv[t]); yo.append(do)
if len(yr) >= 4:
    rho_r, p_r = stats.spearmanr(xr, yr)
    rho_o, p_o = stats.spearmanr(xo, yo)
    print(f"\nRSV-only |delta| vs delta_K_ictv_species (5x depth): rho={rho_r:.4f} p={p_r:.4f} n={len(yr)}")
    print(f"REO-only |delta| vs delta_K_ictv_species (5x depth): rho={rho_o:.4f} p={p_o:.4f} n={len(yo)}")
    print(f"RSV_1e6 series (5x depth): {rsv_1e6}")
    print(f"REO_1e6 series (5x depth): {reo_1e6}")
