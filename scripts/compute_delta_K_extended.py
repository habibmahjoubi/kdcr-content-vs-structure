#!/usr/bin/env python3
"""Extend delta_K (species-level method) back to MSL30 (2015), giving 11 real transitions
instead of 4. ICTV's "Taxon Counts" sheet (the aggregate method) only exists from MSL38
onward, so the aggregate method stays limited to the original 4 transitions (MSL37-41) --
this is a genuine data-availability constraint of ICTV's own spreadsheets, not a choice.
The species-level method ("MSL of Last Change" column) is present in every version back to
MSL30 and is computed uniformly here for all 11 transitions.
"""
import os
import csv
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "data", "ictv_msl")   # ICTV MSL spreadsheets (https://ictv.global/msl)
OUT_DIR = os.path.join(ROOT, "results")

FILES = {
    29: "ICTV_MSL_2014_MSL29.v4.xls",
    30: "ICTV_MSL_2015_MSL30.v1.xlsx",
    31: "ICTV_MSL_2016_MSL31.v1.3.xlsx",
    32: "ICTV_MSL_2017_MSL32.v1.xlsx",
    33: "ICTV_MSL_2018a_MSL33.v1.xlsx",
    34: "ICTV_MSL_2018b_MSL34.v2.xlsx",
    35: "ICTV_MSL_2019_MSL35.v1.xlsx",
    36: "ICTV_MSL_2020_MSL36.v1.xlsx",
    37: "ICTV_MSL_2021_MSL37.v3.xlsx",
    38: "ICTV_MSL_2022_MSL38.v3.xlsx",
    39: "ICTV_MSL_2023_MSL39.v4.xlsx",
    40: "ICTV_MSL_2024_MSL40.v2.xlsx",
    41: "ICTV_MSL_2025_MSL41.v1.xlsx",
}

TRANSITIONS = [(a, a + 1) for a in range(29, 41)]

NON_DATA_SHEETS = {"Version", "Column Definitions", "Taxon Counts", "Taxa Renamed or Abolished",
                   "Taxa Renamed in MSL38", "Taxa Renamed in MSL39"}


def main_sheet(wb):
    candidates = [s for s in wb.sheetnames if s not in NON_DATA_SHEETS and "Renamed" not in s]
    if len(candidates) != 1:
        # fall back: prefer a sheet literally named "MSL" if present, else the longest name
        if "MSL" in candidates:
            return "MSL"
        candidates.sort(key=len, reverse=True)
    return candidates[0]


def delta_k_species_level(path, target_msl):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sn = main_sheet(wb)
    ws = wb[sn]
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    species_idx = header.index("Species")
    last_change_msl_idx = header.index("MSL of Last Change")

    total_species = 0
    changed_species = 0
    for r in rows:
        if r[species_idx] is None:
            continue
        total_species += 1
        moc = r[last_change_msl_idx]
        try:
            moc_num = int(str(moc).strip())
        except (TypeError, ValueError):
            moc_num = None
        if moc_num == target_msl:
            changed_species += 1
    wb.close()
    return (changed_species / total_species if total_species else 0.0), changed_species, total_species


def main():
    results = []
    for prev_v, new_v in TRANSITIONS:
        path = os.path.join(BASE, FILES[new_v])
        dk_sp, changed_sp, total_sp = delta_k_species_level(path, new_v)
        results.append({
            "transition": f"MSL{prev_v}->MSL{new_v}",
            "delta_K_species": round(dk_sp, 5),
            "species_changed": changed_sp,
            "species_total": total_sp,
        })
        print(f"MSL{prev_v}->MSL{new_v}: delta_K_species={dk_sp:.5f} ({changed_sp}/{total_sp})")

    out_csv = os.path.join(OUT_DIR, "delta_K_species_extended_MSL30_41.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)
    print(f"\nSaved: {out_csv}")


if __name__ == "__main__":
    main()
