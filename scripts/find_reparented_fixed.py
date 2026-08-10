#!/usr/bin/env python3
"""A1: among the 5270 fixed-library taxids (present in ALL 13 MSL snapshots, sequence
never entering/leaving the library), find taxa whose PARENT node changed (reparenting =
a genuine existential/structural taxonomy edit) on some transition, excluding RSV(11250)
and REO(351073) already tested. These are candidates for a panel test of whether
taxonomy-structure-only edits (at constant library/reference composition) ever drive
classification change -- generalizing the single REO anecdote (existential event at
MSL40->41, zero effect) to a proper sample.
"""
import csv
import os

BASE = "/home/hm/kdcr"
MSL_ORDER = [f"msl{n}" for n in range(29, 42)]
EXCLUDE = {11250, 351073}

def load_nodes(nodes_dmp_path):
    d = {}
    with open(nodes_dmp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            if len(parts) < 2:
                continue
            taxid_s = parts[0].strip()
            parent_s = parts[1].strip()
            if taxid_s.isdigit() and parent_s.isdigit():
                d[int(taxid_s)] = int(parent_s)
    return d

def load_names(names_dmp_path, wanted):
    names = {}
    with open(names_dmp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            if len(parts) < 4:
                continue
            taxid_s = parts[0].strip()
            if not taxid_s.isdigit():
                continue
            taxid = int(taxid_s)
            if taxid in wanted and "scientific name" in parts[3]:
                names[taxid] = parts[1].strip()
    return names

def main():
    fixed_taxids = set()
    with open(os.path.join(BASE, "kraken2_dbs_fixed", "fixed_taxids.txt")) as f:
        for line in f:
            line = line.strip()
            if line.isdigit():
                fixed_taxids.add(int(line))
    print(f"Fixed taxids loaded: {len(fixed_taxids)}")

    node_maps = {}
    for msl in MSL_ORDER:
        path = os.path.join(BASE, "kraken2_dbs_fixed", msl, "taxonomy", "nodes.dmp")
        node_maps[msl] = load_nodes(path)

    rows = []
    per_transition_candidates = {}
    for i in range(len(MSL_ORDER) - 1):
        a, b = MSL_ORDER[i], MSL_ORDER[i + 1]
        old, new = node_maps[a], node_maps[b]
        trans = f"MSL{a[3:]}->MSL{b[3:]}"
        reparented = []
        for t in fixed_taxids:
            if t in EXCLUDE:
                continue
            po = old.get(t)
            pn = new.get(t)
            if po is None or pn is None:
                continue  # shouldn't happen since t is in all 13 snapshots by construction
            if po != pn:
                reparented.append(t)
        per_transition_candidates[trans] = reparented
        rows.append({"transition": trans, "n_reparented_fixed_taxa": len(reparented)})
        print(f"{trans}: {len(reparented)} fixed-library taxa reparented")

    out_csv = os.path.join(BASE, "results_fixed", "reparented_candidates_per_transition.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["transition", "n_reparented_fixed_taxa"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved: {out_csv}")

    # Pick the transition with the most candidates for the panel test
    best_trans = max(per_transition_candidates, key=lambda k: len(per_transition_candidates[k]))
    best_list = per_transition_candidates[best_trans]
    print(f"\nBest transition for panel test: {best_trans} ({len(best_list)} candidates)")

    # Save candidate taxid list with names for the best transition (cap at 60 for a panel
    # comparable in size to the 56-taxon positive-control panel)
    wanted = set(best_list)
    names = load_names(os.path.join(BASE, "kraken2_dbs_fixed", MSL_ORDER[0], "taxonomy", "names.dmp"), wanted)
    out_panel = os.path.join(BASE, "results_fixed", f"reparented_panel_{best_trans.replace('>','').replace('-','_')}.csv")
    with open(out_panel, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["taxid", "name", "old_parent", "new_parent"])
        a_idx = MSL_ORDER.index(best_trans.split("->")[0].lower().replace("msl", "msl"))
        for t in sorted(best_list):
            w.writerow([t, names.get(t, "?"), "", ""])
    print(f"Saved panel candidate list: {out_panel}")

if __name__ == "__main__":
    main()
