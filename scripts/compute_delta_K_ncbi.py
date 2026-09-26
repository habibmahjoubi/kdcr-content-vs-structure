#!/usr/bin/env python3
"""Compute delta_K directly from consecutive NCBI taxdump nodes.dmp snapshots: predictor and response should be drawn from the same taxonomic source (NCBI),
since Kraken2 classification is driven by NCBI snapshots, not ICTV MSL records directly
(ICTV releases are implemented in NCBI with a lag).
"""
import csv
import os
import sys

BASE = "/home/hm/kdcr"
MSL_ORDER = [f"msl{n}" for n in range(29, 42)]

def load_nodes(nodes_dmp_path):
    """Return dict taxid(int) -> parent_taxid(int)."""
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

def viral_descendant_set(nodes, viruses_taxid=10239):
    """Return the set of taxids that are taxid 10239 (Viruses) or a descendant of it,
    via iterative memoized parent-walk (no recursion -- full NCBI taxonomy has millions
    of nodes). O(n) amortized: each node's ancestor chain is walked at most once thanks
    to path caching."""
    is_viral = {}
    for start in nodes.keys():
        if start in is_viral:
            continue
        path = []
        t = start
        while True:
            if t in is_viral:
                result = is_viral[t]
                break
            if t == viruses_taxid:
                result = True
                break
            parent = nodes.get(t)
            if parent is None or parent == t:
                result = False
                break
            path.append(t)
            t = parent
        for p in path:
            is_viral[p] = result
        is_viral[start] = result
    return {t for t, v in is_viral.items() if v}

def main():
    print(f"Loading nodes.dmp for {len(MSL_ORDER)} MSL versions...", flush=True)
    node_maps = {}
    for msl in MSL_ORDER:
        path = os.path.join(BASE, "kraken2_dbs", msl, "taxonomy", "nodes.dmp")
        node_maps[msl] = load_nodes(path)
        print(f"  {msl}: {len(node_maps[msl])} taxids", flush=True)

    rows = []
    for i in range(len(MSL_ORDER) - 1):
        a, b = MSL_ORDER[i], MSL_ORDER[i + 1]
        old, new = node_maps[a], node_maps[b]
        old_ids, new_ids = set(old.keys()), set(new.keys())
        added = new_ids - old_ids
        removed = old_ids - new_ids
        common = old_ids & new_ids
        reparented = {t for t in common if old[t] != new[t]}
        union = old_ids | new_ids
        changed = added | removed | reparented
        delta_full = len(changed) / len(union)

        # Viral-only restriction (union of viral-taxon sets in both snapshots)
        viral_old = viral_descendant_set(old)
        viral_new = viral_descendant_set(new)
        viral_union_ids = (old_ids | new_ids) & (viral_old | viral_new)
        changed_viral = changed & viral_union_ids
        delta_viral = (len(changed_viral) / len(viral_union_ids)) if viral_union_ids else float("nan")

        trans = f"MSL{a[3:]}->MSL{b[3:]}"
        rows.append({
            "transition": trans,
            "delta_K_ncbi_full": round(delta_full, 5),
            "delta_K_ncbi_viral": round(delta_viral, 5),
            "n_added": len(added), "n_removed": len(removed), "n_reparented": len(reparented),
            "n_union": len(union), "n_viral_union": len(viral_union_ids),
        })
        print(f"{trans}: full={delta_full:.5f} viral={delta_viral:.5f} "
              f"(added={len(added)} removed={len(removed)} reparented={len(reparented)} "
              f"union={len(union)} viral_union={len(viral_union_ids)})", flush=True)

    out_csv = os.path.join(BASE, "results", "delta_K_ncbi_transitions.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved: {out_csv}")

    # --- Merge with existing classification + ICTV delta_K, compute correlations ---
    from scipy import stats

    rank_matched = {}
    with open(os.path.join(BASE, "results", "classification_rank_matched.csv"), encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            rank_matched[rec["transition"]] = rec

    ictv_dk = {}
    with open(os.path.join(BASE, "results", "delta_K_species_extended_MSL30_41.csv"), encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            ictv_dk[rec["transition"]] = float(rec["delta_K_species"])

    merged = []
    for r in rows:
        t = r["transition"]
        rm = rank_matched.get(t, {})
        merged.append({
            "transition": t,
            "delta_K_species_ICTV": ictv_dk.get(t),
            "delta_K_ncbi_full": r["delta_K_ncbi_full"],
            "delta_K_ncbi_viral": r["delta_K_ncbi_viral"],
            "abs_change_RSV": float(rm.get("abs_change_RSV", "nan")),
            "abs_change_REO": float(rm.get("abs_change_REO_rank_matched", "nan")),
        })

    print("\n" + "=" * 100)
    print(f"{'Transition':<16}{'dK_ICTV':<10}{'dK_ncbi_full':<14}{'dK_ncbi_viral':<14}{'|dRSV|':<8}{'|dREO|':<8}")
    for m in merged:
        print(f"{m['transition']:<16}{m['delta_K_species_ICTV']:<10}{m['delta_K_ncbi_full']:<14}"
              f"{m['delta_K_ncbi_viral']:<14}{m['abs_change_RSV']:<8}{m['abs_change_REO']:<8}")

    def corr(x, y, label):
        rho, p = stats.spearmanr(x, y)
        print(f"{label}: rho={rho:.4f} p={p:.4f} n={len(x)}")
        return rho, p

    dk_ictv = [m["delta_K_species_ICTV"] for m in merged]
    dk_ncbi_full = [m["delta_K_ncbi_full"] for m in merged]
    dk_ncbi_viral = [m["delta_K_ncbi_viral"] for m in merged]
    rsv = [m["abs_change_RSV"] for m in merged]
    reo = [m["abs_change_REO"] for m in merged]

    print("\n--- Correlations ---")
    corr(dk_ncbi_full, rsv, "delta_K_ncbi_full vs RSV")
    corr(dk_ncbi_full, reo, "delta_K_ncbi_full vs REO")
    corr(dk_ncbi_viral, rsv, "delta_K_ncbi_viral vs RSV")
    corr(dk_ncbi_viral, reo, "delta_K_ncbi_viral vs REO")
    print()
    corr(dk_ictv, dk_ncbi_full, "delta_K_ICTV vs delta_K_ncbi_full (agreement in ranking)")
    corr(dk_ictv, dk_ncbi_viral, "delta_K_ICTV vs delta_K_ncbi_viral (agreement in ranking)")

    out_merged = os.path.join(BASE, "results", "delta_K_ncbi_vs_ictv_merged.csv")
    with open(out_merged, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(merged[0].keys()))
        w.writeheader()
        w.writerows(merged)
    print(f"Saved: {out_merged}")

if __name__ == "__main__":
    main()
