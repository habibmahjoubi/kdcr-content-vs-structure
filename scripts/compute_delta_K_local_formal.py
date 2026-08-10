#!/usr/bin/env python3
"""Formalize delta_K_local: a knowledge-drift measure computed directly from the NCBI
taxonomy graph (nodes.dmp parent-child structure), working uniformly across all twelve real
transitions (MSL29-41).

REVISED DEFINITION (v2): since classification counts at a target taxid are CUMULATIVE
(Kraken2 sums the target node plus all its descendants), the drift measure that should
predict a change in that cumulative count is drift in the DESCENDANT set below the target,
not drift in the target's own position relative to its parent. A first version of this
script measured the latter (existential change of tau's own node + sibling churn) and it
did not track observed read-count changes at all -- notably, MSL37->MSL38 scored the
highest "local" value for REO under that definition (tau's genus itself was reclassified
under the family-level megataxonomy reorganization) yet produced zero observed read-count
change, because none of Orthoreovirus's own descendant species changed that year. This
revised version measures exactly that quantity instead:

  descendant_drift(tau, e) = |descendants(tau, v+) ^ descendants(tau, v-)| / |descendants(tau, v-)|

i.e. the proportion of tau's descendant taxa (recursively, including tau itself) that were
added or removed between v- and v+, which is the quantity that can mechanically change how
many reference sequences a classifier can assign to tau's clade.
"""
import os
import csv

BASE = "/home/hm/kdcr/kraken2_dbs"
OUT_CSV = "/home/hm/kdcr/results/delta_K_local_formal.csv"

TRANSITIONS = [(a, a + 1) for a in range(29, 41)]
TARGETS = {"RSV": "11250", "REO": "10882"}


def load_nodes(msl):
    path = os.path.join(BASE, f"msl{msl}", "taxonomy", "nodes.dmp")
    parent_of = {}
    children_of = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = parts[0].strip()
            parent = parts[1].strip()
            parent_of[taxid] = parent
            children_of.setdefault(parent, set()).add(taxid)
    return parent_of, children_of


def descendants(tau, children_of):
    """All descendant taxids of tau (inclusive), via BFS over the children index."""
    seen = {tau}
    frontier = [tau]
    while frontier:
        nxt = []
        for t in frontier:
            for c in children_of.get(t, ()):
                if c not in seen:
                    seen.add(c)
                    nxt.append(c)
        frontier = nxt
    return seen


def main():
    rows = []
    cache = {}
    for prev_v, new_v in TRANSITIONS:
        if prev_v not in cache:
            cache[prev_v] = load_nodes(prev_v)
        if new_v not in cache:
            cache[new_v] = load_nodes(new_v)
        _, children_a = cache[prev_v]
        _, children_b = cache[new_v]

        row = {"transition": f"MSL{prev_v}->MSL{new_v}"}
        for name, tau in TARGETS.items():
            desc_a = descendants(tau, children_a) if tau in children_a or any(True for _ in [tau]) else {tau}
            desc_b = descendants(tau, children_b) if tau in children_b or any(True for _ in [tau]) else {tau}
            # descendants() needs tau present as a node at all; guard for tau missing entirely
            desc_a = descendants(tau, children_a)
            desc_b = descendants(tau, children_b)
            union_diff = desc_a ^ desc_b
            drift = len(union_diff) / len(desc_a) if desc_a else 0.0
            row[f"{name}_n_descendants_before"] = len(desc_a)
            row[f"{name}_n_descendants_after"] = len(desc_b)
            row[f"{name}_delta_K_local"] = round(drift, 4)
        rows.append(row)
        print(row)
        if prev_v - 1 in cache:
            del cache[prev_v - 1]

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved: {OUT_CSV}")


if __name__ == "__main__":
    main()
