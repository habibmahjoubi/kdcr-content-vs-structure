#!/usr/bin/env python3
"""Candidate 1 (node-relative existential + sibling neighborhood), rank-matched taxids."""
import os

BASE = "/home/hm/kdcr/kraken2_dbs"
TRANSITIONS = [(a, a + 1) for a in range(29, 41)]
TARGETS = {"RSV": "11250", "REO": "351073"}


def load_nodes(msl):
    path = os.path.join(BASE, f"msl{msl}", "taxonomy", "nodes.dmp")
    nodes = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            nodes[parts[0].strip()] = parts[1].strip()
    return nodes


def compute(nodes_a, nodes_b, tau):
    in_a, in_b = tau in nodes_a, tau in nodes_b
    if in_a != in_b or (in_a and in_b and nodes_a[tau] != nodes_b[tau]):
        existential = 1
    else:
        existential = 0
    parent_a, parent_b = nodes_a.get(tau), nodes_b.get(tau)
    sib_a = {t for t, p in nodes_a.items() if p == parent_a} if parent_a else set()
    sib_b = {t for t, p in nodes_b.items() if p == parent_b} if parent_b else set()
    sibs = (sib_a | sib_b) - {tau}
    if not sibs:
        return existential, 0.0
    changed = sum(1 for s in sibs if (s in nodes_a) != (s in nodes_b) or
                  (s in nodes_a and s in nodes_b and nodes_a[s] != nodes_b[s]))
    return existential, changed / len(sibs)


cache = {}
for prev_v, new_v in TRANSITIONS:
    if prev_v not in cache:
        cache[prev_v] = load_nodes(prev_v)
    if new_v not in cache:
        cache[new_v] = load_nodes(new_v)
    row = {"transition": f"MSL{prev_v}->MSL{new_v}"}
    for name, tau in TARGETS.items():
        e, n = compute(cache[prev_v], cache[new_v], tau)
        row[f"{name}_existential"] = e
        row[f"{name}_neighborhood"] = round(n, 4)
    print(row)
    if prev_v - 1 in cache:
        del cache[prev_v - 1]
