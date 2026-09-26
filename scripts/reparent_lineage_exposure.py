#!/usr/bin/env python3
"""Reparenting simulations (Design 2, all twelve updates), re-evaluated with the lineage-based definition of
exposure: a read is exposed if, in either version, at least one minimizer is stored at a taxon that is neither
within the clade of the accession's own taxid T nor an ancestor of T (hits stored at ancestors of T score for
every path through T and cannot draw a read out of the clade). Membership changes are split by whether the set
of library sequences under T changed on the update. Streaming, constant memory (as reparent_summary.py)."""
import csv
from collections import Counter, defaultdict

import os
K = "/home/hm/kdcr"; W = f"{K}/work_designs"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
R2T = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
MSL = [f"msl{n}" for n in range(29, 42)]


def parents(m):
    p = {}
    for l in open(f"{K}/kraken2_dbs_fixed/{m}/taxonomy/nodes.dmp"):
        x = l.split("\t|\t", 2); p[int(x[0])] = int(x[1])
    return p


def ancestors_fn(par):
    cache = {}
    def anc(t):   # t and all its ancestors
        if t in cache: return cache[t]
        s, x = [], t
        seen = set()
        while x in par and x not in seen:
            s.append(x); seen.add(x)
            if par[x] == x: break
            x = par[x]
        cache[t] = frozenset(s); return cache[t]
    return anc


def parse(line):
    c = line.rstrip("\n").split("\t")
    hits = set(int(h.split(":")[0]) for h in c[4].replace("|:|", " ").split()
               if ":" in h and h.split(":")[0] not in ("0", "A"))
    return c[1], int(c[2]), hits


man = {r["update"]: set(r["accessions"].split(";")) if r["accessions"] else set()
       for r in csv.DictReader(open(f"{R2T}/reparent_manifest.tsv"), delimiter="\t")}
acc_tax = {}
for l in open(f"{K}/kraken2_dbs_fixed/msl29/seqid2taxid.map"):
    s, t = l.split("\t"); acc_tax[s.split("|")[0]] = int(t)


def under(anc):
    u = defaultdict(set)
    for a, t in acc_tax.items():
        for x in anc(t): u[x].add(a)
    return u


rows = []
for a, b in zip(MSL, MSL[1:]):
    u = f"{a}->{b}"; accs = man[u]
    aa, ab = ancestors_fn(parents(a)), ancestors_fn(parents(b))
    ua, ub = under(aa), under(ab)
    tot = Counter()
    with open(f"{W}/reparent/cls_{a}_output.txt") as fa, open(f"{W}/reparent/cls_{b}_output.txt") as fb:
        for la, lb in zip(fa, fb):
            acc = la.split("\t", 2)[1].split("__")[0]
            if acc not in accs: continue
            ra, ta, ha = parse(la); rb, tb, hb = parse(lb)
            assert ra == rb
            T = acc_tax[acc]
            ca, cb = (ta != 0 and T in aa(ta)), (tb != 0 and T in ab(tb))
            ancTa, ancTb = aa(T) - {T}, ab(T) - {T}
            ex_any = any(T not in aa(h) for h in ha) or any(T not in ab(h) for h in hb)
            ex_lin = any(T not in aa(h) and h not in ancTa for h in ha) or any(T not in ab(h) and h not in ancTb for h in hb)
            same_seqs = ua.get(T, set()) == ub.get(T, set())
            tot["n"] += 1; tot["exposed_any"] += ex_any; tot["exposed_lineage"] += ex_lin
            if ca != cb:
                tot["membership_changes"] += 1
                tot["changes_same_sequences_under_T" if same_seqs else "changes_sequences_under_T_changed"] += 1
                if not ex_lin:
                    tot["unexposed_lineage_changes_same_seqs" if same_seqs else "unexposed_lineage_changes_seqs_changed"] += 1
    rows.append([u, len(accs), tot["n"], tot["exposed_any"], tot["exposed_lineage"], tot["membership_changes"],
                 tot["changes_same_sequences_under_T"], tot["changes_sequences_under_T_changed"],
                 tot["unexposed_lineage_changes_same_seqs"], tot["unexposed_lineage_changes_seqs_changed"]])
    print(rows[-1], flush=True)

with open(f"{OUT}/reparent_lineage_exposure.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["update", "accessions", "read_pairs", "exposed_any_hit_outside_clade", "exposed_lineage", "membership_changes",
                "changes_sequences_under_T_unchanged", "changes_sequences_under_T_changed",
                "changes_without_lineage_exposure_sequences_unchanged", "changes_without_lineage_exposure_sequences_changed"])
    w.writerows(rows)
