#!/usr/bin/env python3
"""Streaming summary of the all-update reparenting simulations (stage s7). The Kraken2 outputs of the
13 Design 2 databases list the simulated reads in the same order, so consecutive versions are read in
parallel line by line (constant memory). For each update, only reads of accessions whose taxon was
reparented on that update are evaluated. Read-level decomposition: same/different taxid x membership
of the clade of the accession's own taxid; 'exposed' = a minimizer hit outside that clade in either version."""
import os, csv
from collections import defaultdict, Counter

K = "/home/hm/kdcr"; W = f"{K}/work_designs"; R2 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
MSL = [f"msl{n}" for n in range(29, 42)]


def parents(m):
    p = {}
    for l in open(f"{K}/kraken2_dbs_fixed/{m}/taxonomy/nodes.dmp"):
        x = l.split("\t|\t", 2); p[int(x[0])] = int(x[1])
    return p


def clade_fn(par):
    cache = {}
    def anc(t):
        if t in cache: return cache[t]
        s, x = set(), t
        while x in par and x not in s:
            s.add(x)
            if par[x] == x: break
            x = par[x]
        cache[t] = frozenset(s); return cache[t]
    return lambda t, T: T in anc(t)


def parse(line):
    c = line.rstrip("\n").split("\t")
    hits = set(int(h.split(":")[0]) for h in c[4].replace("|:|", " ").split()
               if ":" in h and h.split(":")[0] not in ("0", "A"))
    return c[1], int(c[2]), hits


man = {r["update"]: set(r["accessions"].split(";")) if r["accessions"] else set()
       for r in csv.DictReader(open(f"{R2}/reparent_manifest.tsv"), delimiter="\t")}
acc_tax = {}
for l in open(f"{K}/kraken2_dbs_fixed/msl29/seqid2taxid.map"):
    s, t = l.split("\t"); acc_tax[s.split("|")[0]] = int(t)

rows, per_acc = [], []
for a, b in zip(MSL, MSL[1:]):
    u = f"{a}->{b}"; accs = man[u]
    ia, ib = clade_fn(parents(a)), clade_fn(parents(b))
    tot = Counter(); st = defaultdict(Counter)
    with open(f"{W}/reparent/cls_{a}_output.txt") as fa, open(f"{W}/reparent/cls_{b}_output.txt") as fb:
        for la, lb in zip(fa, fb):
            acc = la.split("\t", 2)[1].split("__")[0]
            if acc not in accs: continue
            ra, ta, ha = parse(la); rb, tb, hb = parse(lb)
            assert ra == rb
            T = acc_tax[acc]
            ca, cb = (ta != 0 and ia(ta, T)), (tb != 0 and ib(tb, T))
            exposed = any(not ia(h, T) for h in ha) or any(not ib(h, T) for h in hb)
            cat = ("same_taxid" if ta == tb else "different_taxid") + ("_same_membership" if ca == cb else "_membership_changed")
            tot[cat] += 1; tot["n"] += 1; tot["ca"] += ca; tot["cb"] += cb
            if ca != cb and not exposed: tot["unexposed_membership_change"] += 1
            if exposed: tot["exposed"] += 1
            s = st[acc]; s["n"] += 1; s["ca"] += ca; s["cb"] += cb; s[cat] += 1
    n = tot["n"] or 1
    big = [x for x, s in st.items() if abs(s["ca"] - s["cb"]) / s["n"] > 0.01]
    rows.append([u, len(accs), tot["n"], round(100 * tot["ca"] / n, 2), round(100 * tot["cb"] / n, 2),
                 tot["same_taxid_same_membership"], tot["same_taxid_membership_changed"],
                 tot["different_taxid_same_membership"], tot["different_taxid_membership_changed"],
                 tot["exposed"], tot["unexposed_membership_change"], len(big),
                 round(max((abs(s["ca"] - s["cb"]) / s["n"] * 100 for s in st.values()), default=0), 1)])
    print(rows[-1], flush=True)
    for x in big:
        s = st[x]
        per_acc.append([u, x, acc_tax[x], s["n"], round(100 * s["ca"] / s["n"], 1), round(100 * s["cb"] / s["n"], 1),
                        s["same_taxid_membership_changed"], s["different_taxid_membership_changed"] + s["different_taxid_same_membership"]])

H = ["update", "accessions", "read_pairs", "pct_in_own_clade_before", "pct_in_own_clade_after",
     "same_taxid_same_membership", "same_taxid_membership_changed", "different_taxid_same_membership",
     "different_taxid_membership_changed", "exposed_reads", "unexposed_membership_changes",
     "accessions_changing_more_than_1_point", "max_change_points"]
with open(os.path.join(R2, "reparent_all_updates.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(H); w.writerows(rows)
with open(os.path.join(R2, "reparent_accessions_changing.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["update", "accession", "taxid", "read_pairs", "pct_before", "pct_after",
                "reads_same_taxid_membership_changed", "reads_changing_taxid"]); w.writerows(per_acc)
