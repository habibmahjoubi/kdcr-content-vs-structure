#!/usr/bin/env python3
"""Summaries (run after pipeline_designs_panel.sh). Writes TSV tables to ../results.
Usage: python summarize_designs_panel.py [part ...]   parts: spiked conf panel official reparent phix factorial p56 static
Memory-light: one taxonomy at a time, parent links only."""
import os, sys, csv, re, glob
from collections import defaultdict, Counter

K = "/home/hm/kdcr"; W = f"{K}/work_designs"
R2 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
MSL = [f"msl{n}" for n in range(29, 42)]
RUNS14 = ["SRR26352207", "SRR26352206", "SRR26352204", "SRR26352215", "SRR26352205", "SRR26352203", "SRR26352202",
          "SRR26352217", "SRR26352216", "SRR26352195", "SRR26352184", "SRR26352212", "SRR26352209", "SRR26352208"]
SAMPLE = {"SRR26352207": "RSV neg", "SRR26352206": "RSV 1e6", "SRR26352204": "RSV 1e5", "SRR26352215": "RSV 1e4",
          "SRR26352205": "RSV 1e3", "SRR26352203": "RSV 1e2", "SRR26352202": "RSV 1e1",
          "SRR26352217": "MRV neg", "SRR26352216": "MRV 1e6", "SRR26352195": "MRV 1e5", "SRR26352184": "MRV 1e4",
          "SRR26352212": "MRV 1e3", "SRR26352209": "MRV 1e2", "SRR26352208": "MRV 1e1"}
OFFICIAL = ["20201202", "20210517", "20220607", "20220908", "20221209", "20230314", "20230605", "20231009", "20240112",
            "20240605", "20240904", "20241228", "20250402", "20250714", "20251015", "20260226", "20260626"]
RSV, MRV = 11250, 351073


def tsv(name, header, rows):
    with open(os.path.join(R2, name), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(header); w.writerows(rows)
    print("wrote", name, len(rows), flush=True)


def parents_nodes(path):
    p = {}
    with open(path) as f:
        for l in f:
            x = l.split("\t|\t", 2); p[int(x[0])] = int(x[1])
    return p


def parents_report(path):
    """Parent map reconstructed from a kraken2 --report-zero-counts report (DFS order, 2-space indent)."""
    p, stack = {}, []
    with open(path) as f:
        for l in f:
            c = l.rstrip("\n").split("\t")
            if len(c) < 6 or c[4] == "0": continue
            t = int(c[4]); name = c[5]; depth = (len(name) - len(name.lstrip(" "))) // 2
            while len(stack) > depth: stack.pop()
            p[t] = stack[-1] if stack else t
            stack.append(t)
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
        cache[t] = s; return s
    return lambda t, T: T in anc(t)


def cum(report, taxid):
    if not os.path.exists(report): return None
    with open(report) as f:
        for l in f:
            c = l.split("\t")
            if len(c) > 4 and c[4].strip() == str(taxid): return int(c[1])
    return 0


def outputs(path):
    """read -> (taxid, set of hit taxids)"""
    d = {}
    with open(path, errors="replace") as f:
        for l in f:
            c = l.rstrip("\n").split("\t")
            hits = set(int(h.split(":")[0]) for h in c[4].replace("|:|", " ").split() if ":" in h and h.split(":")[0] not in ("0", "A"))
            d[c[1]] = (int(c[2]), hits)
    return d


# ------------------------------------------------------------------ spiked runs: counts per design/version
def spiked():
    rows = []
    def rep(design, m, r):
        if design == "fixed":
            p = f"{K}/results_fixed/{m}/{r}_report.txt"
            return p if os.path.exists(p) else f"{W}/cls/fixed/{m}/{r}_report.txt"
        if design == "d1":
            p = f"{K}/results/{m}/{r}_report.txt"
            if not os.path.exists(p): p = f"{K}/work/SRR26352204/design1/{m}_report.txt" if r == "SRR26352204" else f"{W}/cls/d1/{m}/{r}_report.txt"
            return p
        if design == "d3": return f"{W}/cls/d3/{m}/{r}_report.txt"
    for design in ("fixed", "d1", "d3"):
        for r in RUNS14:
            rsv = [cum(rep(design, m, r), RSV) for m in MSL]; mrv = [cum(rep(design, m, r), MRV) for m in MSL]
            rows.append([design, r, SAMPLE[r], ",".join(map(str, rsv)), ",".join(map(str, mrv))])
    for d in OFFICIAL:
        for r in RUNS14:
            p = f"{W}/official/cls/{d}/{r}_report.txt"
            rows.append([f"official_{d}", r, SAMPLE[r], cum(p, RSV), cum(p, MRV)])
    for r in RUNS14:
        p = f"{K}/results/current_db/{r}_report.txt"
        if not os.path.exists(p): p = f"{W}/cls/current/{r}_report.txt"
        rows.append(["current", r, SAMPLE[r], cum(p, RSV), cum(p, MRV)])
    tsv("spiked_counts.tsv", ["design", "run", "sample", "RSV_11250_by_version", "MRV_351073_by_version"], rows)


def conf():
    rows = []
    for c in ("0.1", "0.5"):
        for design in ("fixed", "d1", "d3"):
            for r in RUNS14:
                rsv = [cum(f"{W}/conf/{design}/{c}/{m}/{r}_report.txt", RSV) for m in MSL]
                mrv = [cum(f"{W}/conf/{design}/{c}/{m}/{r}_report.txt", MRV) for m in MSL]
                rows.append([c, design, r, SAMPLE[r], ",".join(map(str, rsv)), ",".join(map(str, mrv))])
    tsv("confidence_counts.tsv", ["confidence", "design", "run", "sample", "RSV_by_version", "MRV_by_version"], rows)


# ------------------------------------------------------------------ panel: on-target fraction, changes, at-risk reads
def agent_targets():
    """agent -> {version: target taxid} using the Design 3 labels of the agent's reference sequences."""
    man = {r["agent"]: r["accessions"].split(";") for r in csv.DictReader(open(f"{R2}/panel_manifest.tsv"), delimiter="\t")}
    tg = defaultdict(dict)
    for m in MSL:
        lab = {l.split("\t")[0]: int(l.split("\t")[2]) for l in open(f"{W}/design3/{m}_labels.tsv")}
        for a, accs in man.items():
            ls = sorted(set(lab[x] for x in accs if x in lab))
            tg[a][m] = ls[0] if len(ls) == 1 else (ls if ls else None)
    return man, tg


def panel():
    man, tg = agent_targets()
    rows, chg, risk = [], [], []
    for design in ("fixed", "d1", "d3"):
        prev = None
        for m in MSL:
            par = parents_nodes(f"{K}/kraken2_dbs/{m}/taxonomy/nodes.dmp"); inc = clade_fn(par)
            out = outputs(f"{W}/panel/cls/{design}/{m}_output.txt")
            lib = set(l.split("\t")[0].split("|")[0] for l in open(
                {"fixed": f"{K}/kraken2_dbs_fixed/{m}/seqid2taxid.map", "d1": f"{K}/kraken2_dbs/{m}/seqid2taxid.map",
                 "d3": f"{W}/design3/{m}_labels.tsv"}[design]))
            cur = {}
            for rid, (t, hits) in out.items():
                a = rid.split("__")[0]; T = tg[a].get(m)
                on = isinstance(T, int) and t != 0 and inc(t, T)
                outside = isinstance(T, int) and any(not inc(h, T) for h in hits)
                cur[rid] = (t, on, outside)
            for a in man:
                ids = [r for r in cur if r.startswith(a + "__")]
                present = sum(1 for x in man[a] if x in lib)
                n_on = sum(1 for r in ids if cur[r][1]); n_cls = sum(1 for r in ids if cur[r][0] != 0)
                rows.append([design, m, a, tg[a].get(m), f"{present}/{len(man[a])}", len(ids), n_cls, n_on])
                if prev is not None:
                    changed = [r for r in ids if r in prev and prev[r][0] != cur[r][0]]
                    flip = [r for r in ids if r in prev and prev[r][1] != cur[r][1]]
                    flip_safe = [r for r in flip if not prev[r][2] and not cur[r][2]]
                    chg.append([design, f"{pm}->{m}", a, len(changed), len(flip), len(flip_safe),
                                sum(1 for r in ids if cur[r][2] or (r in prev and prev[r][2]))])
            prev, pm = cur, m
    tsv("panel_by_version.tsv", ["design", "version", "agent", "target_taxid", "reference_in_db", "read_pairs", "classified", "on_target"], rows)
    tsv("panel_changes.tsv", ["design", "update", "agent", "reads_changing_taxid", "reads_changing_target_membership",
                                "membership_changes_among_reads_without_outside_hits", "reads_with_outside_hits_either_version"], chg)


def official():
    man, tg = agent_targets()
    # official databases: target = lowest node of the agent's current lineage present in that database
    cur_par = parents_nodes(f"{K}/kraken2_dbs/current/taxonomy/nodes.dmp")
    lin = {}
    for a, accs in man.items():
        t = tg[a].get("msl41")
        if not isinstance(t, int): continue
        l, x = [], t
        while x in cur_par and cur_par[x] != x: l.append(x); x = cur_par[x]
        lin[a] = l
    rows, chg, prev = [], [], None
    for d in OFFICIAL:
        rp = f"{W}/official/cls/{d}/panel_report.txt"
        if not os.path.exists(rp): continue
        par = parents_report(rp); inc = clade_fn(par)
        out = outputs(f"{W}/official/cls/{d}/panel_output.txt")
        cur = {}
        for a in man:
            T = next((x for x in lin.get(a, []) if x in par), None)
            ids = [r for r in out if r.startswith(a + "__")]
            for r in ids:
                t, hits = out[r]; cur[r] = (t, T is not None and t != 0 and inc(t, T), T is not None and any(not inc(h, T) for h in hits))
            rows.append([d, a, T, len(ids), sum(1 for r in ids if out[r][0] != 0), sum(1 for r in ids if cur[r][1])])
            if prev is not None:
                changed = sum(1 for r in ids if r in prev and prev[r][0] != cur[r][0])
                flip = [r for r in ids if r in prev and prev[r][1] != cur[r][1]]
                chg.append([f"{pd}->{d}", a, changed, len(flip), sum(1 for r in flip if not prev[r][2] and not cur[r][2])])
        prev, pd = cur, d
    tsv("official_panel_by_version.tsv", ["index", "agent", "target_taxid", "read_pairs", "classified", "on_target"], rows)
    tsv("official_panel_changes.tsv", ["update", "agent", "reads_changing_taxid", "reads_changing_target_membership",
                                         "membership_changes_among_reads_without_outside_hits"], chg)
    # spiked runs, read level, consecutive official indexes
    rr = []
    for a, b in zip(OFFICIAL, OFFICIAL[1:]):
        both = disc = 0
        for r in RUNS14:
            pa, pb = f"{W}/official/cls/{a}/{r}_output.txt", f"{W}/official/cls/{b}/{r}_output.txt"
            if not (os.path.exists(pa) and os.path.exists(pb)): continue
            da = {l.split("\t", 3)[1]: l.split("\t", 3)[2] for l in open(pa) if l[0] == "C"}
            db = {l.split("\t", 3)[1]: l.split("\t", 3)[2] for l in open(pb) if l[0] == "C"}
            for k, v in da.items():
                if k in db: both += 1; disc += (v != db[k])
        rr.append([f"{a}->{b}", both, disc, round(100 * disc / both, 3) if both else ""])
    tsv("official_read_reassignment.tsv", ["update", "classified_in_both", "discordant", "pct"], rr)


# ------------------------------------------------------------------ reparenting, all updates, read level
def reparent():
    man = {r["update"]: r["accessions"].split(";") if r["accessions"] else [] for r in csv.DictReader(open(f"{R2}/reparent_manifest.tsv"), delimiter="\t")}
    acc_tax = {}
    for l in open(f"{K}/kraken2_dbs_fixed/msl29/seqid2taxid.map"):
        s, t = l.split("\t"); acc_tax[s.split("|")[0]] = int(t)
    rows, per_acc = [], []
    for u, accs in man.items():
        a, b = u.split("->"); accs = set(accs)
        pa, pb = parents_nodes(f"{K}/kraken2_dbs_fixed/{a}/taxonomy/nodes.dmp"), parents_nodes(f"{K}/kraken2_dbs_fixed/{b}/taxonomy/nodes.dmp")
        ia, ib = clade_fn(pa), clade_fn(pb)
        oa, ob = outputs(f"{W}/reparent/cls_{a}_output.txt"), outputs(f"{W}/reparent/cls_{b}_output.txt")
        tot = Counter(); acc_stats = defaultdict(Counter)
        for rid in oa:
            acc = rid.split("__")[0]
            if acc not in accs or rid not in ob: continue
            T = acc_tax[acc]; (ta, ha), (tb, hb) = oa[rid], ob[rid]
            ca, cb = ta != 0 and ia(ta, T), tb != 0 and ib(tb, T)
            safe = not any(not ia(h, T) for h in ha) and not any(not ib(h, T) for h in hb)
            cat = ("same_taxid" if ta == tb else "different_taxid") + ("_same_membership" if ca == cb else "_membership_changed")
            tot[cat] += 1; tot["n"] += 1; tot["correct_a"] += ca; tot["correct_b"] += cb
            if ca != cb and safe: tot["membership_change_without_outside_hits"] += 1
            acc_stats[acc][cat] += 1; acc_stats[acc]["n"] += 1; acc_stats[acc]["ca"] += ca; acc_stats[acc]["cb"] += cb
        n = tot["n"] or 1
        big = [x for x, s in acc_stats.items() if s["n"] and abs(s["ca"] - s["cb"]) / s["n"] > 0.01]
        rows.append([u, len(accs), tot["n"], round(100 * tot["correct_a"] / n, 2), round(100 * tot["correct_b"] / n, 2),
                     tot["same_taxid_same_membership"], tot["same_taxid_membership_changed"],
                     tot["different_taxid_same_membership"], tot["different_taxid_membership_changed"],
                     tot["membership_change_without_outside_hits"], len(big)])
        for x in big:
            s = acc_stats[x]
            per_acc.append([u, x, acc_tax[x], s["n"], round(100 * s["ca"] / s["n"], 1), round(100 * s["cb"] / s["n"], 1),
                            s["same_taxid_membership_changed"], s["different_taxid_membership_changed"] + s["different_taxid_same_membership"]])
    tsv("reparent_all_updates.tsv", ["update", "accessions", "read_pairs", "pct_in_own_clade_before", "pct_in_own_clade_after",
                                       "same_taxid_same_membership", "same_taxid_membership_changed", "different_taxid_same_membership",
                                       "different_taxid_membership_changed", "membership_change_without_outside_hits",
                                       "accessions_changing_more_than_1_point"], rows)
    tsv("reparent_accessions_changing.tsv", ["update", "accession", "taxid", "read_pairs", "pct_before", "pct_after",
                                               "reads_same_taxid_but_membership_changed", "reads_changing_taxid"], per_acc)


# ------------------------------------------------------------------ phiX-excluded read-level reassignment (fixed design)
def phix():
    runs = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195", "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]
    phx = {r: set(l.strip() for l in open(f"{W}/phix/{r}.phix_ids")) for r in runs}
    rows = []
    for a, b in zip(MSL, MSL[1:]):
        both = disc = both_np = disc_np = 0
        for r in runs:
            da = {l.split("\t", 3)[1]: l.split("\t", 3)[2] for l in open(f"{K}/results_fixed/{a}/{r}_output.txt") if l[0] == "C"}
            db = {l.split("\t", 3)[1]: l.split("\t", 3)[2] for l in open(f"{K}/results_fixed/{b}/{r}_output.txt") if l[0] == "C"}
            for k, v in da.items():
                if k not in db: continue
                both += 1; d = v != db[k]; disc += d
                if k not in phx[r]: both_np += 1; disc_np += d
        rows.append([f"{a}->{b}", both, disc, round(100 * disc / both, 3), both_np, disc_np, round(100 * disc_np / both_np, 4)])
    tsv("reassignment_without_phix.tsv", ["update", "classified_both", "discordant", "pct", "classified_both_non_phix", "discordant_non_phix", "pct_non_phix"], rows)


# ------------------------------------------------------------------ factorial and determinism
def factorial():
    F = f"{W}/factorial/cls"
    dbs = {"msl41 (MSL41 content, MSL41 tree)": f"{K}/results/msl41", "current (MSL41+1 content, current tree)": f"{K}/results/current_db",
           "a (MSL41 content, current tree)": f"{F}/a_msl41content_currenttax", "b (MSL41+1 content, MSL41 tree)": f"{F}/b_msl41plus1_msl41tax",
           "c (MSL41 rebuilt, 2 threads)": f"{F}/c_msl41_rebuild_2threads", "d (current rebuilt, 2 threads)": f"{F}/d_current_rebuild_2threads"}
    cur_par = parents_nodes(f"{K}/kraken2_dbs/current/taxonomy/nodes.dmp"); inc = clade_fn(cur_par)
    rows, sets = [], {}
    for name, base in dbs.items():
        s = set()
        for r in RUNS14:
            p = f"{base}/{r}_output.txt"
            if not os.path.exists(p): p = f"{W}/cls/current/{r}_output.txt" if "current" in base else p
            if not os.path.exists(p): continue
            for l in open(p):
                if l[0] != "C": continue
                c = l.split("\t", 3)
                if inc(int(c[2]), MRV): s.add(c[1])
        sets[name] = s
        rows.append([name, len(s)])
    names = list(sets)
    ov = [[a, b, len(sets[a] & sets[b]), len(sets[a] ^ sets[b])] for i, a in enumerate(names) for b in names[i + 1:]]
    tsv("factorial_mrv_clade_reads.tsv", ["database", "reads_in_MRV_clade_all_14_runs"], rows)
    tsv("factorial_overlap.tsv", ["database_a", "database_b", "shared_reads", "reads_in_only_one"], ov)


# ------------------------------------------------------------------ 56-sequence panel: destination of misassigned reads
def p56():
    rows = []
    prev = None
    for m in ("msl37", "msl38", "msl39"):
        out = {l.split("\t", 3)[1]: int(l.split("\t", 3)[2]) for l in open(f"{W}/p56/d1_{m}_output.txt") if l[0] == "C"}
        if prev is not None:
            common = set(out) & set(prev)
            rows.append([f"{pm}->{m}", len(common), sum(1 for r in common if out[r] != prev[r])])
        prev, pm = out, m
    tsv("p56_misassigned_destination_changes.tsv", ["update", "reads_classified_in_both_(genome_absent)", "reads_changing_destination_taxid"], rows)


# ------------------------------------------------------------------ static invariance check: sequences under each target
def static():
    man, tg = agent_targets()
    rows = []
    for design in ("fixed", "d1", "d3"):
        prev = {}
        for m in MSL:
            par = parents_nodes(f"{K}/kraken2_dbs/{m}/taxonomy/nodes.dmp"); inc = clade_fn(par)
            if design == "d3":
                lab = [(l.split("\t")[0], int(l.split("\t")[2])) for l in open(f"{W}/design3/{m}_labels.tsv")]
            else:
                mp = f"{K}/kraken2_dbs_fixed/{m}/seqid2taxid.map" if design == "fixed" else f"{K}/kraken2_dbs/{m}/seqid2taxid.map"
                lab = [(l.split("\t")[0].split("|")[0], int(l.split("\t")[1])) for l in open(mp)]
            lab = [(a, t) for a, t in lab if t in par]
            targets = {"RSV_target_11250": RSV, "MRV_target_351073": MRV}
            for a in man:
                if isinstance(tg[a].get(m), int): targets[a] = tg[a][m]
            cur = {}
            for name, T in targets.items():
                S = frozenset(a for a, t in lab if inc(t, T)); cur[name] = S
                if m != MSL[0] and name in prev:
                    rows.append([design, f"{pm}->{m}", name, len(prev[name]), len(S), len(prev[name] ^ S)])
            prev, pm = cur, m
    tsv("static_sequences_under_target.tsv", ["design", "update", "target", "sequences_before", "sequences_after", "sequences_entering_or_leaving"], rows)


if __name__ == "__main__":
    for p in (sys.argv[1:] or ["spiked", "conf", "panel", "official", "reparent", "phix", "factorial", "p56", "static"]):
        try:
            globals()[p]()
        except Exception as e:
            print("FAILED", p, type(e).__name__, e, flush=True)
