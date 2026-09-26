#!/usr/bin/env python3
"""Exposure, host-read and prebuilt-database analyses. Run after pipeline_exposure_host.sh.
Usage: python exposure_host_prebuilt.py part [part ...]
parts: targets exposure t1l collect spurious oconf jumps vesi determinism denominators
Memory-light: one taxonomy at a time, parent links only."""
import os, sys, csv, subprocess
from collections import defaultdict, Counter

K = "/home/hm/kdcr"; W2 = f"{K}/work_designs"; W3 = f"{K}/work_exposure"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
R2T = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
KR = "/home/hm/miniforge3/envs/kdcr/bin/kraken2"
MSL = [f"msl{n}" for n in range(29, 42)]
OFFICIAL = ["20201202", "20210517", "20220607", "20220908", "20221209", "20230314", "20230605", "20231009", "20240112",
            "20240605", "20240904", "20241228", "20250402", "20250714", "20251015", "20260226", "20260626"]
RUNS14 = ["SRR26352207", "SRR26352206", "SRR26352204", "SRR26352215", "SRR26352205", "SRR26352203", "SRR26352202",
          "SRR26352217", "SRR26352216", "SRR26352195", "SRR26352184", "SRR26352212", "SRR26352209", "SRR26352208"]
SAMPLE = {"SRR26352207": "RSV neg", "SRR26352206": "RSV 1e6", "SRR26352204": "RSV 1e5", "SRR26352215": "RSV 1e4",
          "SRR26352205": "RSV 1e3", "SRR26352203": "RSV 1e2", "SRR26352202": "RSV 1e1",
          "SRR26352217": "MRV neg", "SRR26352216": "MRV 1e6", "SRR26352195": "MRV 1e5", "SRR26352184": "MRV 1e4",
          "SRR26352212": "MRV 1e3", "SRR26352209": "MRV 1e2", "SRR26352208": "MRV 1e1"}
RSV, MRV = 11250, 351073
FQ = f"{K}/data/fastq"


def tsv(name, header, rows):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(header); w.writerows(rows)
    print("wrote", name, len(rows), flush=True)


def parents_nodes(path):
    p = {}
    with open(path) as f:
        for l in f:
            x = l.split("\t|\t", 2); p[int(x[0])] = int(x[1])
    return p


def names_of(path, wanted):
    n = {}
    with open(path) as f:
        for l in f:
            x = l.split("\t|\t")
            t = int(x[0])
            if t in wanted and x[3].startswith("scientific name"): n[t] = x[1]
    return n


def parents_report(path):
    p, stack, names = {}, [], {}
    with open(path) as f:
        for l in f:
            c = l.rstrip("\n").split("\t")
            if len(c) < 6 or c[4] == "0": continue
            t = int(c[4]); name = c[5]; depth = (len(name) - len(name.lstrip(" "))) // 2
            while len(stack) > depth: stack.pop()
            p[t] = stack[-1] if stack else t
            names[t] = name.strip(); stack.append(t)
    return p, names


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


def outputs(path):
    d = {}
    with open(path, errors="replace") as f:
        for l in f:
            c = l.rstrip("\n").split("\t")
            hits = set(int(h.split(":")[0]) for h in c[4].replace("|:|", " ").split() if ":" in h and h.split(":")[0] not in ("0", "A"))
            d[c[1]] = (int(c[2]), hits)
    return d


def tax_path(design, v):
    if design == "fixed": return f"{K}/kraken2_dbs_fixed/{v}/taxonomy/nodes.dmp"
    if design in ("d1", "d3"): return f"{K}/kraken2_dbs/{v}/taxonomy/nodes.dmp"
    if design == "current": return f"{K}/kraken2_dbs/current/taxonomy/nodes.dmp"


# ------------------------------------------------------------------ panel manifest and Design 3 targets
def manifest():
    return {r["agent"]: r["accessions"].split(";") for r in csv.DictReader(open(f"{R2T}/panel_manifest.tsv"), delimiter="\t")}


def d3_labels(m):
    return {l.split("\t")[0]: (int(l.split("\t")[1]), int(l.split("\t")[2])) for l in open(f"{W2}/design3/{m}_labels.tsv")}


def agent_targets():
    man = manifest(); tg = defaultdict(dict)
    for m in MSL:
        lab = d3_labels(m)
        for a, accs in man.items():
            ls = sorted(set(lab[x][1] for x in accs if x in lab))
            tg[a][m] = ls[0] if len(ls) == 1 else (ls if ls else None)
    return man, tg


def targets():
    """Design 3 target of each agent per version, with the labelling route (current / merged / ancestor)."""
    man = manifest(); rows = []
    for m in MSL:
        lab = d3_labels(m)
        mp = parents_nodes(tax_path("d1", m))
        merged = {}
        for l in open(f"{K}/kraken2_dbs/{m}/taxonomy/merged.dmp") if os.path.exists(f"{K}/kraken2_dbs/{m}/taxonomy/merged.dmp") else []:
            x = l.split("\t|\t"); merged[int(x[0])] = int(x[1].split("\t")[0])
        for a, accs in man.items():
            for x in accs:
                if x not in lab: continue
                cur, L = lab[x]
                route = "current" if L == cur else ("ancestor" if L in mp and _is_anc(L, cur, mp) else "merged_predecessor")
                rows.append([a, x, m, cur, L, route])
    tsv("design3_panel_labels.tsv", ["agent", "accession", "version", "current_taxid", "label_taxid", "route"], rows)


def _is_anc(L, cur, par):
    x, seen = cur, set()
    while x in par and x not in seen:
        if x == L: return True
        seen.add(x)
        if par[x] == x: break
        x = par[x]
    return False


# ------------------------------------------------------------------ exposure bound and confidence 0.1
def status(t, hits, T, inc, par):
    """(on target, exposed = any hit outside the clade of T, lineage-exposed = any hit on a node that is neither
    in the clade of T nor an ancestor of T). Hits stored at ancestors of T score for every path through T."""
    out = [h for h in hits if not inc(h, T)]
    anc, x = set(), T
    while out and x in par and par[x] != x: x = par[x]; anc.add(x)
    return (t != 0 and inc(t, T), bool(out), any(h not in anc for h in out))


def exposure():
    man, tg = agent_targets()
    per_ver, summ = [], []
    data = {}   # (scheme, agent) -> {version: {read: (on, exposed)}}
    for design in ("d3", "fixed", "d1"):
        for m in MSL:
            par = parents_nodes(tax_path(design, m)); inc = clade_fn(par)
            out = outputs(f"{W2}/panel/cls/{design}/{m}_output.txt")
            for rid, (t, hits) in out.items():
                a = rid.split("__")[0]; T = tg[a].get(m)
                if not isinstance(T, int): continue
                data.setdefault((design, a), {}).setdefault(m, {})[rid] = status(t, hits, T, inc, par)
            # MRV type 1 Lang reads (spiked strain, not in RefSeq): target = MRV node
            for rid, (t, hits) in outputs(f"{W3}/t1l/cls/{design}_{m}_output.txt").items():
                data.setdefault((design, "MRV1_Lang"), {}).setdefault(m, {})[rid] = status(t, hits, MRV, inc, par)
    # prebuilt: target = lowest node of the agent's current lineage present in the database
    cur_par = parents_nodes(tax_path("current", None)); lin = {}
    for a in man:
        t = tg[a].get("msl41")
        if isinstance(t, int):
            l, x = [], t
            while x in cur_par and cur_par[x] != x: l.append(x); x = cur_par[x]
            lin[a] = l
    for o in OFFICIAL:
        par, _ = parents_report(f"{W2}/official/cls/{o}/panel_report.txt"); inc = clade_fn(par)
        out = outputs(f"{W2}/official/cls/{o}/panel_output.txt")
        for a in lin:
            T = next((x for x in lin[a] if x in par), None)
            if T is None: continue
            for rid, (t, hits) in out.items():
                if rid.startswith(a + "__"):
                    data.setdefault(("prebuilt", a), {}).setdefault(o, {})[rid] = status(t, hits, T, inc, par)
        pt, _ = parents_report(f"{W3}/t1l/cls/official_{o}_report.txt"); pt = {**par, **pt}; it = clade_fn(pt)
        for rid, (t, hits) in outputs(f"{W3}/t1l/cls/official_{o}_output.txt").items():
            data.setdefault(("prebuilt", "MRV1_Lang"), {}).setdefault(o, {})[rid] = status(t, hits, MRV, it, pt)
    # confidence 0.1: on-target in Design 3 (13 versions), prebuilt (17) and current
    conf = defaultdict(dict)
    for scheme, vers in (("d3", MSL), ("prebuilt", OFFICIAL)):
        for v in vers:
            p = f"{W3}/pconf/{'d3_' + v if scheme == 'd3' else 'official_' + v}_output.txt"
            if not os.path.exists(p): continue
            if scheme == "d3":
                inc = clade_fn(parents_nodes(tax_path("d3", v))); T_of = lambda a: tg[a].get(v)
            else:
                par, _ = parents_report(f"{W3}/pconf/official_{v}_report.txt"); inc = clade_fn(par)
                par0, _ = parents_report(f"{W2}/official/cls/{v}/panel_report.txt")
                T_of = lambda a: next((x for x in lin.get(a, []) if x in par0), None)
                inc0 = clade_fn(par0)
            n_on = Counter()
            for l in open(p):
                c = l.split("\t", 4); a = c[1].split("__")[0]; T = T_of(a); t = int(c[2])
                if not isinstance(T, int) or t == 0: continue
                if (inc(t, T) if scheme == "d3" else (t in par0 and inc0(t, T))): n_on[a] += 1
            for a in man: conf[(scheme, a)][v] = n_on[a]
    for (scheme, a), vers in sorted(data.items()):
        vs = [v for v in (MSL if scheme != "prebuilt" else OFFICIAL) if v in vers]
        n = len(next(iter(vers.values())))
        on = {v: sum(1 for r in vers[v].values() if r[0]) for v in vs}
        ex = {v: sum(1 for r in vers[v].values() if r[1]) for v in vs}
        lx = {v: sum(1 for r in vers[v].values() if r[2]) for v in vs}
        for v in vs: per_ver.append([scheme, a, v, n, on[v], ex[v], lx[v]])
        # tight bound: membership changes between any two versions vs reads lineage-exposed in either version
        lmax_flip, lbound, lviol = 0, 0, 0
        for i, v1 in enumerate(vs):
            for v2 in vs[i + 1:]:
                common = [r for r in vers[v1] if r in vers[v2]]
                fl = [r for r in common if vers[v1][r][0] != vers[v2][r][0]]
                lviol += sum(1 for r in fl if not vers[v1][r][2] and not vers[v2][r][2])
                if len(fl) >= lmax_flip:
                    lmax_flip = len(fl); lbound = sum(1 for r in common if vers[v1][r][2] or vers[v2][r][2])
        lx_any = len(set(r for v in vs for r, x in vers[v].items() if x[2]))
        ex_any = set(r for v in vs for r, x in vers[v].items() if x[1])
        # bound on every pair of versions: reads changing membership vs reads exposed in either version
        worst_ratio, max_flip, max_bound_pair, violations = 0, 0, None, 0
        for i, v1 in enumerate(vs):
            for v2 in vs[i + 1:]:
                flips = [r for r in vers[v1] if r in vers[v2] and vers[v1][r][0] != vers[v2][r][0]]
                unexposed = [r for r in flips if not vers[v1][r][1] and not vers[v2][r][1]]
                violations += len(unexposed)
                if len(flips) > max_flip:
                    max_flip = len(flips)
                    max_bound_pair = (v1, v2, sum(1 for r in vers[v1] if r in vers[v2] and (vers[v1][r][1] or vers[v2][r][1])))
        consec = [abs(on[b] - on[a_]) for a_, b in zip(vs, vs[1:])]
        c01 = conf.get((scheme, a), {})
        summ.append([scheme, a, n, min(on.values()), max(on.values()), round(100 * min(ex.values()) / n, 2), round(100 * max(ex.values()) / n, 2),
                     len(ex_any), max(on.values()) - min(on.values()), max(consec) if consec else 0, max_flip,
                     max_bound_pair[2] if max_bound_pair else 0, violations,
                     (min(c01.values()) if c01 else ""), (max(c01.values()) if c01 else ""),
                     round(100 * min(lx.values()) / n, 2), round(100 * max(lx.values()) / n, 2), lx_any, lbound, lviol])
    tsv("panel_exposure_by_version.tsv", ["scheme", "agent", "version", "read_pairs", "on_target", "exposed", "lineage_exposed"], per_ver)
    tsv("panel_exposure_summary.tsv", ["scheme", "agent", "read_pairs", "on_target_min", "on_target_max", "pct_exposed_min", "pct_exposed_max",
                                          "reads_exposed_in_any_version", "on_target_range_all_versions", "largest_change_consecutive",
                                          "max_membership_changes_any_pair", "exposed_either_version_for_that_pair",
                                          "unexposed_membership_changes_all_pairs", "on_target_conf0.1_min", "on_target_conf0.1_max",
                                          "pct_lineage_exposed_min", "pct_lineage_exposed_max", "reads_lineage_exposed_in_any_version",
                                          "lineage_exposed_either_version_for_max_pair", "membership_changes_without_lineage_exposure_all_pairs"], summ)


# ------------------------------------------------------------------ MRV type 1 Lang simulated reads
def t1l():
    rows = []
    def count(out, inc, par_ok=None):
        n = cls = mrv = 0; dest = Counter()
        for l in open(out):
            c = l.split("\t", 4); n += 1; t = int(c[2])
            if t: cls += 1; dest[t] += 1
            if t and (par_ok is None or t in par_ok) and inc(t, MRV): mrv += 1
        return n, cls, mrv, dest
    for design in ("fixed", "d1", "d3"):
        for m in MSL:
            inc = clade_fn(parents_nodes(tax_path(design, m)))
            n, cls, mrv, dest = count(f"{W3}/t1l/cls/{design}_{m}_output.txt", inc)
            rows.append([design, m, n, cls, mrv, ";".join(f"{t}:{c}" for t, c in dest.most_common(4))])
    for o in OFFICIAL:
        par, _ = parents_report(f"{W3}/t1l/cls/official_{o}_report.txt")
        n, cls, mrv, dest = count(f"{W3}/t1l/cls/official_{o}_output.txt", clade_fn(par), par)
        rows.append(["prebuilt", o, n, cls, mrv, ";".join(f"{t}:{c}" for t, c in dest.most_common(4))])
    inc = clade_fn(parents_nodes(tax_path("current", None)))
    for tag, f in (("current", "current"), ("d3_msl41_conf0.1", "d3_msl41_c0.1")):
        n, cls, mrv, dest = count(f"{W3}/t1l/cls/{f}_output.txt", inc)
        rows.append([tag, "", n, cls, mrv, ";".join(f"{t}:{c}" for t, c in dest.most_common(4))])
    # per segment, Design 3 MSL41
    seg = defaultdict(lambda: [0, 0])
    for l in open(f"{W3}/t1l/cls/d3_msl41_output.txt"):
        c = l.split("\t", 4); s = c[1].rsplit("_", 5)[0]; seg[s][0] += 1; seg[s][1] += int(c[2]) != 0 and inc(int(c[2]), MRV)
    for s, (n, k) in sorted(seg.items()): rows.append(["d3_msl41_by_segment", s, n, "", k, ""])
    tsv("mrv_type1_lang_simulated.tsv", ["database", "version", "read_pairs", "classified", "in_MRV_clade_351073", "top_destinations"], rows)


# ------------------------------------------------------------------ reads at the RSV / MRV nodes in every database
def spiked_output(design, v, r):
    cands = {"fixed": [f"{K}/results_fixed/{v}/{r}_output.txt", f"{W2}/cls/fixed/{v}/{r}_output.txt"],
             "d1": [f"{K}/results/{v}/{r}_output.txt", f"{W2}/cls/d1/{v}/{r}_output.txt"],
             "d3": [f"{W2}/cls/d3/{v}/{r}_output.txt"],
             "prebuilt": [f"{W2}/official/cls/{v}/{r}_output.txt"],
             "current": [f"{K}/results/current_db/{r}_output.txt", f"{W2}/cls/current/{r}_output.txt"],
             "factorial": [f"{W2}/factorial/cls/{v}/{r}_output.txt"]}[design]
    return next((p for p in cands if os.path.exists(p)), None)


def clade_members(par, T):
    ch = defaultdict(list)
    for c, p in par.items():
        if c != p: ch[p].append(c)
    s, st = set(), [T]
    while st:
        x = st.pop()
        if x in s: continue
        s.add(x); st.extend(ch.get(x, []))
    return s


def collect():
    hits = defaultdict(set)   # (run, read) -> {"design:version:node"}
    missing = []
    jobs = [(d, v) for d in ("fixed", "d1", "d3") for v in MSL] + [("prebuilt", o) for o in OFFICIAL] + [("current", "current")] + \
           [("factorial", f) for f in ("a_msl41content_currenttax", "b_msl41plus1_msl41tax", "c_msl41_rebuild_2threads", "d_current_rebuild_2threads")]
    for design, v in jobs:
        if design == "prebuilt":
            par = {}
            for r in RUNS14:
                p, _ = parents_report(f"{W2}/official/cls/{v}/{r}_report.txt"); par.update(p)
        elif design == "factorial":
            par = parents_nodes(tax_path("current", None) if v[0] in "ad" else tax_path("d1", "msl41"))
        else:
            par = parents_nodes(tax_path(design, v))
        sets = {"RSV": clade_members(par, RSV), "MRV": clade_members(par, MRV)}; del par
        for r in RUNS14:
            p = spiked_output(design, v, r)
            if p is None: missing.append((design, v, r)); continue
            res = subprocess.run(["awk", '$1=="C"{print $2"\t"$3}', p], capture_output=True, text=True).stdout
            for l in res.splitlines():
                rid, t = l.split("\t"); t = int(t)
                for node, S in sets.items():
                    if t in S: hits[(r, rid)].add(f"{design}:{v}:{node}")
        print("done", design, v, flush=True)
    rows = [[r, rid, SAMPLE[r], len(s), ";".join(sorted(s))] for (r, rid), s in sorted(hits.items())]
    tsv("reads_at_target_nodes_all_databases.tsv", ["run", "read", "sample", "n_databases", "databases"], rows)
    tsv("missing_outputs.tsv", ["design", "version", "run"], missing)
    # FASTA of both mates for alignment
    os.makedirs(f"{W3}/spur", exist_ok=True)
    byrun = defaultdict(set)
    for r, rid in hits: byrun[r].add(rid)
    with open(f"{W3}/spur/reads.fa", "w") as fo:
        for r, ids in byrun.items():
            for mate in (1, 2):
                keep = False
                with open(f"{FQ}/{r}_{mate}.fastq") as f:
                    for i, l in enumerate(f):
                        if i % 4 == 0:
                            rid = l[1:].split()[0].split("/")[0]; keep = rid in ids
                            if keep: name = f"{rid}/{mate}"
                        elif i % 4 == 1 and keep: fo.write(f">{name}\n{l}")


def spurious():
    """Merge the alignment results (viral megablast, blastn-short to RSV/MRV, minimap2 to hamster genomes)."""
    S = f"{W3}/spur"
    reads = list(csv.DictReader(open(f"{OUT}/reads_at_target_nodes_all_databases.tsv"), delimiter="\t"))
    vir = defaultdict(list)
    for l in open(f"{S}/viral_megablast.tsv"):
        q, s, pid, length, ev, title = l.rstrip("\n").split("\t")[:6]
        vir[q.split("/")[0]].append((float(pid), int(length), s, title))
    lem = defaultdict(int)
    for l in open(f"{S}/short_vs_targets.tsv"):
        c = l.split("\t"); q = c[0].split("/")[0]
        if float(c[2]) == 100.0: lem[q] = max(lem[q], int(c[3]))
    host = {}
    for g in ("picr", "chok1"):
        best = defaultdict(lambda: (0.0, 0.0, 0))
        for l in open(f"{S}/host_{g}.paf"):
            c = l.split("\t"); q = c[0]; ql = int(c[1]); aln = int(c[10]); m = int(c[9])
            idn, cov = m / aln, (int(c[3]) - int(c[2])) / ql
            nm = next((int(x[5:]) for x in c[12:] if x.startswith("NM:i:")), -1)
            if cov * idn > best[q][0] * best[q][1]: best[q] = (idn, cov, nm)
        host[g] = best
    rows = []
    for r in reads:
        rid = r["read"]
        v = sorted(vir.get(rid, []), key=lambda x: -x[1] * x[0])
        vtop = f"{v[0][2]} {v[0][0]:.1f}% over {v[0][1]} nt ({v[0][3][:60]})" if v else "none"
        hb = []
        for g in ("picr", "chok1"):
            for mate in ("1", "2"):
                x = host[g].get(f"{rid}/{mate}")
                hb.append(f"{100 * x[0]:.1f}/{100 * x[1]:.0f}/{x[2]}" if x else "-")
        rows.append([r["run"], rid, r["sample"], r["n_databases"], vtop, lem.get(rid, 0)] + hb + [r["databases"]])
    tsv("reads_at_target_nodes_alignment.tsv",
        ["run", "read", "sample", "n_databases", "best_viral_megablast_hit", "longest_exact_match_RSV_or_MRV_genome_nt",
         "PICR_mate1_identity/coverage/NM", "PICR_mate2", "CHOK1_mate1", "CHOK1_mate2", "databases"], rows)


# ------------------------------------------------------------------ prebuilt databases at confidence 0.1
def oconf():
    rows = []
    for o in OFFICIAL:
        for r in RUNS14:
            rp = f"{W3}/oconf/{o}/{r}_report.txt"
            if not os.path.exists(rp): rows.append([o, r, SAMPLE[r], "", "", ""]); continue
            par, _ = parents_report(rp); inc = clade_fn(par)
            ids = {"RSV": [], "MRV": []}
            for l in open(f"{W3}/oconf/{o}/{r}_output.txt"):
                if l[0] != "C": continue
                c = l.split("\t", 3); t = int(c[2])
                if inc(t, RSV): ids["RSV"].append(c[1])
                if inc(t, MRV): ids["MRV"].append(c[1])
            rows.append([o, r, SAMPLE[r], len(ids["RSV"]), len(ids["MRV"]), ";".join(ids["RSV"] + ids["MRV"])])
    tsv("prebuilt_confidence0.1_target_reads.tsv", ["index", "run", "sample", "RSV_clade", "MRV_clade", "reads"], rows)


# ------------------------------------------------------------------ classified-read jumps between prebuilt databases
def jumps():
    tot = {}; direct = {}; names = {}
    for o in OFFICIAL:
        c = Counter(); n = 0
        for r in RUNS14:
            for l in open(f"{W2}/official/cls/{o}/{r}_report.txt"):
                x = l.rstrip("\n").split("\t")
                if len(x) < 6: continue
                t = int(x[4]); c[t] += int(x[2]); names.setdefault(t, x[5].strip())
                if t == 0: n_uncl = int(x[1])
        direct[o] = c; tot[o] = sum(v for t, v in c.items() if t != 0)
    rows = [[o, tot[o]] for o in OFFICIAL]
    tsv("prebuilt_classified_totals.tsv", ["index", "reads_classified_14_runs"], rows)
    rows = []
    for a, b in zip(OFFICIAL, OFFICIAL[1:]):
        diff = Counter({t: direct[b][t] - direct[a][t] for t in set(direct[a]) | set(direct[b]) if t != 0})
        for t, d in sorted(diff.items(), key=lambda x: -abs(x[1]))[:8]:
            rows.append([f"{a}->{b}", tot[b] - tot[a], t, names.get(t, ""), direct[a][t], direct[b][t], d])
    tsv("prebuilt_classified_jumps_top_taxa.tsv", ["update", "change_in_classified", "taxid", "name", "direct_reads_before", "direct_reads_after", "difference"], rows)


# ------------------------------------------------------------------ vesivirus 2117
def vesi():
    rows = []
    for design in ("d3", "d1", "fixed"):
        for m in MSL:
            dest = Counter()
            for l in open(f"{W2}/panel/cls/{design}/{m}_output.txt"):
                if l.startswith("C\tVesivirus2117__"): dest[int(l.split("\t", 3)[2])] += 1
            nm = names_of(f"{K}/kraken2_dbs/{m}/taxonomy/names.dmp", set(dest))
            rows.append([design, m, sum(dest.values()), "; ".join(f"{t} {nm.get(t, '?')}: {c}" for t, c in dest.most_common())])
    for o in OFFICIAL:
        dest = Counter()
        for l in open(f"{W2}/official/cls/{o}/panel_output.txt"):
            if l.startswith("C\tVesivirus2117__"): dest[int(l.split("\t", 3)[2])] += 1
        _, nm = parents_report(f"{W2}/official/cls/{o}/panel_report.txt")
        rows.append(["prebuilt", o, sum(dest.values()), "; ".join(f"{t} {nm.get(t, '?')}: {c}" for t, c in dest.most_common())])
    tsv("vesivirus2117_destinations.tsv", ["scheme", "version", "classified", "destinations"], rows)
    # the references carrying the destination taxa, and their Design 3 labels per version
    want = set()
    for r in rows:
        for part in r[3].split("; "):
            if part: want.add(int(part.split()[0]))
    lab41 = d3_labels("msl41")
    refs = [a for a, (cur, L) in lab41.items() if cur in want or L in want]
    out = []
    for m in MSL:
        lab = d3_labels(m); par = parents_nodes(tax_path("d1", m))
        nm = names_of(f"{K}/kraken2_dbs/{m}/taxonomy/names.dmp", set(lab[a][1] for a in refs) | set(lab[a][0] for a in refs))
        for a in refs:
            cur, L = lab[a]
            out.append([a, m, cur, cur in par, L, nm.get(L, "?"), "current" if cur == L else ("ancestor" if _is_anc(L, cur, par) else "merged_predecessor")])
    tsv("vesivirus2117_reference_labels.tsv", ["accession", "version", "current_taxid", "current_taxid_in_snapshot", "design3_label", "label_name", "route"], out)


# ------------------------------------------------------------------ build determinism over the 14 runs (minor)
def determinism():
    rows = []
    extra = f"{W3}/extra"; os.makedirs(extra, exist_ok=True)
    pairs = (("MSL41 Design 1 vs rebuilt with 2 threads", f"{K}/kraken2_dbs/msl41", "d1", "msl41", "c_msl41_rebuild_2threads"),
             ("Current vs rebuilt with 2 threads", f"{K}/kraken2_dbs/current", "current", "current", "d_current_rebuild_2threads"))
    for label, db, design, v, fac in pairs:
        for r in RUNS14:
            a = spiked_output(design, v, r)
            if a is None:
                a = f"{extra}/{design}_{r}_output.txt"
                if not os.path.exists(a):
                    subprocess.run([KR, "--db", db, "--memory-mapping", "--paired", "--threads", "6", "--output", a, "--report", a.replace("_output", "_report"),
                                    f"{FQ}/{r}_1.fastq", f"{FQ}/{r}_2.fastq"], capture_output=True)
            b = f"{W2}/factorial/cls/{fac}/{r}_output.txt"
            n = diff_t = diff_line = 0
            with open(a) as fa, open(b) as fb:
                for la, lb in zip(fa, fb):
                    n += 1
                    if la != lb:
                        diff_line += 1
                        if la.split("\t", 3)[2] != lb.split("\t", 3)[2]: diff_t += 1
            rows.append([label, r, SAMPLE[r], n, diff_t, diff_line])
    tsv("build_determinism_14_runs.tsv", ["comparison", "run", "sample", "reads", "reads_with_different_taxid", "reads_with_any_output_difference"], rows)


# ------------------------------------------------------------------ Design 2 / Design 3 classified totals per version (minor)
def denominators():
    rows = []
    runs9 = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195", "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]
    for design in ("fixed", "d3"):
        prev = None
        for m in MSL:
            ids = set()
            for r in runs9:
                p = spiked_output(design, m, r)
                for l in open(p):
                    if l[0] == "C":
                        c = l.split("\t", 3)
                        if c[2] != "0": ids.add(c[1])
            rows.append([design, m, len(ids), len(ids - prev) if prev is not None else "", len(prev - ids) if prev is not None else ""])
            prev = ids
    tsv("classified_reads_per_version_9runs.tsv", ["design", "version", "classified_reads", "newly_classified_vs_previous", "no_longer_classified_vs_previous"], rows)


if __name__ == "__main__":
    for p in sys.argv[1:]:
        try:
            globals()[p]()
        except Exception as e:
            import traceback; traceback.print_exc()
            print("FAILED", p, type(e).__name__, e, flush=True)
