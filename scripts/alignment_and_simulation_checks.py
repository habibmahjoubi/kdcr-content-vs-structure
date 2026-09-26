#!/usr/bin/env python3
"""Alignment and simulation checks; tables written to ../results. Run inside WSL:
    python alignment_and_simulation_checks.py [A B C D E]
A  origin of reads reassigned on MSL30->31 and MSL31->32 (phiX174?)  -> reassigned_reads_origin.tsv, reassigned_top_taxid_pairs.tsv, phix_in_fixed_libraries.tsv
B  target-node reads: BLAST confirmation + cross-negatives  -> target_reads_blast.tsv, cross_negatives.tsv
C  simulated RSV / MRV reads against all 26 (+current) databases  -> simulated_targets_counts.tsv, simulated_targets_read_changes.tsv
D  MRV 10^6 run: current vs MSL41 database, where did the reads go  -> mrv_current_vs_msl41.tsv
E  read-level reassignment with the alternative snapshot pairing  -> reassignment_alternative_pairing.tsv
"""
import os, sys, random, subprocess, csv
from collections import Counter, defaultdict

K = "/home/hm/kdcr"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
RUNS9 = ["SRR26352217", "SRR26352207", "SRR26352216", "SRR26352206", "SRR26352195",
         "SRR26352212", "SRR26352205", "SRR26352208", "SRR26352202"]
VIRUS = {"SRR26352217": "MRV neg", "SRR26352207": "RSV neg", "SRR26352216": "MRV 1e6",
         "SRR26352206": "RSV 1e6", "SRR26352195": "MRV 1e5", "SRR26352212": "MRV 1e3",
         "SRR26352205": "RSV 1e3", "SRR26352208": "MRV 1e1", "SRR26352202": "RSV 1e1",
         "SRR26352204": "RSV 1e5"}
MSL = [f"msl{n}" for n in range(29, 42)]
RSV, MRV = 11250, 351073
PHIX = "NC_001422.1"
BIN = "/home/hm/miniconda3/envs/viral-score/bin"
KRAKEN = "/home/hm/miniforge3/envs/kdcr/bin/kraken2"  # v2.17.1, as used to build and classify in the study
LIB_CURRENT = f"{K}/kraken2_dbs/current/library/added/BsVO5XdUTo.fna"
WORK = f"{K}/work"  # scratch on the Linux filesystem (no spaces)
BLASTDB = os.path.join(WORK, "blastdb", "viral_current")


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


def write_tsv(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(header); w.writerows(rows)
    print("  wrote", name, len(rows), "rows")


# ---------------------------------------------------------------- taxonomy helpers
_tax = {}
def taxonomy(taxdir):
    if taxdir not in _tax:
        parent, rank, name = {}, {}, {}
        with open(os.path.join(taxdir, "nodes.dmp")) as f:
            for l in f:
                p = l.split("\t|\t"); t = int(p[0]); parent[t] = int(p[1]); rank[t] = p[2]
        with open(os.path.join(taxdir, "names.dmp"), errors="replace") as f:
            for l in f:
                p = l.split("\t|\t")
                if "scientific name" in p[3]:
                    name[int(p[0])] = p[1]
        _tax[taxdir] = (parent, rank, name)
    return _tax[taxdir]


def lineage(taxdir, t):
    parent, rank, name = taxonomy(taxdir)
    out = []
    while t in parent and t != 1:
        out.append(f"{name.get(t, '?')} [{t},{rank.get(t)}]")
        if parent[t] == t: break
        t = parent[t]
    return " < ".join(out)


def in_clade(taxdir, t, target):
    parent, _, _ = taxonomy(taxdir)
    while t in parent:
        if t == target: return True
        if parent[t] == t or t == 1: return False
        t = parent[t]
    return False


def load_output(path, classified_only=True):
    d = {}
    with open(path, errors="replace") as f:
        for l in f:
            if classified_only and l[0] != "C": continue
            p = l.split("\t", 3)
            d[p[1]] = int(p[2])
    return d


def fetch_reads(run, ids, fqdir=f"{K}/data/fastq", mates=(1, 2)):
    """return {(read_id, mate): seq} for wanted ids."""
    ids = set(ids); got = {}
    for m in mates:
        with open(f"{fqdir}/{run}_{m}.fastq") as f:
            while True:
                h = f.readline()
                if not h: break
                s = f.readline().strip(); f.readline(); f.readline()
                rid = h[1:].split()[0].split("/")[0]
                if rid in ids: got[(rid, m)] = s
    return got


# ---------------------------------------------------------------- BLAST helpers
def ensure_blastdb():
    if os.path.exists(BLASTDB + ".nsq") or os.path.exists(BLASTDB + ".00.nsq"): return
    os.makedirs(os.path.dirname(BLASTDB), exist_ok=True)
    fa = BLASTDB + ".fna"
    with open(LIB_CURRENT) as fi, open(fa, "w") as fo:
        for l in fi:
            if l[0] == ">":
                desc = l.split(" ", 1)[1].strip() if " " in l else ""
                parts = l[1:].split()[0].split("|")
                fo.write(f">{parts[0]} taxid={parts[2] if len(parts) > 2 else '?'} {desc}\n")
            else:
                fo.write(l)
    sh(f"makeblastdb -in {fa} -dbtype nucl -out {BLASTDB} -title viral_current > /dev/null")


def blast_best(seqs, tag):
    """seqs: {qid: seq}. returns {qid: (sacc, stitle, pident, alen, qlen, evalue)} best bitscore."""
    ensure_blastdb()
    os.makedirs(WORK, exist_ok=True)
    q = os.path.join(WORK, f"tmp_{tag}.fa"); o = os.path.join(WORK, f"tmp_{tag}.b6")
    with open(q, "w") as f:
        for k, s in seqs.items(): f.write(f">{k}\n{s}\n")
    sh(f"blastn -task megablast -query {q} -db {BLASTDB} -evalue 1e-5 -max_target_seqs 5 "
       f"-num_threads 8 -outfmt '6 qseqid sseqid pident length qlen evalue bitscore stitle' > {o}")
    best = {}
    with open(o) as f:
        for l in f:
            p = l.rstrip("\n").split("\t")
            b = float(p[6])
            if p[0] not in best or b > best[p[0]][-1]:
                best[p[0]] = (p[1], p[7], float(p[2]), int(p[3]), int(p[4]), p[5], b)
    os.remove(q); os.remove(o)
    return best


# ================================================================= A
def part_A():
    print("A: origin of reassigned reads")
    tax_b = {"msl31": f"{K}/kraken2_dbs_fixed/msl31/taxonomy", "msl32": f"{K}/kraken2_dbs_fixed/msl32/taxonomy"}
    summary, pairs_rows = [], []
    rng = random.Random(1)
    for a, b in (("msl30", "msl31"), ("msl31", "msl32")):
        disc, allcls = {}, {}
        for run in RUNS9:
            da = load_output(f"{K}/results_fixed/{a}/{run}_output.txt")
            db = load_output(f"{K}/results_fixed/{b}/{run}_output.txt")
            for rid, ta in da.items():
                tb = db.get(rid)
                if tb is None: continue
                allcls[(run, rid)] = (ta, tb)
                if ta != tb: disc[(run, rid)] = (ta, tb)
        # BLAST every discordant read (mate 1) and a random sample of 5,000 concordant classified reads
        conc = [k for k in allcls if k not in disc]
        sample = rng.sample(conc, min(5000, len(conc)))
        for label, keys in (("discordant", list(disc)), ("concordant_sample", sample)):
            seqs = {}
            byrun = defaultdict(list)
            for run, rid in keys: byrun[run].append(rid)
            for run, ids in byrun.items():
                for (rid, m), s in fetch_reads(run, ids, mates=(1,)).items():
                    seqs[f"{run}|{rid}"] = s
            best = blast_best(seqs, f"A_{a}_{label}")
            n = len(keys)
            hit = sum(1 for k in seqs if k in best)
            phix = sum(1 for k in seqs if k in best and best[k][0] == PHIX)
            phix97 = sum(1 for k in seqs if k in best and best[k][0] == PHIX and best[k][2] >= 97)
            top = Counter(best[k][1][:70] for k in seqs if k in best).most_common(5)
            summary.append([f"{a}->{b}", label, n, hit, phix, phix97,
                            round(100 * phix / n, 2) if n else "", round(100 * phix97 / n, 2) if n else "",
                            "; ".join(f"{t} ({c})" for t, c in top)])
            print("  ", summary[-1][:8])
        pc = Counter(disc.values())
        _, _, nm = taxonomy(tax_b[b])
        for (ta, tb), c in pc.most_common(10):
            pairs_rows.append([f"{a}->{b}", ta, nm.get(ta, "?"), tb, nm.get(tb, "?"), c])
    write_tsv("reassigned_reads_origin.tsv",
              ["update", "read_set", "n_reads", "n_with_blast_hit", "best_hit_phiX174",
               "best_hit_phiX174_pid>=97", "pct_phiX174", "pct_phiX174_pid>=97", "top_best_hits"], summary)
    write_tsv("reassigned_top_taxid_pairs.tsv", ["update", "taxid_before", "name_in_newer_tree", "taxid_after",
                                        "name_after", "n_reads"], pairs_rows)
    # is phiX174 in the fixed libraries?
    rows = []
    for lib in (f"{K}/kraken2_dbs_fixed/library_fixed.fna",):
        found = subprocess.run(f"grep -c '^>{PHIX}' {lib}", shell=True, capture_output=True, text=True).stdout.strip()
        rows.append([lib, found])
    for d in ("kraken2_dbs_fixed_corrected",):
        for m in ("msl29", "msl41"):
            p = f"{K}/{d}/{m}/seqid2taxid.map"
            if os.path.exists(p):
                found = subprocess.run(f"grep -c '^{PHIX}' {p}", shell=True, capture_output=True, text=True).stdout.strip()
                rows.append([p, found])
    write_tsv("phix_in_fixed_libraries.tsv", ["library", "count_NC_001422.1"], rows)


# ================================================================= B
def part_B():
    print("B: target-node reads")
    tx = f"{K}/kraken2_dbs_fixed/msl41/taxonomy"
    rows, cross = [], []
    sources = [(run, f"{K}/results_fixed/msl41/{run}_output.txt", "fixed MSL41") for run in RUNS9]
    alt204 = f"{K}/results_altpairing/msl41/SRR26352204_output.txt"
    if os.path.exists(alt204): sources.append(("SRR26352204", alt204, "fixed MSL41, alt pairing"))
    for run, path, dbl in sources:
        d = load_output(path)
        rsv = [r for r, t in d.items() if in_clade(tx, t, RSV)]
        mrv = [r for r, t in d.items() if in_clade(tx, t, MRV)]
        cross.append([run, VIRUS[run], dbl, len(rsv), len(mrv)])
        want = [(r, "RSV node") for r in rsv] + [(r, "MRV node") for r in mrv]
        if not want: continue
        seqs = fetch_reads(run, [r for r, _ in want])
        best = blast_best({f"{rid}|{m}": s for (rid, m), s in seqs.items()}, f"B_{run}")
        for rid, node in want:
            for m in (1, 2):
                b = best.get(f"{rid}|{m}")
                rows.append([run, VIRUS[run], node, rid, d[rid], m,
                             b[0] if b else "no hit", b[1][:80] if b else "", b[2] if b else "",
                             b[3] if b else "", b[4] if b else ""])
    # current database: RSV-named or MRV-node reads in every run (the name-based summary flagged SRR26352217)
    txc = f"{K}/kraken2_dbs/current/taxonomy"
    for run in RUNS9:
        p = f"{K}/results/current_db/{run}_output.txt"
        if not os.path.exists(p): continue
        d = load_output(p)
        rsv = [r for r, t in d.items() if in_clade(txc, t, RSV)]
        mrv = [r for r, t in d.items() if in_clade(txc, t, MRV)]
        _, _, nm = taxonomy(txc)
        syn = [r for r, t in d.items() if "syncytial" in nm.get(t, "").lower() and r not in rsv]
        cross.append([run, VIRUS[run], "current (Design 1 library, Aug 2026 tree)", len(rsv), len(mrv)])
        extra = [(r, "RSV node (current)") for r in rsv if VIRUS[run].startswith("MRV")] + \
                [(r, "MRV node (current)") for r in mrv if VIRUS[run].startswith("RSV") or "neg" in VIRUS[run]] + \
                [(r, f"other syncytial name: {nm.get(d[r])}") for r in syn]
        if extra:
            seqs = fetch_reads(run, [r for r, _ in extra])
            best = blast_best({f"{rid}|{m}": s for (rid, m), s in seqs.items()}, f"Bc_{run}")
            for rid, node in extra:
                for m in (1, 2):
                    b = best.get(f"{rid}|{m}")
                    rows.append([run, VIRUS[run], node, rid, d[rid], m,
                                 b[0] if b else "no hit", b[1][:80] if b else "", b[2] if b else "",
                                 b[3] if b else "", b[4] if b else ""])
    write_tsv("cross_negatives.tsv", ["run", "sample", "database", "reads_in_RSV_clade_11250",
                                        "reads_in_MRV_clade_351073"], cross)
    write_tsv("target_reads_blast.tsv", ["run", "sample", "node", "read_id", "kraken_taxid", "mate",
                                           "best_hit", "best_hit_title", "pident", "aln_len", "read_len"], rows)


# ================================================================= C
SIM_TARGETS = {
    "RSV_A_NC_038235": (["NC_038235.1"], RSV),
    "RSV_A_NC_001803": (["NC_001803.1"], RSV),
    "RSV_B_NC_001781": (["NC_001781.1"], RSV),
    "MRV3_538123": ([f"NC_0132{n}.1" for n in range(25, 35)], MRV),
    "MRV3_T3D_10886": ([f"NC_0778{n}.1" for n in range(37, 47)], MRV),
}


def part_C():
    print("C: simulated targets")
    sd = os.path.join(WORK, "sim"); os.makedirs(sd, exist_ok=True)
    r1, r2 = os.path.join(sd, "sim_1.fq"), os.path.join(sd, "sim_2.fq")
    if not os.path.exists(r1):
        seqs, cur = {}, None
        with open(LIB_CURRENT) as f:
            for l in f:
                if l[0] == ">":
                    cur = l[1:].split("|")[0].split()[0]; seqs[cur] = []
                else:
                    seqs[cur].append(l.strip())
        with open(r1, "w") as o1, open(r2, "w") as o2:
            for i, (lab, (accs, _)) in enumerate(SIM_TARGETS.items()):
                fa = os.path.join(sd, f"{lab}.fa")
                with open(fa, "w") as f:
                    for a in accs: f.write(f">{a}\n{''.join(seqs[a])}\n")
                sh(f"wgsim -N 5000 -1 150 -2 150 -e 0.01 -r 0 -R 0 -S {101 + i} {fa} {sd}/t1.fq {sd}/t2.fq > /dev/null")
                for src, dst in ((f"{sd}/t1.fq", o1), (f"{sd}/t2.fq", o2)):
                    with open(src) as f:
                        for j, l in enumerate(f):
                            if j % 4 == 0: l = f"@{lab}__{l[1:]}"
                            dst.write(l)
    dbs = [(f"fixed_{m}", f"{K}/kraken2_dbs_fixed/{m}") for m in MSL] + \
          [(f"design1_{m}", f"{K}/kraken2_dbs/{m}") for m in MSL] +           ([("current", f"{K}/kraken2_dbs/current")] if os.environ.get("WITH_CURRENT") else [])
    counts, perread = [], {}
    for name, db in dbs:
        out = os.path.join(sd, f"{name}.out")
        if not os.path.exists(out):
            sh(f"{KRAKEN} --db {db} --memory-mapping --paired --threads 8 --output {out} {r1} {r2} > /dev/null 2>&1")
        # parent links only, not cached: 27 full taxonomies with names exhaust WSL memory
        par = {}
        with open(os.path.join(db, "taxonomy", "nodes.dmp")) as f:
            for l in f:
                p = l.split("	|	", 2); par[int(p[0])] = int(p[1])
        def clade(t, target):
            while t in par:
                if t == target: return True
                if par[t] == t: return False
                t = par[t]
            return False
        d = load_output(out, classified_only=False)
        tally = defaultdict(lambda: [0, 0, 0])  # in clade, classified elsewhere, unclassified
        for rid, t in d.items():
            lab = rid.split("__")[0]; target = SIM_TARGETS[lab][1]
            if t == 0: tally[lab][2] += 1
            elif clade(t, target): tally[lab][0] += 1
            else: tally[lab][1] += 1
        del par
        perread[name] = d
        for lab in SIM_TARGETS:
            c = tally[lab]; counts.append([name, lab, c[0], c[1], c[2], round(100 * c[0] / 5000, 2)])
    write_tsv("simulated_targets_counts.tsv", ["database", "target", "reads_in_target_clade",
                                                "classified_outside_clade", "unclassified", "pct_in_clade"], counts)
    conc = []
    for design in ("fixed", "design1"):
        for i in range(12):
            a, b = f"{design}_{MSL[i]}", f"{design}_{MSL[i + 1]}"
            for lab in SIM_TARGETS:
                ids = [r for r in perread[a] if r.startswith(lab + "__")]
                changed = sum(1 for r in ids if perread[a][r] != perread[b].get(r))
                conc.append([design, f"{MSL[i]}->{MSL[i + 1]}", lab, len(ids), changed])
    write_tsv("simulated_targets_read_changes.tsv", ["design", "update", "target", "n_read_pairs",
                                                     "n_reads_changing_taxid"], conc)


# ================================================================= D
def part_D():
    print("D: MRV 1e6, current vs MSL41")
    run = "SRR26352216"
    txc, tx41 = f"{K}/kraken2_dbs/current/taxonomy", f"{K}/kraken2_dbs/msl41/taxonomy"
    dc = load_output(f"{K}/results/current_db/{run}_output.txt")
    d41 = load_output(f"{K}/results/msl41/{run}_output.txt")
    ids = [r for r, t in dc.items() if in_clade(txc, t, MRV)] + [r for r, t in d41.items() if in_clade(tx41, t, MRV)]
    rows = []
    for r in sorted(set(ids)):
        tc, t41 = dc.get(r, 0), d41.get(r, 0)
        rows.append([r, tc, lineage(txc, tc)[:300], t41, lineage(tx41, t41)[:300],
                     "same taxid" if tc == t41 else "different taxid"])
    # which library sequences carry the taxids involved, and their lineage in both trees
    involved = {x for row in rows for x in (row[1], row[3])}
    m = defaultdict(list)
    with open(f"{K}/kraken2_dbs/current/seqid2taxid.map") as f:
        for l in f:
            s, t = l.split()[:2]
            if int(t) in involved: m[int(t)].append(s.split("|")[0])
    for t in sorted(involved):
        rows.append(["TAXID", t, lineage(txc, t)[:300], t, lineage(tx41, t)[:300], ",".join(m[t][:12])])
    write_tsv("mrv_current_vs_msl41.tsv", ["read_id", "taxid_current", "lineage_current",
                                             "taxid_msl41", "lineage_msl41", "note"], rows)


# ================================================================= E
def part_E():
    print("E: read-level reassignment, alternative pairing")
    base = f"{K}/results_altpairing"
    runs = sorted({f.split("_")[0] for f in os.listdir(f"{base}/msl30") if f.endswith("_output.txt")})
    lib = {}
    for m in ("msl29", "msl41"):
        p = f"{K}/kraken2_dbs_altpairing/{m}/seqid2taxid.map"
        lib[m] = sum(1 for _ in open(p)) if os.path.exists(p) else "?"
    rows = []
    for i in range(12):
        a, b = MSL[i], MSL[i + 1]
        both = disc = 0
        for run in runs:
            if run not in RUNS9: continue
            da, db = load_output(f"{base}/{a}/{run}_output.txt"), load_output(f"{base}/{b}/{run}_output.txt")
            for r, t in da.items():
                u = db.get(r)
                if u is None: continue
                both += 1; disc += (t != u)
        rows.append([f"{a}->{b}", both, disc, round(100 * disc / both, 3) if both else ""])
        print("  ", rows[-1])
    rows.append([f"library size msl29/msl41: {lib['msl29']}/{lib['msl41']}", "", "", ""])
    write_tsv("reassignment_alternative_pairing.tsv", ["update", "classified_in_both", "discordant", "pct"], rows)


if __name__ == "__main__":
    parts = sys.argv[1:] or list("ABCDE")
    for p in parts:
        {"A": part_A, "B": part_B, "C": part_C, "D": part_D, "E": part_E}[p]()
