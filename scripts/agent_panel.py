#!/usr/bin/env python3
"""Panel of adventitious agents relevant to CHO-derived biologics, simulated from their RefSeq genomes
(vesivirus 2117 from GenBank AY343325.2, absent from RefSeq). 5,000 read pairs per agent, wgsim,
150 bp, 1% substitution error, no mutation (-r 0 -R 0), insert 300 +/- 30 bp, fixed seeds."""
import os, sys, subprocess

K = "/home/hm/kdcr"
R2 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
W = f"{K}/work_designs/panel"
SRC = f"{K}/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna"
# agent: (list of accessions, or taxid to take all sequences of that taxid)
AGENTS = {
    "RSV_A":        {"acc": ["NC_038235.1"]},
    "MRV3":         {"taxid": "538123"},
    "MVM":          {"acc": ["NC_001510.1"]},
    "PPV":          {"acc": ["NC_001718.1"]},
    "PCV1":         {"acc": ["NC_001792.2"]},
    "PCV2":         {"acc": ["NC_005148.1"]},
    "EHDV":         {"taxid": "449133"},
    "CacheValley":  {"taxid": "80935"},
    "BVDV1":        {"acc": ["NC_001461.1"]},
    "MoMuLV":       {"acc": ["NC_001501.1"]},
    "LCMV":         {"acc": ["NC_004294.1", "NC_004291.1"]},
    "HAdV_C":       {"acc": ["NC_001405.1"]},
    "SV40":         {"acc": ["NC_001669.1"]},
    "Vesivirus2117": {"file": os.path.join(R2, "..", "sequences", "vesivirus2117_AY343325.2.fna")},
}


def read_fasta(path):
    seqs, cur = {}, None
    with open(path) as f:
        for l in f:
            if l[0] == ">":
                h = l[1:].split()[0]; cur = h; seqs[cur] = []
            else:
                seqs[cur].append(l.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def simulate():
    os.makedirs(W, exist_ok=True)
    lib = read_fasta(SRC)                       # keys: acc|kraken:taxid|tid
    by_acc = {k.split("|")[0]: (k.split("|")[2], s) for k, s in lib.items()}
    manifest = []
    o1, o2 = open(f"{W}/panel_1.fq", "w"), open(f"{W}/panel_2.fq", "w")
    for i, (name, spec) in enumerate(AGENTS.items()):
        if "file" in spec:
            recs = read_fasta(spec["file"]).items()
        elif "taxid" in spec:
            recs = [(a, s) for a, (t, s) in by_acc.items() if t == spec["taxid"]]
        else:
            recs = []
            for a in spec["acc"]:
                if a in by_acc: recs.append((a, by_acc[a][1]))
                else: print("missing", name, a, flush=True)
        fa = f"{W}/{name}.fa"
        with open(fa, "w") as f:
            for a, s in recs: f.write(f">{a}\n{s}\n")
        subprocess.run(f"wgsim -N 5000 -1 150 -2 150 -e 0.01 -r 0 -R 0 -d 300 -s 30 -S {201 + i} "
                       f"{fa} {W}/t1.fq {W}/t2.fq > /dev/null 2>&1", shell=True, check=True)
        for src, dst in ((f"{W}/t1.fq", o1), (f"{W}/t2.fq", o2)):
            with open(src) as f:
                for j, l in enumerate(f):
                    if j % 4 == 0: l = f"@{name}__{l[1:]}"
                    dst.write(l)
        manifest.append((name, ";".join(a for a, _ in recs), sum(len(s) for _, s in recs), 201 + i))
    o1.close(); o2.close()
    with open(f"{R2}/panel_manifest.tsv", "w") as f:
        f.write("agent\taccessions\ttotal_length\twgsim_seed\n")
        for r in manifest: f.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "simulate": simulate()
