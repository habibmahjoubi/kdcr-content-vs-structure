#!/usr/bin/env python3
"""A1 panel test: does a genuine taxonomic STRUCTURAL edit (reparenting), at CONSTANT
library/reference composition (sequence present throughout), ever produce a measurable
Kraken2 classification change? Uses the same MSL39->MSL40 transition as the existing
56-taxon positive-control (creation-event) panel, so the two panels are a direct,
same-transition contrast: creation-event effect vs. reparenting-event effect.

Protocol mirrors the positive-control panel exactly (wgsim 5000 paired reads/accession,
150bp, 1% error, kraken2 --paired against the fixed-library before/after DBs), but reads
for all candidates are batched into two combined FASTQ pairs (one submitted against MSL39,
one against MSL40) so each database's hash table is loaded only once instead of once per
candidate -- the classification cost is then 2 kraken2 invocations, not ~2N.
"""
import csv
import os
import re
import subprocess
import sys

BASE = "/home/hm/kdcr"
FIXED_LIB = os.path.join(BASE, "kraken2_dbs_fixed", "library_fixed.fna")
WORKDIR = os.path.join(BASE, "reparenting_panel")
TRANS_BEFORE = "msl39"
TRANS_AFTER = "msl40"

def load_nodes(nodes_dmp_path):
    d = {}
    with open(nodes_dmp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            if len(parts) < 2:
                continue
            taxid_s, parent_s = parts[0].strip(), parts[1].strip()
            if taxid_s.isdigit() and parent_s.isdigit():
                d[int(taxid_s)] = int(parent_s)
    return d

def step1_find_candidates():
    fixed_taxids = set()
    with open(os.path.join(BASE, "kraken2_dbs_fixed", "fixed_taxids.txt")) as f:
        for line in f:
            line = line.strip()
            if line.isdigit():
                fixed_taxids.add(int(line))

    old = load_nodes(os.path.join(BASE, "kraken2_dbs_fixed", TRANS_BEFORE, "taxonomy", "nodes.dmp"))
    new = load_nodes(os.path.join(BASE, "kraken2_dbs_fixed", TRANS_AFTER, "taxonomy", "nodes.dmp"))
    exclude = {11250, 351073}  # RSV, REO -- already tested elsewhere
    reparented = sorted(t for t in fixed_taxids if t not in exclude
                         and t in old and t in new and old[t] != new[t])
    print(f"Reparented fixed-library taxa, {TRANS_BEFORE}->{TRANS_AFTER}: {len(reparented)}")
    return set(reparented)

def step2_extract_sequences(wanted_taxids):
    os.makedirs(WORKDIR, exist_ok=True)
    header_re = re.compile(r"^>(\S+)\|kraken:taxid\|(\d+)\s*(.*)$")
    manifest = []  # (accession, taxid, description)
    seqs = {}
    cur_acc = None
    cur_lines = []
    def flush():
        if cur_acc is not None:
            seqs[cur_acc] = cur_lines[:]
    with open(FIXED_LIB, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                flush()
                cur_lines = []
                m = header_re.match(line.strip())
                if m and int(m.group(2)) in wanted_taxids:
                    cur_acc = m.group(1)
                    manifest.append((m.group(1), int(m.group(2)), m.group(3)))
                else:
                    cur_acc = None
            else:
                if cur_acc is not None:
                    cur_lines.append(line)
        flush()
    print(f"Sequences extracted for {len(manifest)} accessions "
          f"(covering {len(set(t for _, t, _ in manifest))} distinct taxids)")

    fasta_dir = os.path.join(WORKDIR, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    for acc, taxid, desc in manifest:
        with open(os.path.join(fasta_dir, f"{acc}.fna"), "w") as out:
            out.write(f">{acc}\n")
            out.writelines(seqs[acc])

    manifest_csv = os.path.join(WORKDIR, "manifest.csv")
    with open(manifest_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["accession", "taxid", "description"])
        w.writerows(manifest)
    print(f"Saved manifest: {manifest_csv}")
    return manifest

def step3_simulate_and_combine(manifest):
    fasta_dir = os.path.join(WORKDIR, "fasta")
    combined_1 = os.path.join(WORKDIR, "combined_1.fastq")
    combined_2 = os.path.join(WORKDIR, "combined_2.fastq")
    n_ok = 0
    with open(combined_1, "w") as out1, open(combined_2, "w") as out2:
        for acc, taxid, desc in manifest:
            fna = os.path.join(fasta_dir, f"{acc}.fna")
            r1 = os.path.join(WORKDIR, f"{acc}_1.fastq")
            r2 = os.path.join(WORKDIR, f"{acc}_2.fastq")
            with open(os.devnull, "w") as devnull:
                rc = subprocess.call(
                    ["wgsim", "-N", "5000", "-1", "150", "-2", "150", "-e", "0.01", "-r", "0",
                     fna, r1, r2],
                    stdout=devnull, stderr=devnull)
            if rc != 0 or not os.path.exists(r1) or os.path.getsize(r1) == 0:
                print(f"  WARNING: wgsim failed for {acc} (taxid {taxid}), skipping")
                continue
            # Re-tag every read ID with the source accession so we can trace it back
            # after batched classification (e.g. @accession__originalreadid/1)
            for src, out in ((r1, out1), (r2, out2)):
                with open(src) as f:
                    for i, line in enumerate(f):
                        if i % 4 == 0:
                            out.write(f"@{acc}__{line[1:]}")
                        else:
                            out.write(line)
            os.remove(r1)
            os.remove(r2)
            n_ok += 1
    print(f"Simulated + tagged reads for {n_ok}/{len(manifest)} accessions")
    return n_ok

def step4_classify():
    for db, label in ((TRANS_BEFORE, "before"), (TRANS_AFTER, "after")):
        out_txt = os.path.join(WORKDIR, f"combined_output_{label}.txt")
        report_txt = os.path.join(WORKDIR, f"combined_report_{label}.txt")
        cmd = ["kraken2", "--db", os.path.join(BASE, "kraken2_dbs_fixed", db),
               "--paired", "--threads", "4",
               "--report", report_txt, "--output", out_txt,
               os.path.join(WORKDIR, "combined_1.fastq"), os.path.join(WORKDIR, "combined_2.fastq")]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

def step5_score(manifest):
    acc2taxid = {acc: taxid for acc, taxid, _ in manifest}
    taxid2acc = {}
    for acc, taxid, _ in manifest:
        taxid2acc.setdefault(taxid, []).append(acc)

    def parse_output(path):
        # per-accession: [n_reads, n_classified, n_correct]
        stats = {acc: [0, 0, 0] for acc in acc2taxid}
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                status, read_id, assigned = parts[0], parts[1], parts[2]
                if "__" not in read_id:
                    continue
                acc = read_id.split("__", 1)[0]
                if acc not in stats:
                    continue
                stats[acc][0] += 1
                if status == "C":
                    stats[acc][1] += 1
                    try:
                        assigned_id = int(assigned)
                    except ValueError:
                        assigned_id = None
                    if assigned_id == acc2taxid[acc]:
                        stats[acc][2] += 1
        return stats

    before = parse_output(os.path.join(WORKDIR, "combined_output_before.txt"))
    after = parse_output(os.path.join(WORKDIR, "combined_output_after.txt"))

    rows = []
    for acc, taxid, desc in manifest:
        nb, cb, kb = before.get(acc, [0, 0, 0])
        na, ca, ka = after.get(acc, [0, 0, 0])
        rows.append({
            "accession": acc, "taxid": taxid, "description": desc,
            "n_reads_before": nb, "pct_before_any": round(100 * cb / nb, 2) if nb else None,
            "pct_before_correct": round(100 * kb / nb, 2) if nb else None,
            "n_reads_after": na, "pct_after_any": round(100 * ca / na, 2) if na else None,
            "pct_after_correct": round(100 * ka / na, 2) if na else None,
        })

    out_csv = os.path.join(BASE, "results_fixed", "reparenting_panel_MSL39_MSL40.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Saved: {out_csv}")

    valid = [r for r in rows if r["pct_before_correct"] is not None and r["pct_after_correct"] is not None]
    deltas = [r["pct_after_correct"] - r["pct_before_correct"] for r in valid]
    print(f"\nn={len(valid)} accessions scored")
    print(f"mean pct_before_correct = {sum(r['pct_before_correct'] for r in valid)/len(valid):.3f}")
    print(f"mean pct_after_correct  = {sum(r['pct_after_correct'] for r in valid)/len(valid):.3f}")
    print(f"mean delta (after-before) = {sum(deltas)/len(deltas):.3f}")
    print(f"max |delta| = {max(abs(d) for d in deltas):.3f}")
    n_moved = sum(1 for d in deltas if abs(d) > 1.0)
    print(f"accessions with |delta| > 1 percentage point: {n_moved}/{len(valid)}")

if __name__ == "__main__":
    cands = step1_find_candidates()
    manifest = step2_extract_sequences(cands)
    step3_simulate_and_combine(manifest)
    step4_classify()
    step5_score(manifest)
