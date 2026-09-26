#!/usr/bin/env python3
"""Resimulate the 56 positive-control candidates at matched effective coverage (30x)
instead of the original fixed 5,000 read pairs regardless of genome length, decoupling
genome length from per-base coverage (previously perfectly anti-correlated, rho=-1.00)
so the genome-length-vs-misattribution-risk correlation can be tested
without that confound. Classifies against the existing (unmodified, naive-tagged) MSL39
('before') and MSL40 ('after') Kraken2 databases, matching the original panel's design.
"""
import csv
import os
import re
import subprocess

BASE = "/home/hm/kdcr"
FASTA_DIR = f"{BASE}/positive_control_56"
OUTDIR = f"{BASE}/positive_control_56_covmatched"
POS56 = f"{BASE}/results/positive_control_56taxa.csv"
COVERAGE = 30
READ_LEN = 150
DB_BEFORE = f"{BASE}/kraken2_dbs/msl39"
DB_AFTER = f"{BASE}/kraken2_dbs/msl40"

os.makedirs(OUTDIR, exist_ok=True)


def genome_length(fna_path):
    total = 0
    with open(fna_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def classify(db, r1, r2, report, output, taxid):
    subprocess.run(
        ["kraken2", "--db", db, "--paired", "--threads", "4",
         "--report", report, "--output", output, r1, r2],
        check=True, capture_output=True,
    )
    n_total = n_correct = n_any = 0
    with open(output) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            n_total += 1
            status, taxid_field = parts[0], parts[2]
            m = re.search(r"\(taxid (\d+)\)", taxid_field) or re.match(r"^(\d+)$", taxid_field)
            assigned = int(m.group(1)) if m else None
            if status != "U":
                n_any += 1
            if assigned == taxid:
                n_correct += 1
    if n_total == 0:
        return 0, 0.0, 0.0
    return n_total, 100.0 * n_correct / n_total, 100.0 * n_any / n_total


def main():
    rows = list(csv.DictReader(open(POS56)))
    out_rows = []
    for r in rows:
        acc, taxid, desc = r["accession"], int(r["taxid"]), r["description"]
        fna = f"{FASTA_DIR}/{acc}.fna"
        gl = genome_length(fna)
        n_pairs = max(50, round(COVERAGE * gl / (2 * READ_LEN)))

        r1 = f"{OUTDIR}/{acc}_1.fastq"
        r2 = f"{OUTDIR}/{acc}_2.fastq"
        log = f"{OUTDIR}/{acc}_wgsim.log"
        subprocess.run(
            ["wgsim", "-N", str(n_pairs), "-1", str(READ_LEN), "-2", str(READ_LEN),
             "-e", "0.01", "-r", "0", fna, r1, r2],
            check=True, stdout=open(log, "w"), stderr=subprocess.STDOUT,
        )

        n_before, pct_before_correct, pct_before_any = classify(
            DB_BEFORE, r1, r2, f"{OUTDIR}/{acc}_before_report.txt",
            f"{OUTDIR}/{acc}_before_output.txt", taxid)
        n_after, pct_after_correct, pct_after_any = classify(
            DB_AFTER, r1, r2, f"{OUTDIR}/{acc}_after_report.txt",
            f"{OUTDIR}/{acc}_after_output.txt", taxid)

        print(f"{acc} taxid={taxid} gl={gl} n_pairs={n_pairs} cov={COVERAGE}x  "
              f"before: any={pct_before_any:.2f}% correct={pct_before_correct:.2f}%  "
              f"after: correct={pct_after_correct:.2f}%")

        out_rows.append({
            "accession": acc, "taxid": taxid, "description": desc,
            "genome_length": gl, "n_pairs_simulated": n_pairs, "coverage_x": COVERAGE,
            "pct_before_any": round(pct_before_any, 2),
            "pct_before_correct": round(pct_before_correct, 2),
            "pct_after_any": round(pct_after_any, 2),
            "pct_after_correct": round(pct_after_correct, 2),
        })

    out_csv = f"{BASE}/results/positive_control_56taxa_covmatched.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nSaved: {out_csv}")


if __name__ == "__main__":
    main()
