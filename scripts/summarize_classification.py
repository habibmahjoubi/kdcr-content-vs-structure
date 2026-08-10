#!/usr/bin/env python3
"""Summarize RSV/REO read counts from Kraken2 reports for the 9 KDCR case-study runs,
robust to a taxon being entirely absent from a given report (unlike the shell version)."""
import glob
import os

REPORT_DIR = "/home/hm/kdcr/results/current_db"

MANIFEST = {
    "SRR26352217": ("REO", "neg"),
    "SRR26352207": ("RSV", "neg"),
    "SRR26352216": ("REO", "1e6"),
    "SRR26352206": ("RSV", "1e6"),
    "SRR26352195": ("REO", "1e5"),
    "SRR26352212": ("REO", "1e3"),
    "SRR26352205": ("RSV", "1e3"),
    "SRR26352208": ("REO", "1e1"),
    "SRR26352202": ("RSV", "1e1"),
}

def count_matching(report_path, keywords):
    total = 0
    with open(report_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            name = cols[5].strip().lower()
            if any(k in name for k in keywords):
                total += int(cols[2])  # column 3 = reads assigned directly at this taxon
    return total

def main():
    print(f"{'Run':<14}{'Virus':<6}{'Dilution':<10}{'RSV_reads':<12}{'REO_reads':<12}")
    rows = []
    for run, (virus, dilution) in MANIFEST.items():
        path = os.path.join(REPORT_DIR, f"{run}_report.txt")
        rsv = count_matching(path, ["syncytial"])
        reo = count_matching(path, ["orthoreovirus", "reoviridae"])
        rows.append((run, virus, dilution, rsv, reo))
        print(f"{run:<14}{virus:<6}{dilution:<10}{rsv:<12}{reo:<12}")

    out_csv = "/home/hm/kdcr/results/rsv_reo_summary_current_db.csv"
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("run,virus,dilution,rsv_reads,reo_reads\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    print(f"\nSaved: {out_csv}")

if __name__ == "__main__":
    main()
