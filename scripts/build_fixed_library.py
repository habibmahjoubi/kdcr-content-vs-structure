#!/usr/bin/env python3
import re, sys

BASE = "/home/hm/kdcr"
MSLS = [f"msl{n}" for n in range(29, 42)]
TAGGED_FASTA = f"{BASE}/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna"
HEADER_RE = re.compile(r"^>(\S+)\|kraken:taxid\|(\d+)")

def load_valid_taxids(nodes_dmp_path):
    valid = set()
    with open(nodes_dmp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            taxid = line.split("\t|\t", 1)[0].strip()
            if taxid.isdigit():
                valid.add(int(taxid))
    return valid

def main():
    tagged_taxids = set()
    acc_by_taxid = {}
    with open(TAGGED_FASTA, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                m = HEADER_RE.match(line.strip())
                if m:
                    acc, taxid = m.group(1), int(m.group(2))
                    tagged_taxids.add(taxid)
    print(f"Tagged taxids in bulk FASTA: {len(tagged_taxids)}")

    per_version_valid = {}
    for msl in MSLS:
        nodes_path = f"{BASE}/kraken2_dbs/{msl}/taxonomy/nodes.dmp"
        valid = load_valid_taxids(nodes_path)
        restricted = valid & tagged_taxids
        per_version_valid[msl] = restricted
        print(f"{msl}: {len(restricted)} tagged taxids valid")

    fixed = set(tagged_taxids)
    for msl in MSLS:
        fixed &= per_version_valid[msl]
    print(f"\nFIXED_TAXIDS (13-way intersection): {len(fixed)}")

    msl29_set = per_version_valid["msl29"]
    disappeared = msl29_set - fixed
    print(f"MSL29 count (paper's stated 6653): {len(msl29_set)}")
    print(f"Taxids present at MSL29 but absent from fixed intersection (disappeared later): {len(disappeared)}")
    if disappeared:
        print(f"  examples: {list(disappeared)[:20]}")

    with open(f"{BASE}/kraken2_dbs_fixed/fixed_taxids.txt", "w") as f:
        for t in sorted(fixed):
            f.write(f"{t}\n")
    print(f"\nSaved {len(fixed)} fixed taxids to kraken2_dbs_fixed/fixed_taxids.txt")

if __name__ == "__main__":
    main()
