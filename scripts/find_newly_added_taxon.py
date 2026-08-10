#!/usr/bin/env python3
"""Find a viral genome sequence (from the tagged bulk RefSeq FASTA) whose taxid did NOT exist
in an older MSL-year taxonomy snapshot but DOES exist in a more recent one -- a real
"newly recognized taxon" positive control for testing whether KDCR's classification pipeline
correctly detects drift-driven change (contrasting with the RSV/REO stability result)."""
import re
import sys

def load_valid_taxids(nodes_dmp_path):
    valid = set()
    with open(nodes_dmp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            taxid = line.split("\t|\t", 1)[0].strip()
            if taxid.isdigit():
                valid.add(int(taxid))
    return valid

def main():
    old_nodes, new_nodes, fasta = sys.argv[1], sys.argv[2], sys.argv[3]
    old_valid = load_valid_taxids(old_nodes)
    new_valid = load_valid_taxids(new_nodes)
    newly_added = new_valid - old_valid
    print(f"Taxids newly present in {new_nodes} vs {old_nodes}: {len(newly_added)}")

    header_re = re.compile(r"^>(\S+)\|kraken:taxid\|(\d+)\s*(.*)$")
    candidates = []
    with open(fasta, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                m = header_re.match(line.strip())
                if m:
                    acc, taxid, desc = m.group(1), int(m.group(2)), m.group(3)
                    if taxid in newly_added:
                        candidates.append((acc, taxid, desc))

    print(f"Sequences in bulk FASTA with a newly-added taxid: {len(candidates)}")
    for acc, taxid, desc in candidates[:20]:
        print(f"  {acc}\t{taxid}\t{desc}")

if __name__ == "__main__":
    main()
