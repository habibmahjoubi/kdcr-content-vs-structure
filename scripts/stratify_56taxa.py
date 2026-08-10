#!/usr/bin/env python3
"""B5: rebuild the 56-taxon misattribution-risk analysis correctly (the previous
'_extended' CSV had a corrupted pct_before_any column, and its generating script no
longer exists), add taxonomic-group (phage vs non-phage) stratification, and test the
coverage-depth confound explicitly (5,000 reads simulated regardless of genome length,
so effective per-base coverage varies ~30-fold across the panel)."""
import csv
import os
import re

BASE = "/home/hm/kdcr"
FASTA_DIR = os.path.join(BASE, "positive_control_56")

def genome_length(fna_path):
    total = 0
    with open(fna_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                total += len(line.strip())
    return total

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

# --- Load base 56-taxon data (correct pct columns) ---
base_rows = []
with open(os.path.join(BASE, "results", "positive_control_56taxa.csv"), encoding="utf-8") as f:
    for rec in csv.DictReader(f):
        base_rows.append(rec)
print(f"Loaded {len(base_rows)} base rows")

# --- Genome length, computed directly from the actual simulated-from FASTA ---
for r in base_rows:
    fna = os.path.join(FASTA_DIR, f"{r['accession']}.fna")
    r["genome_length"] = genome_length(fna) if os.path.exists(fna) else None

missing = [r["accession"] for r in base_rows if r["genome_length"] is None]
if missing:
    print("WARNING: missing FASTA for:", missing)

# --- Sibling count: taxa already present under this taxon's (current) parent at MSL39 ---
# Use MSL40's tree to find each new taxon's parent (since it doesn't exist in MSL39's tree),
# then count how many taxa were already children of that same parent in MSL39.
nodes_40 = load_nodes(os.path.join(BASE, "kraken2_dbs", "msl40", "taxonomy", "nodes.dmp"))
nodes_39 = load_nodes(os.path.join(BASE, "kraken2_dbs", "msl39", "taxonomy", "nodes.dmp"))
children_39 = {}
for taxid, parent in nodes_39.items():
    children_39.setdefault(parent, set()).add(taxid)

for r in base_rows:
    taxid = int(r["taxid"])
    parent = nodes_40.get(taxid)
    r["sibling_count_msl39"] = len(children_39.get(parent, set())) if parent is not None else None

# --- Taxonomic group: phage vs non-phage, from the description text (reliable NCBI convention) ---
for r in base_rows:
    desc = r.get("description", "")
    r["group"] = "phage" if re.search(r"\bphage\b", desc, re.IGNORECASE) else "non-phage virus"

n_phage = sum(1 for r in base_rows if r["group"] == "phage")
print(f"Group split: {n_phage} phage, {len(base_rows) - n_phage} non-phage virus")

# --- Effective per-base coverage: 5000 read pairs * 2 * 150bp / genome_length ---
for r in base_rows:
    if r["genome_length"]:
        r["effective_coverage"] = (5000 * 2 * 150) / r["genome_length"]
    else:
        r["effective_coverage"] = None

# --- Save the corrected, annotated CSV ---
out_csv = os.path.join(BASE, "results", "positive_control_56taxa_stratified.csv")
fieldnames = ["accession", "taxid", "description", "group", "genome_length", "effective_coverage",
              "sibling_count_msl39", "pct_before_any", "pct_before_correct", "pct_after_correct"]
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in base_rows:
        w.writerow({k: r.get(k) for k in fieldnames})
print(f"Saved: {out_csv}")

# --- Statistics ---
from scipy import stats
import statistics as st

def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

valid = [r for r in base_rows if r["genome_length"] and to_float(r["pct_before_any"]) is not None
         and r["sibling_count_msl39"] is not None]
print(f"\nn valid = {len(valid)}")

gl = [r["genome_length"] for r in valid]
sib = [r["sibling_count_msl39"] for r in valid]
cov = [r["effective_coverage"] for r in valid]
pba = [to_float(r["pct_before_any"]) for r in valid]
pbc = [to_float(r["pct_before_correct"]) for r in valid]
pac = [to_float(r["pct_after_correct"]) for r in valid]

print("\n=== Pooled (all 56, both groups) ===")
for label, x in [("genome_length", gl), ("effective_coverage", cov), ("sibling_count", sib)]:
    rho, p = stats.spearmanr(x, pba)
    print(f"{label} vs pct_before_any: rho={rho:.4f} p={p:.4g} n={len(x)}")

print("\n=== Stratified by group ===")
for group in ("phage", "non-phage virus"):
    sub = [r for r in valid if r["group"] == group]
    if len(sub) < 5:
        print(f"{group}: n={len(sub)}, too small for a stable correlation")
        continue
    gl_s = [r["genome_length"] for r in sub]
    cov_s = [r["effective_coverage"] for r in sub]
    sib_s = [r["sibling_count_msl39"] for r in sub]
    pba_s = [to_float(r["pct_before_any"]) for r in sub]
    print(f"\n--- {group} (n={len(sub)}) ---")
    for label, x in [("genome_length", gl_s), ("effective_coverage", cov_s), ("sibling_count", sib_s)]:
        rho, p = stats.spearmanr(x, pba_s)
        print(f"{label} vs pct_before_any: rho={rho:.4f} p={p:.4g} n={len(x)}")

# genome_length vs group (is length itself confounded with group membership?)
gl_phage = [r["genome_length"] for r in valid if r["group"] == "phage"]
gl_nonphage = [r["genome_length"] for r in valid if r["group"] == "non-phage virus"]
print(f"\nGenome length by group: phage mean={st.mean(gl_phage):.0f} (n={len(gl_phage)}), "
      f"non-phage mean={st.mean(gl_nonphage):.0f} (n={len(gl_nonphage)})")
try:
    u, p_mw = stats.mannwhitneyu(gl_phage, gl_nonphage)
    print(f"Mann-Whitney U test, genome length by group: p={p_mw:.4g}")
except Exception as e:
    print("Mann-Whitney failed:", e)

# effective_coverage vs genome_length -- should be near-perfectly anti-correlated (mechanical)
rho_cov_gl, p_cov_gl = stats.spearmanr(gl, cov)
print(f"\ngenome_length vs effective_coverage (mechanical check): rho={rho_cov_gl:.4f} p={p_cov_gl:.4g}")

# Does coverage explain MORE or LESS variance than raw length? partial-ish check:
# correlate residual of pba after regressing out coverage, against genome_length, and vice versa
print(f"\ngenome_length range: {min(gl)}-{max(gl)} bp")
print(f"effective_coverage range: {min(cov):.1f}x-{max(cov):.1f}x")

# three 100%-misattributed taxa mentioned in the manuscript
hundred = [r for r in valid if to_float(r["pct_before_any"]) == 100.0]
print(f"\nTaxa with 100% before-state misattribution (n={len(hundred)}):")
for r in hundred:
    print(f"  {r['accession']} {r['taxid']} [{r['group']}] {r['description'][:60]} "
          f"len={r['genome_length']} cov={r['effective_coverage']:.1f}x siblings={r['sibling_count_msl39']}")
