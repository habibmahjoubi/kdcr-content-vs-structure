#!/usr/bin/env python3
"""M5: joint permutation null for max|rho| across the family of delta_K_local candidate
predictors tested against RSV's rank-matched (pre-fixed-library) classification-change
response in Section 5.3, replacing the prose-only 'exploratory, uncorrected' framing with
a quantitative family-wise null.

Three predictors have an actually-computed Spearman rho against RSV at n=12 in the
manuscript: the aggregate/species-level delta_K (P1's primary test), the descendant-drift
candidate (delta_K_local_formal), and the NCBI-snapshot-based candidate (delta_K_ncbi_viral).
REO is excluded: its response has zero variance in this dataset (rank-matched), so no
Spearman correlation is definable for it, exactly as reported in the manuscript.
"""
import csv
import numpy as np
from scipy import stats

BASE = "/home/hm/kdcr"

def col(path, colname):
    with open(f"{BASE}/{path}") as f:
        return [float(r[colname]) for r in csv.DictReader(f)]

rsv_response = col("results/classification_rank_matched.csv", "abs_change_RSV")
pred_aggregate = col("results/classification_rank_matched.csv", "delta_K_species")
pred_descendant = col("results/delta_K_local_formal.csv", "RSV_delta_K_local")
pred_ncbi_viral = col("results/delta_K_ncbi_vs_ictv_merged.csv", "delta_K_ncbi_viral")

predictors = {
    "aggregate_delta_K (P1 primary)": pred_aggregate,
    "descendant_drift_candidate": pred_descendant,
    "delta_K_ncbi_viral": pred_ncbi_viral,
}

print(f"n = {len(rsv_response)}")
print(f"RSV response vector: {rsv_response}")

observed = {}
for name, pred in predictors.items():
    rho, p = stats.spearmanr(pred, rsv_response)
    observed[name] = rho
    print(f"{name}: rho={rho:.4f} p={p:.4g} (uncorrected)")

obs_max_abs_rho = max(abs(v) for v in observed.values())
print(f"\nObserved max|rho| across the 3-test family: {obs_max_abs_rho:.4f} "
      f"({[k for k,v in observed.items() if abs(v)==obs_max_abs_rho][0]})")

rng = np.random.default_rng(20260808)
B = 200000
pred_matrix = np.array(list(predictors.values()))  # 3 x 12
n = len(rsv_response)
resp = np.array(rsv_response)

null_max_abs_rho = np.empty(B)
idx = np.arange(n)
for b in range(B):
    perm = rng.permutation(idx)
    resp_perm = resp[perm]
    rhos = []
    for p in pred_matrix:
        rho, _ = stats.spearmanr(p, resp_perm)
        rhos.append(0.0 if np.isnan(rho) else rho)
    null_max_abs_rho[b] = max(abs(r) for r in rhos)

pval_joint = np.mean(null_max_abs_rho >= obs_max_abs_rho)
print(f"\nPermutation-based joint (family-wise) null of max|rho| over {B} permutations:")
print(f"  mean={null_max_abs_rho.mean():.4f} sd={null_max_abs_rho.std():.4f}")
print(f"  95th pct={np.percentile(null_max_abs_rho,95):.4f}  99th pct={np.percentile(null_max_abs_rho,99):.4f}")
print(f"  P(null max|rho| >= observed {obs_max_abs_rho:.4f}) = {pval_joint:.4f}")
