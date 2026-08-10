#!/usr/bin/env python3
"""M5 follow-up: the full-permutation null (m5_permutation_joint_null.py) treats the 12
transitions as exchangeable, which ignores that they are cumulative and temporally ordered
(later transitions build on earlier taxonomy states). A circular-shift (cyclic) permutation
null preserves the response vector's internal temporal structure (autocorrelation, run
lengths) while still destroying its alignment with each predictor, giving a null that is
robust to, rather than assuming away, temporal dependence.
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
    "aggregate_delta_K": pred_aggregate,
    "descendant_drift_candidate": pred_descendant,
    "delta_K_ncbi_viral": pred_ncbi_viral,
}

n = len(rsv_response)
resp = np.array(rsv_response)
pred_matrix = np.array(list(predictors.values()))

observed = {}
for name, p in predictors.items():
    rho, _ = stats.spearmanr(p, resp)
    observed[name] = rho
obs_max_abs_rho = max(abs(v) for v in observed.values())
print(f"n={n}, observed max|rho|={obs_max_abs_rho:.4f}")

# All n-1 non-trivial cyclic shifts (shift=0 is the observed data itself, excluded)
cyclic_max_abs_rho = []
for shift in range(1, n):
    resp_shifted = np.roll(resp, shift)
    rhos = []
    for p in pred_matrix:
        rho, _ = stats.spearmanr(p, resp_shifted)
        rhos.append(0.0 if np.isnan(rho) else rho)
    cyclic_max_abs_rho.append(max(abs(r) for r in rhos))
cyclic_max_abs_rho = np.array(cyclic_max_abs_rho)

print(f"Cyclic-shift null (n-1={n-1} non-trivial shifts):")
print(f"  values: {np.round(cyclic_max_abs_rho, 4).tolist()}")
p_cyclic = np.mean(cyclic_max_abs_rho >= obs_max_abs_rho)
print(f"  P(cyclic max|rho| >= observed {obs_max_abs_rho:.4f}) = {p_cyclic:.4f} "
      f"({int(np.sum(cyclic_max_abs_rho >= obs_max_abs_rho))}/{n-1})")

# Block bootstrap alternative: circular block permutation with block size 3 (preserves
# short-range autocorrelation / run structure better than single-point cyclic shift alone)
rng = np.random.default_rng(20260809)
block_size = 3
B = 200000
block_max_abs_rho = np.empty(B)
idx = np.arange(n)
n_blocks = int(np.ceil(n / block_size))
for b in range(B):
    start = rng.integers(0, n)
    perm_idx = [(start + i) % n for i in range(n)]
    # circular block shuffle: cut into blocks of block_size after a random circular start,
    # then randomly permute the block order (Politis-Romano circular block bootstrap, applied
    # as a permutation of block order rather than resampling-with-replacement)
    blocks = [perm_idx[i:i+block_size] for i in range(0, n, block_size)]
    rng.shuffle(blocks)
    shuffled_idx = [i for block in blocks for i in block]
    resp_perm = resp[shuffled_idx]
    rhos = []
    for p in pred_matrix:
        rho, _ = stats.spearmanr(p, resp_perm)
        rhos.append(0.0 if np.isnan(rho) else rho)
    block_max_abs_rho[b] = max(abs(r) for r in rhos)

p_block = np.mean(block_max_abs_rho >= obs_max_abs_rho)
print(f"\nCircular block-permutation null (block size={block_size}, {B} draws):")
print(f"  mean={block_max_abs_rho.mean():.4f} sd={block_max_abs_rho.std():.4f}")
print(f"  P(block max|rho| >= observed {obs_max_abs_rho:.4f}) = {p_block:.4f}")
