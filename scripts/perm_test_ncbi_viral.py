import csv
import numpy as np
from scipy import stats

rows = list(csv.DictReader(open("/home/hm/kdcr/results/delta_K_ncbi_vs_ictv_merged.csv")))
dk = np.array([float(r["delta_K_ncbi_viral"]) for r in rows])
rsv = np.array([float(r["abs_change_RSV"]) for r in rows])

rho_obs, p_asym = stats.spearmanr(dk, rsv)
print(f"Observed: rho={rho_obs:.4f}, asymptotic p={p_asym:.4f}")

rng = np.random.default_rng(12345)
n_perm = 200000
count = 0
for _ in range(n_perm):
    perm = rng.permutation(rsv)
    rho, _ = stats.spearmanr(dk, perm)
    if abs(rho) >= abs(rho_obs) - 1e-12:
        count += 1
p_perm = count / n_perm
print(f"Permutation p (two-sided, {n_perm} resamples): {p_perm:.5f}")
