import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("per_codon_table.tsv", sep="\t")

# Observed 2x2
def compute_or(sub):
    a = int(sub[(sub["region"]=="connecting") & (sub["is_slow"]==1)].shape[0])
    b = int(sub[(sub["region"]=="connecting") & (sub["is_slow"]==0)].shape[0])
    c = int(sub[(sub["region"]=="nonconnecting") & (sub["is_slow"]==1)].shape[0])
    d = int(sub[(sub["region"]=="nonconnecting") & (sub["is_slow"]==0)].shape[0])
    # Haldane-Anscombe correction if any zero
    if min(a,b,c,d) == 0:
        a,b,c,d = a+0.5, b+0.5, c+0.5, d+0.5
    return (a*d)/(b*c)

OR_obs = compute_or(df)

# --- Permutation design ---
# Within each accession:
# keep the number of connecting codons fixed, randomly assign which positions are "connecting"
# This preserves gene-level composition and region size, breaks positional association.

rng = np.random.default_rng(42)
B = 5000  # permutations

ors = np.zeros(B, dtype=float)

# Pre-split data by accession for speed
groups = {acc: g.copy() for acc, g in df.groupby("accession")}

for b in range(B):
    perm_parts = []
    for acc, g in groups.items():
        n = g.shape[0]
        k = int((g["region"]=="connecting").sum())  # number of connecting codons in that gene
        # assign k positions as connecting at random
        idx = np.arange(n)
        rng.shuffle(idx)
        conn_idx = set(idx[:k])
        # create permuted region labels
        reg = np.array(["nonconnecting"]*n, dtype=object)
        reg[list(conn_idx)] = "connecting"
        gg = g.copy()
        gg["region"] = reg
        perm_parts.append(gg)
    perm_df = pd.concat(perm_parts, ignore_index=True)
    ors[b] = compute_or(perm_df)

# One-sided p-value for H1: OR_obs > OR_perm
p_greater = (np.sum(ors >= OR_obs) + 1) / (B + 1)
# One-sided p-value for H1: OR_obs < OR_perm
p_less = (np.sum(ors <= OR_obs) + 1) / (B + 1)
# Two-sided (distance from 1 on log scale)
logdist_obs = abs(np.log(OR_obs))
logdist_perm = abs(np.log(ors))
p_two = (np.sum(logdist_perm >= logdist_obs) + 1) / (B + 1)

print(f"Observed OR = {OR_obs:.6f}")
print(f"Permutation B = {B}")
print(f"Permutation p (greater; OR>expected) = {p_greater:.6g}")
print(f"Permutation p (less; OR<expected)    = {p_less:.6g}")
print(f"Permutation p (two-sided)            = {p_two:.6g}")

with open("permutation_test_or.txt", "w", encoding="utf-8") as f:
    f.write(f"Observed OR = {OR_obs:.6f}\n")
    f.write(f"B = {B}\n")
    f.write(f"Permutation p (greater) = {p_greater:.6g}\n")
    f.write(f"Permutation p (less)    = {p_less:.6g}\n")
    f.write(f"Permutation p (two-sided)= {p_two:.6g}\n")

# Figure: OR distribution
plt.figure()
plt.hist(ors, bins=40, alpha=0.85)
plt.axvline(OR_obs, linewidth=2)
plt.xlabel("Odds Ratio (OR) under permutation")
plt.ylabel("Count")
plt.title("Figure 3. Permutation distribution of OR (within-gene shuffling)")
plt.savefig("Figure3_permutation_OR_distribution.png", dpi=300, bbox_inches="tight")
plt.close()

print("Saved: permutation_test_or.txt, Figure3_permutation_OR_distribution.png")
