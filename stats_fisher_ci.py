import pandas as pd
import math
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2

df = pd.read_csv("per_codon_table.tsv", sep="\t")

# 2x2 counts
# rows: connecting / nonconnecting
# cols: slow / not slow
a = int(df[(df["region"]=="connecting") & (df["is_slow"]==1)].shape[0])
b = int(df[(df["region"]=="connecting") & (df["is_slow"]==0)].shape[0])
c = int(df[(df["region"]=="nonconnecting") & (df["is_slow"]==1)].shape[0])
d = int(df[(df["region"]=="nonconnecting") & (df["is_slow"]==0)].shape[0])

table = [[a, b],
         [c, d]]

# Fisher exact:
# two-sided
oddsratio_two, p_two = fisher_exact(table, alternative="two-sided")
# one-sided greater: H1 = connecting has MORE slow codons (OR > 1)
oddsratio_greater, p_greater = fisher_exact(table, alternative="greater")
# one-sided less: H1 = connecting has FEWER slow codons (OR < 1)
oddsratio_less, p_less = fisher_exact(table, alternative="less")

# 95% CI for OR (Woolf / log OR approx via statsmodels)
t = Table2x2(table)
or_hat = t.oddsratio
ci_low, ci_high = t.oddsratio_confint(alpha=0.05, method="normal")  # log OR normal approx

# proportions
pct_conn = 100 * a / (a+b)
pct_non = 100 * c / (c+d)

out = []
out.append("=== 2x2 Table (rows: region, cols: slow/not slow) ===")
out.append(f"Connecting:    slow={a}  not_slow={b}  (%slow={pct_conn:.2f}%)")
out.append(f"Non-connecting slow={c}  not_slow={d}  (%slow={pct_non:.2f}%)")
out.append("")
out.append(f"Odds Ratio (Table2x2) = {or_hat:.6f}")
out.append(f"95% CI for OR        = [{ci_low:.6f}, {ci_high:.6f}]")
out.append("")
out.append("=== Fisher's Exact Test p-values ===")
out.append(f"Two-sided  p = {p_two:.6g}")
out.append(f"One-sided (greater; H1 OR>1) p = {p_greater:.6g}")
out.append(f"One-sided (less;    H1 OR<1) p = {p_less:.6g}")

text = "\n".join(out)
print(text)

with open("stats_fisher_ci.txt", "w", encoding="utf-8") as f:
    f.write(text + "\n")

print("\nSaved: stats_fisher_ci.txt")
