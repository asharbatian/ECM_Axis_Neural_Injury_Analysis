import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys, os; sys.path.insert(0, os.path.dirname(__file__)); import run_damp_diseases as D
D2 = pd.read_csv("geo_validation/results/damp_D2_coupling.csv"); D1 = pd.read_csv("geo_validation/results/damp_D1_module_effects.csv")
fig, ax = plt.subplots(1, 2, figsize=(13, 6.2), gridspec_kw={"width_ratios": [1.1, 1]})
# A) module effect z per cohort
mods = ["FRAG","SENSE","DAMP","SYNTH","ADAPT","LEUK"]; H = D1.pivot(index="cohort", columns="module", values="z_perm")[mods]
a = ax[0]; im = a.imshow(H.values, cmap="RdBu_r", vmin=-8, vmax=8, aspect="auto"); a.set_xticks(range(len(mods))); a.set_xticklabels(mods); a.set_yticks(range(len(H))); a.set_yticklabels(H.index, fontsize=8)
for i in range(H.shape[0]):
    for j in range(H.shape[1]): v = H.values[i, j]; a.text(j, i, "" if np.isnan(v) else f"{v:.1f}", ha="center", va="center", fontsize=7)
a.set_title("A  Disease vs control: module effect (permutation z)"); plt.colorbar(im, ax=a, shrink=.6)
# B) forest: FRAG-DAMP partial r vs null
s = D2[D2.module == "FRAG"].reset_index(drop=True); y = np.arange(len(s))[::-1]; b = ax[1]
b.errorbar(s.null_mean, y, xerr=1.96*s.null_sd, fmt="none", ecolor="#bbbbbb", lw=6, capsize=0, label="expression-matched null (mean ± 1.96 SD)")
b.scatter(s.partial_r, y, c=["#1B9E77" if p < .05 else "#555555" for p in s.p_up], zorder=3, label="observed FRAG~DAMP | leukocytes (green p<0.05)")
b.set_yticks(y); b.set_yticklabels(s.cohort, fontsize=8); b.axvline(0, c="k", lw=.6); b.set_xlabel("partial Spearman r (within-group residuals)")
pm = s[~s.cohort.str.startswith("RA_")]  # primary set: RA shares controls with OA, excluded from pooling
meta = D.dl_meta(pm.partial_r.values, pm.n.values); b.axvline(meta["r_re"], ls="--", c="#B2182B"); b.set_title(f"B  FRAG~DAMP, pooled r={meta['r_re']:.2f} [{meta['ci_lo']:.2f},{meta['ci_hi']:.2f}]")
b.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1, frameon=False); plt.tight_layout(); plt.savefig("geo_validation/results/fig4_damp_coupling.png", dpi=200)
