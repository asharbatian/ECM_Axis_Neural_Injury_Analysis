#!/usr/bin/env python3
"""Figures for the GEO validation (run after run_geo_validation.py)."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = pd.read_csv("geo_validation/results/axis_enrichment_all_contrasts.csv")
order = ["Hyaluronan","HA_receptor_DAMP","HA_synthesis","HA_degradation","Provisional_Matrix","Proteases_Regulators",
         "Crosslinking_Fibrosis","Basement_Membrane","PNN_CSPG"]
R["lab"] = R.cohort + " " + R.contrast.str.replace("_vs_", " vs ")
rows = list(dict.fromkeys(R.lab))   # keep chronological order from the analysis script
M = R.pivot(index="lab", columns="axis", values="z_perm").reindex(index=rows, columns=order)
fig, ax = plt.subplots(figsize=(9, 5)); im = ax.imshow(M.values, cmap="RdBu_r", vmin=-8, vmax=8, aspect="auto")
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=40, ha="right"); ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=8)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M.values[i, j]; ax.text(j, i, "" if np.isnan(v) else f"{v:.1f}", ha="center", va="center", fontsize=7)
plt.colorbar(im, label="competitive permutation z (expression-matched null)"); ax.set_title("ECM axis activation in public neural-injury cohorts")
plt.tight_layout(); plt.savefig("geo_validation/results/fig1_axis_heatmap.png", dpi=200); plt.close()
H = R[R.axis == "Hyaluronan"].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(6.5, 4.5)); y = np.arange(len(H))[::-1]
ax.barh(y, H.z_perm, color=["#1B9E77" if p < 0.05 else "#bbbbbb" for p in H.p_up]); ax.set_yticks(y); ax.set_yticklabels(H.lab, fontsize=8)
ax.axvline(1.645, ls="--", c="k", lw=.8); ax.set_xlabel("Hyaluronan-axis permutation z (dashed = one-sided p=0.05)"); ax.set_title("HA axis by cohort / timepoint")
plt.tight_layout(); plt.savefig("geo_validation/results/fig2_HA_axis_by_contrast.png", dpi=200)
