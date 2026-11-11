import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
S = pd.read_csv("geo_validation/results/atlas_mouse_GSE10246_module_scores.csv", index_col=0)
CNS = {"amygdala","cerebellum","cerebral_cortex","cerebral_cortex_prefrontal","dorsal_striatum","hippocampus","hypothalamus","nucleus_accumbens","olfactory_bulb","spinal_cord"}
S = S.assign(d=S.Stabilization - S.Fragmentation).sort_values("d", ascending=False)
fig, ax = plt.subplots(figsize=(6.5, 10)); im = ax.imshow(S.drop(columns="d").values, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(S.columns[:4], rotation=30, ha="right"); ax.set_yticks(range(len(S))); ax.set_yticklabels(S.index, fontsize=7)
for t, l in zip(ax.get_yticklabels(), S.index):
    if l in CNS: t.set_color("#B2182B"); t.set_fontweight("bold")
plt.colorbar(im, label="module z (mean gene z across tissues)", shrink=.5); ax.set_title("Mouse tissue atlas (GSE10246)\nCNS tissues in red; sorted by Stabilization - Fragmentation")
plt.tight_layout(); plt.savefig("geo_validation/results/fig3_tissue_atlas_modules.png", dpi=200)
