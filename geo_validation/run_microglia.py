#!/usr/bin/env python3
"""GSE96055 sorted-microglia test. Implements PREREGISTRATION_microglia.md."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import geo_io as g
import run_geo_validation as R
g.C = os.path.join("geo_validation", "cache"); OUT = "geo_validation/results"
d = pd.read_csv(g.fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96055/suppl/GSE96055_Microglia_normcount_table.txt.gz", "GSE96055_normcount.txt.gz"), sep="\t", index_col=0)
d = d[d.mean(1) >= 5]; X = np.log2(d + 1); X.index = X.index.str.upper(); X = X[~X.index.duplicated()]
mu = X.mean(1)
def grp(c): 
    k = c.split("_")[0]; return "NI" if k == "MGNI" else k[4:]            # 72h / 1w / 2w (injury code stripped)
def inj(c): return c.split("_")[0][2:4]
contrasts = []
for t in ["72h", "1w", "2w"]:
    contrasts.append((f"injured_{t}", [c for c in X.columns if grp(c) == t], [c for c in X.columns if grp(c) == "NI"]))
for t in ["72h", "1w", "2w"]:
    for ty, nm in [("FT", "full"), ("HS", "hemi")]:
        contrasts.append((f"{nm}_{t}", [c for c in X.columns if grp(c) == t and inj(c) == ty], [c for c in X.columns if grp(c) == "NI"]))
res, genes = [], []
for name, a, b in contrasts:
    gt = R.hedges_g(X, a, b); gt["q"] = R.bh(gt["p"].values)
    t = R.competitive_test(gt, mu, R.AXES_UP); t.insert(0, "contrast", name); t["n_inj"], t["n_ctrl"] = len(a), len(b); res.append(t)
    sel = gt.loc[gt.index.isin(sum(R.AXES_UP.values(), []))].reset_index(names="gene"); sel.insert(0, "contrast", name); genes.append(sel)
Rz = pd.concat(res); Rz.to_csv(f"{OUT}/microglia_GSE96055_axis_tests.csv", index=False)
G = pd.concat(genes); G.to_csv(f"{OUT}/microglia_GSE96055_gene_level.csv", index=False)
pd.set_option("display.width", 220)
print(Rz[Rz.contrast.str.startswith("injured")].pivot(index="axis", columns="contrast", values="z_perm").round(2))
print(Rz[Rz.contrast.str.startswith("injured")].pivot(index="axis", columns="contrast", values="p_up").round(4))
print(Rz[Rz.contrast.str.startswith(("full", "hemi"))].pivot(index="axis", columns="contrast", values="z_perm").round(2))
print(G[(G.contrast == "injured_72h") & G.gene.isin(R.AXES_UP["Hyaluronan"])][["gene", "logFC", "q"]].round(3).to_string(index=False))
