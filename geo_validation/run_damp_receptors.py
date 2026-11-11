#!/usr/bin/env python3
"""Per-gene breakdown of the pre-registered SENSE module (CD44, CD14, TLR2, TLR4) -- EXPLORATORY decomposition."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import run_damp_diseases as D, run_geo_validation as R
D.NNULL = 1000
GENES = ["CD44", "CD14", "TLR2", "TLR4"]
PRIM = ["IPF_lung_GSE53845","NASH_liver_GSE48452","OA_synovium_GSE55235","atheroma_GSE43292","psoriasis_skin_GSE13355","Alzheimer_brain_GSE5281","bladder_cancer_GSE13507","SCI_rat_GSE45006_D3","SCI_mouse_GSE47681_D3","TBI_rat_GSE45997"]
cohorts = {}
for name, fn in [("IPF_lung_GSE53845", D.c_GSE53845), ("NASH_liver_GSE48452", D.c_GSE48452), ("psoriasis_skin_GSE13355", D.c_GSE13355),
                 ("atheroma_GSE43292", D.c_GSE43292), ("Alzheimer_brain_GSE5281", D.c_GSE5281), ("bladder_cancer_GSE13507", D.c_GSE13507)]:
    cohorts[name] = fn()
X, m, ids, s = D.c_GSE55235(); gr = ["case" if "osteo" in v else "ctrl" if "healthy" in v else "x" for v in s]; cohorts["OA_synovium_GSE55235"] = D.mk(X, ids, gr)
for nm, fn, a, b in [("SCI_rat_GSE45006_D3", R.cohort_GSE45006, "D3", "sham"), ("SCI_mouse_GSE47681_D3", R.cohort_GSE47681, "D3", "sham"), ("TBI_rat_GSE45997", R.cohort_GSE45997, "TBI_ipsi", "naive")]:
    Xr, grp, _, _ = fn(); gg = pd.Series({k: ("case" if v == a else "ctrl" if v == b else "x") for k, v in grp.items()}); cohorts[nm] = D.mk(Xr, list(gg.index), list(gg.values))
rows = []
for name in PRIM:
    X, d = cohorts[name]; mu = X.mean(1); gt, _ = D.effects(X, d); gt["q"] = R.bh(gt["p"].values)
    a, b = d.index[d.grp == "case"], d.index[d.grp == "ctrl"]
    if "pair" in d:   # paired logFC
        pa, pb = d.loc[a, "pair"], d.loc[b, "pair"]; com = sorted(set(pa) & set(pb))
        lfc = (X[[pa[pa == p].index[0] for p in com]].values - X[[pb[pb == p].index[0] for p in com]].values).mean(1); lfc = pd.Series(lfc, index=X.index)
    else: lfc = X[a].mean(1) - X[b].mean(1)
    Zr = D.resid(D.zrows(X), d); damp = Zr.loc[[x for x in D.MOD["DAMP"] if x in Zr.index]].mean(0).values
    leuk = Zr.loc[[x for x in D.MOD["LEUK"] if x in Zr.index]].mean(0).values
    ok = Zr.dropna().index; bins = pd.qcut(mu.loc[ok], 10, labels=False, duplicates="drop")
    for gname in GENES:
        if gname not in X.index: rows.append(dict(cohort=name, gene=gname)); continue
        r = D.pcor(Zr.loc[gname].values, damp, leuk); b0 = bins.get(gname)
        pool = ok[(bins.values == b0)]; pool = pool[~pool.isin(D.ALLMOD)]
        null = np.array([D.pcor(Zr.loc[RNG_pick].values, damp, leuk) for RNG_pick in np.random.default_rng(7).choice(pool, D.NNULL)])
        rows.append(dict(cohort=name, gene=gname, n=X.shape[1], log2FC=lfc[gname], g=gt.loc[gname, "g"], q=gt.loc[gname, "q"], partial_r=r,
                         null_mean=null.mean(), z_null=(r-null.mean())/null.std(ddof=1), p_up=(1+np.sum(null >= r))/(1+len(null))))
T = pd.DataFrame(rows); T.to_csv("geo_validation/results/damp_receptor_genes.csv", index=False)
pd.set_option("display.width", 200)
print("log2FC (disease vs control); * = BH q<0.05"); L = T.pivot(index="cohort", columns="gene", values="log2FC")[GENES].round(2).astype(str)
Q = T.pivot(index="cohort", columns="gene", values="q")[GENES]; print((L + np.where(Q < 0.05, "*", " ")).reindex(PRIM).to_string())
print("\npartial r with DAMP | leukocytes (residualised)"); print(T.pivot(index="cohort", columns="gene", values="partial_r")[GENES].reindex(PRIM).round(2).to_string())
out = []
for gname in GENES:
    s2 = T[(T.gene == gname) & T.partial_r.notna()]; r = D.dl_meta(s2.partial_r.values, s2.n.values); zn = s2.z_null.values
    zz = T[(T.gene == gname) & T.g.notna()]; zs = (zz.g.rank() * 0).sum()
    out.append(dict(gene=gname, k=r["k"], pooled_partial_r=r["r_re"], lo=r["ci_lo"], hi=r["ci_hi"], I2=r["I2"], null_z_stouffer=zn.sum()/np.sqrt(len(zn)),
                    cohorts_up_q05=int(((zz.log2FC > 0) & (zz.q < .05)).sum()), cohorts_down_q05=int(((zz.log2FC < 0) & (zz.q < .05)).sum()), k_eff=len(zz)))
print("\nPOOLED per gene:"); P = pd.DataFrame(out); print(P.round(2).to_string(index=False)); P.to_csv("geo_validation/results/damp_receptor_genes_pooled.csv", index=False)
