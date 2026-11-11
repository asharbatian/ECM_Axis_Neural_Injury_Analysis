#!/usr/bin/env python3
"""Tissue-atlas test of HA-module predictions. Implements PREREGISTRATION_tissue_atlas.md exactly."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import geo_io as g
from run_geo_validation import gene_matrix, bh
g.C = os.path.join("geo_validation", "cache"); OUT = "geo_validation/results"
RNG = np.random.default_rng(1133); NP = 100000
MOD = {"Synthesis": ["HAS1","HAS2","HAS3"], "Stabilization": ["HAPLN1","HAPLN2","HAPLN3","HAPLN4"],
       "Fragmentation": ["HYAL1","HYAL2","HYAL4"], "Clearance": ["CEMIP","TMEM2","SPAM1"]}
GRP = {"CNS": ["amygdala","cerebellum","cerebral_cortex","cerebral_cortex_prefrontal","dorsal_striatum","hippocampus","hypothalamus","nucleus_accumbens","olfactory_bulb","spinal_cord"],
       "GI": ["stomach","intestine_small","intestine_large"], "SpleenLiver": ["spleen","liver"],
       "Adipose": ["adipose_white","adipose_brown"], "Bladder": ["bladder"], "Lung": ["lung"]}
EXCL = ("3T3","Baf3","C2C12","C3H","NK_","RAW","T-cells","B-cells","follicular","dendritic","granul","macrophage","mast_","mega_","min6","neuro2a","nih_3T3",
        "osteo","embryonic","stem_","thymocyte","mIMCD","common_myeloid","microglia","bone_marrow")

def module_scores(X):          # X: genes x tissues (log2)
    Z = X.sub(X.mean(1), axis=0).div(X.std(1, ddof=1), axis=0)
    out, cov = {}, {}
    for m, genes in MOD.items():
        pres = [x for x in genes if x in Z.index]; cov[m] = f"{len(pres)}/{len(genes)}:{','.join(pres)}"
        out[m] = Z.loc[pres].mean(0) if pres else pd.Series(np.nan, index=Z.columns)
    return pd.DataFrame(out), cov

def pct_rank(s):               # 1 = highest
    return s.rank(ascending=False, method="average")

def test(score, members, high=True, label=""):
    sc = score.dropna(); mem = [t for t in members if t in sc.index]
    if not mem: return None
    r = pct_rank(sc); obs = r.loc[mem].mean(); n = len(sc); k = len(mem)
    null = np.array([r.values[RNG.choice(n, k, replace=False)].mean() for _ in range(NP//10)])  # 10k perms
    p = (1+np.sum(null <= obs))/(len(null)+1) if high else (1+np.sum(null >= obs))/(len(null)+1)
    return dict(test=label, n_tissues=n, members=",".join(mem), mean_rank=obs, expected_rank=(n+1)/2, p=p,
                ranks=";".join(f"{t}:{int(r[t])}" for t in mem))

def run(acc, X, groups, name):
    S, cov = module_scores(X); rows = []; print(f"\n=== {name}: {X.shape[1]} tissues; coverage {cov}")
    S.to_csv(f"{OUT}/atlas_{name}_module_scores.csv")
    P = [("P1 CNS high Stabilization","CNS","Stabilization",True),("P2a CNS low Fragmentation","CNS","Fragmentation",False),
         ("P2b CNS low Clearance","CNS","Clearance",False),("P4 SpleenLiver high Fragmentation","SpleenLiver","Fragmentation",True),
         ("P5 GI high Clearance","GI","Clearance",True),("P6 Adipose high Synthesis","Adipose","Synthesis",True)]
    for lab, grp, mod, hi in P:
        if grp in groups:
            r = test(S[mod], groups[grp], hi, lab); 
            if r: r["module"] = mod; rows.append(r)
    if "Bladder" in groups and groups["Bladder"][0] in S.index:
        s = S["Synthesis"].dropna(); rk = int(pct_rank(s)[groups["Bladder"][0]])
        rows.append(dict(test="P3 Bladder top Synthesis", module="Synthesis", n_tissues=len(s), members="bladder", mean_rank=rk, expected_rank=(len(s)+1)/2,
                         p=rk/len(s) if rk>1 else 1/len(s), ranks=f"bladder:{rk}  (top-1 hit={rk==1}; top-3 hit={rk<=3})"))
    if "CNS" in groups:   # P7: stabilization minus fragmentation
        D = (S["Stabilization"] - S["Fragmentation"]); r = test(D, groups["CNS"], True, "P7 CNS high Stabilization-Fragmentation")
        if r: r["module"] = "Stab-Frag"; rows.append(r)
    R = pd.DataFrame(rows); R["q_BH"] = bh(R["p"].values); R.insert(0, "atlas", name)
    print(R[["test","n_tissues","mean_rank","expected_rank","p","q_BH","ranks"]].round(4).to_string(index=False, max_colwidth=70)); return R, S, cov

# ---- mouse
e, m = g.series_matrix("GSE10246"); lab = pd.Series(m["Sample_source_name_ch1"][0], index=m["Sample_geo_accession"][0])
keep = lab[~lab.str.startswith(EXCL) & ~lab.isin(["epidermis_cornea"]) ] if False else lab[[not any(l.startswith(x) for x in EXCL) for l in lab]]
X = gene_matrix(e[keep.index], "GPL1261"); X.columns = keep.values; X = X.T.groupby(level=0).mean().T
# keep only tissues (drop a few non-tissue leftovers)
Rm, Sm, cm = run("GSE10246", X, GRP, "mouse_GSE10246")
# ---- human (CNS / liver / lung only)
e, m = g.series_matrix("GSE1133", "GSE1133-GPL96_series_matrix.txt.gz")
lab = pd.Series([s.split("_", 1)[-1].strip().lower().replace(" ", "") for s in m["Sample_source_name_ch1"][0]], index=m["Sample_geo_accession"][0])
Xh = gene_matrix(e, "GPL96"); Xh.columns = lab.values
HCNS = ["spinalcord","fetalbrain","brain,caudatenucleus","cerebellum","brainthalamus","wholebrain","cerebellumpeduncles","parietallobe","occipitallobe","prefrontalcortex","cingulatecortex","hypothalamus","olfactorybulb","medullaoblongata","pons","temporallobe","brainamygdala","globuspallidus","subthalamicnucleus"]
Xh = Xh.T.groupby(level=0).mean().T
Xh = Xh[[c for c in Xh.columns if not any(k in c for k in ["leukemia","lymphoma","bm_","pb_","hbec","cultured","fetal","placenta","adenocarc","721_"]) or c in HCNS]]
Rh, Sh, ch = run("GSE1133", Xh, {"CNS": HCNS, "SpleenLiver": ["liver"], "Lung": ["lung"]}, "human_GSE1133")
pd.concat([Rm, Rh]).to_csv(f"{OUT}/atlas_prediction_tests.csv", index=False)
