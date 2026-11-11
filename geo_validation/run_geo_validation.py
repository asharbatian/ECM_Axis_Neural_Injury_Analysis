#!/usr/bin/env python3
"""
Independent public-GEO validation of the ECM-axis hypothesis in neural injury.

Design (pre-specified, see README_geo_validation.md):
  * ECM axis gene sets: copied unchanged from scripts/02_annotation_mapping.R
  * Per cohort/contrast: gene-level variance-moderated effect size g (injured vs control) on log2 expression
  * Axis statistic: mean gene-level g over axis genes
  * Competitive permutation test: observed statistic vs 10,000 random gene sets of the
    same size, matched on mean-expression decile (guards against global shifts / low-expr bias)
  * Cross-cohort meta-analysis: Stouffer on the PRE-SPECIFIED primary contrasts only
Run:  python3 geo_validation/run_geo_validation.py   (from repo root)
"""
import os, sys, json, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
import geo_io as g
g.C = os.path.join("geo_validation", "cache")      # raw downloads (git-ignored)
os.makedirs(g.C, exist_ok=True)
OUT = os.path.join("geo_validation", "results"); os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(20261111)
NPERM = 10000

AXES = {
 "Hyaluronan": ["Has1","Has2","Has3","Hyal1","Hyal2","Hyal3","Cemip","Tmem2","Spam1","Cd44","Hmmr","Tlr2","Tlr4","Cd14"],
 "Provisional_Matrix": ["Fn1","Eda","Tnc","Itgav","Itgb1","Itgb3","Itga5","Spp1","Thbs1","Thbs2"],
 "PNN_CSPG": ["Acan","Vcan","Ncan","Bcan","Cspg4","Cspg5","Tnr","Hapln1","Hapln4","Ptprz1"],
 "Basement_Membrane": ["Lama1","Lama2","Lama4","Lama5","Lamb1","Lamb2","Lamc1","Col4a1","Col4a2","Col4a3","Nid1","Nid2","Hspg2","Agrn"],
 "Proteases_Regulators": ["Mmp2","Mmp3","Mmp9","Mmp12","Mmp14","Adamts1","Adamts4","Adamts5","Adamts9","Timp1","Timp2","Timp3"],
 "Crosslinking_Fibrosis": ["Lox","Loxl1","Loxl2","Loxl3","Loxl4","Tgm1","Tgm2","Col1a1","Col1a2","Col3a1","Col5a1","Col6a1","Acta2","Tagln"],
}
HA_SUB = {  # sub-modules of the HA axis (same genes as above, just split)
 "HA_synthesis": ["Has1","Has2","Has3"],
 "HA_degradation": ["Hyal1","Hyal2","Hyal3","Cemip","Tmem2","Spam1"],
 "HA_receptor_DAMP": ["Cd44","Hmmr","Tlr2","Tlr4","Cd14"],
}
AXES_UP = {k: [x.upper() for x in v] for k, v in {**AXES, **HA_SUB}.items()}

def gene_matrix(expr, gpl):
    raw = g.platform(gpl)
    import re as _re
    sc = next((c for c in raw.columns if _re.fullmatch(r"(?i)gene[ _]?symbol|symbol", c)), None)
    if sc is None and "gene_assignment" in raw.columns:          # e.g. Affymetrix ST arrays
        raw["Gene symbol"] = raw["gene_assignment"].dropna().map(
            lambda v: "///".join(x.split("//")[1].strip() for x in v.split("///") if len(x.split("//")) > 1))
        sc = "Gene symbol"
    if sc is None: raise RuntimeError(f"no gene-symbol column in {gpl}: {list(raw.columns)[:12]}")
    pl = raw[["ID", sc]].rename(columns={sc: "Gene symbol"}).dropna()
    pl["Gene symbol"] = pl["Gene symbol"].str.split("///")
    pl = pl.explode("Gene symbol"); pl["Gene symbol"] = pl["Gene symbol"].str.strip().str.upper()
    pl = pl[pl["Gene symbol"] != ""]
    e = expr.copy()
    if np.nanmax(e.values) > 40: e = np.log2(e.clip(lower=1))      # raw-scale safeguard
    e = e.dropna(how="all"); e.index = e.index.astype(str)
    m = e.mean(axis=1)
    pl = pl[pl["ID"].isin(e.index)].assign(mu=lambda d: d["ID"].map(m))
    best = pl.sort_values("mu", ascending=False).drop_duplicates("Gene symbol")
    out = e.loc[best["ID"]]; out.index = best["Gene symbol"].values
    return out

def hedges_g(X, a, b):
    """gene-wise variance-moderated Hedges-type effect size g (a vs b) + Welch p. X: genes x samples (DataFrame)."""
    A, B = X[a].values, X[b].values
    na, nb = A.shape[1], B.shape[1]
    sp = np.sqrt(((na-1)*A.var(1, ddof=1) + (nb-1)*B.var(1, ddof=1)) / (na+nb-2))
    # SAM-style variance moderation: add s0 = median pooled SD, so near-zero-variance
    # genes cannot inflate the effect size (essential with n=3-4 per group).
    s0 = np.nanmedian(sp)
    d = (A.mean(1) - B.mean(1)) / (sp + s0)
    J = 1 - 3/(4*(na+nb)-9)
    t, p = stats.ttest_ind(A, B, axis=1, equal_var=False)
    return pd.DataFrame({"logFC": A.mean(1)-B.mean(1), "g": d*J, "p": p}, index=X.index)

def bh(p):
    p = np.asarray(p, float); n = np.isfinite(p).sum(); o = np.argsort(np.where(np.isfinite(p), p, 2))
    q = np.full_like(p, np.nan); ranked = p[o]*n/np.arange(1, len(p)+1)
    ranked = np.minimum.accumulate(np.where(np.isfinite(ranked), ranked, np.inf)[::-1])[::-1]
    q[o] = np.minimum(ranked, 1); q[~np.isfinite(p)] = np.nan; return q

def competitive_test(gtab, mu, sets):
    ok = gtab["g"].notna()
    gg = gtab.loc[ok, "g"]; bins = pd.qcut(mu.loc[gg.index], 10, labels=False, duplicates="drop")
    pool = {b: gg.index[bins.values == b].values for b in np.unique(bins)}
    res = []
    for name, genes in sets.items():
        present = [x for x in genes if x in gg.index]
        if len(present) < 2: res.append(dict(axis=name, n_genes=len(present))); continue
        obs = gg.loc[present].mean()
        nb = bins.loc[present].values
        null = np.empty(NPERM)
        for i in range(NPERM):
            pick = [RNG.choice(pool[b]) for b in nb]
            null[i] = gg.loc[pick].mean()
        z = (obs - null.mean()) / null.std(ddof=1)
        res.append(dict(axis=name, n_genes=len(present), n_total=len(genes), mean_g=obs, z_perm=z,
                        p_up=(1+np.sum(null >= obs))/(NPERM+1), p_two=min(1, 2*min((1+np.sum(null >= obs))/(NPERM+1), (1+np.sum(null <= obs))/(NPERM+1))),
                        genes_present=",".join(present)))
    return pd.DataFrame(res)

# ------------------------------------------------------------------ cohorts
def cohort_GSE45006():
    e, m = g.series_matrix("GSE45006"); t = m["Sample_title"][0]; gsm = m["Sample_geo_accession"][0]
    grp = {}
    for s, ti in zip(gsm, t):
        if "sham" in ti: grp[s] = "sham"
        else:
            tp = ti.split("injured,")[1].split(",")[0].strip()
            grp[s] = {"1 day": "D1", "3 days": "D3", "1 week": "D7", "2 weeks": "D14", "8 weeks": "D56"}[tp]
    return gene_matrix(e, "GPL1355"), grp, [(k, "sham") for k in ["D1","D3","D7","D14","D56"]], dict(organism="rat", tissue="spinal cord", model="clip impact-compression", qc="")

def cohort_GSE47681():
    e, m = g.series_matrix("GSE47681")
    gsm = m["Sample_geo_accession"][0]; title = m["Sample_title"][0]
    ch = {k: v for k, v in m.items() if k == "Sample_characteristics_ch1"}["Sample_characteristics_ch1"]
    geno = ch[2]; tim = ch[5]; grp = {}; dropped = []
    for s, ti, ge, tm in zip(gsm, title, geno, tim):
        if "WT" not in ge: continue
        title_time = "Sham" if "SHAM" in ti else ti.split("WT")[0].replace("DAY", "D")
        meta_time = {"control": "Sham", "Day 1": "D1", "Day 3": "D3", "Day 7": "D7"}[tm.replace("time: ", "")]
        if title_time != meta_time: dropped.append(f"{s}({ti}|{tm})"); continue   # label conflict -> exclude
        grp[s] = "sham" if meta_time == "Sham" else meta_time
    return gene_matrix(e, "GPL1261"), grp, [("D1","sham"),("D3","sham"),("D7","sham")], dict(organism="mouse", tissue="thoracic spinal cord", model="contusion, WT only", qc="excluded label-conflict samples: " + "; ".join(dropped))

def cohort_GSE45997():
    e, m = g.series_matrix("GSE45997"); gsm = m["Sample_geo_accession"][0]; tr = m["Sample_characteristics_ch1"][4]
    grp = {}
    for s, t in zip(gsm, tr):
        grp[s] = "TBI_ipsi" if t.endswith("traumatic brain injury") and "contra" not in t else ("contra" if "contra" in t else "naive")
    return gene_matrix(e, "GPL1355"), grp, [("TBI_ipsi", "naive"), ("TBI_ipsi", "contra")], dict(organism="rat", tissue="brain (ipsilateral)", model="unilateral TBI", qc="")

COHORTS = {"GSE45006": cohort_GSE45006, "GSE47681": cohort_GSE47681, "GSE45997": cohort_GSE45997}
PRIMARY = [("GSE45006","D3","sham"), ("GSE47681","D3","sham"), ("GSE45997","TBI_ipsi","naive")]   # pre-specified

def main():
    allres, genetabs, info = [], [], {}
    for acc, fn in COHORTS.items():
        X, grp, contrasts, meta = fn(); info[acc] = dict(meta, n_genes=int(X.shape[0]), n_samples_used=len(grp),
                                                         groups=pd.Series(grp).value_counts().to_dict())
        print(f"\n{acc}: {X.shape[0]} genes, groups {info[acc]['groups']}  {meta['qc']}")
        mu = X.mean(axis=1)
        for a, b in contrasts:
            sa = [s for s, v in grp.items() if v == a]; sb = [s for s, v in grp.items() if v == b]
            gt = hedges_g(X, sa, sb); gt["q"] = bh(gt["p"].values)
            res = competitive_test(gt, mu, {**AXES_UP})
            res.insert(0, "contrast", f"{a}_vs_{b}"); res.insert(0, "cohort", acc)
            res["n_a"], res["n_b"] = len(sa), len(sb)
            allres.append(res)
            sel = gt.loc[gt.index.isin(sum(AXES_UP.values(), []))].copy()
            sel["axis"] = [next(k for k in AXES if i in AXES_UP[k]) for i in sel.index]
            sel.insert(0, "contrast", f"{a}_vs_{b}"); sel.insert(0, "cohort", acc); genetabs.append(sel.reset_index(names="gene"))
            h = res[res.axis == "Hyaluronan"].iloc[0]
            print(f"  {a:>8} vs {b:<5} n={len(sa)}v{len(sb)}  HA mean_g={h.mean_g:+.2f} z={h.z_perm:+.2f} p_up={h.p_up:.4f}")
    R = pd.concat(allres, ignore_index=True)
    R["q_axes_6"] = np.nan
    main_axes = R.axis.isin(AXES)
    for (c, k), idx in R[main_axes].groupby(["cohort", "contrast"]).groups.items():
        R.loc[idx, "q_axes_6"] = bh(R.loc[idx, "p_up"].values)
    R.to_csv(f"{OUT}/axis_enrichment_all_contrasts.csv", index=False)
    pd.concat(genetabs).to_csv(f"{OUT}/axis_gene_level_results.csv", index=False)
    json.dump(info, open(f"{OUT}/cohort_info.json", "w"), indent=2, default=str)
    # ---- pre-specified meta-analysis
    meta_rows = []
    P = R[[ (r.cohort, r.contrast.split("_vs_")[0], r.contrast.split("_vs_")[1]) in PRIMARY for r in R.itertuples()]]
    for ax, d in P.groupby("axis"):
        z = d["z_perm"].dropna().values
        zs = z.sum()/np.sqrt(len(z)); meta_rows.append(dict(axis=ax, k_cohorts=len(z), stouffer_z=zs, p_one_sided=stats.norm.sf(zs),
                                                              mean_g=d["mean_g"].mean(), per_cohort_z=";".join(f"{x:.2f}" for x in z)))
    M = pd.DataFrame(meta_rows); M["q_BH"] = bh(M["p_one_sided"].values); M = M.sort_values("stouffer_z", ascending=False)
    M.to_csv(f"{OUT}/meta_primary_contrasts.csv", index=False); print("\nPRE-SPECIFIED META (k=3):\n", M.round(4).to_string(index=False))

if __name__ == "__main__": main()
