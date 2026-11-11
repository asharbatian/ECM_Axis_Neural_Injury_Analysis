#!/usr/bin/env python3
"""HA-fragmentation machinery vs innate DAMP response across diseases. Implements PREREGISTRATION_damp_diseases.md."""
import os, sys, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
import geo_io as g
import run_geo_validation as R
g.C = os.path.join("geo_validation", "cache"); OUT = "geo_validation/results"
RNG = np.random.default_rng(555); NNULL = 2000
MOD = {
 "FRAG": ["HYAL1","HYAL2","CEMIP","TMEM2"], "SYNTH": ["HAS1","HAS2","HAS3"], "SENSE": ["CD44","TLR2","TLR4","CD14","HMMR"],
 "DAMP": ["MYD88","IRAK4","TRAF6","NFKB1","NFKBIA","TNFAIP3","TNF","IL1B","IL1A","IL6","CXCL8","CCL2","CXCL1","CXCL2","CXCL3","CCL20","PTGS2","NLRP3","ICAM1","SOD2"],
 "LEUK": ["PTPRC","LCP1","CORO1A","CD53","LAPTM5","TYROBP","FCER1G"], "ADAPT": ["CD3D","CD3E","CD2","LCK","CD8A","CD19","MS4A1","CD79A"]}
ALLMOD = set(sum(MOD.values(), []))

def ch(m, i): return [x.split(": ", 1)[-1].strip().replace("\xa0", "") if ": " in x else x.strip() for x in m["Sample_characteristics_ch1"][i]]

def load(acc, gpl):
    fn = g.list_matrix(acc)[0]; e, m = g.series_matrix(acc, fn)
    return R.gene_matrix(e, gpl), m, m["Sample_geo_accession"][0]

def mk(X, ids, grp, **kw):
    d = pd.DataFrame({"grp": grp}, index=ids); d = d[d.grp.isin(["case", "ctrl"])]
    for k, v in kw.items(): d[k] = pd.Series(v, index=ids).loc[d.index]
    return X[d.index], d

def c_GSE53845():
    X, m, ids = load("GSE53845", "GPL6480"); gr = ["case" if "pulmonary" in s else "ctrl" for s in ch(m, 0)]; return mk(X, ids, gr)
def c_GSE48452():
    X, m, ids = load("GSE48452", "GPL11532"); gr = [{"Nash": "case", "Control": "ctrl"}.get(s, "x") for s in ch(m, 1)]; return mk(X, ids, gr)
def c_GSE55235():
    X, m, ids = load("GSE55235", "GPL96"); s = [x.lower() for x in ch(m, 0)]
    return X, m, ids, s
def c_GSE13355():
    X, m, ids = load("GSE13355", "GPL570"); gr = ["case" if s.startswith("involved") else "ctrl" if s.startswith("normal") else "x" for s in ch(m, 0)]; return mk(X, ids, gr)
def c_GSE43292():
    X, m, ids = load("GSE43292", "GPL6244"); gr = ["case" if "Atheroma" in s else "ctrl" for s in ch(m, 2)]; return mk(X, ids, gr, pair=ch(m, 0))
def c_GSE5281():
    X, m, ids = load("GSE5281", "GPL570"); gr = ["case" if "alzheimer" in s.lower() else "ctrl" if "normal" in s.lower() else "x" for s in ch(m, 8)]
    return mk(X, ids, gr, region=[s.lower() for s in ch(m, 4)])
def c_GSE13507():
    X, m, ids = load("GSE13507", "GPL6102"); gr = ["case" if "primary bladder cancer" in s else "ctrl" if "normal bladder" in s else "x" for s in ch(m, 0)]; return mk(X, ids, gr)

def effects(X, d):
    """gene-level effect case vs ctrl (variance-moderated), handling pairs / region."""
    X = X.copy()
    if "region" in d: X = X.sub(X.T.groupby(d["region"].values).transform("mean").T)
    a, b = d.index[d.grp == "case"], d.index[d.grp == "ctrl"]
    if "pair" in d:
        pa = d.loc[a, "pair"]; pb = d.loc[b, "pair"]; common = sorted(set(pa) & set(pb))
        A = X[[pa[pa == p].index[0] for p in common]].values; B = X[[pb[pb == p].index[0] for p in common]].values
        diff = A - B; sd = diff.std(1, ddof=1); s0 = np.nanmedian(sd)
        t, p = stats.ttest_1samp(diff, 0, axis=1)
        return pd.DataFrame({"g": diff.mean(1) / (sd + s0), "p": p}, index=X.index), len(common)
    return R.hedges_g(X, list(a), list(b)), (len(a), len(b))

def zrows(X): return X.sub(X.mean(1), axis=0).div(X.std(1, ddof=1).replace(0, np.nan), axis=0)
def spear(x, y): return stats.spearmanr(x, y)[0]
def pcor(x, y, z):
    rxy, rxz, ryz = spear(x, y), spear(x, z), spear(y, z); return (rxy - rxz*ryz) / np.sqrt(max(1e-9, (1-rxz**2)*(1-ryz**2)))

def resid(Z, d):
    Z = Z.copy()
    if "pair" in d: Z = Z.sub(Z.T.groupby(d["pair"].values).transform("mean").T)
    if "region" in d: Z = Z.sub(Z.T.groupby(d["region"].values).transform("mean").T)
    return Z.sub(Z.T.groupby(d["grp"].values).transform("mean").T)

def coupling(X, d, mu):
    Zr = resid(zrows(X), d); Zraw = zrows(X); out = []
    sc = lambda Z, gs: Z.loc[[x for x in gs if x in Z.index]].mean(0).values
    damp, leuk = sc(Zr, MOD["DAMP"]), sc(Zr, MOD["LEUK"])
    ok = Zr.dropna().index
    bins = pd.qcut(mu.loc[ok], 10, labels=False, duplicates="drop"); pool = {b: ok[(bins.values == b)] for b in np.unique(bins)}
    excl = ALLMOD
    for mod in ["FRAG", "SENSE", "SYNTH", "ADAPT"]:
        pres = [x for x in MOD[mod] if x in X.index]
        if len(pres) < 2: continue
        obs = pcor(sc(Zr, pres), damp, leuk); raw = spear(sc(Zraw, pres), sc(Zraw, MOD["DAMP"]))
        nb = bins.reindex(pres).dropna().astype(int).values; null = []
        for _ in range(NNULL):
            pick = []
            for b in nb:
                c = pool[b][~pool[b].isin(excl)]; pick.append(RNG.choice(c))
            null.append(pcor(Zr.loc[pick].mean(0).values, damp, leuk))
        null = np.array(null)
        out.append(dict(module=mod, n_genes=len(pres), n=X.shape[1], partial_r=obs, raw_r=raw, null_mean=null.mean(), null_sd=null.std(ddof=1),
                        z_null=(obs-null.mean())/null.std(ddof=1), p_up=(1+np.sum(null >= obs))/(1+NNULL)))
    return pd.DataFrame(out)

def dl_meta(r, n):
    z = np.arctanh(np.clip(r, -.99, .99)); v = 1/(np.asarray(n)-4); w = 1/v; zf = (w*z).sum()/w.sum(); Q = (w*(z-zf)**2).sum(); k = len(z)
    tau2 = max(0, (Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); ws = 1/(v+tau2); zr = (ws*z).sum()/ws.sum(); se = np.sqrt(1/ws.sum())
    return dict(k=k, r_re=np.tanh(zr), ci_lo=np.tanh(zr-1.96*se), ci_hi=np.tanh(zr+1.96*se), p=2*stats.norm.sf(abs(zr/se)), I2=max(0, (Q-(k-1))/Q)*100 if Q > 0 else 0)

def main():
    cohorts = {}
    for name, fn in [("IPF_lung_GSE53845", c_GSE53845), ("NASH_liver_GSE48452", c_GSE48452), ("psoriasis_skin_GSE13355", c_GSE13355),
                     ("atheroma_GSE43292", c_GSE43292), ("Alzheimer_brain_GSE5281", c_GSE5281), ("bladder_cancer_GSE13507", c_GSE13507)]:
        try: cohorts[name] = fn(); print(f"{name}: genes={cohorts[name][0].shape[0]} groups={cohorts[name][1].grp.value_counts().to_dict()}")
        except Exception as ex: print(f"DROPPED {name}: {type(ex).__name__} {ex}")
    # synovium: OA and RA vs healthy
    try:
        X, m, ids, s = c_GSE55235()
        for nm, key in [("OA_synovium_GSE55235", "osteo"), ("RA_synovium_GSE55235", "rheum")]:
            gr = ["case" if key in v else "ctrl" if "healthy" in v else "x" for v in s]; cohorts[nm] = mk(X, ids, gr); print(nm, cohorts[nm][1].grp.value_counts().to_dict())
    except Exception as ex: print("DROPPED synovium", ex)
    # injury cohorts (rodent)
    for nm, fn, a, b in [("SCI_rat_GSE45006_D3", R.cohort_GSE45006, "D3", "sham"), ("SCI_mouse_GSE47681_D3", R.cohort_GSE47681, "D3", "sham"), ("TBI_rat_GSE45997", R.cohort_GSE45997, "TBI_ipsi", "naive")]:
        X, grp, _, _ = fn(); gr = pd.Series({k: ("case" if v == a else "ctrl" if v == b else "x") for k, v in grp.items()})
        cohorts[nm] = mk(X, list(gr.index), list(gr.values))
    D1, D2, cov = [], [], {}
    for name, (X, d) in cohorts.items():
        mu = X.mean(1); gt, n = effects(X, d); gt["q"] = R.bh(gt["p"].values)
        t = R.competitive_test(gt, mu, {k: [x for x in v] for k, v in MOD.items()}); t = t.rename(columns={"axis": "module"}); t.insert(0, "cohort", name); t["n"] = str(n); D1.append(t)
        c = coupling(X, d, mu); c.insert(0, "cohort", name); D2.append(c)
        cov[name] = {k: f"{sum(x in X.index for x in v)}/{len(v)}" for k, v in MOD.items()}
        h = t.set_index("module")
        print(f"{name:26s} n={n}  g: FRAG {h.loc['FRAG','mean_g']:+.2f}(z{h.loc['FRAG','z_perm']:+.1f})  SENSE {h.loc['SENSE','mean_g']:+.2f}(z{h.loc['SENSE','z_perm']:+.1f})  DAMP {h.loc['DAMP','mean_g']:+.2f}(z{h.loc['DAMP','z_perm']:+.1f})")
    D1 = pd.concat(D1); D2 = pd.concat(D2); D1.to_csv(f"{OUT}/damp_D1_module_effects.csv", index=False); D2.to_csv(f"{OUT}/damp_D2_coupling.csv", index=False)
    pd.DataFrame(cov).T.to_csv(f"{OUT}/damp_gene_coverage.csv")
    pd.set_option("display.width", 220)
    print("\nD2 partial r (module~DAMP | LEUK, residualised) and z vs expression-matched null:"); print(D2.pivot(index="cohort", columns="module", values="partial_r").round(2)); print(D2.pivot(index="cohort", columns="module", values="z_null").round(1)); print(D2.pivot(index="cohort", columns="module", values="p_up").round(3))
    primary = ["IPF_lung_GSE53845", "NASH_liver_GSE48452", "OA_synovium_GSE55235", "atheroma_GSE43292", "psoriasis_skin_GSE13355", "Alzheimer_brain_GSE5281", "bladder_cancer_GSE13507", "SCI_rat_GSE45006_D3", "SCI_mouse_GSE47681_D3", "TBI_rat_GSE45997"]
    primary = [c for c in primary if c in cohorts]; rows = []
    for mod in ["FRAG", "SENSE", "DAMP", "SYNTH", "ADAPT", "LEUK"]:
        z = D1[(D1.module == mod) & D1.cohort.isin(primary)].z_perm.dropna().values; zs = z.sum()/np.sqrt(len(z))
        rows.append(dict(analysis="D1_effect", module=mod, k=len(z), stouffer_z=zs, p_one_sided=stats.norm.sf(zs), cohorts_positive=int((z > 0).sum())))
    for mod in ["FRAG", "SENSE", "SYNTH", "ADAPT"]:
        s = D2[(D2.module == mod) & D2.cohort.isin(primary)]; r = dl_meta(s.partial_r.values, s.n.values)
        zn = s.z_null.values; rows.append(dict(analysis="D2_coupling_FisherRE", module=mod, **r, stouffer_z_vs_null=zn.sum()/np.sqrt(len(zn)), cohorts_p_lt_0_05=int((s.p_up < .05).sum())))
    M = pd.DataFrame(rows); M.to_csv(f"{OUT}/damp_D3_meta.csv", index=False); print("\nMETA:\n", M.round(3).to_string(index=False))
    e = D1[D1.cohort.isin(primary)].pivot(index="cohort", columns="module", values="mean_g")
    print("\nD5 ecological Spearman (module effect, across primary contrasts):", {m: round(stats.spearmanr(e["FRAG"], e[m])[0], 2) for m in ["DAMP", "SENSE", "SYNTH", "ADAPT"]})
    e.to_csv(f"{OUT}/damp_D5_effects.csv")
if __name__ == "__main__": main()
