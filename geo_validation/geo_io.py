import requests,gzip,io,os,re,pandas as pd
H="https://ftp.ncbi.nlm.nih.gov/geo/"
C="geo_cache"; os.makedirs(C,exist_ok=True)
def _stub(a):
    pre = "GSE" if a.startswith("GSE") else "GPL"; return pre + a[3:-3] + "nnn"
def fetch(url,fn):
    p=os.path.join(C,fn)
    if not os.path.exists(p):
        r=requests.get(url,timeout=300); r.raise_for_status(); open(p,"wb").write(r.content)
    return p
def list_matrix(a):
    cached = sorted(f for f in os.listdir(C) if f.startswith(a) and f.endswith("series_matrix.txt.gz"))
    if cached: return cached                      # already downloaded: no network needed
    for i in range(4):
        try:
            r=requests.get(f"{H}series/{_stub(a)}/{a}/matrix/",timeout=90); r.raise_for_status()
            return re.findall(r'href="([^"]+series_matrix\.txt\.gz)"',r.text)
        except Exception:
            if i == 3: raise
def series_matrix(a,fname=None):
    fname=fname or list_matrix(a)[0]
    p=fetch(f"{H}series/{_stub(a)}/{a}/matrix/{fname}",fname)
    meta={};rows=[];tab=[]
    with gzip.open(p,"rt",errors="replace") as f:
        intab=False
        for l in f:
            l=l.rstrip("\n")
            if l.startswith("!series_matrix_table_begin"): intab=True;continue
            if l.startswith("!series_matrix_table_end"): break
            if intab: tab.append(l)
            elif l.startswith("!Sample_"):
                k,*v=l.split("\t"); meta.setdefault(k[1:],[]).append([x.strip('"') for x in v])
    expr=pd.read_csv(io.StringIO("\n".join(tab)),sep="\t",index_col=0).apply(pd.to_numeric,errors="coerce")
    expr.columns=[c.strip('"') for c in expr.columns]
    return expr,meta
def pheno(meta,acc_cols):
    df=pd.DataFrame({"gsm":meta["Sample_geo_accession"][0],"title":meta["Sample_title"][0]})
    for i,ch in enumerate(meta.get("Sample_characteristics_ch1",[])):
        df[f"char{i}"]=ch
    return df
def platform(gpl):
    stub=_stub(gpl)
    for kind,fn in [("annot",f"{gpl}.annot.gz"),("soft",f"{gpl}_family.soft.gz")]:
        try:
            url=f"{H}platforms/{stub}/{gpl}/{kind}/{fn}"
            p=fetch(url,fn); break
        except Exception as e: continue
    else: raise RuntimeError("no platform "+gpl)
    lines=[];on=False
    with gzip.open(p,"rt",errors="replace") as f:
        for l in f:
            if l.startswith("!platform_table_begin"): on=True;continue
            if l.startswith("!platform_table_end"): break
            if on: lines.append(l)
    return pd.read_csv(io.StringIO("".join(lines)),sep="\t",dtype=str,low_memory=False)
