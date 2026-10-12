import numpy as np, pandas as pd
def read_gmt(path, minn=8, maxn=500):
    S={}
    for l in open(path):
        f=l.rstrip("\n").split("\t")
        g=[x for x in f[2:] if x]
        if minn<=len(g)<=maxn: S[f[0]]=set(g)
    return S
def enrich(stat, sets, nperm=2000, seed=0, minn=8):
    """stat: Series gen->valor. Compara la media del conjunto con conjuntos aleatorios
    del mismo tamano extraidos del MISMO universo medido."""
    rng=np.random.default_rng(seed); vals=stat.values; idx=stat.index
    pos={g:i for i,g in enumerate(idx)}
    rows=[]
    for name,gs in sets.items():
        ii=[pos[g] for g in gs if g in pos]
        if len(ii)<minn: continue
        obs=vals[ii].mean()
        null=np.array([vals[rng.integers(0,len(vals),len(ii))].mean() for _ in range(nperm)])
        p=(np.sum(np.abs(null-vals.mean())>=abs(obs-vals.mean()))+1)/(nperm+1)
        rows.append(dict(set=name,n=len(ii),mean=obs,z=(obs-null.mean())/(null.std()+1e-12),p=p))
    d=pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    m=len(d)
    # Benjamini-Hochberg: minimo acumulado DESDE EL FINAL sobre p ordenados de menor a mayor.
    # Un cummin hacia delante propaga el FDR del p mas pequeno y produce FDR no monotonos.
    q=(d.p.values*m/np.arange(1,m+1))
    d["fdr"]=np.minimum.accumulate(q[::-1])[::-1].clip(max=1)
    return d
