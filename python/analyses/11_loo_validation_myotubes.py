import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/lib"); import phase as ph
X=pd.read_pickle("results/new/GSE182117_log2cpm.pkl"); md=pd.read_pickle("results/new/GSE182117_meta.pkl").set_index("gsm")
md["ind"]=md.c4.str.replace("individual: ",""); md["trt"]=md.c1.str.replace("treatment: ","")
md["t"]=md.c3.str.replace("time: ","").str.replace("h","").astype(float); md=md.reindex(X.columns)
Xc=X.sub(X.mean(0),axis=1); g=[x for x in ph.CLK if x in Xc.index]
def cosinor(y,t,per=24.0):
    w=2*np.pi/per; D=np.column_stack([np.ones(len(t)),np.cos(w*t),np.sin(w*t)])
    b,*_=np.linalg.lstsq(D,y,rcond=None); return np.hypot(b[1],b[2])
rows=[]
for (ind,trt),q in md.groupby(["ind","trt"]):
    cols=list(q.index); t=q.t.values
    if len(cols)<6: continue
    for k,held in enumerate(cols):
        keep=[c for c in cols if c!=held]
        # cosinor del donante SIN el punto retenido
        amp_cos=np.mean([cosinor(np.asarray(Xc.loc[gg,keep],float),t[[cols.index(c) for c in keep]]) for gg in g])
        # amplitud inferida del punto retenido, con la base construida SIN ese donante
        others=[c for c in Xc.columns if c not in cols]
        Z=Xc.loc[g,others]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
        ang=ph.unsup(Z.values)
        tpl=pd.DataFrame({"a":[np.mean(Z.loc[gg].values*np.cos(ang)) for gg in g],
                          "b":[np.mean(Z.loc[gg].values*np.sin(ang)) for gg in g]},index=g)
        mu=Xc.loc[g,others].mean(1); sd=Xc.loc[g,others].std(1)+1e-9
        z=(np.asarray(Xc.loc[g,held],float)-mu.values)/sd.values
        th=np.arctan2(np.dot(z,tpl.b.values),np.dot(z,tpl.a.values))
        amp_inf=np.dot(z,tpl.a.values*np.cos(th)+tpl.b.values*np.sin(th))/np.sqrt(len(g))
        rows.append(dict(ind=ind,trt=trt,held=held,amp_inf=amp_inf,amp_cos=amp_cos))
R=pd.DataFrame(rows); R.to_csv("results/new/GSE182117_loo_validation.tsv",sep="\t",index=False)
print(f"n puntos retenidos: {len(R)}")
for trt,q in R.groupby("trt"):
    r=stats.spearmanr(q.amp_inf,q.amp_cos); print(f"  {trt}: rho {r.statistic:+.2f} (P {r.pvalue:.1e}, n={len(q)})")
r=stats.spearmanr(R.amp_inf,R.amp_cos); print(f"  juntas: rho {r.statistic:+.2f} (P {r.pvalue:.1e})")
M=R.groupby(["ind","trt"])[["amp_inf","amp_cos"]].mean().reset_index()
for trt,q in M.groupby("trt"):
    r=stats.spearmanr(q.amp_inf,q.amp_cos); print(f"  promediado por donante, {trt}: rho {r.statistic:+.2f} (n={len(q)})")
