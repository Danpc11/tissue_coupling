import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0,"python/lib"); import phase as ph
X=pd.read_pickle("results/new/GSE129843_log2cpm.pkl"); md=pd.read_pickle("results/new/GSE129843_meta.pkl").set_index("gsm")
md["ind"]=md.title.str.extract(r"Participant: (\d+)")[0]
md["cond"]=np.where(md.c0.str.contains("Time restricted"),"TRF","EXF")
md["t"]=md.c1.str.replace("time point: ","").str.split(":").str[0].astype(float)
md=md.reindex(X.columns); Xc=X.sub(X.mean(0),axis=1); g=[x for x in ph.CLK if x in Xc.index]
def cosinor(y,t,per=24.0):
    w=2*np.pi/per; D=np.column_stack([np.ones(len(t)),np.cos(w*t),np.sin(w*t)])
    b,*_=np.linalg.lstsq(D,y,rcond=None); return np.hypot(b[1],b[2])
rows=[]
for (ind,cond),q in md.groupby(["ind","cond"]):
    cols=list(q.index); t=q.t.values
    if len(cols)<5: continue
    for held in cols:
        keep=[c for c in cols if c!=held]; tk=t[[cols.index(c) for c in keep]]
        amp_cos=np.mean([cosinor(np.asarray(Xc.loc[gg,keep],float),tk) for gg in g])
        others=[c for c in Xc.columns if c not in cols]
        Z=Xc.loc[g,others]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
        ang=ph.unsup(Z.values)
        a=np.array([np.mean(Z.loc[gg].values*np.cos(ang)) for gg in g])
        b=np.array([np.mean(Z.loc[gg].values*np.sin(ang)) for gg in g])
        mu=Xc.loc[g,others].mean(1).values; sd=Xc.loc[g,others].std(1).values+1e-9
        z=(np.asarray(Xc.loc[g,held],float)-mu)/sd
        th=np.arctan2(np.dot(z,b),np.dot(z,a))
        rows.append(dict(ind=ind,cond=cond,amp_inf=np.dot(z,a*np.cos(th)+b*np.sin(th))/np.sqrt(len(g)),amp_cos=amp_cos))
R=pd.DataFrame(rows); R.to_csv("results/new/GSE129843_loo_validation.tsv",sep="\t",index=False)
r=stats.spearmanr(R.amp_inf,R.amp_cos); print(f"por muestra retenida: rho {r.statistic:+.2f} (P {r.pvalue:.1e}, n={len(R)})")
M=R.groupby(["ind","cond"])[["amp_inf","amp_cos"]].mean().reset_index()
r=stats.spearmanr(M.amp_inf,M.amp_cos); print(f"promediado por serie:  rho {r.statistic:+.2f} (P {r.pvalue:.1e}, n={len(M)})")
