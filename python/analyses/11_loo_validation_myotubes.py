import os as _os, sys as _sys
_RAW = _os.environ.get("T2D_RAW", "data/raw"); _OUT = _os.environ.get("T2D_OUT", "results")   # las mismas rutas que run_all.sh
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../lib"))
"""Validacion de la amplitud inferida de muestra unica contra el cosinor del mismo donante, en miotubos sincronizados (GSE182117).

Para cada punto de cada serie (un donante en una condicion):
  - la referencia es la amplitud por cosinor de esa serie SIN el punto retenido;
  - la amplitud inferida del punto retenido usa una base (ordenamiento, plantilla y parametros de
    estandarizacion) construida sin esa serie.
El panel de genes es fijo y esta en el registro; no se elige con estos datos.

DOS EXCLUSIONES. `amp_inf` construye la base sin la serie del punto retenido, que es como se
calcularon los valores publicados: las muestras del mismo donante en la OTRA condicion siguen en la
base. `amp_inf_donor_out` la construye sin ninguna muestra de ese donante, que es lo que el texto
describe ("the whole donor was excluded"). Se escriben las dos y se informan las dos.
"""
import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
import phase as ph
X=pd.read_pickle(f"{_OUT}/new/GSE182117_log2cpm.pkl"); md=pd.read_pickle(f"{_OUT}/new/GSE182117_meta.pkl").set_index("gsm")
md["ind"]=md.c4.str.replace("individual: ",""); md["trt"]=md.c1.str.replace("treatment: ","")
md["t"]=md.c3.str.replace("time: ","").str.replace("h","").astype(float); md=md.reindex(X.columns)
Xc=X.sub(X.mean(0),axis=1); g=[x for x in ph.CLK if x in Xc.index]
def cosinor(y,t,per=24.0):
    w=2*np.pi/per; D=np.column_stack([np.ones(len(t)),np.cos(w*t),np.sin(w*t)])
    b,*_=np.linalg.lstsq(D,y,rcond=None); return np.hypot(b[1],b[2])
def base(others):
    """Ordenamiento, plantilla y parametros de estandarizacion, solo con las muestras `others`."""
    Z=Xc.loc[g,others]; Z=((Z.T-Z.mean(1))/(Z.std(1)+1e-9)).T
    ang=ph.unsup(Z.values)
    a=np.array([np.mean(Z.loc[gg].values*np.cos(ang)) for gg in g])
    b=np.array([np.mean(Z.loc[gg].values*np.sin(ang)) for gg in g])
    return a,b,Xc.loc[g,others].mean(1).values,Xc.loc[g,others].std(1).values+1e-9
def infer(held,B):
    a,b,mu,sd=B
    z=(np.asarray(Xc.loc[g,held],float)-mu)/sd
    th=np.arctan2(np.dot(z,b),np.dot(z,a))
    return np.dot(z,a*np.cos(th)+b*np.sin(th))/np.sqrt(len(g))
rows=[]
for (ind,trt),q in md.groupby(["ind","trt"]):
    cols=list(q.index); t=q.t.values
    if len(cols)<6: continue
    B_series=base([c for c in Xc.columns if c not in cols])            # sin esta serie
    B_donor=base([c for c in Xc.columns if md.ind[c]!=ind])            # sin ninguna muestra del donante
    for k,held in enumerate(cols):
        keep=[c for c in cols if c!=held]
        # cosinor del donante SIN el punto retenido
        amp_cos=np.mean([cosinor(np.asarray(Xc.loc[gg,keep],float),t[[cols.index(c) for c in keep]]) for gg in g])
        rows.append(dict(ind=ind,trt=trt,held=held,amp_inf=infer(held,B_series),amp_cos=amp_cos,
                         amp_inf_donor_out=infer(held,B_donor)))
R=pd.DataFrame(rows); R.to_csv(f"{_OUT}/new/GSE182117_loo_validation.tsv",sep="\t",index=False)
print(f"n puntos retenidos: {len(R)}")
for col,lab in (("amp_inf","base sin la serie del punto retenido"),("amp_inf_donor_out","base sin el donante entero")):
    print(f"{lab}:")
    for trt,q in R.groupby("trt"):
        r=stats.spearmanr(q[col],q.amp_cos); print(f"  {trt}: rho {r.statistic:+.2f} (P {r.pvalue:.1e}, n={len(q)})")
    r=stats.spearmanr(R[col],R.amp_cos); print(f"  juntas: rho {r.statistic:+.2f} (P {r.pvalue:.1e})")
    M=R.groupby(["ind","trt"])[[col,"amp_cos"]].mean().reset_index()
    for trt,q in M.groupby("trt"):
        r=stats.spearmanr(q[col],q.amp_cos); print(f"  promediado por donante, {trt}: rho {r.statistic:+.2f} (n={len(q)})")
