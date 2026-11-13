import os as _os, sys as _sys
_RAW = _os.environ.get("T2D_RAW", "data/raw"); _OUT = _os.environ.get("T2D_OUT", "results")   # las mismas rutas que run_all.sh
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../lib"))
"""Validacion de la amplitud inferida de muestra unica contra el cosinor del mismo donante, en musculo esqueletico in vivo (GSE129843).

Para cada punto de cada serie (un donante en una condicion):
  - la referencia es la amplitud por cosinor de esa serie SIN el punto retenido;
  - la amplitud inferida del punto retenido usa una base (ordenamiento, plantilla y parametros de
    estandarizacion) construida sin esa serie.
El panel de genes es fijo y esta en el registro; no se elige con estos datos.

DOS GEOMETRIAS. `amp_inf` y `amp_inf_donor_out` usan las proyecciones directas sobre los vectores
de la forma de onda, que es como se calcularon los valores publicados; esos vectores no tienen la
misma norma, y la ganancia de esa amplitud cambia con la fase. `amp_ls` y `amp_ls_donor_out` usan
la proyeccion por minimos cuadrados, que no tiene ese defecto (python/lib/phase.py).

DOS EXCLUSIONES. `amp_inf` construye la base sin la serie del punto retenido, que es como se
calcularon los valores publicados: las muestras del mismo donante en la OTRA condicion siguen en la
base. `amp_inf_donor_out` la construye sin ninguna muestra de ese donante, que es lo que el texto
describe ("the whole donor was excluded"). Se escriben las dos y se informan las dos.
"""
import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings("ignore")
from scipy import stats
import phase as ph
X=pd.read_pickle(f"{_OUT}/new/GSE129843_log2cpm.pkl"); md=pd.read_pickle(f"{_OUT}/new/GSE129843_meta.pkl").set_index("gsm")
md["ind"]=md.title.str.extract(r"Participant: (\d+)")[0]
md["cond"]=np.where(md.c0.str.contains("Time restricted"),"TRF","EXF")
md["t"]=md.c1.str.replace("time point: ","").str.split(":").str[0].astype(float)
md=md.reindex(X.columns); Xc=X.sub(X.mean(0),axis=1); g=[x for x in ph.CLK if x in Xc.index]
def cosinor(y,t,per=24.0):
    w=2*np.pi/per; D=np.column_stack([np.ones(len(t)),np.cos(w*t),np.sin(w*t)])
    b,*_=np.linalg.lstsq(D,y,rcond=None); return np.hypot(b[1],b[2])
def base(others):
    """Referencia fija (python/lib/phase.py) construida solo con las muestras `others`."""
    return ph.fit_reference(Xc[others], g)
def infer(held,B,geometry="as published"):
    return ph.infer_from_reference(Xc[held], B, geometry=geometry)[1]
rows=[]
for (ind,cond),q in md.groupby(["ind","cond"]):
    cols=list(q.index); t=q.t.values
    if len(cols)<5: continue
    B_series=base([c for c in Xc.columns if c not in cols])            # sin esta serie
    B_donor=base([c for c in Xc.columns if md.ind[c]!=ind])            # sin ninguna muestra del participante
    for held in cols:
        keep=[c for c in cols if c!=held]; tk=t[[cols.index(c) for c in keep]]
        amp_cos=np.mean([cosinor(np.asarray(Xc.loc[gg,keep],float),tk) for gg in g])
        rows.append(dict(ind=ind,cond=cond,amp_inf=infer(held,B_series),amp_cos=amp_cos,
                         amp_inf_donor_out=infer(held,B_donor),
                         amp_ls=infer(held,B_series,"least squares"),amp_ls_donor_out=infer(held,B_donor,"least squares")))
R=pd.DataFrame(rows); R.to_csv(f"{_OUT}/new/GSE129843_loo_validation.tsv",sep="\t",index=False)
for col,lab in (("amp_inf","base sin la serie del punto retenido"),("amp_ls","lo mismo, geometria por minimos cuadrados"),
                ("amp_ls_donor_out","minimos cuadrados y base sin el donante entero"),("amp_inf_donor_out","base sin el participante entero")):
    r=stats.spearmanr(R[col],R.amp_cos); print(f"{lab}:\n  por muestra retenida: rho {r.statistic:+.2f} (P {r.pvalue:.1e}, n={len(R)})")
    M=R.groupby(["ind","cond"])[[col,"amp_cos"]].mean().reset_index()
    r=stats.spearmanr(M[col],M.amp_cos); print(f"  promediado por serie:  rho {r.statistic:+.2f} (P {r.pvalue:.1e}, n={len(M)})")
