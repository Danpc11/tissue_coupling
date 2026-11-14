import os as _os, sys as _sys
_RAW = _os.environ.get("T2D_RAW", "data/raw"); _OUT = _os.environ.get("T2D_OUT", "results")   # las mismas rutas que run_all.sh
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "../lib"))
import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
from scipy import stats
R=pd.read_csv(f"{_OUT}/new/meta_per_cohort.tsv",sep="\t")
q=R[R.contraste=="NR3C1~reloj"].copy()
q["post_mortem"]=q.cohorte.str.contains("GTEx").astype(float)
q["log_n"]=np.log10(q.n)
q["tipo"]=np.where(q.tejido.str.contains("adiposo"),"adiposo",np.where(q.tejido.str.contains("islotes"),"islotes","músculo"))   # tres tipos de tejido
print("METARREGRESION: cada moderador por separado, ponderada por 1/(se^2+tau2)")
def tau2_resid(X,z,se):
    """tau2 residual de la metarregresion por momentos: (Q_E-(k-p))/tr(P), con pesos 1/se^2."""
    w=1/se**2; W=np.diag(w); k,p=X.shape
    b=np.linalg.solve(X.T@W@X,X.T@W@z); QE=(w*(z-X@b)**2).sum()
    P=W-W@X@np.linalg.inv(X.T@W@X)@X.T@W
    return max(0,(QE-(k-p))/np.trace(P))
def r2_tau2(X):
    """R2 estandar de la metarregresion: 1 - tau2 residual / tau2 base (ambos por momentos)."""
    t=tau2_resid(X,q.z.values,q.se.values); return 0 if tau2<=0 else max(0,1-t/tau2)
def dl_tau2(z,se):
    w=1/se**2; m=(w*z).sum()/w.sum(); Q=(w*(z-m)**2).sum(); k=len(z)
    C=w.sum()-(w**2).sum()/w.sum(); return max(0,(Q-(k-1))/C), Q, k
tau2,Q,k=dl_tau2(q.z.values,q.se.values)
print(f"{len(q)} estimaciones, {q.cohorte.nunique()} cohortes, I2 global {max(0,(Q-(k-1))/Q)*100:.0f}%")
print(f"  tau2 base {tau2:.4f}  Q {Q:.1f} con {k-1} gl  (P {stats.chi2.sf(Q,k-1):.1e})\n")
for mod in ("post_mortem","log_n"):
    X=np.column_stack([np.ones(len(q)),q[mod].values]); w=1/(q.se.values**2+tau2)
    W=np.diag(w); b=np.linalg.solve(X.T@W@X, X.T@W@q.z.values)
    V=np.linalg.inv(X.T@W@X); se_b=np.sqrt(np.diag(V))
    zt=b[1]/se_b[1]; p=2*stats.t.sf(abs(zt),k-2)        # t con k-2 grados de libertad: con trece estimaciones la normal da P demasiado pequenos
    res=q.z.values-X@b; Qr=(w*res**2).sum()
    R2=max(0,1-Qr/Q)*100
    print(f"  {mod:12s} beta {b[1]:+.3f} (se {se_b[1]:.3f}) P {p:.3f}  Q residual {Qr:.1f}  1-Qres/Q {R2:.0f}%  R2 (reduccion de tau2) {r2_tau2(X)*100:.0f}%")
# tejido como factor (tres tipos: adiposo, musculo, islotes)
D=pd.get_dummies(q.tipo,drop_first=True).astype(float)
X=np.column_stack([np.ones(len(q)),D.values]); w=1/(q.se.values**2+tau2); W=np.diag(w)
b=np.linalg.solve(X.T@W@X,X.T@W@q.z.values); res=q.z.values-X@b; Qr=(w*res**2).sum()
print(f"  {'tejido':12s} {D.shape[1]+1} niveles       Q residual {Qr:.1f}  1-Qres/Q {max(0,1-Qr/Q)*100:.0f}%  R2 (reduccion de tau2) {r2_tau2(X)*100:.0f}%")
print(f"\n  por tejido:")
for t,g in q.groupby("tipo"):
    tt,_,kk=dl_tau2(g.z.values,g.se.values)
    print(f"    {t:12s} k={kk}  rho {np.tanh(np.average(g.z,weights=1/g.se**2)):+.2f}  tau2 {tt:.4f}")
