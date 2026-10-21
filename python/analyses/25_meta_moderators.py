import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy import stats
R=pd.read_csv("results/new/meta_per_cohort.tsv",sep="\t")
q=R[R.contraste=="NR3C1~reloj"].copy()
q["post_mortem"]=q.cohorte.str.contains("GTEx").astype(float)
q["log_n"]=np.log10(q.n)
print(f"{len(q)} estimaciones, {q.cohorte.nunique()} cohortes, I2 global 84%\n")
print("METARREGRESION: cada moderador por separado, ponderada por 1/(se^2+tau2)")
def dl_tau2(z,se):
    w=1/se**2; m=(w*z).sum()/w.sum(); Q=(w*(z-m)**2).sum(); k=len(z)
    C=w.sum()-(w**2).sum()/w.sum(); return max(0,(Q-(k-1))/C), Q, k
tau2,Q,k=dl_tau2(q.z.values,q.se.values)
print(f"  tau2 base {tau2:.4f}  Q {Q:.1f} con {k-1} gl  (P {1-stats.chi2.cdf(Q,k-1):.1e})\n")
for mod in ("post_mortem","log_n"):
    X=np.column_stack([np.ones(len(q)),q[mod].values]); w=1/(q.se.values**2+tau2)
    W=np.diag(w); b=np.linalg.solve(X.T@W@X, X.T@W@q.z.values)
    V=np.linalg.inv(X.T@W@X); se_b=np.sqrt(np.diag(V))
    zt=b[1]/se_b[1]; p=2*(1-stats.norm.cdf(abs(zt)))
    res=q.z.values-X@b; Qr=(w*res**2).sum()
    R2=max(0,1-Qr/Q)*100
    print(f"  {mod:12s} beta {b[1]:+.3f} (se {se_b[1]:.3f}) P {p:.3f}  Q residual {Qr:.1f}  R2 {R2:.0f}%")
# tejido como factor
D=pd.get_dummies(q.tejido,drop_first=True).astype(float)
X=np.column_stack([np.ones(len(q)),D.values]); w=1/(q.se.values**2+tau2); W=np.diag(w)
b=np.linalg.solve(X.T@W@X,X.T@W@q.z.values); res=q.z.values-X@b; Qr=(w*res**2).sum()
print(f"  {'tejido':12s} {D.shape[1]} niveles       Q residual {Qr:.1f}  R2 {max(0,1-Qr/Q)*100:.0f}%")
print(f"\n  por tejido:")
for t,g in q.groupby("tejido"):
    tt,_,kk=dl_tau2(g.z.values,g.se.values)
    print(f"    {t:12s} k={kk}  rho {np.tanh(np.average(g.z,weights=1/g.se**2)):+.2f}  tau2 {tt:.4f}")
