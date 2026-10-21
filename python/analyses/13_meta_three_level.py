#!/usr/bin/env python3
"""Metaanalisis de tres niveles: estimados anidados en cohorte (A4 de la revision).

DerSimonian-Laird asume estimados independientes. Aqui no lo son: GTEx aporta tres tejidos de
donantes solapados, y AAGMEx, GSE13070, GSE41168 y GSE159984 aportan dos cada uno. Se ajusta un
modelo de efectos aleatorios con dos componentes de varianza (tau2 entre cohortes, sigma2 entre
estimados dentro de cohorte) por maxima verosimilitud restringida, y se compara con el de dos
niveles. Tambien se reporta varianza robusta por clúster (RVE) como comprobacion.
"""
import os, numpy as np, pandas as pd
from scipy import stats, optimize
OUT = os.environ.get("T2D_OUT", "results")

def dl(z, se):
    w = 1/se**2; m = (w*z).sum()/w.sum(); Q = (w*(z-m)**2).sum(); k = len(z)
    tau2 = max(0, (Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()))
    w2 = 1/(se**2+tau2); m2 = (w2*z).sum()/w2.sum(); se2 = np.sqrt(1/w2.sum())
    I2 = max(0, (Q-(k-1))/Q*100) if Q > 0 else 0
    return m2, se2, tau2, I2, Q, k

def reml_3level(z, se, cl):
    """Dos componentes de varianza: tau2 entre cohortes, sg2 entre estimados dentro de cohorte."""
    cl = np.asarray(cl); idx = {c: i for i, c in enumerate(pd.unique(cl))}
    C = np.array([idx[c] for c in cl]); n = len(z)
    def nll(p):
        tau2, sg2 = np.exp(p)
        V = np.diag(se**2 + sg2) + tau2*(C[:, None] == C[None, :])
        try: Vi = np.linalg.inv(V)
        except np.linalg.LinAlgError: return 1e9
        X = np.ones((n, 1)); XtViX = X.T@Vi@X
        b = np.linalg.solve(XtViX, X.T@Vi@z); r = z - (X@b).ravel()
        s, ld = np.linalg.slogdet(V)
        return .5*(ld + np.log(XtViX[0, 0]) + r@Vi@r)
    best = min((optimize.minimize(nll, s0, method="Nelder-Mead",
                                  options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
                for s0 in ([-4, -4], [-2, -2], [-6, -2], [-2, -6])), key=lambda r: r.fun)
    tau2, sg2 = np.exp(best.x)
    V = np.diag(se**2 + sg2) + tau2*(C[:, None] == C[None, :]); Vi = np.linalg.inv(V)
    X = np.ones((n, 1)); XtViX = X.T@Vi@X
    m = float(np.linalg.solve(XtViX, X.T@Vi@z).ravel()[0]); sem = float(np.sqrt(1/XtViX[0, 0]))
    return m, sem, tau2, sg2

def rve(z, se, cl, tau2):
    """Varianza robusta por clúster con correccion de muestra pequeña (Hedges)."""
    cl = np.asarray(cl); w = 1/(se**2+tau2); m = (w*z).sum()/w.sum()
    num = sum((w[cl == c]*(z[cl == c]-m)).sum()**2 for c in pd.unique(cl))
    M = len(pd.unique(cl)); V = num/w.sum()**2 * M/(M-1)
    return m, np.sqrt(V), M

if __name__ == "__main__":
    R = pd.read_csv(f"{OUT}/new/meta_per_cohort.tsv", sep="\t")
    for contrast in ("NR3C1~reloj", "GC~reloj"):
        q = R[R.contraste == contrast].copy()
        z, se, cl = q.z.values, q.se.values, q.cohorte.values
        m2, se2, tau2, I2, Q, k = dl(z, se)
        m3, se3, t3, s3 = reml_3level(z, se, cl)
        mr, ser, M = rve(z, se, cl, t3)
        print(f"\n=== {contrast}  ({k} estimados, {M} cohortes)")
        print(f"  dos niveles (DerSimonian-Laird): rho {np.tanh(m2):+.3f} "
              f"(IC95 {np.tanh(m2-1.96*se2):+.3f} a {np.tanh(m2+1.96*se2):+.3f})  "
              f"tau2 {tau2:.4f}  I2 {I2:.0f}%")
        print(f"  tres niveles (REML, anidado):    rho {np.tanh(m3):+.3f} "
              f"(IC95 {np.tanh(m3-1.96*se3):+.3f} a {np.tanh(m3+1.96*se3):+.3f})  "
              f"tau2 entre cohortes {t3:.4f}  sigma2 dentro {s3:.4f}")
        print(f"  varianza robusta por clúster:    rho {np.tanh(mr):+.3f} "
              f"(IC95 {np.tanh(mr-1.96*ser):+.3f} a {np.tanh(mr+1.96*ser):+.3f})")
        print(f"  anchura del IC: 2 niveles {np.tanh(m2+1.96*se2)-np.tanh(m2-1.96*se2):.3f} | "
              f"3 niveles {np.tanh(m3+1.96*se3)-np.tanh(m3-1.96*se3):.3f} | "
              f"RVE {np.tanh(mr+1.96*ser)-np.tanh(mr-1.96*ser):.3f}")
