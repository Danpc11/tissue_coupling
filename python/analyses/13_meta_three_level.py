#!/usr/bin/env python3
"""Metaanalisis de tres niveles: estimados anidados en cohorte, y cuanto depende de las decisiones.

DerSimonian-Laird asume estimados independientes. Aqui no lo son: GTEx aporta tres tejidos de
donantes solapados, y AAGMEx, GSE13070, GSE41168 y GSE159984 aportan dos cada uno. Se ajusta un
modelo de efectos aleatorios con dos componentes de varianza (tau2 entre cohortes, sigma2 entre
estimados dentro de cohorte) por maxima verosimilitud restringida, y se compara con el de dos
niveles y con la varianza robusta por cluster.

INTERVALOS. Con siete cohortes, el cuantil normal 1.96 da intervalos demasiado estrechos. Se
informan los dos: el normal, que es el que se cito, y el de la t con (cohortes - 1) grados de
libertad, que es el que corresponde.

SENSIBILIDAD, toda en la misma tabla:
  - dejando fuera cada cohorte;
  - una sola estimacion por cohorte (la de mayor n), que elimina la dependencia de raiz;
  - con el error estandar por bootstrap de individuos en lugar de la aproximacion de Fisher;
  - con la expresion sin centrar por muestra;
  - con el grupo diagnostico o experimental como covariable, donde la cohorte lo tiene;
  - el intervalo de PREDICCION: donde cabe esperar la estimacion de una cohorte nueva. Con mucha
    heterogeneidad es mucho mas ancho que el intervalo de confianza de la media, y es el que dice
    si el efecto es positivo en todas partes o solo en promedio.

Salida: results/new/meta_three_level.tsv
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

def summary(z, se, cl):
    """Los tres estimadores con sus intervalos, normal y t, y el intervalo de prediccion."""
    M = len(pd.unique(cl)); tq = stats.t.ppf(.975, M - 1) if M > 1 else np.nan
    m2, se2, tau2, I2, Q, k = dl(z, se)
    m3, se3, t3, s3 = reml_3level(z, se, cl)
    mr, ser, _ = rve(z, se, cl, t3)
    ci = lambda m, s_, q: (np.tanh(m - q*s_), np.tanh(m + q*s_))
    pq = stats.t.ppf(.975, M - 2) if M > 2 else np.nan
    pi = pq*np.sqrt(se3**2 + t3 + s3)
    return dict(k=k, cohorts=M, I2=I2,
                rho_2level=np.tanh(m2), lo_2level=ci(m2, se2, 1.96)[0], hi_2level=ci(m2, se2, 1.96)[1],
                rho_3level=np.tanh(m3), lo_3level=ci(m3, se3, 1.96)[0], hi_3level=ci(m3, se3, 1.96)[1],
                lo_3level_t=ci(m3, se3, tq)[0], hi_3level_t=ci(m3, se3, tq)[1],
                tau2_between=t3, sigma2_within=s3,
                rho_rve=np.tanh(mr), lo_rve=ci(mr, ser, 1.96)[0], hi_rve=ci(mr, ser, 1.96)[1],
                lo_rve_t=ci(mr, ser, tq)[0], hi_rve_t=ci(mr, ser, tq)[1],
                p_rve_t=2*stats.t.sf(abs(mr/ser), M - 1) if M > 1 else np.nan,
                pred_lo=np.tanh(m3 - pi), pred_hi=np.tanh(m3 + pi))


if __name__ == "__main__":
    R = pd.read_csv(f"{OUT}/new/meta_per_cohort.tsv", sep="\t")
    rows = []
    for contrast in ("NR3C1~reloj", "GC~reloj"):
        q = R[R.contraste == contrast].copy()
        z, se, cl = q.z.values, q.se.values, q.cohorte.values
        S = summary(z, se, cl); rows.append(dict(contrast=contrast, analysis="main", **S))
        print(f"\n=== {contrast}  ({S['k']} estimados, {S['cohorts']} cohortes)"
              + ("" if q.get("complete", pd.Series([True])).all() else "   INCOMPLETO: faltan estimaciones del articulo"))
        print(f"  dos niveles (DerSimonian-Laird): rho {S['rho_2level']:+.3f} (IC95 {S['lo_2level']:+.3f} a {S['hi_2level']:+.3f})  I2 {S['I2']:.0f}%")
        print(f"  tres niveles (REML, anidado):    rho {S['rho_3level']:+.3f} (IC95 {S['lo_3level']:+.3f} a {S['hi_3level']:+.3f}; "
              f"con t de {S['cohorts']-1} gl {S['lo_3level_t']:+.3f} a {S['hi_3level_t']:+.3f})  "
              f"tau2 entre cohortes {S['tau2_between']:.4f}  sigma2 dentro {S['sigma2_within']:.4f}")
        print(f"  varianza robusta por clúster:    rho {S['rho_rve']:+.3f} (IC95 {S['lo_rve']:+.3f} a {S['hi_rve']:+.3f}; "
              f"con t {S['lo_rve_t']:+.3f} a {S['hi_rve_t']:+.3f}, P = {S['p_rve_t']:.4f})")
        print(f"  intervalo de prediccion al 95 % para una cohorte nueva: {S['pred_lo']:+.2f} a {S['pred_hi']:+.2f}")
        # --- sensibilidad ---
        loco = []
        for c in pd.unique(cl):
            keep = cl != c
            if len(pd.unique(cl[keep])) < 3: continue
            T = summary(z[keep], se[keep], cl[keep]); rows.append(dict(contrast=contrast, analysis=f"without {c}", **T))
            loco.append(f"sin {c} {T['rho_3level']:+.2f}")
        print("  dejando una cohorte fuera (tres niveles): " + "; ".join(loco))
        one = q.sort_values("n", ascending=False).drop_duplicates("cohorte")
        m1, s1, _, I1, _, k1 = dl(one.z.values, one.se.values); t1 = stats.t.ppf(.975, k1 - 1)
        rows.append(dict(contrast=contrast, analysis="one estimate per cohort (largest n)", k=k1, cohorts=k1, I2=I1,
                         rho_2level=np.tanh(m1), lo_2level=np.tanh(m1 - t1*s1), hi_2level=np.tanh(m1 + t1*s1)))
        print(f"  una estimacion por cohorte (la de mayor n, k = {k1}): rho {np.tanh(m1):+.3f} (IC95 con t {np.tanh(m1 - t1*s1):+.3f} a {np.tanh(m1 + t1*s1):+.3f})")
        if "se_boot" in q and q.se_boot.notna().all():
            T = summary(z, q.se_boot.values, cl); rows.append(dict(contrast=contrast, analysis="bootstrap standard errors", **T))
            print(f"  con error estandar por bootstrap: rho {T['rho_3level']:+.3f} (IC95 con t {T['lo_3level_t']:+.3f} a {T['hi_3level_t']:+.3f}); "
                  f"razon mediana bootstrap / Fisher {np.median(q.se_boot/q.se):.2f}")
        if "rho_uncentred" in q and q.rho_uncentred.notna().all():
            zu = np.arctanh(np.clip(q.rho_uncentred.values, -.999, .999))
            T = summary(zu, se, cl); rows.append(dict(contrast=contrast, analysis="expression not centred per sample", **T))
            print(f"  sin centrar por muestra: rho {T['rho_3level']:+.3f} (IC95 con t {T['lo_3level_t']:+.3f} a {T['hi_3level_t']:+.3f})")
        if "rho_adjusted_for_group" in q and q.rho_adjusted_for_group.notna().any():
            ra = q.rho_adjusted_for_group.where(q.rho_adjusted_for_group.notna(), q.rho)
            T = summary(np.arctanh(np.clip(ra.values, -.999, .999)), se, cl)
            rows.append(dict(contrast=contrast, analysis="adjusted for diagnostic or study group where the cohort has one", **T))
            print(f"  con el grupo diagnostico como covariable en las {int(q.rho_adjusted_for_group.notna().sum())} estimaciones que lo tienen: "
                  f"rho {T['rho_3level']:+.3f} (IC95 con t {T['lo_3level_t']:+.3f} a {T['hi_3level_t']:+.3f})")
    pd.DataFrame(rows).to_csv(f"{OUT}/new/meta_three_level.tsv", sep="\t", index=False)
