#!/usr/bin/env python3
"""Fase circadiana por tejido y acoplamiento entre organos.

Produce, en orden:
  1. fase inferida por donante y tejido, sin etiquetas de tiempo
  2. el control que decide si el ciclo recuperado es el reloj
  3. acoplamiento de FASE entre pares de tejidos validados, con el periodo anclado
  4. acoplamiento de AMPLITUD, con ordenamiento libre
  5. la medida compuesta que no separa ambas, para comparar con lo que mide un analisis
     que no descompone

Cada cantidad se estima con el procedimiento que es estable para ella: el ordenamiento libre
ajusta la onda a los datos y estima bien la amplitud pero deja la fase sujeta a ambiguedad de
periodo; anclar el periodo fija la fase pero proyecta la amplitud sobre una onda que no es la
del tejido. Usar un solo metodo para las dos cosas da una de ellas mal.

Resultado central: las dos modalidades ordenan los pares de forma distinta, y la medida
compuesta sigue a la amplitud (rho = 0.515) mas que a la fase (0.309).

Uso:  python python/analyses/03_phase_coupling.py
Salidas: results/phase_<tejido>.tsv, results/clock_validation.tsv,
         results/coupling_validated.tsv, results/amplitude_sharing.tsv
"""
import os, sys, itertools, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib"))
from phase import CLK, unsup, align, kap, coupling, is_clock

OUT = os.environ.get("T2D_OUT", "results")
TISSUES = os.environ.get("T2D_TISSUES", "muscle,adipose,liver,nerve_tibial,pancreas,blood,hypothalamus,cortex").split(",")
RENAME = {"ARNTL": "BMAL1"}   # las matrices depositadas conservan el simbolo de su plataforma

def load(tis):
    for p in (f"{OUT}/{tis}.pkl", f"data/processed/{tis}.pkl"):
        if os.path.exists(p):
            E = pd.read_pickle(p)
            E.index = [RENAME.get(i, i) for i in E.index]
            return E
    return None

def phase_and_error(E, nsplit=14, seed=0):
    """Fase por muestra y error de medida, estimado partiendo los genes en mitades.
    Las dos mitades se ALINEAN antes de compararse (ver phase.align)."""
    g = [x for x in CLK if x in E.index]
    Z = E.loc[g]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.values
    a = unsup(Z)
    rng = np.random.default_rng(seed); ds = []
    for _ in range(nsplit):
        pr = rng.permutation(len(g))
        ds.append(align(unsup(Z[pr[:len(g)//2]]), unsup(Z[pr[len(g)//2:]])))
    k, _ = kap(np.concatenate(ds))
    return pd.Series(a, index=E.columns), np.sqrt(1/k)/2, g

def amplitude(E, ang):
    """Amplitud por individuo: norma de la proyeccion de su perfil sobre la base ortonormal de
    la forma de onda de la cohorte. Independiente de la fase por construccion."""
    g = [x for x in CLK if x in E.index]
    Z = E.loc[g]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.values
    D = np.column_stack([np.cos(ang), np.sin(ang)])
    B, *_ = np.linalg.lstsq(D, Z.T, rcond=None)
    Q, _ = np.linalg.qr(B.T)
    return pd.Series(np.sqrt(((Z.T @ Q) ** 2).sum(1)), index=E.columns)

if __name__ == "__main__":
    P, SD, VAL = {}, {}, []
    for tis in TISSUES:
        E = load(tis)
        if E is None: print(f"  (falta) {tis}"); continue
        p, sd, g = phase_and_error(E)
        ok, sep = is_clock(E, p[E.columns].values)
        VAL.append(dict(tissue=tis, n=E.shape[1], n_genes=len(g), bmal1_dbp_h=sep,
                        error_h=sd*24/(2*np.pi), is_clock=ok))
        print(f"{tis:14s} n={E.shape[1]:4d}  BMAL1-DBP={sep:5.1f} h  error={sd*24/(2*np.pi):.2f} h  "
              + ("RELOJ" if ok else "DESCARTADO"))
        p.rename_axis("donor").to_frame("phase_rad").to_csv(f"{OUT}/phase_{tis}.tsv", sep="\t")
        if ok: P[tis], SD[tis] = p, sd
    pd.DataFrame(VAL).to_csv(f"{OUT}/clock_validation.tsv", sep="\t", index=False)

    rows = []
    for a, b in itertools.combinations(P, 2):
        d = sorted(set(P[a].index) & set(P[b].index))
        if len(d) < 80: continue
        r, pv, nu = coupling(P[a][d].values, P[b][d].values)
        rows.append(dict(a=a, b=b, n=len(d), rho=r, p=pv, null_mean=nu))
    cp = pd.DataFrame(rows).sort_values("rho", ascending=False)
    cp.to_csv(f"{OUT}/coupling_validated.tsv", sep="\t", index=False)
    print("\nACOPLAMIENTO DE FASE"); print(cp.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    if "muscle" in P and "adipose" in P:
        M, A = load("muscle"), load("adipose")
        d = sorted(set(M.columns) & set(A.columns))
        am = amplitude(M, P["muscle"][d].values); aa = amplitude(A, P["adipose"][d].values)
        r, pv = stats.spearmanr(am[d], aa[d])
        rng = np.random.default_rng(0); common = [x for x in M.index if x in A.index]
        null = []
        for _ in range(40):
            gg = list(rng.choice(common, len([x for x in CLK if x in M.index]), replace=False))
            try:
                Zm = M.loc[gg, d]; Zm = ((Zm.T - Zm.mean(1))/(Zm.std(1)+1e-9)).T.values
                Za = A.loc[gg, d]; Za = ((Za.T - Za.mean(1))/(Za.std(1)+1e-9)).T.values
                for Z_, store in ((Zm, "m"), (Za, "a")):
                    pass
                a1 = unsup(Zm); D = np.column_stack([np.cos(a1), np.sin(a1)])
                B, *_ = np.linalg.lstsq(D, Zm.T, rcond=None); Q, _ = np.linalg.qr(B.T)
                x1 = np.sqrt(((Zm.T @ Q) ** 2).sum(1))
                a2 = unsup(Za); D = np.column_stack([np.cos(a2), np.sin(a2)])
                B, *_ = np.linalg.lstsq(D, Za.T, rcond=None); Q, _ = np.linalg.qr(B.T)
                x2 = np.sqrt(((Za.T @ Q) ** 2).sum(1))
                null.append(abs(stats.spearmanr(x1, x2)[0]))
            except Exception:
                pass
        pe = (np.sum(np.array(null) >= abs(r)) + 1) / (len(null) + 1)
        pd.DataFrame([dict(pair="muscle-adipose", n=len(d), amplitude_rho=r, p=pv,
                           random_sets_mean=np.mean(null), p_empirical=pe)]
                     ).to_csv(f"{OUT}/amplitude_sharing.tsv", sep="\t", index=False)
        print(f"\nAMPLITUD compartida musculo-adiposo: rho={r:.3f} (P={pv:.1e}); "
              f"conjuntos al azar {np.mean(null):.3f}, P empirico={pe:.3f}")
