#!/usr/bin/env python3
"""Fase y amplitud circadianas por tejido, y acoplamiento entre organos.

Produce, en orden:
  1. filtro: que tejidos tienen un eje ciclico que ES el reloj (antifase BMAL1-DBP)
  2. FASE con el periodo anclado y AMPLITUD con ordenamiento libre, solo en los validados
  3. acoplamiento por pares en las dos modalidades, mas la medida compuesta que no las separa

POR QUE DOS PROCEDIMIENTOS. El ordenamiento libre ajusta la onda a los datos, de modo que la base
de proyeccion es la del tejido y la AMPLITUD sale bien, pero no fija el periodo: si el ciclo
dominante es el armonico de 12 h en un tejido y el fundamental de 24 h en el otro, las dos fases
quedan en escalas distintas y su relacion aparece como un factor de dos. Eso produce un patron
axial (modulo pi) que la correlacion circular ordinaria, sensible solo a pendiente +1, lee como
ausencia de acoplamiento. Anclar la onda a una plantilla de fases relativas obtenida en una cohorte
con horas conocidas fija el periodo y resuelve la FASE, pero entonces la base de proyeccion no es
la del tejido y la amplitud se degrada.

Rangos de estabilidad medidos (musculo-adiposo, plantillas de 13 a 200 genes):
    fase      anclado 0.32-0.47  |  libre 0.21-0.62
    amplitud  libre 0.30-0.37    |  anclado 0.07-0.38
Usar un solo metodo para las dos cosas da una de ellas mal.

Uso:  python python/analyses/03_phase_coupling.py
Salidas: results/clock_validation.tsv, results/phase_<tejido>.tsv, results/coupling_map.tsv
"""
import os
import sys
import itertools
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib"))
from phase import (CLK, unsup, is_clock, coupling, order_anchored, select_template,
                   amplitude, double_angle)

OUT = os.environ.get("T2D_OUT", "results")
# Tamano de la plantilla que ancla la fase. 60 es el valor con el que se calcularon las cifras
# publicadas; la estimacion es estable entre 13 y 200 genes (ESM Fig. 1D) pero NO con la plantilla
# completa de miles de genes.
N_TEMPLATE = int(os.environ.get("T2D_N_TEMPLATE", "60"))
TISSUES = os.environ.get(
    "T2D_TISSUES",
    "muscle,adipose,liver,nerve_tibial,pancreas,blood,hypothalamus,cortex").split(",")
RENAME = {"ARNTL": "BMAL1"}   # las matrices depositadas conservan el simbolo de su plataforma


def load(tis):
    """Carga la matriz de un tejido y normaliza el simbolo de BMAL1."""
    for p in (f"{OUT}/{tis}.pkl", f"data/processed/{tis}.pkl"):
        if os.path.exists(p):
            E = pd.read_pickle(p)
            E.index = [RENAME.get(i, i) for i in E.index]
            return E
    return None


def composite(Ea, Eb, donors, genes):
    """Medida que NO separa fase de amplitud: correlacion media, gen a gen, del perfil entre los
    dos tejidos. Es lo que calcula un analisis que no descompone, y por eso sirve de comparacion:
    resulta seguir a la amplitud mas que a la fase."""
    g = [x for x in genes if x in Ea.index and x in Eb.index]
    Za = Ea.loc[g, donors]
    Za = ((Za.T - Za.mean(1)) / (Za.std(1) + 1e-9)).T.values
    Zb = Eb.loc[g, donors]
    Zb = ((Zb.T - Zb.mean(1)) / (Zb.std(1) + 1e-9)).T.values
    return float(np.mean([abs(stats.spearmanr(Za[i], Zb[i])[0]) for i in range(len(g))]))


if __name__ == "__main__":
    template = None
    tpath = os.path.join(OUT, "muscle_rhythmic_genes.tsv")
    if os.path.exists(tpath):
        template = pd.read_csv(tpath, sep="\t").set_index("gene")
    else:
        print("Falta results/muscle_rhythmic_genes.tsv, la plantilla de fases relativas obtenida")
        print("por cosinor en una cohorte con horas conocidas. Sin ella la fase no puede anclarse")
        print("y se reporta con ordenamiento libre, sujeta a ambiguedad de periodo.")

    if template is not None:
        sel = select_template(template, N_TEMPLATE)
        print(f"Plantilla de fase: {len(sel)} genes (de {len(template)}), los mas ritmicos")
        print(f"  {', '.join(sel[:10])}{' ...' if len(sel) > 10 else ''}")
    print("\nFiltro: la antifase BMAL1-DBP debe quedar entre 8 y 16 h")
    E_, FREE, VAL = {}, {}, []
    for tis in TISSUES:
        E = load(tis)
        if E is None:
            print(f"  (falta) {tis}")
            continue
        g = [x for x in CLK if x in E.index]
        Z = E.loc[g]
        Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.values
        a = unsup(Z)
        ok, sep = is_clock(E, a)
        E_[tis] = E
        FREE[tis] = pd.Series(a, index=E.columns)
        VAL.append(dict(tissue=tis, n=E.shape[1], n_genes=len(g),
                        bmal1_dbp_h=sep, is_clock=bool(ok)))
        print(f"  {tis:14s} n={E.shape[1]:4d}  BMAL1-DBP={sep:5.1f} h  "
              + ("RELOJ" if ok else "descartado"))
    val = pd.DataFrame(VAL)
    val.to_csv(f"{OUT}/clock_validation.tsv", sep="\t", index=False)
    good = list(val.tissue[val.is_clock])
    if len(good) < 2:
        print("\nMenos de dos tejidos validados; no hay pares que comparar.")
        sys.exit(0)

    PH, AM = {}, {}
    for t in good:
        AM[t] = amplitude(E_[t], FREE[t][E_[t].columns].values)
        if template is not None:
            PH[t] = order_anchored(E_[t], template, n_template=N_TEMPLATE)[0]
        else:
            PH[t] = FREE[t]
        PH[t].rename_axis("donor").to_frame("phase_rad").to_csv(f"{OUT}/phase_{t}.tsv", sep="\t")

    rng = np.random.default_rng(0)
    rows = []
    for a, b in itertools.combinations(good, 2):
        d = sorted(set(PH[a].index) & set(PH[b].index))
        if len(d) < 80:
            continue
        rp, pp, _ = coupling(PH[a][d].values, PH[b][d].values)
        ra = abs(stats.spearmanr(AM[a][d].values, AM[b][d].values)[0])
        null = [abs(stats.spearmanr(AM[a][d].values,
                                    AM[b][d].values[rng.permutation(len(d))])[0])
                for _ in range(800)]
        rows.append(dict(a=a, b=b, n=len(d),
                         phase=rp, p_phase=pp,
                         amplitude=ra, p_amp=(np.sum(np.array(null) >= ra) + 1) / 801,
                         composite=composite(E_[a], E_[b], d, CLK),
                         double_angle=double_angle(PH[a][d].values, PH[b][d].values)))
    cp = pd.DataFrame(rows).sort_values("composite", ascending=False)
    cp.to_csv(f"{OUT}/coupling_map.tsv", sep="\t", index=False)
    print("\nACOPLAMIENTO: fase anclada, amplitud libre")
    print(cp.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    if len(cp) > 2:
        ca = stats.spearmanr(cp.composite, cp.amplitude)[0]
        cf = stats.spearmanr(cp.composite, cp.phase)[0]
        print(f"\nLa medida compuesta sigue a la amplitud (rho = {ca:.2f}) mas que a la fase "
              f"({cf:.2f}).")
        print("Un analisis que no descompone las dos modalidades mide sobre todo la primera.")
        hi = cp.double_angle.max()
        if hi > 0.12:
            print(f"AVISO: double_angle llega a {hi:.3f}; por encima de ~0.12 queda ambiguedad de")
            print("periodo y la fase de ese par no debe interpretarse.")

    # --- 4. cuanta varianza de una puntuacion depende de la fase ------------------------------
    # Separa el componente dependiente de fase (que NO se comparte entre tejidos) del nivel
    # (que si se comparte). Se calcula con la fase ANCLADA, no con la libre.
    spath = os.path.join(OUT, "axis_score_donors.tsv")
    if os.path.exists(spath) and template is not None:
        sc = pd.read_csv(spath, sep="\t", index_col=0).iloc[:, 0]
        rows = []
        for t in good:
            don = [d for d in PH[t].index if d in sc.index]
            if len(don) < 80:
                continue
            a = PH[t][don].values
            y = sc[don].values
            D = np.column_stack([np.ones(len(a)), np.cos(a), np.sin(a)])
            b, *_ = np.linalg.lstsq(D, y, rcond=None)
            r2 = 1 - ((y - D @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            null = []
            for _ in range(2000):
                ap = a[rng.permutation(len(a))]
                Dp = np.column_stack([np.ones(len(ap)), np.cos(ap), np.sin(ap)])
                bp, *_ = np.linalg.lstsq(Dp, y, rcond=None)
                null.append(1 - ((y - Dp @ bp) ** 2).sum() / ((y - y.mean()) ** 2).sum())
            rows.append(dict(tissue=t, n=len(don), r2_phase_dependent=r2,
                             amplitude=np.sqrt(b[1] ** 2 + b[2] ** 2),
                             peak_phase_h=(np.arctan2(b[2], b[1]) * 24 / (2 * np.pi)) % 24,
                             p_perm=(np.sum(np.array(null) >= r2) + 1) / 2001,
                             null_mean_r2=float(np.mean(null))))
        if rows:
            v = pd.DataFrame(rows)
            v.to_csv(f"{OUT}/axis_variance_phase.tsv", sep="\t", index=False)
            print("\nVARIANZA DEL EJE DEPENDIENTE DE LA FASE (fase anclada)")
            print(v.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
