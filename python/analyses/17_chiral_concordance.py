#!/usr/bin/env python3
"""Concordancia con CHIRAL: el acoplamiento de fase entre tejidos no depende del estimador (A5).

CHIRAL (Talamanca, Gobet y Naef, Science 2023) infiere fase circadiana de muestra unica. Se calcula
con el, y con el estimador de este trabajo, la misma cantidad: R, la concentracion de la diferencia
de fase entre dos tejidos del mismo donante.

    d_i = arg( exp(i (phi_musculo,i - phi_adiposo,i)) )      diferencia envuelta al circulo
    R   = | mean_i exp(i d_i) |                              concentracion, 0 sin acoplar, 1 rigido

La significancia se evalua permutando la correspondencia de donantes entre tejidos (2.000 veces).
El argumento no es que un estimador sea mejor, sino que el resultado es el mismo con los dos.

AVISO del modulo: la documentacion de scRITMO dice que las muestras van en columnas, pero el codigo
las espera en FILAS. Se transpone antes de llamar.
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
import phase as ph
OUT = os.environ.get("T2D_OUT", "results")

CLK = ["ARNTL", "NPAS2", "CLOCK", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "DBP", "TEF", "HLF",
       "CIART", "CRY1", "CRY2", "BHLHE40", "NFIL3", "RORA", "PDK4", "CEBPD", "FKBP5", "TXNIP", "KLF15"]

def concentration(a, b, seed=0, B=2000):
    """R de la diferencia de fase y su P por permutacion de donantes."""
    d = np.angle(np.exp(1j*(np.asarray(a) - np.asarray(b))))
    R = np.abs(np.mean(np.exp(1j*d)))
    rng = np.random.default_rng(seed)
    null = np.array([np.abs(np.mean(np.exp(1j*(a - np.asarray(b)[rng.permutation(len(b))]))))
                     for _ in range(B)])
    return R, (np.sum(null >= R) + 1)/(B + 1)

def phase_ours(E, genes=None):
    """Fase con el periodo ANCLADO a la plantilla, no con ordenamiento libre.

    El ordenamiento libre da bien la amplitud pero no fija el periodo, y la diferencia de fase
    entre tejidos queda en escalas distintas (ver cabecera de 03_phase_coupling.py). Para fase va
    el anclado, con los 60 genes mas ritmicos de la plantilla.
    """
    T = pd.read_csv(f"{OUT}/template_2harm.tsv", sep="\t", index_col=0)
    # order_anchored espera columnas r2, phase_h y amp; la plantilla guarda armonicos
    T = T.assign(phase_h=(np.arctan2(T.b1, T.a1) % (2*np.pi))*24/(2*np.pi), amp=T.amp1)
    return ph.order_anchored(E, T)[0].values   # devuelve (fases, genes usados)

def phase_chiral(E, genes):
    """CHIRAL vía scRITMO. Se transpone: el codigo espera muestras en filas (ver aviso)."""
    from chiral_backend import infer_phase
    return infer_phase(E.loc[genes], genes)["phi"]

def run(Em, Ea, donors, label, use_chiral=True):
    g = [x for x in CLK if x in Em.index and x in Ea.index]
    rows = []
    methods = [("this work", lambda E, gg: phase_ours(E))] + ([("CHIRAL", phase_chiral)] if use_chiral else [])
    for name, f in methods:
        try:
            R, p = concentration(f(Em[donors], g), f(Ea[donors], g))
            rows.append(dict(cohort=label, method=name, n_donors=len(donors), n_genes=len(g), R=R, P=p))
            print(f"  {label:10s} {name:10s} n={len(donors):4d}  R = {R:.3f}  P = {p:.4f}")
        except ImportError:
            print(f"  {label:10s} {name:10s} no disponible: falta python/lib/chiral_backend.py")
    return rows

if __name__ == "__main__":
    out = []
    # AAGMEx: musculo y adiposo del mismo donante
    Em = pd.read_pickle(f"{OUT}/new/AAGMEx_muscle_log2.pkl"); Ea = pd.read_pickle(f"{OUT}/new/AAGMEx_adipose_log2.pkl")
    don = sorted(set(Em.columns) & set(Ea.columns))
    out += run(Em.sub(Em.mean(0), axis=1), Ea.sub(Ea.mean(0), axis=1), don, "AAGMEx")
    # GTEx: musculo y adiposo subcutaneo del mismo donante
    Em = pd.read_pickle(f"{OUT}/new/GTEx_muscle.pkl"); Ea = pd.read_pickle(f"{OUT}/new/GTEx_adipose_sc.pkl")
    don = sorted(set(Em.columns) & set(Ea.columns))
    out += run(Em.sub(Em.mean(0), axis=1), Ea.sub(Ea.mean(0), axis=1), don, "GTEx")
    pd.DataFrame(out).to_csv(f"{OUT}/new/chiral_concordance.tsv", sep="\t", index=False)
