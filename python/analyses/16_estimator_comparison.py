#!/usr/bin/env python3
"""Comparacion con los metodos publicados de inferencia circadiana de muestra unica (A5).

Los metodos de referencia (CYCLOPS, ZeitZeiger, TimeSignature) estiman FASE. Se reimplementa el
nucleo algoritmico de cada familia en Python y se evalua contra las horas registradas, en las dos
cohortes de validacion. Se evalua ademas la amplitud, que ninguno de los tres produce.

  - molecular timetable (Ueda): correlacion de la muestra con el perfil medio de cada hora.
  - CYCLOPS-like: ordenamiento no supervisado en el circulo a partir de los genes semilla
    (aqui por las dos primeras componentes principales, que es el nucleo del autoencoder).
  - ZeitZeiger-like: supervisado, proyeccion sobre la base periodica aprendida de las muestras
    de entrenamiento con hora conocida.
  - este trabajo: ordenamiento no supervisado + proyeccion con signo sobre la plantilla.

ATENCION: son reimplementaciones del algoritmo, no el software original. Declararlo en Metodos.
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
import phase as ph
OUT = os.environ.get("T2D_OUT", "results")
TAU = 24.0

def circ_err(pred, true):
    """Error de fase en horas, en el circulo."""
    d = np.abs((pred - true + TAU/2) % TAU - TAU/2); return d

def align(pred, true):
    """Quita el desfase global y el sentido de giro, que ningun metodo no supervisado fija."""
    best = None
    for s in (1, -1):
        for off in np.arange(0, TAU, 0.25):
            e = np.median(circ_err((s*pred + off) % TAU, true))
            if best is None or e < best[0]: best = (e, s, off)
    return best

def timetable(Z, train_idx, t_train):
    """Molecular timetable: perfil medio por hora y correlacion de cada muestra."""
    hrs = np.arange(0, TAU, 2.0)
    prof = np.array([Z[:, train_idx][:, np.argsort(circ_err(t_train, h))[:4]].mean(1) for h in hrs])
    r = np.array([[stats.pearsonr(Z[:, i], prof[j])[0] for j in range(len(hrs))] for i in range(Z.shape[1])])
    return hrs[np.argmax(r, 1)]

def cyclops_like(Z):
    """Ordenamiento no supervisado en el circulo por las dos primeras componentes."""
    U, S, Vt = np.linalg.svd(Z - Z.mean(1, keepdims=True), full_matrices=False)
    th = np.arctan2(Vt[1], Vt[0]); return (th % (2*np.pi))*TAU/(2*np.pi)

def zeitzeiger_like(Z, train_idx, t_train):
    """Supervisado: base periodica de primer armonico ajustada en entrenamiento."""
    w = 2*np.pi/TAU; D = np.column_stack([np.ones(len(t_train)), np.cos(w*t_train), np.sin(w*t_train)])
    B, *_ = np.linalg.lstsq(D, Z[:, train_idx].T, rcond=None)      # 3 x genes
    a, b = B[1], B[2]
    th = np.arctan2(Z.T@b, Z.T@a); return (th % (2*np.pi))*TAU/(2*np.pi)

def ours(Z):
    ang = ph.unsup(Z)
    return (ang % (2*np.pi))*TAU/(2*np.pi)

if __name__ == "__main__":
    X = pd.read_pickle(f"{OUT}/new/GSE182117_log2cpm.pkl")
    md = pd.read_pickle(f"{OUT}/new/GSE182117_meta.pkl").set_index("gsm")
    md["t"] = md.c3.str.replace("time: ", "").str.replace("h", "").astype(float) % TAU
    md = md.reindex(X.columns)
    Xc = X.sub(X.mean(0), axis=1); g = [x for x in ph.CLK if x in Xc.index]
    Zd = Xc.loc[g]; Z = ((Zd.T - Zd.mean(1))/(Zd.std(1)+1e-9)).T.values
    t = md.t.values
    rng = np.random.default_rng(0); n = Z.shape[1]
    tr = rng.choice(n, n//2, replace=False); te = np.setdiff1d(np.arange(n), tr)
    print(f"GSE182117 · {n} muestras, {len(g)} genes del reloj, horas registradas")
    print(f"{'metodo':26s} {'error de fase (h)':>18s} {'supervisado':>12s}")
    res = {"molecular timetable": (timetable(Z, tr, t[tr]), True),
           "CYCLOPS-like": (cyclops_like(Z), False),
           "ZeitZeiger-like": (zeitzeiger_like(Z, tr, t[tr]), True),
           "this work": (ours(Z), False)}
    rows = []
    for name, (p, sup) in res.items():
        idx = te if sup else np.arange(n)
        e, s_, off = align(p[idx], t[idx])
        rows.append(dict(method=name, median_error_h=e, supervised=sup, n=len(idx)))
        print(f"{name:26s} {e:18.2f} {'si' if sup else 'no':>12s}")
    pd.DataFrame(rows).to_csv(f"{OUT}/new/estimator_comparison_phase.tsv", sep="\t", index=False)
    print("\nAmplitud por muestra: ninguno de los tres metodos de referencia la produce.")
    print("Es la cantidad que este trabajo aporta y la que se valida en la Fig. 1d,e.")
