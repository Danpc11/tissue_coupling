#!/usr/bin/env python3
"""Concordancia con CHIRAL: la coherencia de fase entre tejidos no depende del estimador.

CHIRAL (Talamanca, Gobet y Naef, Science 2023) infiere fase circadiana de muestra unica. Se calcula
con el, y con el estimador de este trabajo, la misma cantidad: R, la concentracion de la diferencia
de fase entre dos tejidos del mismo donante.

    d_i = arg( exp(i (phi_musculo,i - phi_adiposo,i)) )      diferencia envuelta al circulo
    R   = | mean_i exp(i d_i) |                              0 sin relacion de fase, 1 rigida

La significancia se evalua permutando la correspondencia de donantes entre tejidos (2.000 veces).
El argumento no es que un estimador sea mejor, sino que el resultado es el mismo con los dos.

QUE MIDE R Y QUE NO. R mide coherencia de fase entre los dos tejidos de una persona: que la
diferencia de fase se parezca de un donante a otro. No distingue una interaccion directa entre
tejidos de una senal sistemica que ambos siguen; las dos dan R alto (tests/test_phase_statistic.py
lo muestra con fases simuladas). La permutacion contrasta contra la asignacion al azar de
donantes, no contra una hipotesis causal.

SENTIDO DE GIRO. CHIRAL no tiene ancla: fija la fase salvo rotacion y sentido, por separado en
cada tejido, y el sentido que sale cambia con el orden de las muestras. R, tal como se define
arriba, supone que los dos tejidos giran en el mismo sentido; si no, da cerca de cero aunque esten
alineados. La tabla lleva por eso una segunda columna, R_reflection_aware, que no lo supone. Cual
de las dos se cita para CHIRAL es una decision del articulo; para el estimador propio, que si
esta anclado a la plantilla, las dos coinciden.

CHIRAL ES OPCIONAL Y SE DICE. Si la implementacion no esta (ver python/lib/chiral_backend.py), el
paso calcula solo el estimador propio, lo escribe en la tabla con la fila de CHIRAL marcada como
no ejecutada y lo anuncia al final. Con T2D_REQUIRE_CHIRAL=1 la ausencia detiene la corrida.

AVISO del modulo: CHIRAL espera las muestras en FILAS. La envoltura transpone antes de llamar.
"""
import os, sys, warnings, collections
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
import phase as ph
OUT = os.environ.get("T2D_OUT", "results")

from gene_sets import panel as _panel
CLK = _panel("CHIRAL input")          # registro unico; ver alli por que conserva ARNTL

def concentration(a, b, seed=0, B=2000):
    """R de la diferencia de fase y su P por permutacion de donantes."""
    d = np.angle(np.exp(1j*(np.asarray(a) - np.asarray(b))))
    R = np.abs(np.mean(np.exp(1j*d)))
    rng = np.random.default_rng(seed)
    null = np.array([np.abs(np.mean(np.exp(1j*(a - np.asarray(b)[rng.permutation(len(b))]))))
                     for _ in range(B)])
    return R, (np.sum(null >= R) + 1)/(B + 1)

def concentration_reflection_aware(a, b, seed=0, B=2000):
    """La misma concentracion cuando el SENTIDO DE GIRO de cada serie de fases es arbitrario.

    Un estimador sin ancla, como CHIRAL, fija la fase salvo rotacion y sentido, y lo hace por
    separado en cada tejido. Si en un tejido el ciclo queda recorrido al reves, la diferencia
    phi_a - phi_b deja de concentrarse aunque los dos tejidos esten perfectamente alineados, y lo
    que se concentra es la suma. El sentido que sale depende de detalles sin significado, como el
    orden de las muestras. Aqui se toma el mayor de los dos y el nulo usa el mismo maximo, de modo
    que el P sigue siendo valido."""
    a, b = np.asarray(a), np.asarray(b)
    stat = lambda x, y: max(np.abs(np.mean(np.exp(1j*(x - y)))), np.abs(np.mean(np.exp(1j*(x + y)))))
    R = stat(a, b)
    rng = np.random.default_rng(seed)
    null = np.array([stat(a, b[rng.permutation(len(b))]) for _ in range(B)])
    return R, (np.sum(null >= R) + 1)/(B + 1)

def phase_ours(E, genes=None):
    """Fase con el periodo ANCLADO a la plantilla, no con ordenamiento libre.

    El ordenamiento libre da bien la amplitud pero no fija el periodo, y la diferencia de fase
    entre tejidos queda en escalas distintas. Para fase va
    el anclado, con los 60 genes mas ritmicos de la plantilla.
    """
    T = pd.read_csv(f"{OUT}/template_2harm.tsv", sep="\t", index_col=0)
    # order_anchored espera columnas r2, phase_h y amp; la plantilla guarda armonicos
    T = T.assign(phase_h=(np.arctan2(T.b1, T.a1) % (2*np.pi))*24/(2*np.pi), amp=T.amp1)
    return ph.order_anchored(E, T)[0].values   # devuelve (fases, genes usados)

CHIRAL_WARNINGS = collections.Counter()

def phase_chiral(E, genes):
    """CHIRAL. La envoltura fija la orientacion (genes x muestras de entrada). Los avisos que emita
    no se silencian: se cuentan por tipo y se informan al final."""
    from chiral_backend import infer_phase
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        phi = infer_phase(E.loc[genes], genes)["phi"]
    CHIRAL_WARNINGS.update(f"{x.category.__name__}: {str(x.message)[:70]}" for x in w)
    return phi

def run(Em, Ea, donors, label, use_chiral=True):
    g = [x for x in CLK if x in Em.index and x in Ea.index]
    rows = []
    methods = [("this work", lambda E, gg: phase_ours(E))] + ([("CHIRAL", phase_chiral)] if use_chiral else [])
    for name, f in methods:
        base = dict(cohort=label, method=name, n_donors=len(donors), n_genes=len(g))
        try:
            pm, pa = f(Em[donors], g), f(Ea[donors], g)
            R, p = concentration(pm, pa)
            R2, p2 = concentration_reflection_aware(pm, pa)
            rows.append(dict(base, R=R, P=p, R_reflection_aware=R2, P_reflection_aware=p2, status="ok"))
            print(f"  {label:10s} {name:10s} n={len(donors):4d}  R = {R:.3f}  P = {p:.4f}"
                  f"   | sin fijar el sentido de giro: R = {R2:.3f}  P = {p2:.4f}")
        except ImportError as e:
            rows.append(dict(base, R=np.nan, P=np.nan, R_reflection_aware=np.nan, P_reflection_aware=np.nan,
                             status=f"not run: {e}"))
            print(f"  {label:10s} {name:10s} NO EJECUTADO: {e}")
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

    import chiral_backend
    impl = chiral_backend.implementation()
    T = pd.DataFrame(out)
    T["implementation"] = np.where(T.method == "CHIRAL", impl or "", "python/lib/phase.py")
    T.to_csv(f"{OUT}/new/chiral_concordance.tsv", sep="\t", index=False)

    ran = int(((T.method == "CHIRAL") & (T.status == "ok")).sum())
    if ran:
        print(f"\nCHIRAL: EJECUTADO en {ran} cohortes con {impl}")
        for k, n in CHIRAL_WARNINGS.most_common(5):
            print(f"  aviso de CHIRAL x{n}: {k}")
    else:
        print("\nCHIRAL: NO EJECUTADO. La tabla lleva solo el estimador de este trabajo; la concordancia "
              "con CHIRAL que cita el articulo no se ha regenerado en esta corrida.")
        if os.environ.get("T2D_REQUIRE_CHIRAL") == "1":
            raise SystemExit("T2D_REQUIRE_CHIRAL=1 y CHIRAL no esta disponible")
