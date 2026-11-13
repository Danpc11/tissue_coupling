#!/usr/bin/env python3
"""Cuatro mecanismos compitiendo contra los mismos datos.

El modelo del articulo lleva dentro lo que quiere explicar: el forzamiento efectivo es F*K, el
ritmo de cortisol por la disponibilidad del receptor. Que al bajar K baje la salida del reloj no
prueba nada, porque asi esta escrito. Para saber que aporta ese producto hay que ponerlo a competir
con mecanismos mas simples, con el MISMO oscilador, el mismo ruido, el mismo muestreo y el mismo
calculo de los observables, cambiando solo de que depende el forzamiento de cada tejido:

    M0  sin forzamiento externo            0
    M1  solo el ritmo de cortisol          F_i                (el receptor varia pero no actua)
    M2  cortisol por receptor              F_i * K_it         (el del articulo)
    M3  solo el receptor                   Fmedio * K_it      (el mismo ritmo en todas las personas)

F_i es de la persona y lo comparten sus tejidos; K_it es de cada tejido. Cada persona aporta dos
tejidos y una hora de biopsia. De cada cohorte simulada se calculan los mismos observables que de
los datos, sobre la puntuacion de los genes de salida a la hora de la biopsia y sobre una medida
ruidosa del receptor, que es lo que da el ARNm:

    O1  receptor ~ salida, dentro del tejido
    O2  salida de un tejido ~ salida del otro, misma persona
    O3  receptor de un tejido ~ receptor del otro            (se usa para fijar el modelo, no para puntuarlo)
    O4  receptor de un tejido ~ salida del OTRO
    O5  receptor ~ salida del mismo tejido, descontado el receptor del otro
    O6  salida ~ salida entre tejidos, descontados los dos receptores

QUE DISTINGUE A CADA MODELO. M0 no produce ninguna asociacion. M1 produce concordancia entre tejidos
(O2) pero ninguna relacion con el receptor (O1 = 0). M3 produce O1, y concordancia entre tejidos
solo en la medida en que el receptor de los dos tejidos se parezca. M2 produce las dos cosas.
La pregunta que los datos pueden contestar no es "es cierto M2" sino "que modelos NO pueden producir
lo observado, con ningun valor razonable de lo que no se conoce" (el ruido de la medida del
receptor, cuanto varia el cortisol entre personas, cuanto se parece el receptor entre tejidos).
Por eso cada modelo se simula en una rejilla de esos valores y se informa el intervalo de cada
observable, no un punto.

LO QUE ESTO NO HACE. No ajusta ningun parametro a los datos y no elige un modelo por verosimilitud.
Dice que rasgos cualitativos de los datos puede y no puede generar cada mecanismo. Dos mecanismos
que generan los mismos rasgos quedan sin distinguir, y eso se informa asi.

Entradas (si existen): results/new/AAGMEx_*_log2.pkl, GTEx_*.pkl, GTEx_pheno.pkl y la matriz de serie
                       de AAGMEx. Sin ellas se escriben solo las simulaciones.
Salidas:  results/new/competing_models_observed.tsv     los observables en cada cohorte con dos tejidos
          results/new/competing_models_simulated.tsv    lo mismo en cada modelo y cada punto de la rejilla
          results/new/competing_models_verdict.tsv      que rasgos observados puede producir cada modelo
"""
import os, sys, gzip, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "../lib"))
import forced_clock as fc
from gene_sets import SETS, resolve
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results")

MODELS = {"M0 no forcing": lambda F, K, Fm, Km: 0.0 * F,
          "M1 cortisol only": lambda F, K, Fm, Km: F * Km,
          "M2 cortisol x receptor": lambda F, K, Fm, Km: F * K,
          "M3 receptor only": lambda F, K, Fm, Km: Fm * K}
OBS = ["O1 receptor~output", "O2 output~output", "O3 receptor~receptor", "O4 receptor~other output",
       "O5 receptor~output | other receptor", "O6 output~output | receptors"]
N_DONORS = int(os.environ.get("T2D_MODEL_N", 240))
F_MEAN, K_MEAN = 0.04, 1.0                      # los del articulo (make_fig7_data.py)
# lo que no se conoce: se recorre, no se ajusta
GRID = [dict(sd_F=sf, rho_K=rk, proxy_noise=pn) for sf in (0.25, 0.5) for rk in (0.0, 0.3, 0.6) for pn in (0.5, 1.0, 2.0)]
OUTPUT_GENES = list(SETS["Clock output"]["genes"])        # registro unico


def resid(y, C):
    b, *_ = np.linalg.lstsq(C, y, rcond=None); return y - C @ b


def observables(rA, rM, oA, oM, C=None):
    """Los seis observables de una cohorte con dos tejidos por persona. rA, rM: receptor en cada
    tejido; oA, oM: puntuacion de salida. C: covariables (con intercepto) que se descuentan antes."""
    n = len(rA); C = np.ones((n, 1)) if C is None else C
    rA, rM, oA, oM = (resid(np.asarray(v, float), C) for v in (rA, rM, oA, oM))
    sp = lambda a, b: stats.spearmanr(a, b).statistic
    par = lambda a, b, Z: sp(resid(a, np.column_stack([np.ones(n), Z])), resid(b, np.column_stack([np.ones(n), Z])))
    return {"O1 receptor~output": (sp(rA, oA) + sp(rM, oM)) / 2,
            "O2 output~output": sp(oA, oM),
            "O3 receptor~receptor": sp(rA, rM),
            "O4 receptor~other output": (sp(rA, oM) + sp(rM, oA)) / 2,
            "O5 receptor~output | other receptor": (par(rA, oA, rM) + par(rM, oM, rA)) / 2,
            "O6 output~output | receptors": par(oA, oM, np.column_stack([rA, rM])), "n": n}


# ------------------------------------------------------------------ los datos
def score(X, genes):
    g = [r for r in (resolve(x, X.index) for x in genes) if r]
    Z = X.loc[g]; return ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.mean(0)


def observed():
    rows = []; CL = SETS["Clock output"]["genes"]
    def add(label, window, A, M, C):
        A = A.sub(A.mean(0), axis=1); M = M.sub(M.mean(0), axis=1)          # centrado por muestra, como en el metaanalisis
        o = observables(A.loc["NR3C1"].values, M.loc["NR3C1"].values, score(A, CL).values, score(M, CL).values, C)
        rows.append(dict(cohort=label, window=window, **o))
    fa, fm = f"{OUT}/new/AAGMEx_adipose_log2.pkl", f"{OUT}/new/AAGMEx_muscle_log2.pkl"
    if os.path.exists(fa) and os.path.exists(fm):
        A, M = pd.read_pickle(fa), pd.read_pickle(fm); ids = [c for c in A.columns if c in set(M.columns)]
        t = sx = ag = None
        with gzip.open(f"{RAW}/GSE95674_series_matrix.txt.gz", "rt") as fh:
            for l in fh:
                v = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
                if l.startswith("!Sample_title"): t = v
                if l.startswith("!Sample_characteristics") and "gender" in l: sx = [x.split(": ")[1] for x in v]
                if l.startswith("!Sample_characteristics") and "age" in l: ag = [float(x.split(": ")[1]) for x in v]
        cov = pd.DataFrame({"sex": sx, "age": ag}, index=t).reindex(ids)
        add("AAGMEx adipose-muscle", "morning biopsy", A[ids], M[ids],
            np.column_stack([np.ones(len(ids)), (cov.sex == "M").astype(float), cov.age]))
    fp = f"{OUT}/new/GTEx_pheno.pkl"
    if os.path.exists(fp):
        P = pd.read_pickle(fp); agemid = {"20-29": 25, "30-39": 35, "40-49": 45, "50-59": 55, "60-69": 65, "70-79": 75}
        G = {k: pd.read_pickle(f"{OUT}/new/GTEx_{k}.pkl") for k in ("adipose_sc", "adipose_visc", "muscle")
             if os.path.exists(f"{OUT}/new/GTEx_{k}.pkl")}
        for a, b in (("adipose_sc", "muscle"), ("adipose_visc", "muscle"), ("adipose_sc", "adipose_visc")):
            if a in G and b in G:
                ids = sorted(set(G[a].columns) & set(G[b].columns) & set(P.index)); p = P.loc[ids]
                H = pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int), drop_first=True).values.astype(float)
                C = np.column_stack([np.ones(len(ids)), (p.SEX == 1).astype(float), p.AGE.map(agemid).astype(float), H])
                add(f"GTEx {a}-{b}", "post mortem, any hour", G[a][ids], G[b][ids], C)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ las simulaciones
def cohort(model, sd_F, rho_K, seed, n=N_DONORS, sd_K=0.35):
    """Una cohorte de n personas con dos tejidos. Devuelve por persona y tejido la amplitud y la fase
    del ensamble; no depende de la hora de biopsia ni del ruido de medida, que se aplican despues."""
    rng = np.random.default_rng(seed)
    F = np.clip(F_MEAN * (1 + sd_F * rng.normal(0, 1, n)), 0.0, None)
    shared = rng.normal(0, 1, n)                                             # parte del receptor comun a los dos tejidos
    K = {t: np.clip(K_MEAN + sd_K * (np.sqrt(rho_K) * shared + np.sqrt(1 - rho_K) * rng.normal(0, 1, n)), 0.05, None) for t in "AM"}
    amp, ph = {}, {}
    for t in "AM":
        fk = MODELS[model](F, K[t], F_MEAN, K_MEAN)
        res = [fc.ensemble(*fc.simulate_population(1.0, float(x), n_cells=40, days=12, seed=int(rng.integers(1e6)))) for x in fk]
        amp[t] = np.array([r[0] for r in res]); ph[t] = np.array([r[1] for r in res])
    return F, K, amp, ph


def measure(K, amp, ph, hour, proxy_noise, expr_noise, rng):
    """Lo que se mide: una puntuacion de salida a la hora de biopsia y un sustituto ruidoso del receptor."""
    out = {}
    for t in "AM":
        d = pd.DataFrame(dict(subject=np.arange(len(hour)), hour=hour, amp_ens=amp[t], phase_ens=ph[t]))
        E = fc.synthetic_expression(d, noise=expr_noise, seed=int(rng.integers(1e6)))
        Z = E.loc[OUTPUT_GENES]; out[t] = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.mean(0).values
    rec = {t: K[t] + proxy_noise * np.std(K[t]) * rng.normal(0, 1, len(hour)) for t in "AM"}
    return rec, out


def simulate():
    rows = []
    sims = {}
    for model in MODELS:
        for sd_F in sorted({g["sd_F"] for g in GRID}):
            for rho_K in sorted({g["rho_K"] for g in GRID}):
                if model == "M0 no forcing" and (sd_F, rho_K) != (0.25, 0.0): continue      # sin forzamiento nada de eso actua
                sims[(model, sd_F, rho_K)] = cohort(model, sd_F, rho_K, seed=hash((sd_F, rho_K)) % 1000)
        print(f"  simulado {model}", flush=True)
    rng = np.random.default_rng(0)
    for (model, sd_F, rho_K), (F, K, amp, ph) in sims.items():
        for pn in sorted({g["proxy_noise"] for g in GRID}):
            for window, hour in (("morning biopsy", rng.normal(8.0, 1.5, len(F)) % 24), ("any hour", rng.uniform(0, 24, len(F)))):
                rec, out = measure(K, amp, ph, hour, pn, 0.6, rng)
                rows.append(dict(model=model, sd_F=sd_F, rho_K=rho_K, proxy_noise=pn, window=window,
                                 **observables(rec["A"], rec["M"], out["A"], out["M"])))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ el veredicto
def verdict(O, S):
    """Para cada rasgo observado, puede el modelo producirlo en algun punto de la rejilla?"""
    crit = lambda n: 1.96 / np.sqrt(max(n - 3, 1))                           # |rho| que seria significativa con ese n
    feats = []
    for win, wsim in (("morning biopsy", "morning biopsy"), ("post mortem, any hour", "any hour")):
        o = O[O.window == win]
        if not len(o): continue
        for k in OBS:
            if k.startswith("O3"): continue
            feats.append(dict(window=win, sim_window=wsim, observable=k, observed_median=o[k].median(),
                              observed_min=o[k].min(), observed_max=o[k].max(),
                              present=bool((o[k] > o.n.map(crit)).mean() >= 0.5)))
    rows = []
    for f in feats:
        for model in MODELS:
            s = S[(S.model == model) & (S.window == f["sim_window"])][f["observable"]]
            thr = crit(N_DONORS)
            can_present = bool((s > thr).any()); can_absent = bool((s.abs() <= thr).any())
            ok = can_present if f["present"] else can_absent
            rows.append(dict(window=f["window"], observable=f["observable"], observed_median=f["observed_median"],
                             observed_present=f["present"], model=model, simulated_min=s.min(), simulated_max=s.max(),
                             model_can_reproduce=ok))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    pd.set_option("display.width", 220)
    O = observed()
    if len(O):
        O.to_csv(f"{OUT}/new/competing_models_observed.tsv", sep="\t", index=False)
        print("OBSERVADO, cohortes con dos tejidos por persona:"); print(O.round(3).to_string(index=False))
    print(f"\nSIMULANDO {len(MODELS)} modelos, {N_DONORS} personas y dos tejidos cada uno:")
    S = simulate(); S.to_csv(f"{OUT}/new/competing_models_simulated.tsv", sep="\t", index=False)
    print("\nSIMULADO, intervalo de cada observable al recorrer lo que no se conoce:")
    for win in ("morning biopsy", "any hour"):
        print(f"  [{win}]")
        for model in MODELS:
            q = S[(S.model == model) & (S.window == win)]
            print(f"    {model:24s} " + "  ".join(f"{k.split()[0]} {q[k].min():+.2f}..{q[k].max():+.2f}" for k in OBS))
    if len(O):
        V = verdict(O, S); V.to_csv(f"{OUT}/new/competing_models_verdict.tsv", sep="\t", index=False)
        print("\nVEREDICTO: rasgos observados que cada modelo NO puede producir en ningun punto de la rejilla")
        for model in MODELS:
            bad = V[(V.model == model) & ~V.model_can_reproduce]
            print(f"  {model:24s} " + ("los produce todos" if not len(bad) else
                  "; ".join(f"{r.observable.split()[0]} en {r.window.split(',')[0]} (observado {r.observed_median:+.2f}, modelo {r.simulated_min:+.2f}..{r.simulated_max:+.2f})"
                            for r in bad.itertuples())))
