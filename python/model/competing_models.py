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

COMPARACION EMPAREJADA. Los cuatro modelos ven exactamente a las mismas personas (el mismo cortisol,
el mismo receptor), a las mismas horas de biopsia y con el mismo ruido de medida: esas realizaciones
se sortean una vez por escenario y semilla y se reutilizan en M0, M1, M2 y M3, de modo que lo unico
que difiere entre modelos es el mecanismo. Todo se repite con varias semillas.

CUANTAS VECES, NO SI ALGUNA VEZ. Preguntar si un modelo puede producir lo observado "en algun punto
de la rejilla" es demasiado facil: con suficientes semillas el ruido de muestreo acaba rozando
cualquier intervalo. La cifra principal es por eso una FRACCION: de todas las simulaciones de un
modelo (rejilla por semillas), en cuantas cae cada observable dentro del intervalo de confianza
observado, y en cuantas caen TODOS a la vez (competing_models_joint.tsv). Un mecanismo que no
reproduce el conjunto en ninguna simulacion queda descartado; dos que lo reproducen con frecuencia
parecida quedan sin distinguir.

TRES ESTADOS, NO DOS. Cada observable es positivo, negativo o compatible con cero. En los datos se
decide con el intervalo de confianza al 95 % por bootstrap de personas; en las simulaciones, con el
umbral de muestreo del tamano de cohorte simulado. Un modelo reproduce un rasgo si puede generar el
mismo estado Y si el intervalo de lo que simula se solapa con el intervalo observado.

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
N_DONORS = int(os.environ.get("T2D_MODEL_N", 200))
SEEDS = list(range(int(os.environ.get("T2D_MODEL_SEEDS", 3))))
NBOOT = int(os.environ.get("T2D_MODEL_BOOT", 500))
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
        v = (A.loc["NR3C1"].values, M.loc["NR3C1"].values, score(A, CL).values, score(M, CL).values)
        o = observables(*v, C); rng = np.random.default_rng(0); n = len(v[0]); B = {k: [] for k in OBS}
        for _ in range(NBOOT):                               # se remuestrean personas: los dos tejidos van juntos
            i = rng.integers(0, n, n)
            if np.linalg.matrix_rank(C[i]) < C.shape[1]: continue
            b = observables(*(x[i] for x in v), C[i])
            for k in OBS: B[k].append(b[k])
        ci = {}
        for k in OBS:
            ci[f"{k} lo"], ci[f"{k} hi"] = np.percentile(B[k], [2.5, 97.5])
        rows.append(dict(cohort=label, window=window, **o, **ci))
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
def people(sd_F, rho_K, seed, n=N_DONORS, sd_K=0.35):
    """Las personas de un escenario: su cortisol y su receptor en cada tejido. Las mismas para todos
    los modelos."""
    rng = np.random.default_rng([seed, int(sd_F * 100), int(rho_K * 100)])
    F = np.clip(F_MEAN * (1 + sd_F * rng.normal(0, 1, n)), 0.0, None)
    shared = rng.normal(0, 1, n)                                             # parte del receptor comun a los dos tejidos
    K = {t: np.clip(K_MEAN + sd_K * (np.sqrt(rho_K) * shared + np.sqrt(1 - rho_K) * rng.normal(0, 1, n)), 0.05, None) for t in "AM"}
    cell_seed = {t: rng.integers(1e6, size=n) for t in "AM"}                 # la misma realizacion celular en cada modelo
    return F, K, cell_seed


def ensemble_of(task):
    """Amplitud y fase del ensamble de cada persona y tejido bajo un modelo. No depende de la hora de
    biopsia ni del ruido de medida, que se aplican despues."""
    model, sd_F, rho_K, seed = task
    F, K, cell_seed = people(sd_F, rho_K, seed)
    amp, ph = {}, {}
    for t in "AM":
        fk = MODELS[model](F, K[t], F_MEAN, K_MEAN)
        res = [fc.ensemble(*fc.simulate_population(1.0, float(x), n_cells=40, days=12, seed=int(sd))) for x, sd in zip(fk, cell_seed[t])]
        amp[t] = np.array([r[0] for r in res]); ph[t] = np.array([r[1] for r in res])
    return task, K, amp, ph


def measure(K, amp, ph, hour, proxy_eps, expr_seed, proxy_noise, expr_noise=0.6):
    """Lo que se mide: una puntuacion de salida a la hora de biopsia y un sustituto ruidoso del receptor.
    La hora, el ruido de expresion y el del receptor llegan ya sorteados, iguales para todos los modelos."""
    out = {}
    for t in "AM":
        d = pd.DataFrame(dict(subject=np.arange(len(hour)), hour=hour, amp_ens=amp[t], phase_ens=ph[t]))
        E = fc.synthetic_expression(d, noise=expr_noise, seed=int(expr_seed[t]))
        Z = E.loc[OUTPUT_GENES]; out[t] = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.mean(0).values
    rec = {t: K[t] + proxy_noise * np.std(K[t]) * proxy_eps[t] for t in "AM"}
    return rec, out


def simulate():
    from concurrent.futures import ProcessPoolExecutor
    scen = sorted({(g["sd_F"], g["rho_K"]) for g in GRID})
    tasks = [(m, sf, rk, sd) for m in MODELS for (sf, rk) in scen for sd in SEEDS
             if not (m == "M0 no forcing" and (sf, rk) != scen[0])]          # sin forzamiento nada de eso actua
    workers = max(1, min(int(os.environ.get("T2D_MODEL_WORKERS", os.cpu_count() or 1)), 16))
    print(f"  {len(tasks)} cohortes de {N_DONORS} personas, {workers} procesos", flush=True)
    with ProcessPoolExecutor(workers) as ex:
        sims = list(ex.map(ensemble_of, tasks))
    rows = []
    for (model, sd_F, rho_K, seed), K, amp, ph in sims:
        # las mismas horas y el mismo ruido para todos los modelos de este escenario y esta semilla
        rng = np.random.default_rng([1000 + seed, int(sd_F * 100), int(rho_K * 100)])
        n = len(amp["A"])
        hours = {"morning biopsy": rng.normal(8.0, 1.5, n) % 24, "any hour": rng.uniform(0, 24, n)}
        eps = {t: rng.normal(0, 1, n) for t in "AM"}
        eseed = {w: {t: rng.integers(1e6) for t in "AM"} for w in hours}
        for pn in sorted({g["proxy_noise"] for g in GRID}):
            for window, hour in hours.items():
                rec, out = measure(K, amp, ph, hour, eps, eseed[window], pn)
                rows.append(dict(model=model, sd_F=sd_F, rho_K=rho_K, proxy_noise=pn, seed=seed, window=window,
                                 **observables(rec["A"], rec["M"], out["A"], out["M"])))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ el veredicto
def state(lo, hi):
    return "positive" if lo > 0 else "negative" if hi < 0 else "zero"


def verdict(O, S):
    """Para cada rasgo observado: en que estado esta (positivo, negativo o compatible con cero, segun
    su intervalo de confianza) y si cada modelo puede generar ese estado con un intervalo que se
    solape con el observado, en algun punto de la rejilla y alguna semilla."""
    thr = 1.96 / np.sqrt(N_DONORS - 3)                                       # |rho| distinguible de cero con el n simulado
    rows = []
    for win, wsim in (("morning biopsy", "morning biopsy"), ("post mortem, any hour", "any hour")):
        o = O[O.window == win]
        if not len(o): continue
        for k in OBS:
            if k.startswith("O3"): continue                                  # se usa para fijar el modelo, no para puntuarlo
            states = [state(r[f"{k} lo"], r[f"{k} hi"]) for _, r in o.iterrows()]
            obs_state = states[0] if len(set(states)) == 1 else "mixed: " + ", ".join(sorted(set(states)))
            lo, hi = o[f"{k} lo"].min(), o[f"{k} hi"].max()
            for model in MODELS:
                v = S[(S.model == model) & (S.window == wsim)][k].values
                can = {"positive": bool((v > thr).any()), "negative": bool((v < -thr).any()), "zero": bool((np.abs(v) <= thr).any())}
                state_ok = all(can[x] for x in set(states))
                overlap = bool(v.max() >= lo and v.min() <= hi)
                inside = float(np.mean((v >= lo) & (v <= hi)))
                rows.append(dict(window=win, observable=k, observed_state=obs_state, observed_median=o[k].median(),
                                 fraction_of_simulations_within_observed_ci=inside,
                                 observed_ci_lo=lo, observed_ci_hi=hi, model=model, simulated_min=v.min(), simulated_max=v.max(),
                                 model_can_be_positive=can["positive"], model_can_be_negative=can["negative"], model_can_be_zero=can["zero"],
                                 state_reproduced=state_ok, interval_overlaps=overlap, model_can_reproduce=state_ok and overlap))
    return pd.DataFrame(rows)


def joint(O, S):
    """Fraccion de las simulaciones de cada modelo en la que TODOS los observables puntuados caen a la
    vez dentro del intervalo de confianza observado."""
    keys = [k for k in OBS if not k.startswith("O3")]; rows = []
    for win, wsim in (("morning biopsy", "morning biopsy"), ("post mortem, any hour", "any hour")):
        o = O[O.window == win]
        if not len(o): continue
        for model in MODELS:
            q = S[(S.model == model) & (S.window == wsim)]
            ok = np.all([(q[k].values >= o[f"{k} lo"].min()) & (q[k].values <= o[f"{k} hi"].max()) for k in keys], axis=0)
            rows.append(dict(window=win, model=model, n_simulations=len(q), n_all_within_ci=int(ok.sum()),
                             fraction_all_within_ci=float(ok.mean()),
                             **{f"fraction {k.split()[0]}": float(np.mean((q[k].values >= o[f"{k} lo"].min()) & (q[k].values <= o[f"{k} hi"].max()))) for k in keys}))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    pd.set_option("display.width", 220)
    O = observed()
    if len(O):
        O.to_csv(f"{OUT}/new/competing_models_observed.tsv", sep="\t", index=False)
        print("OBSERVADO, cohortes con dos tejidos por persona (IC 95 % por bootstrap de personas):")
        for _, r in O.iterrows():
            print(f"  {r.cohort} (n = {int(r.n)}, {r.window}): " + "; ".join(f"{k.split()[0]} {r[k]:+.2f} [{r[k + ' lo']:+.2f}, {r[k + ' hi']:+.2f}]" for k in OBS))
    print(f"\nSIMULANDO {len(MODELS)} modelos, {len(SEEDS)} semillas:")
    fs = f"{OUT}/new/competing_models_simulated.tsv"
    if os.environ.get("T2D_MODEL_REUSE") == "1" and os.path.exists(fs):
        S = pd.read_csv(fs, sep="\t"); print(f"  (se reutilizan las {len(S)} simulaciones de {fs})")
    else:
        S = simulate(); S.to_csv(fs, sep="\t", index=False)
    print("\nSIMULADO, intervalo de cada observable al recorrer lo que no se conoce:")
    for win in ("morning biopsy", "any hour"):
        print(f"  [{win}]")
        for model in MODELS:
            q = S[(S.model == model) & (S.window == win)]
            print(f"    {model:24s} " + "  ".join(f"{k.split()[0]} {q[k].min():+.2f}..{q[k].max():+.2f}" for k in OBS))
    if len(O):
        V = verdict(O, S); V.to_csv(f"{OUT}/new/competing_models_verdict.tsv", sep="\t", index=False)
        J = joint(O, S); J.to_csv(f"{OUT}/new/competing_models_joint.tsv", sep="\t", index=False)
        print("\nFRACCION de las simulaciones de cada modelo dentro del intervalo de confianza observado")
        for win in J.window.unique():
            print(f"  [{win}]")
            for r in J[J.window == win].itertuples():
                per = "  ".join(f"{c.split()[1]} {J.loc[r.Index, c]:.2f}" for c in J.columns if c.startswith("fraction O"))
                print(f"    {r.model:24s} {per}   | todos a la vez: {r.n_all_within_ci} de {r.n_simulations} ({r.fraction_all_within_ci:.2f})")
        print("\nESTADO: rasgos observados cuyo signo o intervalo cada modelo no alcanza en ninguna simulacion")
        for win in V.window.unique():
            print(f"  [{win}]")
            for model in MODELS:
                bad = V[(V.model == model) & (V.window == win) & ~V.model_can_reproduce]
                print(f"    {model:24s} " + ("los produce todos" if not len(bad) else
                      "; ".join(f"{r.observable.split()[0]} (observado {r.observed_state}, IC {r.observed_ci_lo:+.2f} a {r.observed_ci_hi:+.2f}; "
                                f"modelo {r.simulated_min:+.2f} a {r.simulated_max:+.2f})" for r in bad.itertuples())))
