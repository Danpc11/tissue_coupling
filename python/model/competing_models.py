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
    O3  receptor de un tejido ~ receptor del otro            (un dato que el modelo tambien tiene que cumplir)
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

EL CRITERIO CONJUNTO. Un modelo reproduce una cohorte si genera los SEIS observables a la vez, y eso
se decide asi:

  - Cohorte por cohorte. Cada simulacion se compara con una cohorte concreta, con sus estimaciones
    y su incertidumbre. No se mezclan cohortes: el limite inferior de una y el superior de otra no
    forman un intervalo de nada.
  - Con una region conjunta, no con seis intervalos sueltos. Los observables estan correlacionados
    (el receptor frente a la salida, con y sin descontar el otro receptor, son casi la misma
    cantidad), y seis intervalos marginales al 95 % no cubren al 95 % en conjunto. La distancia de
    cada simulacion a lo observado se mide con la covarianza de los seis, estimada por bootstrap de
    personas en la escala z de Fisher, a la que se suma la del muestreo de la propia simulacion
    (la misma covarianza, escalada por n observado / n simulado). Una simulacion es compatible si
    esa distancia queda dentro del cuantil 95 % de una ji cuadrado con seis grados de libertad.
  - Con la concordancia del receptor entre tejidos (O3) dentro. Es un dato que el modelo tambien
    tiene que cumplir: una simulacion en la que el receptor de los dos tejidos no se parece nada no
    describe una cohorte en la que si se parecen.
  - La prediccion de un modelo en un punto de la rejilla es la MEDIA de sus semillas, no cada cohorte
    simulada por separado. Una cohorte simulada de 200 personas es una realizacion ruidosa de lo
    que el modelo predice; promediar las semillas acerca a la prediccion, y la covarianza de la
    simulacion se reduce en consecuencia (n simulado = personas por semillas). Asi, anadir semillas
    hace el criterio mas exigente y no mas laxo: al reves que preguntar si alguna simulacion roza
    lo observado.
  - Contando. La cifra es la fraccion de los puntos de la rejilla de cada modelo compatibles con
    cada cohorte, con la distancia minima alcanzada y el punto que la alcanza.

Como comprobacion mas facil de leer se da tambien la fraccion de puntos con los seis observables
dentro de intervalos SIMULTANEOS al 95 % (los del maximo estudentizado del bootstrap), que no
incluyen el ruido de muestreo de la simulacion.

TRES ESTADOS POR OBSERVABLE. competing_models_verdict.tsv da, por cohorte y observable, si lo
observado es positivo, negativo o compatible con cero segun su intervalo, y que fraccion de los
puntos de la rejilla de cada modelo predice un valor dentro de ese intervalo.

LO QUE ESTO NO HACE. No ajusta ningun parametro a los datos y no elige un modelo por verosimilitud.
Dice que rasgos cualitativos de los datos puede y no puede generar cada mecanismo. Dos mecanismos
que generan los mismos rasgos quedan sin distinguir, y eso se informa asi.

Entradas (si existen): results/new/AAGMEx_*_log2.pkl, GTEx_*.pkl, GTEx_pheno.pkl y la matriz de serie
                       de AAGMEx. Sin ellas se escriben solo las simulaciones.
Salidas:  results/new/competing_models_observed.tsv     los observables en cada cohorte con dos tejidos
          results/new/competing_models_simulated.tsv    lo mismo en cada modelo y cada punto de la rejilla
          results/new/competing_models_joint.tsv        el criterio conjunto: simulaciones compatibles con cada cohorte
          results/new/competing_models_verdict.tsv      observable por observable, estado y fraccion dentro del intervalo
          results/new/competing_models_observed_covariance.tsv  covarianza de los seis observables de cada cohorte
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


def fz(x):
    return np.arctanh(np.clip(np.asarray(x, float), -0.999, 0.999))


def observed():
    """Los seis observables de cada cohorte con dos tejidos por persona, con su incertidumbre por
    bootstrap de personas. Devuelve la tabla y, por cohorte, la covarianza de los seis en z de Fisher."""
    rows = []; COV = {}; CL = SETS["Clock output"]["genes"]
    def add(label, window, A, M, C):
        A = A.sub(A.mean(0), axis=1); M = M.sub(M.mean(0), axis=1)          # centrado por muestra, como en el metaanalisis
        v = (A.loc["NR3C1"].values, M.loc["NR3C1"].values, score(A, CL).values, score(M, CL).values)
        o = observables(*v, C); rng = np.random.default_rng(0); n = len(v[0]); B = []
        for _ in range(NBOOT):                               # se remuestrean personas: los dos tejidos van juntos
            i = rng.integers(0, n, n)
            if np.linalg.matrix_rank(C[i]) < C.shape[1]: continue
            b = observables(*(x[i] for x in v), C[i]); B.append([b[k] for k in OBS])
        B = np.array(B); est = np.array([o[k] for k in OBS])
        Z = fz(B); zc = fz(est); se = Z.std(0, ddof=1)
        crit = np.percentile(np.abs((Z - zc) / se).max(1), 95)                # maximo estudentizado: intervalos simultaneos
        ci = {}
        for j, k in enumerate(OBS):
            ci[f"{k} lo"], ci[f"{k} hi"] = np.percentile(B[:, j], [2.5, 97.5])
            ci[f"{k} sim_lo"], ci[f"{k} sim_hi"] = np.tanh(zc[j] - crit * se[j]), np.tanh(zc[j] + crit * se[j])
        COV[label] = np.cov(Z.T)
        rows.append(dict(cohort=label, window=window, **o, **ci, simultaneous_critical_value=crit, n_bootstrap=len(B)))
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
    return pd.DataFrame(rows), COV


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
SIM_WINDOW = {"morning biopsy": "morning biopsy", "post mortem, any hour": "any hour"}
POINT = ["sd_F", "rho_K", "proxy_noise"]


def predictions(S, model, window):
    """La prediccion del modelo en cada punto de la rejilla: la media de sus semillas, en z de Fisher.
    Devuelve una fila por punto, con los observables ya de vuelta en la escala de rho, cuantas
    semillas se promediaron y cuantas personas suman."""
    q = S[(S.model == model) & (S.window == window)].copy()
    q[OBS] = fz(q[OBS].values)
    g = q.groupby(POINT, as_index=False)
    P = g[OBS].mean(); P[OBS] = np.tanh(P[OBS].values)
    P["n_seeds"] = g.size()["size"].values
    P["n_simulated"] = g["n"].sum()["n"].values if "n" in q.columns else P.n_seeds * N_DONORS
    return P


def state(lo, hi):
    return "positive" if lo > 0 else "negative" if hi < 0 else "zero"


def verdict(O, S):
    """Observable por observable y cohorte por cohorte: el estado de lo observado (positivo, negativo
    o compatible con cero, por su intervalo al 95 %) y la fraccion de las simulaciones de cada modelo
    que cae dentro de ese intervalo. Un observable aislado no decide nada; decide el conjunto (joint)."""
    rows = []
    for _, o in O.iterrows():
        for k in OBS:
            lo, hi = o[f"{k} lo"], o[f"{k} hi"]; st = state(lo, hi)
            for model in MODELS:
                P = predictions(S, model, SIM_WINDOW[o.window]); v = P[k].values
                t = 1.96 / np.sqrt(max(P.n_simulated.mean() - 3, 1))            # umbral con las personas que suma cada punto
                can = {"positive": bool((v > t).any()), "negative": bool((v < -t).any()), "zero": bool((np.abs(v) <= t).any())}
                inside = float(np.mean((v >= lo) & (v <= hi)))
                rows.append(dict(cohort=o.cohort, window=o.window, observable=k, observed=o[k], observed_ci_lo=lo, observed_ci_hi=hi,
                                 observed_state=st, model=model, n_grid_points=len(v), simulated_min=v.min(),
                                 simulated_median=float(np.median(v)), simulated_max=v.max(),
                                 model_can_be_positive=can["positive"], model_can_be_negative=can["negative"], model_can_be_zero=can["zero"],
                                 state_reproduced=can[st], fraction_of_grid_points_within_ci=inside,
                                 never_within_ci=inside == 0.0))
    return pd.DataFrame(rows)


def joint(O, COV, S):
    """El criterio conjunto, cohorte por cohorte: que fraccion de las simulaciones de cada modelo es
    compatible con los seis observables de la cohorte a la vez (ver la cabecera)."""
    thr = float(stats.chi2.ppf(0.95, len(OBS))); rows = []
    for _, o in O.iterrows():
        zc = fz([o[k] for k in OBS]); C = np.asarray(COV[o.cohort])
        for model in MODELS:
            P = predictions(S, model, SIM_WINDOW[o.window])
            D = fz(P[OBS].values) - zc
            # incertidumbre de lo observado mas la de la prediccion, que depende de cuantas personas suma cada punto
            d2 = np.array([d @ np.linalg.pinv(C * (1.0 + o.n / ns)) @ d for d, ns in zip(D, P.n_simulated.values)])
            okj = d2 <= thr
            oks = np.all([(P[k].values >= o[f"{k} sim_lo"]) & (P[k].values <= o[f"{k} sim_hi"]) for k in OBS], axis=0)
            best = P.iloc[int(np.argmin(d2))]
            rows.append(dict(cohort=o.cohort, window=o.window, n_observed=int(o.n), model=model, n_grid_points=len(P),
                             seeds_per_point=float(P.n_seeds.mean()), people_per_point=float(P.n_simulated.mean()),
                             n_compatible=int(okj.sum()), fraction_compatible=float(okj.mean()),
                             min_distance=float(d2.min()), median_distance=float(np.median(d2)), distance_threshold=thr,
                             n_all_within_simultaneous_ci=int(oks.sum()), fraction_all_within_simultaneous_ci=float(oks.mean()),
                             closest_sd_F=best.sd_F, closest_rho_K=best.rho_K, closest_proxy_noise=best.proxy_noise,
                             **{f"closest {k.split()[0]}": float(best[k]) for k in OBS}))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    pd.set_option("display.width", 220)
    O, COV = observed()
    if len(O):
        O.to_csv(f"{OUT}/new/competing_models_observed.tsv", sep="\t", index=False)
        pd.DataFrame([dict(cohort=c, observable_1=OBS[i], observable_2=OBS[j], covariance_fisher_z=M[i, j])
                      for c, M in COV.items() for i in range(len(OBS)) for j in range(len(OBS))]
                     ).to_csv(f"{OUT}/new/competing_models_observed_covariance.tsv", sep="\t", index=False)
        print("OBSERVADO, cohortes con dos tejidos por persona (IC 95 % por bootstrap de personas; entre parentesis, el simultaneo):")
        for _, r in O.iterrows():
            print(f"  {r.cohort} (n = {int(r.n)}, {r.window}):")
            print("     " + "; ".join(f"{k.split()[0]} {r[k]:+.2f} [{r[k + ' lo']:+.2f}, {r[k + ' hi']:+.2f}] ({r[k + ' sim_lo']:+.2f}, {r[k + ' sim_hi']:+.2f})" for k in OBS))
    fs = f"{OUT}/new/competing_models_simulated.tsv"
    if os.environ.get("T2D_MODEL_REUSE") == "1" and os.path.exists(fs):
        S = pd.read_csv(fs, sep="\t"); print(f"\n(se reutilizan las {len(S)} simulaciones de {fs})")
    else:
        print(f"\nSIMULANDO {len(MODELS)} modelos, {len(SEEDS)} semillas:")
        S = simulate(); S.to_csv(fs, sep="\t", index=False)
    print("\nSIMULADO, intervalo de cada observable al recorrer lo que no se conoce:")
    for win in ("morning biopsy", "any hour"):
        print(f"  [{win}]")
        for model in MODELS:
            q = S[(S.model == model) & (S.window == win)]
            print(f"    {model:24s} " + "  ".join(f"{k.split()[0]} {q[k].min():+.2f}..{q[k].max():+.2f}" for k in OBS))
    if len(O):
        V = verdict(O, S); V.to_csv(f"{OUT}/new/competing_models_verdict.tsv", sep="\t", index=False)
        J = joint(O, COV, S); J.to_csv(f"{OUT}/new/competing_models_joint.tsv", sep="\t", index=False)
        print(f"\nCRITERIO CONJUNTO, cohorte por cohorte: puntos de la rejilla cuya prediccion (media de {J.seeds_per_point.max():.0f} semillas, "
              f"{J.people_per_point.max():.0f} personas) es compatible con los seis observables a la vez (distancia <= {J.distance_threshold.iloc[0]:.1f})")
        for c in J.cohort.unique():
            print(f"  [{c}]")
            for r in J[J.cohort == c].itertuples():
                print(f"    {r.model:24s} {r.n_compatible:2d} de {r.n_grid_points:2d} ({r.fraction_compatible:.2f})   distancia minima {r.min_distance:8.1f}"
                      f"   | dentro de los seis intervalos simultaneos: {r.n_all_within_simultaneous_ci} de {r.n_grid_points}")
        print("\nPOR OBSERVABLE: los que ningun punto de la rejilla de un modelo predice dentro del intervalo observado")
        for c in V.cohort.unique():
            print(f"  [{c}]")
            for model in MODELS:
                bad = V[(V.cohort == c) & (V.model == model) & V.never_within_ci]
                print(f"    {model:24s} " + ("ninguno" if not len(bad) else
                      "; ".join(f"{r.observable.split()[0]} (observado {r.observed:+.2f}, {r.observed_state}; modelo {r.simulated_min:+.2f} a {r.simulated_max:+.2f})" for r in bad.itertuples())))
