#!/usr/bin/env python3
"""Cuatro mecanismos compitiendo contra los mismos datos.

El modelo del articulo lleva dentro lo que quiere explicar: el forzamiento efectivo es F*K, el
ritmo de cortisol por la disponibilidad del receptor. Para saber que aporta ese producto se pone a
competir con mecanismos mas simples, con el MISMO oscilador, el mismo ruido, el mismo muestreo y el
mismo calculo de los observables, cambiando solo de que depende el forzamiento de cada tejido:

    M0  sin forzamiento externo            0
    M1  solo el ritmo de cortisol          F_i                (el receptor varia pero no actua)
    M2  cortisol por receptor              F_i * K_it         (el del articulo)
    M3  solo el receptor                   Fmedio * K_it      (el mismo ritmo en todas las personas)

F_i es de la persona y lo comparten sus tejidos; K_it es de cada tejido. Cada persona aporta dos
tejidos y una hora de biopsia. De cada cohorte simulada se calculan los mismos seis observables que
de los datos, sobre la puntuacion de los genes de salida a la hora de la biopsia y sobre una medida
ruidosa del receptor, que es lo que da el ARNm:

    O1  receptor ~ salida, dentro del tejido
    O2  salida de un tejido ~ salida del otro, misma persona
    O3  receptor de un tejido ~ receptor del otro            (un dato que el modelo tambien tiene que cumplir)
    O4  receptor de un tejido ~ salida del OTRO
    O5  receptor ~ salida del mismo tejido, descontado el receptor del otro
    O6  salida ~ salida entre tejidos, descontados los dos receptores

LO QUE NO SE CONOCE SE RECORRE, no se ajusta: cuanto varia el cortisol entre personas, cuanto se
parece el receptor de dos tejidos, cuanto ruido tiene el ARNm como medida del receptor y, donde el
diseno no lo fija, A QUE HORAS se tomaron las muestras. Los cuatro modelos se evaluan sobre la
misma rejilla completa, tambien el que no tiene forzamiento, y ven a las mismas personas, a las
mismas horas y con el mismo ruido de medida.

GTEx NO ES UNA COHORTE. La circunstancia de la muerte (escala de Hardy) decide como y cuando se
obtuvo el tejido: en los donantes con ventilador la obtencion es una cirugia programada, y en la
muerte rapida ocurre a la hora que sea. Son dos regimenes de muestreo distintos y el modelo predice
cosas distintas para cada uno, asi que GTEx se puntua entero y tambien por estratos.

LA VENTANA DE MUESTREO. La puntuacion de salida es un nivel a una hora. Si todas las biopsias se
toman a la misma hora, quien tiene mas amplitud tiene mas nivel y el receptor se asocia con la
salida; si se toman a cualquier hora, ese nivel sube o baja segun la hora y la asociacion se
borra, mientras que dos tejidos de la misma persona, tomados a la vez, coinciden casi del todo. En
una cohorte de biopsias de manana la ventana es un dato del diseno (08:00, con 1.5 h de
dispersion). En una cohorte post mortem no lo es: las horas de obtencion no son publicas y no hay
razon para suponerlas repartidas por igual en las 24 h. Por eso alli la ventana es una incognita
mas de la rejilla (centro y dispersion), y el resultado dice que ventana necesitaria cada modelo.
Se da tambien, aparte, el resultado bajo el supuesto de horas uniformes.

EL CRITERIO CONJUNTO. Un modelo reproduce una cohorte si genera los seis observables a la vez:

  - Cohorte por cohorte, sin mezclar intervalos de cohortes distintas.
  - La prediccion de un modelo en un punto de la rejilla es la MEDIA de sus semillas; anadir
    semillas hace el criterio mas exigente. competing_models_seed_stability.tsv repite todo con 3,
    10, 30... semillas para ver si las conclusiones dependen de cuantas se usen.
  - Con una region conjunta: la distancia de la prediccion a lo observado se mide con la covarianza
    de los seis observables (bootstrap de personas, escala z de Fisher) mas la de la simulacion.
    Los seis no son seis dimensiones independientes (O1 y O5, O2 y O6, casi coinciden), asi que
    se informa el rango efectivo de esa covarianza y se usa ese rango, no seis.
  - CALIBRADO, no supuesto. Que esa distancia siga una ji cuadrado es una aproximacion. Se comprueba
    generando, del propio modelo, cohortes del tamano exacto de la observada y midiendo cuantas
    quedan dentro de la region nominal del 95 %. Esas cohortes salen de personas integradas aparte,
    con semillas que ninguna prediccion usa, y de varias poblaciones independientes, para ver cuanto
    cambia el umbral de una a otra (competing_models_calibration.tsv). Se usa el mayor.

LA CACHE. T2D_MODEL_REUSE=1 reutiliza las integraciones y la tabla simulada de una corrida anterior,
pero solo si su manifiesto (competing_models_cache_manifest.json: huella del codigo que las genera
y de todos sus parametros) coincide con el actual. Si no coincide se vuelven a generar.

POR OBSERVABLE. competing_models_verdict.tsv da, por cohorte y observable, el estado de lo observado
(positivo, negativo o compatible con cero, por su intervalo) y la fraccion de puntos de cada modelo
cuya prediccion tiene ese mismo estado Y cae dentro del intervalo, las dos cosas en el mismo punto.
El estado de una prediccion se decide con su error estandar entre semillas.

LO QUE ESTO NO HACE. No ajusta ningun parametro a los datos y no elige un modelo por verosimilitud.
Dice que mecanismos no pueden generar lo observado. Dos que si pueden quedan sin distinguir, y eso
se informa asi.

Entradas (si existen): results/new/AAGMEx_*_log2.pkl, GTEx_*.pkl, GTEx_pheno.pkl y la matriz de serie
                       de AAGMEx. Sin ellas se escriben solo las simulaciones.
Salidas (results/new/): competing_models_observed.tsv, _observed_covariance.tsv, _simulated.tsv,
          _joint.tsv, _verdict.tsv, _calibration.tsv, _seed_stability.tsv
"""
import os, sys, gzip, pickle, warnings; warnings.filterwarnings("ignore")
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
SEEDS = list(range(int(os.environ.get("T2D_MODEL_SEEDS", 10))))
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


HARDY = {"ventilator (Hardy 0)": [0], "fast death (Hardy 1-2)": [1, 2], "intermediate or slow death (Hardy 3-4)": [3, 4]}


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
                # GTEx no es una cohorte sino varias: como murio el donante decide como y cuando se obtuvo el
                # tejido. Con ventilador (Hardy 0) la obtencion es una cirugia programada; en la muerte rapida
                # (1 y 2) ocurre a la hora que sea. Cada estrato se puntua por separado.
                for lab, ks in HARDY.items():
                    m = p.DTHHRDY.isin(ks).values
                    if m.sum() < 60: continue
                    sub_ids = [i for i, k in zip(ids, m) if k]; q = p.loc[sub_ids]
                    Cs = [np.ones(len(sub_ids)), (q.SEX == 1).astype(float), q.AGE.map(agemid).astype(float)]
                    Cs += [(q.DTHHRDY == k).astype(float) for k in ks[1:]]
                    add(f"GTEx {a}-{b}, {lab}", "post mortem, any hour", G[a][sub_ids], G[b][sub_ids], np.column_stack(Cs))
    return pd.DataFrame(rows), COV


# ------------------------------------------------------------------ las simulaciones
SD_F, RHO_K, PROXY = sorted({g["sd_F"] for g in GRID}), sorted({g["rho_K"] for g in GRID}), sorted({g["proxy_noise"] for g in GRID})
MORNING = "centre 08 h, sd 1.5 h"                                            # el diseno de una cohorte de biopsias de manana
WINDOWS = {f"centre {c:02d} h, sd {sd} h": (float(c), float(sd)) for sd in (1.5, 3.0, 4.5, 6.0) for c in (2, 8, 14, 20)}
WINDOWS["uniform over 24 h"] = (None, None)
# que ventanas admite cada tipo de cohorte: la de manana esta fijada por el diseno; la post mortem no se conoce
ALLOWED = {"morning biopsy": [MORNING], "post mortem, any hour": list(WINDOWS)}
POINT = ["sd_F", "rho_K", "proxy_noise", "window"]


def people(sd_F, rho_K, seed, n=N_DONORS, sd_K=0.35):
    """Las personas de un escenario: su cortisol y su receptor en cada tejido. Las mismas para todos los
    modelos. El cortisol depende solo de (semilla, sd_F) y el receptor solo de (semilla, rho_K)."""
    rf = np.random.default_rng([seed, 1, int(sd_F * 100)]); rk = np.random.default_rng([seed, 2, int(rho_K * 100)])
    F = np.clip(F_MEAN * (1 + sd_F * rf.normal(0, 1, n)), 0.0, None)
    shared = rk.normal(0, 1, n)                                              # parte del receptor comun a los dos tejidos
    K = {t: np.clip(K_MEAN + sd_K * (np.sqrt(rho_K) * shared + np.sqrt(1 - rho_K) * rk.normal(0, 1, n)), 0.05, None) for t in "AM"}
    return F, K


def ensemble_of(task):
    """Amplitud y fase del ensamble de cada persona y tejido bajo un modelo. No depende de la hora de
    biopsia ni del ruido de medida, que se aplican despues. La realizacion celular (periodos y fases
    iniciales) es la misma para los cuatro modelos."""
    model, sd_F, rho_K, seed = task
    F, K = people(sd_F, rho_K, seed)
    out = {}
    for i, t in enumerate("AM"):
        a, p, _ = fc.simulate_cohort(MODELS[model](F, K[t], F_MEAN, K_MEAN), n_cells=40, days=12, seed=[seed, 3, i])
        out[t] = (a, p)
    return task, out


def simulate_ensembles(seeds):
    """Todas las cohortes, sobre la rejilla completa para los cuatro modelos. Donde un parametro no
    entra en el forzamiento de un modelo la integracion es la misma y no se repite (sin forzamiento
    no entra ninguno; con solo cortisol no entra el receptor; con solo receptor no entra cuanto varia
    el cortisol), pero el punto de la rejilla se evalua igual, porque el receptor medido si cambia."""
    from concurrent.futures import ProcessPoolExecutor
    key = {"M0 no forcing": lambda sf, rk: (SD_F[0], RHO_K[0]), "M1 cortisol only": lambda sf, rk: (sf, RHO_K[0]),
           "M2 cortisol x receptor": lambda sf, rk: (sf, rk), "M3 receptor only": lambda sf, rk: (SD_F[0], rk)}
    unique = sorted({(m, *key[m](sf, rk), sd) for m in MODELS for sf in SD_F for rk in RHO_K for sd in seeds})
    workers = max(1, min(int(os.environ.get("T2D_MODEL_WORKERS", os.cpu_count() or 1)), 32))
    print(f"  {len(unique)} integraciones de {N_DONORS} personas y dos tejidos, {len(seeds)} semillas, {workers} procesos", flush=True)
    with ProcessPoolExecutor(workers) as ex:
        done = dict(ex.map(ensemble_of, unique, chunksize=max(1, len(unique) // (workers * 4))))
    return {(m, sf, rk, sd): done[(m, *key[m](sf, rk), sd)] for m in MODELS for sf in SD_F for rk in RHO_K for sd in seeds}


def sampling(seed, window, n=N_DONORS):
    """Horas de biopsia y ruido de medida de una semilla y una ventana: los mismos para todos los modelos
    y todos los puntos de la rejilla."""
    rng = np.random.default_rng([seed, 4, list(WINDOWS).index(window)])
    c, sd = WINDOWS[window]
    hour = rng.uniform(0, 24, n) if c is None else rng.normal(c, sd, n) % 24
    return hour, {t: rng.normal(0, 1, n) for t in "AM"}, {t: int(rng.integers(1e6)) for t in "AM"}


def measured(ENS, model, sd_F, rho_K, seed, window, proxy_noise, expr_noise=0.6):
    """Lo que se mide en cada persona: la puntuacion de salida de cada tejido a la hora de la biopsia y
    un sustituto ruidoso del receptor. Devuelve (receptor A, receptor M, salida A, salida M)."""
    _, K = people(sd_F, rho_K, seed); ens = ENS[(model, sd_F, rho_K, seed)]
    hour, eps, eseed = sampling(seed, window)
    out = {}
    for t in "AM":
        d = pd.DataFrame(dict(subject=np.arange(len(hour)), hour=hour, amp_ens=ens[t][0], phase_ens=ens[t][1]))
        E = fc.synthetic_expression(d, noise=expr_noise, seed=eseed[t])
        Z = E.loc[OUTPUT_GENES]; out[t] = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.mean(0).values
    rec = {t: K[t] + proxy_noise * np.std(K[t]) * eps[t] for t in "AM"}
    return rec["A"], rec["M"], out["A"], out["M"]


def table(ENS, seeds):
    rows = []
    for (model, sd_F, rho_K, seed) in ENS:
        if seed not in seeds: continue
        for w in WINDOWS:
            for pn in PROXY:
                rows.append(dict(model=model, sd_F=sd_F, rho_K=rho_K, proxy_noise=pn, window=w, seed=seed,
                                 **observables(*measured(ENS, model, sd_F, rho_K, seed, w, pn))))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ el criterio
def predictions(S, model, windows, nseeds=None):
    """La prediccion del modelo en cada punto de la rejilla y cada ventana admitida: la media de sus
    semillas en z de Fisher, con su error estandar entre semillas."""
    q = S[(S.model == model) & S.window.isin(windows)]
    if nseeds is not None: q = q[q.seed < nseeds]
    q = q.copy(); q[OBS] = fz(q[OBS].values)
    g = q.groupby(POINT, as_index=False, sort=True)
    P = g[OBS].mean(); k = g.size()["size"].values
    SE = g[OBS].std(ddof=1)[OBS].values / np.sqrt(k)[:, None] if (k > 1).all() else np.full((len(P), len(OBS)), np.nan)
    Z = P[OBS].values.copy(); P[OBS] = np.tanh(Z)
    P["n_seeds"] = k; P["n_simulated"] = g["n"].sum()["n"].values
    return P, Z, SE


def region(C, rtol=1e-3):
    """Inversa de la covarianza restringida a sus direcciones informativas, cuantas son, y el numero
    de condicion entre ellas. Dos observables que casi coinciden dan un autovalor casi nulo:
    invertirlo tal cual convertiria ruido numerico en distancias enormes, y contar seis grados de
    libertad seria contar de mas.

    Se niega a seguir, con un mensaje que dice por que, si la matriz no sirve para definir una
    region: valores no finitos, forma que no es cuadrada, ninguna direccion con varianza, o un
    autovalor claramente negativo (no es una covarianza). Un autovalor negativo del tamano del error
    de redondeo se trata como cero."""
    C = np.asarray(C, float)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise ValueError(f"la covarianza debe ser cuadrada y tiene forma {C.shape}")
    if not np.isfinite(C).all():
        raise ValueError("la covarianza tiene valores no finitos: algun observable no se pudo calcular en el bootstrap")
    C = (C + C.T) / 2
    w, V = np.linalg.eigh(C)
    top = w.max()
    if not top > 0:
        raise ValueError("la covarianza no tiene ninguna direccion con varianza: no define una region")
    if w.min() < -1e-8 * top:
        raise ValueError(f"la matriz no es una covarianza: autovalor {w.min():.3g} frente a un maximo de {top:.3g}")
    keep = w > top * rtol
    return (V[:, keep] / w[keep]) @ V[:, keep].T, int(keep.sum()), float(top / w[keep].min())


def distance(Z, zc, C, n_obs, n_sim):
    """Distancia de cada prediccion a lo observado, con la incertidumbre de lo observado mas la de la
    prediccion (la misma covarianza, escalada por las personas que suma cada punto)."""
    Ci, rank, cond = region(C)
    D = Z - zc
    return np.array([d @ (Ci / (1.0 + n_obs / ns)) @ d for d, ns in zip(D, n_sim)]), rank, cond


def state(lo, hi):
    return "positive" if lo > 0 else "negative" if hi < 0 else "zero"


def joint(O, COV, S, threshold=None, nseeds=None):
    """El criterio conjunto, cohorte por cohorte. `threshold` es un diccionario cohorte -> umbral
    calibrado; sin el se usa el cuantil de la ji cuadrado con el rango efectivo de la covarianza."""
    rows = []
    for _, o in O.iterrows():
        zc = fz([o[k] for k in OBS]); C = COV[o.cohort]
        for model in MODELS:
            for scope, wins in (("windows allowed by the design", ALLOWED[o.window]),) + \
                               ((("uniform hours only", ["uniform over 24 h"]),) if len(ALLOWED[o.window]) > 1 else ()):
                P, Z, _ = predictions(S, model, wins, nseeds)
                if not len(P): continue                                      # esa ventana no se simulo
                d2, rank, cond = distance(Z, zc, C, o.n, P.n_simulated.values)
                thr = (threshold or {}).get(o.cohort, float(stats.chi2.ppf(0.95, rank)))
                okj = d2 <= thr
                oks = np.all([(P[k].values >= o[f"{k} sim_lo"]) & (P[k].values <= o[f"{k} sim_hi"]) for k in OBS], axis=0)
                best = P.iloc[int(np.argmin(d2))]
                rows.append(dict(cohort=o.cohort, n_observed=int(o.n), scope=scope, model=model, n_points=len(P),
                                 seeds_per_point=float(P.n_seeds.mean()), people_per_point=float(P.n_simulated.mean()),
                                 n_compatible=int(okj.sum()), fraction_compatible=float(okj.mean()),
                                 min_distance=float(d2.min()), median_distance=float(np.median(d2)), distance_threshold=thr,
                                 threshold_is_calibrated=bool(threshold and o.cohort in threshold),
                                 covariance_effective_rank=rank, covariance_condition_number=cond,
                                 n_all_within_simultaneous_ci=int(oks.sum()),
                                 closest_window=best.window, closest_sd_F=best.sd_F, closest_rho_K=best.rho_K,
                                 closest_proxy_noise=best.proxy_noise, **{f"closest {k.split()[0]}": float(best[k]) for k in OBS}))
    return pd.DataFrame(rows)


def verdict(O, S, nseeds=None):
    """Observable por observable y cohorte por cohorte. El estado de una prediccion se decide con su
    error estandar entre semillas, y un punto cuenta solo si su prediccion tiene el mismo estado que lo
    observado Y cae dentro del intervalo observado: las dos cosas en el mismo punto."""
    rows = []
    for _, o in O.iterrows():
        for model in MODELS:
            P, Z, SE = predictions(S, model, ALLOWED[o.window], nseeds)
            for j, k in enumerate(OBS):
                lo, hi = o[f"{k} lo"], o[f"{k} hi"]; st = state(lo, hi); v = P[k].values
                pst = np.array([state(np.tanh(z - 1.96 * e), np.tanh(z + 1.96 * e)) if np.isfinite(e) else "undetermined"
                                for z, e in zip(Z[:, j], SE[:, j])])
                inside = (v >= lo) & (v <= hi)
                rows.append(dict(cohort=o.cohort, observable=k, observed=o[k], observed_ci_lo=lo, observed_ci_hi=hi,
                                 observed_state=st, model=model, n_points=len(v), predicted_min=v.min(),
                                 predicted_median=float(np.median(v)), predicted_max=v.max(),
                                 fraction_predicted_positive=float(np.mean(pst == "positive")),
                                 fraction_predicted_negative=float(np.mean(pst == "negative")),
                                 fraction_predicted_zero=float(np.mean(pst == "zero")),
                                 fraction_within_ci=float(inside.mean()),
                                 fraction_same_state_and_within_ci=float(np.mean(inside & (pst == st))),
                                 never_reproduced=not bool((inside & (pst == st)).any())))
    return pd.DataFrame(rows)


CALIBRATION_SEED = 10_000        # las semillas de calibracion no coinciden con ninguna de las que definen las predicciones
_POOLS = {}


def independent_population(model, point, n, pool):
    """n personas nuevas del modelo en un punto de la rejilla: su receptor medido y su salida en los dos
    tejidos. Se integran aparte, con semillas que no usa ninguna prediccion, de modo que lo que sirve
    para calibrar no ha intervenido en lo que se calibra. Se guarda por si otra cohorte pide lo mismo."""
    sd_F, rho_K, pn, w = point
    key = (model, sd_F, rho_K, w, pn, pool)
    if key in _POOLS and len(_POOLS[key][0]) >= n:
        return [x[:n] for x in _POOLS[key]]
    seed = CALIBRATION_SEED + pool
    F, K = people(sd_F, rho_K, seed, n=n)
    hour, eps, eseed = sampling(seed, w, n=n)
    out = {}
    for i, t in enumerate("AM"):
        a, p, _ = fc.simulate_cohort(MODELS[model](F, K[t], F_MEAN, K_MEAN), n_cells=40, days=12, seed=[seed, 3, i])
        d = pd.DataFrame(dict(subject=np.arange(n), hour=hour, amp_ens=a, phase_ens=p))
        E = fc.synthetic_expression(d, noise=0.6, seed=eseed[t])
        Z = E.loc[OUTPUT_GENES]; out[t] = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.mean(0).values
    rec = {t: K[t] + pn * np.std(K[t]) * eps[t] for t in "AM"}
    _POOLS[key] = (rec["A"], rec["M"], out["A"], out["M"])
    return list(_POOLS[key])


def calibrate(o, ENS, model, point, seeds, reps=100, boot=150, pools=2, pool_factor=3):
    """Se cumple que la region nominal del 95 % cubre el 95 %? Se toma un punto de un modelo como si fuera
    la verdad y se le aplica a cohortes suyas el mismo procedimiento que a los datos: sus seis
    observables, su covarianza por bootstrap y su distancia a la prediccion del modelo.

    Dos condiciones para que la respuesta valga:
      - Las cohortes de calibracion son INDEPENDIENTES de la prediccion. La prediccion es la media de
        las semillas de la rejilla, igual que en joint(); las cohortes salen de personas integradas
        aparte, con otras semillas.
      - Tienen EXACTAMENTE el tamano de la cohorte observada. La poblacion de la que se extraen tiene
        `pool_factor` veces ese tamano, sea cual sea.
    Se repite con `pools` poblaciones independientes, y se informa de cuanto cambian la cobertura y el
    umbral de una a otra. Las replicas de una misma poblacion comparten personas; las de poblaciones
    distintas, ninguna."""
    def vec(*v):
        ob = observables(*v); return [ob[k] for k in OBS]
    sd_F, rho_K, pn, w = point
    zp = np.array([fz(vec(*measured(ENS, model, sd_F, rho_K, s, w, pn))) for s in seeds]).mean(0)
    n = int(o.n); n_sim = len(seeds) * N_DONORS
    cov, thr, med, ranks = [], [], [], []
    for pool in range(pools):
        pop = independent_population(model, point, pool_factor * n, pool)
        assert len(pop[0]) >= pool_factor * n
        rng = np.random.default_rng([CALIBRATION_SEED, pool]); d2 = []
        for _ in range(reps):
            i = rng.choice(len(pop[0]), n, replace=False); v = [x[i] for x in pop]
            assert len(v[0]) == n
            B = np.array([vec(*(x[j] for x in v)) for j in (rng.integers(0, n, n) for _ in range(boot))])
            d, rank, _ = distance(zp[None, :], fz(vec(*v)), np.cov(fz(B).T), n, [n_sim]); d2.append(d[0]); ranks.append(rank)
        d2 = np.array(d2); nominal = float(stats.chi2.ppf(0.95, int(np.median(ranks))))
        cov.append(float(np.mean(d2 <= nominal))); thr.append(float(np.percentile(d2, 95))); med.append(float(np.median(d2)))
    rank = int(np.median(ranks)); nominal = float(stats.chi2.ppf(0.95, rank))
    return dict(cohort=o.cohort, n_observed=int(o.n), n_calibration_cohort=n, population_size=pool_factor * n,
                independent_populations=pools, replicates_per_population=reps, bootstrap_per_replicate=boot,
                truth_model=model, truth_sd_F=sd_F, truth_rho_K=rho_K, truth_proxy_noise=pn, truth_window=w,
                prediction_seeds=len(seeds), effective_rank=rank, nominal_threshold=nominal,
                coverage_of_nominal_threshold=float(np.mean(cov)), coverage_min=min(cov), coverage_max=max(cov),
                calibrated_threshold=float(np.mean(thr)), calibrated_threshold_min=min(thr), calibrated_threshold_max=max(thr),
                median_distance_under_truth=float(np.mean(med)))


# ------------------------------------------------------------------ la cache y su manifiesto
def manifest(seeds):
    """Que codigo y que parametros generan las integraciones y la tabla simulada. Una cache solo se
    reutiliza si su manifiesto es identico a este: cambiar una ecuacion, el numero de personas, la
    rejilla o las semillas la invalida sin que nadie tenga que acordarse de borrarla."""
    import hashlib, inspect
    def h(*objs):
        m = hashlib.sha256()
        for x in objs: m.update((x if isinstance(x, str) else inspect.getsource(x)).encode())
        return m.hexdigest()
    params = dict(n_donors=N_DONORS, seeds=list(seeds), f_mean=F_MEAN, k_mean=K_MEAN, sd_F=list(SD_F), rho_K=list(RHO_K),
                  models=list(MODELS), n_cells=40, days=12)
    gen = h(fc.simulate_cohort, people, ensemble_of, simulate_ensembles, *[inspect.getsource(v) for v in MODELS.values()], repr(sorted(params.items())))
    tab = h(gen, sampling, measured, table, observables, resid, fc.synthetic_expression,
            repr(sorted(WINDOWS.items(), key=str)), repr(PROXY), repr(OBS), repr(OUTPUT_GENES))
    return dict(parameters=params, integrations_sha256=gen, simulated_table_sha256=tab)


def cached(path, key, current):
    """True si hay cache, se pidio reutilizarla y su manifiesto coincide con el codigo y los parametros
    actuales en la parte `key`. Si no coincide lo dice y devuelve False."""
    import json
    if os.environ.get("T2D_MODEL_REUSE") != "1" or not os.path.exists(path): return False
    mf = f"{OUT}/new/competing_models_cache_manifest.json"
    old = json.load(open(mf)) if os.path.exists(mf) else {}
    if old.get(key) == current[key]: return True
    print(f"  la cache {os.path.basename(path)} no corresponde al codigo o a los parametros actuales: se vuelve a generar", flush=True)
    return False


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    pd.set_option("display.width", 220)
    O, COV = observed()
    if len(O):
        O.to_csv(f"{OUT}/new/competing_models_observed.tsv", sep="\t", index=False)
        pd.DataFrame([dict(cohort=c, observable_1=OBS[i], observable_2=OBS[j], covariance_fisher_z=M[i, j])
                      for c, M in COV.items() for i in range(len(OBS)) for j in range(len(OBS))]
                     ).to_csv(f"{OUT}/new/competing_models_observed_covariance.tsv", sep="\t", index=False)
        print("OBSERVADO, cohortes con dos tejidos por persona (IC 95 % por bootstrap de personas):")
        for _, r in O.iterrows():
            _, rank, cond = region(COV[r.cohort])
            print(f"  {r.cohort} (n = {int(r.n)}, {r.window}); rango efectivo de la covarianza {rank} de {len(OBS)}")
            print("     " + "; ".join(f"{k.split()[0]} {r[k]:+.2f} [{r[k + ' lo']:+.2f}, {r[k + ' hi']:+.2f}]" for k in OBS))
    import json
    fe = f"{OUT}/new/competing_models_ensembles.pkl"; seeds = SEEDS; MAN = manifest(seeds)
    if cached(fe, "integrations_sha256", MAN):
        ENS = pickle.load(open(fe, "rb"))
        print(f"\n(se reutilizan las integraciones de {fe}: {len(seeds)} semillas; su manifiesto coincide)")
    else:
        print(f"\nSIMULANDO {len(MODELS)} modelos sobre la rejilla completa:")
        ENS = simulate_ensembles(seeds); pickle.dump(ENS, open(fe, "wb"))
    fs = f"{OUT}/new/competing_models_simulated.tsv"
    if cached(fs, "simulated_table_sha256", MAN):
        S = pd.read_csv(fs, sep="\t")
    else:
        S = table(ENS, seeds); S.to_csv(fs, sep="\t", index=False)
    json.dump(MAN, open(f"{OUT}/new/competing_models_cache_manifest.json", "w"), indent=1)
    print(f"\nPREDICHO con la ventana de manana ({MORNING}) y con horas uniformes, intervalo sobre la rejilla (media de {len(seeds)} semillas):")
    for w in (MORNING, "uniform over 24 h"):
        print(f"  [{w}]")
        for model in MODELS:
            P, _, _ = predictions(S, model, [w])
            print(f"    {model:24s} " + "  ".join(f"{k.split()[0]} {P[k].min():+.2f}..{P[k].max():+.2f}" for k in OBS))
    if len(O):
        # 1. criterio con el umbral nominal, para saber que punto de que modelo se toma como verdad al calibrar
        J0 = joint(O, COV, S)
        cal = []; thr = {}
        print("\nCALIBRANDO la region conjunta, cohorte por cohorte:", flush=True)
        for _, o in O.iterrows():
            q = J0[(J0.cohort == o.cohort) & (J0.scope == "windows allowed by the design")].sort_values("min_distance").iloc[0]
            c = calibrate(o, ENS, q.model, (q.closest_sd_F, q.closest_rho_K, q.closest_proxy_noise, q.closest_window), seeds,
                          reps=int(os.environ.get("T2D_MODEL_CALIB", 100)), boot=int(os.environ.get("T2D_MODEL_CALIB_BOOT", 150)),
                          pools=int(os.environ.get("T2D_MODEL_CALIB_POOLS", 2)))
            # el umbral que se usa: el mayor de los calibrados y el nominal, para no excluir un modelo por un umbral corto
            cal.append(c); thr[o.cohort] = max(c["calibrated_threshold_max"], c["nominal_threshold"])
            print(f"  calibrada {o.cohort} (n = {c['n_calibration_cohort']}): cobertura del umbral nominal {c['coverage_min'] * 100:.0f} a "
                  f"{c['coverage_max'] * 100:.0f} % entre poblaciones", flush=True)
        CAL = pd.DataFrame(cal); CAL.to_csv(f"{OUT}/new/competing_models_calibration.tsv", sep="\t", index=False)
        print("\nCALIBRACION de la region conjunta: cohortes generadas del propio modelo, del tamano de la observada")
        for r in CAL.itertuples():
            print(f"  {r.cohort}: cohortes de {r.n_calibration_cohort} personas; rango efectivo {r.effective_rank}; el umbral nominal ({r.nominal_threshold:.1f}) "
                  f"cubre {r.coverage_of_nominal_threshold * 100:.0f} % ({r.coverage_min * 100:.0f} a {r.coverage_max * 100:.0f}); "
                  f"percentil 95 de las distancias {r.calibrated_threshold_min:.1f} a {r.calibrated_threshold_max:.1f}; se usa {thr[r.cohort]:.1f}")
        # 2. criterio con el umbral calibrado
        J = joint(O, COV, S, thr); J.to_csv(f"{OUT}/new/competing_models_joint.tsv", sep="\t", index=False)
        V = verdict(O, S); V.to_csv(f"{OUT}/new/competing_models_verdict.tsv", sep="\t", index=False)
        print("\nCRITERIO CONJUNTO, cohorte por cohorte: puntos cuya prediccion es compatible con los seis observables a la vez")
        for c in J.cohort.unique():
            for scope in J[J.cohort == c].scope.unique():
                print(f"  [{c}; {scope}]")
                for r in J[(J.cohort == c) & (J.scope == scope)].itertuples():
                    print(f"    {r.model:24s} {r.n_compatible:3d} de {r.n_points:3d} ({r.fraction_compatible:.2f})  distancia minima {r.min_distance:7.1f} (umbral {r.distance_threshold:.1f})"
                          + (f"  | mas cerca: {r.closest_window}" if scope.startswith("windows") and r.n_points > 18 else ""))
        # 3. depende de cuantas semillas se usen?
        st = []
        for k in sorted({x for x in (3, 10, 30, 50, len(seeds)) if x <= len(seeds)}):
            Jk = joint(O, COV, S, thr, nseeds=k); Jk = Jk[Jk.scope == "windows allowed by the design"]
            st.append(Jk.assign(seeds_used=k)[["seeds_used", "cohort", "model", "n_points", "n_compatible", "fraction_compatible", "min_distance", "distance_threshold"]])
        ST = pd.concat(st); ST.to_csv(f"{OUT}/new/competing_models_seed_stability.tsv", sep="\t", index=False)
        print("\nESTABILIDAD con el numero de semillas (compatibles; distancia minima), primera cohorte:")
        c0 = O.cohort.iloc[0]
        for k in ST.seeds_used.unique():
            print(f"  {k:3d} semillas: " + "; ".join(f"{r.model.split()[0]} {r.n_compatible}/{r.n_points} ({r.min_distance:.1f})"
                                                     for r in ST[(ST.seeds_used == k) & (ST.cohort == c0)].itertuples()))
