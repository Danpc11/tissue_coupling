#!/usr/bin/env python3
"""El criterio con el que se comparan los modelos competidores, comprobado a mano.

Sin simulaciones ni datos: se construyen lo observado y lo predicho.

  1. Tres estados por observable, y las dos condiciones en el MISMO punto de la rejilla: que la
     prediccion tenga el estado de lo observado y que caiga dentro de su intervalo.
  2. El criterio conjunto va cohorte por cohorte: una prediccion que acierta un observable de una
     cohorte y otro de otra no es compatible con ninguna.
  3. Es una region conjunta que respeta la correlacion entre observables.
  4. La concordancia del receptor entre tejidos cuenta.
  5. Se cuenta cuantos puntos, no si alguno.
  6. La prediccion de un punto es la media de sus semillas: mas semillas exigen mas.
  7. Una covarianza con direcciones redundantes no se invierte a ciegas: se usa su rango efectivo.
  8. Donde el diseno no fija la ventana de muestreo, se recorren todas; donde la fija, solo esa.
  9. Una covarianza que no define una region (no finita, nula, no definida positiva) se rechaza.
 10. La calibracion usa personas independientes de las que definen la prediccion, del tamano exacto
     de la cohorte observada, y de mas de una poblacion.
 11. La cache de simulaciones solo se reutiliza si su manifiesto coincide con el codigo y los
     parametros actuales.
"""
import os, sys, importlib.util, tempfile, json
os.environ.setdefault("T2D_MODEL_N", "40")                  # cohortes pequenas: aqui se comprueba el criterio, no el resultado
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("cm", os.path.join(ROOT, "python/model/competing_models.py"))
cm = importlib.util.module_from_spec(spec); spec.loader.exec_module(cm)
K = cm.OBS; names = list(cm.MODELS)

fails = []
def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: fails.append(name)

def cohort_row(label, values, half=0.08, window="morning biopsy", n=200):
    r = dict(cohort=label, window=window, n=n)
    for k, v in zip(K, values):
        r[k] = v; r[f"{k} lo"], r[f"{k} hi"] = v - half, v + half; r[f"{k} sim_lo"], r[f"{k} sim_hi"] = v - 1.4 * half, v + 1.4 * half
    return r

def preds(vectors_by_model, window=None, n=10**6, jitter=1e-4):
    """Cada vector es la prediccion de un punto distinto de la rejilla; dos semillas casi iguales y
    muchas personas, para que el ruido de la simulacion no cuente."""
    rows = []
    for model, vectors in vectors_by_model.items():
        for i, v in enumerate(vectors):
            for sd in (0, 1):
                rows.append(dict(model=model, sd_F=0.25, rho_K=0.3, proxy_noise=float(i), window=window or cm.MORNING, seed=sd, n=n,
                                 **{k: x + (jitter if sd else -jitter) for k, x in zip(K, v)}))
    return pd.DataFrame(rows)

se = 0.08 / 1.96; DIAG = np.eye(len(K)) * se ** 2
a = [0.20, 0.40, 0.25, 0.10, 0.20, 0.35]; b = [0.50, 0.10, 0.25, 0.30, 0.45, 0.05]

print("1. Tres estados, y las dos condiciones en el mismo punto")
O = pd.DataFrame([cohort_row("A", [-0.4] * 6)])
S = preds({names[0]: [[0.0] * 6], names[1]: [[0.4] * 6], names[2]: [[-0.42] * 6], names[3]: [[0.0] * 6, [-0.9] * 6]})
V = cm.verdict(O, S).set_index(["model", "observable"])
check("lo negativo se clasifica como negativo", V.loc[(names[0], K[0]), "observed_state"] == "negative")
check("un modelo que solo predice cero no lo reproduce", bool(V.loc[(names[0], K[0]), "never_reproduced"]))
check("uno que solo predice positivo tampoco", bool(V.loc[(names[1], K[0]), "never_reproduced"]))
check("uno que predice negativo dentro del intervalo si", V.loc[(names[2], K[0]), "fraction_same_state_and_within_ci"] == 1.0)
check("el estado en un punto y el intervalo en ninguno no suman",
      bool(V.loc[(names[3], K[0]), "never_reproduced"]) and V.loc[(names[3], K[0]), "fraction_predicted_negative"] == 0.5,
      "un punto predice el estado y ninguno cae en el intervalo")

print("2. Cohorte por cohorte, sin mezclar")
O = pd.DataFrame([cohort_row("A", a), cohort_row("B", b)]); COV = {"A": DIAG, "B": DIAG}
mix = [a[0], b[1], a[2], b[3], a[4], b[5]]
J = cm.joint(O, COV, preds({names[0]: [a], names[1]: [b], names[2]: [mix], names[3]: [[0.0] * 6]})).set_index(["cohort", "model"])
check("la prediccion igual a la cohorte A es compatible con A y no con B",
      J.loc[("A", names[0]), "n_compatible"] == 1 and J.loc[("B", names[0]), "n_compatible"] == 0)
check("la mezcla de las dos no es compatible con ninguna", J.loc[("A", names[2]), "n_compatible"] == 0 and J.loc[("B", names[2]), "n_compatible"] == 0)

print("3. Region conjunta, no seis intervalos")
C = DIAG.copy(); C[0, 4] = C[4, 0] = 0.9 * se ** 2
O = pd.DataFrame([cohort_row("A", a)]); COV = {"A": C}
along = list(a); along[0] += 0.05; along[4] += 0.05
against = list(a); against[0] += 0.05; against[4] -= 0.05
J = cm.joint(O, COV, preds({names[0]: [along], names[1]: [against], names[2]: [a], names[3]: [[0.0] * 6]})).set_index("model")
check("moverse a favor de la correlacion es compatible", J.loc[names[0], "n_compatible"] == 1, f"distancia {J.loc[names[0], 'min_distance']:.1f}")
check("moverse en contra no, aunque cada valor este dentro de su intervalo",
      J.loc[names[1], "n_compatible"] == 0 and J.loc[names[1], "n_all_within_simultaneous_ci"] == 1, f"distancia {J.loc[names[1], 'min_distance']:.1f}")

print("4. La concordancia del receptor entre tejidos cuenta")
O = pd.DataFrame([cohort_row("A", a)]); COV = {"A": DIAG}
no_o3 = list(a); no_o3[2] = 0.0
J = cm.joint(O, COV, preds({names[0]: [no_o3], names[1]: [a], names[2]: [a], names[3]: [a]})).set_index("model")
check("una prediccion sin esa concordancia no es compatible", J.loc[names[0], "n_compatible"] == 0 and J.loc[names[1], "n_compatible"] == 1)

print("5. Cuantos puntos, no si alguno")
far = [0.9] * 6
J = cm.joint(O, COV, preds({names[0]: [a] * 19 + [far], names[1]: [far] * 19 + [a], names[2]: [far] * 20, names[3]: [a] * 20})).set_index("model")
check("diecinueve de veinte, uno de veinte y ninguno",
      J.loc[names[0], "fraction_compatible"] == 0.95 and J.loc[names[1], "fraction_compatible"] == 0.05 and J.loc[names[2], "n_compatible"] == 0)

print("6. Mas semillas exigen mas")
rng = np.random.default_rng(0); truth = np.array([0.05, 0.40, 0.25, 0.10, 0.05, 0.35])
def seeds(k):
    rows = [dict(model=names[0], sd_F=0.25, rho_K=0.3, proxy_noise=1.0, window=cm.MORNING, seed=sd, n=200,
                 **dict(zip(K, np.tanh(np.arctanh(truth) + rng.normal(0, 1 / np.sqrt(197), 6))))) for sd in range(k)]
    return pd.DataFrame(rows + [dict(model=m, sd_F=0.25, rho_K=0.3, proxy_noise=1.0, window=cm.MORNING, seed=sd, n=200, **dict(zip(K, [0.0] * 6)))
                                for m in names[1:] for sd in (0, 1)])
O = pd.DataFrame([cohort_row("A", a, half=0.06, n=1000)]); COV = {"A": np.eye(6) * (0.06 / 1.96) ** 2}
d2, d60 = (cm.joint(O, COV, seeds(k)).set_index("model").loc[names[0]] for k in (2, 60))
check("un punto, sean cuantas sean las semillas", d2.n_points == 1 and d60.n_points == 1 and d60.seeds_per_point == 60)
check("con muchas semillas el modelo equivocado queda descartado", d60.n_compatible == 0, f"distancia {d60.min_distance:.1f} con 60 semillas, {d2.min_distance:.1f} con 2")
J3 = cm.joint(O, COV, seeds(60), nseeds=3).set_index("model").loc[names[0]]
check("se puede repetir con las primeras semillas", J3.seeds_per_point == 3)

print("7. Covarianza con una direccion redundante")
C = DIAG.copy(); C[4] = C[0]; C[:, 4] = C[:, 0]; C[4, 4] = C[0, 0]            # el quinto observable es una copia del primero
Ci, rank, cond = cm.region(C)
check("el rango efectivo es cinco, no seis", rank == 5, f"numero de condicion {cond:.1e}")
O = pd.DataFrame([cohort_row("A", a)]); J = cm.joint(O, {"A": C}, preds({m: [a] for m in names})).set_index("model")
check("y el umbral es el de cinco grados de libertad", abs(J.loc[names[0], "distance_threshold"] - 11.07) < 0.01 and J.loc[names[0], "covariance_effective_rank"] == 5)

print("8. La ventana de muestreo")
w2 = [w for w in cm.WINDOWS if w != cm.MORNING][0]
S = preds({m: [far] for m in names}); S2 = preds({names[0]: [a], names[1]: [far], names[2]: [far], names[3]: [far]}, window=w2)
both = pd.concat([S, S2])
Jm = cm.joint(pd.DataFrame([cohort_row("A", a, window="morning biopsy")]), {"A": DIAG}, both).set_index(["scope", "model"])
Jp = cm.joint(pd.DataFrame([cohort_row("A", a, window="post mortem, any hour")]), {"A": DIAG}, both).set_index(["scope", "model"])
check("una cohorte de manana solo admite la ventana de manana", Jm.loc[("windows allowed by the design", names[0]), "n_compatible"] == 0)
check("una cohorte post mortem admite cualquiera y dice cual", Jp.loc[("windows allowed by the design", names[0]), "n_compatible"] == 1
      and Jp.loc[("windows allowed by the design", names[0]), "closest_window"] == w2)

print("9. Una matriz que no define una region se rechaza, y dice por que")
def refuses(C):
    try: cm.region(C); return False
    except ValueError: return True
check("valores no finitos", refuses(np.array([[1.0, np.nan], [np.nan, 1.0]])))
check("ninguna direccion con varianza", refuses(np.zeros((3, 3))))
check("no definida positiva", refuses(np.array([[1.0, 2.0], [2.0, 1.0]])))
check("no cuadrada", refuses(np.ones((2, 3))))
near = np.eye(3) * 1e-2; near[2, 2] = 1e-9                   # una direccion casi sin varianza
Ci, rank, cond = cm.region(near)
check("casi singular: se descarta esa direccion y el numero de condicion es el de las que quedan", rank == 2 and abs(cond - 1.0) < 1e-9)
tiny = np.eye(3); tiny[2, 2] = -1e-14                        # negativo solo por redondeo
check("un autovalor negativo del tamano del redondeo se trata como cero", cm.region(tiny)[1] == 2)

print("10. La calibracion no usa lo que calibra, y respeta el tamano de la cohorte")
seeds2 = [0, 1]; pt = (cm.SD_F[0], cm.RHO_K[1], cm.PROXY[1], cm.MORNING)
ENS = dict(cm.ensemble_of((names[2], pt[0], pt[1], sd)) for sd in seeds2)
c = cm.calibrate(pd.Series(dict(cohort="A", n=57)), ENS, names[2], pt, seeds2, reps=6, boot=25, pools=2)
check("las cohortes de calibracion tienen el tamano de la observada", c["n_calibration_cohort"] == 57 and c["n_observed"] == 57)
check("salen de una poblacion varias veces mayor", c["population_size"] == 3 * 57)
check("de dos poblaciones independientes, con su intervalo", c["independent_populations"] == 2 and c["calibrated_threshold_min"] <= c["calibrated_threshold_max"])
check("con semillas que ninguna prediccion usa", cm.CALIBRATION_SEED > max(cm.SEEDS) and all(k[5] in (0, 1) for k in cm._POOLS))
a1 = cm.independent_population(names[2], pt, 30, 0); a2 = cm.independent_population(names[2], pt, 30, 1)
m0 = cm.measured(ENS, names[2], pt[0], pt[1], 0, pt[3], pt[2])
check("y que no repiten a las personas de la prediccion ni entre si",
      not np.allclose(a1[0][:30], m0[0][:30]) and not np.allclose(a1[0], a2[0]))

print("11. La cache solo vale si su manifiesto coincide")
M1 = cm.manifest([0, 1, 2]); M2 = cm.manifest([0, 1, 2, 3])
check("cambiar las semillas cambia la huella", M1["integrations_sha256"] != M2["integrations_sha256"] and M1["simulated_table_sha256"] != M2["simulated_table_sha256"])
check("la misma configuracion da la misma huella", M1 == cm.manifest([0, 1, 2]))
with tempfile.TemporaryDirectory() as d:
    os.makedirs(f"{d}/new"); old_out, cm.OUT = cm.OUT, d; f = f"{d}/new/x.pkl"; open(f, "w").write("x")
    os.environ["T2D_MODEL_REUSE"] = "1"
    no_manifest = cm.cached(f, "integrations_sha256", M1)
    json.dump(M2, open(f"{d}/new/competing_models_cache_manifest.json", "w")); other = cm.cached(f, "integrations_sha256", M1)
    json.dump(M1, open(f"{d}/new/competing_models_cache_manifest.json", "w")); same = cm.cached(f, "integrations_sha256", M1)
    os.environ["T2D_MODEL_REUSE"] = "0"; not_asked = cm.cached(f, "integrations_sha256", M1)
    cm.OUT = old_out; os.environ.pop("T2D_MODEL_REUSE")
check("sin manifiesto no se reutiliza", not no_manifest)
check("con un manifiesto de otra configuracion tampoco", not other)
check("con el mismo si, y solo si se pide", same and not not_asked)

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
