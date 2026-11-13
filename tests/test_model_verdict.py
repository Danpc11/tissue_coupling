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
"""
import os, sys, importlib.util
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

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
