#!/usr/bin/env python3
"""El criterio con el que se comparan los modelos competidores, comprobado a mano.

Sin simulaciones ni datos: se construyen lo observado y lo simulado.

  1. Tres estados por observable: una asociacion negativa no cuenta como ausente.
  2. El criterio conjunto va cohorte por cohorte. Una simulacion que acierta un observable de una
     cohorte y otro de otra no es compatible con ninguna de las dos.
  3. Es una region conjunta: respeta la correlacion entre observables. Una simulacion dentro de
     cada intervalo por separado puede quedar fuera del conjunto.
  4. La concordancia del receptor entre tejidos cuenta.
  5. Se cuenta cuantas veces, no si alguna vez.
  6. La prediccion de un punto de la rejilla es la media de sus semillas: mas semillas hacen el
     criterio mas exigente, no mas laxo.
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

def sims(vectors_by_model, window="morning biopsy"):
    """Cada vector es la prediccion de un punto distinto de la rejilla, con muchas personas para que
    el ruido de la simulacion no cuente en estas comprobaciones."""
    rows = []
    for model, vectors in vectors_by_model.items():
        for i, v in enumerate(vectors):
            rows.append(dict(model=model, window=window, sd_F=0.25, rho_K=0.3, proxy_noise=float(i), seed=0, n=10**6, **dict(zip(K, v))))
    return pd.DataFrame(rows)

se = 0.08 / 1.96; DIAG = np.eye(len(K)) * se ** 2
blank = {m: [[0.0] * 6] for m in names}

print("1. Tres estados por observable")
O = pd.DataFrame([cohort_row("A", [-0.4] * 6)])
S = sims({names[0]: [[0.0] * 6], names[1]: [[0.4] * 6], names[2]: [[-0.42] * 6], names[3]: [[0.0] * 6]})
V = cm.verdict(O, S).set_index(["model", "observable"])
check("lo negativo se clasifica como negativo", V.loc[(names[0], K[0]), "observed_state"] == "negative")
check("un modelo que solo da cero no lo situa en el intervalo", V.loc[(names[0], K[0]), "never_within_ci"] and not V.loc[(names[0], K[0]), "state_reproduced"])
check("uno que solo da positivo tampoco", bool(V.loc[(names[1], K[0]), "never_within_ci"]))
check("uno que da negativo si", V.loc[(names[2], K[0]), "fraction_of_grid_points_within_ci"] == 1.0)

print("2. Cohorte por cohorte, sin mezclar")
a = [0.20, 0.40, 0.25, 0.10, 0.20, 0.35]; b = [0.50, 0.10, 0.25, 0.30, 0.45, 0.05]
O = pd.DataFrame([cohort_row("A", a), cohort_row("B", b)]); COV = {"A": DIAG, "B": DIAG}
mix = [a[0], b[1], a[2], b[3], a[4], b[5]]                 # un observable de cada cohorte
S = sims({names[0]: [a], names[1]: [b], names[2]: [mix], names[3]: [[0.0] * 6]})
J = cm.joint(O, COV, S).set_index(["cohort", "model"])
check("la simulacion igual a la cohorte A es compatible con A y no con B",
      J.loc[("A", names[0]), "n_compatible"] == 1 and J.loc[("B", names[0]), "n_compatible"] == 0)
check("la mezcla de las dos no es compatible con ninguna", J.loc[("A", names[2]), "n_compatible"] == 0 and J.loc[("B", names[2]), "n_compatible"] == 0,
      "con un intervalo envolvente de las dos cohortes habria pasado")

print("3. Region conjunta, no seis intervalos")
C = DIAG.copy(); C[0, 4] = C[4, 0] = 0.98 * se ** 2          # O1 y O5 casi la misma cantidad
O = pd.DataFrame([cohort_row("A", a)]); COV = {"A": C}
along = list(a); along[0] += 0.05; along[4] += 0.05          # se mueven juntos, como manda la correlacion
against = list(a); against[0] += 0.05; against[4] -= 0.05    # dentro de cada intervalo, pero en sentido contrario
S = sims({names[0]: [along], names[1]: [against], names[2]: [a], names[3]: [[0.0] * 6]})
J = cm.joint(O, COV, S).set_index("model")
check("moverse a favor de la correlacion es compatible", J.loc[names[0], "n_compatible"] == 1, f"distancia {J.loc[names[0], 'min_distance']:.1f}")
check("moverse en contra no lo es, aunque cada valor este dentro de su intervalo",
      J.loc[names[1], "n_compatible"] == 0 and J.loc[names[1], "n_all_within_simultaneous_ci"] == 1, f"distancia {J.loc[names[1], 'min_distance']:.1f}")

print("4. La concordancia del receptor entre tejidos cuenta")
O = pd.DataFrame([cohort_row("A", a)]); COV = {"A": DIAG}
no_o3 = list(a); no_o3[2] = 0.0                              # todo igual salvo que el receptor no se parece entre tejidos
J = cm.joint(O, COV, sims({names[0]: [no_o3], names[1]: [a], names[2]: [a], names[3]: [a]})).set_index("model")
check("una simulacion sin esa concordancia no es compatible", J.loc[names[0], "n_compatible"] == 0 and J.loc[names[1], "n_compatible"] == 1)

print("5. Cuantas veces, no si alguna vez")
far = [0.9] * 6
J = cm.joint(O, COV, sims({names[0]: [a] * 19 + [far], names[1]: [far] * 19 + [a], names[2]: [far] * 20, names[3]: [a] * 20})).set_index("model")
check("casi siempre: 0.95", J.loc[names[0], "fraction_compatible"] == 0.95)
check("una de veinte: 0.05", J.loc[names[1], "fraction_compatible"] == 0.05)
check("nunca: 0, y la distancia minima lo dice", J.loc[names[2], "fraction_compatible"] == 0.0 and J.loc[names[2], "min_distance"] > J.loc[names[2], "distance_threshold"])

print("6. La prediccion de un punto es la media de sus semillas, y mas semillas exigen mas")
rng = np.random.default_rng(0); truth = np.array([0.05, 0.40, 0.25, 0.10, 0.05, 0.35])      # un modelo cuyo O1 verdadero es 0.05
def seeds(k):
    rows = []
    for sd in range(k):
        v = np.tanh(np.arctanh(truth) + rng.normal(0, 1 / np.sqrt(200 - 3), 6))             # cohortes simuladas de 200 personas
        rows.append(dict(model=names[0], window="morning biopsy", sd_F=0.25, rho_K=0.3, proxy_noise=1.0, seed=sd, n=200, **dict(zip(K, v))))
    return pd.DataFrame(rows + [dict(model=m, window="morning biopsy", sd_F=0.25, rho_K=0.3, proxy_noise=1.0, seed=0, n=200, **dict(zip(K, [0.0] * 6))) for m in names[1:]])
O = pd.DataFrame([cohort_row("A", a, half=0.06, n=1000)]); COV = {"A": np.eye(6) * (0.06 / 1.96) ** 2}
d_few = cm.joint(O, COV, seeds(2)).set_index("model").loc[names[0]]
d_many = cm.joint(O, COV, seeds(60)).set_index("model").loc[names[0]]
check("un punto, sea cual sea el numero de semillas", d_few.n_grid_points == 1 and d_many.n_grid_points == 1 and d_many.seeds_per_point == 60)
check("con muchas semillas el modelo equivocado queda descartado", d_many.n_compatible == 0, f"distancia {d_many.min_distance:.1f} con 60 semillas, {d_few.min_distance:.1f} con 2")

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
