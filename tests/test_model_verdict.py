#!/usr/bin/env python3
"""El veredicto de los modelos competidores distingue tres estados: positivo, negativo y compatible
con cero. Una asociacion negativa intensa no puede contarse como ausente.

Sin simulaciones ni datos: se construyen a mano lo observado y lo simulado.
"""
import os, sys, importlib.util
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("cm", os.path.join(ROOT, "python/model/competing_models.py"))
cm = importlib.util.module_from_spec(spec); spec.loader.exec_module(cm)

fails = []
def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: fails.append(name)

def observed(value, lo, hi, window="morning biopsy"):
    r = dict(cohort="c", window=window, n=200)
    for k in cm.OBS: r[k], r[f"{k} lo"], r[f"{k} hi"] = value, lo, hi
    return pd.DataFrame([r])

def simulated(values_by_model, window="morning biopsy"):
    rows = []
    for model in cm.MODELS:
        for v in values_by_model[model]:
            rows.append(dict(model=model, window=window, **{k: v for k in cm.OBS}))
    return pd.DataFrame(rows)

names = list(cm.MODELS)
S = simulated({names[0]: [0.0, 0.05], names[1]: [0.4, 0.6], names[2]: [-0.5, -0.3], names[3]: [-0.05, 0.02, 0.5]})

print("1. Lo observado es NEGATIVO")
V = cm.verdict(observed(-0.4, -0.55, -0.25), S).set_index(["model", "observable"])
k = cm.OBS[0]
check("se clasifica como negativo, no como ausente", V.loc[(names[0], k), "observed_state"] == "negative")
check("un modelo que solo da cero no lo reproduce", not V.loc[(names[0], k), "model_can_reproduce"])
check("un modelo que solo da positivo no lo reproduce", not V.loc[(names[1], k), "model_can_reproduce"])
check("un modelo que da negativo en ese intervalo si", bool(V.loc[(names[2], k), "model_can_reproduce"]))

print("2. Lo observado es compatible con CERO")
V = cm.verdict(observed(0.02, -0.12, 0.16), S).set_index(["model", "observable"])
check("se clasifica como cero", V.loc[(names[0], k), "observed_state"] == "zero")
check("un modelo siempre negativo no lo reproduce", not V.loc[(names[2], k), "model_can_reproduce"])
check("un modelo que puede dar cero si", bool(V.loc[(names[0], k), "model_can_reproduce"]) and bool(V.loc[(names[3], k), "model_can_reproduce"]))

print("3. Lo observado es POSITIVO, y el intervalo cuenta")
V = cm.verdict(observed(0.45, 0.35, 0.55), S).set_index(["model", "observable"])
check("un modelo positivo dentro del intervalo lo reproduce", bool(V.loc[(names[1], k), "model_can_reproduce"]))
S2 = simulated({names[0]: [0.0], names[1]: [0.9, 0.95], names[2]: [-0.5], names[3]: [0.2, 0.25]})
V = cm.verdict(observed(0.45, 0.35, 0.55), S2).set_index(["model", "observable"])
check("un modelo positivo pero lejos del intervalo observado no", not V.loc[(names[1], k), "model_can_reproduce"]
      and bool(V.loc[(names[1], k), "state_reproduced"]) and not V.loc[(names[1], k), "interval_overlaps"])

print("4. Cuantas veces, no si alguna vez")
S3 = simulated({names[0]: [0.0] * 20, names[1]: [0.45] * 19 + [-0.9], names[2]: [-0.5] * 19 + [0.45], names[3]: [0.9] * 20})
J = cm.joint(observed(0.45, 0.35, 0.55), S3).set_index("model")
check("un modelo que cae dentro casi siempre", J.loc[names[1], "fraction_all_within_ci"] == 0.95)
check("uno que solo lo roza una vez de veinte", J.loc[names[2], "fraction_all_within_ci"] == 0.05,
      "con el criterio 'alguna vez' los dos pasarian igual")
check("uno que no cae nunca", J.loc[names[3], "fraction_all_within_ci"] == 0.0 and J.loc[names[0], "fraction_all_within_ci"] == 0.0)

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
