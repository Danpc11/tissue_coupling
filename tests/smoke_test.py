#!/usr/bin/env python3
"""Prueba de humo: comprueba que las piezas del analisis funcionan sobre datos sinteticos
con respuesta conocida. No sustituye a verify_package, que comprueba la coherencia de una
corrida real; esta comprueba que el codigo hace lo que dice."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../python/lib"))
import gene_sets as G
from phase import unsup, align, kap, rho_js, CLK
F=[]
def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: F.append(name)

print("1. Registro de conjuntos")
check("se importa", len(G.SETS) >= 6, f"{len(G.SETS)} conjuntos")
check("cada conjunto declara fuente", all(len(v["source"]) > 20 for v in G.SETS.values()))
check("resuelve sinonimos", G.resolve("KDM5D", ["JARID1D"]) == "JARID1D")
check("devuelve None si falta", G.resolve("KDM5D", ["ACTB"]) is None)

print("\n2. Estadistico circular")
rng = np.random.default_rng(0); a = rng.uniform(0, 2*np.pi, 300)
check("invariante a rotacion", abs(abs(rho_js(a, a + 1.3)) - 1) < 1e-9)
check("invariante a reflexion", abs(abs(rho_js(a, -a)) - 1) < 1e-9)
check("nulo con series independientes", abs(rho_js(a, rng.uniform(0, 2*np.pi, 300))) < 0.2)

print("\n3. Ordenamiento no supervisado recupera una fase conocida")
n, ng = 300, 14
true = rng.uniform(0, 2*np.pi, n)
ph = rng.uniform(0, 2*np.pi, ng)
Z = np.array([np.cos(true - p) for p in ph]) + rng.normal(0, 0.35, (ng, n))
Z = (Z - Z.mean(1, keepdims=True)) / Z.std(1, keepdims=True)
est = unsup(Z)
r = abs(rho_js(est, true))
check("recupera la fase simulada", r > 0.85, f"|rho| = {r:.3f}")

print("\n4. align corrige rotacion y reflexion")
# una copia reflejada y rotada de la misma fase: sin alinear la discrepancia es grande,
# alineando debe ser practicamente nula. Es la situacion real al comparar dos ordenaciones.
b = -true + 2.1
d_raw = np.angle(np.exp(1j*(true - b)))
d_al = align(true, b)
def circ_sd_h(d):
    R = np.abs(np.mean(np.exp(1j*d)))
    return np.sqrt(-2*np.log(max(R,1e-12)))*24/(2*np.pi)
sd_raw = circ_sd_h(d_raw - np.angle(np.mean(np.exp(1j*d_raw))))
sd_al  = circ_sd_h(d_al)
check("alinear elimina la indeterminacion", sd_al < 0.1 and sd_raw > 2,
      f"sin alinear {sd_raw:.1f} h -> alineado {sd_al:.2f} h")

print("\n5. kap es estable en los extremos")
k1,R1 = kap(np.zeros(100))
k2,R2 = kap(rng.uniform(0,2*np.pi,2000))
check("R = 1 no produce nan", np.isfinite(k1), f"kappa = {k1:.0f}")
check("R = 0 no produce nan", np.isfinite(k2) and k2 < 1, f"kappa = {k2:.3f}")

print("\n6. la plantilla de fase selecciona los mas ritmicos, no todos")
from phase import select_template, order_anchored
tpl = pd.DataFrame({"r2": np.linspace(0.9, 0.01, 500),
                    "amp": np.linspace(1.2, 0.05, 500),
                    "phase_h": rng.uniform(0, 24, 500)},
                   index=[f"G{i:03d}" for i in range(495)] + ["LOC1", "LINC2", "C1orf3", "LOC4", "G999"])
sel = select_template(tpl, 60)
check("devuelve exactamente n genes", len(sel) == 60, f"{len(sel)}")
check("ordenados por ritmicidad", sel[0] == "G000" and tpl.r2[sel].is_monotonic_decreasing)
check("descarta identificadores sin anotar",
      not any(str(g).startswith(("LOC", "LINC")) or "orf" in str(g).lower() for g in select_template(tpl, 500)))
check("por defecto NO usa la plantilla completa", len(select_template(tpl)) == 60, "n por defecto = 60")
# el anclado debe usar la seleccion, no los 500
E = pd.DataFrame(rng.normal(size=(500, 120)), index=tpl.index)
ang, used = order_anchored(E, tpl)
check("order_anchored usa la seleccion", len(used) == 60, f"{len(used)} genes")

print("\n" + ("TODO CORRECTO" if not F else f"FALLAN {len(F)}: {F}"))
sys.exit(1 if F else 0)
