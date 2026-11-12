#!/usr/bin/env python3
"""Que mide y que no mide la concentracion de la diferencia de fase entre tejidos.

Fases simuladas, sin datos. Cuatro situaciones:

  1. tejidos independientes            -> R cerca de cero, P no significativo
  2. tejidos alineados                 -> R alto
  3. senal comun sin interaccion       -> R alto TAMBIEN: cada tejido sigue por su cuenta una senal
                                          sistemica y no hay ningun termino que los conecte. R mide
                                          coherencia de fase, no acoplamiento directo.
  4. un tejido recorrido al reves      -> R de la diferencia cae a cero aunque esten alineados; la
                                          version que no supone el sentido de giro lo recupera. Es
                                          lo que ocurre con un estimador sin ancla como CHIRAL.

Y dos propiedades: rotar todas las fases de un tejido no cambia R, y la version que no supone el
sentido de giro mantiene el 5 % de falsos positivos.
"""
import os, sys, importlib.util
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python/lib"))
spec = importlib.util.spec_from_file_location("chiral_concordance", os.path.join(ROOT, "python/analyses/17_chiral_concordance.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

fails = []
def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: fails.append(name)

rng = np.random.default_rng(0); n = 300
wrap = lambda x: np.mod(x, 2*np.pi)

print("1. Tejidos independientes")
a, b = rng.uniform(0, 2*np.pi, n), rng.uniform(0, 2*np.pi, n)
R, p = m.concentration(a, b)
check("R bajo y no significativo", R < 0.15 and p > 0.05, f"R = {R:.3f}, P = {p:.3f}")

print("2. Tejidos alineados")
a = rng.uniform(0, 2*np.pi, n); b = wrap(a + 0.8 + rng.normal(0, 0.4, n))
R, p = m.concentration(a, b)
check("R alto y significativo", R > 0.8 and p < 0.01, f"R = {R:.3f}, P = {p:.4f}")

print("3. Senal comun, sin interaccion entre tejidos")
signal = rng.uniform(0, 2*np.pi, n)                      # fase de la senal sistemica de cada persona
a = wrap(signal + 0.3 + rng.normal(0, 0.4, n))           # cada tejido la sigue con su retraso y su ruido
b = wrap(signal + 1.4 + rng.normal(0, 0.4, n))           # no hay termino que una a con b
R, p = m.concentration(a, b)
check("R alto sin acoplamiento directo", R > 0.6 and p < 0.01, f"R = {R:.3f}, P = {p:.4f}")

print("4. Un tejido recorrido en sentido contrario")
a = rng.uniform(0, 2*np.pi, n); b = wrap(-(a + 0.8 + rng.normal(0, 0.4, n)))
R, p = m.concentration(a, b)
R2, p2 = m.concentration_reflection_aware(a, b)
check("la diferencia no lo detecta", R < 0.15, f"R = {R:.3f}, P = {p:.3f}")
check("la version sin sentido fijo si", R2 > 0.8 and p2 < 0.01, f"R = {R2:.3f}, P = {p2:.4f}")
# tomar el maximo de dos estadisticos no debe inflar los falsos positivos: el nulo usa el mismo maximo
r2 = np.random.default_rng(1); ps = []
for k in range(300):
    a, b = r2.uniform(0, 2*np.pi, 120), r2.uniform(0, 2*np.pi, 120)
    ps.append(m.concentration_reflection_aware(a, b, seed=k, B=200)[1])
frac = float(np.mean(np.array(ps) < 0.05))
check("y con tejidos independientes da el 5 % de falsos positivos", 0.02 <= frac <= 0.09, f"{frac:.3f} en 300 replicas")

print("5. Invariancia a la rotacion")
a = rng.uniform(0, 2*np.pi, n); b = wrap(a + rng.normal(0, 0.5, n))
check("rotar un tejido no cambia R", abs(m.concentration(a, b)[0] - m.concentration(a, wrap(b + 2.1))[0]) < 1e-9)

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
