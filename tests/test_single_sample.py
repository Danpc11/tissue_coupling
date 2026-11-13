#!/usr/bin/env python3
"""Que significa "muestra unica" y que no: los estimadores son relativos a una cohorte, y una
muestra nueva se puntua contra una referencia fija.

Cohortes simuladas, sin datos.
  1. order_anchored() se niega a trabajar con una muestra aislada en vez de devolver NaN.
  2. Con una referencia fija, el valor de una muestra no depende de que otras la acompanen.
  3c. La referencia fija recupera la fase exacta aunque los dos vectores de su forma de onda no sean
     ortogonales ni de igual norma. Las fases se comparan por su ERROR circular, no por la
     concentracion de su diferencia, que daria por buena una fase desplazada 12 h.
  3. La referencia es estable: dos referencias ajustadas en mitades distintas de una cohorte dan
     casi la misma amplitud a las mismas muestras nuevas. Lo mismo para la fase anclada.
  4. La amplitud inferida sigue a la amplitud verdadera de cada muestra.
  5. is_clock(): acepta la antifase BMAL1-DBP, rechaza la fase comun, y encuentra BMAL1 bajo su
     simbolo anterior.
"""
import os, sys
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python/lib"))
import phase as ph
from scipy import stats

fails = []
def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: fails.append(name)

COHORT_PHASE = {}
def circ_err_h(a, b):
    """Error circular absoluto, en horas. A diferencia de la concentracion de la diferencia, es
    cero solo si las dos fases COINCIDEN: un desfase constante de 12 h da 12, no 0."""
    return np.abs(np.angle(np.exp(1j * (np.asarray(a) - np.asarray(b))))) * 24 / (2 * np.pi)

def cohort(rng, n, genes, peak, amp_sd=0.35, noise=0.5):
    t = rng.uniform(0, 2*np.pi, n); amp = np.clip(rng.normal(1, amp_sd, n), 0.1, None)
    X = amp[None, :] * np.cos(t[None, :] - peak[:, None]) + rng.normal(0, noise, (len(genes), n)) + 8
    X = pd.DataFrame(X, index=genes, columns=[f"s{i}" for i in range(n)])
    COHORT_PHASE["last"] = t
    return X.sub(X.mean(axis=0), axis=1), amp

rng = np.random.default_rng(0)
genes = list(ph.CLK); peak = rng.uniform(0, 2*np.pi, len(genes))
gi = {g: i for i, g in enumerate(genes)}; peak[gi["DBP"]] = (peak[gi["BMAL1"]] + np.pi) % (2*np.pi)   # antifase
E, amp = cohort(rng, 300, genes, peak); true_phase = COHORT_PHASE["last"]

print("1. Una muestra aislada")
tpl = pd.DataFrame({"r2": 1.0, "phase_h": peak * 24 / (2*np.pi), "amp": 1.0}, index=genes)
try:
    ph.order_anchored(E.iloc[:, :1], tpl, genes=genes); raised = False
except ValueError:
    raised = True
check("order_anchored avisa en vez de devolver valores no numericos", raised)

print("2. Con referencia fija, la muestra no depende de las demas")
ref = ph.fit_reference(E.iloc[:, :200])
new = E.iloc[:, 200:]
alone = np.array([ph.infer_from_reference(new[c], ref)[1] for c in new.columns])
check("amplitud finita para cada muestra nueva, una a una", bool(np.isfinite(alone).all()), f"{len(alone)} muestras")
shuffled = np.array([ph.infer_from_reference(new[c], ref)[1] for c in new.columns[::-1]])[::-1]
check("el orden y la compania no cambian el valor", bool(np.array_equal(alone, shuffled)))

print("3. Estabilidad de la referencia")
ref2 = ph.fit_reference(E.iloc[:, 100:200])
other = np.array([ph.infer_from_reference(new[c], ref2)[1] for c in new.columns])
r = stats.spearmanr(alone, other).statistic
check("dos referencias de cohortes distintas ordenan igual a las muestras nuevas", r > 0.9, f"rho = {r:.3f}")

print("3b. Fase anclada con referencia fija")
aref = ph.fit_anchored_reference(E.iloc[:, :200], tpl, genes=genes)
one = np.array([ph.phase_from_anchored_reference(new[c], aref) for c in new.columns])
check("fase finita para cada muestra nueva, una a una", bool(np.isfinite(one).all()))
cohort_phase = ph.order_anchored(E, tpl, genes=genes)[0].values[200:]          # la misma muestra, dentro de la cohorte
e = circ_err_h(one, cohort_phase)
check("coincide con la fase que recibe dentro de la cohorte", np.median(e) < 0.5, f"diferencia mediana {np.median(e):.2f} h, maxima {e.max():.2f} h")
aref2 = ph.fit_anchored_reference(E.iloc[:, 100:200], tpl, genes=genes)
two = np.array([ph.phase_from_anchored_reference(new[c], aref2) for c in new.columns])
e = circ_err_h(one, two)
check("y no depende de que referencia se use", np.median(e) < 0.5, f"diferencia mediana {np.median(e):.2f} h")
e = circ_err_h(one, true_phase[200:])
check("y es la hora verdadera de la muestra, no solo una fase coherente", np.median(e) < 1.0 and np.mean(e <= 3) > 0.9,
      f"error mediano {np.median(e):.2f} h, {np.mean(e <= 3) * 100:.0f} % dentro de 3 h")
check("un desfase constante de 12 h no pasaria esta prueba", np.median(circ_err_h(one + np.pi, true_phase[200:])) > 10)

print("3c. Geometria de la referencia: los vectores de la forma de onda no son ortonormales")
g2 = genes[:8]; rng2 = np.random.default_rng(3)
fake = dict(genes=g2, mu=np.zeros(8), sd=np.ones(8), a=2.0 * rng2.normal(0, 1, 8), b=0.7 * rng2.normal(0, 1, 8))
thetas = np.linspace(0, 2 * np.pi, 48, endpoint=False)
ls, pub, amp_ls, amp_pub = [], [], [], []
for th in thetas:
    x = pd.Series(1.3 * (fake["a"] * np.cos(th) + fake["b"] * np.sin(th)), index=g2)
    p1, a1 = ph.infer_from_reference(x, fake); p0, a0 = ph.infer_from_reference(x, fake, geometry="as published")
    ls.append(p1); pub.append(p0); amp_ls.append(a1); amp_pub.append(a0)
check("minimos cuadrados recupera la fase exacta", circ_err_h(ls, thetas).max() < 1e-6, f"error maximo {circ_err_h(ls, thetas).max():.1e} h")
check("las proyecciones directas la sesgan", circ_err_h(pub, thetas).max() > 1.0, f"error maximo {circ_err_h(pub, thetas).max():.1f} h")
degenerate = dict(fake, b=fake["a"] * 0.5)
try:
    ph.infer_from_reference(pd.Series(fake["a"], index=g2), degenerate); raised = False
except ValueError:
    raised = True
check("una referencia sin estructura ciclica se rechaza", raised)

print("4. La amplitud inferida sigue a la verdadera")
r = stats.spearmanr(alone, amp[200:]).statistic
check("correlacion con la amplitud simulada de cada muestra", r > 0.5, f"rho = {r:.2f}")

print("5. Control de antifase")
Z = E[genes] if False else E
ang = ph.unsup(((Z.T - Z.mean(axis=1)) / (Z.std(axis=1) + 1e-9)).T.values)
ok, sep = ph.is_clock(E, ang)
check("BMAL1 y DBP en antifase: se acepta", ok and sep > 8, f"separacion {sep:.1f} h")
peak2 = peak.copy(); peak2[gi["DBP"]] = peak2[gi["BMAL1"]]
E2, _ = cohort(rng, 300, genes, peak2)
ang2 = ph.unsup(((E2.T - E2.mean(axis=1)) / (E2.std(axis=1) + 1e-9)).T.values)
ok2, sep2 = ph.is_clock(E2, ang2)
check("BMAL1 y DBP en fase: se rechaza", (not ok2) and sep2 < 4, f"separacion {sep2:.1f} h")
E3 = E.rename(index={"BMAL1": "ARNTL"})
ok3, sep3 = ph.is_clock(E3, ang)
check("BMAL1 anotado como ARNTL: se encuentra", ok3 and abs(sep3 - sep) < 1e-9)

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
