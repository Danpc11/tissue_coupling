#!/usr/bin/env python3
"""El estimador aplicado a cohortes simuladas: recupera el efecto cuando existe y no lo inventa
cuando no existe.

La figura del modelo aplica el estimador a cohortes generadas por un modelo en el que la ganancia
del receptor K controla el forzamiento. Eso solo dice que el estimador ve el efecto cuando esta.
Aqui se anade el caso contrario, con el mismo generador:

  1. con forzamiento (F > 0), la amplitud inferida crece con K;
  2. sin forzamiento (F = 0), K varia igual entre personas pero no actua sobre el oscilador, y la
     amplitud inferida no debe relacionarse con K;
  3. con K barajada entre personas despues de simular, tampoco.

Sin datos. El modelo demuestra que el mecanismo propuesto basta para producir la relacion; no
demuestra que sea el mecanismo real. Este test solo protege contra un estimador que encuentre
relacion en cualquier cosa.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python/model")); sys.path.insert(0, os.path.join(ROOT, "python/lib"))
import forced_clock as fc

fails = []
def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: fails.append(name)

def cohort(F, seed, n=90):
    rng = np.random.default_rng(seed); K = np.clip(rng.normal(1.0, 0.35, n), 0.05, None); hour = rng.uniform(0, 24, n)
    rows = []
    for i, k in enumerate(K):
        phi, r = fc.simulate_population(F, k, n_cells=40, days=12, seed=int(rng.integers(1e6)))
        a, p, R = fc.ensemble(phi, r); rows.append(dict(subject=i, K=k, hour=hour[i], amp_ens=a, phase_ens=p, R=R))
    d = pd.DataFrame(rows)
    E = fc.synthetic_expression(d, noise=0.3, seed=int(rng.integers(1e6)))
    _, amp = fc.run_estimator(E)
    return d, np.asarray(amp)

print("1. Con forzamiento")
d, amp = cohort(0.04, 1)
r = stats.spearmanr(d.K, amp)
check("la amplitud inferida crece con la ganancia del receptor", r.statistic > 0.3 and r.pvalue < 0.01, f"rho = {r.statistic:+.2f}, P = {r.pvalue:.1e}")

print("2. Sin forzamiento: K no actua")
d0, amp0 = cohort(0.0, 2)
r0 = stats.spearmanr(d0.K, amp0)
check("no aparece relacion con K", abs(r0.statistic) < 0.25 and r0.pvalue > 0.01, f"rho = {r0.statistic:+.2f}, P = {r0.pvalue:.2f}")

print("3. K barajada entre personas")
rs = [stats.spearmanr(np.random.default_rng(k).permutation(d.K.values), amp).statistic for k in range(200)]
check("la relacion desaparece al romper la correspondencia", abs(np.mean(rs)) < 0.05 and np.percentile(np.abs(rs), 95) < r.statistic,
      f"media {np.mean(rs):+.3f}; el 95 % por debajo de {np.percentile(np.abs(rs), 95):.2f}")

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
