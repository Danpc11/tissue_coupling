#!/usr/bin/env python3
"""La validacion de la plantilla de fase no usa nada del sujeto que se deja fuera.

Cohortes simuladas, sin datos. Tres comprobaciones:

  1. Independencia: lo que un pliegue aprende (parametros de estandarizacion, orden de los genes)
     no cambia si las muestras del sujeto retenido se sustituyen por ruido.
  2. Sin ritmo no hay acierto: en una cohorte de puro ruido el error dejando sujetos fuera es el
     del azar (mediana de 6 h). El procedimiento anterior a v2.2.2, que elegia los genes y
     estandarizaba con todas las muestras, da ahi un error claramente menor: esa diferencia es la
     fuga, y se mide aqui con una reimplementacion de aquel procedimiento.
  3. Con ritmo, acierta: en una cohorte con genes ritmicos el error es pequeno, tambien en la
     validacion anidada, que elige la configuracion sin ver al sujeto retenido.
"""
import os, sys, importlib.util
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("phase_template", os.path.join(ROOT, "python/analyses/02b_phase_template.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

fails = []
def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: fails.append(name)

def cohort(rng, n_genes, n_rhythmic, noise=0.7):
    subject = np.repeat([f"s{i}" for i in range(10)], 6)
    hours = np.tile([0., 4., 8., 12., 16., 20.], 10)
    t = hours * 2 * np.pi / 24
    X = rng.normal(0, noise, (n_genes, len(t))) + rng.normal(0, 1, (n_genes, 1)) + 8
    peak = rng.uniform(0, 2 * np.pi, n_rhythmic)
    X[:n_rhythmic] += np.cos(t[None, :] - peak[:, None])
    return X, t, subject, np.ones(n_genes, bool)

def old_procedure(X, t, subject, harm, n):
    """Como se hacia hasta v2.2.1: seleccion y estandarizacion con todas las muestras."""
    Z = (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-9)
    sel = np.argsort(-m.fit_cosinor(Z, t, harm)[0])[:n]; e = []
    for s in pd.unique(subject):
        te = subject == s
        B, *_ = np.linalg.lstsq(m.design(t[~te], harm, intercept=False), Z[sel][:, ~te].T, rcond=None)
        e.extend(m.circular_error_h(m.predict(Z[sel][:, te], B, harm), t[te]))
    return np.array(e)

rng = np.random.default_rng(0)

print("1. El pliegue no depende de las muestras del sujeto retenido")
X, t, subject, ok = cohort(rng, 1500, 150)
te = subject == "s3"
a = m.Fold(X, t, ~te, ok)
X2 = X.copy(); X2[:, te] = rng.normal(0, 5, (X.shape[0], int(te.sum())))
b = m.Fold(X2, t, ~te, ok)
check("misma media y desviacion por gen", np.array_equal(a.mu, b.mu) and np.array_equal(a.sd, b.sd))
check("mismo orden de genes, con uno y con dos armonicos",
      a.ranking([24.0]) == b.ranking([24.0]) and a.ranking([24.0, 12.0]) == b.ranking([24.0, 12.0]))

print("2. Cohorte de puro ruido")
X, t, subject, ok = cohort(rng, 3000, 0)
clean = m.loso_by_config(X, t, subject, ok)
c30 = next(c for c, (h, lab, n) in enumerate(m.CONFIGS) if lab == "24 h" and n == 30)
med_clean = float(np.median(np.concatenate(list(clean.values()))))
med_old = float(np.median(old_procedure(X, t, subject, [24.0], 30)))
check("sin fuga, el error es el del azar", 4.5 <= med_clean <= 7.5, f"mediana {med_clean:.2f} h; azar 6 h")
check("el procedimiento anterior parecia acertar sobre ruido", med_old < med_clean - 1.0,
      f"{med_old:.2f} h con 30 genes frente a {np.median(clean[c30]):.2f} h sin fuga")

print("3. Cohorte con genes ritmicos")
X, t, subject, ok = cohort(rng, 2000, 200)
clean = m.loso_by_config(X, t, subject, ok)
best = min(np.median(e) for e in clean.values())
check("sin fuga, error pequeno", best < 1.5, f"mejor configuracion {best:.2f} h")
N, e = m.nested_loso(X, t, subject, ok)
check("validacion anidada: una fila por sujeto y error pequeno", len(N) == 10 and len(e) == 60 and np.median(e) < 1.5,
      f"mediana {np.median(e):.2f} h")
check("la configuracion elegida es una de las doce", set(zip(N.chosen_harmonics, N.chosen_n_genes)) <= {(lab, n) for _, lab, n in m.CONFIGS})

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
