#!/usr/bin/env python3
"""Que mide el estimador de amplitud, y que no puede separar: una simulacion factorial.

Una biopsia es una mezcla de celulas. Cuatro cosas distintas bajan la oscilacion que se ve en ella,
y aqui se varian por separado, de forma independiente entre muestras:

    a      amplitud de la oscilacion de CADA celula
    sigma  dispersion de fase entre las celulas de la muestra. El vector medio de celulas con fases
           dispersas tiene longitud R = exp(-sigma^2 / 2): menos sincronia, menos amplitud colectiva,
           aunque cada celula oscile igual
    c      fraccion de celulas ritmicas. El resto es otro tipo celular, sin ritmo y con otro nivel
           basal de cada gen (composicion del tejido)
    ruido  de medida

La expresion del gen g en la muestra i es

    x_gi = c_i * a_i * R(sigma_i) * cos(theta_i - phi_g)  +  (1 - c_i) * delta_g  +  ruido

con theta_i la fase de la muestra (desconocida, uniforme) y delta_g el desplazamiento basal del
segundo tipo celular. La amplitud colectiva verdadera es el producto A_i = c_i * a_i * R(sigma_i).
A la matriz se le aplica el estimador del articulo (ordenamiento libre y proyeccion) y se pregunta:

  1. Sigue la amplitud inferida a la amplitud colectiva A?
  2. Distingue entre sus tres componentes? Si dos muestras tienen la misma A por caminos distintos
     (celulas mas debiles, o celulas igual de fuertes pero desincronizadas, o menos celulas
     ritmicas), las separa?
  3. La confunde un cambio de composicion que no toca el reloj (delta distinto de cero)?
  4. Que devuelve cuando no hay ritmo?

Sin datos. Salidas: results/new/factorial_estimator.tsv (una fila por escenario) y
results/new/factorial_estimator_samples.tsv (las muestras del escenario de referencia).
"""
import os, sys, warnings; warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
import phase as ph
OUT = os.environ.get("T2D_OUT", "results")
N, SEED = 600, 0


def cohort(rng, noise=0.5, delta_sd=0.0, rhythm=True, n=N, genes=None):
    genes = genes or list(ph.CLK); G = len(genes)
    phi = rng.uniform(0, 2 * np.pi, G)
    a = rng.uniform(0.5, 1.5, n)                         # amplitud de cada celula
    sigma = rng.uniform(0.0, 1.5, n)                     # dispersion de fase entre celulas (rad)
    c = rng.uniform(0.5, 1.0, n)                         # fraccion de celulas ritmicas
    theta = rng.uniform(0, 2 * np.pi, n)
    R = np.exp(-sigma**2 / 2)
    A = c * a * R if rhythm else np.zeros(n)
    delta = rng.normal(0, delta_sd, G)
    X = A[None, :] * np.cos(theta[None, :] - phi[:, None]) + (1 - c)[None, :] * delta[:, None] + rng.normal(0, noise, (G, n))
    E = pd.DataFrame(X, index=genes, columns=[f"s{i}" for i in range(n)])
    return E, pd.DataFrame(dict(a=a, sigma=sigma, R=R, c=c, A=A, theta=theta), index=E.columns)


def infer(E):
    Z = ((E.T - E.mean(1)) / (E.std(1) + 1e-9)).T
    ang = ph.unsup(Z.values)
    return np.asarray(ph.amplitude(E, ang, genes=list(E.index))), ang


def partial(x, y, Z):
    """Correlacion de rangos de x con y una vez descontado Z de los dos."""
    r = lambda v: stats.rankdata(v)
    D = np.column_stack([np.ones(len(x))] + [r(z) for z in np.atleast_2d(Z)])
    res = lambda v: r(v) - D @ np.linalg.lstsq(D, r(v), rcond=None)[0]
    return stats.pearsonr(res(x), res(y))


def summarise(label, E, T, rows, keep=None):
    amp, ang = infer(E)
    sp = lambda u: stats.spearmanr(amp, u).statistic
    row = dict(scenario=label, n=len(T), rho_collective_amplitude=sp(T.A) if T.A.std() > 0 else np.nan,
               rho_cell_amplitude=sp(T.a), rho_synchrony=sp(T.R), rho_rhythmic_fraction=sp(T.c),
               median_inferred_amplitude=float(np.median(amp)))
    if T.A.std() > 0:
        # con la amplitud colectiva descontada, queda algo de cada componente?
        for k, col in (("cell_amplitude", "a"), ("synchrony", "R"), ("rhythmic_fraction", "c")):
            r = partial(amp, T[col].values, T.A.values)
            row[f"partial_{k}_given_collective"] = r.statistic; row[f"p_partial_{k}"] = r.pvalue
        # error de fase, salvo rotacion y sentido
        best = max(abs(np.mean(np.exp(1j * (ang - T.theta.values)))), abs(np.mean(np.exp(1j * (ang + T.theta.values)))))
        row["phase_concordance_R"] = best
    rows.append(row)
    if keep is not None:
        keep.append(T.assign(inferred_amplitude=amp))
    return row


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    rng = np.random.default_rng(SEED); rows = []; samples = []
    print(f"{N} muestras por escenario; amplitud celular, sincronia y fraccion ritmica varian de forma independiente\n")
    r = summarise("reference (noise 0.5)", *cohort(rng, noise=0.5), rows, samples)
    print("1. La amplitud inferida sigue a la amplitud colectiva")
    print(f"   rho con la amplitud colectiva c*a*R: {r['rho_collective_amplitude']:+.2f}")
    print(f"   rho con cada componente por separado: amplitud celular {r['rho_cell_amplitude']:+.2f}, "
          f"sincronia {r['rho_synchrony']:+.2f}, fraccion ritmica {r['rho_rhythmic_fraction']:+.2f}")
    print("2. No distingue entre sus componentes")
    print("   descontada la amplitud colectiva, lo que queda de cada uno: "
          + ", ".join(f"{k.replace('_', ' ')} {r[f'partial_{k}_given_collective']:+.2f} (P = {r[f'p_partial_{k}']:.2f})"
                      for k in ("cell_amplitude", "synchrony", "rhythmic_fraction")))
    print("   -> dos muestras con la misma amplitud colectiva son indistinguibles para el estimador, se llegue a ella\n"
          "      por celulas mas debiles, por celulas desincronizadas o por menos celulas ritmicas.")
    print("3. Ruido de medida")
    for nz in (0.25, 1.0, 2.0):
        q = summarise(f"noise {nz}", *cohort(rng, noise=nz), rows)
        print(f"   ruido {nz:4.2f}: rho con la amplitud colectiva {q['rho_collective_amplitude']:+.2f}; concordancia de fase R = {q['phase_concordance_R']:.2f}")
    print("4. Un cambio de composicion que no toca el reloj (el segundo tipo celular tiene otro nivel basal)")
    for ds in (0.5, 1.0, 2.0):
        q = summarise(f"composition shift sd {ds}", *cohort(rng, noise=0.5, delta_sd=ds), rows)
        print(f"   desplazamiento basal sd {ds:3.1f}: rho con la amplitud colectiva {q['rho_collective_amplitude']:+.2f}; "
              f"fraccion ritmica, descontada la colectiva {q['partial_rhythmic_fraction_given_collective']:+.2f} "
              f"(P = {q['p_partial_rhythmic_fraction']:.3f})")
    print("5. Sin ritmo")
    q0 = summarise("no rhythm", *cohort(rng, noise=0.5, rhythm=False), rows)
    ref = rows[0]["median_inferred_amplitude"]
    E0, T0 = cohort(rng, noise=0.5, rhythm=False)
    Z0 = ((E0.T - E0.mean(1)) / (E0.std(1) + 1e-9)).T
    ok0, sep0 = ph.is_clock(E0, ph.unsup(Z0.values))
    print(f"   amplitud inferida mediana {q0['median_inferred_amplitude']:.2f}, frente a {ref:.2f} con ritmo: el estimador devuelve un "
          f"numero aunque no haya nada que medir")
    print(f"   el control de antifase BMAL1-DBP {'lo deja pasar' if ok0 else 'lo rechaza'} (separacion {sep0:.1f} h)")
    rows[-1]["antiphase_control_passed"] = bool(ok0)
    # con que frecuencia pasa el control una cohorte sin ritmo? La separacion entre dos fases sin
    # relacion es uniforme entre 0 y 12 h, y el control exige 8 h o mas: un tercio de las veces.
    passed = []
    for k in range(200):
        Ek, _ = cohort(np.random.default_rng(1000 + k), noise=0.5, rhythm=False, n=200)
        Zk = ((Ek.T - Ek.mean(1)) / (Ek.std(1) + 1e-9)).T
        passed.append(ph.is_clock(Ek, ph.unsup(Zk.values))[0])
    rate = float(np.mean(passed))
    print(f"   en 200 cohortes sin ritmo, el control de antifase deja pasar {rate*100:.0f} %: por si solo no basta para "
          f"decir que un tejido tiene reloj")
    rows.append(dict(scenario="no rhythm, 200 cohorts", n=200, antiphase_false_pass_rate=rate))
    pd.DataFrame(rows).to_csv(f"{OUT}/new/factorial_estimator.tsv", sep="\t", index=False)
    samples[0].to_csv(f"{OUT}/new/factorial_estimator_samples.tsv", sep="\t")
