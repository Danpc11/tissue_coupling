"""Acoplamiento entre dos tejidos medido dentro de cada individuo, con anclaje de simultaneidad.

EL PROBLEMA. Estimar la fase de cada tejido por separado y restarlas deja una ambiguedad de medio
ciclo: el perfil de verosimilitud tiene un maximo secundario a 12 h, y en una fraccion de las
muestras el maximo cae ahi. En GTEx eso afectaba al 15 % de los donantes, que aparecian con
diferencias de fase proximas a 12 h y acoplamiento negativo. No es desincronia: es lectura
equivocada.

EL ANCLAJE ES FISICO, NO UNA REGULARIZACION. Las dos biopsias proceden del MISMO cuerpo en el
MISMO instante, de modo que ambas lecturas corresponden a la misma hora de reloj. Estimar las dos
fases de forma independiente equivale a fingir que podrian venir de momentos distintos, lo que
contradice el diseno del muestreo. El modelo incorpora esa restriccion buscando la PAREJA de
fases que maximiza la verosimilitud conjunta mas un termino que penaliza diferencias alejadas de
la tipica de la cohorte.

QUE SE ANCLA Y QUE NO. Se ancla la SIMULTANEIDAD, no el acoplamiento. Un individuo cuyo reloj
muscular va genuinamente tres horas por delante del adiposo es admisible y el modelo lo recupera;
lo que se excluye es la lectura de doce horas, que no es desincronia sino error. La diferencia
tipica no se impone: se estima de la cohorte, y en GTEx sale +0.66 h, practicamente cero, como
corresponde a dos tejidos tomados a la vez.

ELECCION DEL PESO. Se fija por el criterio de que desaparezcan las lecturas a 12 h, no por el
valor del acoplamiento resultante. En GTEx con 2000 genes ese criterio da 0.10: los casos cerca de
12 h pasan de 15 % a 0 % y los negativos de 19 % a 3 %. Pesos mayores siguen subiendo el
acoplamiento, pero ahi ya se estaria imponiendo el resultado en lugar de desempatar.

LIMITACION QUE PERMANECE. La correlacion de residuos no es independiente de la diferencia de fase
(rho = -0.57 tras el anclaje, -0.83 sin el). La medida recodifica en parte esa diferencia y no
debe presentarse como informacion separada de ella.
"""
import numpy as np
import pandas as pd
from scipy import stats

GRID = np.linspace(0, 2 * np.pi, 241, endpoint=False)


def waveform(template, genes, two_harmonics=True):
    """Perfil esperado de cada gen en cada fase de la rejilla, estandarizado por fase."""
    a1 = template.a1[genes].values
    b1 = template.b1[genes].values
    if two_harmonics and "a2" in template.columns:
        a2, b2 = template.a2[genes].values, template.b2[genes].values
        S = np.array([a1 * np.cos(t) + b1 * np.sin(t) + a2 * np.cos(2 * t) + b2 * np.sin(2 * t)
                      for t in GRID])
    else:
        S = np.array([a1 * np.cos(t) + b1 * np.sin(t) for t in GRID])
    return (S - S.mean(1, keepdims=True)) / (S.std(1, keepdims=True) + 1e-9)


def _profile(E, template, genes, two_harmonics=True):
    S = waveform(template, genes, two_harmonics)
    Z = E.loc[genes]
    Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.values
    O = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    return S @ O / len(genes), S, O


def joint_phases(E_a, E_b, template, genes, donors=None, weight=0.10, two_harmonics=True):
    """Fases de los dos tejidos estimadas CONJUNTAMENTE bajo la restriccion de simultaneidad.

    weight = 0 recupera la estimacion independiente (cada tejido por su cuenta).
    """
    cols = list(donors) if donors is not None else [c for c in E_a.columns if c in E_b.columns]
    g = [x for x in genes if x in E_a.index and x in E_b.index and x in template.index]
    Pa, Sa, Oa = _profile(E_a[cols], template, g, two_harmonics)
    Pb, Sb, Ob = _profile(E_b[cols], template, g, two_harmonics)
    # diferencia tipica de la cohorte, estimada de los datos con las fases independientes
    d0 = np.angle(np.exp(1j * (GRID[np.argmax(Pa, axis=0)] - GRID[np.argmax(Pb, axis=0)])))
    delta = np.angle(np.mean(np.exp(1j * d0)))
    prior = weight * np.cos(np.angle(np.exp(1j * (GRID[:, None] - GRID[None, :] - delta))))
    ia, ib = [], []
    for i in range(len(cols)):
        J = Pa[:, i][:, None] + Pb[:, i][None, :] + prior
        k = np.unravel_index(np.argmax(J), J.shape)
        ia.append(k[0]); ib.append(k[1])
    return np.array(ia), np.array(ib), (Sa, Oa, Sb, Ob), cols, g, delta


def coupling(E_a, E_b, template, genes, donors=None, weight=0.10, two_harmonics=True, seed=0):
    """Acoplamiento por individuo, con su nulo de personas cruzadas.

    El contraste correcto NO es contra cero: dos personas sin relacion comparten la misma
    plantilla, y eso solo produce un fondo apreciable. El nulo empareja el tejido de una persona
    con el de otra.
    """
    ia, ib, (Sa, Oa, Sb, Ob), cols, g, delta = joint_phases(
        E_a, E_b, template, genes, donors, weight, two_harmonics)
    nrm = lambda M: (M - M.mean(0)) / (M.std(0) + 1e-9)
    Ra = Oa - nrm(Sa[ia].T)
    Rb = Ob - nrm(Sb[ib].T)
    n = len(cols)
    r = np.array([stats.spearmanr(Ra[:, i], Rb[:, i])[0] for i in range(n)])
    rng = np.random.default_rng(seed)
    pr = rng.permutation(n)
    null = np.array([stats.spearmanr(Ra[:, i], Rb[:, pr[i]])[0] for i in range(n)])
    dphase = np.angle(np.exp(1j * (GRID[ia] - GRID[ib]))) * 24 / (2 * np.pi)
    return pd.DataFrame(dict(donor=cols, coupling=r, null=null,
                             phase_a=GRID[ia] * 24 / (2 * np.pi),
                             phase_b=GRID[ib] * 24 / (2 * np.pi),
                             dphase=dphase)), g, delta


def choose_weight(E_a, E_b, template, genes, donors=None, grid=(0, 0.02, 0.05, 0.10, 0.20),
                  two_harmonics=True, tol=0.01):
    """Peso minimo que reduce las lecturas a medio ciclo por debajo de tol.

    El criterio es la desaparicion de las diferencias proximas a 12 h, NO el valor del
    acoplamiento: elegir el peso por el resultado seria imponerlo. La tolerancia evita exigir
    cero exacto, que empujaria a pesos altos por un solo donante y ahi el prior ya no desempata
    sino que impone. Con tol = 0.01 (menos del 1 % de los donantes) GTEx elige 0.10.
    """
    out = []
    for w in grid:
        d, _, _ = coupling(E_a, E_b, template, genes, donors, w, two_harmonics)
        near12 = float(np.mean(np.abs(np.abs(d.dphase) - 12) < 3))
        out.append(dict(weight=w, near_12h=near12, negative=float((d.coupling < 0).mean()),
                        median_coupling=float(d.coupling.median())))
    t = pd.DataFrame(out)
    ok = t[t.near_12h <= tol]
    return (float(ok.weight.min()) if len(ok) else float(t.weight.max())), t


def reliability(E_a, E_b, template, genes, donors=None, weight=0.10):
    """Fiabilidad de la cifra individual: dos mitades disjuntas de genes deben ordenar igual a las
    personas. Decide si la medida sirve como variable de sujeto."""
    g = [x for x in genes if x in E_a.index and x in E_b.index and x in template.index]
    a = coupling(E_a, E_b, template, g[0::2], donors, weight)[0].coupling.values
    b = coupling(E_a, E_b, template, g[1::2], donors, weight)[0].coupling.values
    return float(stats.spearmanr(a, b)[0])
