"""Proyeccion de una muestra sobre la plantilla de fase: la forma de onda esperada de cada gen en una
rejilla de fases, y el perfil de concordancia de cada muestra con ella.

Lo usa 02e_gse87761_benchmark.py. El fichero conserva el nombre de un analisis de acoplamiento entre
tejidos dentro del individuo que ya no forma parte del articulo; las funciones de aquel analisis se
retiraron en v2.8.0 porque ningun paso las llamaba.
"""
import numpy as np

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
