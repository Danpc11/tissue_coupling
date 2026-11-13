#!/usr/bin/env python3
"""Proyeccion de una muestra sobre la plantilla de fase: fase, amplitud y calidad de ajuste.

ESTE FICHERO ES UNA LIBRERIA, NO UN PASO. No calcula nada por si solo ni escribe ninguna tabla.
Sus funciones las usan 05_clamp_within_person.py (clamp y tiazolidindiona, dentro de la persona)
y 00d_reference_tables.py (tablas de amplitud). Las respuestas a hidrocortisona, dexametasona y
sindrome de Cushing se calculan en 00d_reference_tables.py y 23_build_signatures_table.py.

QUE SE MIDE Y COMO. Por muestra se estima una fase y una amplitud contra una plantilla fija. La
comparacion entre condiciones se hace sobre DISTRIBUCIONES, no sobre medias:

  - la fase, por la concentracion circular R del grupo, porque la perturbacion puede agrupar las
    fases sin desplazar su centro, y una prueba de medianas no ve eso;
  - la amplitud, por la diferencia ABSOLUTA entre condiciones, nunca por el cociente: con basal
    pequeno el cociente explota y produce ademas correlacion espuria con el valor basal.

AVISO. Esta proyeccion sobre los 800 genes mas ritmicos de la plantilla NO es el estimador
validado dejando uno fuera (ese es el ordenamiento libre sobre el panel de genes del reloj, pasos
04 y 05 de run_all.sh). Sus resultados dependen de la plantilla que se use: con la plantilla
versionada y con la regenerada el cambio de amplitud tras el clamp tiene distinto signo.
"""
import os
import sys
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib"))

OUT = os.environ.get("T2D_OUT", "results")
GRID = np.linspace(0, 2 * np.pi, 721, endpoint=False)


def waveform(template, genes, two_harmonics=True):
    """Perfil esperado de cada gen en cada fase, normalizado a norma unidad por fase.

    La normalizacion por NORMA (no por desviacion de la muestra) es lo que deja la amplitud en
    unidades de expresion: estandarizar la muestra antes de proyectar la vuelve degenerada, porque
    entonces coincide con el coeficiente de ajuste.
    """
    a1 = template.a1[genes].values
    b1 = template.b1[genes].values
    if two_harmonics and "a2" in template.columns:
        a2, b2 = template.a2[genes].values, template.b2[genes].values
        S = np.array([a1 * np.cos(t) + b1 * np.sin(t) + a2 * np.cos(2 * t) + b2 * np.sin(2 * t)
                      for t in GRID])
    else:
        S = np.array([a1 * np.cos(t) + b1 * np.sin(t) for t in GRID])
    return S / np.linalg.norm(S, axis=1, keepdims=True)


def phase_amplitude(E, template, genes):
    """Fase, amplitud y calidad de ajuste por muestra.

    E se centra por gen DENTRO del conjunto que se compara: el centrado define el cero respecto al
    cual se mide la oscilacion, de modo que cambiarlo cambia el resultado. Restringir primero y
    centrar despues es la regla que se aplica en todo el analisis.
    """
    g = [x for x in genes if x in E.index and x in template.index]
    S = waveform(template, g)
    X = E.loc[g].sub(E.loc[g].mean(1), axis=0).values
    proj = S @ X
    idx = np.argmax(np.abs(proj), axis=0)
    fit = np.abs(S @ (X / np.linalg.norm(X, axis=0, keepdims=True)))
    return pd.DataFrame(dict(
        sample=E.columns,
        phase=GRID[idx] * 24 / (2 * np.pi),
        amplitude=[abs(proj[idx[i], i]) for i in range(X.shape[1])],
        fit=[fit[idx[i], i] for i in range(X.shape[1])]))


def concentration(phase_hours):
    """Concentracion circular R: 0 disperso, 1 agrupado. Es la cantidad que cambia cuando una
    perturbacion agrupa las fases sin desplazar su centro."""
    a = np.asarray(phase_hours) * 2 * np.pi / 24
    return float(np.abs(np.mean(np.exp(1j * a))))


def compare(df, condition, lo=5, hi=95):
    """Resumen por condicion, restringido al rango central de ajuste."""
    q = df[(df.fit >= np.percentile(df.fit, lo)) & (df.fit <= np.percentile(df.fit, hi))]
    rows = []
    for c, sub in q.groupby(condition):
        rows.append(dict(condition=c, n=len(sub),
                         R_phase=concentration(sub.phase.values),
                         median_amplitude=float(sub.amplitude.median()),
                         median_fit=float(sub.fit.median())))
    return pd.DataFrame(rows)


def sweep_fit_ranges(df, condition, ranges=((0, 100), (10, 90), (5, 95), (12, 88), (15, 85))):
    """Comprueba que el efecto no depende del rango de ajuste elegido."""
    out = []
    for lo, hi in ranges:
        t = compare(df, condition, lo, hi)
        t.insert(0, "fit_range", f"{lo}-{hi}")
        out.append(t)
    return pd.concat(out, ignore_index=True)


if __name__ == "__main__":
    print("04_perturbations.py es una libreria de funciones; no es un paso del pipeline y no escribe nada.")
