#!/usr/bin/env python3
"""Respuesta de la fase y la amplitud a perturbaciones metabolicas.

Tres perturbaciones, con la misma plantilla y el mismo estimador en las tres:

  GSE107537  retraso de 10 h del periodo de sueno habitual, en laboratorio controlado, sangre.
             Es el control POSITIVO: las horas de muestreo estan registradas y el resultado esta
             publicado por otro grupo, de modo que sirve para comprobar que el estimador asigna la
             perturbacion a la componente correcta.
  GSE13070   clamp de insulina (agudo) y tiazolidindiona (sostenida), musculo y adiposo.

QUE SE MIDE Y COMO. Por muestra se estima una fase y una amplitud contra una plantilla fija. La
comparacion entre condiciones se hace sobre DISTRIBUCIONES, no sobre medias:

  - la fase, por la concentracion circular R del grupo, porque la perturbacion puede agrupar las
    fases sin desplazar su centro, y una prueba de medianas no ve eso;
  - la amplitud, por la diferencia ABSOLUTA entre condiciones, nunca por el cociente: con basal
    pequeno el cociente explota (se observaron cambios de +1000 % que son division por casi cero)
    y produce ademas correlacion espuria con el valor basal.

FILTRO DE AJUSTE. Las muestras se restringen a los percentiles 5-95 del ajuste a la onda. El
criterio se fija por regla y no por resultado, y el efecto debe comprobarse en varios rangos: en
GSE13070 la concentracion de fase tras el clamp aparece igual sin filtro y con cualquiera de los
rangos ensayados.

Uso:  python python/analyses/04_perturbations.py
Salidas: results/perturbation_phase_amplitude.tsv
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
    tpath = os.path.join(OUT, "template_2harm.tsv")
    if not os.path.exists(tpath):
        print(f"Falta {tpath}. Ejecuta antes 02b_phase_template.py con los dos armonicos.")
        sys.exit(1)
    template = pd.read_csv(tpath, sep="\t").set_index("gene")
    genes = list(template.sort_values("r2", ascending=False).index[:800])
    print(f"plantilla: {len(template)} genes, se usan {len(genes)}")
    print("Este paso necesita las matrices por cohorte preparadas en 00_prepare_data.py;")
    print("las funciones de arriba son las que aplica el analisis publicado.")
