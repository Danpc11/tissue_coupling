"""Acceso a CHIRAL (paquete scRITMO) para inferir fase circadiana por muestra.

POR QUE UN METODO EXTERNO. CHIRAL esta publicado y validado (Talamanca et al., y el paquete
scRITMO de Salati et al.): infiere la fase de cada muestra por esperanza-maximizacion sobre un
modelo armonico multivariado, ajustando a la vez las fases de las muestras y los parametros de
los genes. Reimplementarlo fue un error nuestro: la version propia reproducia el metodo pero sin
su validacion ni sus controles.

INSTALACION. El paquete completo arrastra scanpy y pydeseq2, que no hacen falta para CHIRAL. Se
clona el repositorio y se carga el subpaquete pychiral directamente:

    git clone --depth 1 https://github.com/AndreaSalati/scRITMO.git
    pip install numpy pandas scipy statsmodels tqdm anndata

ORIENTACION DE LA MATRIZ. La documentacion de CHIRAL dice que las muestras van en columnas; en
realidad las espera en FILAS (muestras x genes). Pasarla al reves falla con un error de forma que
no lo explica. Esta envoltura fija la orientacion para que no vuelva a ocurrir.
"""
import os
import sys
import types

import numpy as np

SCRITMO_SRC = os.environ.get("SCRITMO_SRC", "/tmp/scRITMO/src")


PYCHIRAL_SRC = os.environ.get("PYCHIRAL_SRC", "")


def _load():
    """Devuelve la funcion CHIRAL. Tres vias, en este orden (rama repro-hugo):
      1. paquete `pychiral` instalado o en PYCHIRAL_SRC (github.com/AndreaSalati/pyCHIRAL, src/)
      2. subpaquete `scritmo.pychiral` de un clon de scRITMO en SCRITMO_SRC
    Las dos copias son la misma implementacion de Andrea Salati; ninguna esta en PyPI."""
    if PYCHIRAL_SRC and os.path.isdir(PYCHIRAL_SRC) and PYCHIRAL_SRC not in sys.path:
        sys.path.insert(0, PYCHIRAL_SRC)
    try:
        from pychiral.chiral import CHIRAL
        return CHIRAL
    except ImportError:
        pass
    if "scritmo" not in sys.modules:
        if not os.path.isdir(SCRITMO_SRC):
            raise ImportError(
                f"No encuentro pychiral (PYCHIRAL_SRC={PYCHIRAL_SRC!r}) ni scRITMO en {SCRITMO_SRC}. "
                "Clona uno de los dos: git clone --depth 1 https://github.com/AndreaSalati/pyCHIRAL.git "
                "(y PYCHIRAL_SRC=<clon>/src) o https://github.com/AndreaSalati/scRITMO.git")
        sys.path.insert(0, SCRITMO_SRC)
        pkg = types.ModuleType("scritmo")
        pkg.__path__ = [os.path.join(SCRITMO_SRC, "scritmo")]
        sys.modules["scritmo"] = pkg
    from scritmo.pychiral.chiral import CHIRAL
    return CHIRAL


def infer_phase(E, genes, standardize=True, **kw):
    """Fase por muestra a partir de una matriz genes x muestras.

    E      : array o DataFrame de genes x muestras (la orientacion habitual en expresion)
    genes  : nombres de los genes, en el orden de las filas de E

    Devuelve un diccionario con phi (fase por muestra, en radianes), sigma (incertidumbre),
    params_g (parametros por gen) y weights (peso de cada gen en el ajuste).

    La fase queda definida salvo rotacion y sentido de giro: lo interpretable son diferencias,
    no valores absolutos.
    """
    CHIRAL = _load()
    X = np.asarray(E, dtype=float)
    if X.shape[0] != len(genes):
        raise ValueError(f"E tiene {X.shape[0]} filas y se dieron {len(genes)} genes; "
                         "E debe ser genes x muestras")
    np.random.seed(0)   # pyCHIRAL usa el generador global sin semilla
    out = CHIRAL(X.T, clockgenes=list(genes), standardize=standardize, **kw)   # CHIRAL quiere muestras x genes; con ndarray ya debe venir restringida a los genes del reloj
    out["phi"] = np.asarray(out["phi"]).ravel()
    return out
