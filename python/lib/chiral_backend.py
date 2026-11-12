"""Acceso a CHIRAL (paquete scRITMO) para inferir fase circadiana por muestra.

POR QUE UN METODO EXTERNO. CHIRAL esta publicado y validado (Talamanca et al., y el paquete
scRITMO de Salati et al.): infiere la fase de cada muestra por esperanza-maximizacion sobre un
modelo armonico multivariado, ajustando a la vez las fases de las muestras y los parametros de
los genes. Reimplementarlo fue un error nuestro: la version propia reproducia el metodo pero sin
su validacion ni sus controles.

INSTALACION. CHIRAL no esta en PyPI. Hay un unico camino documentado: `bash download_public.sh`
clona el repositorio del metodo, github.com/AndreaSalati/pyCHIRAL, en data/raw/pychiral, en la
revision fijada en ese script, y este modulo lo busca ahi. Necesita ademas tqdm y anndata, que
estan en requirements.txt como opcionales.

Dos variables de entorno permiten usar otra copia: PYCHIRAL_SRC (directorio src de un clon de
pyCHIRAL) y SCRITMO_SRC (directorio src de un clon de scRITMO, que incluye la misma
implementacion como subpaquete). Las dos difieren solo en el nombre de un argumento y de algunas
claves del resultado; por eso los genes se pasan por posicion.

ORIENTACION DE LA MATRIZ. La documentacion de CHIRAL dice que las muestras van en columnas; en
realidad las espera en FILAS (muestras x genes). Pasarla al reves falla con un error de forma que
no lo explica. Esta envoltura fija la orientacion para que no vuelva a ocurrir.
"""
import os
import sys
import types

import numpy as np

RAW = os.environ.get("T2D_RAW", "data/raw")
PYCHIRAL_SRC = os.environ.get("PYCHIRAL_SRC") or os.path.join(RAW, "pychiral", "src")
SCRITMO_SRC = os.environ.get("SCRITMO_SRC") or os.path.join(RAW, "scritmo", "src")


def _load():
    """Devuelve la funcion CHIRAL. Dos vias, en este orden:
      1. paquete `pychiral` instalado o en PYCHIRAL_SRC (por defecto data/raw/pychiral/src)
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
                f"no encuentro pyCHIRAL en {PYCHIRAL_SRC} ni scRITMO en {SCRITMO_SRC}; "
                "`bash download_public.sh` clona pyCHIRAL en data/raw/pychiral")
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
    # genes por posicion: el argumento se llama `clockgenes` en pyCHIRAL y `genes` en scRITMO
    out = CHIRAL(X.T, list(genes), standardize=standardize, **kw)   # CHIRAL quiere muestras x genes
    out["phi"] = np.asarray(out["phi"]).ravel()
    return out


def implementation():
    """De donde sale CHIRAL en esta ejecucion y en que revision, para dejarlo escrito junto al
    resultado. Devuelve None si no esta disponible."""
    import subprocess
    try:
        _load()
    except ImportError:
        return None
    mod = sys.modules.get("pychiral") or sys.modules.get("scritmo.pychiral") or sys.modules.get("scritmo")
    name = "pyCHIRAL" if "pychiral" in sys.modules else "scRITMO"
    src = PYCHIRAL_SRC if name == "pyCHIRAL" else SCRITMO_SRC
    try:
        rev = subprocess.run(["git", "-C", src, "rev-parse", "--short", "HEAD"], capture_output=True,
                             text=True, timeout=10).stdout.strip()
    except Exception:
        rev = ""
    return f"{name} {rev}".strip()
