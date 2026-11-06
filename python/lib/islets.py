"""GSE159984: que GSM son islotes de donantes (GPL9115) y cuales preparaciones cultivadas
(GPL16791). Se lee de las dos matrices de serie en lugar de un umbral de GSM escrito a mano."""
import gzip, os

RAW = os.environ.get("T2D_RAW", "data/raw")


def _gsms(fn):
    with gzip.open(fn, "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!Sample_geo_accession"):
                return [x.strip('"') for x in l.rstrip("\n").split("\t")[1:]]
    return []


def groups(raw=None):
    """Devuelve (donantes, cultivos) como listas de GSM."""
    raw = raw or RAW
    return _gsms(f"{raw}/GSE159984-GPL9115_series_matrix.txt.gz"), _gsms(f"{raw}/GSE159984-GPL16791_series_matrix.txt.gz")


def log2cpm_symbols(counts_file, cols, annot, min_count=5, min_samples=10):
    """Cuentas NCBI (GeneID) -> simbolo -> log2 CPM, con el tamano de libreria calculado ANTES
    del filtro de genes para que no dependa del universo de genes."""
    import pandas as pd, numpy as np
    ann = pd.read_csv(annot, sep="\t", usecols=["GeneID", "Symbol"], dtype=str); sym = dict(zip(ann.GeneID, ann.Symbol))
    C = pd.read_csv(counts_file, sep="\t", index_col=0); C = C[[c for c in C.columns if c in set(cols)]]
    C.index = [sym.get(str(i), str(i)) for i in C.index]; C = C[~C.index.str.isdigit()].groupby(level=0).sum()
    lib = C.sum(); C = C[(C > min_count).sum(1) >= min_samples]
    return np.log2(C / lib * 1e6 + 1)
