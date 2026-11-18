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


def diagnosis(raw):
    """Serie GSM -> 'T2D' o 'ND' de los donantes (GPL9115), leida de la matriz de serie."""
    import gzip, pandas as pd, numpy as np
    gsm = dx = None
    with gzip.open(f"{raw}/GSE159984-GPL9115_series_matrix.txt.gz", "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!Sample_geo_accession"): gsm = [x.strip('"') for x in l.rstrip("\n").split("\t")[1:]]
            if l.startswith("!Sample_characteristics_ch1") and "diagnosis" in l:
                dx = [x.strip('"').split(": ", 1)[1] for x in l.rstrip("\n").split("\t")[1:]]
    return pd.Series(np.where(pd.Series(dx).str.contains("type 2"), "T2D", "ND"), index=gsm)


def counts_source(raw):
    """De donde salen las cuentas de islotes que va a leer el pipeline. Hay dos rutas y no son
    equivalentes: las de NCBI (HISAT2 y featureCounts sobre RefSeq), que son las de las cifras del
    articulo, y las de ARCHS4 (kallisto sobre GENCODE), que escribe 00c_archs4_counts.py. Se decide por
    el fichero al que apunta GSE159984_raw_counts.tsv.gz y, si el nombre no lo dice, por sus
    identificadores. Devuelve tambien la huella del fichero y de la anotacion."""
    import os, gzip, hashlib
    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for c in iter(lambda: fh.read(1 << 22), b""): h.update(c)
        return h.hexdigest()
    path = f"{raw}/GSE159984_raw_counts.tsv.gz"; ann = f"{raw}/Human_GRCh38_p13_annot_tsv.gz"
    if not os.path.exists(path):
        return dict(source="missing", is_the_reported_source=False, file=None)
    real = os.path.basename(os.path.realpath(path))
    with gzip.open(path, "rt") as fh:
        next(fh); ids = [next(fh).split("\t", 1)[0].strip('"') for _ in range(200)]
    numeric = sum(x.isdigit() for x in ids) / len(ids)
    if "ARCHS4" in real:
        src, quant, ok = "ARCHS4", "kallisto against GENCODE", False
    elif "NCBI" in real and numeric > 0.9:
        src, quant, ok = "NCBI", "HISAT2 and featureCounts against RefSeq (NCBI-generated counts)", True
    elif numeric > 0.9:
        src, quant, ok = "NCBI by its identifiers (GeneID); the file name does not say", "HISAT2 and featureCounts against RefSeq, if it is the NCBI file", True
    else:
        src, quant, ok = "unrecognised: gene symbols as identifiers, so not the NCBI file as downloaded", "unknown", False
    # la copia versionada en el repositorio es la de las cifras del articulo: se dice si es esa misma
    ref = os.path.join(os.environ.get("T2D_REF", "data/reference"), "ncbi_counts", "GSE159984_raw_counts_GRCh38_p13_NCBI.tsv.gz")
    same = (sha(ref) == sha(path)) if os.path.exists(ref) else None
    if same is False and ok:
        src += " (but not the file versioned in data/reference/ncbi_counts)"
    return dict(source=src, quantification=quant, is_the_reported_source=ok, is_the_versioned_file=same, file=real, sha256=sha(path),
                annotation_file=os.path.basename(ann) if os.path.exists(ann) else None,
                annotation_sha256=sha(ann) if os.path.exists(ann) else None)
