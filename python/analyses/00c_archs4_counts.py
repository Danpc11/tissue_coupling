#!/usr/bin/env python3
"""Cuentas de RNA-seq por serie desde el HDF5 de ARCHS4 (human_gene_v2).

Sustituto de las cuentas generadas por NCBI (`*_raw_counts_GRCh38.p13_NCBI.tsv.gz`), que solo se
sirven desde www.ncbi.nlm.nih.gov y no son accesibles desde este servidor. ARCHS4 cuantifica con
kallisto contra GENCODE; NCBI con HISAT2 + featureCounts contra RefSeq. Las cuentas no son
identicas, pero para log2 CPM y scores por gen la diferencia es menor que la de un experimento.

Escribe `data/raw/<acc>_raw_counts.tsv.gz` (indice = simbolo, columnas = GSM), que es lo que leen
09_hpo_by_tissue.py, 22_pc_sensitivity_core.py y fig3_eje.py para GSE159984. Esos scripts
aceptan indice de simbolo porque solo traducen las filas que son numericas (GeneID).

Uso:  python python/analyses/00c_archs4_counts.py GSE159984 [GSE...]
      T2D_ARCHS4=data/raw/archs4/human_gene_v2.latest.h5 (por defecto)
"""
import os, sys
import numpy as np, pandas as pd, h5py

RAW = os.environ.get("T2D_RAW", "data/raw")
H5 = os.environ.get("T2D_ARCHS4", f"{RAW}/archs4/human_gene_v2.latest.h5")


URL = "https://s3.k8s.maayanlab.cloud/archs4/files/human_gene_v2.latest.h5"


def open_h5(h5):
    """Fichero local si existe; si no, lectura por rangos HTTP (fsspec). El dataset esta
    troceado por muestra (chunks 2000 genes x 1 muestra), asi que leer una serie de ~100
    muestras cuesta unos miles de peticiones pequenas, no los 62 GB."""
    if os.path.exists(h5):
        return h5py.File(h5, "r")
    import fsspec
    print(f"  {h5} no existe: leyendo por HTTP desde {URL}")
    return h5py.File(fsspec.open(URL, mode="rb", block_size=1 << 20, cache_type="readahead").open(), "r")


def series_counts(h5, acc):
    with open_h5(h5) as f:
        series = f["meta/samples/series_id"][:].astype(str)
        gsm = f["meta/samples/geo_accession"][:].astype(str)
        sym = f["meta/genes/symbol"][:].astype(str)
        idx = np.where(np.array([acc in s.split(",") for s in series]))[0]
        if len(idx) == 0:
            raise SystemExit(f"{acc}: ninguna muestra en {h5}")
        X = f["data/expression"]
        # el dataset es genes x muestras (v2); lee columnas en bloque ordenado
        idx = np.sort(idx)
        M = np.column_stack([X[:, i] for i in idx]) if X.shape[0] == len(sym) else X[idx, :].T
    C = pd.DataFrame(M, index=sym, columns=gsm[idx])
    return C.groupby(level=0).sum().astype(int)


if __name__ == "__main__":
    for acc in sys.argv[1:] or ["GSE159984"]:
        C = series_counts(H5, acc)
        out = f"{RAW}/{acc}_raw_counts.tsv.gz"
        C.to_csv(out, sep="\t")
        print(f"{acc}: {C.shape[0]} genes x {C.shape[1]} muestras -> {out}")
