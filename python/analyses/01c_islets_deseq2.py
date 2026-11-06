#!/usr/bin/env python3
"""Expresion diferencial en islotes de donantes (GSE159984, GPL9115): T2D frente a no diabeticos.

Produce results/new/GSE159984_islets_T2D_vs_ND.tsv, que leen la Fig. 1b-c (fila de islotes) y
01b_fig1_signature_test.py. Methods: cuentas crudas, genes con mas de 5 cuentas en al menos 10
muestras, PyDESeq2 con las covariables de la cohorte en el diseno y log fold changes encogidos.
La matriz de serie de los donantes solo declara tejido y diagnostico, asi que el diseno es
~ diagnostico. Las 44 preparaciones cultivadas (GPL16791) no entran.

Columnas: las de PyDESeq2 (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) mas t, logFC, P
y FDR con los nombres que usa 01b para las demas cohortes.

Uso:  python python/analyses/01c_islets_deseq2.py
"""
import os, sys, gzip, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

RAW = os.environ.get("T2D_RAW", "data/raw")
OUT = os.environ.get("T2D_OUT", "results")
CULTURE_FROM = 4852119      # GSM a partir del cual son cultivos (criterio de 22_pc_sensitivity_core.py)


def donor_metadata():
    gsm = dx = None
    with gzip.open(f"{RAW}/GSE159984-GPL9115_series_matrix.txt.gz", "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!Sample_geo_accession"): gsm = [x.strip('"') for x in l.rstrip("\n").split("\t")[1:]]
            if l.startswith("!Sample_characteristics_ch1") and "diagnosis" in l:
                dx = [x.strip('"').split(": ", 1)[1] for x in l.rstrip("\n").split("\t")[1:]]
    md = pd.DataFrame({"gsm": gsm, "dx": dx}).set_index("gsm")
    md["group"] = np.where(md.dx.str.contains("type 2"), "T2D", "ND")
    return md


if __name__ == "__main__":
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    ann = pd.read_csv(f"{RAW}/Human_GRCh38_p13_annot_tsv.gz", sep="\t", usecols=["GeneID", "Symbol"], dtype=str)
    sym = dict(zip(ann.GeneID, ann.Symbol))
    C = pd.read_csv(f"{RAW}/GSE159984_raw_counts.tsv.gz", sep="\t", index_col=0)
    C = C[[g for g in C.columns if int(g[3:]) < CULTURE_FROM]]
    C.index = [sym.get(str(i), str(i)) for i in C.index]
    C = C[~C.index.str.isdigit()].groupby(level=0).sum()
    C = C[(C > 5).sum(1) >= 10]
    md = donor_metadata().reindex(C.columns)
    print(f"GSE159984 donantes: {C.shape[1]} muestras ({(md.group=='T2D').sum()} T2D / {(md.group=='ND').sum()} ND), {C.shape[0]} genes")
    dds = DeseqDataSet(counts=C.T.astype(int), metadata=md, design="~group", quiet=True,
                       n_cpus=int(os.environ.get("OMP_NUM_THREADS", "8")))
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "T2D", "ND"], quiet=True)
    st.summary()
    try:
        st.lfc_shrink(coeff="group[T.T2D]")
    except Exception as ex:
        print(f"  (aviso) sin encogimiento de LFC: {ex}")
    r = st.results_df.copy()
    r["t"] = r.stat; r["logFC"] = r.log2FoldChange; r["P"] = r.pvalue; r["FDR"] = r.padj
    r.index.name = "gene"
    os.makedirs(f"{OUT}/new", exist_ok=True)
    r.to_csv(f"{OUT}/new/GSE159984_islets_T2D_vs_ND.tsv", sep="\t")
    print(f"genes con FDR < 0.1: {int((r.padj < 0.1).sum())}  -> {OUT}/new/GSE159984_islets_T2D_vs_ND.tsv")
