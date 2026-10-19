"""Anotación de GPL570 por Entrez, no por el símbolo del array.

La anotación de la plataforma es de octubre de 2014: el 8,3 % de sus sondas con un solo Entrez
llevan un símbolo que ya no es el vigente (ARNTL→BMAL1, ATP5A1→ATP5F1A, FAM122C→PABIR3, …).
Casar los conjuntos de genes contra el símbolo del array pierde esos genes en silencio.

Procedimiento: sonda → Entrez (de la anotación de la plataforma) → símbolo actual (anotación de
NCBI GRCh38.p13). Las sondas multi-gen se asignan a cada símbolo que no tenga sonda propia, igual
que en 05_clamp_within_person; las sondas de control y las sin Entrez se descartan.
"""
import os
import numpy as np
import pandas as pd

RAW = os.environ.get("T2D_RAW", "data/raw")


def probe_map(path=None, ncbi=None):
    """Devuelve DataFrame con ID de sonda, Entrez y símbolo actual.

    Vale igual para GPL570 (columna ENTREZ_GENE_ID) y para GPL15207/PrimeView (columna
    "Entrez Gene"); ambas anotaciones de Affymetrix comparten estructura y separador "///".
    """
    g = pd.read_csv(path or f"{RAW}/GPL570_annot.txt", sep="\t", comment="#",
                    low_memory=False, dtype=str)
    col = "ENTREZ_GENE_ID" if "ENTREZ_GENE_ID" in g.columns else "Entrez Gene"
    g = g[["ID", "Gene Symbol", col]].rename(columns={col: "ENTREZ_GENE_ID"})
    ann = pd.read_csv(ncbi or f"{RAW}/Human_GRCh38_p13_annot_tsv.gz", sep="\t",
                      usecols=["GeneID", "Symbol"], dtype=str)
    cur = dict(zip(ann.GeneID, ann.Symbol))
    out = []
    for pid, sym, eid in zip(g.ID, g["Gene Symbol"].fillna(""), g.ENTREZ_GENE_ID.fillna("")):
        es = [x.strip() for x in eid.split("///") if x.strip()]
        new = [cur[e] for e in es if e in cur]
        out.append(dict(ID=pid, old=";".join(x.strip() for x in sym.split("///") if x.strip()),
                        n_genes=len(es), new=";".join(new)))
    return pd.DataFrame(out)


def annotate(expr, path=None, ncbi=None):
    """Colapsa una matriz sonda × muestra a gen × muestra con símbolos actuales.

    Para cada gen toma la sonda de mayor expresión media entre las que lo miden en exclusiva; si
    ninguna lo mide en exclusiva, usa la mejor sonda multi-gen que lo incluya.
    """
    M = probe_map(path, ncbi)
    M = M[(M.new != "") & M.ID.isin(expr.index)]
    mean = expr.mean(1)
    best = {}
    for multi in (False, True):
        sub = M[(M.n_genes > 1) == multi]
        for pid, new in zip(sub.ID, sub.new):
            for gene in new.split(";"):
                if gene in best and not multi:
                    if mean[pid] <= mean[best[gene]]:
                        continue
                elif gene in best:
                    continue
                best[gene] = pid
    genes = sorted(best)
    out = expr.loc[[best[g] for g in genes]].copy()
    out.index = genes
    return out
