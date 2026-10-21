#!/usr/bin/env python3
"""Ensambla signatures_log2_common_genes.tsv, la tabla que leen los paneles H e I de la Figura 3.

Cada columna es el cambio en log2 de un contraste y las filas son los genes medidos en todos, de
modo que las correlaciones entre columnas son comparables. Procedencia de cada columna:

  hydrocortisone  delta por gen de hydrocortisone() en 02_meta_and_validation.py: brazos leidos de
                  la matriz de serie y media tratamiento menos media placebo sobre el fichero
                  procesado del deposito.
  Cushing         mediana Cushing menos mediana control sobre log2 CPM, 5 frente a 11.
  TZD             median_log2fc del contraste tzd_adipose de clamp_genes_change.tsv.
  IR              contraste resistente frente a sensible de GSE13070. Esa comparacion esta
                  desbalanceada (35 frente a 6 en adiposo) y no se pudo regenerar de forma
                  consistente, de modo que se versiona como dato de entrada en
                  data/reference/GSE13070_IR_vs_sensitive_log2.tsv y se declara asi en Metodos.

AVISO de convencion: esta tabla NO usa el centrado por muestra. El centrado corresponde a las
puntuaciones entre individuos; aqui se comparan diferencias entre grupos y centrar cambiaria las
correlaciones. Mezclar las dos convenciones intercambia valores entre columnas.
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
from gene_sets import SETS
OUT = os.environ.get("T2D_OUT", "results"); RAW = os.environ.get("T2D_RAW", "data/raw")
REF = os.environ.get("T2D_REF", "data/reference")

def col_hydrocortisone():
    sys.path.insert(0, os.path.join(HERE))
    src = open(os.path.join(HERE, "02_meta_and_validation.py")).read().split("if __name__")[0]
    g = {"__file__": os.path.join(HERE, "02_meta_and_validation.py"), "__name__": "m2"}
    exec(src, g)
    f = f"{RAW}/GSE148640_processed.txt.gz"
    d = pd.read_csv(f, sep="\t", skiprows=9, low_memory=False)
    d = d.rename(columns={d.columns[1]: "gene"})
    d = d[d.gene.notna() & (d.gene != "---")]
    vals = d.iloc[:, 5:25].apply(pd.to_numeric, errors="coerce")
    pla, trt = g["_hydrocortisone_arms"](raw=RAW, ncol=vals.shape[1])
    d["delta"] = vals.iloc[:, trt].mean(1) - vals.iloc[:, pla].mean(1)
    return d.groupby("gene").delta.mean().rename("hydrocortisone adipose (GSE148640)")

def col_cushing():
    X = pd.read_pickle(f"{OUT}/new/GSE66446_log2cpm.pkl")
    X.index = [{"ARNTL": "BMAL1"}.get(i, i) for i in X.index]
    md = pd.read_pickle(f"{OUT}/new/GSE66446_meta.pkl").set_index("gsm").reindex(X.columns)
    cu = md.c0.str.contains("Cushing").values
    return (X.loc[:, cu].median(1) - X.loc[:, ~cu].median(1)).rename("Cushing adipose (GSE66446)")

def col_tzd():
    W = pd.read_csv(f"{OUT}/clamp_genes_change.tsv", sep="\t")
    w = W[W.contrast == "tzd_adipose"].drop_duplicates("gene").set_index("gene")
    return w.median_log2fc.rename("TZD 3mo adipose (GSE13070)")

def reference():
    """Tabla de referencia versionada: las doce columnas del analisis, incluidas las que no se
    regeneran aqui (contrastes de musculo, manana, dexametasona y miotubos)."""
    return pd.read_csv(f"{REF}/signatures_reference.tsv", sep="\t", index_col=0)

if __name__ == "__main__":
    R = reference()
    # Se regeneran las tres columnas reproducibles y se comprueban contra la referencia; el resto
    # se arrastra de la referencia, que es dato de entrada versionado.
    regen = {c.name: c for c in (col_tzd(), col_hydrocortisone(), col_cushing())}
    for name, v in regen.items():
        com = [g for g in R.index if g in v.index]
        r = np.corrcoef(R.loc[com, name], v[com])[0, 1]
        flag = "ok" if r > 0.95 else "DISCREPANCIA"
        print(f"  {name[:38]:40s} r con la referencia {r:.4f}  {flag}")
        if r <= 0.95:
            raise SystemExit("una columna regenerada no coincide con la referencia: revisar antes de seguir")
        R.loc[com, name] = v[com].values
    S = R.dropna()
    S.to_csv(f"{OUT}/new/signatures_log2_common_genes.tsv", sep="\t")
    clk = [g for g in SETS["Clock output"]["genes"] if g in S.index]
    from scipy import stats
    ir = S.loc[clk, "IR adipose (GSE13070, 35 vs 6)"]
    print(f"{S.shape[0]} genes medidos en los {S.shape[1]} contrastes; {len(clk)} del brazo de salida")
    for c, lab in (("hydrocortisone adipose (GSE148640)", "hydrocortisone"),
                   ("Cushing adipose (GSE66446)", "Cushing"),
                   ("TZD 3mo adipose (GSE13070)", "TZD")):
        r = stats.spearmanr(ir, S.loc[clk, c])
        print(f"  IR vs {lab:15s} rho {r.statistic:+.2f} (P {r.pvalue:.3f})")
    tz = S.loc[clk, "TZD 3mo adipose (GSE13070)"]
    dep = [g for g in clk if ir[g] < 0]
    print(f"  pendiente TZD sobre IR {np.polyfit(ir, tz, 1)[0]:+.3f}; "
          f"{len(dep)}/{len(clk)} genes deprimidos; "
          f"fraccion recuperada mediana {np.median([tz[g]/(-ir[g]) for g in dep])*100:.1f}%")
