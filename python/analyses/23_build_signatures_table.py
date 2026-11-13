#!/usr/bin/env python3
"""Ensambla signatures_log2_common_genes.tsv, la tabla que leen los paneles H e I de la Figura 3,
y calcula la semejanza entre la firma de resistencia a la insulina y las de exposicion a
glucocorticoides.

Cada columna es el cambio en log2 de un contraste y las filas son los genes medidos en todos, de
modo que las correlaciones entre columnas son comparables.

DE DONDE SALE CADA COLUMNA, sin ambiguedad. La tabla que usan las figuras parte de la copia
VERSIONADA (data/reference/signatures_reference.tsv) y sustituye tres columnas por su regeneracion
de esta corrida, comprobando antes que coinciden (r > 0.95): hidrocortisona, Cushing y
tiazolidindiona. Las otras nueve, entre ellas la de resistencia a la insulina, son las de la copia
versionada. El paso 02d regenera las doce desde los datos publicos y mide cuanto coincide cada una
con la versionada (reference_tables_check.tsv); la de resistencia a la insulina coincide exactamente.

Para que el resultado no dependa de esa mezcla, las mismas estadisticas se calculan ademas sobre
la tabla regenerada entera (results/new/signatures_reference.tsv) cuando existe, y las dos van a
signature_correlations.tsv.

LA CORRELACION ENTRE FIRMAS SE CALCULA SOBRE OCHO GENES. Tres comprobaciones la acompanan:
  - frente a conjuntos de ocho genes al azar de la misma tabla: dos perfiles de expresion pueden
    correlacionar en cualquier conjunto de genes, y eso es el fondo contra el que hay que medir;
  - la correlacion de los dos perfiles en todos los genes de la tabla;
  - el intervalo al quitar cada gen, porque con ocho puntos uno solo puede sostener el resultado.
El P de Spearman sobre ocho genes trata los genes como independientes y no lo son.

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

IR = "IR adipose (GSE13070, 35 vs 6)"
COMPARE = (("hydrocortisone adipose (GSE148640)", "hydrocortisone"), ("Cushing adipose (GSE66446)", "Cushing"),
           ("TZD 3mo adipose (GSE13070)", "TZD"))

def signature_stats(S, table, nnull=5000, seed=0):
    """Correlacion de la firma de resistencia con cada perfil, sobre los genes de salida del reloj,
    con su fondo (conjuntos al azar y todos los genes) y su estabilidad al quitar un gen."""
    from scipy import stats
    clk = [g for g in SETS["Clock output"]["genes"] if g in S.index]
    ir = S.loc[clk, IR]; rng = np.random.default_rng(seed); rows = []
    print(f"[{table}] {S.shape[0]} genes medidos en los {S.shape[1]} contrastes; {len(clk)} del brazo de salida: {', '.join(clk)}")
    for c, lab in COMPARE:
        if c not in S.columns: continue
        r = stats.spearmanr(ir, S.loc[clk, c])
        A, B = S[IR].values, S[c].values
        null = np.array([stats.spearmanr(A[i], B[i]).statistic
                         for i in (rng.choice(len(S), len(clk), replace=False) for _ in range(nnull))])
        p_null = (np.sum(null >= r.statistic) + 1) / (nnull + 1)
        allg = stats.spearmanr(A, B).statistic
        loo = [stats.spearmanr(ir.drop(g), S.loc[[x for x in clk if x != g], c]).statistic for g in clk]
        rows.append(dict(table=table, comparison=f"IR vs {lab}", n_genes=len(clk), rho=r.statistic, p_spearman=r.pvalue,
                         rho_all_genes=allg, null_median=np.median(null), null_p95=np.percentile(null, 95), p_vs_random_sets=p_null,
                         rho_leave_one_gene_min=min(loo), rho_leave_one_gene_max=max(loo),
                         gene_whose_removal_lowers_most=clk[int(np.argmin(loo))]))
        print(f"  IR vs {lab:15s} rho {r.statistic:+.2f} (P {r.pvalue:.3f}) | todos los genes {allg:+.2f} | "
              f"conjuntos al azar: mediana {np.median(null):+.2f}, percentil 95 {np.percentile(null, 95):+.2f}, P = {p_null:.4f} | "
              f"quitando un gen: {min(loo):+.2f} a {max(loo):+.2f} (el que mas baja: {clk[int(np.argmin(loo))]})")
    tzd = "TZD 3mo adipose (GSE13070)"
    if tzd in S.columns:
        tz = S.loc[clk, tzd]; dep = [g for g in clk if ir[g] < 0]
        slope = np.polyfit(ir, tz, 1)[0]; rec = np.median([tz[g]/(-ir[g]) for g in dep])*100 if dep else np.nan
        rows.append(dict(table=table, comparison="TZD slope on IR", n_genes=len(clk), rho=slope, genes_depressed=len(dep), recovered_fraction_pct=rec))
        print(f"  pendiente TZD sobre IR {slope:+.3f}; {len(dep)}/{len(clk)} genes deprimidos; fraccion recuperada mediana {rec:.1f}%")
    return rows


if __name__ == "__main__":
    R = reference()
    # Se regeneran tres columnas y se comprueban contra la referencia; el resto se toma de la
    # referencia versionada (ver la cabecera).
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
    rows = signature_stats(S, "versioned reference, three columns regenerated")
    full = f"{OUT}/new/signatures_reference.tsv"
    if os.path.exists(full):
        F = pd.read_csv(full, sep="\t", index_col=0)
        if IR in F.columns:
            rows += signature_stats(F.dropna(), "all twelve columns regenerated")
        else:
            print(f"  (aviso) {full} no lleva la columna de resistencia a la insulina; no se calcula la variante regenerada")
    else:
        print("  (aviso) falta la tabla regenerada por el paso 02d: solo se calcula sobre la referencia versionada")
    pd.DataFrame(rows).to_csv(f"{OUT}/new/signature_correlations.tsv", sep="\t", index=False)
