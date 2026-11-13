#!/usr/bin/env python3
"""Ensambla signatures_log2_common_genes.tsv, la tabla que leen los paneles H e I de la Figura 3,
y calcula la semejanza entre la firma de resistencia a la insulina y las de exposicion a
glucocorticoides.

Cada columna es el cambio en log2 de un contraste y las filas son los genes medidos en todos, de
modo que las correlaciones entre columnas son comparables.

DOS TABLAS, CON NOMBRE, porque reproducir lo publicado y recalcular desde los datos publicos son dos
cosas distintas y las dos valen:

  signatures_published_reference.tsv   la copia VERSIONADA (data/reference/signatures_reference.tsv)
        con tres de sus doce columnas sustituidas por su regeneracion, tras comprobar que coinciden
        (r > 0.95): hidrocortisona, Cushing y tiazolidindiona. Las otras nueve, entre ellas la de
        resistencia a la insulina, son las versionadas. Con esta tabla se obtuvieron los valores
        publicados.
  signatures_fully_regenerated.tsv     las doce columnas recalculadas desde los datos publicos por
        el paso 02d.

Las mismas estadisticas se calculan sobre las dos (signature_correlations.tsv), se escribe la
diferencia gen a gen (signature_table_differences.tsv) y se dice si llevan a la misma conclusion.
Las figuras leen la regenerada entera, de modo que todo lo que se dibuja sale de los datos
publicos; T2D_SIGNATURES=published hace que lean la de referencia publicada.

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
    # Tabla A, "published reference": la copia versionada con tres columnas sustituidas por su
    # regeneracion, tras comprobar que coinciden. Es con la que se calcularon los valores publicados.
    regen = {c.name: c for c in (col_tzd(), col_hydrocortisone(), col_cushing())}
    for name, v in regen.items():
        com = [g for g in R.index if g in v.index]
        r = np.corrcoef(R.loc[com, name], v[com])[0, 1]
        flag = "ok" if r > 0.95 else "DISCREPANCIA"
        print(f"  {name[:38]:40s} r con la referencia {r:.4f}  {flag}")
        if r <= 0.95:
            raise SystemExit("una columna regenerada no coincide con la referencia: revisar antes de seguir")
        R.loc[com, name] = v[com].values
    A = R.dropna()
    A.to_csv(f"{OUT}/new/signatures_published_reference.tsv", sep="\t")
    rows = signature_stats(A, "published reference (nine versioned columns, three regenerated)")
    # Tabla B, "fully regenerated": las doce columnas desde los datos publicos (paso 02d).
    full = f"{OUT}/new/signatures_reference.tsv"; B = None
    if os.path.exists(full):
        F = pd.read_csv(full, sep="\t", index_col=0)
        if IR in F.columns:
            B = F.dropna(); B.to_csv(f"{OUT}/new/signatures_fully_regenerated.tsv", sep="\t")
            rows += signature_stats(B, "fully regenerated (twelve columns from the public data)")
    if B is None:
        print("  (aviso) falta la tabla regenerada por el paso 02d: solo existe la de referencia publicada")
    pd.DataFrame(rows).to_csv(f"{OUT}/new/signature_correlations.tsv", sep="\t", index=False)

    # Las dos tablas, gen a gen, en los genes y contrastes que entran en la estadistica
    if B is not None:
        cols = [IR] + [c for c, _ in COMPARE]
        genes = [g for g in SETS["Clock output"]["genes"] if g in A.index or g in B.index]
        D = pd.concat({"published_reference": A.reindex(genes)[cols], "fully_regenerated": B.reindex(genes)[cols]}, axis=1)
        for c in cols:
            D[("difference", c)] = D[("fully_regenerated", c)] - D[("published_reference", c)]
        D.to_csv(f"{OUT}/new/signature_table_differences.tsv", sep="\t")
        T = pd.DataFrame(rows); q = T[T.comparison == "IR vs hydrocortisone"].set_index("table")
        same = bool((np.sign(q.rho).nunique() == 1) and ((q.p_vs_random_sets < 0.05).nunique() == 1))
        print(f"\nLas dos tablas llevan a la misma conclusion sobre la firma de hidrocortisona: {'SI' if same else 'NO'} "
              f"(rho {q.rho.iloc[0]:+.2f} y {q.rho.iloc[1]:+.2f}; frente a conjuntos al azar P = {q.p_vs_random_sets.iloc[0]:.3f} y {q.p_vs_random_sets.iloc[1]:.3f})")

    # Cual leen las figuras. Por defecto la regenerada entera, para que todo lo que se dibuja salga
    # de los datos publicos; T2D_SIGNATURES=published devuelve la de referencia publicada.
    use = os.environ.get("T2D_SIGNATURES", "regenerated")
    S = B if (use == "regenerated" and B is not None) else A
    S.to_csv(f"{OUT}/new/signatures_log2_common_genes.tsv", sep="\t")
    print(f"las figuras 3h y 3i leen la tabla: {'fully regenerated' if S is B else 'published reference'} "
          f"({S.shape[0]} genes, {sum(g in S.index for g in SETS['Clock output']['genes'])} del brazo de salida)")
