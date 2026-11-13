#!/usr/bin/env python3
"""Es especifica de NR3C1 la asociacion con la salida del reloj?

Una correlacion de 0.4 entre un gen y una puntuacion no dice nada por si sola: en estos datos muchos
genes correlacionan con cualquier puntuacion. Aqui se mide NR3C1 contra tres fondos, con la misma
correlacion parcial, las mismas covariables y el mismo metaanalisis que el resultado principal:

  1. TODOS los genes medidos en todas las estimaciones: que percentil ocupa NR3C1 y cuantos genes
     quedan por encima. Los genes de la propia puntuacion se excluyen del fondo.
  2. Los demas RECEPTORES NUCLEARES (panel del registro): si la asociacion fuera propia de esa
     familia, o de los factores de transcripcion en general, NR3C1 no destacaria entre ellos.
  3. Por estimacion: el percentil de NR3C1 en cada cohorte y tejido, para ver que no depende de una.

Entradas: results/new/meta_per_cohort.tsv, meta_all_genes_rho.pkl (las escribe 06_meta_analysis.py)
Salidas:  results/new/receptor_specificity_genes.tsv              rho agrupada de cada gen y su puesto
          results/new/receptor_specificity_nuclear_receptors.tsv  lo mismo para los receptores nucleares
          results/new/receptor_specificity.tsv                    resumen
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
from gene_sets import SETS, panel, SYNONYMS
OUT = os.environ.get("T2D_OUT", "results")
ANCHOR = "NR3C1"


def pooled_random_effects(Z, se):
    """DerSimonian-Laird por filas: Z es genes x estimaciones (z de Fisher), se el error de cada estimacion."""
    w = 1 / se**2; k = Z.shape[1]
    m = (Z * w).sum(1) / w.sum(); Q = (w * (Z - m[:, None])**2).sum(1)
    tau2 = np.clip((Q - (k - 1)) / (w.sum() - (w**2).sum() / w.sum()), 0, None)
    w2 = 1 / (se[None, :]**2 + tau2[:, None])
    return np.tanh((Z * w2).sum(1) / w2.sum(1))


if __name__ == "__main__":
    R = pd.read_csv(f"{OUT}/new/meta_per_cohort.tsv", sep="\t")
    R = R[R.contraste == "NR3C1~reloj"].copy(); R["key"] = R.cohorte + " | " + R.tejido
    A = pd.read_pickle(f"{OUT}/new/meta_all_genes_rho.pkl")
    R = R[R.key.isin(A.columns)]; A = A[list(R.key)]
    incomplete = "complete" in R.columns and not bool(R.complete.all())
    C = A.dropna()                                                   # genes medidos en todas las estimaciones
    score = set(SETS["Clock output"]["genes"]) | {a for g in SETS["Clock output"]["genes"] for a in SYNONYMS.get(g, [])}
    pooled = pd.Series(pooled_random_effects(np.arctanh(C.clip(-.999, .999).values.astype(float)), R.se.values), index=C.index)
    bg = pooled.drop([g for g in pooled.index if g in score or g == ANCHOR])
    above = bg[bg > pooled[ANCHOR]].sort_values(ascending=False)
    pct = float((bg < pooled[ANCHOR]).mean() * 100)
    G = pooled.rename("pooled_rho").to_frame(); G["in_clock_output_score"] = [g in score for g in G.index]
    G["rank"] = G.pooled_rho.rank(ascending=False).astype(int)
    G.sort_values("pooled_rho", ascending=False).to_csv(f"{OUT}/new/receptor_specificity_genes.tsv", sep="\t")
    print(f"{len(R)} estimaciones" + (" (INCOMPLETO: faltan estimaciones del articulo)" if incomplete else "")
          + f", {len(C)} genes medidos en todas")
    print(f"NR3C1: rho agrupada {pooled[ANCHOR]:+.3f}; percentil {pct:.2f} entre los {len(bg)} genes del fondo; "
          f"{len(above)} genes por encima: " + ", ".join(f"{g} {v:+.2f}" for g, v in above.head(10).items()))

    nr = [g for g in panel("Nuclear receptors") if g in pooled.index and g not in score]
    N = pooled[nr].sort_values(ascending=False).rename("pooled_rho").to_frame(); N["rank"] = range(1, len(N) + 1)
    N.to_csv(f"{OUT}/new/receptor_specificity_nuclear_receptors.tsv", sep="\t")
    print(f"receptores nucleares medidos en todas las estimaciones ({len(N)} de {len(panel('Nuclear receptors'))}): NR3C1 es el "
          f"{int(N.loc[ANCHOR, 'rank'])}.o; " + ", ".join(f"{g} {v:+.2f}" for g, v in N.pooled_rho.head(8).items()))
    print("percentil de NR3C1 entre todos los genes, por estimacion: "
          + "; ".join(f"{k} {v:.1f}" for k, v in zip(R.key, R.percentile_among_all_genes)))
    pd.DataFrame([dict(n_estimates=len(R), complete=not incomplete, n_common_genes=len(C), pooled_rho_NR3C1=pooled[ANCHOR],
                       percentile_among_all_genes=pct, genes_above=len(above), genes_above_names="; ".join(above.index[:20]),
                       n_nuclear_receptors=len(N), rank_among_nuclear_receptors=int(N.loc[ANCHOR, "rank"]),
                       top_nuclear_receptor=N.index[0], percentile_per_estimate_median=R.percentile_among_all_genes.median(),
                       percentile_per_estimate_min=R.percentile_among_all_genes.min())]).to_csv(
        f"{OUT}/new/receptor_specificity.tsv", sep="\t", index=False)
