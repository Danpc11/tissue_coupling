#!/usr/bin/env python3
"""Metaanalisis entre cohortes, reto de glucocorticoides y lectura entre tejidos.

Combina las estimaciones por cohorte del paso 21 con efectos aleatorios, aplica la puntuacion
CONGELADA a una perturbacion glucocorticoide controlada, y mide cuanto puede leerse el eje de un
tejido desde otro sin fuga de informacion.

Tres cautelas implementadas de forma explicita:
  - La puntuacion NO se reajusta sobre los datos del reto: se aplica tal cual se definio.
  - En la prediccion entre tejidos, TODO el aprendizaje (ajuste de covariables, estandarizacion y
    modelo) ocurre dentro de los pliegues de entrenamiento, y los dos tejidos de un donante van
    siempre al mismo pliegue.
  - Predictor y resultado usan genes DISJUNTOS, de modo que ningun gen aparece a ambos lados.

Salidas: results/gtex/meta_random_effects.tsv, hydrocortisone_response.tsv, crosstissue_prediction.tsv
"""
import os, sys, gzip, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/gtex")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib"))
from gene_sets import genes as _g
GC = _g("Glucocorticoid")          # registro unico: python/lib/gene_sets.py
INS = _g("Insulin core")
TARGET = ["MT1X", "MT2A", "CEBPB", "TSC22D1", "SESN3", "KLF9", "DUSP1", "ERRFI1"]
assert not (set(GC) & set(TARGET)), "predictor y resultado deben ser disjuntos"

def dersimonian_laird(r, se):
    z = np.arctanh(r); sez = se / (1 - r ** 2)
    w = 1 / sez ** 2; Q = np.sum(w * (z - np.sum(w * z) / np.sum(w)) ** 2); k = len(z)
    Cc = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / Cc) if Cc > 0 else 0.0
    wr = 1 / (sez ** 2 + tau2); zm = np.sum(wr * z) / np.sum(wr); sem = np.sqrt(1 / np.sum(wr))
    return np.tanh(zm), np.tanh(zm - 1.96 * sem), np.tanh(zm + 1.96 * sem), tau2, \
           (max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0.0), Q, k

def meta(files):
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files if os.path.exists(f)], ignore_index=True)
    rows = []
    for s in d.set.unique():
        sub = d[d.set == s]
        if len(sub) < 2: continue
        r, lo, hi, tau2, I2, Q, k = dersimonian_laird(sub.r.values, sub.se.values)
        rows.append(dict(set=s, r=r, lo=lo, hi=hi, tau2=tau2, I2=I2, Q=Q, k=k))
    return pd.DataFrame(rows), d

def _hydrocortisone_arms(raw=RAW, ncol=20):
    """Qué columna del fichero procesado es placebo y cuál hidrocortisona.

    La matriz de serie declara el tratamiento de cada muestra. Se lee de ahí en lugar de asumir
    que las columnas alternan placebo/tratamiento: si GEO reordena el fichero o alguien lo
    reprocesa, asumir el orden invertiría el signo del resultado sin ningún aviso.
    Si la matriz de serie no está, se vuelve al orden alternante y se avisa.
    """
    f = f"{raw}/GSE148640_series_matrix.txt.gz"
    if not os.path.exists(f):
        print("  (aviso) sin matriz de serie: se asume orden placebo/tratamiento alternante")
        return list(range(0, ncol, 2)), list(range(1, ncol, 2))
    import gzip
    titles = None
    with gzip.open(f, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!Sample_characteristics_ch1") and "treatment:" in line:
                titles = [x.strip('"').split(":", 1)[1].strip()
                          for x in line.rstrip().split("\t")[1:]]
                break
    if titles is None or len(titles) != ncol:
        print("  (aviso) no pude leer el tratamiento: se asume orden alternante")
        return list(range(0, ncol, 2)), list(range(1, ncol, 2))
    pla = [i for i, t in enumerate(titles) if t.lower().startswith("placebo")]
    trt = [i for i, t in enumerate(titles) if not t.lower().startswith("placebo")]
    print(f"  brazos leídos de la matriz de serie: {len(pla)} placebo, {len(trt)} hidrocortisona")
    return pla, trt


def hydrocortisone():
    """GSE148640: hidrocortisona frente a placebo, intrasujeto, 10 pacientes. Puntuacion congelada.

    La puntuación NO se reajusta aquí: los genes vienen definidos de las cohortes pareadas y este
    experimento sólo los aplica. Reajustarla sobre estos datos haría la prueba circular.
    """
    f = f"{RAW}/GSE148640_processed.txt.gz"
    if not os.path.exists(f): print("  (omitido) falta", f); return None
    d = pd.read_csv(f, sep="\t", skiprows=9, low_memory=False)
    d = d.rename(columns={d.columns[1]: "gene", d.columns[2]: "p", d.columns[3]: "FC"})
    d = d[d.gene.notna() & (d.gene != "---")]
    vals = d.iloc[:, 5:25].apply(pd.to_numeric, errors="coerce")
    pla, trt = _hydrocortisone_arms(ncol=vals.shape[1])
    d["delta"] = vals.iloc[:, trt].mean(1) - vals.iloc[:, pla].mean(1)
    d = d.groupby("gene").delta.mean()
    rng = np.random.default_rng(0); rows = []
    for nm, gs in [("Glucocorticoid", GC), ("Insulin core", INS)]:
        g = [x for x in gs if x in d.index]
        obs = d[g].mean()
        null = np.array([d.sample(len(g), random_state=int(rng.integers(1e9))).mean() for _ in range(5000)])
        rows.append(dict(set=nm, n=len(g), delta=obs, p=(np.sum(np.abs(null) >= abs(obs)) + 1) / 5001))
        print(f"  {nm:16s} delta={obs:+.3f}  p={rows[-1]['p']:.4f}")
    return pd.DataFrame(rows)

def crosstissue_prediction(raw=RAW, out=OUT, nfold=5, seed=0):
    """Predice la puntuacion del eje en musculo a partir de la de adiposo, en GTEx.

    Sin fuga de informacion, de forma explicita:
      - TODO el aprendizaje (ajuste de covariables tecnicas, estandarizacion y modelo) ocurre
        dentro de los pliegues de entrenamiento.
      - Los dos tejidos de un donante van SIEMPRE al mismo pliegue.
      - Predictor y resultado usan genes DISJUNTOS (GC en adiposo, TARGET en musculo).
    Exporta ademas las predicciones fuera de pliegue por donante, para que sean auditables.
    """
    from sklearn.model_selection import KFold
    from sklearn.linear_model import RidgeCV
    TN = {"muscle": "Muscle - Skeletal", "adipose": "Adipose - Subcutaneous", "blood": "Whole Blood"}
    E = {k: pd.read_pickle(f"{out}/{k}.pkl") for k in TN if os.path.exists(f"{out}/{k}.pkl")}
    if "muscle" not in E or "adipose" not in E:
        print("  (omitida) faltan las matrices de GTEx"); return None
    ph = pd.read_csv(f"{raw}/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt", sep="\t").set_index("SUBJID")
    sa = pd.read_csv(f"{raw}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt", sep="\t", low_memory=False)
    sa["donor"] = sa.SAMPID.str.split("-").str[:2].str.join("-")
    rows = []
    for src in [k for k in ("adipose", "blood") if k in E]:
        donors = sorted(set(E[src].columns) & set(E["muscle"].columns) & set(ph.index))
        if len(donors) < 60: continue
        def mat(k, cols):
            g = [x for x in cols if x in E[k].index]
            t = sa[(sa.SMTSD == TN[k]) & (sa.SMAFRZE == "RNASEQ")].drop_duplicates("donor").set_index("donor").reindex(donors)
            Cv = np.column_stack([np.ones(len(donors)),
                                  t.SMRIN.fillna(t.SMRIN.median()).values,
                                  t.SMTSISCH.fillna(t.SMTSISCH.median()).values])
            return E[k].loc[g, donors].to_numpy().T, Cv, g
        Yx, Cx, gx = mat(src, GC); Yy, Cy, gy = mat("muscle", TARGET)
        p = ph.loc[donors]
        hardy = pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int), drop_first=True).values.astype(float)
        age = p.AGE.str[:2].astype(float).values; sex = (p.SEX.values == 2).astype(float)
        for label, extra, covars_only in [("raw", None, False), ("death-adjusted", hardy, False),
                                          ("covariates only", hardy, True)]:
            kf = KFold(nfold, shuffle=True, random_state=seed)
            pred = np.zeros(len(donors)); obs = np.zeros(len(donors)); fold = np.zeros(len(donors), int)
            for fi, (tr, te) in enumerate(kf.split(donors)):
                def prep(Y, Cv):
                    C = Cv if extra is None else np.column_stack([Cv, extra])
                    b, *_ = np.linalg.lstsq(C[tr], Y[tr], rcond=None)
                    R = Y - C @ b
                    return (R - R[tr].mean(0)) / (R[tr].std(0) + 1e-9)
                y = prep(Yy, Cy).mean(1)
                Z = np.column_stack([age, sex, hardy]) if covars_only else prep(Yx, Cx)
                m = RidgeCV(alphas=np.logspace(-2, 3, 20)).fit(Z[tr], y[tr])
                pred[te] = m.predict(Z[te]); obs[te] = y[te]; fold[te] = fi
            r2 = 1 - ((obs - pred) ** 2).sum() / ((obs - obs.mean()) ** 2).sum()
            rows.append(dict(predictor=src, outcome="muscle", setting=label, n=len(donors),
                             n_pred_genes=len(gx), n_out_genes=len(gy), folds=nfold, seed=seed, r2=r2))
            print(f"  {src:7s} -> muscle, {label:16s} n={len(donors)}  R2={r2:+.3f}")
            if label == "raw":
                pd.DataFrame(dict(donor=donors, fold=fold, observed=obs, predicted=pred)).to_csv(
                    f"{out}/oof_predictions_{src}_to_muscle.tsv", sep="\t", index=False)
    d = pd.DataFrame(rows); d.to_csv(f"{out}/crosstissue_prediction.tsv", sep="\t", index=False)
    return d

if __name__ == "__main__":
    pooled, per_cohort = meta([f"{OUT}/ci_AAGMEx.tsv", f"{OUT}/ci_GSE41168.tsv"])
    per_cohort.to_csv(f"{OUT}/ci_all_cohorts.tsv", sep="\t", index=False)
    pooled.to_csv(f"{OUT}/meta_random_effects.tsv", sep="\t", index=False)
    print("METAANALISIS de efectos aleatorios")
    for _, r in pooled.iterrows():
        print(f"  {r.set:16s} r={r.r:+.3f} [{r.lo:+.3f}, {r.hi:+.3f}]  I2={r.I2:5.1f}%  Q={r.Q:.2f}")
    print("\nRETO DE HIDROCORTISONA (membresia congelada; escala log2 del deposito,")
    print("no la escala estandarizada del eje en GTEx)")
    h = hydrocortisone()
    if h is not None: h.to_csv(f"{OUT}/hydrocortisone_response.tsv", sep="\t", index=False)
    print("\nPREDICCION ENTRE TEJIDOS (genes disjuntos, aprendizaje dentro de pliegue)")
    crosstissue_prediction()
