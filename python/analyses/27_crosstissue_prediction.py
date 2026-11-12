#!/usr/bin/env python3
"""Prediccion entre tejidos: la senal glucocorticoide de un tejido predice la respuesta en otro
del mismo donante, en GTEx.

Este paso faltaba en run_all.sh: la funcion vivia en 02_meta_and_validation.py, que no es un paso
del pipeline, de modo que la afirmacion del resumen no tenia tabla ni panel que la respaldara.

Sin fuga de informacion, de forma explicita:
  - predictor y resultado usan genes DISJUNTOS (conjunto glucocorticoide en el tejido predictor,
    dianas del receptor en el tejido predicho), de modo que ningun gen aparece a ambos lados;
  - los dos tejidos de un donante van SIEMPRE al mismo pliegue;
  - TODO el aprendizaje (ajuste de covariables tecnicas, estandarizacion y modelo) ocurre dentro
    del pliegue de entrenamiento.
El control "covariates only" usa edad, sexo y escala de Hardy en lugar de la expresion del tejido
predictor: separa la senal del tejido de la demografia y del estado terminal.

Salidas: results/new/crosstissue_prediction.tsv  (una fila por ajuste)
         results/new/oof_predictions_adipose_to_muscle.tsv  (prediccion fuera de pliegue por donante)
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
from sklearn.model_selection import KFold
from sklearn.linear_model import RidgeCV
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
from gene_sets import genes as _g, panel as _panel
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results")

GC = _g("Glucocorticoid")
# dianas del receptor fuera del conjunto predictor: ningun gen a ambos lados
TARGET = _panel("Receptor targets, disjoint")
assert not (set(GC) & set(TARGET)), "predictor y resultado deben ser disjuntos"
TN = {"muscle": "Muscle - Skeletal", "adipose": "Adipose - Subcutaneous"}
AGEMID = {"20-29": 25, "30-39": 35, "40-49": 45, "50-59": 55, "60-69": 65, "70-79": 75}

def run(nfold=5, seed=0):
    need = [f"{OUT}/new/GTEx_adipose_sc.pkl", f"{OUT}/new/GTEx_muscle.pkl", f"{OUT}/new/GTEx_pheno.pkl"]
    if not all(os.path.exists(p) for p in need):
        print("  (omitido) faltan las matrices de GTEx"); return None
    E = {"adipose": pd.read_pickle(need[0]), "muscle": pd.read_pickle(need[1])}
    PH = pd.read_pickle(need[2])
    sa = next((p for p in (f"{RAW}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",
                           f"{RAW}/gtex/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt")
               if os.path.exists(p)), None)
    if sa is None:
        print("  (omitido) falta SampleAttributesDS"); return None
    SA = pd.read_csv(sa, sep="\t", low_memory=False)
    SA["donor"] = SA.SAMPID.str.split("-").str[:2].str.join("-")
    donors = sorted(set(E["adipose"].columns) & set(E["muscle"].columns) & set(PH.index))
    def mat(k, cols):
        g = [x for x in cols if x in E[k].index]
        t = (SA[(SA.SMTSD == TN[k]) & (SA.SMAFRZE == "RNASEQ")]
             .drop_duplicates("donor").set_index("donor").reindex(donors))
        Cv = np.column_stack([np.ones(len(donors)),
                              t.SMRIN.fillna(t.SMRIN.median()).values,
                              t.SMTSISCH.fillna(t.SMTSISCH.median()).values])
        return E[k].loc[g, donors].to_numpy().T, Cv, g
    Yx, Cx, gx = mat("adipose", GC); Yy, Cy, gy = mat("muscle", TARGET)
    p = PH.loc[donors]
    hardy = pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int), drop_first=True).values.astype(float)
    age = p.AGE.map(AGEMID).values.astype(float); sex = (p.SEX.values == 2).astype(float)
    rows = []; oof = None
    for label, extra, covonly in (("raw", None, False), ("death-adjusted", hardy, False),
                                  ("covariates only", hardy, True)):
        kf = KFold(nfold, shuffle=True, random_state=seed)
        pred = np.zeros(len(donors)); obs = np.zeros(len(donors)); fold = np.zeros(len(donors), int)
        for fi, (tr, te) in enumerate(kf.split(donors)):
            def prep(Y, Cv):
                C = Cv if extra is None else np.column_stack([Cv, extra])
                b, *_ = np.linalg.lstsq(C[tr], Y[tr], rcond=None); R = Y - C @ b
                return (R - R[tr].mean(0)) / (R[tr].std(0) + 1e-9)
            y = prep(Yy, Cy).mean(1)
            Z = np.column_stack([age, sex, hardy]) if covonly else prep(Yx, Cx)
            m = RidgeCV(alphas=np.logspace(-2, 3, 20)).fit(Z[tr], y[tr])
            pred[te] = m.predict(Z[te]); obs[te] = y[te]; fold[te] = fi
        r2 = 1 - ((obs - pred) ** 2).sum() / ((obs - obs.mean()) ** 2).sum()
        rho = stats.spearmanr(obs, pred)
        rows.append(dict(predictor="adipose", outcome="muscle", setting=label, n=len(donors),
                         n_pred_genes=len(gx), n_out_genes=len(gy), folds=nfold, seed=seed,
                         r2=r2, rho=rho.statistic, p=rho.pvalue))
        print(f"  adipose -> muscle, {label:16s} n={len(donors)}  R2={r2:+.3f}  rho={rho.statistic:+.3f}")
        if label == "death-adjusted":
            oof = pd.DataFrame(dict(donor=donors, fold=fold, observed=obs, predicted=pred))
    d = pd.DataFrame(rows)
    d.to_csv(f"{OUT}/new/crosstissue_prediction.tsv", sep="\t", index=False)
    oof.to_csv(f"{OUT}/new/oof_predictions_adipose_to_muscle.tsv", sep="\t", index=False)
    return d

if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    print("PREDICCION ENTRE TEJIDOS (genes disjuntos, aprendizaje dentro de pliegue)")
    run()
