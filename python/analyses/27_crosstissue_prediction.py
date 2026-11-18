#!/usr/bin/env python3
"""Prediccion entre tejidos: la senal glucocorticoide de un tejido predice la respuesta en otro
del mismo donante, en GTEx.

Sin fuga de informacion, de forma explicita:
  - predictor y resultado usan genes DISJUNTOS (conjunto glucocorticoide en el tejido predictor,
    dianas del receptor en el tejido predicho), de modo que ningun gen aparece a ambos lados;
  - los dos tejidos de un donante van SIEMPRE al mismo pliegue;
  - TODO el aprendizaje (ajuste de covariables tecnicas, estandarizacion y modelo) ocurre dentro
    del pliegue de entrenamiento.

TRES CONTROLES, porque dos tejidos de una misma persona comparten mucho mas que una hormona:
  1. "covariates only": edad, sexo y escala de Hardy en lugar de la expresion del tejido predictor.
     Separa la senal del tejido de la demografia y del estado terminal.
  2. Predictor al azar: el mismo calculo con conjuntos de genes de adiposo elegidos al azar, del
     mismo tamano que el glucocorticoide. Mide cuanto predice CUALQUIER conjunto de genes, que no
     es cero: hay un componente de donante generico. La senal especifica es el exceso sobre eso.
  3. Resultado al azar: el mismo predictor frente a conjuntos de genes de musculo al azar.

RIN y tiempo isquemico son los de la muestra cuya expresion esta en la matriz (00b guarda cual en
<tejido>_samples.tsv). Se informa ademas la variante con expresion centrada por muestra, que es la
convencion del resto de las puntuaciones entre individuos; las filas sin centrar son el calculo
con el que se obtuvo el valor publicado.

Salidas: results/new/crosstissue_prediction.tsv         una fila por ajuste
         results/new/crosstissue_prediction_null.tsv    R2 de cada conjunto al azar
         results/new/oof_predictions_adipose_to_muscle.tsv  prediccion fuera de pliegue por donante
         results/new/oof_predictions_all_settings.tsv   lo mismo para las seis configuraciones
         results/new/crosstissue_prediction_design.tsv  el diseno de las seis: predictores, resultado y covariables
         results/new/crosstissue_prediction_genes.tsv   los genes de cada lado
         results/new/crosstissue_prediction_folds.tsv   en cada pliegue, que donantes entrenan y cuales se predicen
"""
import os, sys, warnings; warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
import numpy as np, pandas as pd
from scipy import stats
from sklearn.model_selection import KFold
from sklearn.linear_model import RidgeCV
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
from gene_sets import genes as _g, panel as _panel
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results")
NNULL = int(os.environ.get("T2D_NULL", 300))

GC = _g("Glucocorticoid")
# dianas del receptor fuera del conjunto predictor: ningun gen a ambos lados
TARGET = _panel("Receptor targets, disjoint")
assert not (set(GC) & set(TARGET)), "predictor y resultado deben ser disjuntos"
TN = {"muscle": ("GTEx_muscle", "Muscle - Skeletal"), "adipose": ("GTEx_adipose_sc", "Adipose - Subcutaneous")}
AGEMID = {"20-29": 25, "30-39": 35, "40-49": 45, "50-59": 55, "60-69": 65, "70-79": 75}


def load():
    need = [f"{OUT}/new/GTEx_adipose_sc.pkl", f"{OUT}/new/GTEx_muscle.pkl", f"{OUT}/new/GTEx_pheno.pkl"]
    if not all(os.path.exists(p) for p in need):
        print("  (omitido) faltan las matrices de GTEx"); return None
    sa = next((p for p in (f"{RAW}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",
                           f"{RAW}/gtex/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt")
               if os.path.exists(p)), None)
    if sa is None:
        print("  (omitido) falta SampleAttributesDS"); return None
    E = {"adipose": pd.read_pickle(need[0]), "muscle": pd.read_pickle(need[1])}
    PH = pd.read_pickle(need[2])
    SA = pd.read_csv(sa, sep="\t", low_memory=False)
    SA["donor"] = SA.SAMPID.str.split("-").str[:2].str.join("-")
    donors = sorted(set(E["adipose"].columns) & set(E["muscle"].columns) & set(PH.index))
    tech = {}
    for k, (fn, smtsd) in TN.items():
        sm = f"{OUT}/new/{fn}_samples.tsv"
        if os.path.exists(sm):                       # la muestra cuya expresion esta en la matriz
            used = pd.read_csv(sm, sep="\t", index_col=0).iloc[:, 0].reindex(donors)
            t = SA.set_index("SAMPID").reindex(used.values)
        else:                                        # sin ese fichero, la primera que lista la anotacion
            t = (SA[(SA.SMTSD == smtsd) & (SA.SMAFRZE == "RNASEQ")]
                 .drop_duplicates("donor").set_index("donor").reindex(donors))
        tech[k] = np.column_stack([np.ones(len(donors)), t.SMRIN.fillna(t.SMRIN.median()).values,
                                   t.SMTSISCH.fillna(t.SMTSISCH.median()).values])
    p = PH.loc[donors]
    hardy = pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int), drop_first=True).values.astype(float)
    age = p.AGE.map(AGEMID).values.astype(float); sex = (p.SEX.values == 2).astype(float)
    return {k: E[k][donors] for k in E}, donors, tech, hardy, np.column_stack([age, sex, hardy])


def oof(Yx, Yy, Cx, Cy, extra, Z_cov, nfold, seed):
    """Prediccion fuera de pliegue. Con Z_cov el predictor son las covariables y no la expresion.

    Dos precisiones sobre la escala:
      - el resultado se estandariza en cada pliegue con la media y la desviacion de SU entrenamiento
        (no hay otra forma de hacerlo sin fuga), de modo que cada pliegue queda en su propia escala.
        Por eso, ademas del R2 con todas las predicciones juntas, se devuelve el R2 de cada pliegue,
        que no depende de esa escala;
      - los predictores, sean genes o covariables, se estandarizan igual, con el entrenamiento, antes
        de la regresion con penalizacion: asi la penalizacion pesa lo mismo en los dos modelos que
        se comparan.
      - para tener ademas un R2 con todas las muestras en UNA escala, el mismo modelo se ajusta al
        resultado sin dividir por la desviacion del pliegue: la media, entre los genes diana, de la
        expresion residual en log2. Esa es la cifra de referencia (`r2_common_scale`).
    Devuelve (R2 conjunto, observado, predicho, pliegue, R2 por pliegue, R2 en escala comun)."""
    n = Yy.shape[0]; pred = np.zeros(n); obs = np.zeros(n); fold = np.zeros(n, int); r2f = []
    pred_c = np.zeros(n); obs_c = np.zeros(n)
    for fi, (tr, te) in enumerate(KFold(nfold, shuffle=True, random_state=seed).split(np.arange(n))):
        def prep(Y, Cv, scale=True):
            C = Cv if extra is None else np.column_stack([Cv, extra])
            b, *_ = np.linalg.lstsq(C[tr], Y[tr], rcond=None); R = Y - C @ b
            return (R - R[tr].mean(0)) / (R[tr].std(0) + 1e-9) if scale else R - R[tr].mean(0)
        y = prep(Yy, Cy).mean(1); yc = prep(Yy, Cy, scale=False).mean(1)
        if Z_cov is not None:
            keep = Z_cov[tr].std(0) > 0                      # sin columnas constantes en el entrenamiento
            Z = (Z_cov[:, keep] - Z_cov[tr][:, keep].mean(0)) / Z_cov[tr][:, keep].std(0)
        else:
            Z = prep(Yx, Cx)
        m = RidgeCV(alphas=np.logspace(-2, 3, 20)).fit(Z[tr], y[tr])
        pred[te] = m.predict(Z[te]); obs[te] = y[te]; fold[te] = fi
        pred_c[te] = RidgeCV(alphas=np.logspace(-2, 3, 20)).fit(Z[tr], yc[tr]).predict(Z[te]); obs_c[te] = yc[te]
        r2f.append(1 - ((obs[te] - pred[te]) ** 2).sum() / ((obs[te] - obs[te].mean()) ** 2).sum())
    r2 = 1 - ((obs - pred) ** 2).sum() / ((obs - obs.mean()) ** 2).sum()
    r2c = 1 - ((obs_c - pred_c) ** 2).sum() / ((obs_c - obs_c.mean()) ** 2).sum()
    return r2, obs, pred, fold, np.array(r2f), r2c


def run(nfold=5, seed=0):
    L = load()
    if L is None: return None
    X, donors, tech, hardy, cov = L
    rows = []; out_oof = None; observed_r2 = None; all_oof = []; XY = {}
    for centring in ("none", "per sample"):
        E = X if centring == "none" else {k: v.sub(v.mean(0), axis=1) for k, v in X.items()}
        gx = [g for g in GC if g in E["adipose"].index]; gy = [g for g in TARGET if g in E["muscle"].index]
        Yx = E["adipose"].loc[gx].to_numpy().T; Yy = E["muscle"].loc[gy].to_numpy().T
        for label, extra, zc in (("raw", None, None), ("death-adjusted", hardy, None), ("covariates only", hardy, cov)):
            r2, obs, pred, fold, r2f, r2c = oof(Yx, Yy, tech["adipose"], tech["muscle"], extra, zc, nfold, seed)
            rho = stats.spearmanr(obs, pred)
            rows.append(dict(predictor="adipose", outcome="muscle", setting=label, centring=centring, n=len(donors),
                             n_pred_genes=len(gx), n_out_genes=len(gy), folds=nfold, seed=seed,
                             r2=r2, rho=rho.statistic, p=rho.pvalue,
                             r2_fold_mean=r2f.mean(), r2_fold_min=r2f.min(), r2_fold_max=r2f.max(), r2_common_scale=r2c))
            print(f"  adipose -> muscle, {label:16s} centrado: {centring:10s} n={len(donors)}  R2={r2:+.3f}  rho={rho.statistic:+.3f}"
                  f"   por pliegue: media {r2f.mean():+.3f}, de {r2f.min():+.3f} a {r2f.max():+.3f}   en escala comun: {r2c:+.3f}")
            all_oof.append(pd.DataFrame(dict(setting=label, centring=centring, donor=donors, fold=fold, observed=obs, predicted=pred)))
            XY[centring] = (Yx, Yy, gx, gy)
            if label == "death-adjusted" and centring == "none":
                out_oof = pd.DataFrame(dict(donor=donors, fold=fold, observed=obs, predicted=pred))
                observed_r2 = r2c                 # los nulos se comparan en la escala comun

    # especificidad: el mismo calculo (ajustado por causa de muerte, sin centrar) con genes al azar
    E = X; rng = np.random.default_rng(seed)
    gx = [g for g in GC if g in E["adipose"].index]; gy = [g for g in TARGET if g in E["muscle"].index]
    Yx = E["adipose"].loc[gx].to_numpy().T; Yy = E["muscle"].loc[gy].to_numpy().T
    excl = set(GC) | set(TARGET)
    pool_a = [g for g in E["adipose"].index if g not in excl]; pool_m = [g for g in E["muscle"].index if g not in excl]
    null = []
    for kind in ("random predictor genes", "random outcome genes"):
        v = []
        for _ in range(NNULL):
            if kind == "random predictor genes":
                Ya = E["adipose"].loc[list(rng.choice(pool_a, len(gx), replace=False))].to_numpy().T
                v.append(oof(Ya, Yy, tech["adipose"], tech["muscle"], hardy, None, nfold, seed)[5])
            else:
                Yb = E["muscle"].loc[list(rng.choice(pool_m, len(gy), replace=False))].to_numpy().T
                v.append(oof(Yx, Yb, tech["adipose"], tech["muscle"], hardy, None, nfold, seed)[5])
        v = np.array(v); pval = (np.sum(v >= observed_r2) + 1) / (NNULL + 1)
        null += [dict(null=kind, draw=i, r2=x) for i, x in enumerate(v)]
        rows.append(dict(predictor="adipose", outcome="muscle", setting=f"null: {kind}", centring="none", n=len(donors),
                         n_pred_genes=len(gx), n_out_genes=len(gy), folds=nfold, seed=seed,
                         r2=np.median(v), r2_p95=np.percentile(v, 95), r2_max=v.max(), n_null=NNULL, p_vs_observed=pval))
        print(f"  nulo, {kind:22s} ({NNULL} conjuntos): R2 mediana {np.median(v):+.3f}, percentil 95 {np.percentile(v, 95):+.3f}, "
              f"maximo {v.max():+.3f}; observado {observed_r2:+.3f}, P = {pval:.4f}")
    d = pd.DataFrame(rows)
    d.to_csv(f"{OUT}/new/crosstissue_prediction.tsv", sep="\t", index=False)
    pd.DataFrame(null).to_csv(f"{OUT}/new/crosstissue_prediction_null.tsv", sep="\t", index=False)
    out_oof.to_csv(f"{OUT}/new/oof_predictions_adipose_to_muscle.tsv", sep="\t", index=False)
    # Las seis configuraciones y el diseno con el que se ajustaron, para poder recalcular cada una desde fuera:
    # predictores y resultado sin centrar (x:, y:) y centrados por muestra (xc:, yc:), covariables tecnicas de
    # cada tejido (cx:, cy:), circunstancia de la muerte (h:) y las covariables del modelo de comparacion (z:).
    pd.concat(all_oof, ignore_index=True).to_csv(f"{OUT}/new/oof_predictions_all_settings.tsv", sep="\t", index=False)
    block = lambda M, pre, names=None: pd.DataFrame(M, columns=[f"{pre}:{n}" for n in (names if names is not None else range(M.shape[1]))])
    design = pd.concat([pd.DataFrame({"donor": donors}),
                        block(XY["none"][0], "x", XY["none"][2]), block(XY["none"][1], "y", XY["none"][3]),
                        block(XY["per sample"][0], "xc", XY["per sample"][2]), block(XY["per sample"][1], "yc", XY["per sample"][3]),
                        block(tech["adipose"], "cx"), block(tech["muscle"], "cy"), block(hardy, "h"), block(cov, "z")], axis=1)
    design.to_csv(f"{OUT}/new/crosstissue_prediction_design.tsv", sep="\t", index=False)
    # Lo que hace falta para comprobar desde fuera que no hay fuga: que genes hay a cada lado y, en
    # cada pliegue, que donantes entrenan y cuales se predicen.
    pd.DataFrame([dict(role="predictor", tissue="adipose", gene=g) for g in gx] + [dict(role="outcome", tissue="muscle", gene=g) for g in gy]
                 ).to_csv(f"{OUT}/new/crosstissue_prediction_genes.tsv", sep="\t", index=False)
    parts = []
    for fi, (tr, te) in enumerate(KFold(nfold, shuffle=True, random_state=seed).split(np.arange(len(donors)))):
        parts += [dict(fold=fi, donor=donors[i], role="train") for i in tr] + [dict(fold=fi, donor=donors[i], role="test") for i in te]
    pd.DataFrame(parts).to_csv(f"{OUT}/new/crosstissue_prediction_folds.tsv", sep="\t", index=False)
    return d


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    print("PREDICCION ENTRE TEJIDOS (genes disjuntos, aprendizaje dentro de pliegue)")
    run()
