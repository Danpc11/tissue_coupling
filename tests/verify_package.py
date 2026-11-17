#!/usr/bin/env python3
"""Verificacion del paquete antes de enviar: coherencia entre codigo, resultados y tablas.

Comprueba lo que fallo en la revision de octubre de 2026 y que no detecta ninguna prueba unitaria:
ficheros de corridas distintas mezclados en el mismo paquete, definiciones de conjuntos que
divergen entre el codigo y las tablas, y cifras citadas que no coinciden con su export.

Uso:  python tests/verify_package.py [directorio_de_resultados]
Devuelve codigo 1 si alguna comprobacion falla.
"""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../python/lib"))
RES = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("T2D_OUT", "results")
FAIL = []

def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: FAIL.append(name)

print("1. Registro de conjuntos de genes")
try:
    import gene_sets as G
    check("el registro se importa", True, f"{len(G.SETS)} conjuntos")
    dup = [k for k, v in G.SETS.items() if len(v["genes"]) != len(set(v["genes"]))]
    check("sin genes repetidos dentro de un conjunto", not dup, str(dup))
    check("cada conjunto declara su fuente", all(len(v["source"]) > 20 for v in G.SETS.values()))
except Exception as e:
    check("el registro se importa", False, str(e)); G = None

print("\n2. Procedencia de los ficheros de concordancia")
for f, coh, plat in [(f"{RES}/crosstissue_concordance_AAGMEx.tsv", "AAGMEx", "Illumina"),
                     (f"{RES}/crosstissue_concordance_GSE41168.tsv", "GSE41168", "Affymetrix")]:
    if not os.path.exists(f):
        check(f"existe {os.path.basename(f)}", False); continue
    d = pd.read_csv(f, sep="\t", nrows=5)
    check(f"{os.path.basename(f)} declara su cohorte", "cohort" in d.columns and d.cohort.iloc[0] == coh,
          d.cohort.iloc[0] if "cohort" in d.columns else "sin columna cohort")
    check(f"{os.path.basename(f)} declara su plataforma", "platform" in d.columns and plat in d.platform.iloc[0],
          d.platform.iloc[0] if "platform" in d.columns else "sin columna platform")
    ident = d.iloc[:, 0].astype(str)
    affy = ident.str.contains("_at").any()
    check(f"{os.path.basename(f)}: identificadores coherentes con la plataforma",
          (affy and plat == "Affymetrix") or (not affy and plat == "Illumina"),
          "sondas _at" if affy else "simbolos")

print("\n3. Conjuntos del codigo frente a los exports")
for f in glob.glob(f"{RES}/ci_*.tsv"):
    d = pd.read_csv(f, sep="\t")
    if G is None or "n_declared" not in d.columns: continue
    bad = [r.set for r in d.itertuples() if r.set in G.SETS and r.n_declared != len(G.SETS[r.set]["genes"])]
    check(f"{os.path.basename(f)}: n declarado coincide con el registro", not bad, str(bad))
    check(f"{os.path.basename(f)}: registra que genes se midieron", "genes_measured" in d.columns)

print("\n4. Aritmetica del metaanalisis")
fci, fme = f"{RES}/new/meta_per_cohort.tsv", f"{RES}/new/meta_pooled.tsv"
if os.path.exists(fci) and os.path.exists(fme):
    ci = pd.read_csv(fci, sep="\t"); me = pd.read_csv(fme, sep="\t")
    me = me[me.grupo == "todos"].set_index("contraste")          # el agrupado de todas las estimaciones
    worst = 0.0
    for s in me.index:
        sub = ci[ci.contraste == s]
        if len(sub) < 2: continue
        z, se = sub.z.values, sub.se.values                         # se ya esta en la escala z de Fisher
        w = 1 / se ** 2; Q = np.sum(w * (z - np.sum(w * z) / np.sum(w)) ** 2); k = len(z)
        C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
        tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
        wr = 1 / (se ** 2 + tau2)
        worst = max(worst, abs(np.tanh(np.sum(wr * z) / np.sum(wr)) - float(me.loc[s, "rho"])))
        check(f"{s}: el numero de estimaciones agrupadas es el de la tabla por cohorte", int(me.loc[s, "k"]) == k, f"k = {k}")
    check("el agrupado se reproduce desde la tabla por cohorte", worst < 1e-9, f"max dif {worst:.2e}")
    check("z es la transformada de Fisher de rho en cada fila", bool(np.allclose(np.arctanh(ci.rho.clip(-.999, .999)), ci.z, atol=1e-9)))
    if "complete" in ci.columns:
        check("estan las trece estimaciones del articulo", bool(ci.complete.all()) and ci[ci.contraste == "NR3C1~reloj"].shape[0] == 13,
              f"{ci[ci.contraste == 'NR3C1~reloj'].shape[0]} estimaciones")
else:
    check("existen los exports del metaanalisis", False,
          "corra el pipeline: pasos 10 y 11 escriben meta_per_cohort.tsv y meta_pooled.tsv")

print("\n5. Benjamini-Hochberg en la expresion diferencial")
# Las tablas gen a gen de casos frente a controles (paso 03) llevan P y FDR. Se recalcula el FDR desde P.
# La de islotes no entra: su valor ajustado es el de PyDESeq2, que filtra genes antes de ajustar.
de = sorted(f for f in glob.glob(f"{RES}/fig1_de_GSE*.tsv"))
check("hay tablas de expresion diferencial que comprobar", len(de) >= 4, f"{len(de)} tablas")
for f in de:
    d = pd.read_csv(f, sep="\t", index_col=0)
    if not {"P", "FDR"} <= set(d.columns):
        check(f"{os.path.basename(f)} lleva P y FDR", False, str(list(d.columns))); continue
    d = d.dropna(subset=["P", "FDR"]).sort_values("P"); m = len(d)
    q = np.minimum.accumulate((d.P.values * m / np.arange(1, m + 1))[::-1])[::-1].clip(max=1)
    check(f"{os.path.basename(f)}: el FDR es el de Benjamini-Hochberg sobre sus P", bool(np.allclose(q, d.FDR.values, atol=1e-9)),
          f"{m} genes, max dif {np.max(np.abs(q - d.FDR.values)):.1e}")

print("\n6. Prediccion sin fuga de informacion")
fp = f"{RES}/new/crosstissue_prediction.tsv"
if os.path.exists(fp):
    d = pd.read_csv(fp, sep="\t")
    check("la prediccion declara n, pliegues y semilla", {"n", "folds", "seed"} <= set(d.columns))
    check("existen predicciones fuera de pliegue auditables",
          len(glob.glob(f"{RES}/new/oof_predictions_*.tsv")) > 0)
    fg, ff = f"{RES}/new/crosstissue_prediction_genes.tsv", f"{RES}/new/crosstissue_prediction_folds.tsv"
    if os.path.exists(fg) and os.path.exists(ff):
        g = pd.read_csv(fg, sep="\t"); F = pd.read_csv(ff, sep="\t")
        pred, outc = set(g[g.role == "predictor"].gene), set(g[g.role == "outcome"].gene)
        check("predictor y resultado no comparten ningun gen", len(pred) > 0 and len(outc) > 0 and not (pred & outc),
              f"{len(pred)} y {len(outc)} genes; en comun: {sorted(pred & outc) or 'ninguno'}")
        if G is not None:
            reg_p, reg_o = set(G.genes("Glucocorticoid")), set(G.panel("Receptor targets, disjoint"))
            check("y son los del registro", pred <= reg_p and outc <= reg_o)
        leak = [int(k) for k, q in F.groupby("fold") if set(q[q.role == "train"].donor) & set(q[q.role == "test"].donor)]
        check("en ningun pliegue un donante entrena y se predice a la vez", not leak, f"pliegues con solape: {leak or 'ninguno'}")
        tests = F[F.role == "test"]
        check("cada donante se predice exactamente una vez", bool((tests.groupby("donor").size() == 1).all()) and tests.donor.nunique() == F.donor.nunique(),
              f"{tests.donor.nunique()} donantes, {F.fold.nunique()} pliegues")
        oof = pd.read_csv(glob.glob(f"{RES}/new/oof_predictions_*.tsv")[0], sep="\t")
        m = oof.merge(tests[["donor", "fold"]], on="donor", suffixes=("", "_declared"))
        check("las predicciones guardadas son las de esos pliegues", len(m) == len(oof) and bool((m.fold == m.fold_declared).all()))
    else:
        check("la prediccion exporta sus genes y sus pliegues", False, "lo escribe el paso 11b")
else:
    check("existe el export de la prediccion", False,
          "corra el pipeline: el paso 11b escribe results/new/crosstissue_prediction.tsv")

print("\n7. Un solo origen: todo lo que hay en resultados lo escribio esta corrida")
import json
fm = f"{RES}/new/run_manifest.json"; fo = f"{RES}/new/output_files.tsv"
if os.path.exists(fm) and os.path.exists(fo):
    M0 = json.load(open(fm)); OF = pd.read_csv(fo, sep="\t")
    check("los ficheros de salida llevan el identificador de la corrida", bool((OF.run_id == M0["run_id"]).all()), f"corrida {M0['run_id']}, {len(OF)} ficheros")
    old = OF[~OF.written_by_this_run & ~OF.declared_cache]
    check("ningun fichero es anterior al inicio de la corrida", not len(old),
          ("anteriores: " + ", ".join(old.file.head(6))) if len(old) else f"{int(OF.declared_cache.sum())} son cache declarada y reutilizada")
    now = {os.path.relpath(os.path.join(d, f), RES) for d, _, fs in os.walk(RES) for f in fs if not f.startswith(".")}
    extra = sorted(now - set(OF.file) - {"new/output_files.tsv", "new/run_manifest.json"})
    check("no hay ficheros que el inventario de la corrida no conozca", not extra, ", ".join(extra[:6]))
else:
    check("existen el manifiesto y el inventario de salidas de la corrida", False, "los escriben los pasos 01a y 29z")

print("\n8. Que corrida es esta")
if os.path.exists(fm):
    M = json.load(open(fm))
    check("la corrida registra su perfil, su version y sus entradas", all(k in M for k in ("profile", "release", "input_files_sha256", "settings")),
          f"perfil {M.get('profile')}, version {M.get('release')}, {M.get('input_files')} ficheros de entrada")
    check("las cuentas de islotes son las de las cifras del articulo (NCBI)", bool(M["islet_counts"]["is_the_reported_source"]), M["islet_counts"]["source"])
    if M["scientific_reproduction"]:
        print("  Es una REPRODUCCION COMPLETA: perfil full con los ajustes por defecto.")
    else:
        print(f"  AVISO: NO es una reproduccion completa (perfil {M['profile']}). Sirve para comprobar que el pipeline funciona; "
              "sus cifras no son las del articulo.")
        for k, v in M["settings_that_differ_from_a_full_reproduction"].items():
            print(f"         {k} = {v['this_run']} (reproduccion completa: {v['full_reproduction']})")
else:
    check("existe el manifiesto de la corrida (run_manifest.json)", False, "lo escribe el paso 01a")

print("\n9. Comparacion de mecanismos")
CM = {k: f"{RES}/new/competing_models_{k}.tsv" for k in ("observed", "simulated", "joint", "verdict", "calibration", "confusion", "seed_stability")}
missing = [k for k, f in CM.items() if not os.path.exists(f)]
check("estan las siete tablas del paso 20c", not missing, "faltan: " + ", ".join(missing) if missing else "")
if not missing:
    J = pd.read_csv(CM["joint"], sep="\t"); CAL = pd.read_csv(CM["calibration"], sep="\t"); CF = pd.read_csv(CM["confusion"], sep="\t")
    S = pd.read_csv(CM["simulated"], sep="\t", usecols=["model", "seed"]); O = pd.read_csv(CM["observed"], sep="\t")
    check("los cuatro mecanismos se evaluan sobre el mismo numero de puntos en cada cohorte",
          bool((J.groupby(["cohort", "scope"]).n_points.nunique() == 1).all()) and J.model.nunique() == 4)
    check("todas las cohortes observadas tienen criterio conjunto y calibracion", set(O.cohort) == set(J.cohort) == set(CAL.cohort))
    check("las cohortes de calibracion tienen el tamano de la observada", bool((CAL.n_calibration_cohort == CAL.n_observed).all()))
    check("el umbral usado no es menor que el nominal ni que ninguno calibrado",
          bool(np.all([J[J.cohort == c].distance_threshold.min() >= CAL[CAL.cohort == c][[x for x in ("nominal_threshold", "calibrated_threshold_max", "calibrated_threshold_fresh_cohorts") if x in CAL.columns]].max().max() - 1e-9
                       for c in CAL.cohort.unique()])))
    check("la matriz de confusion cubre los cuatro mecanismos como generadores y como candidatos",
          CF.generating_model.nunique() == 4 and CF.candidate_model.nunique() == 4 and bool(((CF.accepted_fraction >= 0) & (CF.accepted_fraction <= 1)).all()),
          f"{CF.scenario.nunique() if 'scenario' in CF.columns else 1} escenarios, {int(CF.iloc[0].get('independent_cohorts', CF.iloc[0].get('replicates', 0)))} cohortes por celda")
    check("las semillas de la tabla simulada son las del criterio", int(J.seeds_per_point.max()) == S.seed.nunique(), f"{S.seed.nunique()} semillas")
    check("existe el manifiesto de la cache de simulaciones", os.path.exists(f"{RES}/new/competing_models_cache_manifest.json"))

print("\n" + ("TODO CORRECTO" if not FAIL else f"FALLAN {len(FAIL)}: {FAIL}"))
sys.exit(1 if FAIL else 0)
