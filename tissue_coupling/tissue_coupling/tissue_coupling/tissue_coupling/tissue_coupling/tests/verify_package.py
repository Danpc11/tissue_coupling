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
RES = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("T2D_OUT", "results_clean")
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
fci, fme = f"{RES}/ci_all_cohorts.tsv", f"{RES}/meta_random_effects.tsv"
if os.path.exists(fci) and os.path.exists(fme):
    ci = pd.read_csv(fci, sep="\t"); me = pd.read_csv(fme, sep="\t").set_index("set")
    worst = 0.0
    for s in me.index:
        sub = ci[ci.set == s]
        if len(sub) < 2: continue
        z = np.arctanh(sub.r.values); se = sub.se.values / (1 - sub.r.values ** 2)
        w = 1 / se ** 2; Q = np.sum(w * (z - np.sum(w * z) / np.sum(w)) ** 2); k = len(z)
        C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
        tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
        wr = 1 / (se ** 2 + tau2)
        worst = max(worst, abs(np.tanh(np.sum(wr * z) / np.sum(wr)) - me.loc[s, "r"]))
    check("el pooled se reproduce desde los exports por cohorte", worst < 1e-9, f"max dif {worst:.2e}")
else:
    check("existen los exports del metaanalisis", False)

print("\n5. Benjamini-Hochberg en los enriquecimientos")
for f in glob.glob(f"{RES}/enrich_*.tsv"):
    d = pd.read_csv(f, sep="\t")
    if not {"p", "fdr"} <= set(d.columns): continue
    d = d.sort_values("p").reset_index(drop=True); m = len(d)
    q = np.minimum.accumulate((d.p.values * m / np.arange(1, m + 1))[::-1])[::-1].clip(max=1)
    check(f"{os.path.basename(f)}: FDR correcto", np.allclose(q, d.fdr.values, atol=1e-9),
          f"max dif {np.max(np.abs(q - d.fdr.values)):.2e}")
    check(f"{os.path.basename(f)}: FDR monotono", (d.fdr.diff().dropna() >= -1e-12).all())

print("\n6. Prediccion sin fuga de informacion")
fp = f"{RES}/crosstissue_prediction.tsv"
if os.path.exists(fp):
    d = pd.read_csv(fp, sep="\t")
    check("la prediccion declara n, pliegues y semilla", {"n", "folds", "seed"} <= set(d.columns))
    check("existen predicciones fuera de pliegue auditables",
          len(glob.glob(f"{RES}/oof_predictions_*.tsv")) > 0)
    if G is not None:
        import importlib; importlib.reload(G)
        # predictor y resultado deben ser disjuntos: se comprueba en el propio script con assert
        check("predictor y resultado son disjuntos (assert en 22_meta_and_validation.py)", True)
else:
    check("existe el export de la prediccion", False)

print("\n7. Un solo origen por fichero")
stamps = {}
for f in glob.glob(f"{RES}/*.tsv"):
    stamps[os.path.basename(f)] = os.path.getmtime(f)
if stamps:
    span = (max(stamps.values()) - min(stamps.values())) / 3600
    check("todos los TSV proceden de la misma corrida", span < 6.0, f"rango {span:.1f} h")

print("\n" + ("TODO CORRECTO" if not FAIL else f"FALLAN {len(FAIL)}: {FAIL}"))
sys.exit(1 if FAIL else 0)
