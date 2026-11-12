#!/usr/bin/env python3
"""Prueba de integracion del paso 14b: lo ejecuta de verdad sobre una serie pequena con respuesta
conocida y comprueba el contenido de lo que escribe, no solo que termine.

La serie simulada tiene el mismo formato que GSE40589 en GEO y lleva dentro lo que el paso debe
encontrar:
  - el brazo de salida del reloj cae entre el minuto 0 y el 240 en los tres brazos (hora);
  - el conjunto de respuesta a insulina cae solo con los dos clamps (insulina);
  - el primer par de salina esta formado, a proposito, por biopsias de dos personas distintas;
  - hay arreglos sin pareja y brazos GIP, como en la serie real.

No necesita datos ni red. Complementa a verify_inputs.py, que es un analisis estatico del codigo y
no puede saber si un paso produce el contenido correcto.
"""
import gzip, os, subprocess, sys, tempfile
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python/lib"))
from gene_sets import SETS

fails = []
def check(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FALLA'} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok: fails.append(name)


def write_series(raw):
    rng = np.random.default_rng(5)
    genes = sorted({g for s in SETS.values() for g in s["genes"]}) + [f"G{i}" for i in range(1200)]
    # la plataforma es de 2008: dos genes van con su simbolo de entonces
    genes = ["C1orf51" if g == "CIART" else g for g in genes]
    probes = [f"A_23_P{100000 + i}" for i in range(len(genes))]
    person = {k: rng.normal(0, .6, len(genes)) for k in range(40)}
    titles = ["SCAT_0_NACL_EU_2", "SCAT_240_GIP_EU_3", "SCAT_240_NACL_EU_4", "SCAT_0_NACL_INF_9"]
    who = [30, 31, 32, 33]
    n, k = 13, 0
    for rep in range(7):
        for arm in ("NACL_EU", "GIP_EU", "NACL_INF", "GIP_INF", "NACL_HC", "GIP_HC"):
            titles += [f"SCAT_0_{arm}_{n}", f"SCAT_240_{arm}_{n + 1}"]
            w = k % 17
            who += [w, 39 if (arm == "NACL_INF" and rep == 0) else w]     # par 17/18: dos personas
            n += 2; k += 1
    out = set(SETS["Clock output"]["genes"]) | {"C1orf51"}; ins = set(SETS["Insulin core"]["genes"])
    X = np.column_stack([person[w] + rng.normal(0, .25, len(genes)) for w in who])
    for j, t in enumerate(titles):
        if "_240_" in t:
            for i, g in enumerate(genes):
                if g in out: X[i, j] -= 0.6
                if g in ins and ("_EU_" in t or "_HC_" in t): X[i, j] -= 0.6
    acc = [f"GSM{997377 + i}" for i in range(len(titles))]
    with gzip.open(os.path.join(raw, "GSE40589_series_matrix.txt.gz"), "wt") as f:
        f.write("!Sample_title\t" + "\t".join(f'"{t}"' for t in titles) + "\n")
        f.write("!Sample_geo_accession\t" + "\t".join(f'"{a}"' for a in acc) + "\n")
        f.write("!series_matrix_table_begin\n" + '"ID_REF"\t' + "\t".join(f'"{a}"' for a in acc) + "\n")
        for p, r in zip(probes, X):
            f.write(f'"{p}"\t' + "\t".join(f"{v:.5f}" for v in r) + "\n")
        f.write("!series_matrix_table_end\n")
    with open(os.path.join(raw, "GPL6848_annot.txt"), "w") as f:
        f.write("ID\tSPOT_ID\tCONTROL_TYPE\tREFSEQ\tGB_ACC\tGENE\tGENE_SYMBOL\tGENE_NAME\n")
        for p, g in zip(probes, genes):
            f.write(f"{p}\t\tFALSE\t\t\t1\t{g}\tname\n")
    return len(titles)


with tempfile.TemporaryDirectory() as tmp:
    raw, out = os.path.join(tmp, "raw"), os.path.join(tmp, "out")
    os.makedirs(raw)
    n_arrays = write_series(raw)
    env = dict(os.environ, T2D_RAW=raw, T2D_OUT=out)
    r = subprocess.run([sys.executable, os.path.join(ROOT, "python/analyses/28_clamp_vs_saline.py")],
                       env=env, capture_output=True, text=True)
    print("1. El paso se ejecuta")
    check("termina sin error", r.returncode == 0, r.stderr.strip().split("\n")[-1][:120] if r.returncode else f"{n_arrays} arreglos")
    if r.returncode == 0:
        new = os.path.join(out, "new")
        read = lambda name: pd.read_csv(os.path.join(new, name), sep="\t")
        P, S, U, G, Q = (read(f"clamp_vs_saline_{k}.tsv") for k in ("pairs", "scores", "unpaired", "genes", "pairing_pairs"))

        print("2. Pareo y regla de identidad")
        check("21 pares en los tres brazos del analisis", len(P) == 21 and set(P.arm) == {"NACL_INF", "NACL_EU", "NACL_HC"}, f"{len(P)}")
        check("el par formado por dos personas queda excluido", list(P.pair[P.excluded]) == ["17/18"], str(list(P.pair[P.excluded])))
        check("ningun par bien formado queda excluido", int(P.excluded.sum()) == 1)
        check("el excluido tiene la correlacion de identidad mas baja", Q.set_index("pair").r_identity.idxmin() == "17/18")

        print("3. Contenido de las tablas")
        cols = {"pair", "arm", "excluded", "d_Clock output", "d_Insulin core", "d_Glucocorticoid", "d_Core oscillator"}
        check("columnas de la tabla por par", cols <= set(P.columns), str(sorted(cols - set(P.columns))))
        check("la tabla de resumen lleva los dos analisis", set(S.analysis) == {"primary", "all pairs"})
        check("el principal usa 20 pares", int(S[(S.analysis == "primary") & (S.set == "Clock output") & S.arm.isin(["saline", "euglycaemic clamp", "hyperglycaemic clamp"])].n.sum()) == 20)
        check("gen a gen: los nueve genes de salida, con el simbolo antiguo resuelto", (G.set == "Clock output").sum() == 9 and "C1orf51" in set(G.gene))

        print("4. Recupera lo que se simulo")
        prim = S[S.analysis == "primary"].set_index(["set", "arm"])
        co = prim.loc["Clock output"]
        check("la salida del reloj cae en los tres brazos", all(co.loc[a, "median"] < -0.3 and co.loc[a, "n_down"] == co.loc[a, "n"]
                                                               for a in ("saline", "euglycaemic clamp", "hyperglycaemic clamp")))
        check("y no difiere entre clamp y salina", co.loc["euglycaemic clamp - saline", "p_mannwhitney"] > 0.05,
              f"P = {co.loc['euglycaemic clamp - saline', 'p_mannwhitney']:.3f}")
        ic = prim.loc["Insulin core"]
        check("el conjunto de insulina no cae con salina", abs(ic.loc["saline", "median"]) < 0.5, f"{ic.loc['saline', 'median']:+.2f}")
        check("y si con los dos clamps", all(ic.loc[f"{a} - saline", "p_mannwhitney"] < 0.01 and ic.loc[f"{a} - saline", "median"] < 0
                                            for a in ("euglycaemic clamp", "hyperglycaemic clamp")))
        u = U[(U.set == "Clock output") & U.arm.isin(["saline", "euglycaemic clamp", "hyperglycaemic clamp"])]
        check("sin usar el pareo: la misma caida, con los arreglos sueltos", bool((u.change < -0.3).all()) and int(u.n0.sum() + u.n240.sum()) == 45,
              f"{int(u.n0.sum() + u.n240.sum())} arreglos")

print("\nTODO CORRECTO" if not fails else f"\nFALLAN {len(fails)}: {fails}")
sys.exit(1 if fails else 0)
