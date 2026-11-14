#!/usr/bin/env python3
"""Clamp de insulina frente a infusion salina a las mismas horas, en adiposo humano (GSE40589).

Es el control horario del clamp: Tuvia et al. (Diabetes 2021;70:1985) tomaron biopsias de tejido
adiposo subcutaneo a los -40 y a los 240 min de una infusion que empieza a las 8:00, con salina
sola, con clamp hiperinsulinemico euglucemico o con clamp hiperglucemico. Los tres brazos
comparten las horas, de modo que el cambio bajo salina es el paso del tiempo y la diferencia
entre un clamp y la salina es lo que anade la insulina.

Titulos de GEO: SCAT_<minuto>_<infusion>_<condicion>_<numero>
    NACL_INF  salina sola                      -> tiempo
    NACL_EU   salina + clamp euglucemico       -> tiempo + insulina
    NACL_HC   salina + clamp hiperglucemico    -> tiempo + insulina + glucosa
    GIP_*     brazos con GIP; solo se usan para comprobar el pareo

PAREO Y REGLA DE IDENTIDAD. La serie no declara sujeto. El par antes/despues se toma de la
numeracion consecutiva del titulo (n = minuto 0, n + 1 = minuto 240, misma infusion y condicion)
y se comprueba con genes que identifican a un individuo y no dependen del lote ni de la hora
(delecion de GSTM1 y GSTT1, ERAP2 y HLA, del registro de conjuntos): las dos biopsias de un par
deben correlacionar en esos genes mas que dos biopsias de personas distintas. La regla, fijada
antes de mirar el resultado, es:

    un par se excluye del analisis principal si su correlacion de identidad queda por debajo
    del percentil 90 de las correlaciones entre pares distintos de toda la serie.

Se usa el valor y no el puesto porque el ensayo es cruzado: una persona puede estar en varios
brazos y su pareja compite con otras biopsias suyas igual de parecidas.

Se reportan tres analisis: el principal (pares que pasan la regla), todos los pares, y uno que no
usa el pareo (minuto 240 frente a minuto 0 como dos grupos, con todos los arreglos de cada brazo,
incluidos los que no tienen pareja).

Los valores depositados son log2 normalizados por gen, asi que la diferencia 240 - 0 dentro de un
par es un log2 del cambio.

Salidas en results/new/:
    clamp_vs_saline_pairs.tsv           cambio por par y por conjunto, con la marca de exclusion
    clamp_vs_saline_scores.tsv          resumen por brazo y contraste clamp - salina (principal y todos)
    clamp_vs_saline_unpaired.tsv        el mismo contraste sin usar el pareo
    clamp_vs_saline_genes.tsv           gen a gen para el reloj (analisis principal)
    clamp_vs_saline_pairing.tsv         resumen de la comprobacion del pareo
    clamp_vs_saline_pairing_pairs.tsv   correlacion de identidad de cada par
    clamp_vs_saline_pairing_matrix.tsv  matriz de correlaciones de identidad, antes x despues
"""
import os, sys, gzip, io, warnings; warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
from gene_sets import SETS, resolve
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results")

ACC = "GSE40589"
ARMS = {"NACL_INF": "saline", "NACL_EU": "euglycaemic clamp", "NACL_HC": "hyperglycaemic clamp"}
SETS_USED = ["Clock output", "Core oscillator", "Glucocorticoid", "Insulin core"]
IDENTITY_PERCENTILE = 90          # la regla de exclusion, ver la cabecera
SEED, NBOOT, NPERM = 0, 10000, 10000


def read_series():
    f = f"{RAW}/{ACC}_series_matrix.txt.gz"
    if not os.path.exists(f):
        raise SystemExit(f"falta {f} (lo baja download_public.sh)")
    with gzip.open(f, "rt", errors="ignore") as fh:
        lines = fh.read().split("\n")
    row = lambda key: [x.strip().strip('"').strip() for x in
                       next(l for l in lines if l.startswith(key)).split("\t")[1:]]
    acc, tit = row("!Sample_geo_accession"), row("!Sample_title")
    i = next(k for k, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    j = next(k for k, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    X = pd.read_csv(io.StringIO("\n".join(lines[i + 1:j])), sep="\t", index_col=0)
    X.columns = [c.strip('"') for c in X.columns]; X.index = [str(k).strip('"') for k in X.index]
    M = pd.DataFrame({"gsm": acc, "title": tit})
    p = M.title.str.extract(r"SCAT_(\d+)_([A-Z]+)_([A-Z]+)_(\d+)")
    M["minute"], M["arm"], M["num"] = p[0].astype(float), p[1] + "_" + p[2], p[3].astype(float)
    return X.apply(pd.to_numeric, errors="coerce"), M.dropna(subset=["num"]).reset_index(drop=True)


def annotate(X):
    """Sonda Agilent -> simbolo, con la tabla de la plataforma (GPL6848)."""
    f = next((p for p in (f"{RAW}/GPL6848_annot.txt", f"{RAW}/GPL6848_annot.txt.gz") if os.path.exists(p)), None)
    if f is None:
        raise SystemExit("falta data/raw/GPL6848_annot.txt (lo baja download_public.sh)")
    a = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    col = next((c for c in a.columns if c.upper().replace(" ", "_") in ("GENE_SYMBOL", "SYMBOL")), None)
    if col is None:
        raise SystemExit(f"GPL6848: no encuentro la columna de simbolo; columnas: {list(a.columns)}")
    sym = a.dropna(subset=[col]).set_index(a.columns[0])[col]
    keep = [i for i in X.index if i in sym.index]
    G = X.loc[keep].groupby(sym.reindex(keep).values).mean()
    print(f"anotacion: {len(keep)} de {X.shape[0]} sondas -> {G.shape[0]} genes")
    return G


def pairs(M):
    """(n, n + 1) con minuto 0 y 240 en el mismo brazo. Devuelve pares y arreglos sueltos."""
    out, used = [], set()
    by = {int(r.num): r for r in M.itertuples()}
    for n, r in sorted(by.items()):
        s = by.get(n + 1)
        if r.minute == 0 and s is not None and s.minute == 240 and s.arm == r.arm:
            out.append(dict(pair=f"{n}/{n + 1}", arm=r.arm, pre=r.gsm, post=s.gsm)); used |= {n, n + 1}
    loose = M[~M.num.astype(int).isin(used) & M.arm.isin(ARMS)]
    return pd.DataFrame(out), loose


def check_pairing(G, P_all, rng):
    """Aplica la regla de identidad de la cabecera a todos los pares completos de la serie y
    devuelve la lista de pares que no la pasan. Como referencia reporta tambien los 500 genes mas
    variables una vez quitada la hora, que arrastran lote y por eso no deciden nada."""
    pre, post = list(P_all.pre), list(P_all.post)
    Gc = pd.concat([G[pre].sub(G[pre].mean(1), axis=0), G[post].sub(G[post].mean(1), axis=0)], axis=1)
    Gz = Gc.div(Gc.std(1) + 1e-9, axis=0)
    ident = [r for r in (resolve(x, G.index) for x in SETS["Genetic"]["genes"] + SETS["HLA"]["genes"]) if r]
    top = list(Gc.var(1).sort_values(ascending=False).index[:500])
    rows = []; detail = P_all[["pair", "arm"]].copy()
    print(f"\n== comprobacion del pareo ({len(P_all)} pares completos en toda la serie)")
    for key, label, genes in (("identity", "identity genes (registry)", ident), ("top500", "500 most variable genes", top)):
        A, B = Gz.loc[genes, pre].values, Gz.loc[genes, post].values
        n = A.shape[1]
        C = np.corrcoef(A.T, B.T)[:n, n:]
        diag = np.diag(C); off = C[~np.eye(n, dtype=bool)]
        null = np.array([C[np.arange(n), rng.permutation(n)].mean() for _ in range(NPERM)])
        pval = (np.sum(null >= diag.mean()) + 1) / (NPERM + 1)
        thr = np.percentile(off, IDENTITY_PERCENTILE)
        detail[f"r_{key}"] = diag
        detail[f"others_as_similar_{key}"] = [int(np.sum(np.delete(C[i], i) >= C[i, i])) for i in range(n)]
        print(f"  {label} ({len(genes)}): r dentro del par, mediana {np.median(diag):+.2f}; entre pares distintos, "
              f"mediana {np.median(off):+.2f} y percentil {IDENTITY_PERCENTILE} {thr:+.2f}; permutacion P = {pval:.4f}")
        rows.append(dict(gene_set=label, n_genes=len(genes), n_pairs=n, r_within_median=np.median(diag),
                         r_between_median=np.median(off), r_between_threshold=thr, p_permutation=pval,
                         pairs_below_threshold=int((diag < thr).sum())))
        if key == "identity":
            thr_id = thr
            C_id = pd.DataFrame(C, index=[f"pre {x}" for x in P_all.pair], columns=[f"post {x}" for x in P_all.pair])
    detail["excluded"] = detail.r_identity < thr_id
    used = detail[detail.arm.isin(ARMS)]
    print("  pares del analisis, genes de identidad (r dentro del par; otras biopsias igual o mas parecidas):")
    for a in ARMS:
        q = used[used.arm == a]
        print(f"    {ARMS[a]:22s} " + "  ".join(f"{r.pair} {r.r_identity:+.2f} ({r.others_as_similar_identity})"
                                                  + ("*" if r.excluded else "") for r in q.itertuples()))
    bad = list(used.pair[used.excluded])
    print(f"  excluidos del analisis principal (r por debajo de {thr_id:+.2f}, marcados con *): {bad if bad else 'ninguno'}")
    pd.DataFrame(rows).to_csv(f"{OUT}/new/clamp_vs_saline_pairing.tsv", sep="\t", index=False)
    detail.to_csv(f"{OUT}/new/clamp_vs_saline_pairing_pairs.tsv", sep="\t", index=False)
    C_id.to_csv(f"{OUT}/new/clamp_vs_saline_pairing_matrix.tsv", sep="\t")
    return bad


def boot_diff(a, b, rng):
    """Diferencia de medianas a - b con IC 95 % por bootstrap, y cociente de medianas b / a."""
    d = np.array([np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(NBOOT)])
    q = np.array([np.median(rng.choice(b, len(b))) / np.median(rng.choice(a, len(a))) for _ in range(NBOOT)])
    q = q[np.isfinite(q)]
    return (np.median(a) - np.median(b), *np.percentile(d, [2.5, 97.5]),
            np.median(b) / np.median(a), *np.percentile(q, [2.5, 97.5]))


def paired_summary(per_pair, measured, rng, analysis, verbose=True):
    """Cambio por par en cada brazo y contraste de cada clamp con la salina."""
    rows = []
    for s in SETS_USED:
        v = {a: per_pair.loc[per_pair.arm == a, f"d_{s}"].values for a in ARMS}
        if verbose:
            print(f"\n== {s}  ({len(measured[s])} genes: {', '.join(measured[s])})")
        else:
            print(f"  {s}: " + "; ".join(f"{ARMS[a]} {np.median(v[a]):+.2f}, bajan {int((v[a] < 0).sum())}/{len(v[a])}"
                                           for a in ARMS))
        for a, lab in ARMS.items():
            x = v[a]; w = stats.wilcoxon(x).pvalue if len(x) >= 3 else np.nan
            if verbose:
                print(f"  {lab:22s} n={len(x)}  cambio por par: {np.round(np.sort(x), 2).tolist()}")
                print(f"  {'':22s} mediana {np.median(x):+.2f}; bajan {int((x < 0).sum())}/{len(x)}; Wilcoxon P = {w:.4f}")
            rows.append(dict(analysis=analysis, set=s, arm=lab, n=len(x), median=np.median(x),
                             n_down=int((x < 0).sum()), p_wilcoxon=w))
        txt = []
        for a in ("NACL_EU", "NACL_HC"):
            mw = stats.mannwhitneyu(v[a], v["NACL_INF"]).pvalue
            d, lo, hi, q, qlo, qhi = boot_diff(v[a], v["NACL_INF"], rng)
            # el cociente solo tiene sentido si el clamp mueve el conjunto; si no, se deja vacio
            if not (len(v[a]) >= 3 and stats.wilcoxon(v[a]).pvalue < 0.05):
                q = qlo = qhi = np.nan
            if verbose:
                print(f"  {ARMS[a]} - saline: {d:+.2f} (IC95 {lo:+.2f} a {hi:+.2f}); Mann-Whitney P = {mw:.4f}")
                if np.isfinite(q):
                    print(f"  fraccion del cambio del clamp que ya aparece con salina: {q:.2f} (IC95 {qlo:.2f} a {qhi:.2f})")
            else:
                txt.append(f"{ARMS[a]} - saline P = {mw:.4f}")
            rows.append(dict(analysis=analysis, set=s, arm=f"{ARMS[a]} - saline", n=len(v[a]) + len(v["NACL_INF"]),
                             median=d, ci_lo=lo, ci_hi=hi, p_mannwhitney=mw,
                             fraction_under_saline=q, fraction_lo=qlo, fraction_hi=qhi))
        if not verbose:
            print("      " + "; ".join(txt))
    return rows


def unpaired(G, M, measured, rng):
    """El mismo contraste sin usar el pareo: todos los arreglos de cada brazo, minuto 240 frente
    a minuto 0 como dos grupos, incluidos los arreglos sin pareja. La comparacion entre brazos es
    una diferencia de diferencias de medianas; su P sale de permutar la etiqueta de brazo dentro
    de cada minuto."""
    Mu = M[M.arm.isin(ARMS)].copy()
    Z = G[list(Mu.gsm)]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T
    rows = []
    print("\n== sin usar el pareo: minuto 240 frente a minuto 0, todos los arreglos de cada brazo")
    for s in SETS_USED:
        Mu["score"] = Z.loc[measured[s]].mean(0).values
        line = []
        for a, lab in ARMS.items():
            x0 = Mu.score[(Mu.arm == a) & (Mu.minute == 0)].values; x1 = Mu.score[(Mu.arm == a) & (Mu.minute == 240)].values
            chg = np.median(x1) - np.median(x0); p = stats.mannwhitneyu(x1, x0).pvalue
            line.append(f"{lab} {chg:+.2f} (n = {len(x0)} y {len(x1)}, P = {p:.4f})")
            rows.append(dict(set=s, arm=lab, n0=len(x0), n240=len(x1), change=chg, p_mannwhitney=p))
        print(f"  {s}: " + "; ".join(line))
        for a in ("NACL_EU", "NACL_HC"):
            sub = Mu[Mu.arm.isin([a, "NACL_INF"])]
            minute, score = sub.minute.values, sub.score.values
            def did(arm):
                d = {k: np.median(score[(arm == k) & (minute == 240)]) - np.median(score[(arm == k) & (minute == 0)])
                     for k in (a, "NACL_INF")}
                return d[a] - d["NACL_INF"]
            obs = did(sub.arm.values); null = np.empty(NPERM)
            for i in range(NPERM):
                arm = sub.arm.values.copy()
                for m in (0, 240):
                    idx = np.where(minute == m)[0]; arm[idx] = rng.permutation(arm[idx])
                null[i] = did(arm)
            p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)
            print(f"      {ARMS[a]} - saline: {obs:+.2f}, permutacion P = {p:.4f}")
            rows.append(dict(set=s, arm=f"{ARMS[a]} - saline", change=obs, p_permutation=p))
    pd.DataFrame(rows).to_csv(f"{OUT}/new/clamp_vs_saline_unpaired.tsv", sep="\t", index=False)


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    X, M = read_series()
    G = annotate(X)
    P_all, loose = pairs(M)
    P_all = P_all.reset_index(drop=True)
    P = P_all[P_all.arm.isin(ARMS)].reset_index(drop=True)
    print(f"{ACC}: {X.shape[1]} arreglos; pares completos por brazo: {P.arm.value_counts().to_dict()}")
    if len(loose):
        print("arreglos sin pareja (solo entran en el analisis sin pareo):", ", ".join(loose.title))

    rng = np.random.default_rng(SEED)
    excluded = check_pairing(G, P_all, rng)

    # z por gen sobre los arreglos pareados de los tres brazos, para que los conjuntos pesen igual
    cols = list(P.pre) + list(P.post)
    Z = G[cols]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T
    D_log2 = pd.DataFrame(G[list(P.post)].values - G[list(P.pre)].values, index=G.index, columns=P.pair)
    D_z = pd.DataFrame(Z[list(P.post)].values - Z[list(P.pre)].values, index=G.index, columns=P.pair)

    per_pair = P[["pair", "arm"]].copy(); measured = {}
    per_pair["excluded"] = per_pair.pair.isin(excluded)
    for s in SETS_USED:
        g = [r for r in (resolve(x, G.index) for x in SETS[s]["genes"]) if r]
        measured[s] = g
        per_pair[f"d_{s}"] = D_z.loc[g].mean(0).values
        per_pair[f"dlog2_{s}"] = D_log2.loc[g].mean(0).values
    per_pair.to_csv(f"{OUT}/new/clamp_vs_saline_pairs.tsv", sep="\t", index=False)

    primary = per_pair[~per_pair.excluded]
    print(f"\n#### ANALISIS PRINCIPAL: {len(primary)} pares que pasan la regla de identidad")
    rows = paired_summary(primary, measured, rng, "primary")
    if excluded:
        print(f"\n#### SENSIBILIDAD: los {len(per_pair)} pares, incluidos los excluidos {excluded}")
        rows += paired_summary(per_pair, measured, rng, "all pairs", verbose=False)
    pd.DataFrame(rows).to_csv(f"{OUT}/new/clamp_vs_saline_scores.tsv", sep="\t", index=False)

    print("\n#### SENSIBILIDAD", end="")
    unpaired(G, M, measured, rng)

    # gen a gen, en log2, para los genes del reloj, con los pares del analisis principal
    ok = ~P.pair.isin(excluded).values
    grow = []
    print("\n== gen a gen (log2 del cambio 240 - 0; mediana y pares que bajan; analisis principal)")
    print(f"  {'gen':8s} {'salina':>16s} {'clamp eu':>16s} {'clamp - salina':>15s} {'P':>7s}")
    for s in ("Clock output", "Core oscillator"):
        for g in measured[s]:
            a = D_log2.loc[g, (P.arm == "NACL_INF").values & ok].values
            b = D_log2.loc[g, (P.arm == "NACL_EU").values & ok].values
            p = stats.mannwhitneyu(b, a).pvalue
            print(f"  {g:8s} {np.median(a):+7.2f} ({int((a < 0).sum())}/{len(a)}) {np.median(b):+7.2f} ({int((b < 0).sum())}/{len(b)})"
                  f" {np.median(b) - np.median(a):+14.2f} {p:7.3f}")
            grow.append(dict(set=s, gene=g, saline_median=np.median(a), saline_down=int((a < 0).sum()), saline_n=len(a),
                             clamp_median=np.median(b), clamp_down=int((b < 0).sum()), clamp_n=len(b),
                             difference=np.median(b) - np.median(a), p_mannwhitney=p))
    pd.DataFrame(grow).to_csv(f"{OUT}/new/clamp_vs_saline_genes.tsv", sep="\t", index=False)
    print(f"\n-> {OUT}/new/clamp_vs_saline_*.tsv")
