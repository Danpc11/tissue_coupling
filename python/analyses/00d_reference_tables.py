#!/usr/bin/env python3
"""Productores de las tablas de `data/reference/` que el repositorio publicado no regeneraba.

Cada funcion recalcula una tabla desde los datos publicos y las
matrices de 00b_build_matrices.py, la escribe en results/new/ con el mismo nombre y la compara
con la version versionada en data/reference/ (filas comunes, correlacion por columna numerica).
El resumen va a results/new/reference_tables_check.tsv.

Convenciones (Methods): las diferencias dentro de persona o entre grupos se calculan sin centrar
por muestra; los scores entre individuos, centrados. Donde la tabla de referencia se calculo sin
centrar, la comparacion se hace sin centrar y la tabla escrita sigue la politica de la rama
(columnas adicionales o version centrada), segun se indica en cada funcion.

Uso:  python python/analyses/00d_reference_tables.py [tabla ...]
"""
import os, sys, gzip, glob, re, warnings, importlib.util
warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../lib"))
from gene_sets import SETS, resolve
import phase as ph
import gpl570, islets

RAW = os.environ.get("T2D_RAW", "data/raw")
OUT = os.environ.get("T2D_OUT", "results")
REF = os.environ.get("T2D_REF", "data/reference")
NEW = f"{OUT}/new"
CENTER = os.environ.get("T2D_CENTER", "1") != "0"
CLK = SETS["Clock output"]["genes"]
CHECK = []


def mod(name, fn):
    sp = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m


def sc(X, g):
    gg = [r for r in (resolve(x, X.index) for x in g) if r]
    Z = X.loc[gg]; return ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.mean(0)


def center(X):
    return X.sub(X.mean(0), axis=1)


def pres(y, cv):
    D = np.column_stack([np.ones(len(y))] + list(cv)) if len(cv) else np.ones((len(y), 1))
    b, *_ = np.linalg.lstsq(D, np.asarray(y, float), rcond=None); return np.asarray(y, float) - D @ b


def aagmex_cov():
    A = pd.read_pickle(f"{NEW}/AAGMEx_adipose_log2.pkl")
    with gzip.open(f"{RAW}/GSE95674_series_matrix.txt.gz", "rt") as f:
        for l in f:
            if l.startswith("!Sample_title"): t = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
            if l.startswith("!Sample_characteristics") and "gender" in l: sx = [x.strip('"').split(": ")[1] for x in l.rstrip().split("\t")[1:]]
            if l.startswith("!Sample_characteristics") and "age" in l: ag = [float(x.strip('"').split(": ")[1]) for x in l.rstrip().split("\t")[1:]]
    cov = pd.DataFrame({"sex": sx, "age": ag}, index=t).reindex(A.columns)
    return A, [(cov.sex == "M").astype(float).values, cov.age.values]


# ----------------------------------------------------------------------------- comparacion
def read_ref(name):
    f = f"{REF}/{name}"
    first = open(f).readline().split("\t")[0]
    return pd.read_csv(f, sep="\t", index_col=0 if first == "" else None)


def compare(name, new, key=None, note=""):
    """Compara `new` con la tabla versionada. key: columnas que identifican la fila (None = indice)."""
    ref = read_ref(name)
    if key:
        r = ref.set_index(key); n = new.set_index(key)
    else:
        r, n = ref, new
    common = r.index.intersection(n.index)
    cols = [c for c in r.columns if c in n.columns and pd.api.types.is_numeric_dtype(r[c])]
    stats_ = []
    for c in cols:
        a = pd.to_numeric(r.loc[common, c], errors="coerce"); b = pd.to_numeric(n.loc[common, c], errors="coerce")
        m = a.notna() & b.notna()
        if m.sum() >= 3 and a[m].std() > 0 and b[m].std() > 0:
            if str(c).lower().startswith("phase"):
                # una fase es circular: la correlacion lineal no la mide. Se compara la concentracion de
                # la diferencia, permitiendo los dos sentidos de giro y un desfase global, y el error
                # que queda una vez quitado ese desfase (el azar da 6 h).
                x, y = a[m].values * 2 * np.pi / 24, b[m].values * 2 * np.pi / 24
                R, way = max((abs(np.mean(np.exp(1j * (x - y)))), "same direction"), (abs(np.mean(np.exp(1j * (x + y)))), "reversed"))
                d = x - y if way == "same direction" else x + y
                e = np.abs(np.angle(np.exp(1j * (d - np.angle(np.mean(np.exp(1j * d))))))) * 24 / (2 * np.pi)
                stats_.append((f"{c} (circular: concentration of the difference, {way}; median error after offset)", R, float(np.median(e))))
            else:
                stats_.append((c, np.corrcoef(a[m], b[m])[0, 1], float((a[m] - b[m]).abs().median())))
    rmin = min([s[1] for s in stats_], default=np.nan)
    row = dict(tabla=name, filas_ref=len(r), filas_nueva=len(n), filas_comunes=len(common),
               columnas_comparadas=len(stats_), r_minimo=rmin,
               detalle="; ".join(f"{c}: r={rr:.3f}, |d| med={d:.3g}" for c, rr, d in stats_), nota=note)
    CHECK.append(row)
    print(f"  [{name}] ref {len(r)} / nueva {len(n)} / comunes {len(common)}; r min {rmin:.3f}" if stats_ else
          f"  [{name}] ref {len(r)} / nueva {len(n)} / comunes {len(common)}; sin columnas numericas comparables")
    for c, rr, d in stats_: print(f"      {c:40s} r = {rr:+.3f}   |diferencia| mediana {d:.3g}")
    new.to_csv(f"{NEW}/{name}", sep="\t", index=key is None)
    return row


# ----------------------------------------------------------------------------- tablas
def t_within_donor():
    """within_donor_decomposition.tsv (Fig. 2d). rho = NR3C1 ~ score de salida del reloj sin
    ajustar; rho_sin_nivel_donante = lo mismo tras restar a NR3C1 la media del donante en los
    tejidos que tiene. Referencia sin centrar; se escriben las dos versiones."""
    rows = []
    for cen in (False, True):
        G = {t: pd.read_pickle(f"{NEW}/{fn}.pkl") for t, fn in (("adiposo sc", "GTEx_adipose_sc"), ("adiposo visceral", "GTEx_adipose_visc"), ("músculo", "GTEx_muscle"))}
        if cen: G = {t: center(X) for t, X in G.items()}
        nr = pd.DataFrame({t: X.loc["NR3C1"] for t, X in G.items()})
        multi = nr.notna().sum(1) >= 2
        for t, X in G.items():
            co = sc(X, CLK); x = X.loc["NR3C1"]
            d = nr.index[multi & nr[t].notna()]
            resid = pd.Series(pres(nr.loc[d, t], [nr.loc[d].mean(1).values]), index=d)   # residuo de la regresion sobre la media del donante
            rows.append(dict(cohorte="GTEx", tejido=t, rho=stats.spearmanr(x, co).statistic,
                             rho_sin_nivel_donante=stats.spearmanr(resid, co.reindex(d)).statistic, centrado=cen))
        Ea = pd.read_pickle(f"{NEW}/GSE41168_adipose.pkl"); Em = pd.read_pickle(f"{NEW}/GSE41168_muscle.pkl")
        if cen: Ea, Em = center(Ea), center(Em)
        pair = [d for d in Ea.columns if d in Em.columns]
        m = (Ea.loc["NR3C1", pair] + Em.loc["NR3C1", pair]) / 2
        for t, X in (("adiposo", Ea), ("músculo", Em)):
            co = sc(X, CLK)
            rows.append(dict(cohorte="GSE41168", tejido=t, rho=stats.spearmanr(X.loc["NR3C1"], co).statistic,
                             rho_sin_nivel_donante=stats.spearmanr(pres(X.loc["NR3C1", pair], [m.values]), co[pair]).statistic, centrado=cen))
    D = pd.DataFrame(rows).round(2)
    compare("within_donor_decomposition.tsv", D[~D.centrado].drop(columns="centrado"), key=["cohorte", "tejido"],
            note="referencia sin centrar y con 28 parejas de GSE41168; aqui residuo de la regresion sobre la media del donante")
    out = D[D.centrado].drop(columns="centrado") if CENTER else D[~D.centrado].drop(columns="centrado")
    out.to_csv(f"{NEW}/within_donor_decomposition.tsv", sep="\t", index=False)


def t_gcreceptor_adiposity():
    """AAGMEx_GCreceptor_vs_adiposity.tsv (Fig. 4d): rho parcial (sexo, edad) de cada gen con LEP
    y con la salida del reloj, AAGMEx adiposo."""
    A, cv = aagmex_cov()
    if CENTER: A = center(A)
    genes = ["NR3C1", "FKBP5", "TSC22D3", "ZBTB16", "KLF15", "HSD11B1", "PER3", "DBP", "HLF"]   # Fig. 4d
    co = pres(sc(A, CLK), cv); lep = pres(A.loc["LEP"], cv)
    rows = [dict(gene=g, rho_vs_LEP=stats.spearmanr(pres(A.loc[g], cv), lep).statistic,
                 rho_vs_clock=stats.spearmanr(pres(A.loc[g], cv), co).statistic) for g in genes if g in A.index]
    compare("AAGMEx_GCreceptor_vs_adiposity.tsv", pd.DataFrame(rows).round(3), key=["gene"])


def rank_vs_nr3c1():
    A, cv = aagmex_cov()
    Xc = center(A)
    D = np.column_stack([np.ones(Xc.shape[1])] + cv + [Xc.loc["LEP"].values, Xc.loc["CD68"].values])
    Q, _ = np.linalg.qr(D); V = Xc.values.T; Vr = V - Q @ (Q.T @ V)
    rk = np.apply_along_axis(stats.rankdata, 0, Vr)
    nr = rk[:, list(Xc.index).index("NR3C1")]
    z = (rk - rk.mean(0)) / rk.std(0); nz = (nr - nr.mean()) / nr.std()
    return pd.Series((z * nz[:, None]).mean(0), index=Xc.index).dropna().rename("rho")


def t_rank_vs_nr3c1():
    """AAGMEx_rank_vs_NR3C1.tsv (Fig. 4e): rho de cada gen con NR3C1 ajustado por sexo, edad, LEP y
    CD68 (el mismo calculo de 09_hpo_by_tissue.py, que no lo guardaba)."""
    R = rank_vs_nr3c1().sort_values(ascending=False)
    compare("AAGMEx_rank_vs_NR3C1.tsv", R.to_frame())


def gmt(path):
    out = {}
    for l in open(path):
        p = l.rstrip("\n").split("\t"); out[p[0]] = set(p[2:])
    return out


LIBS = {"h_all": ("h.all", None), "c2_cp_reactome": ("c2.cp.reactome", None), "c2_cp_wikipathways": ("c2.cp.wikipathways", None),
        "c2_cp_kegg_medicus": ("c2.cp.kegg_medicus", None), "c2_cp_pid": ("c2.cp.pid", None), "c2_cgp": ("c2.cgp", None),
        "c3_tft_gtrd": ("c3.tft.gtrd", None), "c3_mir_mirdb": ("c3.mir.mirdb", None), "c4_all": ("c4.all", None),
        "c5_go_bp": ("c5.go.bp", None), "c5_hpo": ("c5.hpo", None), "c7_all": ("c7.all", None), "c8_all": ("c8.all", None)}


def enrichment(R, libs, nperm=2000, seed=0, min_n=10, anchor="NR3C1"):
    """z de la media de rho de cada conjunto frente a 2,000 conjuntos aleatorios del mismo tamano.

    El ranking es la correlacion de cada gen con `anchor`, y la del ancla consigo misma vale 1. Un
    conjunto que contiene al ancla recibe ese 1 en su media: es circular, y en conjuntos pequenos
    pesa (hasta dos unidades de z). La columna `z` excluye al ancla de los conjuntos que la
    contienen; `z_including_anchor` es el calculo anterior, tal como se publico, y
    `contains_anchor` dice a cuales afecta. En los demas conjuntos las dos columnas coinciden."""
    rng = np.random.default_rng(seed); vals = R.values; idx = {g: i for i, g in enumerate(R.index)}
    ia = idx.get(anchor); rng2 = np.random.default_rng(seed + 1); null2 = {}
    pool = np.array([i for i in range(len(vals)) if i != ia])            # el ancla tampoco entra en el nulo corregido
    null = {}
    rows = []
    for lib in libs:
        fn = f"{RAW}/msigdb/{LIBS[lib][0]}.v2026.2.Hs.symbols.gmt"
        if not os.path.exists(fn): print(f"    (sin {fn})"); continue
        for name, genes in gmt(fn).items():
            ii = [idx[g] for g in genes if g in idx]
            n = len(ii)
            if n < min_n: continue
            if n not in null:
                null[n] = np.array([vals[rng.choice(len(vals), n, replace=False)].mean() for _ in range(nperm)])
            z0 = (vals[ii].mean() - null[n].mean()) / null[n].std(); z = z0; has = ia is not None and ia in ii
            if has:
                jj = [i for i in ii if i != ia]; m = len(jj)
                if m not in null2:
                    null2[m] = np.array([vals[rng2.choice(pool, m, replace=False)].mean() for _ in range(nperm)])
                z = (vals[jj].mean() - null2[m].mean()) / null2[m].std()
            rows.append(dict(lib=lib, set=name, n=n, z=z, z_including_anchor=z0, contains_anchor=has))
    return pd.DataFrame(rows)


def t_enrichment():
    """AAGMEx_NR3C1_enrichment.tsv y _extra.tsv (Fig. 4b, f, g): enriquecimiento de colecciones
    MSigDB en el ranking de AAGMEx_rank_vs_NR3C1."""
    R = rank_vs_nr3c1()
    E1 = enrichment(R, ["c2_cp_kegg_medicus", "c2_cp_reactome", "c2_cp_wikipathways", "c3_tft_gtrd", "c5_go_bp", "c5_hpo", "c8_all"])
    compare("AAGMEx_NR3C1_enrichment.tsv", E1, key=["lib", "set"], note="c9_all no existe en MSigDB 2026.2")
    E2 = enrichment(R, ["c2_cgp", "c2_cp_pid", "c3_mir_mirdb", "c4_all", "c7_all", "h_all"])
    compare("AAGMEx_NR3C1_enrichment_extra.tsv", E2, key=["lib", "set"])


def cohorts_axis():
    """Las matrices de las diez cohortes-tejido de AUDIT_global_axis.tsv, sin centrar."""
    don, _ = islets.groups(RAW)
    yield "islotes", islets.log2cpm_symbols(f"{RAW}/GSE159984_raw_counts.tsv.gz", don, f"{RAW}/Human_GRCh38_p13_annot_tsv.gz")
    for lab, fn in (("GTEx adiposo sc", "GTEx_adipose_sc"), ("GTEx visceral", "GTEx_adipose_visc"), ("GTEx músculo", "GTEx_muscle"),
                    ("AAGMEx adiposo", "AAGMEx_adipose_log2"), ("AAGMEx músculo", "AAGMEx_muscle_log2"), ("GSE78721 adiposo", "GSE78721_log2"),
                    ("GSE41168 adiposo", "GSE41168_adipose"), ("GSE41168 músculo", "GSE41168_muscle"), ("GSE25462 músculo", "GSE25462_entrez")):
        yield lab, pd.read_pickle(f"{NEW}/{fn}.pkl")


def t_audit_global_axis():
    """AUDIT_global_axis.tsv (Supp. Fig. 3a): rho NR3C1 ~ reloj sin ajustar, mediana del nulo de
    800 genes al azar y percentil de NR3C1, sin centrar y centrado."""
    rows = []
    for lab, X in cohorts_axis():
        rng = np.random.default_rng(0); genes = rng.choice(X.index, 800, replace=False)
        r = dict(cohorte=lab)
        for suf, Y in (("", X), ("_c", center(X))):
            co = sc(Y, CLK); rho = stats.spearmanr(Y.loc["NR3C1"], co).statistic
            null = np.array([stats.spearmanr(Y.loc[g], co).statistic for g in genes])
            r["rho" + ("_centrado" if suf else "")] = round(rho, 3); r["null_med" + suf] = round(float(np.median(null)), 3)
            r["pct" + suf] = round(float(np.mean(null < rho) * 100), 3)
        rows.append(r)
    compare("AUDIT_global_axis.tsv", pd.DataFrame(rows), key=["cohorte"])


def t_msigdb_coverage():
    """msigdb_symbol_coverage.tsv (Supp. Fig. 3): genes de cada matriz presentes en MSigDB con el
    simbolo de la plataforma y tras mapear sinonimos a NCBI."""
    allsym = set()
    for fn in glob.glob(f"{RAW}/msigdb/*.gmt"): [allsym.update(v) for v in gmt(fn).values()]
    syn = pd.read_csv(f"{NEW}/ncbi_synonym_map.tsv", sep="\t", index_col=0)["current"].to_dict()
    rows = []
    for lab, fn in (("AAGMEx adiposo", "AAGMEx_adipose_log2"), ("GTEx adiposo sc", "GTEx_adipose_sc"), ("GSE41168 adiposo", "GSE41168_adipose"),
                    ("GSE182120 músculo", "GSE182120_log2"), ("GSE66446 adiposo", "GSE66446_log2cpm")):
        g = list(pd.read_pickle(f"{NEW}/{fn}.pkl").index)
        raw = sum(x in allsym for x in g); mapped = len({syn.get(x, x) for x in g} & allsym)
        rows.append(dict(cohorte=lab, genes=len(g), en_msigdb_crudo=raw, en_msigdb_mapeado=mapped, ganancia=mapped - raw))
    compare("msigdb_symbol_coverage.tsv", pd.DataFrame(rows), key=["cohorte"],
            note="las matrices de aqui tienen otro universo de genes (GTEx 39k frente a 15.6k)")


def t_synonym_map():
    """ncbi_synonym_map.tsv: la escribe 00b_build_matrices.py (gene_info); aqui solo se compara."""
    D = pd.read_csv(f"{NEW}/ncbi_synonym_map.tsv", sep="\t", index_col=0)
    ref = read_ref("ncbi_synonym_map.tsv"); common = ref.index.intersection(D.index)
    agree = (ref.loc[common, "current"] == D.loc[common, "current"]).mean() if len(common) else np.nan
    print(f"  [ncbi_synonym_map.tsv] ref {len(ref)} / nueva {len(D)} / comunes {len(common)}; coinciden {agree*100:.1f} %")
    CHECK.append(dict(tabla="ncbi_synonym_map.tsv", filas_ref=len(ref), filas_nueva=len(D), filas_comunes=len(common),
                      columnas_comparadas=1, r_minimo=agree, detalle=f"acuerdo en 'current' {agree*100:.1f} %", nota="gene_info actual (p14), no p13"))


def t_gse26637():
    """GSE26637_clamp_adipose_per_person.tsv (Fig. 3c-d): hiperinsulinemia - ayuno por persona."""
    X = pd.read_pickle(f"{NEW}/GSE26637_log2.pkl"); md = pd.read_pickle(f"{NEW}/GSE26637_meta.pkl").set_index("gsm")
    md["subject"] = md.title.str.split(" ").str[-1]
    D = {}
    for s, q in md.groupby("subject"):
        f = q.index[q.stim == "fasting"]; h = q.index[q.stim != "fasting"]
        if len(f) == 1 and len(h) == 1: D[s] = X[h[0]] - X[f[0]]
    D = pd.DataFrame(D)
    compare("GSE26637_clamp_adipose_per_person.tsv", D)


def t_gse148640_stats():
    """GSE148640_HC_response_stats.tsv (Fig. 3e, Supp. Fig. 6): por gen, n de personas en que sube,
    mediana y media de la diferencia, P de Wilcoxon exacto y FDR."""
    H = pd.read_pickle(f"{OUT}/GSE148640_diff_genes.pkl")
    p = H.apply(lambda r: stats.wilcoxon(r.values, method="exact").pvalue if r.notna().all() and (r != 0).any() else np.nan, axis=1)
    from statsmodels.stats.multitest import multipletests
    D = pd.DataFrame(dict(n_up=(H > 0).sum(1), median=H.median(1), mean=H.mean(1), p=p))
    D["fdr"] = np.nan; m = D.p.notna(); D.loc[m, "fdr"] = multipletests(D.p[m], method="fdr_bh")[1]
    compare("GSE148640_HC_response_stats.tsv", D, note="universo: 29,525 simbolos aqui frente a 12,642")


def t_gras():
    """GRAS_signature.tsv (Supp. Fig. 6): genes que cambian en la misma direccion en todas las personas
    con hidrocortisona, con |mediana| >= 0.2.

    Cuantos de esos genes se esperan por azar: con n personas, un gen sin efecto cambia en el mismo
    sentido en todas con probabilidad 2 * 0.5^n, y entre decenas de miles de genes eso son decenas.
    El nulo invierte el signo de la respuesta de personas enteras (todos los patrones), lo que
    conserva la correlacion entre genes, y cuenta cuantos pasan los dos criterios."""
    import itertools
    H = pd.read_pickle(f"{OUT}/GSE148640_diff_genes.pkl")
    n = H.shape[1]
    up = (H > 0).sum(1); med = H.median(1)
    cons = H[((up == n) | (up == 0)) & (med.abs() >= 0.2)]
    D = pd.DataFrame(dict(gene=cons.index, direction=np.where(med[cons.index] > 0, "up", "down"),
                          median_log2=med[cons.index].round(5), n_consistent=n)).sort_values("median_log2", key=abs, ascending=False)
    V = H.dropna().values; null = []
    for signs in itertools.product([1, -1], repeat=n):
        if signs[0] == -1: break                         # la otra mitad es la misma con el signo global invertido
        if all(x == 1 for x in signs): continue
        W = V * np.array(signs)[None, :]; u = (W > 0).sum(1)
        null.append(int((((u == n) | (u == 0)) & (np.abs(np.median(W, axis=1)) >= 0.2)).sum()))
    null = np.array(null); obs = len(cons)
    pd.DataFrame([dict(n_participants=n, genes_tested=V.shape[0], signature_genes=obs, null_median=float(np.median(null)),
                       null_p95=float(np.percentile(null, 95)), expected_false_fraction=float(np.median(null) / max(obs, 1)),
                       p_signature_size=float((np.sum(null >= obs) + 1) / (len(null) + 1)), sign_patterns=len(null))]
                 ).to_csv(f"{NEW}/GRAS_signature_null.tsv", sep="\t", index=False)
    print(f"  firma: {obs} genes; por azar (inversion de signo de personas enteras, {len(null)} patrones) pasan {np.median(null):.0f} "
          f"(percentil 95: {np.percentile(null, 95):.0f}); fraccion esperada de falsos {np.median(null) / max(obs, 1) * 100:.0f} %; "
          f"P del tamano de la firma {(np.sum(null >= obs) + 1) / (len(null) + 1):.4f}")
    compare("GRAS_signature.tsv", D, key=["gene"], note="umbral |mediana| >= 0.2 deducido de la referencia")


def t_foxo3():
    """FOXO3_regulon_ED4.tsv (Supp. Fig. 4): dianas de FOXO3 en DoRothEA A-C con signo y nivel."""
    d = pd.read_csv(f"{RAW}/dorothea_ABC.tsv", sep="\t", low_memory=False,
                    usecols=["source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition", "dorothea_level"])
    d["mor"] = np.where(d.is_stimulation & ~d.is_inhibition, 1, np.where(d.is_inhibition & ~d.is_stimulation, -1, 0))
    d = d[(d.mor != 0) & (d.source_genesymbol == "FOXO3")].drop_duplicates(["source_genesymbol", "target_genesymbol"])
    D = pd.DataFrame(dict(TF="FOXO3", target=d.target_genesymbol.values, sign=d.mor.values,
                          confidence=d.dorothea_level.str.split(";").str[0].values))
    compare("FOXO3_regulon_ED4.tsv", D, key=["TF", "target"])


def gse13070():
    c5 = mod("c5", "05_clamp_within_person.py"); E, d = c5.load(); return E, d


def clamp_contrasts(E, dd):
    """Mediana entre parejas de (post - pre) por gen: clamp 3 h en musculo basal, TZD 3 meses en
    musculo (sin clamp) y en adiposo. Mismo criterio que 05 para los 45 genes que este guarda."""
    def paired(tis, fix, var):
        q = dd[(dd.tissue == tis) & dd[var].isin(["pre", "post"])]
        if fix: q = q[q[fix[0]] == fix[1]]
        D = []
        for s, g in q.groupby("subject"):
            a = g.gsm[g[var] == "pre"]; b_ = g.gsm[g[var] == "post"]
            if len(a) == 1 and len(b_) == 1: D.append(E[b_.iloc[0]] - E[a.iloc[0]])
        return pd.concat(D, axis=1).median(1), len(D)
    out = {}
    for name, args in (("clamp_baseline", ("muscle", ("tzd", "pre"), "clamp")), ("tzd_muscle", ("muscle", ("clamp", "pre"), "tzd")), ("tzd_adipose", ("adipose", None, "tzd"))):
        out[name], out[name + "_n"] = paired(*args)
    return out


def t_gse13070_ir():
    """GSE13070_IR_vs_sensitive_log2.tsv (Fig. 3h-i, 5c) y observed_signature_for_model.tsv:
    media de resistentes menos media de sensibles, basal, musculo (51 vs 18) y adiposo (35 vs 6)."""
    E, d = gse13070()
    mus = d[(d.tissue == "muscle") & (d.tzd == "pre") & (d.clamp == "pre")].drop_duplicates("subject")
    adi = d[(d.tissue == "adipose") & (d.tzd == "pre")].drop_duplicates("subject")
    out = {}
    for lab, q in (("IR muscle (GSE13070, 51 vs 18)", mus), ("IR adipose (GSE13070, 35 vs 6)", adi)):
        r = q.gsm[q.sensitivity == "resistant"]; s = q.gsm[q.sensitivity == "sensitive"]
        print(f"    {lab}: {len(r)} resistentes, {len(s)} sensibles")
        out[lab] = E[list(r)].mean(1) - E[list(s)].mean(1)
    D = pd.DataFrame(out)
    compare("GSE13070_IR_vs_sensitive_log2.tsv", D, note="todos los genes; la referencia solo lleva los 6,108 comunes a los 12 contrastes")
    genes = ["PER3", "DBP", "TEF", "NR1D1", "NR1D2", "PER2", "CIART", "BMAL1"]   # brazo de salida, Fig. 5c
    O = D["IR adipose (GSE13070, 35 vs 6)"].reindex(genes).rename("IR_minus_sensitive_log2").to_frame()
    compare("observed_signature_for_model.tsv", O)
    return D


def gse129843_by_time():
    X = pd.read_pickle(f"{NEW}/GSE129843_log2cpm.pkl"); md = pd.read_pickle(f"{NEW}/GSE129843_meta.pkl").set_index("gsm")
    md["ind"] = md.title.str.extract(r"Participant: (\d+)")[0]
    md["cond"] = np.where(md.c0.str.contains("Time restricted"), "TRF", "EXF")
    md["t"] = md.c1.str.replace("time point: ", "").str.split(":").str[0].astype(float)
    return X, md.reindex(X.columns)


def gse108539_by_time():
    p = mod("b2", "02b_phase_template.py"); md, _ = p.read_series("GSE108539"); X = p.load_counts("GSE108539", md)
    md["clock_h"] = md["time"].astype(str).str.split(":").str[0].astype(float); md["subject"] = md.title.str.split("_").str[0]
    return X, md.set_index("gsm").reindex(X.columns)


def morning(X, md, t0, t1, ind="ind", tcol="t"):
    """Media sobre personas de (expresion a t1 - expresion a t0), personas con las dos horas."""
    D = []
    for s, q in md.groupby(ind):
        a = q.index[q[tcol] == t0]; b = q.index[q[tcol] == t1]
        if len(a) == 1 and len(b) == 1: D.append(X[b[0]] - X[a[0]])
    return pd.concat(D, axis=1).mean(1), len(D)


def t_time_calibration():
    """time_calibration_vs_clamp.tsv (Fig. 3b): cambio medio por gen entre dos horas de la manana
    en las cohortes con hora, frente al cambio con el clamp de 3 h."""
    X, md = gse129843_by_time()
    exf, n1 = morning(X, md[md.cond == "EXF"], 7, 11); trf, n2 = morning(X, md[md.cond == "TRF"], 7, 11)
    X8, md8 = gse108539_by_time(); m8, n3 = morning(X8, md8, 8, 12, ind="subject", tcol="clock_h")
    E, dd = gse13070(); cc = clamp_contrasts(E, dd)
    genes = ["BMAL1", "CDKN1A", "CIART", "CRY1", "DBP", "DDIT4", "FBXO32", "FKBP5", "GADD45A", "HK2", "HLF", "KLF15", "NR1D1", "NR1D2",
             "PDK4", "PER1", "PER3", "SLC2A4", "TEF", "TRIM63", "TSC22D3", "TXNIP"]   # Fig. 3b y Supp. Table 7
    D = pd.DataFrame({f"GSE129843 EXF 11h-07h (n={n1})": exf, f"GSE129843 TRF 11h-07h (n={n2})": trf,
                      f"GSE108539 12h-08h (n={n3})": m8, f"GSE13070 clamp 3h (n={cc['clamp_baseline_n']})": cc["clamp_baseline"]}).reindex(genes).round(2)
    ref = read_ref("time_calibration_vs_clamp.tsv")
    if list(D.columns) != list(ref.columns): print(f"    (aviso) rotulos distintos de la referencia: {list(D.columns)}")
    D.columns = list(ref.columns) if D.shape[1] == ref.shape[1] else D.columns
    compare("time_calibration_vs_clamp.tsv", D, note="clamp: mediana de log2FC por pareja (clamp_contrasts)")
    return exf, trf, m8, cc


def t_template():
    """template_2harm.tsv: cosinor de dos armonicos por gen en GSE108539 (datos estandarizados por
    gen), con r2, coeficientes, amplitudes y asimetria."""
    p = mod("b2", "02b_phase_template.py")
    X, md = gse108539_by_time()
    t = md.clock_h.values * 2 * np.pi / 24
    Z = ((X.T - X.mean(1)) / (X.std(1) + 1e-9)).T
    D = p.design(t, [24.0, 12.0]); B, *_ = np.linalg.lstsq(D, Z.values.T, rcond=None)
    fit = (D @ B).T; r2 = 1 - ((Z.values - fit) ** 2).sum(1) / (((Z.values - Z.values.mean(1, keepdims=True)) ** 2).sum(1) + 1e-12)
    a1, b1, a2, b2 = B[1], B[2], B[3], B[4]
    T = pd.DataFrame(dict(gene=Z.index, r2=r2, a1=a1, b1=b1, a2=a2, b2=b2, amp1=np.hypot(a1, b1), amp2=np.hypot(a2, b2)))
    T["asym"] = T.amp2 / T.amp1
    T = T.sort_values("r2", ascending=False)
    compare("template_2harm.tsv", T, key=["gene"], note="GSE108539 desde las cuentas de NCBI; 57 muestras")
    T.to_csv(f"{OUT}/template_2harm.tsv", sep="\t", index=False)   # 04, 05 y 17 la leen de results/


def gse17992():
    """GSE17992 (Agilent 4x44K, dos colores): log2(dex/control) por sujeto y tejido, promediando
    las parejas con intercambio de colorante. Devuelve (musculo genes x sujeto, adiposo)."""
    b = mod("b0", "00b_build_matrices.py"); p, e = b.read_series("GSE17992")
    e = e.mask(e.abs() > 50)   # GSM446418 (BS) lleva -1e5 como valor faltante en 9,415 sondas
    ann = pd.read_csv(f"{RAW}/GPL4133_annot.txt", sep="\t", dtype=str, low_memory=False)
    col = next(c for c in ann.columns if c.upper() in ("GENE_SYMBOL", "GENE SYMBOL"))
    symb = ann.set_index("ID")[col].dropna()
    p["subj"] = p.title.str.extract(r"(\b[A-Z]{2})( rep)?$")[0]; p["rep"] = p.title.str.endswith("rep")
    p["tissue"] = np.where(p.c2.str.contains("muscle"), "muscle", "adipose")
    # orientacion: valor = log(ch1/ch2); dex - control = +valor si ch1 es dex, -valor si ch1 es control
    ch1dex = p.c0.str.contains("dexa")
    out = {}
    for tis in ("muscle", "adipose"):
        q = p[p.tissue == tis]; D = {}
        for s, qq in q.groupby("subj"):
            v = [e[g] * (1 if ch1dex[i] else -1) for i, g in zip(qq.index, qq.gsm)]
            D[s] = pd.concat(v, axis=1).mean(1)
        D = pd.DataFrame(D)
        D = D.loc[D.index.intersection(symb.index)]; D.index = symb[D.index].values
        D = D.groupby(level=0).mean()
        syn = pd.read_csv(f"{NEW}/ncbi_synonym_map.tsv", sep="\t", index_col=0)["current"].to_dict(); have = set(D.index)
        D.index = [syn[g] if (g in syn and syn[g] not in have) else g for g in D.index]   # simbolo vigente (ARNTL -> BMAL1)
        out[tis] = D.groupby(level=0).mean()
    return out, p


def t_gse17992():
    """GSE17992_dex_log2ratio_muscle.tsv (Fig. 3f)."""
    out, p = gse17992()
    D = out["muscle"]
    ref = read_ref("GSE17992_dex_log2ratio_muscle.tsv"); common = D.index.intersection(ref.index)
    # escala: Agilent deposita log10 o log2; se elige la que mejor ajusta a la referencia
    with np.errstate(divide="ignore", invalid="ignore"):          # 0/0 donde las dos tablas valen cero: se descartan justo despues
        ratio = (ref.loc[common].values / D.loc[common, ref.columns].values).ravel()
    k = np.median(ratio[np.isfinite(ratio)])
    print(f"    razon mediana referencia/nueva = {k:.3f} (log10 -> log2 seria 3.32)")
    compare("GSE17992_dex_log2ratio_muscle.tsv", D, note=f"razon mediana ref/nueva {k:.2f}; {p.tissue.value_counts().to_dict()}")
    return out


def t_signatures(IR, exf, m8, dex, cc):
    """signatures_reference.tsv: las doce columnas, genes medidos en todas. "T2D myotubes" se define
    como la mediana de T2D menos la mediana de NGT en el brazo control de GSE182117 (todas las
    horas); es una definicion deducida, el manuscrito no la da."""
    c23 = mod("c23", "23_build_signatures_table.py")
    col = lambda c: cc[c]
    X = pd.read_pickle(f"{NEW}/GSE182117_log2cpm.pkl"); md = pd.read_pickle(f"{NEW}/GSE182117_meta.pkl").set_index("gsm").reindex(X.columns)
    ctl = md.c1.str.contains("CTL"); t2d = md.c0.str.contains("T2D")
    myo = X.loc[:, ctl & t2d].median(1) - X.loc[:, ctl & ~t2d].median(1)
    S = pd.DataFrame({IR.columns[0]: IR.iloc[:, 0], IR.columns[1]: IR.iloc[:, 1],
                      "clamp 3h muscle (GSE13070)": col("clamp_baseline"), "TZD 3mo muscle (GSE13070)": col("tzd_muscle"),
                      "TZD 3mo adipose (GSE13070)": col("tzd_adipose"), "morning 07->11h muscle (GSE129843 EXF)": exf,
                      "morning 08->12h muscle (GSE108539)": m8, "dex 4d muscle (GSE17992)": dex["muscle"].mean(1),
                      "dex 4d adipose (GSE17992, n=4)": dex["adipose"].mean(1), "hydrocortisone adipose (GSE148640)": c23.col_hydrocortisone(),
                      "Cushing adipose (GSE66446)": c23.col_cushing(), "T2D myotubes (GSE182117, cell-autonomous)": myo})
    print('    genes con valor por columna:', S.notna().sum().to_dict())
    S = S.dropna()
    compare("signatures_reference.tsv", S, note="todas las columnas regeneradas; TZD/clamp desde 05")


# --- estimadores de amplitud por muestra -------------------------------------------------------
def template800():
    T = pd.read_csv(f"{NEW}/template_2harm.tsv", sep="\t").set_index("gene")
    return T, list(T.sort_values("r2", ascending=False).index[:800])


def amp_free(E, genes=None):
    """phase.amplitude con ordenamiento libre sobre los genes del reloj (restringido y luego centrado por gen)."""
    g = [x for x in (genes or ph.CLK) if x in E.index]
    Z = E.loc[g]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T
    ang = ph.unsup(Z.values)
    return ph.amplitude(E, ang, g), ang


def t_amplitude_two_cohorts():
    """amplitude_phase_distributions_two_cohorts.tsv (dispersion, Supp. Fig. 2): por muestra, fase y
    amplitud con la plantilla de 800 genes (amp_800) y con los genes del reloj (amp_clock), y la
    amplitud libre (amp_free), en GSE13070 musculo basal y GSE182120."""
    p4 = mod("p4", "04_perturbations.py"); T, g800 = template800()
    E, d = gse13070(); mus = d[(d.tissue == "muscle") & (d.tzd == "pre") & (d.clamp == "pre")].drop_duplicates("subject").set_index("gsm")
    Ee = pd.read_pickle(f"{NEW}/GSE13070_entrez.pkl")   # anotacion por Entrez (Methods); 05 usa el simbolo de plataforma
    X1 = Ee[list(mus.index)]
    X2 = pd.read_pickle(f"{NEW}/GSE182120_log2.pkl"); md2 = pd.read_pickle(f"{NEW}/GSE182120_meta.pkl").set_index("gsm").reindex(X2.columns)
    if CENTER: X1, X2 = center(X1), center(X2)   # centrado por muestra antes de cualquier score (Methods); T2D_CENTER=0 lo quita
    rows = []
    for X, grp, coh, M, lab in ((X1, mus.sensitivity, "GSE13070 músculo basal", None, None),
                                (X2, pd.Series(np.where(md2.disease == "T2D", "resistant", "sensitive"), index=X2.columns), "GSE182120 músculo", md2.mvalue, md2.lab)):
        a8 = p4.phase_amplitude(X, T, g800).set_index("sample"); ac = p4.phase_amplitude(X, T, [g for g in ph.CLK if g in T.index]).set_index("sample")
        af, _ = amp_free(X)
        D = pd.DataFrame(dict(amp_800=a8.amplitude, phase_800=a8.phase, amp_clock=ac.amplitude, phase_clock=ac.phase, amp_free=af,
                              group=grp.reindex(X.columns).values, cohort=coh, M=M.reindex(X.columns).values if M is not None else "",
                              lab=lab.reindex(X.columns).values if lab is not None else ""))
        rows.append(D)
    D = pd.concat(rows)
    compare("amplitude_phase_distributions_two_cohorts.tsv", D, note=f"amp_free = phase.amplitude sobre genes del reloj; centrado={CENTER}, matriz Entrez")


def cosinor_amp(y, t, per=24.0):
    w = 2 * np.pi / per; D = np.column_stack([np.ones(len(t)), np.cos(w * t), np.sin(w * t)])
    b, *_ = np.linalg.lstsq(D, y, rcond=None); return np.hypot(b[1], b[2]), 1 - ((y - D @ b) ** 2).sum() / (((y - y.mean()) ** 2).sum() + 1e-12)


def t_gse182117_inferred():
    """GSE182117_inferred_free_vs_measured.tsv (Fig. 1f): por donante y tratamiento, amplitud
    medida (cosinor sobre los 8 puntos, media de los genes del reloj) y amplitudes inferidas."""
    p4 = mod("p4", "04_perturbations.py"); T, g800 = template800()
    X = pd.read_pickle(f"{NEW}/GSE182117_log2cpm.pkl"); md = pd.read_pickle(f"{NEW}/GSE182117_meta.pkl").set_index("gsm").reindex(X.columns)
    md["ind"] = md.c4.str.replace("individual: ", ""); md["trt"] = md.c1.str.replace("treatment: ", ""); md["disease"] = md.c0.str.replace("disease: ", "")
    md["t"] = md.c3.str.replace("time: ", "").str.replace("h", "").astype(float)
    Xc = center(X); g = [x for x in ph.CLK if x in Xc.index]; core = [x for x in SETS["Core oscillator"]["genes"] if resolve(x, Xc.index)]
    a8 = p4.phase_amplitude(Xc, T, g800).set_index("sample"); af, _ = amp_free(Xc)
    rows = []
    for (ind, trt), q in md.groupby(["ind", "trt"]):
        cols = list(q.index); t = q.t.values
        ca = [cosinor_amp(np.asarray(Xc.loc[gg, cols], float), t) for gg in g]
        cc = [cosinor_amp(np.asarray(Xc.loc[resolve(gg, Xc.index), cols], float), t)[0] for gg in core]
        rows.append(dict(ind=ind, trt=trt, disease=q.disease.iloc[0], amp_inferred=a8.amplitude[cols].mean(), fit=a8.fit[cols].mean(),
                         amp_mean_log2=np.mean([c[0] for c in ca]), r2_mean=np.mean([c[1] for c in ca]), amp_core=np.mean(cc), amp_free=af[cols].mean()))
    compare("GSE182117_inferred_free_vs_measured.tsv", pd.DataFrame(rows), key=["ind", "trt"])


def t_gse129843_amp():
    """GSE129843_amplitude_per_person.tsv (Fig. 1g): por persona y condicion, amplitud medida
    (cosinor sobre las 6 biopsias) y las inferidas."""
    p4 = mod("p4", "04_perturbations.py"); T, g800 = template800()
    X, md = gse129843_by_time(); Xc = center(X); g = [x for x in ph.CLK if x in Xc.index]
    a8 = p4.phase_amplitude(Xc, T, g800).set_index("sample"); ac = p4.phase_amplitude(Xc, T, [x for x in g if x in T.index]).set_index("sample")
    af, _ = amp_free(Xc)
    rows = []
    for (ind, cond), q in md.groupby(["ind", "cond"]):
        cols = list(q.index); t = q.t.values
        if len(cols) < 5: continue
        rows.append(dict(ind=ind, cond=cond, n=len(cols), amp_measured=np.mean([cosinor_amp(np.asarray(Xc.loc[gg, cols], float), t)[0] for gg in g]),
                         amp_inf_800=a8.amplitude[cols].mean(), amp_inf_clk=ac.amplitude[cols].mean(), amp_free=af[cols].mean()))
    D = pd.DataFrame(rows); D["ind"] = D.ind.astype(int)
    compare("GSE129843_amplitude_per_person.tsv", D, key=["ind", "cond"])


def t_symbol_audit():
    """symbol_audit_renamed_genes.tsv (Supp. Fig. 3b-c): simbolos de plataforma que ya no son el
    vigente, por plataforma, y si el gen esta en los conjuntos del articulo."""
    syn = pd.read_csv(f"{NEW}/ncbi_synonym_map.tsv", sep="\t", index_col=0)["current"].to_dict()
    ann = pd.read_csv(f"{RAW}/Human_GRCh38_p13_annot_tsv.gz", sep="\t", usecols=["GeneID", "Symbol"], dtype=str); cur = set(ann.Symbol)
    insets = {g for s in SETS.values() for g in s["genes"]} | set(ph.CLK)
    rows = []
    def add(platform, old_syms):
        for o in sorted({x for x in old_syms if isinstance(x, str)}):
            if o in cur: continue
            c = syn.get(o)
            if c and c != o: rows.append(dict(platform=platform, old=o, current=c, in_sets=c in insets))
    b = mod("b0", "00b_build_matrices.py")
    _, e = b.read_series("GSE95674"); add("GPL10904 Illumina HT-12 v4", e.index)
    for plat, fn, col in (("GPL570 Affymetrix U133 Plus 2.0", "GPL570_annot.txt", "Gene Symbol"), ("GPL17586 HTA-2.0", "GPL17586_annot.txt", None),
                          ("GPL4133 Agilent 4x44K", "GPL4133_annot.txt", "GENE_SYMBOL"), ("GPL15331 Agilent", "GPL15331_annot.txt", "GENE_SYMBOL")):
        if not os.path.exists(f"{RAW}/{fn}"): print(f"    (sin {fn})"); continue
        a = pd.read_csv(f"{RAW}/{fn}", sep="\t", dtype=str, low_memory=False)
        if col is None:
            s = a.gene_assignment.dropna().str.split(" // ").str[1]
        else:
            c2 = next((c for c in a.columns if c.upper() == col.upper().replace("_", " ") or c.upper() == col.upper()), None)
            s = a[c2].dropna().str.split(" /// ").str[0] if c2 else pd.Series(dtype=str)
        add(plat, s)
    D = pd.DataFrame(rows)
    ref = read_ref("symbol_audit_renamed_genes.tsv"); k = ["platform", "old"]
    common = len(ref.set_index(k).index.intersection(D.set_index(k).index)) if len(D) else 0
    print(f"  [symbol_audit_renamed_genes.tsv] ref {len(ref)} / nueva {len(D)} / comunes {common}")
    CHECK.append(dict(tabla="symbol_audit_renamed_genes.tsv", filas_ref=len(ref), filas_nueva=len(D), filas_comunes=common, columnas_comparadas=0,
                      r_minimo=np.nan, detalle=f"{common} de {len(ref)} pares (plataforma, simbolo antiguo) de la referencia estan en la nueva",
                      nota="la referencia es un subconjunto de 428 pares; criterio de seleccion no documentado"))
    D.to_csv(f"{NEW}/symbol_audit_renamed_genes.tsv", sep="\t", index=False)


# orden: la plantilla se regenera antes de las tablas de amplitud que la usan
STEPS = {"template": t_template, "within_donor": t_within_donor, "gcreceptor": t_gcreceptor_adiposity, "rank": t_rank_vs_nr3c1,
         "synonyms": t_synonym_map, "coverage": t_msigdb_coverage, "gse26637": t_gse26637, "gse148640": t_gse148640_stats,
         "gras": t_gras, "foxo3": t_foxo3, "symbol_audit": t_symbol_audit, "audit": t_audit_global_axis,
         "amplitude": t_amplitude_two_cohorts, "gse182117": t_gse182117_inferred, "gse129843": t_gse129843_amp,
         "enrichment": t_enrichment}


def run_signature_group():
    """Las tablas de firmas comparten calculos: IR, calibracion, dexametasona y la tabla de 12 columnas."""
    IR = t_gse13070_ir(); exf, trf, m8, cc = t_time_calibration(); dex = t_gse17992(); t_signatures(IR, exf, m8, dex, cc)


if __name__ == "__main__":
    want = sys.argv[1:] or (list(STEPS) + ["signatures"])
    for k in want:
        print(f"\n=== {k}")
        try:
            (run_signature_group if k == "signatures" else STEPS[k])()
        except Exception as ex:
            import traceback; traceback.print_exc(); print(f"  FALLO {k}: {type(ex).__name__}: {ex}")
            CHECK.append(dict(tabla=k, filas_ref=np.nan, filas_nueva=np.nan, filas_comunes=np.nan, columnas_comparadas=0, r_minimo=np.nan, detalle="", nota=f"FALLO: {ex}"))
    pd.DataFrame(CHECK).to_csv(f"{NEW}/reference_tables_check.tsv", sep="\t", index=False)
    print(f"\nresumen -> {NEW}/reference_tables_check.tsv")
    failed = [c["tabla"] for c in CHECK if str(c.get("nota", "")).startswith("FALLO")]
    if failed: sys.exit(f"productores fallidos: {', '.join(failed)}")
