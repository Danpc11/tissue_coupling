#!/usr/bin/env python3
"""GSE13070: el reloj y el programa de ayuno DENTRO de cada persona, bajo insulina y bajo TZD.

Pregunta. El modelo del articulo revisado dice que la senalizacion de insulina concentra y entrena
los relojes perifericos y que en resistencia a la insulina domina el programa de ayuno
glucocorticoide/FOXO3 (TXNIP, PDK4, ANGPTL4, FBXO32...). Si es asi, una infusion aguda de insulina
debe (i) bajar ese programa en el musculo de la misma persona y (ii) la caida por persona debe
acompanar la ganancia de amplitud o de alineamiento de fase de esa persona. Es la prueba dentro
del individuo que falta al articulo.

Diseno real de GSE13070 (364 biopsias, 72 sujetos; en el titulo de cada muestra):
  musculo:  preTZD/postTZD x preClamp/postClamp;  adiposo: preTZD/postTZD, SIN clamp.
  grupos:   NRd (Rd normal, sensibles); IRd, IRd+ (respondedor a TZD), IRd- (no respondedor).
Contrastes pareados que se evaluan, cada uno estimado JUNTO (restringir primero, centrar despues):
  clamp_baseline   musculo preTZD: preClamp -> postClamp        (PRIMARIO)
  clamp_postTZD    musculo postTZD: preClamp -> postClamp       (replicacion interna)
  tzd_muscle       musculo preClamp: preTZD -> postTZD
  tzd_adipose      adiposo: preTZD -> postTZD

Estimadores. Si existe results/template_2harm.tsv se usa la proyeccion de 04_perturbations.py con
el argmax sobre la proyeccion CON SIGNO (en 04 es sobre |proj|, que empata fase y antifase; ver
phase_amplitude_signed). Si no existe, se ancla la fase con la plantilla de genes del reloj que
lleva docs/phase.js (20 genes, fases de la cohorte con horas conocidas) mediante
phase.order_anchored, y la amplitud se estima con ordenamiento libre (phase.unsup + phase.amplitude),
que es la combinacion que el articulo declara para cada cantidad.

Bloques de resultado por contraste:
  A. R de fase y amplitud mediana pre/post, por grupo.          B. cambio por persona de los
  scores de conjunto y de sus genes.                              C. acoplamiento dentro del
  individuo: rho(cambio de score, cambio de amplitud | ganancia de alineamiento de fase), con
  P por permutacion del emparejamiento. gain_i = cos(phi_post,i - mu_post) - cos(phi_pre,i - mu_pre).
  D. control positivo: media circular de los desplazamientos de fase por persona; en el clamp
  debe aproximarse a su duracion (la biopsia post se toma al final de la infusion).

Salidas: results/clamp_samples.tsv, clamp_groups.tsv, clamp_scores_change.tsv,
         clamp_genes_change.tsv, clamp_within_person.tsv, clamp_per_person_<contraste>.tsv
Uso:     python python/analyses/05_clamp_within_person.py
"""
import os
import re
import sys
import glob
import gzip
import warnings
import importlib.util

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../lib"))
from gene_sets import SETS, resolve  # noqa: E402
import phase as ph  # noqa: E402

RAW = os.environ.get("T2D_RAW", "data/raw")
OUT = os.environ.get("T2D_OUT", "results")
ACC = "GSE13070"
NPERM = int(os.environ.get("T2D_PERM", 5000))
NTEMPLATE = int(os.environ.get("T2D_CLAMP_TEMPLATE", 800))

# Programa de ayuno / dianas FOXO-glucocorticoide en musculo, pre-especificado. Si el entorno tiene
# el regulon FOXO3 de DoRothEA usado en la Fig. 5, sustituirlo aqui y decirlo en Metodos.
FASTING = ["TXNIP", "PDK4", "ANGPTL4", "FBXO32", "TRIM63", "CDKN1A", "GADD45A", "SESN1",
           "BNIP3L", "KLF15", "DDIT4"]
UPTAKE = ["HK2", "SLC2A4"]
SCORE_SETS = {
    "Glucocorticoid": SETS["Glucocorticoid"]["genes"],
    "Insulin core": SETS["Insulin core"]["genes"],
    "Clock output": SETS["Clock output"]["genes"],
    "Core oscillator": SETS["Core oscillator"]["genes"],
    "Oxidative": SETS["Oxidative"]["genes"],
    "Fasting/FOXO3": FASTING,
    "Glucose uptake": UPTAKE,
}
# Plantilla de genes del reloj de docs/phase.js: (fase del maximo en h, amplitud relativa)
CLOCK_TEMPLATE = {
    "BMAL1": (13.6, 1.07), "NPAS2": (14.3, 1.10), "CLOCK": (15.2, 0.62), "PER1": (9.2, 0.88),
    "PER2": (4.0, 0.96), "PER3": (12.0, 0.92), "NR1D1": (20.2, 0.95), "NR1D2": (12.0, 0.72),
    "DBP": (23.7, 1.09), "TEF": (12.0, 0.83), "HLF": (12.0, 0.64), "CIART": (12.0, 0.92),
    "CRY1": (0.0, 0.84), "CRY2": (12.0, 0.66), "BHLHE40": (0.0, 1.11), "BHLHE41": (1.0, 1.11),
    "NFIL3": (2.0, 0.90), "RORA": (6.0, 0.55)}
RENAME = {"ARNTL": "BMAL1", "C1ORF51": "CIART", "C1orf51": "CIART"}


# --------------------------------------------------------------------------- lectura
def read_series(acc):
    meta, tab, intab, nchar = {}, [], False, 0
    with gzip.open(f"{RAW}/{acc}_series_matrix.txt.gz", "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!series_matrix_table_begin"):
                intab = True
                continue
            if l.startswith("!series_matrix_table_end"):
                break
            if intab:
                tab.append(l.rstrip("\n").split("\t"))
                continue
            for k, t in [("gsm", "!Sample_geo_accession"), ("title", "!Sample_title")]:
                if l.startswith(t) and k not in meta:
                    meta[k] = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
    p = pd.DataFrame(meta)
    e = pd.DataFrame(tab[1:], columns=[c.strip('"') for c in tab[0]]).set_index(tab[0][0].strip('"'))
    e.index = [i.strip('"') for i in e.index]
    return p, e.apply(pd.to_numeric, errors="coerce").dropna()


def annotate_gpl570(e):
    """Sonda -> simbolo (columna Gene Symbol), media entre sondas. Las sondas multi-gen
    ("MIR6883 /// PER1", "NR1D1 /// THRA") se asignan a cada simbolo sin sonda propia; sin esto
    PER1 y NR1D1 se declaran ausentes en GPL570."""
    files = sorted(glob.glob(f"{RAW}/GPL570*"))
    if not files:
        sys.exit(f"Falta la anotacion GPL570 en {RAW}.")
    fn = files[0]
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt", errors="ignore") as fh:
        lines = [l for l in fh if not l.startswith(("#", "^", "!"))]
    a = pd.read_csv(pd.io.common.StringIO("".join(lines)), sep="\t", dtype=str, low_memory=False)
    sym = next(c for c in a.columns if c.lower().replace("_", " ") == "gene symbol")
    idc = next(c for c in a.columns if c.upper() == "ID")
    m = a.set_index(idc)[sym].dropna()
    single = m[~m.str.contains("///")]
    e1 = e.loc[e.index.intersection(single.index)]
    e1.index = single[e1.index].values
    G = e1.groupby(level=0).mean()
    rows = {}
    for probe, syms in m[m.str.contains("///")].items():
        if probe in e.index:
            for g in [x.strip() for x in syms.split("///")]:
                if g and g not in G.index and not g.startswith(("MIR", "LOC", "SNOR")):
                    rows.setdefault(g, []).append(e.loc[probe].values)
    if rows:
        G = pd.concat([G, pd.DataFrame({g: np.mean(v, axis=0) for g, v in rows.items()},
                                       index=e.columns).T])
    G.index = [RENAME.get(i, i) for i in G.index]
    return G.groupby(level=0).mean()


def load():
    cache = f"{OUT}/{ACC}_genes.pkl"
    if os.path.exists(cache):
        E, p = pd.read_pickle(cache), pd.read_pickle(f"{OUT}/{ACC}_meta.pkl")
    else:
        p, e = read_series(ACC)
        E = annotate_gpl570(e)
        E = np.log2(E.clip(lower=1)) if E.values.max() > 60 else E
        E.to_pickle(cache); p.to_pickle(f"{OUT}/{ACC}_meta.pkl")
    E.index = [RENAME.get(i, i) for i in E.index]
    E = E.groupby(level=0).mean()
    # diseno desde el titulo: muscle_preTZD_preClamp_IRd+_s1 / adipose_preTZD_IRd+_s1
    t = p.title.str.replace('"', "")
    d = pd.DataFrame(dict(gsm=p.gsm, title=t))
    d["tissue"] = t.str.split("_").str[0]
    d["tzd"] = t.str.extract(r"(pre|post)TZD")[0]
    d["clamp"] = t.str.extract(r"(pre|post)Clamp")[0]
    d["group"] = t.str.extract(r"_(NRd|IRd[+-]?)_")[0]
    d["sensitivity"] = np.where(d.group == "NRd", "sensitive", "resistant")
    d["subject"] = t.str.extract(r"_s(\d+)$")[0]
    return E, d[d.gsm.isin(E.columns)]


# --------------------------------------------------------------------------- estimadores
def phase_amplitude_signed(E, template, genes, pert):
    """Como pert.phase_amplitude pero con argmax sobre la proyeccion CON SIGNO. En 04 se toma
    sobre |proj|: con un armonico S(t+pi) = -S(t) y el empate fase/antifase lo rompe el ruido
    (52 % de antifases en simulacion sin 12 h real; 0-1 % con 12 h real)."""
    g = [x for x in genes if x in E.index and x in template.index]
    S = pert.waveform(template, g)
    X = E.loc[g].sub(E.loc[g].mean(1), axis=0).values
    proj = S @ X
    idx = np.argmax(proj, axis=0); idx_abs = np.argmax(np.abs(proj), axis=0)
    fit = S @ (X / np.linalg.norm(X, axis=0, keepdims=True))
    dh = np.abs(np.angle(np.exp(1j * (pert.GRID[idx] - pert.GRID[idx_abs])))) * 24 / (2 * np.pi)
    return pd.DataFrame(dict(sample=E.columns, phase=pert.GRID[idx] * 24 / (2 * np.pi),
                             amplitude=[proj[idx[i], i] for i in range(X.shape[1])],
                             fit=[fit[idx[i], i] for i in range(X.shape[1])],
                             antiphase_flag=(dh > 9).astype(int)))


def estimate(E, estimator, basis_cols=None):
    """Fase (h), amplitud y ajuste por muestra; E ya restringida al conjunto que se compara."""
    kind, obj = estimator
    if kind == "2harm":
        template, genes, pert = obj
        pa = phase_amplitude_signed(E, template, genes, pert)
        return pa.set_index("sample")
    template = obj
    genes = [g for g in template.index if g in E.index]
    ang, used = ph.order_anchored(E, template, genes=genes)
    # ajuste: correlacion del perfil con la onda en la fase asignada
    Z = E.loc[used]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T.values
    T = np.array([template.amp[used].values * np.cos(a - template.phase_h[used].values * 2 * np.pi / 24)
                  for a in ang.values])
    T = (T - T.mean(1, keepdims=True)) / (T.std(1, keepdims=True) + 1e-9)
    O = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
    fit = np.array([T[i] @ O[:, i] / len(used) for i in range(O.shape[1])])
    # amplitud por ordenamiento libre. La base de la onda se ajusta SOLO en las muestras basales
    # (basis_cols): si se ajusta con pre y post juntas, el ciclo dominante que recupera unsup es el
    # efecto de la perturbacion y no el reloj (BMAL1-DBP a 6.5 h en el clamp, fuera de 8-16 h).
    # La estandarizacion por gen sigue siendo la del conjunto completo (restringir, luego centrar).
    Zall = ((E.loc[used].T - E.loc[used].mean(1)) / (E.loc[used].std(1) + 1e-9)).T
    bcols = basis_cols if basis_cols is not None else list(E.columns)
    free = ph.unsup(Zall[bcols].values)
    ok, sep = ph.is_clock(E[bcols], free, genes=used)
    D = np.column_stack([np.cos(free), np.sin(free)])
    B, *_ = np.linalg.lstsq(D, Zall[bcols].values.T, rcond=None)
    Qb, _ = np.linalg.qr(B.T)
    amp = np.sqrt(((Zall.values.T @ Qb) ** 2).sum(1))
    return pd.DataFrame(dict(phase=ang.values * 24 / (2 * np.pi), amplitude=amp,
                             fit=fit, antiphase_flag=0, clock_ok=ok, bmal1_dbp_sep=sep),
                        index=E.columns)


# --------------------------------------------------------------------------- estadisticos
def circ_R(h):
    return float(np.abs(np.mean(np.exp(1j * np.asarray(h) * 2 * np.pi / 24))))


def circ_mean_rad(h):
    return float(np.angle(np.mean(np.exp(1j * np.asarray(h) * 2 * np.pi / 24))))


def perm_pvalue(stat_fn, pre, post, rng, nperm=NPERM):
    obs = stat_fn(pre, post); null = []
    for _ in range(nperm):
        flip = rng.random(len(pre)) < 0.5
        null.append(stat_fn(np.where(flip, post, pre), np.where(flip, pre, post)))
    null = np.array(null)
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (nperm + 1)


def perm_spearman(x, y, rng, nperm=NPERM):
    r = stats.spearmanr(x, y).statistic
    null = np.array([stats.spearmanr(x, rng.permutation(y)).statistic for _ in range(nperm)])
    return r, (np.sum(np.abs(null) >= abs(r)) + 1) / (nperm + 1)


def set_scores(E, genes):
    g = [r for r in (resolve(x, E.index) for x in genes) if r is not None]
    if len(g) < 2:
        return None, g
    Z = E.loc[g]; Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T
    return Z.mean(0), g


def wilcox(x):
    x = np.asarray(x)
    return stats.wilcoxon(x).pvalue if np.any(x != 0) and len(x) >= 6 else np.nan


# --------------------------------------------------------------------------- un contraste
def run_contrast(E, pre, post, label, estimator, rng, rows):
    """pre/post: DataFrames indexados por sujeto con columnas gsm y group, mismos sujetos."""
    ids = list(pre.index)
    gpre, gpost = list(pre.gsm), list(post.gsm)
    grp = pre.sensitivity
    Et = E[gpre + gpost]                                   # restringir primero
    pa = estimate(Et, estimator, basis_cols=gpre)          # centrar dentro; base de onda en pre
    pa["subject"] = ids + ids; pa["condition"] = ["pre"] * len(ids) + ["post"] * len(ids)
    pa["contrast"] = label; pa["sensitivity"] = list(grp) * 2
    rows["samples"].append(pa.reset_index(names="gsm"))
    P, Q = pa.iloc[:len(ids)], pa.iloc[len(ids):]
    if "clock_ok" in pa:
        print(f"  reloj recuperado por ordenamiento libre en las basales: BMAL1-DBP {pa.bmal1_dbp_sep.iloc[0]:.1f} h "
              f"({'OK' if pa.clock_ok.iloc[0] else 'FUERA de 8-16 h'})")
    if pa.antiphase_flag.mean() > 0:
        print(f"  antifase |proj| vs proj: {pa.antiphase_flag.mean()*100:.0f} % de muestras")

    # A ----------------------------------------------------------------------------------
    for g in ["all", "sensitive", "resistant"]:
        m = np.ones(len(ids), bool) if g == "all" else (grp == g).values
        if m.sum() < 4:
            continue
        for c, sub in (("pre", P[m]), ("post", Q[m])):
            rows["groups"].append(dict(contrast=label, group=g, condition=c, n=int(m.sum()),
                                       R_phase=circ_R(sub.phase), median_amplitude=float(sub.amplitude.median())))
    dR, pR = perm_pvalue(lambda a, b: circ_R(b) - circ_R(a), P.phase.values, Q.phase.values, rng)
    dA = (Q.amplitude.values - P.amplitude.values)
    dphi = Q.phase.values - P.phase.values
    rot = circ_mean_rad(dphi) * 24 / (2 * np.pi)
    hist = np.histogram(((dphi + 12) % 24) - 12, bins=np.arange(-12, 13, 3))[0]
    print(f"  desplazamientos por persona (h, bins de 3 desde -12): {hist.tolist()}")
    a_pre, a_post = P.phase.values * 2 * np.pi / 24, Q.phase.values * 2 * np.pi / 24
    gain = np.cos(a_post - circ_mean_rad(Q.phase)) - np.cos(a_pre - circ_mean_rad(P.phase))
    shift = np.angle(np.exp(1j * (a_post - a_pre))) * 24 / (2 * np.pi)
    print(f"  {len(ids)} pares | R fase {circ_R(P.phase):.3f} -> {circ_R(Q.phase):.3f} "
          f"(dR {dR:+.3f}, P {pR:.3g}) | amplitud mediana {P.amplitude.median():.2f} -> "
          f"{Q.amplitude.median():.2f} (Wilcoxon P {wilcox(dA):.3g}) | rotacion media {rot:+.2f} h, "
          f"R de los desplazamientos {circ_R(dphi):.2f}")
    W = rows["within"]
    W += [dict(contrast=label, n=len(ids), quantity="phase_R", pre=circ_R(P.phase), post=circ_R(Q.phase),
               delta=dR, p=pR, test="perm sign-flip"),
          dict(contrast=label, n=len(ids), quantity="amplitude_median", pre=P.amplitude.median(),
               post=Q.amplitude.median(), delta=float(np.median(dA)), p=wilcox(dA), test="wilcoxon"),
          dict(contrast=label, n=len(ids), quantity="mean_rotation_h", pre=np.nan, post=np.nan, delta=rot,
               p=np.nan, test=f"R_shift={circ_R(dphi):.2f}; control positivo")]
    for g in ("sensitive", "resistant"):
        m = (grp == g).values
        if m.sum() >= 4:
            dRg, pRg = perm_pvalue(lambda a, b: circ_R(b) - circ_R(a), P.phase.values[m], Q.phase.values[m], rng)
            W += [dict(contrast=label, n=int(m.sum()), quantity=f"phase_R [{g}]", pre=circ_R(P.phase[m]),
                       post=circ_R(Q.phase[m]), delta=dRg, p=pRg, test="perm sign-flip"),
                  dict(contrast=label, n=int(m.sum()), quantity=f"amplitude_median [{g}]",
                       pre=P.amplitude[m].median(), post=Q.amplitude[m].median(),
                       delta=float(np.median(dA[m])), p=wilcox(dA[m]), test="wilcoxon")]
    m = (grp == "sensitive").values
    if m.sum() >= 4 and (~m).sum() >= 4:
        for lab, y in (("delta_amplitude", dA), ("alignment_gain", gain)):
            W.append(dict(contrast=label, n=len(ids), quantity=f"{lab}: sensitive vs resistant",
                          pre=float(np.median(y[m])), post=float(np.median(y[~m])),
                          delta=float(np.median(y[m]) - np.median(y[~m])),
                          p=stats.mannwhitneyu(y[m], y[~m]).pvalue, test="mann-whitney"))

    # B y C ------------------------------------------------------------------------------
    per = pd.DataFrame(dict(subject=ids, sensitivity=list(grp), delta_amplitude=dA, alignment_gain=gain,
                            phase_shift_h=shift, amp_pre=P.amplitude.values, phase_pre=P.phase.values))
    for name, gl in SCORE_SETS.items():
        sc, used = set_scores(Et, gl)
        if sc is None:
            continue
        dS = sc[gpost].values - sc[gpre].values
        per[f"d_{name}"] = dS; per[f"pre_{name}"] = sc[gpre].values
        rows["scores"].append(dict(contrast=label, score=name, n_genes=len(used), n_pairs=len(ids),
                                   median_delta=float(np.median(dS)), frac_down=float(np.mean(dS < 0)),
                                   p_wilcoxon=wilcox(dS), genes=";".join(used)))
        ms = (grp == "sensitive").values
        if ms.sum() >= 4 and (~ms).sum() >= 4:
            W.append(dict(contrast=label, n=len(ids), quantity=f"d{name}: sensitive vs resistant",
                          pre=float(np.median(dS[ms])), post=float(np.median(dS[~ms])),
                          delta=float(np.median(dS[ms]) - np.median(dS[~ms])),
                          p=stats.mannwhitneyu(dS[ms], dS[~ms]).pvalue, test="mann-whitney"))
        for g in used:
            dg = Et.loc[g, gpost].values - Et.loc[g, gpre].values
            rows["genes"].append(dict(contrast=label, score=name, gene=g, median_log2fc=float(np.median(dg)),
                                      frac_down=float(np.mean(dg < 0)), p_wilcoxon=wilcox(dg)))
        for lab, y in (("delta_amplitude", dA), ("alignment_gain", gain), ("abs_phase_shift_h", np.abs(shift))):
            r, pp = perm_spearman(dS, y, rng)
            W.append(dict(contrast=label, n=len(ids), quantity=f"rho[d{name} vs {lab}]", pre=np.nan,
                          post=np.nan, delta=r, p=pp, test="perm spearman"))
        r0, p0 = perm_spearman(sc[gpre].values, P.amplitude.values, rng)
        W.append(dict(contrast=label, n=len(ids), quantity=f"rho[baseline {name} vs amplitude]",
                      pre=np.nan, post=np.nan, delta=r0, p=p0, test="perm spearman"))
    per.to_csv(f"{OUT}/clamp_per_person_{label}.tsv", sep="\t", index=False)


def pairs(d, m_pre, m_post):
    a = d[m_pre].drop_duplicates("subject").set_index("subject")
    b = d[m_post].drop_duplicates("subject").set_index("subject")
    ids = sorted(set(a.index) & set(b.index), key=int)
    return a.reindex(ids), b.reindex(ids)


# --------------------------------------------------------------------------- principal
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(0)
    E, d = load()
    print(f"{ACC}: {E.shape[0]} genes, {E.shape[1]} muestras, {d.subject.nunique()} sujetos")
    print(pd.crosstab([d.tissue, d.tzd], d.clamp.fillna("-")))

    tpath = f"{OUT}/template_2harm.tsv"
    if os.path.exists(tpath):
        spec = importlib.util.spec_from_file_location("pert04", os.path.join(HERE, "04_perturbations.py"))
        pert = importlib.util.module_from_spec(spec); spec.loader.exec_module(pert)
        template = pd.read_csv(tpath, sep="\t").set_index("gene")
        genes = list(template.sort_values("r2", ascending=False).index[:NTEMPLATE])
        estimator = ("2harm", (template, genes, pert))
        print(f"estimador: plantilla de dos armonicos, {len(genes)} genes, proyeccion con signo")
    else:
        template = pd.DataFrame({g: dict(phase_h=v[0], amp=v[1], r2=1.0) for g, v in CLOCK_TEMPLATE.items()}).T
        estimator = ("clock", template)
        print(f"estimador: SIN template_2harm.tsv -> fase anclada con {len(template)} genes del reloj "
              f"(docs/phase.js), amplitud por ordenamiento libre. Rehacer con la plantilla completa.")

    rows = dict(samples=[], groups=[], scores=[], genes=[], within=[])
    mus, adi = d.tissue == "muscle", d.tissue == "adipose"
    contrasts = [
        ("clamp_baseline", mus & (d.tzd == "pre") & (d.clamp == "pre"), mus & (d.tzd == "pre") & (d.clamp == "post")),
        ("clamp_postTZD", mus & (d.tzd == "post") & (d.clamp == "pre"), mus & (d.tzd == "post") & (d.clamp == "post")),
        ("tzd_muscle", mus & (d.tzd == "pre") & (d.clamp == "pre"), mus & (d.tzd == "post") & (d.clamp == "pre")),
        ("tzd_adipose", adi & (d.tzd == "pre"), adi & (d.tzd == "post")),
    ]
    for label, m1, m2 in contrasts:
        pre, post = pairs(d, m1, m2)
        print(f"\n== {label}: {len(pre)} pares ({(pre.sensitivity == 'sensitive').sum()} sensibles, "
              f"{(pre.sensitivity == 'resistant').sum()} resistentes)")
        if len(pre) < 6:
            print("  insuficiente"); continue
        run_contrast(E, pre, post, label, estimator, rng, rows)

    pd.concat(rows["samples"]).to_csv(f"{OUT}/clamp_samples.tsv", sep="\t", index=False)
    pd.DataFrame(rows["groups"]).to_csv(f"{OUT}/clamp_groups.tsv", sep="\t", index=False)
    pd.DataFrame(rows["scores"]).to_csv(f"{OUT}/clamp_scores_change.tsv", sep="\t", index=False)
    pd.DataFrame(rows["genes"]).to_csv(f"{OUT}/clamp_genes_change.tsv", sep="\t", index=False)
    W = pd.DataFrame(rows["within"]); W.to_csv(f"{OUT}/clamp_within_person.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200); pd.set_option("display.max_rows", 500)
    print("\n== B. cambio de scores post - pre ==")
    print(pd.DataFrame(rows["scores"]).drop(columns="genes").to_string(index=False))
    print("\n== B'. respuesta por grupo (mediana del cambio, sensibles frente a resistentes) ==")
    print(W[W.quantity.str.contains("sensitive vs resistant")].to_string(index=False))
    print("\n== C. acoplamiento dentro del individuo (solo |rho| >= 0.2 o P < 0.1) ==")
    Wc = W[W.quantity.str.startswith("rho")]
    print(Wc[(Wc.delta.abs() >= 0.2) | (Wc.p < 0.1)].to_string(index=False))
