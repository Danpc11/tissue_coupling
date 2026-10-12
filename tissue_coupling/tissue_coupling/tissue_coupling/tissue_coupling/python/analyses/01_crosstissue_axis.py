#!/usr/bin/env python3
"""Eje individual compartido entre musculo y tejido adiposo: concordancia, metaanalisis y redes.

Analisis principal del articulo. Pregunta que genes situan a una persona de forma consistente
respecto a las demas en SUS DOS tejidos a la vez, cantidad distinta de la correlacion entre genes
dentro de un tejido.

Pasos:
  1. Concordancia entre tejidos: para cada gen, correlacion de Spearman ENTRE INDIVIDUOS de su
     expresion en musculo y en adiposo. El fondo es la media de esa cantidad sobre genes al azar.
  2. Significacion por conjunto frente a conjuntos aleatorios del mismo tamano y, donde se indica,
     equiparados ademas en expresion media y desviacion tipica (los 80 genes mas proximos en ese
     espacio bidimensional).
  3. Intervalos por bootstrap SOBRE INDIVIDUOS, no sobre genes: el individuo es la unidad de
     replicacion.
  4. Metaanalisis de efectos aleatorios (DerSimonian-Laird) sobre la z de Fisher, con I2 y Q.
  5. Dos redes: una sobre la concordancia entre tejidos, otra DENTRO de cada tejido con el mismo
     exponente y los mismos participantes, de modo que la conectividad sea comparable.

Los conjuntos de genes proceden de la literatura y NUNCA se derivan de los datos que se ponen a
prueba; se definen aqui de forma explicita para que queden congelados.

RESULTADO DE LOS DOS NULOS (AAGMEx, n = 232). El control equiparado es notablemente mas severo:

  conjunto          r      p (tamano)   p (expresion y variabilidad)
  Genetic        +0.808      0.001              0.002
  Clock output   +0.269      0.001              0.014
  Glucocorticoid +0.168      0.001              0.100
  Insulin core   +0.190      0.007              0.116
  Oscillator     +0.127      0.055              0.260
  Oxidative      +0.030      0.533              0.810

Solo el componente genetico y la salida del reloj superan ambos controles. Los glucocorticoides y
el nucleo de insulina superan el de tamano pero no el de expresion. Caben dos lecturas y los datos
no permiten elegir: que parte de su concordancia se deba a ser genes muy expresados y variables,
que se miden mejor y por tanto correlacionan mejor entre tejidos; o que el control sea demasiado
conservador, porque los genes equiparados a FKBP5 en expresion y variabilidad son ellos mismos
candidatos a estar bajo influencia hormonal. Debe reportarse el par de valores, no solo el primero.

Entradas: GSE95674 (adiposo) y GSE95675 (musculo) de AAGMEx; GSE41168 como replicacion.
Salidas: results/gtex/live_crosstissue_r.tsv, ci_*.tsv, meta_random_effects.tsv,
         crosstissue_modules.tsv, connectivity_AAGMEx.tsv
"""
import os, sys, gzip, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/gtex"); os.makedirs(OUT, exist_ok=True)
NBOOT = int(os.environ.get("T2D_BOOT", "400")); NPERM = int(os.environ.get("T2D_PERM", "3000"))

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib"))
from gene_sets import SETS as _REG, qc_sex_pairing, resolve as _resolve, missing as _missing
SETS = {k: v["genes"] for k, v in _REG.items()}   # registro unico: python/lib/gene_sets.py

def read_series(acc):
    meta, tab, intab = {}, [], False
    with gzip.open(f"{RAW}/{acc}_series_matrix.txt.gz", "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!series_matrix_table_begin"): intab = True; continue
            if l.startswith("!series_matrix_table_end"): break
            if intab: tab.append(l.rstrip("\n").split("\t")); continue
            for k, t in [("gsm", "!Sample_geo_accession"), ("title", "!Sample_title")]:
                if l.startswith(t) and k not in meta: meta[k] = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
    e = pd.DataFrame(tab[1:], columns=[c.strip('"') for c in tab[0]]).set_index(tab[0][0].strip('"'))
    e.index = [i.strip('"') for i in e.index]
    return pd.DataFrame(meta), e.apply(pd.to_numeric, errors="coerce").dropna()

def paired_aagmex():
    pa, ea = read_series("GSE95674"); pm, em = read_series("GSE95675")
    pa["s"] = pa.title.str.replace("-M", "", regex=False); pm["s"] = pm.title.str.replace("-M", "", regex=False)
    ids = sorted(set(pa.s) & set(pm.s))
    A = np.log2(ea[[pa.set_index("s").gsm[s] for s in ids]].clip(lower=1)); A.columns = ids
    M = np.log2(em[[pm.set_index("s").gsm[s] for s in ids]].clip(lower=1)); M.columns = ids
    g = sorted(set(A.index) & set(M.index))
    return M.loc[g], A.loc[g], ids

def concordance(M, A):
    return pd.Series([stats.spearmanr(M.loc[x], A.loc[x])[0] for x in M.index], index=M.index).dropna()

def matched_null(target, mu, sd, rng, pool):
    """Conjunto nulo equiparado en expresion media y variabilidad, no solo en tamano."""
    out = []
    for g in target:
        d = np.abs(mu - mu[g]) / mu.std() + np.abs(sd - sd[g]) / sd.std()
        out.append(rng.choice(d.sort_values().index[1:81]))
    return out

def boot_ci(M, A, genes, ids, rng, nboot=NBOOT):
    n = len(ids); vals = []
    for _ in range(nboot):
        i = rng.integers(0, n, n)
        vals.append(np.mean([stats.spearmanr(M.loc[x].values[i], A.loc[x].values[i])[0] for x in genes]))
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5), float(np.std(vals))

def dersimonian_laird(r, se):
    """Efectos aleatorios sobre z de Fisher; devuelve r combinada, IC, tau2, I2, Q."""
    z = np.arctanh(r); sez = se / (1 - r ** 2)
    w = 1 / sez ** 2; zbar = np.sum(w * z) / np.sum(w)
    Q = np.sum(w * (z - zbar) ** 2); k = len(z)
    Cc = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / Cc) if Cc > 0 else 0.0
    wr = 1 / (sez ** 2 + tau2); zm = np.sum(wr * z) / np.sum(wr); sem = np.sqrt(1 / np.sum(wr))
    I2 = max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0.0
    return np.tanh(zm), np.tanh(zm - 1.96 * sem), np.tanh(zm + 1.96 * sem), tau2, I2, Q, k

def crosstissue_network(M, A, con, thr=0.30, k=8):
    """Red SOBRE la concordancia entre tejidos: agrupa genes que ademas covarian entre individuos."""
    hi = con[con > thr].index.tolist()
    Zm = ((M.loc[hi].T - M.loc[hi].mean(1)) / (M.loc[hi].std(1) + 1e-9)).T
    Za = ((A.loc[hi].T - A.loc[hi].mean(1)) / (A.loc[hi].std(1) + 1e-9)).T
    Z = (Zm + Za) / 2
    Cm = np.corrcoef(Z.values); D = 1 - np.abs(Cm); np.fill_diagonal(D, 0); D = (D + D.T) / 2
    lab = fcluster(linkage(squareform(D, checks=False), method="average"), t=k, criterion="maxclust")
    return pd.Series(lab, index=hi)

def within_tissue_connectivity(M, A, beta=6, ntop=4000):
    """Conectividad dentro de cada tejido con el MISMO exponente y los MISMOS participantes."""
    v = (A.std(1).rank(pct=True) + M.std(1).rank(pct=True)) / 2
    sel = list(v.sort_values(ascending=False).index[:ntop])
    for s in set(sum(SETS.values(), [])):
        r = _resolve(s, A.index)
        if r is not None and r not in sel: sel.append(r)
    out = {}
    for nm, E in [("adipose", A), ("muscle", M)]:
        X = E.loc[sel]; Cm = np.corrcoef((X.T - X.mean(1)).T.values)
        Adj = np.abs(Cm) ** beta; np.fill_diagonal(Adj, 0)
        out[nm] = pd.Series(Adj.sum(1), index=sel)
    return pd.DataFrame(out)

PROBES_U133 = {"FKBP5": "204560_at", "TSC22D3": "207001_x_at", "ZBTB16": "205883_at",
 "KLF15": "219836_at", "DDIT4": "202887_s_at", "SESN1": "218346_s_at", "GLUL": "200648_s_at",
 "SGK1": "201739_at", "CEBPD": "203973_s_at", "DBP": "209782_s_at", "PER1": "202861_at",
 "PER2": "205251_at", "PER3": "1554462_a_at", "NR1D1": "204760_s_at", "NR1D2": "209750_at",
 "CIART": "229600_at", "TEF": "204018_x_at", "HLF": "204755_x_at", "ARNTL": "209747_at",
 "CLOCK": "204980_at", "CRY1": "209674_at", "CRY2": "203729_at", "NPAS2": "205125_at",
 "TXNIP": "201010_s_at", "PDK4": "205960_at", "HK2": "202934_at", "IRS2": "209184_s_at",
 "ANGPTL4": "221009_s_at", "RRAD": "204802_at", "PPARGC1A": "219195_at", "NDUFA4": "201483_s_at",
 "NDUFB8": "203189_s_at", "SDHB": "202675_at", "UQCRC1": "201903_at", "COX5A": "203663_s_at",
 "ATP5F1A": "213738_s_at", "CYCS": "208905_at", "ESRRA": "203193_at", "TFAM": "203177_x_at",
 "GSTM1": "204550_x_at", "GSTT1": "203815_at", "ERAP2": "219759_at", "HLA-DQA1": "203290_at",
 "HLA-DRB4": "209728_at"}

def paired_gse41168():
    """Replicacion independiente: 28 mujeres con musculo y adiposo, Affymetrix U133 Plus 2.0.
    Los identificadores son sondas, de modo que los conjuntos se traducen con PROBES_U133."""
    p, e = read_series("GSE41168")
    e = np.log2(e.clip(lower=1)) if e.values.max() > 60 else e
    p["s"] = p.title.str.split(",").str[0].str.strip()
    p["t"] = np.where(p.title.str.upper().str.contains("MUSCLE"), "m", "a")
    b = p[p.title.str.contains("Before") & p.gsm.isin(e.columns)]
    mu = b[b.t == "m"].drop_duplicates("s").set_index("s").gsm
    ad = b[b.t == "a"].drop_duplicates("s").set_index("s").gsm
    ids = sorted(set(mu.index) & set(ad.index))
    M = e[[mu[s] for s in ids]]; A = e[[ad[s] for s in ids]]
    M.columns = ids; A.columns = ids
    return M, A, ids

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    M, A, ids = paired_aagmex()
    print(f"AAGMEx: {len(ids)} individuos con ambos tejidos, {len(M)} genes comunes")
    qc = qc_sex_pairing(M, A, ids)
    if qc is not None:
        bad = qc[~qc.agree]
        print(f"control de emparejamiento por sexo: {qc.agree.sum()}/{len(qc)} pares coherentes"
              + (f" | DISCORDANTES: {list(bad.index)}" if len(bad) else ""))
        qc.rename_axis("individual").to_csv(f"{OUT}/qc_sex_pairing_AAGMEx.tsv", sep="\t")
    con = concordance(M, A)
    out = con.to_frame("r").rename_axis("gene")
    out.insert(0, "cohort", "AAGMEx"); out.insert(1, "platform", "Illumina HumanHT-12 v4")
    out.insert(2, "n_individuals", len(ids))
    out.to_csv(f"{OUT}/crosstissue_concordance_AAGMEx.tsv", sep="\t")
    bg = con.sample(min(3000, len(con)), random_state=0).mean()
    print(f"fondo: {bg:+.3f} | genes con r>0.5: {(con > 0.5).mean()*100:.1f}% | >0.7: {(con > 0.7).mean()*100:.2f}%")
    mu, sd = A.mean(1), A.std(1)
    rows = []
    for nm, gs in SETS.items():
        g = [r for r in (_resolve(x, con.index) for x in gs) if r is not None]
        miss = _missing(nm, con.index)
        if len(g) < 3: continue
        obs = con[g].mean()
        null_sz = [con.sample(len(g), random_state=int(rng.integers(1e9))).mean() for _ in range(NPERM)]
        null_mt = [con[matched_null(g, mu, sd, rng, con.index)].mean() for _ in range(min(NPERM, 500))]
        lo, hi, se = boot_ci(M, A, g, ids, rng)
        rows.append(dict(cohort="AAGMEx", set=nm, n=len(ids), n_declared=len(gs), n_genes=len(g),
                         genes_measured="; ".join(g), genes_missing="; ".join(sorted(miss)),
                         r=obs, lo=lo, hi=hi, se=se,
                         p_size=(np.sum(np.array(null_sz) >= obs) + 1) / (NPERM + 1),
                         p_matched=(np.sum(np.array(null_mt) >= obs) + 1) / (len(null_mt) + 1)))
        print(f"  {nm:16s} r={obs:+.3f} [{lo:+.3f}, {hi:+.3f}]  p={rows[-1]['p_size']:.4f} (matched {rows[-1]['p_matched']:.4f})")
    ci = pd.DataFrame(rows); ci.to_csv(f"{OUT}/ci_AAGMEx.tsv", sep="\t", index=False)
    mods = crosstissue_network(M, A, con)
    mods.to_frame("module").rename_axis("gene").to_csv(f"{OUT}/crosstissue_modules.tsv", sep="\t")
    print(f"\nred de concordancia: {mods.nunique()} modulos sobre {len(mods)} genes")
    K = within_tissue_connectivity(M, A); K.to_csv(f"{OUT}/connectivity_AAGMEx.tsv", sep="\t")
    print("conectividad dentro de tejido (percentil medio por conjunto):")
    for nm, gs in SETS.items():
        g = [x for x in gs if x in K.index]
        if len(g) < 3: continue
        pa_, pm_ = K.adipose.rank(pct=True)[g].mean(), K.muscle.rank(pct=True)[g].mean()
        print(f"  {nm:16s} adiposo {pa_:.3f} | musculo {pm_:.3f}")

    # ---- replicacion independiente: GSE41168 ----
    try:
        M2, A2, ids2 = paired_gse41168()
        print(f"\nGSE41168: {len(ids2)} mujeres con ambos tejidos")
        con2 = concordance(M2, A2)
        out2 = con2.to_frame("r").rename_axis("probe")
        out2.insert(0, "cohort", "GSE41168"); out2.insert(1, "platform", "Affymetrix U133 Plus 2.0")
        out2.insert(2, "n_individuals", len(ids2))
        out2.to_csv(f"{OUT}/crosstissue_concordance_GSE41168.tsv", sep="\t")
        mu2, sd2 = A2.mean(1), A2.std(1)
        rows2 = []
        for nm, gs in SETS.items():
            g = [PROBES_U133[x] for x in gs if x in PROBES_U133 and PROBES_U133[x] in con2.index]
            miss = [x for x in gs if x not in PROBES_U133 or PROBES_U133[x] not in con2.index]
            if len(g) < 3: continue
            obs = con2[g].mean()
            null_sz = [con2.sample(len(g), random_state=int(rng.integers(1e9))).mean() for _ in range(NPERM)]
            null_mt = [con2[matched_null(g, mu2, sd2, rng, con2.index)].mean() for _ in range(min(NPERM, 500))]
            lo, hi, se = boot_ci(M2, A2, g, ids2, rng)
            rows2.append(dict(cohort="GSE41168", set=nm, n=len(ids2), n_declared=len(gs), n_genes=len(g),
                              genes_measured="; ".join(g), genes_missing="; ".join(miss),
                              r=obs, lo=lo, hi=hi, se=se,
                              p_size=(np.sum(np.array(null_sz) >= obs) + 1) / (NPERM + 1),
                              p_matched=(np.sum(np.array(null_mt) >= obs) + 1) / (len(null_mt) + 1)))
            print(f"  {nm:16s} r={obs:+.3f} [{lo:+.3f}, {hi:+.3f}]  p={rows2[-1]['p_size']:.4f} (matched {rows2[-1]['p_matched']:.4f})")
        pd.DataFrame(rows2).to_csv(f"{OUT}/ci_GSE41168.tsv", sep="\t", index=False)
    except FileNotFoundError as ex:
        print("  (omitida) replicacion GSE41168:", ex)
