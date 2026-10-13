#!/usr/bin/env python3
"""Plantilla de fases relativas desde cohortes con hora de muestreo conocida.

Este paso produce `results/muscle_rhythmic_genes.tsv`, que el paso 03 necesita para anclar el
periodo. Sin el, la fase solo puede estimarse con ordenamiento libre, sujeta a la ambiguedad entre
el fundamental de 24 h y su armonico de 12 h.

Que hace:
  1. ajusta por gen un cosinor de 24 h sobre una cohorte con horas conocidas (biopsias seriadas)
  2. valida el modelo dejando SUJETOS ENTEROS fuera, no muestras, porque las muestras de una misma
     persona no son independientes
  3. compara 24 h solo frente a 24 + 12 h, y barre el numero de genes
  4. exporta la plantilla: para cada gen su R2, amplitud y fase del maximo

Cohorte de referencia: GSE108539, diez voluntarios con vasto lateral biopsiado cada 4 h durante
24 h alternando piernas. El muestreo cada 4 h deja el armonico de 12 h dentro del limite de
Nyquist; el de 8 h ya no, y por eso no se ajusta.

Uso:  python python/analyses/02b_phase_template.py
Salidas: results/muscle_rhythmic_genes.tsv, results/phase_model_validation.tsv
"""
import os
import sys
import gzip
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

RAW = os.environ.get("T2D_RAW", "data/raw")
OUT = os.environ.get("T2D_OUT", "results")
ACC = os.environ.get("T2D_CLOCK_COHORT", "GSE108539")
HARMONICS = [24.0, 12.0]          # periodos ajustados, en horas


def read_series(acc):
    """Lee la matriz de serie de GEO y los metadatos de muestra."""
    meta, tab, intab = {}, [], False
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
            if l.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
                kk = v[0].split(":")[0].strip()
                if kk not in meta:
                    meta[kk] = [x.split(":", 1)[1].strip() if ":" in x else "" for x in v]
    p = pd.DataFrame({k: v for k, v in meta.items() if len(v) == len(meta["gsm"])})
    if len(tab) > 1:
        e = pd.DataFrame(tab[1:], columns=[c.strip('"') for c in tab[0]]).set_index(tab[0][0].strip('"'))
        e.index = [i.strip('"') for i in e.index]
        e = e.apply(pd.to_numeric, errors="coerce").dropna()
    else:
        e = None
    return p, e


def load_counts(acc, p):
    """Matriz de cuentas en log2 CPM. Si la serie no trae tabla, busca el fichero de cuentas."""
    for fn in (f"{RAW}/{acc}_raw_counts_GRCh38_p13_NCBI.tsv.gz", f"{RAW}/{acc}_counts.tsv.gz"):
        if os.path.exists(fn):
            c = pd.read_csv(fn, sep="\t", index_col=0)
            ann = f"{RAW}/Human_GRCh38_p13_annot.tsv.gz"
            if os.path.exists(ann) and str(c.index[0]).isdigit():
                a = pd.read_csv(ann, sep="\t", usecols=["GeneID", "Symbol"], index_col=0)
                c.index = a.Symbol.reindex(c.index).values
                c = c[pd.notna(c.index)].groupby(level=0).sum()
            X = np.log2(c / c.sum() * 1e6 + 1)
            return X[(X > 1).mean(1) > 0.5]
    return None


def design(t_rad, harmonics=HARMONICS, intercept=True):
    cols = [np.ones(len(t_rad))] if intercept else []
    for h in harmonics:
        cols += [np.cos(t_rad * 24.0 / h), np.sin(t_rad * 24.0 / h)]
    return np.column_stack(cols)


def fit_cosinor(Z, t_rad, harmonics=HARMONICS):
    """Devuelve R2, amplitud y fase del maximo por gen, para el fundamental de 24 h."""
    D = design(t_rad, harmonics)
    B, *_ = np.linalg.lstsq(D, Z.T, rcond=None)
    fit = (D @ B).T
    ss_res = ((Z - fit) ** 2).sum(1)
    ss_tot = ((Z - Z.mean(1, keepdims=True)) ** 2).sum(1)
    r2 = 1 - ss_res / (ss_tot + 1e-12)
    a, b = B[1], B[2]                      # terminos del fundamental
    return r2, np.sqrt(a ** 2 + b ** 2), (np.arctan2(b, a) * 24 / (2 * np.pi)) % 24


def predict(Zs, B, harmonics=HARMONICS, grid=721):
    g = np.linspace(0, 2 * np.pi, grid, endpoint=False)
    E = design(g, harmonics, intercept=False) @ B
    E = (E - E.mean(1, keepdims=True)) / (E.std(1, keepdims=True) + 1e-9)
    O = (Zs - Zs.mean(0)) / (Zs.std(0) + 1e-9)
    return g[np.argmax(E @ O, axis=0)]


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    p, e = read_series(ACC)
    X = e if e is not None else load_counts(ACC, p)
    if X is None:
        print(f"No encuentro la matriz de {ACC}. Ejecuta primero 00_prepare_data.py.")
        sys.exit(1)
    tcol = next((c for c in p.columns if c.lower() in ("time", "time point", "timepoint")), None)
    if tcol is None:
        print(f"{ACC} no declara la hora de muestreo; no puede servir de plantilla.")
        sys.exit(1)
    p["clock_h"] = p[tcol].astype(str).str.split(":").str[0].astype(float)
    p["subject"] = p.title.astype(str).str.split("_").str[0]
    cols = [c for c in X.columns if c in set(p.gsm)]
    q = p.set_index("gsm").reindex(cols)
    t = q.clock_h.values * 2 * np.pi / 24
    Z = X[cols]
    Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T
    print(f"{ACC}: {len(cols)} muestras, {q.subject.nunique()} sujetos, "
          f"horas {sorted(set(q.clock_h))}")

    # --- validacion dejando SUJETOS fuera, y barrido de numero de genes y armonicos ----------
    rows = []
    for harm, lab in [([24.0], "24 h"), ([24.0, 12.0], "24 + 12 h")]:
        r2_all, _, _ = fit_cosinor(Z.values, t, harm)
        cand = [g for g in pd.Series(r2_all, index=Z.index).sort_values(ascending=False).index
                if not str(g).startswith(("LOC", "LINC")) and "orf" not in str(g).lower()]
        for n in (13, 30, 60, 120, 250, 500):
            sel = cand[:n]
            Zs = X.loc[sel, cols]
            Zs = ((Zs.T - Zs.mean(1)) / (Zs.std(1) + 1e-9)).T.values
            err = []
            for s in q.subject.unique():
                te = (q.subject == s).values
                tr = ~te
                B, *_ = np.linalg.lstsq(design(t[tr], harm, intercept=False), Zs[:, tr].T, rcond=None)
                pred = predict(Zs[:, te], B, harm)
                err.extend(np.angle(np.exp(1j * (pred - t[te]))) * 24 / (2 * np.pi))
            err = np.abs(np.array(err))
            rows.append(dict(harmonics=lab, n_genes=n, median_error_h=np.median(err),
                             within_3h=float(np.mean(err <= 3))))
            print(f"  {lab:10s} n={n:4d}  error mediano {np.median(err):.2f} h  "
                  f"±3 h {np.mean(err) * 0 + np.mean(err <= 3) * 100:.0f} %")
    v = pd.DataFrame(rows)
    v.to_csv(f"{OUT}/phase_model_validation.tsv", sep="\t", index=False)

    # --- plantilla final: cosinor con los dos armonicos sobre todos los genes ----------------
    r2, amp, ph = fit_cosinor(Z.values, t, HARMONICS)
    tpl = pd.DataFrame({"gene": Z.index, "r2": r2, "amp": amp, "phase_h": ph})
    tpl = tpl.sort_values("r2", ascending=False)
    tpl.insert(1, "cohort", ACC)
    tpl.insert(2, "n_subjects", q.subject.nunique())
    tpl.to_csv(f"{OUT}/muscle_rhythmic_genes.tsv", sep="\t", index=False)
    best = v.loc[v.median_error_h.idxmin()]
    print(f"\nplantilla: {len(tpl)} genes -> {OUT}/muscle_rhythmic_genes.tsv")
    print(f"mejor configuracion en validacion: {best.harmonics}, {int(best.n_genes)} genes, "
          f"{best.median_error_h:.2f} h")
