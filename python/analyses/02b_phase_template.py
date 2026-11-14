#!/usr/bin/env python3
"""Plantilla de fases relativas desde una cohorte con hora de muestreo conocida, y cuanto acierta.

Cohorte de referencia: GSE108539, diez voluntarios con vasto lateral biopsiado cada 4 h durante
24 h alternando piernas. El muestreo cada 4 h deja el armonico de 12 h dentro del limite de
Nyquist; el de 8 h ya no, y por eso no se ajusta.

Que hace:
  1. ajusta por gen un cosinor de 24 h (y de 24 + 12 h) sobre las horas registradas;
  2. mide el error de fase dejando SUJETOS ENTEROS fuera, porque las muestras de una misma
     persona no son independientes;
  3. exporta la plantilla ajustada con todos los sujetos: R2, amplitud y fase del maximo por gen.

LA VALIDACION NO USA NADA DEL SUJETO QUE SE DEJA FUERA. En cada pliegue, con las muestras de
entrenamiento y solo con ellas:
  - se calculan la media y la desviacion de cada gen, y con esos mismos parametros se
    estandarizan despues las muestras del sujeto retenido;
  - se ordenan los genes por R2 y se eligen los n primeros;
  - se ajustan los coeficientes de la plantilla. Se ajusta sin intercepto porque cada gen ya
    esta centrado con su media de entrenamiento; el R2 de la seleccion si lleva intercepto.
Hasta v2.2.1 la seleccion de genes y la estandarizacion se hacian con todas las muestras, incluidas
las del sujeto de prueba, y el error que se anunciaba era el de la mejor de doce configuraciones:
dos formas de optimismo.

DOS CIFRAS, Y CUAL ES LA QUE VALE.
  - phase_model_validation.tsv: el error de cada configuracion (armonicos x numero de genes), sin
    fuga. Sirve para ver como depende el error de la configuracion. El minimo de esa tabla NO es
    una estimacion del error: se ha elegido mirando el resultado.
  - phase_model_validation_nested.tsv: validacion anidada. Para cada sujeto retenido, la
    configuracion se elige con una segunda validacion por sujetos DENTRO del entrenamiento, y se
    aplica al sujeto retenido. Es la estimacion sin sesgo del error del procedimiento completo.

Uso:  python python/analyses/02b_phase_template.py
Salidas: results/muscle_rhythmic_genes.tsv, results/phase_model_validation.tsv,
         results/phase_model_validation_nested.tsv
"""
import os
import sys
import gzip
import warnings

warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
import numpy as np
import pandas as pd

RAW = os.environ.get("T2D_RAW", "data/raw")
OUT = os.environ.get("T2D_OUT", "results")
ACC = os.environ.get("T2D_CLOCK_COHORT", "GSE108539")
HARMONICS = [24.0, 12.0]          # periodos ajustados, en horas


def read_series(acc):
    """Lee la matriz de serie de GEO y los metadatos de muestra."""
    meta, tab, intab = {}, [], False
    with gzip.open(f"{RAW}/{acc}_series_matrix.txt.gz", "rt", errors="replace") as fh:
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


CONFIGS = [(harm, lab, n) for harm, lab in (([24.0], "24 h"), ([24.0, 12.0], "24 + 12 h"))
           for n in (13, 30, 60, 120, 250, 500)]


def annotated(names):
    """Mascara de genes con nombre propio: se descartan LOC, LINC y marcos de lectura sin nombre."""
    return np.array([not str(g).startswith(("LOC", "LINC")) and "orf" not in str(g).lower() for g in names])


def circular_error_h(pred, true):
    return np.abs(np.angle(np.exp(1j * (pred - true)))) * 24 / (2 * np.pi)


class Fold:
    """Todo lo que un pliegue aprende, y lo aprende solo de las muestras de entrenamiento `tr`."""

    def __init__(self, X, t, tr, ok):
        self.X, self.t, self.tr, self.ok = X, t, tr, ok
        self.mu = X[:, tr].mean(1, keepdims=True)
        self.sd = X[:, tr].std(1, keepdims=True) + 1e-9
        self.Ztr = (X[:, tr] - self.mu) / self.sd
        self._rank = {}

    def ranking(self, harm):
        """Genes ordenados por R2 del cosinor en el entrenamiento."""
        key = tuple(harm)
        if key not in self._rank:
            r2 = fit_cosinor(self.Ztr, self.t[self.tr], harm)[0]
            self._rank[key] = [i for i in np.argsort(-r2) if self.ok[i]]
        return self._rank[key]

    def predict(self, te, harm, n):
        """Fase de las muestras `te`, estandarizadas con los parametros de entrenamiento."""
        sel = self.ranking(harm)[:n]
        B, *_ = np.linalg.lstsq(design(self.t[self.tr], harm, intercept=False), self.Ztr[sel].T, rcond=None)
        Zte = (self.X[sel][:, te] - self.mu[sel]) / self.sd[sel]
        return predict(Zte, B, harm)


def loso_by_config(X, t, subject, ok):
    """Error de cada configuracion dejando sujetos fuera, sin fuga. Devuelve {config: errores}."""
    errs = {c: [] for c in range(len(CONFIGS))}
    for s in pd.unique(subject):
        te = subject == s
        fold = Fold(X, t, ~te, ok)
        for c, (harm, _, n) in enumerate(CONFIGS):
            errs[c].extend(circular_error_h(fold.predict(te, harm, n), t[te]))
    return {c: np.array(e) for c, e in errs.items()}


def nested_loso(X, t, subject, ok):
    """Validacion anidada: la configuracion se elige dentro del entrenamiento de cada pliegue."""
    rows, errs = [], []
    for s in pd.unique(subject):
        te = subject == s
        tr_idx = np.where(~te)[0]
        inner = loso_by_config(X[:, tr_idx], t[tr_idx], subject[tr_idx], ok)
        best = min(range(len(CONFIGS)), key=lambda c: (np.median(inner[c]), c))
        harm, lab, n = CONFIGS[best]
        e = circular_error_h(Fold(X, t, ~te, ok).predict(te, harm, n), t[te])
        errs.extend(e)
        rows.append(dict(held_out_subject=s, n_samples=int(te.sum()), chosen_harmonics=lab, chosen_n_genes=n,
                         inner_median_error_h=np.median(inner[best]), median_error_h=np.median(e)))
    return pd.DataFrame(rows), np.array(errs)


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
    p["clock_h"] = pd.to_numeric(p[tcol].astype(str).str.split(":").str[0], errors="coerce")
    p["subject"] = p.title.astype(str).str.split("_").str[0]
    cols = [c for c in X.columns if c in set(p.gsm)]
    q = p.set_index("gsm").reindex(cols)
    # los metadatos se leen del titulo y de un campo de texto: se comprueba que salio lo esperado
    problems = []
    if q.clock_h.isna().any():
        problems.append(f"{int(q.clock_h.isna().sum())} muestras sin hora legible en '{tcol}'")
    elif not q.clock_h.between(0, 24).all():
        problems.append("horas fuera de 0-24")
    per_subject = q.groupby("subject").size()
    if len(per_subject) < 3 or (per_subject == 1).mean() > 0.5:
        problems.append(f"el sujeto no se lee del titulo (p. ej. '{q.title.iloc[0]}'): salen {len(per_subject)} "
                        f"sujetos, {int((per_subject == 1).sum())} con una sola muestra")
    if problems:
        print(f"{ACC}: metadatos no interpretables: " + "; ".join(problems))
        sys.exit(1)
    t = q.clock_h.values * 2 * np.pi / 24
    subject = q.subject.values
    Xv = X[cols].values.astype(float)
    ok = annotated(X.index)
    Z = X[cols]
    Z = ((Z.T - Z.mean(1)) / (Z.std(1) + 1e-9)).T            # solo para la plantilla final, con todos
    print(f"{ACC}: {len(cols)} muestras, {q.subject.nunique()} sujetos "
          f"(muestras por sujeto: {sorted(per_subject.values.tolist())}), horas {sorted(set(q.clock_h))}")

    # --- error de cada configuracion, dejando SUJETOS fuera y sin usar nada del sujeto retenido ---
    by_config = loso_by_config(Xv, t, subject, ok)
    rows = []
    for c, (harm, lab, n) in enumerate(CONFIGS):
        e = by_config[c]
        rows.append(dict(harmonics=lab, n_genes=n, median_error_h=np.median(e), within_3h=float(np.mean(e <= 3))))
        print(f"  {lab:10s} n={n:4d}  error mediano {np.median(e):.2f} h  ±3 h {np.mean(e <= 3) * 100:.0f} %")
    v = pd.DataFrame(rows)
    v.to_csv(f"{OUT}/phase_model_validation.tsv", sep="\t", index=False)

    # --- validacion anidada: la configuracion se elige dentro del entrenamiento ------------------
    N, e = nested_loso(Xv, t, subject, ok)
    N.to_csv(f"{OUT}/phase_model_validation_nested.tsv", sep="\t", index=False)
    chosen = (N.chosen_harmonics + ", " + N.chosen_n_genes.astype(str)).value_counts()
    print("\nvalidacion anidada (la configuracion se elige sin ver al sujeto retenido):")
    print(f"  error mediano {np.median(e):.2f} h, medio {np.mean(e):.2f} h, dentro de ±3 h {np.mean(e <= 3) * 100:.0f} %, "
          f"{len(e)} muestras de {len(N)} sujetos")
    print("  configuraciones elegidas: " + "; ".join(f"{k} genes, en {n} de {len(N)} pliegues" for k, n in chosen.items()))

    # --- plantilla final: cosinor con los dos armonicos sobre todos los genes ----------------
    r2, amp, ph = fit_cosinor(Z.values, t, HARMONICS)
    tpl = pd.DataFrame({"gene": Z.index, "r2": r2, "amp": amp, "phase_h": ph})
    tpl = tpl.sort_values("r2", ascending=False)
    tpl.insert(1, "cohort", ACC)
    tpl.insert(2, "n_subjects", q.subject.nunique())
    tpl.to_csv(f"{OUT}/muscle_rhythmic_genes.tsv", sep="\t", index=False)
    print(f"\nplantilla: {len(tpl)} genes -> {OUT}/muscle_rhythmic_genes.tsv")
    print("el error que describe al estimador es el de la validacion anidada; el minimo de la tabla por "
          "configuracion no lo es, porque se elige mirando el resultado.")
