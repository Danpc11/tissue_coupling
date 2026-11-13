#!/usr/bin/env python3
"""Benchmark de la plantilla de musculo en adiposo con hora de muestreo registrada (GSE87761).

GSE87761 (Christou et al., Agilent GPL15331): 7 varones en rutina constante, 5 biopsias de adiposo
subcutaneo cada 6 h, 33 muestras con hora de reloj y fase circadiana respecto al inicio de
melatonina en los metadatos. La plantilla se construyo en otro tejido y en otra cohorte, de modo
que ni los genes ni sus fases se han elegido con estos datos.

Que hace:
  1. mapea las sondas a simbolo con la tabla de la plataforma y promedia por gen
  2. infiere la fase de cada muestra con la plantilla de dos armonicos de results/template_2harm.tsv
     (misma proyeccion firmada que usa el resto del pipeline: lib/within_person_coupling._profile)
  3. mide el error circular absoluto frente a la hora registrada

EL DESFASE. La plantilla fija la fase salvo un cero arbitrario: hay que estimar UN desfase entre el
cero de la plantilla y la hora de reloj, y ese parametro se aprende de datos con hora conocida. Hay
dos maneras de hacerlo y dan cifras distintas:
  - con todas las muestras, las mismas que luego se evaluan: es optimista, porque el desfase se
    elige minimizando el error que se informa;
  - con los demas sujetos, dejando fuera al que se evalua: es la estimacion honesta, y es la cifra
    PRINCIPAL de este paso (columna `primary`).
Para una prediccion del todo independiente el desfase tendria que venir de una tercera cohorte.
El azar se mide con el mismo procedimiento: fases al azar y desfase aprendido de los demas sujetos.

Se informan varias listas de genes. La principal es la de genes del reloj del registro, fijada de
antemano; las demas muestran cuanto depende el error de la lista, y su minimo no es una estimacion.

Uso:  python python/analyses/02e_gse87761_benchmark.py
Salidas: results/new/gse87761_benchmark.tsv (resumen por lista de genes y forma de fijar el desfase),
         results/new/gse87761_phases.tsv (fase inferida por muestra, lista de genes reloj)
"""
import os, sys, importlib.util, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../lib"))
import within_person_coupling as wpc
import phase as phase_lib

RAW = os.environ.get("T2D_RAW", "data/raw")
OUT = os.environ.get("T2D_OUT", "results")
from gene_sets import SETS as _SETS, panel as _panel
CLOCK = _panel("Benchmark clock")              # registro unico: python/lib/gene_sets.py
OUTPUT = list(_SETS["Clock output"]["genes"])


def mod(name, fn):
    sp = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m


def load():
    b = mod("b0", "00b_build_matrices.py")
    p, e = b.read_series("GSE87761")
    chars = {c: p[c].iloc[0].split(": ")[0] for c in p.columns if c.startswith("c")}
    col = {v: k for k, v in chars.items()}
    md = pd.DataFrame(index=p.gsm.values)
    md["subject"] = p[col["subject"]].str.split(": ").str[1].values
    hm = p[col["time sample taken"]].str.split(": ").str[1].str.split(":", expand=True).astype(float)
    md["clock_h"] = (hm[0] + hm[1] / 60).values
    md["circ_phase_h"] = p[col["circadian phase"]].str.split(": ").str[1].astype(float).values % 24
    ann = pd.read_csv(f"{RAW}/GPL15331_annot.txt", sep="\t", dtype=str, low_memory=False)
    symb = ann.set_index("ID")["GENE_SYMBOL"].dropna()
    e = e[md.index]
    e.index = e.index.map(lambda i: symb.get(i, np.nan))
    E = e[e.index.notna()].groupby(level=0).mean()     # ya en log2 (rango 1-18)
    return E, md


def circ_err(phi_h, ref_h, step=0.05):
    """Error absoluto circular tras alinear con el mejor desfase global; devuelve (mae, mediana, desfase, R)."""
    best = None
    for off in np.arange(0, 24, step):
        d = (phi_h + off - ref_h + 12) % 24 - 12; mae = np.abs(d).mean()
        if best is None or mae < best[0]:
            best = (mae, np.median(np.abs(d)), off)
    R = np.abs(np.mean(np.exp(1j * (phi_h - ref_h) * 2 * np.pi / 24)))
    return best + (R,)


def loso_err(phi_h, ref_h, subject):
    """Error absoluto de cada muestra con el desfase aprendido SOLO de los demas sujetos."""
    out = np.empty(len(phi_h))
    for s in pd.unique(subject):
        te = subject == s
        off = circ_err(phi_h[~te], ref_h[~te])[2]
        out[te] = np.abs((phi_h[te] + off - ref_h[te] + 12) % 24 - 12)
    return out


def infer(E, T, genes, two_harmonics=True):
    g = [x for x in genes if x in E.index and x in T.index]
    P, _, _ = wpc._profile(E, T, g, two_harmonics)
    return wpc.GRID[P.argmax(0)] * 24 / (2 * np.pi), g


if __name__ == "__main__":
    E, md = load()
    T = pd.read_csv(f"{OUT}/template_2harm.tsv", sep="\t").set_index("gene")
    print(f"GSE87761: {E.shape[1]} muestras, {md.subject.nunique()} sujetos, {E.shape[0]} genes; plantilla {len(T)} genes")
    lists = [("clock_19_2harm", CLOCK, True), ("clock_19_1harm", CLOCK, False),
             ("lib_phase_CLK_2harm", phase_lib.CLK, True), ("output_9_2harm", OUTPUT, True)]
    byr2 = T.sort_values("r2", ascending=False).index
    for n in (50, 100, 200, 500):
        lists.append((f"top{n}_r2_2harm", list(byr2[:n]), True))
    rows = []; subj = md.subject.values; clock = md.clock_h.values
    ALL, LOSO = "offset fitted on all samples (optimistic)", "offset from the other subjects"
    for name, genes, two in lists:
        phi, g = infer(E, T, genes, two)
        for ref, rname in ((clock, "clock"), (md.circ_phase_h.values, "circadian")):
            mae, med, off, R = circ_err(phi, ref)
            rows.append(dict(genes=name, n_genes=len(g), reference=rname, offset=ALL, primary=False,
                             mae_h=mae, median_h=med, within_3h=float(np.mean(np.abs((phi + off - ref + 12) % 24 - 12) <= 3)), offset_h=off, R=R))
        e = loso_err(phi, clock, subj)
        rows.append(dict(genes=name, n_genes=len(g), reference="clock", offset=LOSO, primary=name == "clock_19_2harm",
                         mae_h=e.mean(), median_h=np.median(e), within_3h=float(np.mean(e <= 3)), offset_h=np.nan, R=np.nan))
        if name == "clock_19_2harm":
            off = circ_err(phi, clock)[2]
            ph = md.copy(); ph["phase_inferred_h"] = (phi + off) % 24
            ph["abs_error_h"] = np.abs((ph.phase_inferred_h - ph.clock_h + 12) % 24 - 12)
            ph["abs_error_h_offset_from_other_subjects"] = e
            per_subj = ph.groupby("subject").abs_error_h_offset_from_other_subjects.median()
    # azar: fases al azar, con las dos formas de fijar el desfase
    rng = np.random.default_rng(0); null_all = []; null_loso = []
    for _ in range(2000):
        r = rng.uniform(0, 24, len(md))
        null_all.append(circ_err(r, clock, step=0.5)[0])
    for _ in range(300):
        null_loso.append(loso_err(rng.uniform(0, 24, len(md)), clock, subj).mean())
    for lab, v in ((ALL, null_all), (LOSO, null_loso)):
        rows.append(dict(genes="null_random_phases", n_genes=0, reference="clock", offset=lab, primary=False,
                         mae_h=np.mean(v), null_p05_mae_h=np.percentile(v, 5), offset_h=np.nan, R=np.nan))
    res = pd.DataFrame(rows)
    os.makedirs(f"{OUT}/new", exist_ok=True)
    res.to_csv(f"{OUT}/new/gse87761_benchmark.tsv", sep="\t", index=False, float_format="%.3f")
    ph.to_csv(f"{OUT}/new/gse87761_phases.tsv", sep="\t", float_format="%.3f")
    with pd.option_context("display.width", 200):
        print(res.round(2).to_string(index=False))
    main = res[res.primary].iloc[0]
    opt = res[(res.genes == "clock_19_2harm") & (res.reference == "clock") & (res.offset == ALL)].iloc[0]
    print(f"\nPRINCIPAL (genes del reloj, dos armonicos, desfase aprendido de los demas sujetos): error absoluto medio "
          f"{main.mae_h:.2f} h, mediana {main.median_h:.2f} h, {main.within_3h*100:.0f} % dentro de 3 h; "
          f"por sujeto, mediana {per_subj.median():.2f} h (de {per_subj.min():.2f} a {per_subj.max():.2f}). "
          f"Azar con el mismo procedimiento: {np.mean(null_loso):.2f} h (percentil 5: {np.percentile(null_loso, 5):.2f} h).")
    print(f"Con el desfase ajustado a las mismas muestras que se evaluan, que es optimista: {opt.mae_h:.2f} h de media, "
          f"{opt.median_h:.2f} h de mediana (azar {np.mean(null_all):.2f} h).")
