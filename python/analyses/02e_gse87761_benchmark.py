#!/usr/bin/env python3
"""Benchmark de la plantilla de musculo en adiposo con hora de muestreo registrada (GSE87761).

Results §2 afirma que la plantilla construida en musculo "predicted the recorded sampling time of
adipose tissue under constant routine (GSE87761) to within 1.7 h". GSE87761 (Christou et al.,
Agilent GPL15331): 7 varones en rutina constante, 5 biopsias de adiposo subcutaneo cada 6 h, 33
muestras con hora de reloj y fase circadiana respecto al inicio de melatonina en los metadatos.

Que hace:
  1. mapea las sondas a simbolo con la tabla de la plataforma y promedia por gen
  2. infiere la fase de cada muestra con la plantilla de dos armonicos de results/template_2harm.tsv
     (misma proyeccion firmada que usa el resto del pipeline: lib/within_person_coupling._profile)
  3. alinea el cero arbitrario de la plantilla con un unico desfase global y mide el error circular
     absoluto frente a la hora registrada y frente a la fase circadiana
  4. repite con varias listas de genes y con un nulo de fases al azar (mismo desfase optimo)

Uso:  python python/analyses/02e_gse87761_benchmark.py
Salidas: results/new/gse87761_benchmark.tsv (resumen por lista de genes),
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
    rows = []
    for name, genes, two in lists:
        phi, g = infer(E, T, genes, two)
        for ref, rname in ((md.clock_h.values, "clock"), (md.circ_phase_h.values, "circadian")):
            mae, med, off, R = circ_err(phi, ref)
            rows.append(dict(genes=name, n_genes=len(g), reference=rname, mae_h=mae, median_h=med, offset_h=off, R=R))
        if name == "clock_19_2harm":
            mae, med, off, R = circ_err(phi, md.clock_h.values)
            ph = md.copy(); ph["phase_inferred_h"] = (phi + off) % 24
            ph["abs_error_h"] = np.abs((ph.phase_inferred_h - ph.clock_h + 12) % 24 - 12)
            # error por sujeto (mediana de sus 4-5 muestras) y dejando cada sujeto fuera del desfase
            per_subj = ph.groupby("subject").abs_error_h.median()
            loso = []
            for s in ph.subject.unique():
                tr = ph.subject != s
                _, _, off_s, _ = circ_err(phi[tr.values], md.clock_h.values[tr.values])
                d = (phi[~tr.values] + off_s - md.clock_h.values[~tr.values] + 12) % 24 - 12
                loso.append(np.abs(d))
            loso = np.concatenate(loso)
            rows.append(dict(genes=name, n_genes=len(g), reference="clock_per_subject_median", mae_h=per_subj.mean(), median_h=per_subj.median(), offset_h=off, R=R))
            rows.append(dict(genes=name, n_genes=len(g), reference="clock_leave_subject_out", mae_h=loso.mean(), median_h=np.median(loso), offset_h=np.nan, R=R))
    rng = np.random.default_rng(0); null = []
    for _ in range(2000):
        null.append(circ_err(rng.uniform(0, 24, len(md)), md.clock_h.values, step=0.5)[0])
    rows.append(dict(genes="null_random_phases", n_genes=0, reference="clock", mae_h=np.mean(null), median_h=np.percentile(null, 5), offset_h=np.nan, R=np.nan))
    res = pd.DataFrame(rows)
    os.makedirs(f"{OUT}/new", exist_ok=True)
    res.to_csv(f"{OUT}/new/gse87761_benchmark.tsv", sep="\t", index=False, float_format="%.3f")
    ph.to_csv(f"{OUT}/new/gse87761_phases.tsv", sep="\t", float_format="%.3f")
    with pd.option_context("display.width", 160):
        print(res.round(2).to_string(index=False))
    main = res[(res.genes == "clock_19_2harm") & (res.reference == "clock")].iloc[0]
    print(f"\nGenes reloj, dos armonicos: error absoluto medio {main.mae_h:.2f} h, mediana {main.median_h:.2f} h, R {main.R:.2f} "
          f"(manuscrito: 'to within 1.7 h'; nulo {np.mean(null):.2f} h)")
