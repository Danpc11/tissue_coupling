#!/usr/bin/env python3
"""El estimador de fase frente a dos metodos publicados, con los mismos datos y los mismos pliegues.

ZeitZeiger (Hughey et al. 2016; Hughey 2017) y tauFisher (Duan, Ngo et al., Nat Commun 2024) predicen
la hora circadiana de una muestra a partir de un conjunto de entrenamiento con horas conocidas. Se
comparan con la plantilla del articulo en dos pruebas:

  1. DENTRO DE LA COHORTE: GSE108539 (musculo, 10 personas, biopsias cada 4 h), dejando a cada
     persona entera fuera. Los tres metodos ven exactamente los mismos pliegues, y en cada uno todo
     lo que el metodo elige (genes, componentes, configuracion) lo elige con el entrenamiento.
  2. TRANSFERENCIA: entrenar en GSE108539 y predecir GSE87761, que es otro tejido (adiposo), otra
     plataforma (microarreglo) y otro laboratorio. Entre musculo y adiposo puede haber un desfase
     real, y la plantilla del articulo fija la fase salvo un cero arbitrario; por eso el error se da
     de dos maneras para los tres metodos: tal cual sale, y con un unico desfase aprendido de los
     demas sujetos de la cohorte de prueba.

No se busca que el metodo propio gane. ZeitZeiger y tauFisher son metodos supervisados maduros; lo
que este paso deja escrito es donde queda la plantilla respecto a ellos con el mismo rasero.

R ES OPCIONAL Y SE DICE. Los dos metodos son paquetes de R que no estan en CRAN (se instalan desde
github.com/hugheylab/zeitzeiger y github.com/micnngo/tauFisher). Si no hay Rscript con los dos
paquetes, el paso escribe solo las filas del estimador propio y marca las otras como no ejecutadas.
T2D_RSCRIPT indica que Rscript usar; T2D_REQUIRE_R=1 convierte la ausencia en un error.

Salidas: results/new/benchmark_methods.tsv, benchmark_methods_predictions_*.tsv
"""
import os, sys, shutil, subprocess, tempfile, importlib.util, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "../lib"))
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results")
ROOT = os.path.abspath(os.path.join(HERE, "../.."))


def mod(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def err(pred_h, true_h):
    return np.abs((np.asarray(pred_h, float) - np.asarray(true_h, float) + 12) % 24 - 12)


def summary(e):
    e = np.asarray(e, float); e = e[np.isfinite(e)]
    return dict(n=len(e), mean_abs_error_h=e.mean(), median_abs_error_h=np.median(e), within_3h=float(np.mean(e <= 3))) if len(e) else dict(n=0)


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    tpl, bm = mod("tpl", "02b_phase_template.py"), mod("bm", "02e_gse87761_benchmark.py")
    # ---- cohorte de entrenamiento, igual que en el paso 02
    p, e = tpl.read_series(tpl.ACC); X = e if e is not None else tpl.load_counts(tpl.ACC, p)
    tcol = next(c for c in p.columns if c.lower() in ("time", "time point", "timepoint"))
    p["clock_h"] = pd.to_numeric(p[tcol].astype(str).str.split(":").str[0], errors="coerce")
    p["subject"] = p.title.astype(str).str.split("_").str[0]
    cols = [c for c in X.columns if c in set(p.gsm)]; q = p.set_index("gsm").reindex(cols)
    ok = tpl.annotated(X.index); Xv = X[cols].values.astype(float)
    t = q.clock_h.values * 2 * np.pi / 24; subject = q.subject.values
    # ---- cohorte de prueba
    E2, md2 = bm.load()
    T = pd.read_csv(f"{OUT}/template_2harm.tsv", sep="\t").set_index("gene")
    rows = []

    # ---- el estimador propio
    N, e_own = tpl.nested_loso(Xv, t, subject, ok)
    rows.append(dict(benchmark="leave-subject-out, GSE108539", method="template (this work, nested)", status="ok", **summary(e_own)))
    phi, _ = bm.infer(E2, T, bm.CLOCK, True)
    e_tr = bm.loso_err(phi, md2.clock_h.values, md2.subject.values)
    rows.append(dict(benchmark="transfer to GSE87761, offset from the other subjects", method="template (this work)", status="ok", **summary(e_tr)))
    rows.append(dict(benchmark="transfer to GSE87761, as predicted", method="template (this work)", status="not defined: the template has an arbitrary zero"))

    # ---- ZeitZeiger y tauFisher, en R
    rscript = os.environ.get("T2D_RSCRIPT") or shutil.which("Rscript")
    if rscript and not (os.path.exists(rscript) or shutil.which(rscript)): rscript = None
    have = False
    if rscript:
        chk = subprocess.run([rscript, "-e", "suppressMessages({library(zeitzeiger); library(tauFisher)})"], capture_output=True, text=True)
        have = chk.returncode == 0
    if not have:
        why = "Rscript not found" if not rscript else "R packages zeitzeiger and tauFisher not installed"
        print(f"ZeitZeiger y tauFisher: NO EJECUTADOS ({why}). La tabla lleva solo el estimador propio.")
        if os.environ.get("T2D_REQUIRE_R") == "1": sys.exit("T2D_REQUIRE_R=1 y R no esta disponible")
        for b in ("leave-subject-out, GSE108539", "transfer to GSE87761, as predicted", "transfer to GSE87761, offset from the other subjects"):
            for m in ("ZeitZeiger", "tauFisher"):
                rows.append(dict(benchmark=b, method=m, status=f"not run: {why}"))
    else:
        with tempfile.TemporaryDirectory() as wd:
            tr = X.loc[ok, cols]; tr = tr.loc[tr.var(1) > 0]
            tr.to_csv(f"{wd}/train_expr.tsv", sep="\t", index_label="gene")
            pd.DataFrame(dict(sample=cols, subject=subject, clock_h=q.clock_h.values)).to_csv(f"{wd}/train_meta.tsv", sep="\t", index=False)
            from gene_sets import panel
            open(f"{wd}/clock_genes.txt", "w").write("\n".join(panel("Benchmark clock") + ["ARNTL"]) + "\n")   # registro unico
            E2.to_csv(f"{wd}/test_expr.tsv", sep="\t", index_label="gene")
            pd.DataFrame(dict(sample=E2.columns, subject=md2.subject.values, clock_h=md2.clock_h.values)).to_csv(f"{wd}/test_meta.tsv", sep="\t", index=False)
            r = subprocess.run([rscript, os.path.join(ROOT, "R/benchmark_supervised.R"), wd], capture_output=True, text=True)
            if r.returncode != 0:
                print(r.stdout[-1500:]); print(r.stderr[-3000:]); sys.exit("el script de R fallo")
            L = pd.read_csv(f"{wd}/predictions_loso.tsv", sep="\t"); TR = pd.read_csv(f"{wd}/predictions_transfer.tsv", sep="\t")
        L.to_csv(f"{OUT}/new/benchmark_methods_predictions_loso.tsv", sep="\t", index=False)
        TR.to_csv(f"{OUT}/new/benchmark_methods_predictions_transfer.tsv", sep="\t", index=False)
        for m, col in (("ZeitZeiger", "zeitzeiger_h"), ("ZeitZeiger, genes standardised within cohort", "zeitzeiger_scaled_h"),
                       ("tauFisher", "taufisher_h")):
            def row(b, e, v):
                okv = np.isfinite(np.asarray(v, float))
                return dict(benchmark=b, method=m, status="ok" if okv.all() else f"{int((~okv).sum())} of {len(okv)} samples not predicted", **summary(e))
            if col in L.columns:
                rows.append(row("leave-subject-out, GSE108539", err(L[col], L.clock_h), L[col]))
            rows.append(row("transfer to GSE87761, as predicted", err(TR[col], TR.clock_h), TR[col]))
            v = TR[col].values.astype(float); g = np.isfinite(v)
            e = np.full(len(v), np.nan)
            if g.sum() > 5: e[g] = bm.loso_err(v[g], TR.clock_h.values[g], TR.subject.values[g])
            rows.append(row("transfer to GSE87761, offset from the other subjects", e, TR[col]))
    # ---- azar, con cada procedimiento
    rng = np.random.default_rng(0)
    rows.append(dict(benchmark="leave-subject-out, GSE108539", method="chance (uniform hour)", status="ok",
                     **summary(np.concatenate([err(rng.uniform(0, 24, len(q)), q.clock_h.values) for _ in range(200)]))))
    rows.append(dict(benchmark="transfer to GSE87761, as predicted", method="chance (uniform hour)", status="ok",
                     **summary(np.concatenate([err(rng.uniform(0, 24, len(md2)), md2.clock_h.values) for _ in range(200)]))))
    rows.append(dict(benchmark="transfer to GSE87761, offset from the other subjects", method="chance (uniform hour)", status="ok",
                     **summary(np.concatenate([bm.loso_err(rng.uniform(0, 24, len(md2)), md2.clock_h.values, md2.subject.values) for _ in range(100)]))))
    R = pd.DataFrame(rows).sort_values(["benchmark", "method"])
    R.to_csv(f"{OUT}/new/benchmark_methods.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(R.round(2).to_string(index=False))
