#!/usr/bin/env python3
"""Descarga y prepara los datos publicos que el analisis necesita.

Ningun dato se redistribuye aqui: todo se obtiene de GEO y del portal de GTEx. Este script
documenta exactamente que se descarga, como se anota y como se emparejan los tejidos, que es la
parte donde es facil equivocarse de cohorte sin darse cuenta.

Dos cautelas implementadas:
  - Cada matriz producida declara su cohorte y su plataforma en columnas propias. Sin eso es
    posible mezclar ficheros de cohortes distintas en un mismo analisis, que es un error real
    que cometimos y que el verificador ahora impide.
  - Los simbolos de gen se resuelven contra el registro con tabla de sinonimos, porque las
    matrices depositadas conservan la anotacion de su plataforma, que puede tener quince anos.
"""
import os, sys, gzip, warnings; warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib"))

RAW = os.environ.get("T2D_RAW", "data/raw")
OUT = os.environ.get("T2D_OUT", "results")

# Cohortes necesarias, agrupadas por el papel que cumplen. Las del primer grupo son
# imprescindibles; sin las demas el orquestador omite los pasos correspondientes.
GEO = {
    # --- analisis principal: dos tejidos de la misma persona ---
    "GSE95674":  "AAGMEx adipose, 256 samples, Illumina HumanHT-12 v4 (primary)",
    "GSE95675":  "AAGMEx muscle, 255 samples; 232 individuals have both tissues (primary)",
    "GSE41168":  "muscle and adipose before treatment, 35 people per tissue, 28 with both, Affymetrix U133 Plus 2.0 (replication)",
    # --- plantilla de fase: cohortes con hora de muestreo conocida ---
    "GSE108539": "serial vastus lateralis biopsies every 4 h over 24 h, n = 10 (phase template)",
    "GSE107537":  "10 h sleep delay under controlled laboratory conditions, blood, recorded hours (external control of the clock-output gene set, by cosinor)",
    # --- validacion farmacologica ---
    "GSE148640": "hydrocortisone vs placebo in adipose tissue, within-person, n = 10",
    # --- definicion del nucleo de respuesta a la insulina ---
    "GSE13070":  "paired muscle and adipose, clamp, Rd-defined groups",
    # --- expresion diferencial clasica (Fig. 1) ---
    "GSE182120": "muscle, NGT vs type 2 diabetes",
    "GSE25462":  "muscle, family history of type 2 diabetes",
    "GSE159984": "pancreatic islets, 27 T2D / 58 non-diabetic donors + 45 cultured preparations",
    "GSE78721":  "adipose, 130 samples of 95 people, type 2 diabetes vs normal; one sample per person is used (Affymetrix PrimeView)",
    "GSE26637":  "adipose, fasting and hyperinsulinaemia in the same 10 people",
    # --- eje receptor-reloj: cohortes anadidas en v1.1 ---
    "GSE66446":  "Cushing syndrome adipose, 5 vs 11",
    "GSE17992":  "dexamethasone 4 days, muscle (n=10) and adipose (n=4), two-colour arrays",
    "GSE129843": "muscle, six biopsies per person per condition, recorded times",
    "GSE182117": "synchronised myotubes, time course in two conditions, 12 donors (validation of the amplitude estimator against cosinor)",
    "GSE87761":  "adipose under constant routine, recorded times",
    "GSE40589":  "adipose, saline infusion vs euglycaemic and hyperglycaemic clamp at the same clock times (time-matched control for the clamp)",
}
GTEX = "https://gtexportal.org/home/datasets/   (gene reads and sample annotations, v11)"
# Ficheros de anotacion que el analisis necesita ademas de las matrices:
ANNOT = {
    "GPL570_annot.txt":   "Affymetrix U133 Plus 2.0 annotation (GEO platform file); probe -> Entrez",
    "GPL15207_annot.txt": "Affymetrix PrimeView annotation; probe -> Entrez",
    "Human_GRCh38_p13_annot_tsv.gz": "NCBI gene annotation; Entrez -> current symbol and synonyms",
    "dorothea_ABC.tsv":   ("OmniPath DoRothEA A-C regulons with sign. Download with:\n"
        "      curl -o $T2D_RAW/dorothea_ABC.tsv 'https://omnipathdb.org/interactions"
        "?datasets=dorothea&dorothea_levels=A,B,C&organisms=9606&genesymbols=yes"
        "&fields=dorothea_level,sources,references'"),
    "msigdb/*.gmt":       "MSigDB v2026.2 collections H, C2, C5, C8 (symbols)",
    "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt": "RIN and ischaemic time per sample",
    "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt": "sex, age bracket, Hardy scale",
}

# Lo que el pipeline lee ademas de las matrices de serie. Si falta algo de aqui, o si un .gz no es un
# gzip (NCBI devuelve a veces una pagina de CAPTCHA con el nombre del fichero), la corrida fallaria
# pasos mas tarde; se dice aqui, antes de calcular nada.
FILES = {
    "cuentas y tablas de las series": [
        "GSE159984_raw_counts.tsv.gz", "GSE159984-GPL9115_series_matrix.txt.gz", "GSE159984-GPL16791_series_matrix.txt.gz",
        "GSE129843_RESTRICT.txt.gz", "GSE129843_README.txt", "GSE182117_counts.tsv.gz", "GSE148640_processed.txt.gz"],
    "recount3 (GSE108539 y GSE66446)": [f"recount3/sra.{k}.{s}.{e}" for s in ("SRP127581", "SRP055749")
                                         for k, e in (("gene_sums", "G026.gz"), ("sra", "MD.gz"), ("recount_qc", "MD.gz"))],
    "GTEx, cuentas por tejido": [f"gene_reads_v11_{t}.gct.gz" for t in ("muscle_skeletal", "adipose_subcutaneous", "adipose_visceral_omentum")],
    "anotacion de plataformas y de genes": ["GPL570_annot.txt", "GPL15207_annot.txt", "GPL15331_annot.txt", "GPL17586_annot.txt",
                                            "GPL4133_annot.txt", "GPL6848_annot.txt", "GPL24324-69524.txt",
                                            "_soft/Homo_sapiens.gene_info.gz", "Human_GRCh38_p13_annot_tsv.gz", "dorothea_ABC.tsv"],
    "MSigDB v2026.2": [f"msigdb/{c}.v2026.2.Hs.symbols.gmt" for c in ("h.all", "c2.cp.reactome", "c2.cp.wikipathways", "c2.cp.kegg_medicus",
                       "c2.cp.pid", "c2.cgp", "c3.tft.gtrd", "c3.mir.mirdb", "c4.all", "c5.go.bp", "c5.hpo", "c7.all", "c8.all")],
}


def state(path):
    """presente, FALTA, o NO ES GZIP para un .gz que no empieza por la firma de gzip."""
    if not os.path.exists(path): return "FALTA"
    if path.endswith(".gz"):
        with open(path, "rb") as fh:
            if fh.read(2) != b"\x1f\x8b": return "NO ES GZIP"
    return "presente"


def report():
    missing = 0
    print("Datos requeridos\n")
    print("GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>):")
    for a, d in GEO.items():
        p = os.path.join(RAW, f"{a}_series_matrix.txt.gz")
        ok = os.path.exists(p); missing += 0 if ok else 1
        print(f"  [{'presente' if ok else 'FALTA   '}] {a}  {d}")
    print(f"\nGTEx: {GTEX}")
    for f in ("GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt",
              "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"):
        p = os.path.join(RAW, f)
        ok = os.path.exists(p); missing += 0 if ok else 1
        print(f"  [{'presente' if ok else 'FALTA   '}] {f}")
    for title, files in FILES.items():
        bad = [(f, state(os.path.join(RAW, f))) for f in files]; bad = [(f, s) for f, s in bad if s != "presente"]
        missing += len(bad)
        print(f"\n{title}: {len(files) - len(bad)} de {len(files)} presentes")
        for f, s in bad:
            print(f"  [{s:10s}] {f}" + ("   <- NCBI devolvio una pagina web en vez del fichero: descarguelo con el navegador (docs/DATA.md)" if s == "NO ES GZIP" else ""))
    for a in GEO:
        s = state(os.path.join(RAW, f"{a}_series_matrix.txt.gz"))
        if s == "NO ES GZIP": missing += 1; print(f"  [NO ES GZIP] {a}_series_matrix.txt.gz")
    print("\nLos ficheros que falten deben descargarse manualmente: GEO y GTEx requieren")
    print("aceptar sus terminos, de modo que este script no los descarga por ti.")
    return missing

if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True); os.makedirs(OUT, exist_ok=True)
    if report():
        raise SystemExit(1)
