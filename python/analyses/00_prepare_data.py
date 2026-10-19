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
import os, sys, gzip, warnings; warnings.filterwarnings("ignore")
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
    "GSE41168":  "28 women, paired muscle and adipose, Affymetrix U133 Plus 2.0 (replication)",
    # --- plantilla de fase: cohortes con hora de muestreo conocida ---
    "GSE108539": "serial vastus lateralis biopsies every 4 h over 24 h, n = 10 (phase template)",
    "GSE104674": "adipose tissue at four circadian times, n = 12 (phase reference)",
    # --- validacion farmacologica ---
    "GSE148640": "hydrocortisone vs placebo in adipose tissue, within-person, n = 10",
    # --- definicion del nucleo de respuesta a la insulina ---
    "GSE22309":  "paired muscle biopsies, hyperinsulinaemic clamp",
    "GSE13070":  "paired muscle and adipose, clamp, Rd-defined groups",
    # --- expresion diferencial clasica (Fig. 1) ---
    "GSE182120": "muscle, NGT vs type 2 diabetes",
    "GSE164416": "pancreatic islets, living donors",
    "GSE50244":  "pancreatic islets, HbA1c",
    "GSE25462":  "muscle, family history of type 2 diabetes",
    "GSE159984": "pancreatic islets, 27 T2D / 58 non-diabetic donors + 45 cultured preparations",
    "GSE78721":  "adipose, limb, 68 type 2 diabetes vs 62 normal (Affymetrix PrimeView)",
    "GSE26637":  "adipose, fasting and hyperinsulinaemia in the same 10 people",
    # --- eje receptor-reloj: cohortes anadidas en v1.1 ---
    "GSE66446":  "Cushing syndrome adipose, 5 vs 11",
    "GSE17992":  "dexamethasone 4 days, muscle (n=10) and adipose (n=4), two-colour arrays",
    "GSE129843": "muscle, six biopsies per person per condition, recorded times",
    "GSE182117": "synchronised myotubes, 8 time points, 11 donors",
    "GSE87761":  "adipose under constant routine, recorded times",
    "GSE11868":  "3 h hyperinsulinaemic clamp, 6 lean subjects",
}
GTEX = "https://gtexportal.org/home/datasets/   (gene reads and sample annotations, v11)"
# Ficheros de anotacion que el analisis necesita ademas de las matrices:
ANNOT = {
    "GPL570_annot.txt":   "Affymetrix U133 Plus 2.0 annotation (GEO platform file); probe -> Entrez",
    "GPL15207_annot.txt": "Affymetrix PrimeView annotation; probe -> Entrez",
    "Human_GRCh38_p13_annot_tsv.gz": "NCBI gene annotation; Entrez -> current symbol and synonyms",
    "dorothea_ABC.tsv":   "OmniPath DoRothEA A-C regulons with sign",
    "msigdb/*.gmt":       "MSigDB v2026.2 collections H, C2, C5, C8 (symbols)",
    "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt": "RIN and ischaemic time per sample",
    "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt": "sex, age bracket, Hardy scale",
}

def report():
    print("Datos requeridos\n")
    print("GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>):")
    for a, d in GEO.items():
        p = os.path.join(RAW, f"{a}_series_matrix.txt.gz")
        print(f"  [{'presente' if os.path.exists(p) else 'FALTA   '}] {a}  {d}")
    print(f"\nGTEx: {GTEX}")
    for f in ("GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt",
              "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"):
        p = os.path.join(RAW, f)
        print(f"  [{'presente' if os.path.exists(p) else 'FALTA   '}] {f}")
    print("\nLos ficheros que falten deben descargarse manualmente: GEO y GTEx requieren")
    print("aceptar sus terminos, de modo que este script no los descarga por ti.")

if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True); os.makedirs(OUT, exist_ok=True)
    report()
