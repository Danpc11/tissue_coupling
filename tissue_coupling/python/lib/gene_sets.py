"""Registro UNICO de los conjuntos de genes del articulo.

Un solo lugar de definicion para evitar que el codigo y las tablas diverjan. Cada conjunto lleva
su fuente. Ningun conjunto se deriva de los datos que se ponen a prueba.

Los scripts y las tablas suplementarias deben importar de aqui; nunca redefinir in situ.
"""

SETS = {
    "Genetic": dict(
        genes=["GSTM1", "GSTT1", "ERAP2"],
        source="Genes silenced by common structural variants: whole-gene deletions of GSTM1 and "
               "GSTT1, and the haplotype silencing ERAP2 (Andres et al., PLoS Genet 2010, "
               "doi:10.1371/journal.pgen.1001157). Pre-specified as the primary genetic set."),
    "Sex": dict(
        genes=["XIST", "RPS4Y1", "RPS4Y2", "USP9Y", "UTY", "KDM5D", "EIF1AY", "ZFY", "DDX3Y",
               "TXLNGY", "NLGN4Y", "TMSB4Y", "PRKY"],
        source="Y-chromosome genes and XIST. Included because chromosomal sex is identical in "
               "every tissue of a person, so these genes set the ceiling of cross-tissue "
               "concordance. They also serve as the sample-pairing control (see qc_sex_pairing)."),
    "HLA": dict(
        genes=["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1",
               "HLA-DRB4", "HLA-DRB5"],
        source="HLA genes whose expression depends on haplotype, including haplotype-dependent "
               "gene presence (DRB4, DRB5). Kept separate from the deletion set: merging the two "
               "raises between-cohort heterogeneity from I2 = 0% to 89%."),
    "Clock output": dict(
        genes=["DBP", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "CIART", "TEF", "HLF"],
        source="Output arm of the circadian oscillator; PAR-bZIP factors and Rev-erb targets."),
    "Core oscillator": dict(
        genes=["ARNTL", "CLOCK", "CRY1", "CRY2", "NPAS2"],
        source="Core transcriptional-translational feedback loop."),
    "Glucocorticoid": dict(
        genes=["FKBP5", "TSC22D3", "ZBTB16", "KLF15", "DDIT4", "SESN1", "GLUL", "SGK1", "CEBPD"],
        source="Canonical glucocorticoid-receptor target genes."),
    "Insulin core": dict(
        genes=["TXNIP", "PDK4", "HK2", "IRS2", "ANGPTL4", "RRAD"],
        source="Genes with the most consistent response to insulin in paired clamp biopsies."),
    "Oxidative": dict(
        genes=["PPARGC1A", "NDUFA4", "NDUFB8", "SDHB", "UQCRC1", "COX5A", "ATP5F1A", "CYCS",
               "ESRRA", "TFAM"],
        source="PGC-1alpha-associated oxidative phosphorylation programme (Mootha 2003, Patti 2003)."),
}

# Sinonimos: los conjuntos usan simbolos HGNC actuales, pero las plataformas antiguas estan
# anotadas con los de su epoca. Sin esta tabla, genes que SI estan medidos se declaran ausentes.
# Comprobado en AAGMEx (Illumina HumanHT-12, anotada en 2010): KDM5D y TXLNGY entran entre los
# treinta primeros de la lista de concordancia bajo sus nombres antiguos.
SYNONYMS = {
    "CIART": ["C1ORF51", "C1orf51", "CHRONO"],
    "KDM5D": ["JARID1D", "SMCY"],
    "TXLNGY": ["CYORF15A", "CYorf15A", "TXLNG2P"],
    "ATP5F1A": ["ATP5A1", "ATP5A"],
    "ATP5F1B": ["ATP5B"],
    "NDUFA4": ["NDUFA4L"],
    "TSC22D3": ["GILZ", "DSIPI"],
    "DDIT4": ["REDD1", "RTP801"],
    "SESN1": ["PA26"],
    "ANGPTL4": ["PGAR", "FIAF"],
    "RRAD": ["REM3"],
    "ARNTL": ["BMAL1", "MOP3"],
    "NR1D1": ["EAR1", "REVERBA"],
    "NR1D2": ["EAR1R", "REVERBB"],
    "FBXO32": ["ATROGIN1", "MAFBX"],
    "TRIM63": ["MURF1"],
    "SQSTM1": ["P62"],
    "MAP1LC3B": ["LC3B"],
    "CXCL8": ["IL8"],
    "ERRFI1": ["MIG6", "GENE33"],
    "HLA-DRB4": ["HLA-DRB4.1"],
}

def resolve(name, index):
    """Devuelve el identificador con el que un gen aparece en `index`, o None si no esta.

    Prueba el simbolo actual y despues sus sinonimos. Es necesario porque las matrices
    depositadas conservan la anotacion de su plataforma, que puede tener quince anos.
    """
    if name in index:
        return name
    for alt in SYNONYMS.get(name, []):
        if alt in index:
            return alt
    return None

def genes(name, index=None):
    """Genes del conjunto. Con `index`, devuelve los identificadores realmente presentes."""
    g = list(SETS[name]["genes"])
    if index is None:
        return g
    return [r for r in (resolve(x, index) for x in g) if r is not None]

def missing(name, index):
    """Genes del conjunto que no estan en `index` ni por sinonimo."""
    return [x for x in SETS[name]["genes"] if resolve(x, index) is None]

def qc_sex_pairing(M, A, ids):
    """Control de emparejamiento de muestras usando genes de sexo.

    Las dos muestras atribuidas a una persona deben venir de esa persona. En un estudio multiorgano
    los tejidos se procesan en lotes y momentos distintos, de modo que una etiqueta mal puesta o un
    intercambio de tubos asigna el musculo de alguien al adiposo de otro. XIST y los genes del
    cromosoma Y lo detectan: si un par esta cruzado entre un hombre y una mujer, el sexo inferido
    en los dos tejidos no coincide.

    Devuelve, por individuo, el sexo inferido en cada tejido y si coinciden.
    """
    import numpy as np, pandas as pd
    Y = [g for g in ("RPS4Y1", "DDX3Y", "USP9Y", "UTY", "EIF1AY", "KDM5D", "ZFY") if g in M.index]
    if not Y or "XIST" not in M.index:
        return None
    def call(E):
        y = E.loc[Y].mean(0); x = E.loc["XIST"]
        zy = (y - y.mean()) / (y.std() + 1e-9); zx = (x - x.mean()) / (x.std() + 1e-9)
        return pd.Series(np.where(zy > zx, "M", "F"), index=E.columns)
    cm, ca = call(M), call(A)
    d = pd.DataFrame({"muscle": cm.reindex(ids), "adipose": ca.reindex(ids)})
    d["agree"] = d.muscle == d.adipose
    return d

def as_table():
    import pandas as pd
    return pd.DataFrame([dict(gene_set=k, n_genes=len(v["genes"]),
                              genes="; ".join(v["genes"]), source=v["source"])
                         for k, v in SETS.items()])
