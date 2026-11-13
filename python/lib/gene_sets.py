"""Registro UNICO de los conjuntos de genes del articulo.

Un solo lugar de definicion para evitar que el codigo y las tablas diverjan. Cada conjunto lleva
su fuente. Ningun conjunto se deriva de los datos que se ponen a prueba.

Hay dos tablas:

  SETS    los conjuntos que se PUNTUAN y se comparan entre si (salida del reloj, oscilador
          central, glucocorticoide, insulina, oxidativo) y los controles de identidad (genetico,
          sexo, HLA). Varios pasos recorren SETS entero, de modo que anadir aqui un conjunto
          anade una fila a esos analisis.
  PANELS  los paneles de genes que usa un estimador, un control o un paso concreto, y que no se
          puntuan junto a los de SETS: los genes con los que se ordena la fase, los que recibe
          CHIRAL, el programa de ayuno del clamp, las dianas disjuntas de la prediccion entre
          tejidos y la firma del modelo.

La regla: todo conjunto que define una puntuacion, una inferencia o un control vive aqui, y los
scripts lo importan; nunca se redefine in situ. Fuera de este fichero solo quedan las listas que
eligen que genes dibuja un panel (en python/figures/) y las tablas de parametros por gen (fase de
pico del modelo, sonda por gen de una plataforma), que tests/test_gene_lists.py enumera una a una.
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

# Paneles de genes de pasos concretos. Se guardan tal como se usan: mismo contenido y mismo orden
# que tenian en el script que los definia, de modo que moverlos aqui no cambia ningun resultado.
PANELS = {
    "Phase ordering": dict(
        genes=["BMAL1", "NPAS2", "CLOCK", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "DBP", "TEF", "HLF", "CIART",
               "CRY1", "CRY2", "BHLHE40", "BHLHE41", "NFIL3", "RORA", "PDK4", "CEBPD", "FKBP5", "KLF15", "TXNIP"],
        used_by="python/lib/phase.py (CLK)",
        source="Clock genes plus clock-controlled metabolic and glucocorticoid-responsive genes with "
               "robust 24-h rhythms in human muscle and adipose tissue. It is a panel for ordering "
               "samples around the cycle, not a circadian signature to be scored: the metabolic genes "
               "are there because they are rhythmic, and the estimator is run with and without them."),
    "CHIRAL input": dict(
        genes=["ARNTL", "NPAS2", "CLOCK", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "DBP", "TEF", "HLF",
               "CIART", "CRY1", "CRY2", "BHLHE40", "NFIL3", "RORA", "PDK4", "CEBPD", "FKBP5", "TXNIP", "KLF15"],
        used_by="python/analyses/17_chiral_concordance.py",
        source="The phase-ordering panel as it was passed to CHIRAL for the published concordance. It "
               "differs from 'Phase ordering' in two ways that are kept on purpose so the published "
               "values can be regenerated: BHLHE41 is absent, and BMAL1 is listed under its former "
               "symbol ARNTL and matched literally, so in matrices annotated with the current symbol "
               "it does not enter. Whether to resolve the synonym is a decision for the manuscript."),
    "Benchmark clock": dict(
        genes=["BMAL1", "CLOCK", "NPAS2", "NR1D1", "NR1D2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "DBP",
               "TEF", "HLF", "CIART", "NFIL3", "BHLHE40", "BHLHE41", "RORA", "RORC"],
        used_by="python/analyses/02e_gse87761_benchmark.py",
        source="Core clock and direct clock-output genes only, with no metabolic genes, used to benchmark "
               "the muscle template against recorded sampling times in adipose tissue."),
    "Fasting/FOXO3": dict(
        genes=["TXNIP", "PDK4", "ANGPTL4", "FBXO32", "TRIM63", "CDKN1A", "GADD45A", "SESN1",
               "BNIP3L", "KLF15", "DDIT4"],
        used_by="python/analyses/05_clamp_within_person.py",
        source="Fasting programme: FOXO and glucocorticoid targets induced by fasting in skeletal muscle, "
               "pre-specified."),
    "Glucose uptake": dict(
        genes=["HK2", "SLC2A4"],
        used_by="python/analyses/05_clamp_within_person.py",
        source="Insulin-stimulated glucose uptake and phosphorylation."),
    "Receptor targets, disjoint": dict(
        genes=["MT1X", "MT2A", "CEBPB", "TSC22D1", "SESN3", "KLF9", "DUSP1", "ERRFI1"],
        used_by="python/analyses/27_crosstissue_prediction.py, 02_meta_and_validation.py",
        source="Glucocorticoid-receptor targets that share no gene with the 'Glucocorticoid' set, used as "
               "the outcome of the cross-tissue prediction so that no gene sits on both sides."),
    "Nuclear receptors": dict(
        genes=["THRA", "THRB", "RARA", "RARB", "RARG", "PPARA", "PPARD", "PPARG", "NR1D1", "NR1D2", "RORA", "RORB", "RORC",
               "NR1H4", "NR1H3", "NR1H2", "VDR", "NR1I2", "NR1I3", "HNF4A", "HNF4G", "RXRA", "RXRB", "RXRG", "NR2C1", "NR2C2",
               "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6", "ESR1", "ESR2", "ESRRA", "ESRRB", "ESRRG", "NR3C1", "NR3C2", "PGR",
               "AR", "NR4A1", "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NR6A1", "NR0B1", "NR0B2"],
        used_by="python/analyses/30_receptor_specificity.py",
        source="The 48 members of the human nuclear receptor superfamily (HGNC gene group 'Nuclear hormone "
               "receptors'). Comparators for NR3C1: whether the association with clock output is a property "
               "of the glucocorticoid receptor or of nuclear receptors, or of expressed genes, in general. "
               "NR1D1 and NR1D2 belong to the clock-output set and are left out of the comparison."),
    "Model signature": dict(
        genes=["PER3", "DBP", "HLF", "TEF", "NR1D1", "NR1D2", "PER2", "CIART"],
        used_by="python/model/make_fig7_data.py, python/figures/supplementary_figures.py",
        source="Output-arm genes for which the forced-oscillator model generates synthetic expression "
               "and that are compared with the observed signature."),
}

def panel(name):
    """Genes de un panel de PANELS, en su orden."""
    return list(PANELS[name]["genes"])

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

def as_table(panels=False):
    """Tabla de los conjuntos puntuados. Con panels=True anade los paneles de PANELS."""
    import pandas as pd
    rows = [dict(gene_set=k, kind="scored set", n_genes=len(v["genes"]),
                 genes="; ".join(v["genes"]), source=v["source"]) for k, v in SETS.items()]
    if panels:
        rows += [dict(gene_set=k, kind="panel", n_genes=len(v["genes"]), genes="; ".join(v["genes"]),
                      source=v["source"], used_by=v["used_by"]) for k, v in PANELS.items()]
    d = pd.DataFrame(rows)
    return d if panels else d.drop(columns="kind")
