#!/usr/bin/env python3
"""Construye las matrices intermedias `results/new/*.pkl` que leen los pasos 02-27.

Reconstruccion (Hugo Tovar, nov 2026): la version publicada del repositorio no incluia el
paso que convierte las descargas publicas en estas matrices. Las convenciones de aqui se dedujeron
de lo que cada script consumidor espera (indice = simbolo actual, columnas = GSM o donante, log2)
y de Methods "Cohorts and preprocessing":

  - microarrays anotados por Entrez -> simbolo actual (lib/gpl570.annotate); sonda exclusiva de
    mayor media por gen, sonda multi-gen solo si no hay exclusiva
  - RNA-seq como cuentas crudas -> log2 CPM, genes con > 5 cuentas en >= 10 muestras
  - los metadatos de cada serie van en `<acc>_meta.pkl` con columnas gsm, title y c0..cN, una por
    linea `!Sample_characteristics_ch1` en el orden del fichero, mas columnas con nombre cuando un
    consumidor las pide (fh, stim, status, dx, lab, disease)

Entradas sustitutas (ver data/raw/README_descarga.md): GSE108539 desde recount3 (SRP127581);
GSE159984 desde ARCHS4 (00c_archs4_counts.py) porque la web de NCBI no es accesible desde el
servidor.

Uso:  python python/analyses/00b_build_matrices.py [cohorte ...]
"""
import os, sys, gzip, glob, re, io, warnings, importlib.util
warnings.filterwarnings("ignore", category=FutureWarning); warnings.filterwarnings("ignore", category=DeprecationWarning)
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../lib"))
import gpl570

RAW = os.environ.get("T2D_RAW", "data/raw")
OUT = os.environ.get("T2D_OUT", "results")
REF = os.environ.get("T2D_REF", "data/reference")
NEW = f"{OUT}/new"
os.makedirs(NEW, exist_ok=True)


# ----------------------------------------------------------------------------- utilidades
def read_series(acc):
    """Matriz de serie de GEO -> (metadatos, expresion o None)."""
    meta, chars, tab, intab = {}, [], [], False
    with gzip.open(f"{RAW}/{acc}_series_matrix.txt.gz", "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!series_matrix_table_begin"):
                intab = True; continue
            if l.startswith("!series_matrix_table_end"):
                break
            if intab:
                tab.append(l.rstrip("\n").split("\t")); continue
            v = [x.strip('"') for x in l.rstrip("\n").split("\t")[1:]]
            if l.startswith("!Sample_geo_accession"): meta["gsm"] = v
            elif l.startswith("!Sample_title"): meta["title"] = v
            elif l.startswith("!Sample_characteristics_ch1"): chars.append(v)
            elif l.startswith("!Sample_description") and "desc" not in meta: meta["desc"] = v
    for i, v in enumerate(chars):
        meta[f"c{i}"] = v
    p = pd.DataFrame(meta)
    if len(tab) > 1:
        e = pd.DataFrame(tab[1:], columns=[c.strip('"') for c in tab[0]]).set_index(tab[0][0].strip('"'))
        e.index = [i.strip('"') for i in e.index]
        e = e.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    else:
        e = None
    return p, e


def val(s):
    """'clave: valor' -> 'valor'."""
    return s.str.split(": ", n=1).str[1].str.strip()


def log2_if_needed(E):
    return np.log2(E.clip(lower=1)) if np.nanmax(E.values) > 60 else E


def log2cpm(C, min_count=5, min_samples=10):
    """log2 CPM con el tamano de libreria calculado antes del filtro de genes."""
    lib = C.sum(); C = C[(C > min_count).sum(1) >= min_samples]
    return np.log2(C / lib * 1e6 + 1)


_ENS = None
def ensembl_to_symbol():
    """ENSG -> simbolo actual, desde dbXrefs de gene_info (NCBI)."""
    global _ENS
    if _ENS is None:
        fn = f"{RAW}/_soft/Homo_sapiens.gene_info.gz"
        g = pd.read_csv(fn, sep="\t", usecols=["Symbol", "dbXrefs"], dtype=str)
        m = {}
        for sym, x in zip(g.Symbol, g.dbXrefs.fillna("")):
            for tok in x.split("|"):
                if tok.startswith("Ensembl:"):
                    m.setdefault(tok[8:], sym)
        _ENS = m
    return _ENS


def by_symbol_from_ensembl(C):
    m = ensembl_to_symbol()
    idx = pd.Index([str(i).split(".")[0] for i in C.index])
    C = C.copy(); C.index = [m.get(i, None) for i in idx]
    C = C[pd.notna(C.index)]
    return C.groupby(level=0).sum()


def synonym_map():
    """results/new/ncbi_synonym_map.tsv: sinonimo (que no es simbolo vigente) -> simbolo vigente,
    desde gene_info de NCBI. Lo leen aagmex(), 09_hpo_by_tissue.py y 00d_reference_tables.py."""
    out = f"{NEW}/ncbi_synonym_map.tsv"
    g = pd.read_csv(f"{RAW}/_soft/Homo_sapiens.gene_info.gz", sep="\t", usecols=["Symbol", "Synonyms"], dtype=str)
    cur = set(g.Symbol); m = {}
    for sym, syn in zip(g.Symbol, g.Synonyms.fillna("-")):
        for x in syn.split("|"):
            if x and x != "-" and x not in cur and x not in m: m[x] = sym
    D = pd.Series(m, name="current").to_frame(); D.index.name = None
    D.to_csv(out, sep="\t"); print(f"  -> {out}  ({len(D)} sinonimos)")
    return D["current"].to_dict()


def collapse(expr, M):
    """Misma regla que gpl570.annotate, con un mapa sonda->simbolos ya construido
    (columnas ID, n_genes, new)."""
    M = M[(M.new != "") & M.ID.isin(expr.index)]
    mean = expr.mean(1); best = {}
    for multi in (False, True):
        sub = M[(M.n_genes > 1) == multi]
        for pid, new in zip(sub.ID, sub.new):
            for gene in new.split(";"):
                if gene in best and not multi:
                    if mean[pid] <= mean[best[gene]]: continue
                elif gene in best:
                    continue
                best[gene] = pid
    genes = sorted(best)
    out = expr.loc[[best[g] for g in genes]].copy(); out.index = genes
    return out


def c5():
    """Modulo 05_clamp_within_person (lector y anotacion por simbolo de plataforma)."""
    sp = importlib.util.spec_from_file_location("c5", os.path.join(HERE, "05_clamp_within_person.py"))
    m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m


def save(obj, name):
    obj.to_pickle(f"{NEW}/{name}.pkl")
    shape = getattr(obj, "shape", None)
    print(f"  -> {NEW}/{name}.pkl  {shape}")


# ----------------------------------------------------------------------------- cohortes
def aagmex():
    """GSE95674 adiposo / GSE95675 musculo, Illumina HT-12 v4 cuantil-normalizado, no log.
    Columnas = identificador de donante (titulo; el de musculo lleva sufijo -M)."""
    fn = f"{NEW}/ncbi_synonym_map.tsv"
    syn = pd.read_csv(fn, sep="\t", index_col=0)["current"].to_dict() if os.path.exists(fn) else synonym_map()
    mats = {}
    for acc, tis in (("GSE95674", "adipose"), ("GSE95675", "muscle")):
        p, e = read_series(acc)
        E = log2_if_needed(e)
        E.columns = [t.replace("-M", "") for t in p.set_index("gsm").title.reindex(E.columns)]
        # simbolo de plataforma (2010) -> actual cuando el actual no tiene fila propia
        have = set(E.index)
        E.index = [syn[g] if (g in syn and syn[g] not in have) else g for g in E.index]
        E = E.groupby(level=0).mean()
        print(f"{acc} {tis}: {E.shape[1]} muestras, {E.shape[0]} genes, rango {E.values.min():.1f}-{E.values.max():.1f}")
        mats[tis] = E
        p["donor"] = p.title.str.replace("-M", "", regex=False)
        save(p, f"AAGMEx_{tis}_meta")
    # 06/07/09 alinean las covariables por posicion entre adiposo y musculo: las dos matrices
    # llevan los mismos donantes (los 232 con ambos tejidos) en el mismo orden
    common = [d for d in mats["adipose"].columns if d in set(mats["muscle"].columns)]
    print(f"AAGMEx: {len(common)} donantes con ambos tejidos")
    for tis, E in mats.items():
        save(E[common], f"AAGMEx_{tis}_log2")


def gse41168():
    """GPL570, RMA log2. Musculo y adiposo de las mismas mujeres, muestras basales ('before').
    Columnas = identificador de sujeto."""
    p, e = read_series("GSE41168")
    E = gpl570.annotate(log2_if_needed(e))
    p["tissue"] = val(p.c0); p["group"] = val(p.c1); p["treatment"] = val(p.c2); p["sex"] = val(p.c3)
    p["subject"] = p.title.str.split(",").str[0].str.split(" ").str[0].str.strip()
    save(p, "GSE41168_meta")
    for tis, key in (("adipose", "adipose"), ("muscle", "muscle")):
        q = p[(p.treatment == "before") & p.tissue.str.contains(key)].drop_duplicates("subject")
        X = E[list(q.gsm)]; X.columns = list(q.subject)
        print(f"GSE41168 {tis}: {X.shape[1]} sujetos basales")
        save(X, f"GSE41168_{tis}")


def gpl570_case_control():
    """GSE25462 (MAS5, no log) y GSE26637 (GC-RMA, log2). Dos anotaciones: por Entrez
    (`_entrez`) y por simbolo de plataforma como en 05 (`_log2`)."""
    m5 = c5()
    for acc in ("GSE25462", "GSE26637"):
        p, e = read_series(acc)
        E = log2_if_needed(e)
        if acc == "GSE25462":
            p["sex"] = val(p.c1); p["age"] = pd.to_numeric(val(p.c2), errors="coerce")
            p["fh"] = val(p.c3)
            sicol = next(c for c in p.columns if c.startswith("c") and p[c].astype(str).str.startswith("si:").all())
            p["si"] = pd.to_numeric(val(p[sicol]), errors="coerce")
        else:
            p["sex"] = val(p.c0); p["stim"] = val(p.c2); p["status"] = val(p.c3)
        save(p, f"{acc}_meta")
        Xe = gpl570.annotate(E); save(Xe, f"{acc}_entrez")
        Xs = m5.annotate_gpl570(E); save(Xs, f"{acc}_log2")
        print(f"{acc}: {E.shape[1]} muestras; entrez {Xe.shape[0]} genes, simbolo {Xs.shape[0]} genes")


def gse78721():
    """PrimeView (GPL15207), RMA log2; 130 muestras de adiposo (limb, subcutaneo, visceral)."""
    p, e = read_series("GSE78721")
    p["sex"] = val(p.c0); p["dx"] = val(p.c1); p["depot"] = val(p.c2).str.replace("Sudcutaneous", "Subcutaneous")
    # 130 muestras de 95 personas: las 35 de cirugia abdominal aportan subcutaneo y visceral (sin identificador
    # de persona en GEO). `one_per_person` marca una muestra por persona: extremidad o subcutaneo.
    p["one_per_person"] = p.depot.isin(["Limb", "Subcutaneous"])
    E = gpl570.annotate(log2_if_needed(e), path=f"{RAW}/GPL15207_annot.txt")
    save(p, "GSE78721_meta"); save(E, "GSE78721_log2")
    print(f"GSE78721: {E.shape}, dx {p.dx.value_counts().to_dict()}")


def gse182120():
    """HTA 2.0 (GPL17586), RMA log2 por transcript cluster. Sonda -> Entrez por la columna
    gene_assignment de la tabla de plataforma, y de ahi al simbolo actual."""
    p, e = read_series("GSE182120")
    p["disease"] = val(p.c0); p["age"] = pd.to_numeric(val(p.c1), errors="coerce")
    p["bmi"] = pd.to_numeric(val(p.c2), errors="coerce"); p["mvalue"] = pd.to_numeric(val(p.c3), errors="coerce")
    p["lab"] = val(p.c4)
    g = pd.read_csv(f"{RAW}/GPL17586_annot.txt", sep="\t", usecols=["ID", "gene_assignment"], dtype=str)
    ann = pd.read_csv(f"{RAW}/Human_GRCh38_p13_annot_tsv.gz", sep="\t", usecols=["GeneID", "Symbol"], dtype=str)
    cur = dict(zip(ann.GeneID, ann.Symbol))
    rows = []
    for pid, ga in zip(g.ID, g.gene_assignment.fillna("")):
        ents = set()
        for tr in ga.split("///"):
            f = [x.strip() for x in tr.split("//")]
            if len(f) >= 5 and f[4].isdigit(): ents.add(f[4])
        new = sorted({cur[x] for x in ents if x in cur})
        rows.append(dict(ID=pid, n_genes=len(ents), new=";".join(new)))
    M = pd.DataFrame(rows)
    E = collapse(log2_if_needed(e), M)
    save(p, "GSE182120_meta"); save(E, "GSE182120_log2")
    print(f"GSE182120: {E.shape}, lab {p.lab.value_counts().to_dict()}")


def gse182117():
    """Miotubos sincronizados: cuentas del autor (Ensembl, columnas = titulo) -> log2 CPM."""
    p, _ = read_series("GSE182117")
    C = pd.read_csv(f"{RAW}/GSE182117_counts.tsv.gz", sep="\t", index_col=0)
    t2g = dict(zip(p.title, p.gsm))
    C = C[[c for c in C.columns if c in t2g]]; C.columns = [t2g[c] for c in C.columns]
    X = log2cpm(by_symbol_from_ensembl(C))
    save(p, "GSE182117_meta"); save(X, "GSE182117_log2cpm")
    print(f"GSE182117: {X.shape}")


def gse129843():
    """Musculo in vivo, 11 hombres x 2 condiciones x 6 horas: CPM del autor (RESTRICT.txt),
    columnas Sample.ID -> GSM por el README de la serie."""
    p, _ = read_series("GSE129843")
    rd = pd.read_csv(f"{RAW}/GSE129843_README.txt", sep="\t", skiprows=4).dropna()
    s2g = dict(zip(rd.iloc[:, 0].str.strip(), rd.iloc[:, 1].str.strip()))
    d = pd.read_csv(f"{RAW}/GSE129843_RESTRICT.txt.gz", sep="\t")
    # simbolo actual por Ensembl (el del autor es de 2018: ARNTL, no BMAL1); el del autor de respaldo
    m = ensembl_to_symbol()
    auth = d.Symbol.where(d.Symbol.notna() & ~d.Symbol.str.startswith("ENSG"), None)
    sym = [m.get(g, None) or (s if s else None) for s, g in zip(auth, d.genes)]
    C = d.drop(columns=["genes", "Symbol", "entrez"]); C.index = sym
    C = C[pd.notna(C.index)].groupby(level=0).sum()
    C = C[[c for c in C.columns if c in s2g]]; C.columns = [s2g[c] for c in C.columns]
    C = C[(C > 1).sum(1) >= 10]
    X = np.log2(C + 1)
    save(p, "GSE129843_meta"); save(X, "GSE129843_log2cpm")
    print(f"GSE129843: {X.shape}")


def recount3_counts(srp):
    """Sumas por gen de recount3 (SRP) agregadas por GSM, Ensembl -> simbolo actual."""
    d = f"{RAW}/recount3"
    md = pd.read_csv(f"{d}/sra.sra.{srp}.MD.gz", sep="\t", usecols=["external_id", "sample_name"])
    C = pd.read_csv(f"{d}/sra.gene_sums.{srp}.G026.gz", sep="\t", comment="#", index_col=0)
    r2g = dict(zip(md.external_id, md.sample_name))
    C = C[[c for c in C.columns if c in r2g]]
    # las gene_sums de recount3 son cobertura por base: se dividen por la longitud media de lectura
    # (star.average_input_read_length de recount_qc) para dejarlas en lecturas
    qc = pd.read_csv(f"{d}/sra.recount_qc.{srp}.MD.gz", sep="\t")
    rl = dict(zip(qc.external_id, pd.to_numeric(qc["star.average_input_read_length"], errors="coerce")))
    C = C / np.array([rl.get(c, np.nan) for c in C.columns])
    C = C.T.groupby([r2g[c] for c in C.columns]).sum().T.round()
    return by_symbol_from_ensembl(C)


def gse66446():
    """Cushing adiposo -> log2 CPM. Preferencia: recount3 (SRP055749, GRCh38), cuya firma
    correlaciona 0.98 con la columna versionada de signatures_reference.tsv; las cuentas htseq
    del autor (GRCh37) solo llegan a 0.91 y hacen abortar 23_build_signatures_table.py."""
    p, _ = read_series("GSE66446")
    if os.path.exists(f"{RAW}/recount3/sra.gene_sums.SRP055749.G026.gz"):
        X = log2cpm(recount3_counts("SRP055749")); print("GSE66446: recount3")
    else:
        C = pd.read_csv(f"{RAW}/GSE66446_Cushings_htseq_gene_counts_GRCh37.74.txt.gz", sep="\t", index_col=0)
        C = C[~C.index.str.startswith("__")]
        t2g = {t.split(" ")[0].lower(): g for t, g in zip(p.title, p.gsm)}
        C.columns = [t2g[c.lower()] for c in C.columns]
        X = log2cpm(by_symbol_from_ensembl(C)); print("GSE66446: cuentas del autor (GRCh37)")
    save(p, "GSE66446_meta"); save(X, "GSE66446_log2cpm")
    print(f"GSE66446: {X.shape}, {p.c0.value_counts().to_dict()}")


def gse13070():
    """GPL570 anotado por Entrez (la version por simbolo la cachea 05 en results/)."""
    p, e = read_series("GSE13070")
    E = gpl570.annotate(log2_if_needed(e))
    save(E, "GSE13070_entrez")
    print(f"GSE13070: {E.shape}")


def gtex():
    """GTEx v11: cuentas por gen -> un donante por columna -> log2 CPM. Fenotipos de la tabla
    de sujetos (SEX, AGE, DTHHRDY)."""
    ph = pd.read_csv(f"{RAW}/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt", sep="\t").set_index("SUBJID")
    PH = ph[["SEX", "AGE", "DTHHRDY"]]
    save(PH, "GTEx_pheno")
    SA = pd.read_csv(f"{RAW}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt", sep="\t", low_memory=False)
    SA["donor"] = SA.SAMPID.str.split("-").str[:2].str.join("-")
    for fn, name, smtsd in (("gene_reads_v11_muscle_skeletal.gct.gz", "GTEx_muscle", "Muscle - Skeletal"),
                            ("gene_reads_v11_adipose_subcutaneous.gct.gz", "GTEx_adipose_sc", "Adipose - Subcutaneous"),
                            ("gene_reads_v11_adipose_visceral_omentum.gct.gz", "GTEx_adipose_visc", "Adipose - Visceral (Omentum)")):
        G = pd.read_csv(f"{RAW}/{fn}", sep="\t", skiprows=2)
        # simbolo vigente por Ensembl (gene_info); si no hay, el de GENCODE de la GCT
        ens = ensembl_to_symbol(); G["Description"] = [ens.get(n.split(".")[0], s) for n, s in zip(G.Name, G.Description)]
        C = G.drop(columns=["Name"]).set_index("Description")
        C = C.groupby(level=0).sum()
        # una muestra por donante: la primera que SampleAttributes lista para ese tejido (misma regla que 08)
        # una muestra por donante: la primera que SampleAttributes lista para ese tejido y que esta en la GCT
        sa = SA[SA.SMTSD == smtsd]; sa = sa[sa.SAMPID.isin(C.columns)]
        first = sa.drop_duplicates("donor").set_index("donor").SAMPID
        keep = [first[d] for d in first.index if d in PH.index]
        C = C[keep]; C.columns = ["-".join(c.split("-")[:2]) for c in keep]
        pd.Series(keep, index=C.columns, name="SAMPID").to_csv(f"{NEW}/{name}_samples.tsv", sep="\t")
        X = log2cpm(C)
        save(X, name)
        print(f"{name}: {X.shape}")


def gse108539_recount3():
    """GSE108539 desde recount3 (SRP127581): suma de SRR por GSM, Ensembl -> simbolo.
    Escribe data/raw/GSE108539_counts.tsv.gz, que es el nombre alternativo que lee 02b."""
    C = recount3_counts("SRP127581").round().astype(int)
    C.to_csv(f"{RAW}/GSE108539_counts.tsv.gz", sep="\t")
    print(f"GSE108539 (recount3): {C.shape} -> {RAW}/GSE108539_counts.tsv.gz")


def gse148640():
    """Hidrocortisona frente a placebo en adiposo (HuGene 2.0 ST, n = 10, pareado): diferencia
    por persona (tratamiento - placebo, log2) de los valores procesados del autor, que ya vienen
    centrados por pareja. results/GSE148640_diff_genes.pkl, que lee la Fig. 3d-e."""
    d = pd.read_csv(f"{RAW}/GSE148640_processed.txt.gz", sep="\t", skiprows=9, low_memory=False, dtype=str)
    d = d.rename(columns={d.columns[1]: "gene"}); d = d[d.gene.notna() & (d.gene != "---")]
    hdr = pd.read_csv(f"{RAW}/GSE148640_processed.txt.gz", sep="\t", header=None, nrows=9, dtype=str)
    trt = hdr[hdr[4] == "treatment"].iloc[0, 5:25].str.lower().values
    pats = pd.read_csv(f"{RAW}/GSE148640_processed.txt.gz", sep="\t", skiprows=9, nrows=0, dtype=str).columns[5:25]
    pats = [c.split(".")[0] for c in pats]
    vals = d.iloc[:, 5:25].apply(pd.to_numeric, errors="coerce"); vals.columns = range(20)
    pla = [i for i in range(20) if trt[i].startswith("placebo")]; hc = [i for i in range(20) if not trt[i].startswith("placebo")]
    assert all(pats[p] == pats[h] for p, h in zip(pla, hc)), "placebo e hidrocortisona no emparejan por paciente"
    diff = pd.DataFrame({pats[h]: vals[h].values - vals[p].values for p, h in zip(pla, hc)}, index=d.gene.values)
    diff = diff.groupby(level=0).mean()
    diff.to_pickle(f"{OUT}/GSE148640_diff_genes.pkl")
    print(f"GSE148640: {diff.shape} -> {OUT}/GSE148640_diff_genes.pkl; PER1 sube en {(diff.loc['PER1']>0).sum()}/{diff.shape[1]}")


STEPS = dict(synonyms=synonym_map, aagmex=aagmex, gse41168=gse41168, gpl570=gpl570_case_control, gse78721=gse78721,
             gse182120=gse182120, gse182117=gse182117, gse129843=gse129843, gse66446=gse66446,
             gse13070=gse13070, gse148640=gse148640, gtex=gtex, gse108539=gse108539_recount3)

if __name__ == "__main__":
    want = sys.argv[1:] or list(STEPS)
    failed = []
    for k in want:
        print(f"\n=== {k}")
        try:
            STEPS[k]()
        except Exception as ex:
            import traceback; traceback.print_exc(); print(f"  FALLO {k}: {type(ex).__name__}: {ex}"); failed.append(k)
    if failed:
        sys.exit(f"pasos fallidos: {', '.join(failed)}")
