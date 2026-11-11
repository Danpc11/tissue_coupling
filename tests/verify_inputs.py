#!/usr/bin/env python3
"""Agente 12: la cadena de dependencias esta cerrada.

Cada fichero que un script de figuras lee tiene que escribirlo otro script del repositorio, y
`run_all.sh` tiene que ejecutarlos en ese orden. Este test falla si alguna entrada aparece de la
nada, que es exactamente lo que rompia `run_all.sh` en un clon limpio.
"""
import os, re, sys, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
READ = re.compile(r'read_(?:csv|pickle)\(\s*f?["\'][^"\']*?(?:\{OUT\}|\{RAW\}|\{REF\})?/?([\w./\-]+\.(?:tsv|csv|pkl|npy|gz))')
WRITE = re.compile(r'to_(?:csv|pickle)\(\s*f?["\'][^"\']*?([\w./\-]+\.(?:tsv|csv|pkl))|np\.save\(\s*f?["\'][^"\']*?([\w./\-]+\.npy)')

def basenames(pattern, text):
    out = set()
    for m in pattern.finditer(text):
        for g in m.groups():
            if g: out.add(os.path.basename(g))
    return out

def main():
    reads, writes = {}, set()
    for f in glob.glob(f"{ROOT}/python/**/*.py", recursive=True):
        t = open(f).read()
        writes |= basenames(WRITE, t)
        if "_legacy" in f: continue
        if "/figures/" in f or "/analyses/" in f or "/model/" in f:
            for b in basenames(READ, t): reads.setdefault(b.replace('.gz',''), []).append(os.path.relpath(f, ROOT))
    # lo que esta versionado como dato de entrada no necesita productor
    versioned = {os.path.basename(p) for p in glob.glob(f"{ROOT}/data/reference/*")}
    # nombres que un script compone con una variable y escribe el mismo (p. ej. clamp_per_person_<contraste>)
    templated = {w.split("{")[0] for w in writes if "{" in w} | {"clamp_per_person_"}
    reads_templated = lambda b: any(b.startswith(t) for t in templated if t)
    raw = {os.path.basename(p) for p in glob.glob(f"{ROOT}/data/raw/**/*", recursive=True)}
    raw |= {b[:-3] for b in raw if b.endswith(".gz")}
    raw |= {b.replace('.','_') for b in raw} | {b.replace('_','.') for b in raw}
    # las matrices de expresion las produce 00_prepare_data.py desde los depositos publicos
    produced_by_prepare = lambda b: b.endswith(".pkl") and (b.startswith(("GTEx_", "GSE", "AAGMEx")))
    # descargas publicas que lista 00_prepare_data.py y el usuario coloca en data/raw
    external = lambda b: ("GRCh38" in b or b.startswith("GTEx_Analysis") or b.endswith(".gct")
                          or b.startswith("GSE") and "raw_counts" in b
                          or b in {".gct", "dorothea_ABC.tsv", "Homo_sapiens.gene_info",
                                   "GSE129843_RESTRICT.txt", "GSE182117_counts.tsv",
                                   "GSE66446_Cushings_htseq_gene_counts_GRCh37.74.txt"})
    missing = {b: v for b, v in reads.items()
               if not b.startswith(".") and b not in writes and b not in versioned and b not in raw
               and not produced_by_prepare(b) and not external(b)
               and not reads_templated(b)}
    for b, v in sorted(missing.items()):
        print(f"  SIN PRODUCTOR  {b}  (leido por {', '.join(sorted(set(v)))})")
    # orden: cada productor debe correr en run_all.sh antes que sus lectores
    order = [l for l in open(f"{ROOT}/run_all.sh") if l.startswith("step ")]
    pos = {}
    for i, l in enumerate(order):
        m = re.search(r"(python/[\w/]+\.py)", l)
        if m: pos[m.group(1)] = i
    writers = {}
    for f in glob.glob(f"{ROOT}/python/**/*.py", recursive=True):
        rel = os.path.relpath(f, ROOT)
        for b in basenames(WRITE, open(f).read()): writers.setdefault(b, []).append(rel)
    bad = 0
    for b, readers in reads.items():
        ws = [w for w in writers.get(b, []) if w in pos]
        for r in readers:
            if r in pos and ws and min(pos[w] for w in ws) > pos[r] and b not in versioned:
                print(f"  ORDEN  {b}: lo lee {r} (paso {pos[r]+1}) antes de que lo escriba {', '.join(ws)}"); bad += 1
    print(f"\n{len(reads)} entradas distintas, {len(writes)} producidas por el repositorio, "
          f"{len(versioned)} versionadas, {len(missing)} sin productor")
    return 1 if (missing or bad) else 0

if __name__ == "__main__":
    sys.exit(main())
