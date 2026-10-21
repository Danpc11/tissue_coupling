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
    raw = {os.path.basename(p) for p in glob.glob(f"{ROOT}/data/raw/*")}
    raw |= {b.replace('.','_') for b in raw} | {b.replace('_','.') for b in raw}
    # las matrices de expresion las produce 00_prepare_data.py desde los depositos publicos
    produced_by_prepare = lambda b: b.endswith(".pkl") and (b.startswith(("GTEx_", "GSE", "AAGMEx")))
    # descargas publicas que lista 00_prepare_data.py y el usuario coloca en data/raw
    external = lambda b: ("GRCh38" in b or b.startswith("GTEx_Analysis") or b.endswith(".gct")
                          or b.startswith("GSE") and "raw_counts" in b or b in {".gct"})
    missing = {b: v for b, v in reads.items()
               if b not in writes and b not in versioned and b not in raw
               and not produced_by_prepare(b) and not external(b)}
    for b, v in sorted(missing.items()):
        print(f"  SIN PRODUCTOR  {b}  (leido por {', '.join(sorted(set(v)))})")
    print(f"\n{len(reads)} entradas distintas, {len(writes)} producidas por el repositorio, "
          f"{len(versioned)} versionadas, {len(missing)} sin productor")
    return 1 if missing else 0

if __name__ == "__main__":
    sys.exit(main())
