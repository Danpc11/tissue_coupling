#!/usr/bin/env python3
"""La cadena de dependencias esta cerrada, hasta donde puede verse sin ejecutar nada.

Cada fichero que un script lee tiene que escribirlo otro script del repositorio, estar versionado
en data/reference/ o ser una descarga publica declarada, y `run_all.sh` tiene que ejecutar al
productor antes que al lector. Este test falla si alguna entrada aparece de la nada, que es lo que
rompia `run_all.sh` en un clon limpio.

ALCANCE. Es un analisis ESTATICO del codigo, y conviene saber que no ve:
  - reconoce lecturas con read_csv / read_pickle y escrituras con to_csv / to_pickle / np.save
    cuya ruta es un literal o un f-string; no ve las que pasan por otra funcion o construyen la
    ruta en una variable (las matrices que escribe save() en 00b_build_matrices.py si las resuelve,
    ver saved_by_helper);
  - compara nombres de fichero, no rutas completas: dos ficheros con el mismo nombre en directorios
    distintos se confunden;
  - no abre ningun fichero: que exista un productor no dice que produzca el contenido correcto.
Lo que cubre esos huecos son las pruebas que SI ejecutan codigo: tests/test_clamp_vs_saline.py
(un paso completo sobre una serie simulada con respuesta conocida), tests/smoke_test.py (el
estimador de fase) y tests/verify_package.py (aritmetica y coherencia de las tablas de una
corrida real, incluido que todas procedan de la misma ejecucion).
"""
import ast, os, re, sys, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
READ = re.compile(r'read_(?:csv|pickle)\(\s*f?["\'][^"\']*?(?:\{OUT\}|\{RAW\}|\{REF\})?/?([\w./\-]+\.(?:tsv|csv|pkl|npy|gz))')
WRITE = re.compile(r'to_(?:csv|pickle)\(\s*f?["\'][^"\']*?([\w./\-]+\.(?:tsv|csv|pkl))|np\.save\(\s*f?["\'][^"\']*?([\w./\-]+\.npy)')

def basenames(pattern, text):
    out = set()
    for m in pattern.finditer(text):
        for g in m.groups():
            if g: out.add(os.path.basename(g))
    return out

def saved_by_helper():
    """Nombres de las matrices que escribe una funcion auxiliar del tipo `save(obj, nombre)`.

    Las expresiones regulares no ven estas escrituras porque el nombre se pasa como argumento.
    Aqui se leen del arbol sintactico: un nombre literal se toma tal cual; en un f-string o una
    variable, cada parte se sustituye por los valores del bucle `for` que la define cuando ese
    bucle recorre una coleccion literal, y por un comodin cuando no puede resolverse.
    Devuelve ({nombre: script que lo escribe}, patrones con comodin).
    """
    exact, patterns = {}, []
    for f in glob.glob(f"{ROOT}/python/**/*.py", recursive=True):
        src = open(f, encoding="utf-8").read()
        if "def save(" not in src or "to_pickle" not in src:
            continue
        tree = ast.parse(src)
        for fn in [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]:
            # valores posibles de cada variable de bucle con iterable literal, dentro de la funcion
            values = {}
            for loop in [n for n in ast.walk(fn) if isinstance(n, ast.For)]:
                if not isinstance(loop.iter, (ast.Tuple, ast.List)):
                    continue
                targets = loop.target.elts if isinstance(loop.target, ast.Tuple) else [loop.target]
                for item in loop.iter.elts:
                    parts = item.elts if isinstance(item, (ast.Tuple, ast.List)) and len(targets) > 1 else [item]
                    for t, v in zip(targets, parts):
                        if isinstance(t, ast.Name) and isinstance(v, ast.Constant) and isinstance(v.value, str):
                            values.setdefault(t.id, set()).add(v.value)
            def expand(node):
                """Lista de cadenas posibles; None en una posicion significa comodin."""
                if isinstance(node, ast.Constant) and isinstance(node.value, str):
                    return [node.value]
                if isinstance(node, ast.Name):
                    return sorted(values[node.id]) if node.id in values else [None]
                if isinstance(node, ast.JoinedStr):
                    outs = [""]
                    for part in node.values:
                        opts = expand(part.value) if isinstance(part, ast.FormattedValue) else expand(part)
                        outs = [None if (o is None or x is None) else o + x for o in outs for x in opts] \
                            if None not in opts else [None]
                        if outs == [None]:
                            # no resoluble: se guarda como patron con comodin
                            pat = "".join(re.escape(p.value) if isinstance(p, ast.Constant) else ".+" for p in node.values)
                            return [("pattern", pat)]
                    return outs
                return [None]
            for call in [n for n in ast.walk(fn) if isinstance(n, ast.Call)]:
                if getattr(call.func, "id", None) == "save" and len(call.args) >= 2:
                    for name in expand(call.args[1]):
                        if isinstance(name, tuple):
                            patterns.append(re.compile(name[1] + r"\.pkl$"))
                        elif name is not None:
                            exact[name + ".pkl"] = os.path.relpath(f, ROOT)
    return exact, patterns

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
    # las matrices de expresion: solo cuentan las que un save() del repositorio escribe de verdad,
    # no cualquier fichero cuyo nombre empiece por GSE, GTEx o AAGMEx
    saved, saved_patterns = saved_by_helper()
    writes |= set(saved)
    produced_by_prepare = lambda b: b in saved or any(p.fullmatch(b) for p in saved_patterns)
    # descargas publicas que lista 00_prepare_data.py y el usuario coloca en data/raw
    external = lambda b: ("GRCh38" in b or b.startswith("GTEx_Analysis") or b.endswith(".gct")
                          or b.startswith("GSE") and "raw_counts" in b
                          or b in {".gct", "dorothea_ABC.tsv", "Homo_sapiens.gene_info",
                                   "GSE129843_RESTRICT.txt", "GSE148640_processed.txt", "GSE182117_counts.tsv",
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
    for b, rel in saved.items(): writers.setdefault(b, []).append(rel)
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
