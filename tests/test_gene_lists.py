#!/usr/bin/env python3
"""La regla del registro de genes, comprobada sobre el codigo.

python/lib/gene_sets.py es el unico lugar donde se define un conjunto de genes que determina una
puntuacion, una inferencia o un control. Este test recorre el arbol sintactico de cada script y
falla si encuentra, fuera del registro, una coleccion literal con cinco o mas simbolos de genes
del registro que no este declarada aqui abajo con su motivo.

Quedan permitidas dos clases de listas, y solo esas:
  - en python/figures/, las que eligen que genes dibuja un panel;
  - las de DECLARED: tablas de parametros por gen y selecciones de filas de una tabla, una a una.
"""
import ast, glob, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python/lib"))
from gene_sets import SETS, PANELS, SYNONYMS

MIN_GENES = 5
DECLARED = {
    ("python/model/forced_clock.py", "genes_phase"):
        "model parameter: peak phase of each gene in the synthetic expression",
    ("python/analyses/05_clamp_within_person.py", "CLOCK_TEMPLATE"):
        "estimator parameter: reference phase and amplitude of each clock gene when no template is available",
    ("python/analyses/01_crosstissue_axis.py", "PROBES_U133"):
        "platform table: the U133 probe set that measures each gene in GSE41168",
    ("python/analyses/00d_reference_tables.py", "genes"):
        "row selection: the genes listed in three versioned tables that the step regenerates",
}

universe = {g for v in list(SETS.values()) + list(PANELS.values()) for g in v["genes"]}
universe |= set(SYNONYMS) | {a for alts in SYNONYMS.values() for a in alts}

def literal_genes(node):
    """Simbolos del registro que aparecen como constantes en una lista, tupla, conjunto o dict."""
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        items = node.elts
    elif isinstance(node, ast.Dict):
        items = [k for k in node.keys if k is not None]
    else:
        return 0
    return sum(isinstance(e, ast.Constant) and isinstance(e.value, str) and e.value in universe for e in items)

def names_of(tree):
    """Nombre de la variable a la que se asigna cada nodo literal (si es una asignacion simple)."""
    owner = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    for sub in ast.walk(node.value):
                        owner.setdefault(id(sub), t.id)
    return owner

found, undeclared, used = [], [], set()
for path in sorted(glob.glob(os.path.join(ROOT, "python", "**", "*.py"), recursive=True)):
    rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
    if rel == "python/lib/gene_sets.py":
        continue
    tree = ast.parse(open(path, encoding="utf-8").read())
    owner = names_of(tree)
    for node in ast.walk(tree):
        n = literal_genes(node)
        if n < MIN_GENES:
            continue
        name = owner.get(id(node), "<expresion>")
        found.append((rel, node.lineno, name, n))
        if rel.startswith("python/figures/"):
            continue
        if (rel, name) in DECLARED:
            used.add((rel, name)); continue
        undeclared.append((rel, node.lineno, name, n))

print(f"colecciones literales con {MIN_GENES} o mas genes del registro, fuera del registro: {len(found)}")
for rel, line, name, n in found:
    kind = ("panel de figura" if rel.startswith("python/figures/") else
            "declarada" if (rel, name) in DECLARED else "SIN DECLARAR")
    print(f"  {kind:16s} {rel}:{line}  {name}  ({n} genes)")
stale = sorted(set(DECLARED) - used)
if stale:
    print("declaraciones que ya no corresponden a nada en el codigo:", stale)
if undeclared or stale:
    print(f"\nFALLA: {len(undeclared)} sin declarar, {len(stale)} declaraciones sobrantes. "
          "Un conjunto que define una puntuacion va en python/lib/gene_sets.py.")
    sys.exit(1)
print("\nTODO CORRECTO")
