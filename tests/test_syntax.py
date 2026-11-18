#!/usr/bin/env python3
"""Todos los scripts compilan, en la version de Python declarada en el entorno de CI.

Un f-string con una barra invertida dentro de la expresion compila en 3.12 y falla antes. Este
test lo detecta sin ejecutar nada, de modo que CI no se entera en el paso de figuras.
"""
import ast, glob, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    bad = []
    for f in sorted(glob.glob(f"{ROOT}/python/**/*.py", recursive=True) +
                    glob.glob(f"{ROOT}/tests/*.py")):
        if "_legacy" in f: continue
        try:
            ast.parse(open(f, encoding="utf-8").read(), filename=f)
        except SyntaxError as e:
            bad.append(f"{os.path.relpath(f, ROOT)}:{e.lineno}  {e.msg}")
    for b in bad: print("  ", b)
    v = f"{sys.version_info.major}.{sys.version_info.minor}"
    print(f"\n{'todos los scripts compilan' if not bad else str(len(bad)) + ' con error de sintaxis'}"
          f" (Python {v})")
    if sys.version_info >= (3, 12):
        print("  aviso: 3.12 admite barras invertidas dentro de expresiones de f-string y las "
              "versiones anteriores no.\n  CI corre 3.11: compruebe ahi antes de empujar.")
    # La cita: su version es la de la primera entrada del CHANGELOG y su fecha no esta en el futuro.
    import datetime, re
    cff = open(f"{ROOT}/CITATION.cff", encoding="utf-8").read(); log = open(f"{ROOT}/CHANGELOG.md", encoding="utf-8").read()
    ver = re.search(r"^version:\s*(\S+)", cff, re.M).group(1); day = re.search(r"^date-released:\s*(\S+)", cff, re.M).group(1)
    top = re.search(r"^## v(\S+?):", log, re.M).group(1)
    meta = []
    if ver != top: meta.append(f"CITATION.cff dice version {ver} y la primera entrada del CHANGELOG es v{top}")
    if datetime.date.fromisoformat(day) > datetime.datetime.now(datetime.timezone.utc).date() + datetime.timedelta(days=1):
        meta.append(f"CITATION.cff tiene una fecha de publicacion futura: {day}")
    # Los datos versionados con huella declarada: que sean lo que dicen ser.
    import hashlib
    for sums in sorted(glob.glob(f"{ROOT}/data/reference/**/SHA256SUMS", recursive=True)):
        for line in open(sums):
            h, name = line.split(); p = os.path.join(os.path.dirname(sums), name)
            if not os.path.exists(p): meta.append(f"falta {os.path.relpath(p, ROOT)}"); continue
            raw = open(p, "rb").read()
            if hashlib.sha256(raw).hexdigest() != h: meta.append(f"{os.path.relpath(p, ROOT)} no coincide con su huella")
            if name.endswith(".gz") and raw[:2] != b"\x1f\x8b": meta.append(f"{os.path.relpath(p, ROOT)} no es un gzip")
    for m in meta: print("  ", m)
    print("la cita coincide con el CHANGELOG, su fecha no es futura y los datos versionados conservan su huella" if not meta
          else f"{len(meta)} problema(s) en la cita o en los datos versionados")
    return 1 if (bad or meta) else 0

if __name__ == "__main__":
    sys.exit(main())
