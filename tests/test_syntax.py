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
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
