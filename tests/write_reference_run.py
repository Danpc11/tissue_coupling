#!/usr/bin/env python3
"""Escribe docs/REFERENCE_RUN.md a partir del manifiesto de una corrida terminada.

El documento decia de que corrida salen las cifras del articulo y se mantenia a mano: quedo atras
(seguia nombrando la corrida de la version 2.13.2 cuando las tablas eran de otra). Ahora se escribe
desde `run_manifest.json`, de modo que no puede decir otra cosa que lo que la corrida registro.

Uso, tras una corrida completa con `run_all.sh`:
    python tests/write_reference_run.py [directorio de resultados]
Se niega a escribir si la corrida no esta cerrada, no paso la verificacion o no es la configuracion
de referencia: una corrida asi no puede ser la que respalda las cifras.
"""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("T2D_OUT", os.path.join(ROOT, "results"))


def text(M):
    c = M.get("components", {}); isl = M.get("islet_counts", {})
    what = ("every result of the release, the external benchmarks included" if M.get("reproduces_everything")
            else "the main results of the release; the external benchmarks listed below as not computed are missing")
    rows = [("release", M["release"]), ("commit", f"`{M['git_commit']}`"), ("run identifier", f"`{M['run_id']}`"),
            ("profile", f"{M['profile']} (default settings)"), ("started", M.get("written", "")), ("finished", M.get("finished", "")),
            ("duration", f"{M.get('duration_minutes', '?')} minutes"),
            ("input files", f"{M['input_files']} (sha256 of their hashes `{M['input_files_sha256']}`)"),
            ("output files", f"{M.get('output_files', '?')} (sha256 of their hashes `{M.get('output_files_sha256', '?')}`)"),
            ("islet counts", str(isl.get("source", "?"))),
            ("uncommitted changes in versioned files", "yes" if M.get("git_uncommitted_changes") else "none"),
            ("CHIRAL (step 06)", c.get("chiral", "?")), ("ZeitZeiger and tauFisher (step 02f)", c.get("supervised_benchmark", "?")),
            ("final check", f"`verify_package.py` clean; the run reproduces {what}")]
    return ("# The run the reported numbers come from\n\n"
            "This file is written by `tests/write_reference_run.py` from the `run_manifest.json` of the run; it is not\n"
            "edited by hand. The numbers of the manuscript, its figures and its Supplementary Tables are those of\n"
            "this one complete run of the pipeline, from the public downloads to the last check.\n\n"
            "| | |\n|---|---|\n" + "".join(f"| {k} | {v} |\n" for k, v in rows) +
            "\nThe record of the run is its `run_manifest.json`, with `input_files.tsv` and `output_files.tsv`: the\n"
            "sha256 of every input and of every output. They are deposited with the archived release, so that a\n"
            "later reproduction can be compared file by file.\n\n"
            "The history of this repository was rewritten on 17 November 2026 (see the entry of v2.13.2 in\n"
            "`CHANGELOG.md`); commits are identified as they are after the rewriting.\n")


if __name__ == "__main__":
    fm = os.path.join(OUT, "new", "run_manifest.json")
    if not os.path.exists(fm): sys.exit(f"no hay manifiesto en {fm}: corra antes run_all.sh")
    M = json.load(open(fm))
    why = [w for ok, w in ((M.get("verification_passed"), "la verificacion no paso o la corrida no se cerro (paso 31)"),
                           (M.get("reference_configuration"), f"no es la configuracion de referencia (perfil {M.get('profile')})"),
                           (M.get("reproduces_main_results"), "los insumos de islotes no son los de las cifras publicadas"),
                           (not M.get("git_uncommitted_changes"), "habia cambios sin confirmar en ficheros versionados")) if not ok]
    if why: sys.exit("esta corrida no puede ser la de referencia: " + "; ".join(why))
    dst = os.path.join(ROOT, "docs", "REFERENCE_RUN.md")
    open(dst, "w", encoding="utf-8").write(text(M))
    print(f"escrito {os.path.relpath(dst, ROOT)}: version {M['release']}, corrida {M['run_id']}")
    print("deposite con la version: " + ", ".join(os.path.join(OUT, "new", f) for f in ("run_manifest.json", "input_files.tsv", "output_files.tsv")))
