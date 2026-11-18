#!/usr/bin/env python3
"""Cierra el manifiesto de la corrida, despues de que la verificacion haya pasado, y dice que se reprodujo.

Que una corrida use la configuracion de referencia no dice que haya terminado bien ni que lo haya
calculado todo. Son tres cosas distintas y el manifiesto las guarda por separado:

  reference_configuration   perfil full con los ajustes por defecto (lo fija el paso 01a).
  verification_passed       verify_package.py termino sin fallos. Solo este paso lo pone a verdadero, y
                            run_all.sh solo llega aqui si el paso anterior, la verificacion, no fallo.
  components                que partes se calcularon en esta corrida:
      main_results              las tablas y figuras del articulo: es lo que comprueba la verificacion.
      chiral                    la concordancia con CHIRAL (paso 06), que depende de un paquete externo.
      supervised_benchmark      la comparacion con ZeitZeiger y tauFisher (paso 02f), que depende de R y
                                de dos paquetes que no estan en CRAN.

De ahi salen dos afirmaciones, que no son la misma:
  reproduces_main_results          configuracion de referencia, verificacion superada y fuente de islotes correcta.
  reproduces_everything            lo anterior y ademas los dos componentes externos calculados.
"""
import os, sys, json, csv, datetime
OUT = os.environ.get("T2D_OUT", "results") if len(sys.argv) < 2 else sys.argv[1]


def ran(path, column, value, status_column, ok="ok"):
    """Hay alguna fila de la tabla con ese metodo y estado ok?"""
    if not os.path.exists(path): return False
    with open(path, newline="") as fh:
        return any(r.get(column) == value and r.get(status_column) == ok for r in csv.DictReader(fh, delimiter="\t"))


if __name__ == "__main__":
    fm = f"{OUT}/new/run_manifest.json"
    M = json.load(open(fm))
    chiral = ran(f"{OUT}/new/chiral_concordance.tsv", "method", "CHIRAL", "status")
    bench = ran(f"{OUT}/new/benchmark_methods.tsv", "method", "ZeitZeiger", "status")
    isl = M["islet_counts"]
    inputs_ok = bool(isl.get("is_the_reported_source")) and isl.get("annotation_is_the_versioned_table") is not False
    M.update(verification_passed=True, closed=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
             components=dict(main_results="reproduced", chiral="computed" if chiral else "not computed: the package was not found",
                             supervised_benchmark="computed" if bench else "not computed: R with ZeitZeiger and tauFisher was not available"),
             reproduces_main_results=bool(M["reference_configuration"] and inputs_ok),
             reproduces_everything=bool(M["reference_configuration"] and inputs_ok and chiral and bench))
    json.dump(M, open(fm, "w"), indent=1)
    print(f"Run {M['run_id']}, release {M['release']}")
    print(f"  configuration: {'reference (profile full, default settings)' if M['reference_configuration'] else 'NOT the reference: profile ' + M['profile']}")
    for k, v in M["settings_that_differ_from_a_full_reproduction"].items():
        print(f"      {k} = {v['this_run']} (reference: {v['full_reproduction']})")
    print("  verification:  passed")
    if not inputs_ok:
        print(f"  islet inputs:  NOT those of the reported numbers ({isl['source']}"
              + ("; another table of gene symbols" if isl.get("annotation_is_the_versioned_table") is False else "") + ")")
    print(f"  CHIRAL (step 06):                     {M['components']['chiral']}")
    print(f"  ZeitZeiger and tauFisher (step 02f):  {M['components']['supervised_benchmark']}")
    if M["reproduces_everything"]:
        print("This run reproduces every result of the release, the external benchmarks included.")
    elif M["reproduces_main_results"]:
        print("This run reproduces the main results of the release. It is NOT a reproduction of everything:\n"
              "the external benchmarks listed above as not computed are missing.")
    else:
        print("This run is NOT a reproduction of the reported numbers: see the lines above.")
