#!/usr/bin/env python3
"""Inventario de lo que ha escrito la corrida, para cerrar su manifiesto.

Comparar fechas entre ficheros no dice si vienen de la misma corrida: dos corridas seguidas caben en
pocas horas y una corrida larga no. Lo que si se puede afirmar es si un fichero se escribio DESPUES de
que esta corrida empezara. El paso 01a apunta ese instante y un identificador; este paso recorre todo
el directorio de resultados, no solo una parte, y anota de cada fichero su tamano, su huella y si es
de esta corrida. Lo unico que puede ser anterior es una cache que el paso que la usa haya declarado y
validado (las simulaciones del paso 20c con T2D_MODEL_REUSE=1).

  results/new/output_files.tsv    file, bytes, sha256, modified, run_id, written_by_this_run, declared_cache
  results/new/run_manifest.json   se le anaden el numero de salidas, su huella conjunta y la hora de fin
"""
import os, sys, json, hashlib, time, datetime
OUT = os.environ.get("T2D_OUT", "results") if len(sys.argv) < 2 else sys.argv[1]
CACHES = {"new/competing_models_ensembles.pkl", "new/competing_models_simulated.tsv", "new/competing_models_cache_manifest.json"}


def sha256(path, block=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""): h.update(chunk)
    return h.hexdigest()


if __name__ == "__main__":
    fm = f"{OUT}/new/run_manifest.json"
    if not os.path.exists(fm):
        sys.exit("falta results/new/run_manifest.json: lo escribe el paso 01a al empezar la corrida")
    M = json.load(open(fm)); t0 = M["started_epoch"] - 2.0          # dos segundos de margen por la resolucion del reloj del disco
    reuse = os.environ.get("T2D_MODEL_REUSE") == "1"
    rows = []
    for d, _, files in os.walk(OUT):
        for f in sorted(files):
            p = os.path.join(d, f); rel = os.path.relpath(p, OUT)
            if f.startswith(".") or rel in ("new/output_files.tsv", "new/run_manifest.json"): continue
            mt = os.path.getmtime(p)
            rows.append(dict(file=rel, bytes=os.path.getsize(p), sha256=sha256(p),
                             modified=datetime.datetime.fromtimestamp(mt, datetime.timezone.utc).isoformat(timespec="seconds"),
                             run_id=M["run_id"], written_by_this_run=mt >= t0, declared_cache=(reuse and rel in CACHES and mt < t0)))
    rows.sort(key=lambda r: r["file"])
    with open(f"{OUT}/new/output_files.tsv", "w") as fh:
        cols = list(rows[0])
        fh.write("\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))
    old = [r["file"] for r in rows if not r["written_by_this_run"] and not r["declared_cache"]]
    M.update(finished=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), duration_minutes=round((time.time() - M["started_epoch"]) / 60, 1),
             output_files=len(rows), output_files_sha256=hashlib.sha256("".join(r["sha256"] for r in rows).encode()).hexdigest(),
             output_files_from_before_this_run=old)
    json.dump(M, open(fm, "w"), indent=1)
    print(f"corrida {M['run_id']}: {len(rows)} ficheros de salida, {M['duration_minutes']} min")
    if old:
        print(f"  {len(old)} ficheros son ANTERIORES a esta corrida y no son cache declarada: " + ", ".join(old[:8]))
