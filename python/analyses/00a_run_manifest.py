#!/usr/bin/env python3
"""Deja escrito que corrida es esta: con que codigo, con que ajustes y con que ficheros de entrada.

Dos corridas del mismo codigo pueden no ser comparables. Una corrida rapida, con menos semillas y
menos replicas, sirve para comprobar que todo funciona y no para citar sus cifras; y unas cuentas
de islotes de otra procedencia dan otros numeros sin que ningun paso falle. Este paso lo registra
antes de calcular nada, y verify_package.py lo lee al final.

  results/new/run_manifest.json   perfil de la corrida (full o quick), ajustes efectivos y cuales
                                  difieren de los de una reproduccion completa; version y commit;
                                  versiones de Python y de las librerias; origen de las cuentas de
                                  islotes
  results/new/input_files.tsv     cada fichero de entrada con su tamano y su huella sha256
"""
import os, sys, json, hashlib, subprocess, platform, datetime
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "../..")); sys.path.insert(0, os.path.join(HERE, "../lib"))
import islets
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results")

# Los ajustes que cambian las cifras, con el valor de una reproduccion completa. El perfil quick de
# run_all.sh los reduce; cualquier otro valor distinto de estos tambien queda anotado.
FULL = {"T2D_MODEL_N": "200", "T2D_MODEL_SEEDS": "10", "T2D_MODEL_BOOT": "500", "T2D_MODEL_CALIB": "100",
        "T2D_MODEL_CALIB_BOOT": "150", "T2D_MODEL_CALIB_POOLS": "2", "T2D_MODEL_CALIB_FRESH": "30",
        "T2D_MODEL_CALIB_COMPARE": "1", "T2D_MODEL_CONFUSION": "50", "T2D_NULL": "300", "T2D_BOOT": "400",
        "T2D_PERM": None, "T2D_BOOT_META": "500", "T2D_CENTER": "1", "T2D_ALLOW_PARTIAL": None,
        "T2D_GSE13070_SYMBOLS": None, "T2D_SIGNATURES": "regenerated", "T2D_LOO_COLUMN": "amp_ls_donor_out",
        "T2D_CLOCK_COHORT": "GSE108539", "T2D_CLAMP_TEMPLATE": "800"}


def sha256(path, block=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""): h.update(chunk)
    return h.hexdigest()


def git(*args):
    try:
        return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, timeout=20).stdout.strip() or None
    except Exception:
        return None


def versions():
    out = {"python": platform.python_version()}
    for m in ("numpy", "pandas", "scipy", "sklearn", "statsmodels", "pydeseq2", "matplotlib"):
        try: out[m] = __import__(m).__version__
        except Exception: out[m] = None
    return out


if __name__ == "__main__":
    os.makedirs(f"{OUT}/new", exist_ok=True)
    profile = os.environ.get("T2D_PROFILE", "full")
    effective = {k: os.environ.get(k, v) for k, v in FULL.items()}
    changed = {k: dict(this_run=effective[k], full_reproduction=v) for k, v in FULL.items() if effective[k] != v}
    release = next((l.split(":", 1)[1].strip() for l in open(os.path.join(ROOT, "CITATION.cff")) if l.startswith("version:")), None)
    rows = []
    for d, _, files in os.walk(RAW, followlinks=True):
        if any(p in d for p in ("pychiral", "archs4", ".git")): continue
        for f in sorted(files):
            p = os.path.join(d, f)
            if f.startswith(".") or not os.path.isfile(p): continue
            rows.append(dict(file=os.path.relpath(p, RAW), bytes=os.path.getsize(p), sha256=sha256(p)))
    rows.sort(key=lambda r: r["file"])
    with open(f"{OUT}/new/input_files.tsv", "w") as fh:
        fh.write("file\tbytes\tsha256\n" + "".join(f"{r['file']}\t{r['bytes']}\t{r['sha256']}\n" for r in rows))
    isl = islets.counts_source(RAW)
    import uuid, time
    commit = git("rev-parse", "HEAD")
    M = dict(run_id=f"{datetime.datetime.now(datetime.timezone.utc):%Y%m%dT%H%M%SZ}-{(commit or 'nogit')[:8]}-{uuid.uuid4().hex[:6]}",
             started_epoch=time.time(), written=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
             profile=profile, scientific_reproduction=(profile == "full" and not changed),
             settings_that_differ_from_a_full_reproduction=changed, settings=effective,
             release=release, git_commit=git("rev-parse", "HEAD"),
             git_uncommitted_changes=bool(git("status", "--porcelain", "--untracked-files=no")),     # cambios en ficheros versionados
             git_untracked_files=len((git("status", "--porcelain", "--untracked-files=all") or "").splitlines())
                                 - len((git("status", "--porcelain", "--untracked-files=no") or "").splitlines()),
             versions=versions(), input_files=len(rows),
             input_files_sha256=hashlib.sha256("".join(r["sha256"] for r in rows).encode()).hexdigest(),
             islet_counts=isl)
    json.dump(M, open(f"{OUT}/new/run_manifest.json", "w"), indent=1)
    print(f"corrida {M['run_id']}")
    print(f"perfil de la corrida: {profile}" + ("" if M["scientific_reproduction"] else "  -> NO es una reproduccion completa: sus cifras no son las del articulo"))
    for k, v in changed.items(): print(f"  {k} = {v['this_run']} (reproduccion completa: {v['full_reproduction']})")
    print(f"version {release}, commit {M['git_commit'] or 'desconocido'}" + (" con CAMBIOS SIN CONFIRMAR en ficheros versionados" if M["git_uncommitted_changes"] else "")
          + (f" ({M['git_untracked_files']} ficheros sin versionar, que no cuentan como cambio)" if M["git_untracked_files"] else ""))
    print(f"{len(rows)} ficheros de entrada con su huella en {OUT}/new/input_files.tsv")
    print(f"cuentas de islotes: {isl['source']}" + ("" if isl["is_the_reported_source"] else "  -> NO es la fuente de las cifras del articulo"))
