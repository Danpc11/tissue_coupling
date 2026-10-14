#!/usr/bin/env bash
# Reproduce todos los numeros del articulo, en orden, desde datos publicos.
#
#   bash run_all.sh
#
# Variables: T2D_RAW (datos crudos, por defecto data/raw)
#            T2D_OUT (resultados, por defecto results)
#            T2D_BOOT, T2D_PERM (remuestreos; bajarlos acelera a costa de precision)
#
# El ultimo paso NO es opcional: verify_package comprueba que los resultados son coherentes
# entre si y proceden de una sola corrida. Un paquete que no lo pase no debe reportarse.
set -euo pipefail
export T2D_RAW="${T2D_RAW:-data/raw}"
export T2D_OUT="${T2D_OUT:-results}"
mkdir -p "$T2D_OUT"
step () { printf '\n=== %s\n' "$1"; shift; "$@"; }

step "00 download and prepare public data"   python python/analyses/00_prepare_data.py
step "01 cross-tissue concordance"           python python/analyses/01_crosstissue_axis.py
step "02 meta-analysis and validation"       python python/analyses/02_meta_and_validation.py
step "02b phase template (known times)"      python python/analyses/02b_phase_template.py
step "03 circadian phase and coupling"       python python/analyses/03_phase_coupling.py
step "04 perturbation response"              python python/analyses/04_perturbations.py
step "04 smoke test"                         python tests/smoke_test.py
step "05 verify the package"                 python tests/verify_package.py "$T2D_OUT"

printf '\nDone. Results in %s\n' "$T2D_OUT"
