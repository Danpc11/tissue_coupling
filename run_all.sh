#!/usr/bin/env bash
# Reproduce todos los numeros del articulo, en orden, desde datos publicos.
#
#   bash run_all.sh
#
# Variables: T2D_RAW (datos crudos, por defecto data/raw)
#            T2D_OUT (resultados, por defecto results)
#            T2D_BOOT, T2D_PERM (remuestreos; bajarlos acelera a costa de precision)
#
# El paso de verificacion NO es opcional: verify_package comprueba que los resultados son
# coherentes entre si y proceden de una sola corrida. Un paquete que no lo pase no se reporta.
set -euo pipefail
export T2D_RAW="${T2D_RAW:-data/raw}"
export T2D_OUT="${T2D_OUT:-results}"
mkdir -p "$T2D_OUT" "$T2D_OUT/new" "$T2D_OUT/figs"
step () { printf '\n=== %s\n' "$1"; shift; "$@"; }

# --- datos y analisis base
step "00 download and prepare public data"   python python/analyses/00_prepare_data.py
step "01 cross-tissue concordance"           python python/analyses/01_crosstissue_axis.py
step "01b conventional differential expression (Fig. 1)" python python/analyses/01b_fig1_signature_test.py
step "02 meta-analysis and validation"       python python/analyses/02_meta_and_validation.py
step "02b phase template (known times)"      python python/analyses/02b_phase_template.py
step "03 circadian phase and coupling"       python python/analyses/03_phase_coupling.py
step "04 perturbation response"              python python/analyses/04_perturbations.py
step "05 within-person clamp and TZD"        python python/analyses/05_clamp_within_person.py

# --- eje receptor-reloj (v1.1)
step "06 random-effects meta-analysis"       python python/analyses/06_meta_analysis.py
step "07 TF activity (DoRothEA A-C)"         python python/analyses/07_tf_activity.py
step "08 GTEx adjusted for ischaemia and RIN" python python/analyses/08_gtex_adjusted.py

# --- modelo de oscilador forzado (lento: ~5 min, cachea)
step "09 entrainment tongue sweep"           python python/model/make_tongue.py
step "10 synthetic single-biopsy cohorts"    python python/model/make_fig7_data.py

# --- figuras
step "11 assemble the signature table"      python python/analyses/23_build_signatures_table.py
step "12 HPO enrichment by tissue"          python python/analyses/09_hpo_by_tissue.py
step "13 paper figures F1-F5"                python python/figures/paper_figures.py

# --- verificacion
step "14 inputs have a producer"             python tests/verify_inputs.py
step "15 smoke test"                         python tests/smoke_test.py
step "16 verify the package"                 python tests/verify_package.py "$T2D_OUT"

printf '\nDone. Results in %s, figures in %s/figs\n' "$T2D_OUT" "$T2D_OUT"
