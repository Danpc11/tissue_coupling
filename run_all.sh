#!/usr/bin/env bash
# Reproduce every result and figure of
#   "Glucocorticoid signalling tracks circadian output across human metabolic tissues"
# from public data. Steps run in order; each writes what the next one reads.
#
#   bash run_all.sh
#
# Expects the public downloads listed in docs/DATA.md under data/raw/, and the GTEx v11
# annotation files for the GTEx steps. Everything else is produced here or versioned in
# data/reference/ (see docs/DATA.md for why).
set -euo pipefail
export T2D_RAW="${T2D_RAW:-data/raw}"
export T2D_OUT="${T2D_OUT:-results}"
export T2D_REF="${T2D_REF:-data/reference}"
mkdir -p "$T2D_OUT/new" "$T2D_OUT/figs"
step () { printf '\n=== %s\n' "$1"; shift; "$@"; }

step "00 every script compiles"                "${PY:-python3}" tests/test_syntax.py

# --- data ------------------------------------------------------------------------------
step "01 prepare the public series"            "${PY:-python3}" python/analyses/00_prepare_data.py

# --- the estimator and its validation (Fig. 1) -------------------------------------------
step "02 phase and amplitude template"         "${PY:-python3}" python/analyses/02b_phase_template.py
step "03 case-control analysis, Fig. 1b-d"     "${PY:-python3}" python/analyses/01b_fig1_signature_test.py
step "04 leave-one-out, myotubes"              "${PY:-python3}" python/analyses/11_loo_validation_myotubes.py
step "05 leave-one-out, muscle in vivo"        "${PY:-python3}" python/analyses/12_loo_validation_invivo.py
step "06 phase concordance with CHIRAL"        "${PY:-python3}" python/analyses/17_chiral_concordance.py
step "07 external control: 10 h sleep delay"   "${PY:-python3}" python/analyses/26_shift_work_control.py

# --- the receptor-clock axis (Fig. 2) -----------------------------------------------------
step "08 cross-tissue axis"                    "${PY:-python3}" python/analyses/01_crosstissue_axis.py
step "09 GTEx adjusted for ischaemia and RIN"  "${PY:-python3}" python/analyses/08_gtex_adjusted.py
step "10 meta-analysis, two-level"             "${PY:-python3}" python/analyses/06_meta_analysis.py
step "11 meta-analysis, three-level"           "${PY:-python3}" python/analyses/13_meta_three_level.py
step "12 heterogeneity moderators"             "${PY:-python3}" python/analyses/25_meta_moderators.py
step "13 sensitivity to principal components"  "${PY:-python3}" python/analyses/22_pc_sensitivity_core.py

# --- perturbations (Fig. 3) ----------------------------------------------------------------
step "14 within-person clamp and TZD"          "${PY:-python3}" python/analyses/05_clamp_within_person.py
step "15 hydrocortisone, dexamethasone, Cushing" "${PY:-python3}" python/analyses/04_perturbations.py
step "16 assemble the signature table"         "${PY:-python3}" python/analyses/23_build_signatures_table.py

# --- transcriptional context (Fig. 4) --------------------------------------------------------
step "17 transcription-factor activity"        "${PY:-python3}" python/analyses/07_tf_activity.py
step "18 phenotype ontology by tissue"         "${PY:-python3}" python/analyses/09_hpo_by_tissue.py

# --- the model and its prediction (Fig. 5) -----------------------------------------------------
step "19 entrainment tongue"                   "${PY:-python3}" python/model/make_tongue.py
step "20 synthetic single-biopsy cohorts"      "${PY:-python3}" python/model/make_fig7_data.py
step "21 amplitude dispersion"                 "${PY:-python3}" python/analyses/10_amplitude_dispersion.py
step "22 estimator noise by group"             "${PY:-python3}" python/analyses/15_dispersion_noise_by_group.py

# --- figures ------------------------------------------------------------------------------------
step "23 Figures 1, 3 and 4"                   "${PY:-python3}" python/figures/paper_figures.py
step "24 Figure 2"                             "${PY:-python3}" python/figures/fig3_eje.py
step "25 Figure 5"                             "${PY:-python3}" python/figures/fig7_modelo.py
step "26 Supplementary Figures S1-S3, S5-S7"   "${PY:-python3}" python/figures/supplementary_figures.py
step "27 Supplementary Figure S4"              "${PY:-python3}" python/figures/ed4_gc_architecture.py

# --- checks ---------------------------------------------------------------------------------------
step "28 every figure input has a producer"    "${PY:-python3}" tests/verify_inputs.py
step "29 smoke test"                           "${PY:-python3}" tests/smoke_test.py
step "30 verify the package"                   "${PY:-python3}" tests/verify_package.py "$T2D_OUT"

printf '\nDone. Figures in %s/figs, tables in %s/new.\n' "$T2D_OUT" "$T2D_OUT"
