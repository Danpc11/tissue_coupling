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
cd "$(dirname "${BASH_SOURCE[0]}")"        # las rutas son relativas a la raiz del repositorio
export T2D_RAW="${T2D_RAW:-data/raw}"
export T2D_OUT="${T2D_OUT:-results}"
export T2D_REF="${T2D_REF:-data/reference}"
# CHIRAL (paso 06) se busca en $T2D_RAW/pychiral, donde lo clona download_public.sh. Si no esta,
# el paso lo dice y sigue con el estimador propio; T2D_REQUIRE_CHIRAL=1 lo vuelve obligatorio.
mkdir -p "$T2D_OUT/new" "$T2D_OUT/figs"
step () { printf '\n=== %s\n' "$1"; shift; "$@"; }

step "00 every script compiles"                "${PY:-python3}" tests/test_syntax.py

# --- data ------------------------------------------------------------------------------
step "01 prepare the public series"            "${PY:-python3}" python/analyses/00_prepare_data.py
step "01b build the expression matrices"       "${PY:-python3}" python/analyses/00b_build_matrices.py

# --- the estimator and its validation (Fig. 1) -------------------------------------------
step "02 phase and amplitude template"         "${PY:-python3}" python/analyses/02b_phase_template.py
step "02c islets, T2D vs non-diabetic (PyDESeq2)" "${PY:-python3}" python/analyses/01c_islets_deseq2.py
step "02d regenerate the versioned reference tables" "${PY:-python3}" python/analyses/00d_reference_tables.py
step "02e template benchmark in adipose, GSE87761" "${PY:-python3}" python/analyses/02e_gse87761_benchmark.py
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
step "11b cross-tissue prediction"            "${PY:-python3}" python/analyses/27_crosstissue_prediction.py
step "13 sensitivity to principal components"  "${PY:-python3}" python/analyses/22_pc_sensitivity_core.py

# --- perturbations (Fig. 3) ----------------------------------------------------------------
step "14 within-person clamp and TZD"          "${PY:-python3}" python/analyses/05_clamp_within_person.py
step "14b clamp vs saline at the same clock times" "${PY:-python3}" python/analyses/28_clamp_vs_saline.py
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
step "27b clamp vs saline figure"              "${PY:-python3}" python/figures/fig_clamp_vs_saline.py

# --- checks ---------------------------------------------------------------------------------------
step "28 every figure input has a producer"    "${PY:-python3}" tests/verify_inputs.py
step "28b no gene set outside the registry"    "${PY:-python3}" tests/test_gene_lists.py
step "29 smoke test"                           "${PY:-python3}" tests/smoke_test.py
step "29b phase-coherence statistic"           "${PY:-python3}" tests/test_phase_statistic.py
step "29c integration test, step 14b"          "${PY:-python3}" tests/test_clamp_vs_saline.py
step "30 verify the package"                   "${PY:-python3}" tests/verify_package.py "$T2D_OUT"

# --- what this run did not do ---------------------------------------------------------------------
printf '\n=== run complete\n'
if grep -q "CHIRAL.*[[:space:]]ok[[:space:]]" "$T2D_OUT/new/chiral_concordance.tsv" 2>/dev/null; then
  echo "CHIRAL concordance (step 06): regenerated."
else
  echo "CHIRAL concordance (step 06): NOT regenerated in this run. The implementation was not found;"
  echo "  bash download_public.sh clones it. Everything else is complete."
fi
