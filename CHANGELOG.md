# Changelog

## v2.0.0 — the version behind the Nature Communications submission

### Added
- `27_crosstissue_prediction.py`, as step 11b. The abstract claimed that the signal of one tissue
  predicts the response in another, and nothing in the pipeline produced it: the function lived in
  `02_meta_and_validation.py`, which is not a step. It now runs and writes
  `crosstissue_prediction.tsv` and the out-of-fold predictions per donor. Across the 647 GTEx
  donors with both tissues, the adipose glucocorticoid signal predicts the muscle target score with
  **R² = 0.44** out of fold (0.52 unadjusted), against **0.03** using age, sex and Hardy scale
  alone. Predictor and outcome genes are disjoint, both tissues of a donor share a fold, and all
  fitting happens inside the training fold.

### Fixed
- **Panels 1F and 1G did not show the leave-one-out validation the text reports.** They read the
  tables without leave-one-out and drew ρ = 0.90, 0.88 and 0.69 while the caption said 0.28, 0.34,
  0.17 and 0.82. They now read `GSE182117_loo_validation.tsv` and `GSE129843_loo_validation.tsv`,
  and the hardcoded "22 series" is gone: the count comes from the table and is 20.

### Verified
- GTEx skeletal muscle gives **ρ = 0.279 with n = 818** in the meta-analysis, reproduced
  independently from the gene-reads GCT. All 818 donors carry phenotype data, and the confidence
  interval half-width is 0.06, as that sample size requires. The ρ ≈ 0.04 with a wide interval seen
  in an earlier submitted figure does not come from this code.


## v1.3.0 — the repository is the code behind the paper

The repository now contains the model and the code that reproduces the published results, and
nothing else. No manuscript drafts, no figure files, no working notes.

### Fixed
- **The "within 4%" agreement of the model was not reproducible.** It came from reading the forcing
  amplitude of the model (`F = 0.04`) as a percentage; the ratio of mean magnitudes was computed in
  the figure script but never printed or used. Recomputed on the output-arm genes: all seven agree
  in direction and the mean magnitude of the model is **within 14%** of the observed one. The panel
  now prints that figure in its own title and the script writes the per-gene comparison to
  `model_vs_observed.tsv`, so the claim can be checked from the repository.
- **BMAL1 was included in the model comparison.** It belongs to the core oscillator, not the output
  arm, and the model overshoots it by a factor of ten; removing it takes the error from 31% to 14%.

### Removed
- `03_phase_coupling.py`: nothing in the paper reads what it produces.
- Manuscript drafts, supplementary documents, figure files and working notes. The repository is
  code and data, and the paper is the paper.

### Added
- `26_shift_work_control.py`: the response of inferred amplitude to a controlled 10-h delay of the
  sleep period (GSE107537), as a step of `run_all.sh`.
- `requirements-lock.txt` with exact library versions.

## v1.2.x — reproducibility chain closed

### Fixed
- The step that assembles the signature table, which three panels read, was missing: `run_all.sh`
  failed on a clean clone. Added as `23_build_signatures_table.py`, which regenerates three of its
  four columns from public data and **checks each against the versioned reference**, aborting on
  disagreement.
- Gene sets were redefined inside two figure panels, each differently. Both now import from
  `python/lib/gene_sets.py`, as that module instructs.
- Absolute paths from a working environment in one figure script.
- `00_prepare_data.py` listed three accessions no other script uses and omitted GSE107537; a
  missing file now exits with an error instead of letting the run fail three steps later.
- A `\uXXXX` escape inside an f-string expression, which Python before 3.12 rejects and CI runs on
  3.11.
- The README described a different paper: title, summary and the perturbation table came from an
  earlier analysis.

### Added
- `tests/verify_inputs.py`: every input a figure reads must be produced by the pipeline, versioned
  in `data/reference/`, or a declared public download.
- `tests/test_syntax.py`, as step 00, so a syntax error stops the run in the first second.
- `docs/DATA.md`: what to download, and why each versioned input is versioned.
- Three-level meta-analysis, heterogeneity moderators, leave-one-out validation,
  principal-component sensitivity, estimator noise by group, and phase concordance with CHIRAL.

## v1.1.0 — glucocorticoid axis, forced-oscillator model, five figures
