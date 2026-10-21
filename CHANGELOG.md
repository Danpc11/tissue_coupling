# Changelog

## v1.2.0 — the repository reproduces the paper

The repository is now the code that produces the published results, and nothing else.

### Removed
- Eleven exploratory scripts that do not feed the manuscript: the phase-error benchmark against
  published estimators, the CHIRAL reproducibility check, the nine-tissue GTEx panel and its two
  audits, a superseded noise check, a receptor-versus-targets comparison, and four figure scripts
  from earlier versions of the paper. Those analyses informed decisions but their outputs are not
  in the manuscript, and keeping them invited questions the paper does not answer.
- The working log and the review checklists, which belong with the authors and not with the code.

### Fixed
- **The README described a different paper.** Title, summary and the perturbation table came from
  an earlier analysis. Rewritten to match the manuscript.
- **The Zenodo badge showed "pending"** instead of doi:10.5281/zenodo.23379059.
- **`run_all.sh` failed on a clean clone**: the step that assembles the signature table, which
  three panels read, was missing. Added as step 15, grouped with the perturbations it belongs to.
- **Absolute paths** from a working environment in one figure script.
- **Gene sets were redefined inside two figure panels**, each differently, adding BMAL1, which
  belongs to the core oscillator and not to the output arm. Both now import from
  `python/lib/gene_sets.py`, as that module instructs.
- **A `\uXXXX` escape inside an f-string expression**, which Python before 3.12 rejects and CI
  runs on 3.11.

### Added
- `tests/verify_inputs.py`: every input a figure reads must be produced by the pipeline, versioned
  in `data/reference/`, or a declared public download. This is the check that would have caught the
  missing step.
- `tests/test_syntax.py`, as step 00, so a syntax error stops the run in the first second.
- `docs/DATA.md`: what to download, and why each versioned input is versioned.
- `requirements-lock.txt` with exact library versions.
- Three-level meta-analysis, heterogeneity moderators, leave-one-out validation, principal-component
  sensitivity, estimator noise by group, and phase concordance with CHIRAL.

### Changed
- `run_all.sh` is grouped by the figure each block feeds, so a reader can find the code behind any
  panel without reading the whole pipeline.
