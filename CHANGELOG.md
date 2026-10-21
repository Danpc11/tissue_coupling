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

## v1.2.3 — the external control reaches the paper

### Added
- `python/analyses/26_shift_work_control.py`: the response of inferred amplitude to a controlled
  10-h delay of the habitual sleep period (GSE107537), as step 07 of `run_all.sh`. Amplitude falls
  in all eight participants, by a median of 29% (P = 0.008), while phase moves by 0.6 h — which is
  the published result for that cohort. This validates the amplitude arm of the estimator against an
  external finding, and is now panel K of Figure 1.
- Kervezee et al. 2018 to the reference list, and GSE107537 to `docs/DATA.md` with the platform
  annotation it needs.

### Fixed
- `00_prepare_data.py` listed three accessions no other script uses, and omitted GSE107537; the list
  now matches the nineteen series of the manuscript. A missing file now exits with an error instead
  of letting the run fail three steps later.
- The transition schematic in Figure 1 was still drawn although it had been removed from the
  caption; its slot now holds the external control.
