# Changelog

## v2.2.2: the validations use nothing from what they are tested on

A second external review found information leaking from the test subject into the validation of
the phase template. It was right. The same question was then put to the two validations whose
numbers are in the manuscript.

### Fixed
- **Phase template (step 02): the gene selection and the standardisation saw the held-out
  subject.** Genes were ranked by R² on every sample, each gene was standardised on every sample,
  and only then was a subject left out; and the error announced was the smallest of twelve
  configurations. On a simulated cohort of pure noise that procedure reported a median phase error
  of 1.1 h, where chance is 6 h. Each fold now learns its standardisation parameters, its gene
  ranking and its template from the training subjects only, and the configuration is chosen by a
  second leave-subject-out inside the training set (nested validation). On the authors' normalised
  matrix of GSE108539 the announced error moves from 0.79 h to **0.97 h** (median; 84% of samples
  within 3 h, against 95%). The definitive figure is the one of a full run on the pipeline's own
  matrix. No number in the manuscript rests on the old value: its text gives 1.0 to 1.2 h for
  supervised estimators.
- The subject and the sampling hour of that cohort are read from free-text fields. The step now
  stops with a message if an hour cannot be read or the subjects do not come out of the titles.
- `run_all.sh` could report the CHIRAL concordance as regenerated when it was not.

### Added
- **Amplitude validation (steps 04 and 05): a second column with the whole donor excluded.** The
  manuscript says the whole donor was excluded from the ordering and the template. The code
  excluded the series of the held-out point, so the same donor's samples in the other condition
  stayed in the base. `amp_inf` is unchanged, to the last bit, and `amp_inf_donor_out` builds the
  base without any sample of that donor. The two agree: in myotubes ρ = 0.28 and 0.34 per sample
  become 0.32 and 0.32, and 0.97 and 0.94 per donor become 0.98 and 0.87; in muscle in vivo 0.17
  and 0.82 become 0.17 and 0.81. The gene panel of these validations is fixed in advance and the
  standardisation already used the base only, so there was nothing else to correct.
- Step 06 adds a row with our estimator restricted to the genes CHIRAL receives. Each estimator in
  its usual configuration uses different genes (60 against 22), so a difference between them could
  come from the genes; between this row and CHIRAL's only the algorithm changes. `n_genes` now
  says how many genes each row used.
- `tests/test_phase_template.py`: what a fold learns does not change when the held-out subject's
  samples are replaced by noise; on a cohort of pure noise the error is that of chance, and the
  previous procedure, reimplemented in the test, is shown to do much better than chance there;
  on a rhythmic cohort the nested validation recovers the phase.

## v2.2.1: answers to an external review of the repository

No reported number changes. The lists that moved into the registry are identical in content and
order to the ones they replace, and the model step, rerun with the registry, returns the published
signature to the third decimal.

### Fixed
- **CHIRAL could not run on a clean clone, for two reasons that hid each other.** v2.1.0 made
  `download_public.sh` clone scRITMO, but the backend looked for it in `/tmp/scRITMO/src` and
  `run_all.sh` pointed `PYCHIRAL_SRC` at a directory in a developer's home. And had it been found,
  the call would have failed: the backend passed the genes as `clockgenes=`, which is the name of
  the argument in pyCHIRAL and not in scRITMO. The genes are now passed by position, which works
  with both; there is one documented install path, pyCHIRAL cloned by `download_public.sh` into
  `data/raw/pychiral` at commit `e1ea76e`, and the backend looks there. Checked on a simulated
  matrix with known phases: both implementations return one phase per sample and recover the
  truth at R = 0.995, so the genes-by-samples orientation of the wrapper is right.
- **Ten analysis steps ignored `T2D_RAW` and `T2D_OUT`** and read `data/raw/...` and
  `results/new/...` relative to the working directory, so the pipeline only worked when run from
  the repository root with the default locations. They now follow the same variables as the rest;
  with the defaults all 35 paths are the same strings as before. `run_all.sh` changes to its own
  directory first. Step 07 was rerun on the real GSE107537 from another directory with both
  variables set elsewhere and returns the same values (median −28.6%, 8 of 8, *P* = 0.0078).
- **The registry was not the only place where gene sets were defined**, although the README said
  so. Seven lists that define a score or an inference lived in the scripts that used them, one of
  them duplicated in two files. They are now in `python/lib/gene_sets.py` under `PANELS`, apart
  from the scored sets of `SETS` so that no existing analysis gains a row, each with its source,
  and the scripts import them.
- **`verify_inputs.py` accepted any matrix whose name began with GSE, GTEx or AAGMEx.** It now
  reads from the syntax tree which matrices `save()` actually writes (28, all resolved by name)
  and includes them in the check of execution order. Its header states what a static analysis
  cannot see.
- The README attributed the data-preparation step to a branch that has been part of `main` since
  v2.0.0, stated Python 3.9 where 3.11 is required, and summarised results as if association were
  mechanism. `docs/DATA.md` described two different ways of installing CHIRAL and said that a
  contrast "was not reproduced consistently" which the pipeline now regenerates exactly.

### Added
- `tests/test_gene_lists.py`: fails if a collection of five or more registry genes is defined
  outside the registry and is not a figure panel or one of four declared parameter tables.
- `tests/test_clamp_vs_saline.py`: an integration test. It runs step 14b on a small simulated
  series with a known answer, including a pair deliberately made of two people, and checks the
  contents of the output: the pair is excluded, the time effect and the insulin effect are
  recovered, the tables have the expected rows and columns.
- `tests/test_phase_statistic.py`: what the phase-coherence statistic measures, on simulated
  phases. Independent tissues give no signal; aligned tissues give a high value; **so do two
  tissues that follow a shared signal with no interaction between them**; and a tissue whose
  cycle is traversed in the opposite direction gives no signal although it is aligned.
- That last case is not hypothetical. **CHIRAL has no anchor, and the direction in which it
  traverses the cycle changes with the order of the samples**: on the same simulated cohort,
  sorting the donors alphabetically turned R from 0.84 into 0.09. Step 06 therefore writes a
  second column, `R_reflection_aware`, that takes the larger of the two directions and uses the
  same maximum in the permutation null (5.7% false positives at the 5% level in 300 simulated
  null cohorts). The original column is unchanged. Which one is reported for CHIRAL is a decision
  for the manuscript; for our own estimator, anchored to the template, the two coincide.
- Step 06 writes whether CHIRAL ran and with which implementation and revision, `run_all.sh`
  ends by saying whether the concordance was regenerated, and `T2D_REQUIRE_CHIRAL=1` makes its
  absence an error. CHIRAL's warnings are counted and reported instead of silenced.
- `docs/DATA.md`: the equivalence table between each versioned table and its regeneration in a
  complete run, including the three that do not reproduce.
- The four no-data tests run in continuous integration, which until now only checked syntax and
  the smoke test.

## v2.2.0: a time-matched control for the clamp

Adds one analysis and the data it needs. It does not change any existing result, and how it
enters the manuscript has not been decided.

### Added
- `python/analyses/28_clamp_vs_saline.py`, step 14b. GSE40589 has adipose biopsies before and
  after four hours of saline alone, of a euglycaemic clamp and of a hyperglycaemic clamp, all
  starting at 08:00, so the change under saline is the passage of time and the difference between
  a clamp and saline is what insulin adds. The step reports, for the clock-output, core-oscillator,
  glucocorticoid and insulin-response sets, the change of every pair, the test within each arm and
  the contrast of each clamp with saline, and the same gene by gene for the clock.
- **An identity rule for the pairing.** The series does not declare the participant, so the pairs
  are taken from the numbering of the sample titles and checked with genes that identify an
  individual. A pair whose identity correlation falls below the 90th percentile of the
  correlations between different pairs is excluded from the primary analysis. Three analyses are
  written: primary, all pairs, and one that does not use the pairing.
- `python/figures/fig_clamp_vs_saline.py`, step 27b: one point per participant, excluded pairs
  drawn hollow, P values written out.
- `download_public.sh` fetches the GSE40589 series matrix and the GPL6848 platform table;
  `00_prepare_data.py` lists the series, so its absence stops the run at step 01; `docs/DATA.md`
  describes it.

### Fixed
- `CITATION.cff` carried the title of an earlier version of the manuscript and version 1.1.0.
- Panel letters asked for one font by name and printed a warning for every panel on machines
  without it; they now follow the font list of `journal.mplstyle`.


## v2.1.0: a reviewer can clone this and run it

Nothing here changes a result. These are the things that stopped a clean clone from running.

### Fixed
- **`requirements.txt` did not list what the pipeline imports.** `pydeseq2` (islets),
  `statsmodels` (heterogeneity moderators) and `anndata` were missing, so a reviewer installing
  from it hit an ImportError three steps in. The file now lists every runtime dependency with
  minimum versions, says that Python 3.11 or later is required (the lock pins numpy 2.4 and
  pandas 3.0, which do not exist for 3.9), and names the optional ones.
- **`download_public.sh` did not fetch two files the pipeline reads.** `sra.recount_qc.*` for
  SRP127581 and SRP055749 carry the average read length that converts recount3 base coverage into
  reads; without them the Cushing and phase-template steps failed. No test caught this, because
  the tests check that every input has a producer, not that the download script fetches everything
  the producers read.
- **NCBI's CAPTCHA aborted the islet step.** The counts for GSE159984 are not optional. The script
  now retries five times with increasing waits, deletes any HTML it receives, and prints the ARCHS4
  fallback if it still fails.
- **CHIRAL silently reported the wrong cause.** Step 06 said `python/lib/chiral_backend.py` was
  missing when the file was there and scRITMO, which is not on PyPI, was not. The script is now
  cloned by `download_public.sh` and the message names the real cause.
- **`verify_package` failed two checks by design**, looking for exports of a script that is not a
  pipeline step. It now reads `meta_per_cohort.tsv`, `meta_pooled.tsv` and
  `new/crosstissue_prediction.tsv`, which the pipeline does produce, and says which step writes
  them when they are absent.


## v2.0.0: the version behind the Nature Communications submission

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


## v1.3.0: the repository is the code behind the paper

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

## v1.2.x: reproducibility chain closed

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

## v1.1.0: glucocorticoid axis, forced-oscillator model, five figures
