# Glucocorticoid signalling tracks circadian output across human metabolic tissues

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23379059.svg)](https://doi.org/10.5281/zenodo.23379059)
[![tests](https://github.com/Danpc11/tissue_coupling/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/tissue_coupling/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Code for the manuscript. `run_all.sh` recomputes its tables and figures from public data.

Recomputing is not the same as reproducing the published numbers, and both are kept apart. Every
figure is drawn from tables computed in the run. The tables shipped with earlier releases stay in
`data/reference/` only so that each can be compared with its recomputation, and `docs/DATA.md`
lists, table by table, how close the two are, including the ones that differ.

## What the code does

It infers circadian phase and amplitude from single transcriptomic samples, validates that
inference against time series from the same donors, and uses it to ask what the output arm of
the peripheral clock follows in human skeletal muscle, adipose tissue and pancreatic islets.
Three kinds of analysis answer that question, and they do not carry the same weight:

| analysis | what it can show | what it cannot |
|---|---|---|
| **Association across individuals** between the level of the glucocorticoid receptor (*NR3C1*) and the clock-output score, pooled over cohorts and tissues | that the two vary together, how strongly, and how robust that is to covariates, tissue and cohort | that receptor availability sets clock output. Cell composition, inflammation, feeding state and the downregulation of the receptor by the hormone itself are alternatives these data do not exclude |
| **Perturbations within individuals** with a time-matched control: hydrocortisone against placebo at a fixed hour, and a clamp against saline over the same hours | which changes are produced by the hormone and which by the passage of time | that the oscillator was reset. An acute change in *PER1* is a transcriptional response; phase, amplitude and sustained synchrony are separate quantities and are estimated separately |
| **A forced-oscillator model** whose effective forcing is the product of the cortisol rhythm and receptor availability | that the proposed mechanism is sufficient to produce the observed signature, and what else it predicts | that the mechanism is the true one. Its parameters are fixed in advance rather than fitted, which is not the same as being uninformed by data |

The numbers are in the manuscript and, after a run, in `results/new/`.

## Which version is the reference

The reference implementation is the latest tagged release on `main`. Releases up to v1.3.1,
including the one first archived on Zenodo, did not contain the step that turns the public
downloads into the expression matrices, so they could not regenerate the results from a clean
download. That step was reconstructed by Hugo Tovar and merged in v2.0.0; from that release on,
`run_all.sh` runs end to end from the public data. `CHANGELOG.md` lists what each release
changed, and `docs/DATA.md` has the equivalence table between the tables shipped with the earlier
releases and their regeneration, including the ones that do not reproduce.

## Reproducing the results

```bash
pip install -r requirements.txt        # Python 3.11 or later; requirements-lock.txt for exact versions
bash download_public.sh                # public inputs into data/raw/ (see docs/DATA.md)
bash run_all.sh                        # every step, in order; figures in results/figs, tables in results/new
```

`run_all.sh` can be called from any directory. `T2D_RAW`, `T2D_OUT` and `T2D_REF` point it at
other locations for the downloads, the outputs and the versioned inputs.

**A complete run and a partial one.** The phase concordance with CHIRAL (step 06) needs an
implementation that is not on PyPI; `download_public.sh` clones it at a fixed revision. If it is
absent the step still runs with our own estimator, writes the CHIRAL rows as not run and says so
at the end of the pipeline. `T2D_REQUIRE_CHIRAL=1` makes its absence an error.

The comparison with ZeitZeiger and tauFisher (step 02f) needs R and those two packages, which are
not on CRAN (`docs/DATA.md` says how to install them). Without them the step writes the rows of
our own estimator, marks the others as not run and the pipeline says so at the end.
`T2D_RSCRIPT` selects the `Rscript` to use and `T2D_REQUIRE_R=1` makes their absence an error.

Stochastic routines use NumPy `default_rng` with seed 0. Continuous integration runs on Python 3.11.

**All thirteen estimates or none.** The meta-analysis lists the thirteen cohort and tissue
estimates of the manuscript and stops if one is missing, instead of pooling whatever is there.
`T2D_ALLOW_PARTIAL=1` lets it continue, and the tables then say that they are incomplete.

## Which script draws which figure

The file names date from earlier numberings of the figures.

| manuscript | script | output in `results/figs/` |
|---|---|---|
| Figures 1, 3 and 4 | `python/figures/paper_figures.py` | `F1`, `F3`, `F4` |
| Figure 2 | `python/figures/fig3_eje.py`, called by `paper_figures.py` | `F2` (also written as `F3_eje`) |
| Figure 5 | `python/figures/fig7_modelo.py`, called by `paper_figures.py` | `F5` (also written as `F7_modelo`) |
| Supplementary Figures 1 to 3 and 5 to 7 | `python/figures/supplementary_figures.py` | `ED1` to `ED3`, `ED5` to `ED7` |
| Supplementary Figure 4 | `python/figures/ed4_gc_architecture.py` | `ED4` |
| not in the manuscript | `python/figures/fig_clamp_vs_saline.py` | `clamp_vs_saline` |

## Layout

| path | contents |
|---|---|
| `python/analyses/` | the analysis steps |
| `python/model/` | the forced-oscillator model and the data it generates; four mechanisms scored against the same observables (`competing_models.py`); what the amplitude estimator can and cannot separate (`factorial_estimator.py`) |
| `R/` | the one R script, which runs ZeitZeiger and tauFisher for step 02f |
| `python/figures/` | the figures |
| `python/lib/gene_sets.py` | the registry of gene sets |
| `data/reference/` | versioned inputs, and how each compares with its regeneration (`docs/DATA.md`) |
| `tests/` | see below |

## Three conventions worth knowing before reading the code

**Gene sets live in one place.** Every set that defines a score, an inference or a control is in
`python/lib/gene_sets.py`, and scripts import it. The registry has two tables: `SETS`, the sets
that are scored and compared with one another, and `PANELS`, the gene panels a particular
estimator or step uses, each with its source and the reason it differs from the others. Two kinds
of list remain outside it: the ones that choose which genes a figure panel draws, and
per-gene parameter tables. `tests/test_gene_lists.py` enforces this on the code.

**Expression is centred per sample for scores, and not for group differences.** Centring removes
a common component that otherwise inflates every correlation in RNA-seq (the median correlation of
a random gene with the clock-output score is 0.44 rather than zero without it). The signature table
compares differences between groups and is built without centring. Mixing the two conventions
changes results.

**"Single sample" means one sample per person, not a sample without a cohort.** The estimators
standardise each gene across the samples they are given and learn the waveform from them, so
they are relative to a cohort. A sample that is not part of it is scored against a reference
fitted beforehand (`fit_reference` and `infer_from_reference` in `python/lib/phase.py`), which is
how the leave-one-out validation works. The amplitude is a transcriptomic proxy: it does not
separate lower amplitude in each cell from wider dispersion of phase between cells.

**"Coupling" is a statistical statement here.** The repository name refers to the relationship
between glucocorticoid signalling and clock output. Where two tissues of one person are compared,
the quantity is the coherence of their phase difference across donors. It is high when both
tissues follow a shared systemic signal, with no direct interaction between them, so it is not
evidence of communication between tissues (`tests/test_phase_statistic.py` shows this with
simulated phases).

## Tests

| test | what it checks | needs data |
|---|---|---|
| `tests/test_syntax.py` | every script compiles | no |
| `tests/verify_inputs.py` | every file a script reads has a producer that runs earlier, is versioned, or is a declared download. Static analysis: see its header for what it cannot see | no |
| `tests/test_gene_lists.py` | no gene set is defined outside the registry | no |
| `tests/test_phase_statistic.py` | what the phase-coherence statistic measures, on simulated phases | no |
| `tests/test_phase_template.py` | the validation of the phase template uses nothing from the held-out subject: on pure noise its error is that of chance | no |
| `tests/test_single_sample.py` | "single sample" means one sample per person scored against a fixed reference: a lone sample is refused, and its value does not depend on which others accompany it | no |
| `tests/test_model_null.py` | on simulated cohorts the estimator recovers the effect of the receptor gain when the model has one and finds none when it has none | no |
| `tests/test_model_verdict.py` | the criterion that compares the competing models: three states per observable; a joint region, cohort by cohort, that respects the correlation between observables; and predictions averaged over seeds, so that more seeds make it stricter | no |
| `tests/test_clamp_vs_saline.py` | runs step 14b on a small simulated series with a known answer and checks the contents of its output | no |
| `tests/smoke_test.py` | the phase estimator recovers a known phase | no |
| `tests/verify_package.py` | arithmetic and consistency of the tables of a real run | yes |

All but the last run in continuous integration and at the end of `run_all.sh`. A passing test suite
means the code does what it says; it does not validate the biological interpretation.

## Citation

See `CITATION.cff`. Cite the Zenodo record of the tagged release you used.
