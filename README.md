# Glucocorticoid signalling tracks circadian output across human metabolic tissues

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23379059.svg)](https://doi.org/10.5281/zenodo.23379059)
[![tests](https://github.com/Danpc11/tissue_coupling/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/tissue_coupling/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Code and data for the manuscript. Everything in the paper is produced by `run_all.sh` from
public data.

## What the paper shows

The circadian clock of peripheral tissues is disrupted in obesity and type 2 diabetes, but the
signal that sets it in humans is unknown. Across 13 cohort–tissue estimates from four tissue
types, the level of the glucocorticoid receptor (NR3C1) predicts the output arm of the clock
(pooled ρ = 0.45, 95% CI 0.31–0.58). Three hours of hyperinsulinaemic clamp reproduce the
transcriptional change of three hours of morning, whereas hydrocortisone against placebo at the
same hour of day moves PER1 in every participant, in the direction opposite to insulin. A
forced-oscillator model in which the effective forcing is the product of the cortisol rhythm and
receptor availability reproduces the signature with no fitted parameter.

## Reproducing the results

```bash
pip install -r requirements.txt        # or requirements-lock.txt for exact versions
# place the public downloads listed in docs/DATA.md under data/raw/
bash run_all.sh
```

`run_all.sh` runs 29 steps in order, grouped by the figure they feed, and ends with three checks:
that every input a figure reads is produced by the pipeline or versioned in `data/reference/`,
a smoke test, and a verification of the output package. Figures land in `results/figs` and tables
in `results/new`.

Analyses use Python 3.9 or later; stochastic routines use NumPy `default_rng` with seed 0.
Continuous integration runs on 3.11.

## Layout

| path | contents |
|---|---|
| `python/analyses/` | the analysis steps, numbered in execution order |
| `python/model/` | the forced-oscillator model and the data it generates |
| `python/figures/` | the five main figures and the seven supplementary figures |
| `python/lib/gene_sets.py` | the single definition of every gene set used anywhere |
| `data/reference/` | inputs that are not regenerated here, and why (`docs/DATA.md`) |
| `tests/` | syntax, dependency chain, smoke test, package verification |
| `manuscript/` | the manuscript and its supplementary information |

## Two conventions worth knowing before reading the code

**Gene sets live in one place.** `python/lib/gene_sets.py` is the only definition; scripts and
supplementary tables import from it and never redefine a list in place. The clock-output set is
the output arm and excludes the core negative-feedback components, which are scored separately.

**Expression is centred per sample for scores, and not for group differences.** Centring removes
a common component that otherwise inflates every correlation in RNA-seq (the median correlation of
a random gene with the clock-output score is 0.44 rather than zero without it). The signature table
compares differences between groups and is built without centring. Mixing the two conventions
changes results.

## Citation

See `CITATION.cff`. The archived release is doi:10.5281/zenodo.23379059.
