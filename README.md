# Glucocorticoid signalling tracks circadian output across human metabolic tissues

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23379059.svg)](https://doi.org/10.5281/zenodo.23379059)
[![tests](https://github.com/Danpc11/tissue_coupling/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/tissue_coupling/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Code for the manuscript. Everything in the paper is produced by `run_all.sh` from
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

`run_all.sh` runs the pipeline in order, grouped by the figure they feed, and ends with three checks:
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

## Branch `repro-hugo` (November 2026)

Reconstruction of the data-preparation step that the published package lacked, so that `run_all.sh` regenerates every result and figure from the public downloads: `download_public.sh`, `python/analyses/00b_build_matrices.py` (expression matrices), `00c_archs4_counts.py` (fallback counts from ARCHS4), `02e_gse87761_benchmark.py` (the muscle template against recorded sampling times in adipose tissue, GSE87761), `00d_reference_tables.py` (the tables formerly versioned in `data/reference/`, now regenerated and compared with the versioned copies) and `01c_islets_deseq2.py` (islet differential expression); `05_clamp_within_person.py` gains block E, the clamp response by insulin sensitivity in muscle and adipose tissue (`results/new/clamp_response_by_sensitivity.tsv`). Sample centring is applied in every cross-individual score, as Methods states. See `docs/DATA.md`. The branch was handed over to the corresponding author in November 2026; the reproduction report that compares each Results statement with the branch's numbers is kept with the manuscript files, not in this repository.
