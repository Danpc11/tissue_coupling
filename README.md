# Independent hormonal and circadian signals coordinate human metabolic tissues

Analysis code for the manuscript of the same title. Everything reported in the paper is produced
by `run_all.sh` from public data.

## What the analysis asks

Among the genes expressed in both skeletal muscle and subcutaneous adipose tissue, which ones place
an individual consistently relative to others in **both** of their organs at once, and through which
inputs? Cross-tissue concordance measures how systemic an input is: the genome is perfectly systemic,
circulating signals partly so, local programmes not at all.

## Reproducing the results

```bash
git clone https://github.com/Danpc11/tissue_coupling.git
cd tissue_coupling
pip install -r requirements.txt
python python/analyses/00_prepare_data.py     # lists what to download and from where
bash run_all.sh                               # runs everything and verifies the output
```

GEO and GTEx require accepting their terms, so the data are not redistributed here and are not
downloaded automatically. `00_prepare_data.py` prints exactly which files are needed and whether
each is present.

```
python/
  lib/gene_sets.py    single registry of every gene set, with sources and a synonym table
  lib/phase.py        circadian phase from a single sample, with its three controls
  lib/enrichment.py   set-level enrichment against matched null models
  analyses/00..03     the four steps, in order
  figures/            journal style and figure helpers
tests/verify_package.py   coherence check over the whole result set
```

## Three design decisions worth knowing before you read the code

**Gene sets live in one place.** `python/lib/gene_sets.py` is imported by every script and by the
supplementary tables, so code and tables cannot diverge. It carries a synonym table because
deposited matrices keep the annotation of their platform: four genes were wrongly declared absent
until we resolved current HGNC symbols against a 2010 array, which moved the sex component from
0.559 to 0.590.

**Two null models, both reported.** Random sets matched for size, and random sets additionally
matched on expression level and variability. They disagree for the glucocorticoid and insulin sets,
and we report both rather than the favourable one. Highly expressed and variable genes are better
measured and therefore correlate better between tissues; the matched null removes that advantage and
some real signal with it.

**Phase is a pseudo-quantification.** The unsupervised ordering defines phase only up to rotation and
reflection, so absolute hours are not interpretable and values are not comparable between cohorts
with different sampling windows. All comparisons use rotation-invariant statistics, and two
orderings are aligned before being compared — omitting that step inflates the error estimate
fivefold and makes it exceed the signal.

## The verification step is not optional

`tests/verify_package.py` checks that every result file declares its cohort and platform, that
identifiers match that platform, that declared set sizes match the registry, that the pooled
meta-analysis reproduces from the per-cohort exports, that Benjamini–Hochberg values are correct and
monotone, and that all files come from one run. It exists because these checks caught real errors:

- a supplementary table drawn from the wrong cohort (Affymetrix probes labelled as the Illumina cohort)
- a cross-tissue prediction inflated by the same genes appearing on both sides (R² 0.14 → 0.47)
- four genes wrongly declared absent through symbol drift
- an apparent treatment effect that random non-rhythmic gene sets reproduced (P = 0.408)

A fifth check lives in `phase.is_clock`: the recovered cycle must reproduce the canonical BMAL1–DBP
antiphase. Blood, hypothalamus and cerebral cortex fail it. Without that filter blood would have been
reported as the best-coupled tissue, which is the opposite of what it is.

## Relation to the earlier repository

This repository supersedes `Danpc11/t2d_landscape`, which contains the earlier analysis of
directional coherence in the muscle response to insulin. That hypothesis depended entirely on one
cohort and did not replicate; it is retained there as a record of what was tried and withdrawn.

## Citation

Repository: https://github.com/Danpc11/tissue_coupling

[Manuscript reference and archived Zenodo DOI to be added on acceptance. Connect the repository to
Zenodo **before** creating the GitHub release, or the release will not be captured.]

## Licence

See `LICENSE`.
