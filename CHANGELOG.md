# Changelog

## v2.10.0: a third audit of the whole repository, code and documentation

Every script and every document read again as a whole, and the pipeline run from the first step to
the last on the real data with numerical warnings no longer silenced.

### Fixed: things that could break a run or read the wrong file
- **The step that builds the matrices stopped the pipeline** if a link to the GTEx annotation
  already existed and pointed nowhere. The link was a shortcut relative to the working directory,
  and where it did resolve it could make two steps read a file other than the one in `T2D_RAW`.
  The same shortcut existed for MSigDB. Both are gone: every input is read from one place.
- **Two routes computed the response to hydrocortisone**, one of them by loading an old module
  through `exec`. They agree to 2e-16; one is kept and the module (`02_meta_and_validation.py`,
  which also carried a superseded copy of the cross-tissue prediction) is removed.

### Fixed: method
- **Supplementary Fig. 4 scored transcription factors against genes of the score itself.** A factor
  whose regulon contains genes of the glucocorticoid score correlates with it by construction. Step
  17 already excluded clock genes for that reason; this figure did not exclude the glucocorticoid
  ones. The androgen receptor was in its top five at 0.35 because of two such targets; without
  them it is 0.25 and leaves. FOXO3, FOXO4, POU2F2 and RFX5 do not change.
- **Supplementary Fig. 1 still took RIN and ischaemic time from another sample** of the donor, the
  defect corrected in step 09 in v2.3.0 and not in the figure that shows the same quantities.
- **Meta-regression *P* values used the normal distribution** with thirteen estimates. They use the
  t distribution with k − 2 degrees of freedom.
- Supplementary Fig. 7 carried two more copies of the model equations. It calls the model; the
  values are identical to the last bit.

### Fixed: what the documentation said and the code did not do
- **Two series were downloaded and required, and no step read them**: GSE104674 and GSE11868. They
  are in Supplementary Table S1 with a role each. The analyses use seventeen of the nineteen series
  of that table, plus GSE40589. They are no longer downloaded, and `docs/DATA.md` says so.
- **Seventeen settings were undocumented.** The README has a table of every environment variable
  the code reads, with its default checked against the code.
- The README said every stochastic routine uses seed 0, which stopped being true with several
  seeds; it said nothing of the DOI badge pointing to v1.0.0, a release that cannot regenerate the
  results; and its layout table omitted `python/lib/`.
- Descriptions of cohorts in the input check were out of date: GSE41168 has 35 people per tissue
  and 28 with both, not 28; GSE182117 has 12 donors, not 11; GSE78721 is used one sample per
  person; and the sleep-delay series was called a control "for the estimator".

### Removed
- Three files that should never have been versioned: two `.DS_Store` and a pickle of GTEx
  phenotypes that nothing read.
- Two declared dependencies that no script imports, `rdata` and `openpyxl`.

### Changed
- **Numerical warnings are no longer silenced.** Thirty-four scripts ignored every warning. They now
  ignore only deprecation notices from the libraries, so a division by zero or a non-finite value
  shows. A complete run emits three kinds, all accounted for: one inside PyDESeq2, the ones CHIRAL
  raises (which step 06 already counted and reported) and a 0/0 in a diagnostic ratio, now handled
  where it happens.

### Documented
- CHIRAL uses NumPy's global random generator, which is seeded before each call; the README said
  every routine used its own generator.
- The optional ARCHS4 route for the islet counts is not equivalent to the NCBI counts the reported
  numbers use, and nothing in the outputs records which was used. `docs/DATA.md` says so.

### Validation
The pipeline was run from step 00 to step 30 on the real data with this release: the thirteen
estimates, eleven figures and the clamp figure, the ten tests and `verify_package` all pass, and
the whole run emits two numerical warnings. Supplementary Fig. 4 draws 16 targets instead of 25
once the genes of the score are excluded.

### Not changed, and needs a decision
- The order of authors in `CITATION.cff` differs from the manuscript in the second and third places.

## v2.9.1: the calibration of the joint criterion, made independent

Four points of a further review of `python/model/competing_models.py`. All four were right.

### Fixed
- **The calibration drew its cohorts from the same simulations that defined the prediction.** The
  two were not independent, so the coverage it reported was not an independent check. Calibration
  cohorts now come from people integrated separately, with seeds that no prediction uses, and from
  two independent populations, so that the table gives the range of the coverage and of the
  calibrated threshold between them. The largest threshold is the one used.
- **Calibration cohorts could be smaller than the observed cohort.** Their size was capped at half
  the pooled simulations. Each population now has three times the size of the observed cohort,
  whatever that is, and the cohorts drawn from it have exactly that size; the table records both.
- **`region()` inverted whatever it was given.** It now refuses, saying why, a matrix with
  non-finite values, with no direction of variance, that is not square or that is not positive
  semidefinite; a negative eigenvalue of rounding size is treated as zero; and the condition number
  it returns is that of the directions it keeps.
- **A cache could be reused with code that no longer produced it.** `T2D_MODEL_REUSE=1` now reuses
  the integrations and the simulated table only if their manifest matches: a hash of the source of
  every function that generates them and of all their parameters (people, seeds, grid, windows,
  models). Otherwise they are generated again, and the step says so.

### The result with the independent calibration (thirty seeds)
Coverage of the nominal 95% threshold, by cohort: 99% in AAGMEx (98 to 100 between populations),
and between 87% and 98% in the twelve GTEx cohorts and strata. Where the calibrated threshold is
above the nominal one it is used (up to 17.1 against 12.6).

The conclusions do not change. AAGMEx: no forcing 0 of 18 grid points (closest distance 63.5),
cortisol alone 0 (22.4), cortisol times receptor 2 (6.3), receptor alone 2 (8.3). GTEx, adipose
tissue and muscle, ventilator donors: cortisol times receptor is compatible in both adipose depots
(2.0 and 3.3); receptor alone is compatible in one (10.6) and just outside in the other (14.0
against a threshold of 13.7). Fast deaths: the three forced models in the subcutaneous depot, and
all but cortisol alone in the visceral one.

`tests/test_model_verdict.py` checks each of the four corrections.

## v2.9.0: the comparison of mechanisms, calibrated; and why GTEx looked incompatible

### GTEx is not one cohort, and the sampling window is not known there
In v2.6.0 to v2.8.1 no model reproduced GTEx. Two assumptions were behind that, and neither holds.

- **The circumstances of death decide how and when tissue was obtained.** In donors on a
  ventilator (Hardy 0) procurement is scheduled surgery; in fast deaths (Hardy 1 and 2) it happens
  at whatever hour. Split that way, the data move in the direction the forced models predict when
  the sampling window widens: receptor against output is 0.35 (0.27 to 0.42) in ventilator donors
  and 0.18 (0.07 to 0.29) in fast deaths, and the concordance of output between adipose tissue and
  muscle is 0.46 (0.38 to 0.55) and 0.71 (0.64 to 0.78). GTEx as a whole is the mixture.
- **The hours of collection are not public, and there is no reason to assume them uniform.** Where
  the design does not fix the sampling window it is now one more unknown of the grid (centre and
  spread), and the result says which window each model would need. Under uniform hours no model is
  compatible with any stratum (closest distances 45 to 650). With the window free, for adipose
  tissue and muscle: ventilator donors are compatible with cortisol times receptor at a narrow
  window (08:00, spread 1.5 h; distance 2.0 and 3.3 in the two adipose depots) and not with no
  forcing or cortisol alone; fast deaths are compatible with the three forced models at a wider
  window (spread 3 to 4.5 h). The model therefore makes a prediction that someone with access to
  the protected GTEx data can test: collection times in ventilator donors are concentrated in a
  few hours of the day.
- Not everything is reproduced. The pair of two adipose depots is compatible with no model: its
  receptor concordance (0.51) is above what the grid reaches. And compatibility is thin: one or two
  grid points of 306 for ventilator donors.
- The phase estimator cannot settle the window: it spreads over the whole day the samples of a
  cohort known to be taken in the morning.

### Fixed, in the comparison of mechanisms
- **The joint region is calibrated instead of assumed.** For each cohort, cohorts of the same size
  are drawn from the closest model and put through the same procedure as the data. The nominal 95%
  threshold covers between 91% and 98% (thirteen cohorts); where the 95th percentile of those
  distances is above the nominal threshold it is used instead. The covariance has full effective
  rank, six of six, in every cohort, and the code now uses the effective rank whatever it is.
- **No forcing is evaluated on the full grid**, like the other three. Its receptor concordance does
  change over the grid, as the review said; it remains compatible with nothing (closest distance
  63.5 in AAGMEx).
- **Thirty seeds, and the answer does not depend on how many.** AAGMEx with 3, 10 and 30 seeds:
  cortisol times receptor 3, 2 and 2 compatible grid points of 18 (closest distance 7.1, 6.3, 6.3);
  receptor alone 4, 2 and 2 (11.0, 9.6, 8.3); cortisol alone 0 (13.9, 22.4, 22.4); no forcing 0.
  The default is ten. `competing_models_seed_stability.tsv` has every cohort.
- **The per-observable verdict uses the same error model as the joint criterion**: the state of a
  prediction comes from its standard error across seeds, and a grid point counts only if the
  prediction has the state of the observed value and falls inside its interval, both at that point.
- The calibration computed each set of observables six times over.

### With all of that, AAGMEx
No forcing and cortisol alone are excluded. Cortisol times receptor and receptor alone are both
compatible, at the same two grid points; the product is closer (6.3 against 8.3) because receptor
alone falls short of the concordance between tissues (0.36 and 0.32 against 0.43 and 0.39 observed).
The data do not separate them.

### Two ways of measuring amplitude, said in one place
Cosinor where a study has a series with recorded hours in each person; the single-sample estimator
where there is one biopsy per person and the tissue passes the antiphase check. The sleep-delay
control in blood is of the first kind, and its step says so instead of reading as a failed
validation: cosinor reproduces the published effect (8 of 8, *P* = 0.008), and blood does not pass
the antiphase check, so the single-sample estimator is not the tool there.

### Changed
- `forced_clock.simulate_cohort()` integrates a whole cohort at once, with the same equations; it
  is six times faster than person by person and gives the same ensemble amplitude and synchrony.
- Step 20c takes about twenty minutes on one core with the default ten seeds; integrations run in
  parallel (`T2D_MODEL_WORKERS`), and `T2D_MODEL_REUSE=1` reuses them.

## v2.8.1: the joint criterion for the competing models

The criterion introduced in v2.8.0 counted the simulations in which every observable fell inside
the observed interval. It had four defects, and the two tables it wrote are recomputed.

### Fixed
- **It mixed cohorts.** For the post-mortem window the "observed interval" of each observable ran
  from the lowest lower bound to the highest upper bound of three cohorts, which is the interval of
  nothing. Every comparison is now with one cohort, with its own estimates and uncertainty.
- **It used six separate intervals as if they were a joint region.** The observables are strongly
  correlated (receptor against output, with and without the other receptor removed, are nearly the
  same quantity), and six marginal 95% intervals do not cover 95% jointly. The distance from a
  prediction to the observed values is now measured with the covariance of the six, from a
  bootstrap of people on Fisher's z scale, against the 95% quantile of a chi-square with six
  degrees of freedom. The fraction inside simultaneous intervals is given beside it.
- **It left the receptor concordance between tissues out.** It was described as used "to fix the
  model", and nothing fixed anything with it: simulations in which the receptor of the two tissues
  is unrelated were scored against a cohort in which it is related. It is now one of the six.
- **It compared single simulated cohorts**, each a noisy realisation of what a model predicts, so
  that sampling noise decided. The prediction of a model at a point of the grid is now the mean of
  its seeds, with the simulation's own covariance reduced accordingly. More seeds make the
  criterion stricter, where before they made it looser.

### The result with the corrected criterion (three seeds, 600 people per grid point)
- **AAGMEx, morning biopsies.** No forcing: 0 of 3 grid points, closest distance 52.5 against a
  threshold of 12.6. Cortisol alone: 0 of 18, closest 15.2. Cortisol times receptor: 4 of 18,
  closest 2.1; at that point the six predicted values are 0.20, 0.42, 0.16, 0.13, 0.19 and 0.39,
  and the observed ones 0.20, 0.43, 0.22, 0.10, 0.19 and 0.39. Receptor alone: 2 of 18, closest
  4.6; it falls short of the concordance between tissues (0.34 and 0.28 against 0.43 and 0.39).
  The first two mechanisms are excluded; the product fits best and receptor alone is not excluded.
- **GTEx, post mortem.** No grid point of any model is compatible with any of the three tissue
  pairs; the closest distances are between 169 and 322.

`competing_models_joint.tsv` now has one row per cohort and model, with the number of compatible
grid points, the closest distance and the grid point that reaches it; `competing_models_verdict.tsv`
has one row per cohort, observable and model. `tests/test_model_verdict.py` checks each of the four
corrections on values built by hand.

## v2.8.0: one annotation rule for the meta-analysis, and a second audit of every script

### Fixed
- **GSE13070 entered the meta-analysis annotated differently from the other cohorts of its
  platform.** GSE25462 and GSE41168 are collapsed to genes by Entrez identifier keeping the probe
  with the highest mean, as the Methods say. GSE13070 came in through the loader of the clamp
  analysis, by platform symbol and averaging probes. Its two estimates, the highest of the
  thirteen, were 0.1 higher for that reason: 0.78 and 0.67, against 0.69 and 0.58 with the rule of
  the others. With one rule for all, **the pooled estimate is 0.415 (0.330 to 0.494) instead of
  0.443, and the three-level estimate 0.457 (0.284 to 0.602 with the t distribution) instead of
  0.484**; heterogeneity falls from 86% to 82% and the prediction interval is −0.09 to 0.79. The
  specificity result does not move: *NR3C1* at the 99.93rd percentile of 14,793 genes, second of
  44 nuclear receptors after *NR3C2*. `T2D_GSE13070_SYMBOLS=1` restores the previous annotation.
- **The waveform of a fixed reference was fitted by simple averages**, which equal the
  least-squares coefficients only when the phases of the reference cohort are spread evenly over
  the cycle. With a cohort covering half the cycle the phase of a new sample was off by up to
  3 h; fitted by least squares with an intercept it is off by 0.02 h. In the references of the real
  data the phases are close to uniform and the effect is small. The averages are kept so that the
  published values can be regenerated (they are, to the last bit).
- **The verdict on the competing models could not see a negative association**, and asked whether a
  model could produce a feature "at some point of the grid", which more seeds make ever easier.
  It now has three states, uses bootstrap confidence intervals of the observed values, and reports
  how often: the fraction of the simulations of each model in which every observable falls inside
  the observed interval at once. The four models now see the same people, the same biopsy hours
  and the same measurement noise, over three seeds.
- The noise analysis of the amplitude carried its own copy of the estimator, with the geometry
  corrected in v2.7.0 still uncorrected there. It uses the library.
- The entrainment sweep carried its own copy of the model equations. It calls the model; its output
  is identical to the last bit.
- The clamp analysis took "the first file whose name begins with GPL570". It takes a named file.
- The library module for within-person coupling still described, with numbers, an analysis that is
  no longer in the paper, and four functions that nothing called. They are gone.
- Fig. 2: two sample sizes were written by hand and one axis label was in Spanish.

### Added
- **Is the association only cases against controls?** For the seven estimates whose cohort has a
  diagnostic or study group, the meta-analysis also gives the correlation with the group as a
  covariate. It hardly moves (0.53 against 0.54 in GSE78721, where it is 0.56 in the normal group
  and 0.57 in the diabetic one), and the pooled estimate with those adjustments is 0.434.
- **Transcription-factor activities against random regulons** of the same size. *FOXO3* stands
  (z = 3.7, *P* = 0.002) and so does NF-κB at the other end (z = −4.3). The regulon of *NR3C1*
  itself does not: fifth of 134 by rank, but z = 1.7 and *P* = 0.08.
- **How many genes of the hydrocortisone signature are expected by chance.** With ten participants
  a gene with no effect changes in the same direction in all of them once in 512 times, and there
  are 29,525 genes. Flipping the sign of whole participants, 27 genes pass both criteria by chance
  (95th percentile 54) against 131 observed: the signature is real as a whole (*P* = 0.004) and
  about one gene in five of it is expected to be chance.

### What the competing models say with the corrected criterion
With morning biopsies (AAGMEx), the fraction of simulations in which all five observables fall
inside the observed intervals at once: no forcing 0 of 9; cortisol alone 0 of 54; cortisol times
receptor 6 of 54; receptor alone 3 of 54. The first two are excluded; the product is favoured over
receptor alone but not separated from it. With post-mortem samples (GTEx) it is 0 for every model.

## v2.7.0: the geometry of the fixed reference, and every figure drawn from the run

Answers to a further external review. One of its points referred to code already changed (Fig. 5c
has had separate axes and no ratio of magnitudes since v2.5.0); the others were right, and looking
into them turned up one more case of the same kind.

### Fixed
- **The fixed reference projected onto vectors that are not orthonormal.** The waveform of a
  reference is z = a cos θ + b sin θ with one vector a and one b over genes. Phase was taken as
  atan2(z·b, z·a) and amplitude as the projection on a cos θ + b sin θ, which is only right if a
  and b are orthogonal and of equal norm. In the references fitted to the real data they are close
  to orthogonal (93° and 83°) but not of equal norm (ratio 1.58 in myotubes, 1.13 in muscle in
  vivo): the phase was biased by up to 1.7 h and the gain of the amplitude changed with phase,
  between 0.76 and 1.21 times its mean, in myotubes. `infer_from_reference()` now solves the two
  coefficients by least squares and takes the amplitude as the norm of the projection onto the
  plane of a and b, the same definition as the cohort-relative `amplitude()`; it rejects a
  degenerate reference. The published geometry is kept as an option so that the published values
  can be regenerated.
  **The validation holds with the corrected estimator**, and is cleaner with the whole donor
  excluded: myotubes ρ = 0.35 and 0.34 per sample and 0.94 and 0.95 per donor (0.32, 0.34, 0.99
  and 0.90 as published, on the same matrix); muscle in vivo ρ = 0.20 per sample (*P* = 0.033) and
  0.81 per series (0.18, *P* = 0.053, and 0.80).
- **Fig. 1f and 1g were drawn from versioned tables**, not from the validation computed in the run.
  They now read the run, and plot the corrected estimator with the whole donor excluded.
- **Fig. 3h and 3i were drawn from a table with nine versioned columns.** Step 16 now writes two
  named tables, `signatures_published_reference.tsv` and `signatures_fully_regenerated.tsv`, their
  difference gene by gene, and whether they lead to the same conclusion (they do: ρ = 0.71 and
  0.70 with hydrocortisone, *P* = 0.026 and 0.016 against random gene sets). The figures read the
  fully regenerated one, with nine output-arm genes instead of eight.
- **The tests compared phases by the concentration of their difference**, which is 1 for two
  phases 12 h apart. They now use the circular error: a fixed-reference phase is within 0.45 h
  (median) of the true hour of a simulated sample, and a constant 12 h offset fails.
- **Phases were compared with a linear correlation** in the equivalence table. They are compared
  as circular quantities, allowing a global offset and either direction. The answer does not
  change and is now unambiguous: for `amplitude_phase_distributions_two_cohorts.tsv` the two builds
  give unrelated phases to the same sample (median error 4.6 to 6.3 h after removing the offset,
  where chance gives 6 h). In cohorts biopsied at one hour of the day the true phases differ
  little and the per-sample phase is mostly noise.
- **Cross-tissue prediction: an R² on one scale.** Besides the R² of each fold, the model is also
  fitted to the outcome without dividing by the fold's standard deviation, so that all samples are
  in log2 units: 0.47 (0.44 with fold-wise standardisation; 0.43 on average per fold), against
  0.03 with covariates alone and a 95th percentile of 0.25 for random predictor genes.

### Changed
- The README says "recomputes" and explains the difference from reproducing published numbers.

## v2.6.0: the three analyses that could refute the paper, and a full run on the real data

v2.5.0 was run here from the public downloads to the last check, on another machine: 46 steps, the
thirteen estimates, every figure, `verify_package` clean. The pooled estimate is 0.443 (0.348 to
0.529), as on the authors' cluster. This release adds three steps, each run on the real data.

### Added
- **Four mechanisms against the same observables** (`python/model/competing_models.py`, step 20c).
  The model of the paper has its conclusion written into it: the forcing is the cortisol rhythm
  times receptor availability. Here it competes with no forcing, cortisol alone and receptor alone,
  with the same oscillator, noise, sampling and estimator, over a grid of what is not known (how
  much cortisol varies between people, how alike the receptor is in two tissues, how noisy mRNA is
  as a measure of it). Six observables are computed in the same way in the simulations and in the
  cohorts with two tissues per person.
  With morning biopsies (AAGMEx): no forcing produces nothing; cortisol alone produces the
  concordance between tissues but no relation with the receptor (−0.08 to 0.13, observed 0.20 in
  that cohort); cortisol times receptor reproduces every feature; receptor alone reproduces the
  signs but not the size of the concordance between tissues once the receptors are removed
  (0.07 to 0.31, observed 0.39; the product gives 0.23 to 0.68).
  **With post-mortem samples (GTEx) no model reproduces the data.** With biopsies at any hour all
  three forced models predict that the relation between receptor and output vanishes (−0.06 to
  0.06) and that the two tissues agree almost perfectly (0.9); the data give 0.28 to 0.32 and 0.56
  to 0.61. Either deaths in GTEx are not spread over the 24 hours, or the association has a
  component of level that does not depend on the hour, which a model that only moves amplitude
  does not contain.
- **The phase template against ZeitZeiger and tauFisher** (`32_benchmark_methods.py` and
  `R/benchmark_supervised.R`, step 02f), on the same data and the same folds.
  Leaving subjects out in GSE108539: template 1.53 h mean and 0.92 h median error; ZeitZeiger 1.95
  and 1.11; tauFisher 2.32 and 2.00; chance 6.05.
  Trained in muscle, tested in adipose tissue on another platform (GSE87761): template 3.03 h mean
  and 2.10 h median with one offset learned from the other subjects; ZeitZeiger 2.93 h with no
  offset at all once each gene is standardised within its cohort, 6.50 h without that step;
  tauFisher 4.98 h. The template is on a par with ZeitZeiger and neither is better.
- **What the amplitude estimator can and cannot separate** (`python/model/factorial_estimator.py`,
  step 20d). Cell amplitude, phase dispersion between cells, fraction of rhythmic cells and noise
  are varied independently. The inferred amplitude follows the collective amplitude (ρ = 0.85) and,
  once that is removed, nothing of its three components is left (partial ρ of −0.02, 0.02 and
  −0.02): two samples with the same collective amplitude are indistinguishable, however they got
  there. A change of composition that does not touch the clock does not bias it. Without any rhythm
  it still returns a number, half the size, and **the BMAL1–DBP antiphase check lets one cohort in
  three through** (32% of 200), as it must when two unrelated phases are compared with a threshold
  of 8 h out of 12.

### Fixed
- **Fig. 2a printed "P = 0e+00"** and the pooled table stored a P of exactly zero: one minus the
  normal distribution function loses all precision for a large z. The tail is used instead and the
  label reads "P < 10⁻¹⁵".
- `verify_inputs.py` now sees the files an R script writes.

### What the full run settles (all in `results/new/` after a run)
- Specificity with the thirteen estimates: *NR3C1* at the 99.93rd percentile of 14,127 genes; ten
  are above it, among them *BHLHE41*, *CRY2* and *NR3C2*; second of 44 nuclear receptors, after
  *NR3C2*.
- The signature table regenerated in full: insulin resistance against hydrocortisone, ρ = 0.70 over
  nine output-arm genes (*P* = 0.016 against random sets of nine), against 0.71 over eight with the
  versioned columns. Against Cushing, 0.67 (*P* = 0.065).
- Nested leave-subject-out error of the template on the pipeline's own matrix: 0.92 h median.

## v2.5.0: answers to a second look by the same reviewer, and one new control

This release contains v2.4.0, which was never published on its own.

Two of the reviewer's points refer to code that was no longer in the repository: the leak in the
validation of the phase template was removed in v2.2.2 (`02b_phase_template.py` has 256 lines and
no selection on all samples), and the fixed reference for a sample outside the cohort was added in
v2.3.0. The other four were right.

### Fixed
- **The benchmark of the template in adipose tissue led with an optimistic error.** The template
  fixes phase up to an arbitrary zero, so one offset has to be learned from samples with a known
  hour. The headline figure fitted it on the very samples it then scored. The primary figure is
  now the one with the offset learned from the other subjects, for every gene list, and chance is
  measured with the same procedure (6.0 h). On the real GSE87761 with the versioned template:
  2.64 h mean and 2.33 h median, 67% within 3 h, against 2.62 h with the optimistic offset.
- **The cross-tissue prediction pooled folds that are on different scales.** The outcome is
  standardised within each fold with its own training set, which is the only way without a leak,
  so the pooled R² mixes scales. The table now also carries the R² of each fold, which does not:
  0.43 on average, from 0.39 to 0.48, against 0.44 pooled. And the covariates of the
  covariates-only model were not standardised before the penalised regression while the genes
  were; they are now, and that row moves from 0.033 to 0.036.
- **The sleep-delay control measured a different amplitude from the estimator it was said to
  validate.** It fits a cosinor on the recorded hours. The step now also applies the single-sample
  estimator under the same design, and says what it finds: the cosinor amplitude falls in 8 of 8
  participants (median −28.6%, *P* = 0.008); the inferred amplitude falls in 3 of 8 (median +8.9%,
  *P* = 0.95) and does not follow the cosinor change across participants (ρ = −0.33). Blood does
  not pass the antiphase check the estimator requires of a tissue (BMAL1 and DBP 2.3 h apart). The
  control supports the gene set, not the single-sample estimator.

### Added
- **Is the association specific to the receptor?** (step 12b). The same partial correlation, with
  the same covariates and the same random-effects pooling, for every gene measured in all the
  estimates. On the eleven estimates that could be built here, *NR3C1* is at the 99.97th
  percentile of 15,041 genes: four are above it, among them the clock gene *BHLHE41* and *NR3C2*,
  the mineralocorticoid receptor, which also binds cortisol. Among the nuclear receptors *NR3C1*
  is second to *NR3C2*, and the third, *RORA*, is at 0.29 against 0.45. Its percentile within each
  estimate has a median of 98 and a minimum of 77.
- The meta-analysis table carries, for each estimate, the correlation of each predictor with clock
  output once the other is removed: the receptor given the glucocorticoid score, and the reverse.
- A fixed reference for the anchored phase as well (`fit_anchored_reference`,
  `phase_from_anchored_reference`), tested: a new sample gets the phase it would get inside the
  cohort, whichever reference is used.

## v2.4.0: the last scripts read, and what the figures draw

The scripts not read in v2.3.0 (perturbations, signature table, case-control analysis, islets,
phenotype terms) and every figure script. Steps were rerun on the real data where the downloads
reach, and all eleven figures were rendered with the edited code.

### Fixed
- **The enrichments scored the receptor against itself.** Fig. 4b, 4c, 4f and 4g rank every gene
  by its correlation with *NR3C1* and score gene sets on that ranking. The correlation of *NR3C1*
  with itself is 1, and it was left in the sets that contain it: 7 of the phenotype terms shown,
  and 189 sets of the Hallmark and phenotype collections. It is circular, and it weighs in small
  sets. With the receptor excluded, "increased urinary cortisol level" (16 genes) goes from z = 3.8
  to 2.0 in adipose tissue and from 1.9 to −0.2 in muscle; "truncal obesity" from 3.7 to 2.8 and
  from 0.9 to 0.0; the GO terms of response to dexamethasone from +1.4 to about zero. Hallmark
  changes little (three sets, by 0.5). `z` now excludes the receptor, `z_including_anchor` (or
  `z_including_NR3C1`) keeps the previous value, and sets without it are unchanged to the last bit.
- **Step 15, "hydrocortisone, dexamethasone, Cushing", computed nothing.** `04_perturbations.py`
  printed two lines and exited; its docstring named an output that no code writes. It is a library
  of functions used by two other steps, says so now, and is no longer a step.
- **The signature table is not fully regenerated, and the documentation said it was.** Nine of its
  twelve columns, among them the insulin-resistance contrast, come from the versioned copy. That
  is now stated in the script and in `docs/DATA.md`, and step 16 computes the same statistics on
  the fully regenerated table too.
- **Fig. 5c compared the model in its own units with observed log2 differences** on one axis, and
  its title gave their ratio as an agreement "within 14%". The magnitude of the model follows from
  an oscillator amplitude fixed at 1 and from the contrast of receptor gain chosen. Each series now
  has its own axis and the title gives what is comparable: the signs, and the correlation of the
  pattern across genes, which with the model's hand-set phases is 0.14.
- **Fig. 5e drew a band "observed with NR3C1 (0.09–0.22)"** whose two numbers no step produces.
  Removed until they have a producer.
- **Numbers written by hand in the figures.** Every count in a panel annotation is now read from
  the data the panel draws: participants, cohorts, genes, transcription factors. One was already
  stale (34,020 genes where the ranking has 34,615) and one was wrong (nine clock genes where
  eight are drawn). Three annotations that stated a conclusion now describe the panel.
- The case-control analysis added each cohort only if its matrix existed, and printed published
  values that no longer hold. Its five cohorts are a manifest, as in the meta-analysis.
- A failed supplementary figure was reported and the step still succeeded.
- The islet donors were selected by a GEO accession threshold in two steps and by reading the
  series matrices in a third; all three read the series matrices now.

### Added
- **The model's prediction without units, with measured phases** (step 20b). The model says that a
  biopsy at a fixed hour sees each gene pulled towards its mean in proportion to its amplitude.
  With the phases and amplitudes measured in the time-course cohort instead of the ones written in
  the model, the predicted sign agrees in 8 of 8 genes and the predicted pattern correlates with
  the observed one at ρ = 0.74 (P = 0.027 against random template genes) at 08:00, the hour of the
  model; 0.31 at 06:00 and 0.86 at 10:00.
- **The signature correlation against its background** (step 16). Over the eight output-arm genes
  the insulin-resistance and hydrocortisone profiles correlate at 0.71; over all 6,108 genes at
  −0.05; random sets of eight genes reach 0.62 at their 95th percentile (P = 0.026); leaving one
  gene out gives 0.57 to 0.79. The Cushing profile does not stand out from random sets (P = 0.16).
- Panel titles that were passed to a function that ignores them, sixteen of them stating results,
  are gone from the code.

## v2.3.0: a full review of every script, run on the real data

The published code was run here on the real downloads, as far as they reach, on a different
machine and with a different source of annotation. GTEx, the cross-tissue axis, the
transcription-factor activities, the clamp scores, the cross-tissue prediction and the sleep-delay
control come out identical to the third decimal; of the eleven cohort estimates of the
meta-analysis that could be built, nine are identical and two differ by 0.002 and 0.003.

### Fixed
- **The sensitivity to principal components ran without covariates in AAGMEx.** It indexed sex
  and age by GEO accession while the matrix columns are donors, found no match, printed a notice
  and went on. Its first row is now 0.345, the value of the meta-analysis, instead of 0.401; the
  conclusion (the association disappears when two components are removed) is the same.
- **The meta-analysis pooled whatever estimates it found.** Four of the thirteen were added only
  if their matrix existed, so a run could silently pool eleven. The thirteen are now a manifest:
  a missing one stops the step, `T2D_ALLOW_PARTIAL=1` lets it continue and marks every table as
  incomplete.
- **GTEx RIN and ischaemic time came from another sample** of the same donor and tissue in 57, 27
  and 81 donors, the first one the annotation lists rather than the one whose expression is in the
  matrix. With the right sample the adjusted estimates are 0.277, 0.352 and 0.277 instead of
  0.276, 0.354 and 0.274.
- **`verify_package` could not have passed on a real run.** Its check of the meta-analysis had
  never been executed against real output: it crashed on the pooled table and treated a standard
  error already on Fisher's scale as if it were not; and it looked for the out-of-fold predictions
  in the wrong directory. Rewritten and run against real tables.
- The step for phenotype terms read MSigDB from a fixed path, and one module was loaded by a path
  relative to the working directory: both would fail with the data in another location.
- A heterogeneity of "84%" was written by hand in a message (it is computed now), two messages
  quoted values from the manuscript, and a docstring carried a table of results.
- `is_clock()` tested a separation "between 8 and 16 h" that can never exceed 12. It now says
  what it does, within 4 h of antiphase, with the same result, and finds BMAL1 under ARNTL.
- The dispersion analysis called a permutation test a bootstrap, could divide by a dispersion of
  zero, and printed Anderson-Darling values of 0.25 that are scipy's upper bound. It writes a table.
- Two library modules that nothing imported are gone. `pandas` is capped below 4, which will
  refuse the positional axis the code passes to `mean()` and `std()`.

### Added
- **"Single sample" has a precise meaning in the code.** The estimators are relative to a cohort;
  a sample outside it is scored against a fixed reference. `fit_reference()` and
  `infer_from_reference()` in `python/lib/phase.py` do that, the leave-one-out validation uses
  them (its output is identical to the last bit), and `order_anchored()` refuses a lone sample
  instead of returning NaN.
- **Meta-analysis, step 11.** Intervals with the t distribution on cohorts minus one degrees of
  freedom beside the normal ones; leaving each cohort out; one estimate per cohort; bootstrap
  standard errors; expression not centred per sample; and the prediction interval for a new
  cohort. Step 10 writes the bootstrap standard error and the uncentred correlation of each
  estimate. All in `meta_three_level.tsv`.
- **Cross-tissue prediction, step 11b.** Two specificity nulls, random predictor genes and random
  outcome genes, and the variant with per-sample centring. The published rows are unchanged.
- Cross-tissue axis: the same concordances without the pair whose inferred sex differs between
  tissues.
- `tests/test_single_sample.py` and `tests/test_model_null.py`. The second applies the estimator
  to simulated cohorts in which the receptor gain varies but does not act: it finds no relation
  there, and finds it when the forcing is on.

### Not changed, and why
The figures still carry hand-written statements in panel titles, and Fig. 5c compares the model in
its own units with observed log2 differences. Both belong to decisions about the manuscript.

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
