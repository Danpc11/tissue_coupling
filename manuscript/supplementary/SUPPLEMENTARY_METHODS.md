# Supplementary Methods

**Glucocorticoid signalling, not insulin, tracks circadian output across human metabolic tissues**

This section gives the full account of every procedure summarised in Methods, with the parameter
values, the preprocessing decisions and the sensitivity analyses. Section numbers match those of
the main Methods.

---

## S2.1 Cohorts and preprocessing

**Annotation.** Microarray platforms were annotated by resolving each probe to its Entrez
identifier and that identifier to the current gene symbol (NCBI GRCh38.p13), rather than using the
symbol supplied with the platform file. The platform annotation of GPL570 is dated October 2014; of
its 41,177 probes with a single Entrez identifier, 3,413 (8.3%) carry a symbol that has since been
retired (ARNTL to BMAL1, ATP5A1 to ATP5F1A, among others). Matching gene sets against the platform
symbol loses those genes silently. Where several probes measure the same gene, the probe with the
highest mean expression among those measuring it exclusively was used; where no probe measures it
exclusively, the best multi-gene probe containing it was used. Two genes of the positive-control
sets, GSTT1 and HLA-DRB4, have no current Entrez identifier because they are absent from the
primary assembly, and were retained by symbol.

**Sample centring.** Every expression matrix was centred per sample, subtracting from each column
its mean across all measured genes, before any score was computed. Without centring, in RNA-seq the
median Spearman correlation of a randomly chosen gene with the clock-output score is 0.44 rather
than zero in pancreatic islets and 0.12 to 0.26 in the three GTEx tissues, because the number of
detected genes per sample correlates 0.85 with the score. Any subsequent correlation is then
inflated by that common component. In quantile-normalised microarray data the effect is negligible
(null median 0.002 in AAGMEx). The effect of centring on each cohort–tissue matrix is in
Supplementary Table S13.

**GTEx covariates.** For GTEx v11, residuals were taken against sex, age (mid-point of the
reported bracket), Hardy scale, ischaemic time (SMTSISCH) and RNA integrity (SMRIN), with
median imputation of the few missing values. Both technical covariates correlate with the
clock-output score (ρ = 0.21 to 0.29 for ischaemic time; −0.19 for RIN in visceral adipose tissue),
which is why the adjustment is reported.

## S2.3 Single-sample estimation of phase and amplitude

**Ordering.** Samples are placed on the unit circle from the standardised clock-gene matrix. Two
variants were used. The *free* variant, used throughout unless stated, fits the waveform to the
data of the cohort itself, so the projection basis is that tissue's own and amplitude is resolved,
but it does not fix the period: if the dominant cycle is the 12-h harmonic in one tissue and the
24-h fundamental in the other, the two phases are on different scales. The *anchored* variant fixes
the period by projecting onto a template of relative phases obtained by cosinor in cohorts with
recorded sampling times, which resolves phase but degrades amplitude. Measured stability ranges
across template sizes of 13 to 200 genes: phase 0.32 to 0.47 anchored and 0.21 to 0.62 free;
amplitude 0.30 to 0.37 free and 0.07 to 0.38 anchored. Using one variant for both quantities gives
one of them wrong. The default template is the 60 most rhythmic genes; with 13 genes a residual
period ambiguity remains (double-angle statistic 0.23 against 0.06 with 60).

**Signed projection.** Phase is assigned by the signed projection onto the circadian basis, not by
the absolute value. With a single harmonic, phase and antiphase give the same absolute projection
and the tie is broken by noise: on real data the absolute rule reassigns 34 to 36% of samples by
more than 9 h.

**Validation.** Validation was leave-one-out at two levels. For each retained sample, (i) the
cosinor fit of its donor was computed from that donor's remaining time points only, and (ii) the
donor as a whole was excluded from the ordering and from the template used to infer the retained
sample. In synchronised myotubes (GSE182117, 11 donors × 2 conditions × 8 time points, 170 retained
samples) the agreement is ρ = 0.28 in control and 0.34 in high-glucose conditions; averaging the
eight independent out-of-sample inferences of each donor raises it to 0.97 and 0.94. In skeletal
muscle in vivo (GSE129843, 117 retained samples) the figures are 0.17 and 0.82. The cohorts
analysed in Results contribute one sample per person, so the single-sample figure is the one that
applies to them.

**Relation to published methods.** CHIRAL, CYCLOPS, ZeitZeiger and TimeSignature all infer
circadian phase from a single sample; none returns a per-sample amplitude. We therefore used
published phase inference as a concordance control rather than as a replacement. Computing the
concentration of the muscle-to-adipose phase difference across donors,
R = |mean_i exp(i(φ_muscle,i − φ_adipose,i))|, gives R = 0.19 with our estimator and 0.15 with
CHIRAL in AAGMEx (232 donors), and 0.46 and 0.29 in GTEx (647 donors), all with P ≤ 0.007 by
permutation of the donor correspondence. The coupling is detected in both cohorts with both
estimators; its magnitude depends on the implementation. CHIRAL's expectation–maximisation
converges to one of two optima depending on initialisation, so the NumPy global seed was fixed
before each call; the values above are for seed 0. A benchmark of phase error against recorded
sampling times, for completeness: supervised methods reach 1.0 to 1.2 h and unsupervised methods,
including ours, 2.1 to 2.9 h. Our estimator is not proposed as a better phase estimator.

## S2.4 Differential expression

**Moderated t.** For microarray cohorts, the prior variance was estimated by moments on log *s*²,
each gene variance was shrunk towards it, and the degrees of freedom were augmented accordingly. In
a null simulation (4,000 genes, 25 samples) the procedure holds the type-I error at 5.1% against
4.9% for the classical t and returns no gene at FDR 0.1; with 200 genes shifted by one log2 unit it
detects 91 with 10 false positives, the same as the classical t at this sample size. Its advantage
is in the designs with 5 to 10 samples per group.

**RNA-seq.** Raw counts, genes with more than 5 counts in at least 10 samples, PyDESeq2 with the
cohort covariates in the design and shrunken log fold changes.

**Within-person designs.** The per-person difference, a Wilcoxon signed-rank test, and the number
of participants changing in each direction, which is what the figures report.

## S2.5 Meta-analysis

Spearman correlations were transformed as *z* = arctanh(ρ) with variance 1/(*n* − 3 − *k*). The
primary estimate is a three-level random-effects model fitted by restricted maximum likelihood,
with estimates nested in cohort: V = diag(se² + σ²) + τ²·J, where J is the block-diagonal indicator
of shared cohort, τ² the between-cohort variance and σ² the within-cohort variance. For the
receptor–clock association this gives ρ = 0.45 (95% CI 0.31 to 0.58), τ² = 0.036 and σ² = 0.011,
against ρ = 0.43 (0.34 to 0.51) by two-level DerSimonian–Laird and ρ = 0.43 (0.29 to 0.55) by
cluster-robust variance with a small-sample correction. The two-level interval is 55% narrower than
the three-level one, which is the expected consequence of treating as independent the three GTEx
tissues from overlapping donors.

## S2.6 Forced-oscillator model

Fixed parameters: relaxation rate λ = 0.05 h⁻¹, limit-cycle amplitude *A* = 1.0, dispersion of
intrinsic periods 1.2 to 1.5 h, 40 to 150 cells per tissue, 15 to 30 days of integration, explicit
Euler with a step of 0.05 h. Amplitude and synchrony are identical to three decimal places for
steps of 0.20, 0.10, 0.05 and 0.02 h. None of these values was fitted to the data analysed here.
Synthetic cohorts draw one coupling value *K* per subject and one sampling hour, integrate the
population, and read the mean of the population at that hour; expression is then generated for each
clock gene as a sinusoid of the ensemble amplitude at its own template phase, plus Gaussian noise.
The level difference between low and high coupling is −0.447 ± 0.028 with a morning window and
+0.020 ± 0.025 with uniform sampling hours, across five seeds, while the amplitude difference is
−0.657 ± 0.024 and −0.658 ± 0.023.

## S2.6b Sensitivity to the removal of principal components

Removing the leading principal components of an expression matrix is a standard control against
between-sample heterogeneity, but it is not neutral here: a systemic signal that moves many genes
together lives in those components by construction, so the adjustment removes the effect it is
meant to test along with any confounder. We report it as a sensitivity analysis. Removing three
components leaves the association unchanged in skeletal muscle (ρ = 0.270 to 0.264 in GTEx,
percentile 98 of the within-tissue null in both cases) and attenuates it in adipose tissue (0.302 to
0.133 in GTEx visceral). The attenuation is therefore compatible both with a compositional component
in adipose tissue and with over-correction of a systemic signal, and these data do not separate the
two. The perturbation results are not affected, because they are within-person comparisons in which
the composition of the donor is the same before and after.

## S2.7 Dispersion of inferred amplitude

Because per-sample amplitude is attenuated, the variance between individuals contains a measurement
component. In the cohort with six biopsies per person (GSE129843) the standard deviation of the
estimator within a single series is 0.481 against 0.264 between series, so only 38% of the
between-person variance is not measurement error. To test whether that error differs between
groups, amplitude was re-inferred 80 times per person from random 60% subsets of the clock genes;
the dispersion across re-inferences is the measurement noise of that person. In GSE13070 skeletal
muscle the noise does not differ between groups (median 0.149 in insulin-resistant against 0.160 in
sensitive participants; Mann–Whitney P = 0.99, Levene P = 0.31), while the variance of the
person-level amplitude differs by a factor of 2.7 and the 5–95 percentile range by 1.70
(P = 0.034 by 3,000 permutations). Constant measurement noise moves the variance ratio away from
one rather than towards it, so it does not generate the difference. The variance of the sensitive
group is nonetheless smaller than its own squared noise, so what can be claimed is that resistant
participants disperse more, not that sensitive participants have a particular dispersion.

## S2.8 Software

Python 3.9.0; NumPy, pandas, SciPy, statsmodels, PyDESeq2 0.5.4. Stochastic routines use NumPy
`default_rng` with seed 0 unless stated. All results are regenerated from public data by
`run_all.sh` in the accompanying repository (github.com/Danpc11/tissue_coupling, archived at
https://doi.org/10.5281/zenodo.23379059).
