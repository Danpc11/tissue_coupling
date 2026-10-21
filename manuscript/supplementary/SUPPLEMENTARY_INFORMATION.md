# Supplementary Information

**Glucocorticoid signalling, not insulin, tracks circadian output across human metabolic tissues**

Mejía-Arratia YE, Escutia-Guadarrama L, Pérez-Calixto M, Cervera A, Vázquez-Victorio G, Tovar H, Pérez-Calixto D

This document contains the supplementary figures with their legends, the legends of the supplementary tables, and the Supplementary Methods. The tables themselves are supplied as a separate Excel workbook, one sheet per table.

---

# Supplementary figures

## Supplementary Fig. S1 | The receptor–clock association in GTEx under technical covariates

![](ED1.png)

**A**, Spearman correlation between NR3C1 and the clock-output score in three GTEx v11 tissues,
adjusted for sex, age and Hardy scale, and additionally for ischaemic time and RIN. **B**,
Correlation of ischaemic time and RIN with the clock-output score in the same tissues. **C**,
Visceral adipose tissue (*n* = 587); dashed line, least-squares fit. Expression is centred per
sample.

---

## Supplementary Fig. S2 | Inferred clock amplitude in two muscle cohorts

![](ED2.png)

**A–C**, Inferred amplitude per person in GSE13070 (51 resistant, 18 sensitive) and GSE182120 (25
T2D, 24 NGT) under three estimators: projection onto an 800-gene template, clock genes with anchored
ordering, and clock genes with free ordering. Grey, sensitive or NGT; orange, resistant or T2D;
bars, medians; *P*, Mann–Whitney. **D**, Difference in medians for each estimator and cohort.
**E**, The two estimators compared across all samples. **F**, Ratio of the 5–95 percentile range
between groups; *P* from 3,000 permutations, with the Levene *P* (median-centred) alongside.
Coloured markers, *P* < 0.05.

---

## Supplementary Fig. S3 | Preprocessing decisions

![](ED3.png)

**A**, Receptor–clock correlation against the median correlation of a random gene with the
clock-output score, in ten cohort–tissue matrices, before and after centring each sample; lines
join the same matrix. **B**, Genes matched in MSigDB v2026.2 (H + C2 + C5) using the platform
symbol and after resolving Entrez identifiers to current symbols; numbers, genes gained. **C**,
Probes carrying a retired symbol, by platform.

---

## Supplementary Fig. S4 | Regulatory architecture of the glucocorticoid axis in adipose tissue

![](ED4.png)

**A**, Chord diagram of the five transcription factors whose inferred activity (DoRothEA A–C
regulons) correlates most strongly with the glucocorticoid signalling score, and their curated
targets with |ρ| > 0.35 against that score. Upper arcs, factors, with the correlation of their
inferred activity given outside; lower arcs, targets, coloured by the sign of their correlation;
chords, curated factor-to-target edges. **B**, Correlation of each of the 38 FOXO3 targets with
the glucocorticoid signalling score; TXNIP and FBXO32, referred to in the text, are in bold. AAGMEx adipose tissue, *n* = 232, sample-centred expression. The
ranking of factor activity against clock output is Fig. 4A and is not repeated here. Full regulon
in Supplementary Table S18.

---

## Supplementary Fig. S5 | Within-person response to the clamp and to thiazolidinedione

![](ED5.png)

**A–D**, Median per-person change of each gene-set score in the four contrasts of GSE13070, with
the Wilcoxon *P* and the percentage of people changing downwards. Coloured bars, *P* < 0.05;
*n*, pairs.

---

## Supplementary Fig. S6 | An empirically derived glucocorticoid signature

![](ED6.png)

**A**, Per-person response to hydrocortisone against placebo across 12,641 genes (GSE148640; 10
people); highlighted, the 68 genes selected for the empirical signature. **B**, Genes by the number
of participants changing in the same direction. **C**, Score difference for the empirical signature
and for the pre-specified nine-gene set in four cohorts not used to build either.

---

## Supplementary Fig. S7 | Model robustness

![](ED7.png)

**A**, Ensemble amplitude against effective forcing for four integration steps. **B**, Level and
amplitude differences between low and high coupling, for a morning window and for uniform sampling
hours; mean ± s.d. of five seeds. **C**, Ensemble synchrony against effective forcing for three
values of the dispersion of intrinsic cellular periods.

---

# Supplementary table legends

**Supplementary Table S1.** Cohorts: accession, tissue, platform, sample number and the role each series plays in the analysis.

**Supplementary Table S2.** Gene sets with their sources. Both sets were fixed before any analysis of these cohorts.

**Supplementary Table S3.** Differential expression summary per cohort: genes reaching FDR < 0.1, group sizes and statistical engine.

**Supplementary Table S4.** Set-level enrichment per cohort: z against 2,000 size-matched random sets, for the four gene sets.

**Supplementary Table S5.** Leave-one-out validation of inferred amplitude, per retained sample and per donor, in both validation cohorts.

**Supplementary Table S6.** Meta-analysis: per-cohort estimates with covariates, Fisher z, standard errors, and the pooled two-level, three-level and cluster-robust estimates.

**Supplementary Table S7.** Calibration of the hyperinsulinaemic clamp against the three cohorts with recorded sampling times, gene by gene.

**Supplementary Table S8.** Clock gene signatures across twelve contrasts, restricted to genes measured in all of them.

**Supplementary Table S9.** Transcription-factor activity against clock output: 134 DoRothEA A–C regulons with regulon size.

**Supplementary Table S10.** Human Phenotype Ontology terms in the receptor axis of adipose tissue, skeletal muscle and pancreatic islets.

**Supplementary Table S11.** Inferred amplitude and phase per person in the two muscle cohorts, under the three estimator variants.

**Supplementary Table S12.** Full enrichment results across the MSigDB collections tested against the receptor axis.

**Supplementary Table S13.** Effect of sample centring on the null distribution in ten cohort–tissue matrices.

**Supplementary Table S14.** Gene matching in MSigDB with platform symbols and after resolving Entrez identifiers.

**Supplementary Table S15.** Probes carrying a retired symbol, by platform.

**Supplementary Table S16.** Per-gene response to hydrocortisone against placebo, with the number of participants changing in each direction.

**Supplementary Table S17.** The empirically derived glucocorticoid signature and the external validation it failed.

**Supplementary Table S18.** The DoRothEA A–C FOXO3 regulon used in Fig. 4a and Supplementary Fig. S4.

**Supplementary Table S19.** Within-person score changes in the four contrasts of the clamp and thiazolidinedione cohort.

**Supplementary Table S20.** Model predictions: coupling against inferred amplitude and phase concentration, at three noise levels.

**Supplementary Table S21.** Sensitivity of the receptor–clock association to the removal of principal components, in the cohorts of the main analysis.

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