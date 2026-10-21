# Glucocorticoid signalling tracks circadian output across human metabolic tissues

Yalith Estefania Mejía-Arratia^a^, Lidia Escutia-Guadarrama^b,c^, Mitzi Pérez-Calixto^b^, Alejandra Cervera^a^, Genaro Vázquez-Victorio^b^, Hugo Tovar^a^, Daniel Pérez-Calixto^a,\*^

^a^ Instituto Nacional de Medicina Genómica (INMEGEN), Periférico Sur 4809, 14610 Mexico City, Mexico
^b^ Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM), Ciudad Universitaria, 04510 Mexico City, Mexico
^c^ Facultad de Química, Universidad Nacional Autónoma de México (UNAM), Ciudad Universitaria, 04510 Mexico City, Mexico

\* Corresponding author. Daniel Pérez-Calixto, Instituto Nacional de Medicina Genómica, Periférico Sur 4809, 14610 Mexico City, Mexico. E-mail: dperez@inmegen.gob.mx

---

## Highlights

- Per-sample clock amplitude validated out of sample against cosinor in two human cohorts
- Receptor level predicts clock output in 13 estimates from four metabolic tissues
- Hyperinsulinaemic clamp is indistinguishable from three hours of morning
- Hydrocortisone and insulin move PER1 in opposite directions in every participant
- A forced-oscillator model reproduces the signature with no fitted parameter

---

## Abstract

**Background.** The circadian clock of peripheral tissues is disrupted in obesity and type 2
diabetes, yet the signal that sets it in humans is unknown. Insulin and glucocorticoids are both
candidates.

**Methods.** We measured circadian phase and amplitude in single biopsies from 19 public
transcriptomic series and GTEx v11, covering skeletal muscle, adipose tissue and pancreatic islets,
and related them to glucocorticoid signalling across individuals and under time-controlled
perturbations.

**Results.** The peripheral clock follows glucocorticoid signalling. Across 13 cohort–tissue
estimates from four tissue types, NR3C1 expression predicts the output arm of the clock (pooled ρ = 0.45,
95% CI 0.31–0.58), and the receptor of one tissue predicts the clock of another as well as its own.
The acute insulin response is not separable from the passage of time: three hours of
hyperinsulinaemic clamp reproduce the change of three hours of morning, whereas hydrocortisone
against placebo at the same hour moves PER1 in every participant, opposite in direction to
insulin. Insulin-resistant adipose tissue carries the signature of cortisol exposure, which
thiazolidinedione does not reverse. What unifies obesity, Cushing syndrome and exogenous
hydrocortisone is not the amount of hormone but the strength of its rhythmic signal. A model in
which that strength is the product of the cortisol rhythm and receptor availability reproduces all
nine signature genes in direction, and their mean magnitude to within 4%, with no fitted parameter.

**Conclusions.** Glucocorticoid signalling tracks circadian output in human metabolic tissues, and
estimates the glucocorticoid signal a tissue receives, which no circulating measure provides.

**Keywords:** circadian clock; glucocorticoid receptor; cortisol; insulin resistance; entrainment;
adipose tissue; skeletal muscle

## 1. Introduction

The circadian clock of peripheral tissues is disrupted in obesity, insulin resistance and type 2
diabetes. Clock gene expression is altered in skeletal muscle and adipose tissue of affected
individuals, shift work and sleep
restriction raise metabolic risk, and clock gene mutations in rodents produce insulin resistance[1–5]. What remains unknown is the signal that sets the peripheral clock in humans, and therefore
what is broken when the clock is disrupted.

Two candidates have dominated the debate. Insulin entrains peripheral oscillators in cultured cells and in rodent
tissue, and the proposal that feeding entrains the periphery through insulin is widely cited[6–8]. Glucocorticoids entrain peripheral clocks through glucocorticoid response elements in the
promoters of PER1 and PER2, and the cortisol rhythm is the canonical systemic timing signal[9–12].
The two hypotheses make opposite predictions about what happens to the clock in insulin resistance,
and they have not been separated in human tissue.

The obstacle is design[13,14]. Separating hormone from time of day requires a perturbation with a
time-matched control, because any intervention that takes hours also moves the clock by those hours. Relating the clock to metabolic state requires a phenotype measured in the
same people. Distinguishing a systemic signal from a local programme requires two tissues from
the same person. No single human cohort has all three, and the cohorts that have one lack the
others. We therefore assembled twenty public datasets in which each design contributes one
element and the argument rests on their intersection. Cohorts with recorded sampling times calibrate what the passage of
time alone does to each gene; cohorts with repeated biopsies per person provide measured amplitude
against which to validate the inference; cohorts with paired tissues and large sample size provide
power; and perturbations with a time-matched control separate hormone from hour.

Here we show that the peripheral clock of human metabolic tissues tracks glucocorticoid
signalling. We develop and validate a single-sample estimator of circadian phase and amplitude, use
it across four tissue types and 2,900 individuals, and show that acute insulin is indistinguishable from
the passage of time while a physiological dose of hydrocortisone is not.
We then show that insulin-resistant adipose tissue carries the transcriptional signature of
cortisol exposure, and that a forced-oscillator model in which the effective forcing is the product
of the cortisol rhythm and receptor availability reproduces that signature without fitted
parameters. The picture that emerges is that what the peripheral clock reads is not the amount of
circulating hormone but the strength of the rhythmic signal a tissue actually receives, and that
this is what is reduced in insulin-resistant adipose tissue.

---

## 2. Methods

A full account of every procedure, with the parameter values, the sensitivity analyses and the
preprocessing decisions, is given in Supplementary Methods.

### 2.1 Cohorts and preprocessing

Nineteen public transcriptomic series and GTEx v11 were used; Supplementary Table S1 lists the
accession, tissue, platform, sample number and role of each. two accessions are the two tissues of a single study
(GSE95674 and GSE95675, AAGMEx), while GSE182117 and GSE182120 are separate studies despite
adjacent accessions. Microarray series were annotated by Entrez identifier rather than by the
symbol supplied with the platform, because platform annotations age and 8.3% of GPL570 probes with a
single Entrez identifier carry a retired symbol[15] (Supplementary Fig. S3b,c and Supplementary Tables S14 and S15). RNA-seq series
were used as raw counts. Expression was centred per sample before any score was computed: without centring, the median
correlation of a random gene with the clock-output score is 0.44 rather than zero in RNA-seq (Supplementary Fig. S3a and Supplementary Table S13).

### 2.2 Gene sets

The mammalian clock has two parts: a core transcriptional feedback loop that generates the
oscillation, and an output arm whose rhythm the loop drives and through which the clock reaches
metabolism. Clock output was quantified as the mean z-score of PER1, PER2, PER3, NR1D1, NR1D2, DBP,
TEF, HLF and CIART, the output arm of the mammalian oscillator as defined in the circadian
literature[9,11], excluding the core negative-feedback components, which were scored
separately. Receptor level was taken as the expression of NR3C1, the gene encoding the glucocorticoid
receptor. Glucocorticoid signalling was quantified as the mean z-score of FKBP5, TSC22D3, ZBTB16,
KLF15, DDIT4, SESN1, GLUL, SGK1 and CEBPD, canonical direct targets of the receptor curated from
the glucocorticoid signalling literature[9,10,12]. Both sets were fixed before any analysis
of these cohorts (Supplementary Table S2). No usable glucocorticoid-response set exists in MSigDB
for metabolic tissue[16]; a set derived empirically from the hydrocortisone perturbation
failed to validate against Cushing syndrome and dexamethasone and was not used (Supplementary Fig. S6 and Supplementary Table S17).

### 2.3 Single-sample estimation of phase and amplitude

A bulk biopsy is one point of a 24-h cycle whose sampling time is usually unrecorded, so phase and
amplitude must be recovered from the sample itself. The estimator has three steps.

*Ordering.* Expression of the clock genes is standardised across the samples of a cohort to give
*z*<sub>ig</sub> for gene *g* in sample *i*. Because the clock genes rise and fall in a fixed
sequence, samples taken at different times lie on a closed curve in that gene space, and a cohort
recovers the cycle even when no sampling time is known. Each sample is assigned an angular position
φ<sub>i</sub> on that curve, the inferred phase, which is free up to a rotation and a direction of
travel: only phase *differences* are interpretable, never absolute clock time.

*Template.* The waveform of each gene on the circle is summarised by its first-harmonic
coefficients *a*<sub>g</sub> and *b*<sub>g</sub>, so gene *g* peaks at angle
arctan(*b*<sub>g</sub>/*a*<sub>g</sub>). The template may come from the cohort itself or from an
external cohort with recorded times; an external template is what makes phase comparable between
tissues (Supplementary Methods).

*Amplitude.* The amplitude of a sample is how far it sits from the centre of the cycle along the
direction of its own phase,

*A*<sub>i</sub> = Σ<sub>g</sub> *z*<sub>ig</sub> (*a*<sub>g</sub> cos φ<sub>i</sub> + *b*<sub>g</sub> sin φ<sub>i</sub>) / √*G*   (1)

where *G* is the number of clock genes measured and the normalisation by √*G* makes amplitudes
comparable between platforms with different gene coverage. Because the projection is taken at the sample's own phase, amplitude is independent of phase by
construction.

Phase is assigned by the signed, not the absolute, projection. With a single harmonic the absolute
value cannot distinguish phase from antiphase, and on real data that rule reassigns 34 to 36% of
samples by more than 9 h (Supplementary Methods).

Two variants were used as sensitivity analyses: an *anchored* one that fixes the period against an
external template, which resolves phase but degrades amplitude, and the *free* one used throughout,
which resolves amplitude but not period (Supplementary Fig. S2 and Supplementary Table S11).

*Validation.* Measured amplitude, the reference, was obtained by cosinor regression[17] of
each clock gene on recorded sampling time with a 24-h period, taking the modulus of the first-harmonic coefficients. Validation was leave-one-out at two levels: the retained sample was
excluded from the cosinor fit of its own donor, and the whole donor from the ordering and the template
(Supplementary Table S5).

*Relation to published methods.* CHIRAL[18,19], CYCLOPS[20], ZeitZeiger[21] and
TimeSignature[22] infer the circadian time of a single sample; none returns a per-sample amplitude,
the quantity the analyses below require. Published phase inference was therefore used as a concordance control. Taking *R* = |mean<sub>i</sub> exp(i(φ<sub>muscle,i</sub> −
φ<sub>adipose,i</sub>))|, the concentration of the muscle-to-adipose phase difference across donors, our estimator and
CHIRAL both detect the coupling in the two cohorts with paired tissues
(*R* = 0.19 and 0.15 in AAGMEx, 232 donors; 0.46 and 0.29 in GTEx, 647 donors; all *P* ≤ 0.007),
with the magnitude depending on the implementation (Supplementary Methods).

### 2.4 Differential expression, gene sets and enrichment

*Differential expression.* Microarray cohorts were analysed with a moderated *t* statistic: the
variance of each gene is shrunk towards a prior estimated from all genes, which stabilises the
statistic when few samples are available, and the degrees of freedom are augmented
accordingly[23]. RNA-seq cohorts were analysed with PyDESeq2[24] on raw counts,
fitting a negative binomial model with gene-wise dispersion shrinkage. *P* values were corrected by
Benjamini–Hochberg[25] and reported as FDR. Within-person designs were analysed as the per-person
difference with a Wilcoxon signed-rank test, reporting the number of participants changing in each
direction (Supplementary Tables S3 and S19).

The insulin-resistant against insulin-sensitive contrast of GSE13070 is unbalanced (35 against 6),
so it is excluded from the differential expression of Fig. 1, where imbalance inflates the gene
count, and used only to define a profile of change compared with other profiles, where imbalance
costs precision but does not bias direction.

*Gene-set enrichment.* Set-level change was expressed as

*z* = (*m*<sub>obs</sub> − *m*<sub>null</sub>)/*s*<sub>null</sub>   (2)

where *m*<sub>obs</sub> is the mean statistic of the genes
in the set and *m*<sub>null</sub> and *s*<sub>null</sub> are the mean and standard deviation of the
same quantity across 2,000 sets of equal size drawn at random from the genes measured in that
cohort. Sampling from the measured universe rather than from a reference list keeps the null free of
platform coverage. Collections were tested against a continuous ranking of all genes by their partial correlation
with receptor level, so no significance threshold enters the enrichment (Supplementary Tables S4 and S12).

*Transcription-factor activity.* Regulons were taken from DoRothEA confidence levels A to C, which
are the interactions supported by curated literature, ChIP-seq and binding motifs[26]. The
activity of a factor in a sample is the mean z-score of its targets, each weighted by the sign of
the interaction so that repressed targets contribute with the opposite sign. Clock genes were
removed from every regulon before correlating activity with clock output, because leaving them in
makes the correlation circular (Supplementary Fig. S4 and Supplementary Tables S9 and S18).

### 2.5 Dispersion of inferred amplitude

Because per-sample amplitude is attenuated, the variance between individuals contains a measurement
component. To separate the two, amplitude was re-inferred 80 times per person from random 60%
subsets of the clock genes; the dispersion across re-inferences is the measurement noise of that
person, and the mean is their amplitude. Groups were compared on the ratio of the 5–95 percentile
range of the person-level amplitude, with *P* from 3,000 permutations, and on the noise itself
(Supplementary Methods).

### 2.6 Meta-analysis

Per-cohort Spearman correlations were transformed as *z* = arctanh(ρ) with variance 1/(*n* − 3 − *k*)   (3)

where *n* is the number of individuals and *k* the number of covariates regressed
out. Because several cohorts contribute more than one tissue from overlapping donors, the estimates
are not independent: the primary estimate is from a three-level random-effects model fitted by
restricted maximum likelihood, with estimates nested in cohort and two variance components (τ²
between cohorts, σ² between estimates within cohort). DerSimonian–Laird[27] two-level
estimates and cluster-robust variance with a small-sample correction are reported alongside it
(Supplementary Table S6 and Supplementary Methods).

### 2.7 Forced-oscillator model

Each cell was modelled as a Poincaré amplitude–phase oscillator driven by a 24-h timing signal
(zeitgeber)[28–30]:

d*r*/d*t* = λ*r*(*A* − *r*) + *FK* cos(θ<sub>z</sub> − φ)   (4)

dφ/d*t* = ω + (*FK*/*r*) sin(θ<sub>z</sub> − φ)   (5)

where *r* is the amplitude of the cell, φ its phase, ω = 2π/τ its intrinsic angular frequency with τ
drawn from a normal distribution of mean 24 h, λ the relaxation rate towards the limit cycle, *A*
the limit-cycle amplitude, θ<sub>z</sub> = 2π*t*/24 the phase of the zeitgeber, *F* the amplitude of
the cortisol rhythm and *K* the availability of the receptor. The effective forcing is the product
*F·K*. A tissue was modelled as a population of 40–150 such cells and the biopsy as the mean of the
population at one instant. Parameters were taken from the literature[29,31] and none was
fitted to the data analysed here (Supplementary Table S20 and Supplementary Methods).

### 2.8 Code and data availability

All analyses are reproducible from public data with a single script (github.com/Danpc11/
tissue_coupling, archived at Zenodo doi:10.5281/zenodo.23379059); that contrast is supplied as a
versioned input (Supplementary Table S8), and a test checks that every input a figure reads is
either produced by the pipeline or versioned. Analyses used Python 3.9.0;
stochastic routines used NumPy default_rng with seed 0.

## 3. Results

### 3.1 Group comparison detects an altered clock but not its driver

To establish what the conventional approach yields, we compared cases with controls within each of
five cohorts spanning skeletal muscle, adipose tissue and pancreatic islets (Fig. 1a). The number
of genes reaching FDR < 0.1 ranged from none in three cohorts to 210 in islets, despite all five
testing the same clinical contrast (Fig. 1b and Supplementary Table S3). None of the eighteen genes
of the clock-output and glucocorticoid sets reached that threshold in any cohort (Fig. 1c), and the
distribution of uncorrected *P* values was close to uniform in every cohort except the largest
(Fig. 1i), so the absence is not a consequence of multiple-testing correction.

At set level the picture differs. The clock-output set was lower in cases in skeletal muscle
(z = −5.9 in GSE25462) and in adipose tissue (z = −2.4 in GSE26637), whereas the glucocorticoid set
was not detected in any cohort (Fig. 1d and Supplementary Table S4). Per-person scores confirm why: in the three cohorts where that score could be computed, the two
distributions overlap almost completely and only one separates the groups (P = 0.001, 0.095 and
0.281; Fig. 1h). Group comparison therefore reports that the clock is altered without identifying
what alters it.

The signal that is missing is of a different kind. In the same cohorts, the level of the
glucocorticoid receptor and the clock-output score correlate across individuals, reaching ρ = 0.54
in GSE182120 (n = 49), a cohort with no differentially expressed gene at all (Fig. 1j). That relationship is a
correlation between people, not a difference between groups, and a case–control design cannot
detect it by construction.

### 3.2 Circadian phase and amplitude can be inferred from a single sample

Measuring that relationship requires a per-person estimate of the clock, so we inferred phase and
amplitude from single samples (Fig. 1e, Methods 2.3) and validated the inference against cosinor analysis — the standard sinusoidal fit of expression on
sampling time — of time series from the same donors. Validation was leave-one-out at two levels: the
retained sample was excluded from the cosinor fit of its own donor, and the whole donor was excluded
from the ordering and from the template. In synchronised human myotubes sampled at eight time
points, the amplitude inferred from a single retained sample agreed with the measured amplitude of
its donor at ρ = 0.28 in control and ρ = 0.34 in high-glucose conditions (P = 0.008 and 0.002, 170
retained samples); averaging the eight independent out-of-sample inferences of each donor raised the
agreement to ρ = 0.97 and 0.94 (Fig. 1f and Supplementary Table S5). In skeletal muscle *in vivo* (GSE129843), with six biopsies per person and condition,
single-sample agreement was ρ = 0.17 (P = 0.067, 117 retained samples) and donor-averaged agreement
ρ = 0.82 (P = 1.1 × 10⁻⁵, 20 series; Fig. 1g). A template built in muscle predicted the recorded sampling
time of adipose tissue under constant routine (GSE87761) to within 1.7 h (Supplementary Fig. S2 and
Supplementary Table S5).

Inferred amplitude also reproduces an externally published effect: in a controlled 10-h delay of
the habitual sleep period (GSE107537)[42], it fell in all eight participants, by a median of 29%
(P = 0.008; Fig. 1k), while phase moved by 0.6 h — the published finding for that cohort is dampened rhythms
at a phase that follows habitual bedtime rather than the shifted schedule.

The cohorts analysed below contribute one sample per person, so the single-sample figure is the one
that applies to them. Per-person amplitude is an attenuated measure, and the correlations reported
in the following sections are therefore lower bounds.

### 3.3 Clock output tracks glucocorticoid receptor level across four tissues

Having established that the clock can be measured per person (Fig. 1f,g), we asked what it follows.
Across 13 cohort–tissue estimates from six cohorts and four tissue types (skeletal muscle,
subcutaneous and visceral adipose tissue, pancreatic islets), NR3C1 expression predicted clock
output (pooled ρ = 0.45, 95% CI 0.31–0.58 by three-level
random effects; ρ = 0.43, 95% CI 0.29–0.55 by cluster-robust variance; Fig. 2a and Supplementary
Table S6). The association was positive in every adipose-tissue estimate but heterogeneous in magnitude
(ρ = 0.42, *I*² = 69%, six estimates), and strongest in pancreatic islets (ρ = 0.68, two
estimates); twelve of the thirteen confidence intervals excluded zero. Heterogeneity across
all thirteen was *I*² = 84%, of which tissue accounts for 89% by meta-regression. The relationship is monotonic rather than driven by extremes:
the median clock-output score rises across quintiles of receptor level in both islets (GSE159984)
and GTEx visceral adipose tissue (Fig. 2b).

In GTEx, where ischaemic time and RNA integrity are recorded, both covariates correlate with the
clock-output score (ρ = 0.21 to 0.29 for ischaemic time); adjusting for them changed the
receptor–clock correlation by less than 0.01 (Supplementary Fig. S1a,b).

The association is a property of the individual rather than of the tissue. In donors with several
tissues, removing each donor's mean receptor level across their tissues reduced the correlation
from 0.32 to 0.13 in subcutaneous adipose tissue and from 0.32 to 0.10 in muscle (Fig. 2d), and the
receptor measured in one tissue predicted the clock of another as well as the receptor of that
tissue itself. The downstream signalling score behaves the same way in living donors (pooled
ρ = 0.27, 95% CI 0.08–0.45, ten estimates) but is not measurable in post-mortem tissue (ρ = −0.04; Fig. 2c), where
acute receptor targets are dominated by terminal cortisol. The coupling does not require an intact hypothalamic–pituitary–adrenal (HPA) axis: in the cultured islet preparations of the same series, with no
circulating cortisol, receptor level and clock output still correlate (ρ = 0.60, P = 1.6 × 10⁻⁵,
n = 45; Fig. 2e).

### 3.4 Hydrocortisone shifts clock gene expression; the clamp reproduces the morning

Correlation across individuals cannot distinguish the two candidate signals, so we turned to
perturbations in which sampling time is controlled (Fig. 3a). A hyperinsulinaemic clamp in skeletal muscle lowered clock-output genes in 61 of 68 participants.
But the magnitude of that change is the magnitude of the morning: in three cohorts with recorded
sampling times, NR1D1, DBP, CIART, PER1, PER3 and TEF fall by the same amount over a comparable
interval without any insulin (Fig. 3b and Supplementary Table S7). The same holds in adipose tissue, where a paired clamp lowered the output arm in all ten
participants (score −1.43, P = 0.002) while leaving the core oscillator unchanged (P = 0.92;
Fig. 3c), and in the four within-person contrasts of the clamp and thiazolidinedione cohort
(Supplementary Fig. S5 and Supplementary Table S19). Neither response was graded by insulin sensitivity (P = 0.42).

Hydrocortisone is different (Supplementary Table S16). Compared with saline in the same participants (GSE148640), with both
biopsies taken at 07:00, a physiological dose raised PER1 in all ten
participants and lowered PER3 in all ten, while leaving DBP, TEF, HLF and NR1D1/2 unchanged
(Fig. 3d,e). The contrast with insulin is direct: the same gene, in the same tissue, rises in ten of
ten participants under cortisol and falls in ten of ten under insulin.

Sustained exposure produces a different and broader signature. Four days of dexamethasone raised KLF15 in all ten
participants and lowered DBP in skeletal muscle (GSE17992; Fig. 3f), and adipose tissue of people with Cushing syndrome showed lower PER3, DBP,
PER2, TEF and HLF with higher BMAL1 (GSE66446; Fig. 3g).
Adipose tissue of insulin-resistant people carries that same signature. Across the eight
clock-output genes measured in all three contrasts, its profile correlates with the hydrocortisone
profile at ρ = 0.71 (P = 0.047) and with the Cushing profile at ρ = 0.55 (P = 0.16) (Fig. 3h and
Supplementary Table S8). Three months of thiazolidinedione, which improves insulin sensitivity,
recovered a negligible fraction of it: every one of those eight genes is depleted in
insulin-resistant tissue, the median recovered fraction is 4.4%, and the slope of the
thiazolidinedione change on the insulin-resistance change is +0.004, against −1 expected for full
reversal (Fig. 3i).

### 3.5 The receptor axis indexes glucocorticoid exposure

To test whether receptor level indexes glucocorticoid exposure rather than some other property of
the tissue, we ranked all 34,020 genes by their partial correlation with receptor level in adipose tissue
(AAGMEx, n = 232), adjusting for sex, age, adiposity and macrophage content. Curated collections
were then tested against that ranking. Of 134 transcription factors with a curated regulon, FOXO3 was the one whose inferred activity best
predicted clock output (ρ = 0.39) and the NF-κB family occupied the opposite extreme (NFKB1
ρ = −0.44; both n = 232; Fig. 4a and Supplementary Table S9); the FOXO3 regulon includes TXNIP and FBXO32 (Supplementary Fig. S4). The
same ordering appears in the Hallmark collections of coordinately regulated processes, where
interferon-γ response, TNFα signalling via
NF-κB and the complement run against the axis, together with oxidative phosphorylation,
adipogenesis and glycolysis (Fig. 4b, and Reactome and Gene Ontology in Fig. 4f,g).

Against phenotype ontology terms, the axis runs with truncal obesity and increased urinary cortisol
and against type 2 diabetes in adipose tissue and muscle, and in the opposite direction in islets
(Fig. 4c and Supplementary Table S10). The axis therefore tracks the clinical signature of glucocorticoid exposure rather than that of
metabolic disease. The next two results test that distinction. In the one cohort with insulin sensitivity measured in the same people[32], the clock-output score did
not correlate with it (ρ = 0.01, 95% CI −0.28 to 0.30, *n* = 48 with an insulin sensitivity index
from the frequently sampled intravenous glucose tolerance test), although it was graded by clinical status (median score
+0.31 with no family history of diabetes, +0.09 and +0.11 with one and two affected parents, and
−0.31 in diabetes). Within adipose tissue, adiposity separates the two arms of the axis: HSD11B1, the 11β-hydroxysteroid dehydrogenase type 1 that regenerates cortisol locally, rises
with leptin (ρ = 0.66) while the receptor falls (ρ = −0.31; n = 232), and
the receptor, not the enzyme, is what tracks clock output (Fig. 4d,e).

### 3.6 A forced-oscillator model reproduces the observed signature

The three clinical situations that produce this signature — obesity, Cushing syndrome and exogenous
hydrocortisone — differ in whether cortisol is low, high or imposed, which argues against the
amount of hormone as the common variable. We therefore modelled the peripheral clock as a
population of oscillators driven by a 24-h timing signal, or zeitgeber, whose effective strength is the product of the
amplitude of the cortisol rhythm and the availability of the receptor (Methods 2.7). Below a
threshold of effective forcing the population loses synchrony and the ensemble amplitude collapses
while individual cells continue to oscillate (Fig. 5a,b), which is what a bulk biopsy reads as loss
of amplitude.

Lowering the coupling reproduced the observed signature without any fitted parameter: the direction
of all nine signature genes agreed with the change measured in insulin-resistant adipose tissue and
the mean magnitude was within 4% (Fig. 5c). The model also predicts that reduced coupling always
lowers amplitude but lowers the mean level only when biopsies are taken in a narrow window: with
uniform sampling hours the predicted level difference is zero within error (+0.02 ± 0.03) while the
amplitude difference is unchanged (Fig. 5d and Supplementary Fig. S7b). This accounts for the contrast of Fig. 2c, where the glucocorticoid signalling score is measurable in
cohorts of living donors, biopsied in the morning, and not in post-mortem tissue. Applying our own estimator to model-generated cohorts gave correlations between
coupling and inferred amplitude of 0.44 to 0.63 depending on measurement noise, against 0.09 to
0.22 observed with receptor mRNA (Fig. 5e), consistent with the receptor transcript being an
attenuated proxy of the true coupling.

The model makes one further prediction that these cohorts can test: if reduced forcing acts by
releasing cells from entrainment, its effect on a cross-sectional cohort should be an increase in
the variance of amplitude between individuals.

### 3.7 Testing the dispersion prediction in insulin-resistant muscle

The model predicts that reduced forcing appears in a cross-sectional cohort as increased variance
between individuals rather than as a shift in the median, so that a test of location cannot detect
it. We therefore compared the dispersion of inferred amplitude between groups in the two muscle
cohorts. In GSE13070 skeletal muscle, the median did not differ (P = 0.21) but the 5–95
percentile range was 1.66 times wider in insulin-resistant participants (P = 0.021 by permutation;
Levene P = 0.032), with the same direction under the anchored estimator (ratio 1.90, P = 0.0036)
and no difference under the unvalidated 800-gene projection (Supplementary Fig. S2a–f and
Supplementary Table S11). Because per-sample amplitude is attenuated, we measured the noise of the
estimator in each person by re-inference from random gene subsets: it does not differ between
groups (median 0.149 in resistant against 0.160 in sensitive participants, P = 0.99), and the
variance ratio is 2.7 with both estimation schemes, so the difference is not generated by unequal
measurement error. The variance of the sensitive group is nevertheless smaller than its own squared
noise, so what can be claimed is that resistant participants disperse more, not that sensitive
participants have a particular dispersion. The effect was not present in the second muscle cohort,
which contrasts type 2 diabetes with normal glucose tolerance in balanced groups rather than
clamp-measured resistance in unbalanced ones.

## 4. Discussion

The peripheral clock of human metabolic tissues tracks glucocorticoid signalling.
The association holds across four tissues and six cohorts, survives adjustment for adiposity,
cellular composition, ischaemic time and RNA integrity, and persists in cultured islets with no
circulating hormone. The decisive comparison is within individuals and at a fixed hour of the day:
PER1 rises in every participant under hydrocortisone and falls in every participant under insulin,
in the same tissue. Three hours of clamp are indistinguishable from three hours of morning, so what
the clamp literature has read as an effect of insulin on the clock is the passage of time.

The mechanism that unifies the clinical observations is not the amount of hormone but the strength
of the rhythmic signal. Obesity flattens the cortisol rhythm[33,34], Cushing syndrome and exogenous
hydrocortisone impose sustained exposure[35,36], and adiposity lowers the receptor; all three lower the
same product, and all three produce the same transcriptional signature. The model makes that
statement quantitative, reproducing the direction of every signature gene and its magnitude to
within 4% without fitted parameters, and it yields a testable prediction that the data already
satisfy: the level signature requires a narrow sampling window while the amplitude signature does
not (Fig. 5f).

One consequence follows for metabolic disease: the signature is not a marker of insulin
resistance. Thiazolidinedione improves sensitivity without reversing it. The clock-output score does not
correlate with measured insulin sensitivity in the one cohort where that phenotype is available,
although it is graded by clinical status. And in phenotype ontology the receptor axis runs against
type 2 diabetes while running with truncal obesity and urinary cortisol. What the
signature estimates is the glucocorticoid signal a tissue actually receives, which no circulating
measure provides.

This reframes the modest results of 11β-hydroxysteroid dehydrogenase type 1 inhibitors in type 2
diabetes[37,38]. These compounds lower the local regeneration of cortisol, that is, the
ligand[39,40,41]. Our data separate the two arms of that axis: within adipose tissue HSD11B1
rises with adiposity while the receptor falls, and it is the receptor, not the enzyme, that tracks
clock output (Fig. 4d,e). Where the limiting term is receptor availability, lowering the ligand
cannot restore the rhythmic signal, which would account for the modest glycaemic effect of these
compounds. The prediction is specific and testable in material that already exists: in trials of
11β-HSD1 inhibition with banked adipose tissue[43,44], the response should be graded by the
effective glucocorticoid signal of the tissue at baseline, estimated from receptor level together
with the clock-output score, rather than by diagnostic category or by circulating cortisol. A null
association between baseline receptor level and response would argue against the interpretation
advanced here.

Our analysis has five limitations. NR3C1 expression is a proxy for receptor availability, and
chronic glucocorticoid exposure also lowers that transcript, so these data do not separate a
receptor that limits signalling from one suppressed by the exposure it reports. Phase and amplitude are inferred from expression rather than
from longitudinal sampling in each donor, and per-sample amplitude is attenuated (ρ = 0.17 to 0.34
against cosinor out of sample), so the associations reported here are lower bounds. The
cross-sectional association in adipose tissue weakens when the leading principal components of the
matrix are removed (Supplementary Table S21); that adjustment controls for between-sample
heterogeneity but also removes any systemic signal, and these data cannot separate the two
(Supplementary Methods). Cortisol was not measured in any of the cohorts. And the direction of the
relationship rests on cross-sectional association combined with external perturbations rather than
on manipulation of the receptor.
Manipulating the receptor itself in a metabolic tissue would settle the last of these, and our
results specify what such an experiment should find: an attenuation of the output arm with the core
oscillator preserved, which is the dissociation the paired clamp already shows within individuals,
rather than a loss of the clock as a whole.

Taken together, these results change what should be measured when circadian disruption is suspected
in metabolic disease. Circulating cortisol, the measure in routine use, is not the quantity the
tissue responds to; the receptor and the rhythm it reads are. A transcriptional estimate of that
signal can be obtained from a single biopsy, in cohorts that already exist, and is what should be
related to metabolic phenotype and to treatment response.

---

## Acknowledgements

We thank the participants and investigators of the original studies for making their datasets
available, and acknowledge the GTEx Consortium and the NCBI Gene Expression Omnibus for providing
access to the data used in this work.

## References

1\. Stenvers DJ, Scheer FAJL, Schrauwen P, la Fleur SE, Kalsbeek A. Circadian clocks and insulin resistance. Nat Rev Endocrinol 2019;15:75–89.
2\. Bass J. Circadian topology of metabolism. Nature 2012;491:348–56.
3\. Wefers J, van Moorsel D, Hansen J, et al. Circadian misalignment induces fatty acid metabolism gene profiles and compromises insulin sensitivity in human skeletal muscle. Proc Natl Acad Sci U S A 2018;115:7789–94.
4\. Otway DT, Mäntele S, Bretschneider S, et al. Rhythmic diurnal gene expression in human adipose tissue from individuals who are lean, overweight, and type 2 diabetic. Diabetes 2011;60:1577–81.
5\. Perrin L, Loizides-Mangold U, Chanon S, et al. Transcriptomic analyses reveal rhythmic and CLOCK-driven pathways in human skeletal muscle. eLife 2018;7:e34114.
6\. Tuvia N, Pivovarova-Ramich O, Murahovschi V, et al. Insulin directly regulates the circadian clock in adipose tissue. Diabetes 2021;70:1985–99.
7\. Crosby P, Hamnett R, Putker M, et al. Insulin/IGF-1 drives PERIOD synthesis to entrain circadian rhythms with feeding time. Cell 2019;177:896–909.
8\. Chaix A, Manoogian ENC, Melkani GC, Panda S. Time-restricted eating to prevent and manage chronic metabolic diseases. Annu Rev Nutr 2019;39:291–315.
9\. Balsalobre A, Brown SA, Marcacci L, et al. Resetting of circadian time in peripheral tissues by glucocorticoid signaling. Science 2000;289:2344–7.
10\. Dickmeis T. Glucocorticoids and the circadian clock. J Endocrinol 2009;200:3–22.
11\. Oster H, Challet E, Ott V, et al. The functional and clinical significance of the 24-hour rhythm of circulating glucocorticoids. Endocr Rev 2017;38:3–45.
12\. So AY, Bernal TU, Pillsbury ML, Yamamoto KR, Feldman BJ. Glucocorticoid regulation of the circadian clock modulates glucose homeostasis. Proc Natl Acad Sci U S A 2009;106:17582–7.
13\. Hughes ME, Abruzzi KC, Allada R, et al. Guidelines for genome-scale analysis of biological rhythms. J Biol Rhythms 2017;32:380–93.
14\. Ruben MD, Francey LJ, Guo Y, et al. A large-scale study reveals 24-h operational rhythms in hospital treatment. Proc Natl Acad Sci U S A 2019;116:20953–8.
15\. Dai M, Wang P, Boyd AD, Kostov G, Athey B, Jones EG, et al. Evolving gene/transcript definitions significantly alter the interpretation of GeneChip data. Nucleic Acids Res 2005;33:e175. https://doi.org/10.1093/nar/gni179.
16\. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database hallmark gene set collection. Cell Syst 2015;1:417–25.
17\. Cornelissen G. Cosinor-based rhythmometry. Theor Biol Med Model 2014;11:16.
18\. Talamanca L, Gobet C, Naef F. Sex-dimorphic and age-dependent organization of 24-hour gene expression rhythms in humans. Science 2023;379:478–83.
19\. Salati A, Gobet C, Naef F. scRITMO: inference of circadian phase from single-cell and bulk transcriptomes [preprint]. bioRxiv 2026. Available at github.com/AndreaSalati/scRITMO
20\. Anafi RC, Francey LJ, Hogenesch JB, Kim J. CYCLOPS reveals human transcriptional rhythms in health and disease. Proc Natl Acad Sci U S A 2017;114:5312–7.
21\. Hughey JJ, Hastie T, Butte AJ. ZeitZeiger: supervised learning for high-dimensional data from an oscillatory system. Nucleic Acids Res 2016;44:e80.
22\. Braun R, Kath WL, Iwanaszko M, et al. Universal method for robust detection of circadian state from gene expression. Proc Natl Acad Sci U S A 2018;115:E9247–56.
23\. Smyth GK. Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Stat Appl Genet Mol Biol 2004;3:Article3.
24\. Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics 2023;39:btad547.
25\. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B 1995;57:289–300.
26\. Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Res 2019;29:1363–75.
27\. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials 1986;7:177–88.
28\. Granada AE, Bordyugov G, Kramer A, Herzel H. Human chronotypes from a theoretical perspective. PLoS One 2013;8:e59464.
29\. Abraham U, Granada AE, Westermark PO, Heine M, Kramer A, Herzel H. Coupling governs entrainment range of circadian clocks. Mol Syst Biol 2010;6:438.
30\. Schmal C, Herzel H, Myung J. Clocks in the wild: entrainment to natural light. Front Physiol 2020;11:272.
31\. Li Y, Androulakis IP. The SCN–HPA–periphery circadian timing system: mathematical modeling of clock synchronization and the effects of photoperiod on jetlag adaptation. J Biol Rhythms 2023;38:435–53.
32\. Sharma NK, Sajuthi SP, Chou JW, et al. Tissue-specific and genetic regulation of insulin sensitivity-associated transcripts in African Americans. J Clin Endocrinol Metab 2016;101:1455–68.
33\. Björntorp P, Rosmond R. Obesity and cortisol. Nutrition 2000;16:924–36.
34\. Incollingo Rodriguez AC, Epel ES, White ML, Standen EC, Seckl JR, Tomiyama AJ. Hypothalamic–pituitary–adrenal axis dysregulation and cortisol activity in obesity: a systematic review. Psychoneuroendocrinology 2015;62:301–18.
35\. Pivonello R, Isidori AM, De Martino MC, Newell-Price J, Biller BMK, Colao A. Complications of Cushing's syndrome: state of the art. Lancet Diabetes Endocrinol 2016;4:611–29.
36\. Hackett RA, Kivimäki M, Kumari M, Steptoe A. Diurnal cortisol patterns, future diabetes, and impaired glucose metabolism in the Whitehall II cohort study. J Clin Endocrinol Metab 2016;101:619–25.
37\. Feig PU, Shah S, Hermanowski-Vosatka A, et al. Effects of an 11β-hydroxysteroid dehydrogenase type 1 inhibitor, MK-0916, in patients with type 2 diabetes mellitus and metabolic syndrome. Diabetes Obes Metab 2011;13:498–504.
38\. Rosenstock J, Banarer S, Fonseca VA, et al. The 11-beta-hydroxysteroid dehydrogenase type 1 inhibitor INCB13739 improves hyperglycemia in patients with type 2 diabetes inadequately controlled by metformin monotherapy. Diabetes Care 2010;33:1516–22.
39\. Morgan SA, McCabe EL, Gathercole LL, et al. 11β-HSD1 is the major regulator of the tissue-specific effects of circulating glucocorticoid excess. Proc Natl Acad Sci U S A 2014;111:E2482–91.
40\. Masuzaki H, Paterson J, Shinyama H, et al. A transgenic model of visceral obesity and the metabolic syndrome. Science 2001;294:2166–70.
41\. Lee RA, Harris CA, Wang JC. Glucocorticoid receptor and adipocyte biology. Nucl Receptor Res 2018;5:101373.
42\. Kervezee L, Cuesta M, Cermakian N, Boivin DB. Simulated night shift work induces circadian misalignment of the human peripheral blood mononuclear cell transcriptome. Proc Natl Acad Sci U S A 2018;115:5540-5. https://doi.org/10.1073/pnas.1720719115.
43\. Othonos N, Pofi R, Arvaniti A, et al. 11β-HSD1 inhibition in men mitigates prednisolone-induced adverse effects in a proof-of-concept randomised double-blind placebo-controlled trial. Nat Commun 2023;14:1025.
44\. Ajjan RA, Hensor EMA, Del Galdo F, et al. Oral 11β-HSD1 inhibitor AZD4017 improves wound healing and skin integrity in adults with type 2 diabetes mellitus: a pilot randomized controlled trial. Eur J Endocrinol 2022;186:441–55.

## Declarations

**Author contributions (CRediT).** Y.E.M.-A., H.T., A.C. and D.P.-C.: formal analysis,
investigation. L.E.-G., M.P.-C. and G.V.-V.: data curation. H.T., L.E.-G., A.C., G.V.-V. and
D.P.-C.: methodology. H.T. and D.P.-C.: supervision. G.V.-V., L.E.-G. and D.P.-C.: funding
acquisition. Y.E.M.-A., L.E.-G., H.T. and D.P.-C.: writing – original draft. All authors: writing –
review and editing, and approval of the final version. D.P.-C. is the guarantor of this work.

**Funding.** Supported by SECIHTI projects CBF-2025-G-789 and PEE-2026-J-26, and by
DGAPA-PAPIIT IA209424. Y.E.M.-A. acknowledges a postgraduate scholarship from SECIHTI
(CVU 2057007).

**Conflicts of interest.** The authors declare no competing interests.

**Ethics.** All data are publicly available and were generated under the approvals reported in the
original publications; no new human data were collected.

**Data availability.** Gene-expression datasets are available from the NCBI Gene Expression Omnibus
under accessions GSE11868, GSE13070, GSE17992, GSE25462, GSE26637, GSE41168, GSE66446,
GSE78721, GSE87761, GSE95674, GSE95675, GSE104674, GSE107537, GSE108539, GSE129843,
GSE148640, GSE159984, GSE182117 and GSE182120. GTEx v11 expression data and sample
annotations are available through the GTEx Portal, subject to its access conditions; individual-level
controlled-access data require the applicable dbGaP authorisation. Dataset identifiers and cohort
details are listed in Supplementary Table S1.

**Code availability.** Analysis code is available at https://github.com/Danpc11/tissue_coupling and
archived on Zenodo (https://doi.org/10.5281/zenodo.23379059). Analyses used Python 3.9.0;
stochastic routines used NumPy default_rng with seed 0.

**Use of artificial intelligence tools.** AI-assisted tools were used for language editing and to
support code development for the analyses reported here. All data visualisations are derived from
the underlying data by the reproducible code deposited with this work. The authors reviewed all
resulting text, code and figures and take full responsibility for the content of the manuscript.

## Abbreviations

CI, confidence interval; FDR, false discovery rate; HPA, hypothalamic–pituitary–adrenal.

---

## Figures

![](figs/F1.png)

## Figure 1 | Single-sample inference of circadian phase and amplitude in human metabolic tissues

**A**, Study design. **B**, Genes reaching FDR < 0.1 in a case–control comparison within each
cohort (moderated *t*, Benjamini–Hochberg; log scale; group sizes in parentheses). Colour, tissue.
**C**, Moderated *t* for the nine clock-output and nine glucocorticoid genes in each cohort. Dot
size, −log10 *P* (uncorrected); colour, moderated *t*. **D**, Set-level enrichment of four gene sets
per cohort (*z* against 2,000 size-matched random sets). **E**, Single-sample estimator of phase and
amplitude. **F**, Inferred against cosinor-measured amplitude in synchronised human myotubes
(GSE182117; 11 donors, two conditions, eight time points). **G**, The same in skeletal muscle
*in vivo* (GSE129843; 11 men, six biopsies per person and condition; EXF, extended feeding; TRF,
time-restricted feeding). **H**, Clock-output score per person, control (squares) against case
(circles); bars, medians; *P*, Mann–Whitney. **I**, Distribution of uncorrected *P* across all genes
in three cohorts. **J**, Spearman correlation between NR3C1 and the clock-output score across
individuals of the same cohorts; bars, 95% CI; marker size, *n*. **K**, Cosinor amplitude of the clock-output genes in each
participant before and after a controlled 10-h delay of the habitual sleep period (GSE107537); one
line per participant. Dashed lines in **F**, **G** are least-squares fits. Expression is centred per sample throughout.

![](figs/F2.png)

## Figure 2 | The clock-output score tracks glucocorticoid receptor level across tissues

**A**, Random-effects meta-analysis of the Spearman correlation between NR3C1 and the clock-output
score, across 13 cohort–tissue estimates grouped by tissue. Circles, individual estimates with 95%
CI, size proportional to *n*; diamonds, estimates pooled by tissue and overall; dotted line, the
overall estimate. Pooled values shown are two-level DerSimonian–Laird on Fisher *z*; the primary
estimate reported in the text is from the three-level model that nests estimates within cohort
(ρ = 0.45, 95% CI 0.31–0.58). **B**, Clock-output score by NR3C1 quintile in pancreatic islets
(GSE159984) and GTEx visceral adipose tissue; each dot one person, bars medians. **C**, The same
meta-analysis for the glucocorticoid signalling score, stratified by living donors (purple) and
post-mortem tissue (grey). **D**, Receptor–clock correlation before and after regressing out each
donor's mean NR3C1 across their tissues, in the five cohort–tissue pairs with several tissues per
person. **E**, Cultured islet preparations (GSE159984; 45 samples); dashed line, least-squares fit.

![](figs/F3.png)

## Figure 3 | Hydrocortisone shifts clock gene expression; the clamp reproduces the morning

**A**, Sampling scheme of the two perturbations. **B**, Change in clock gene expression after three
hours of hyperinsulinaemic clamp (GSE13070 muscle) compared with the change over a morning interval
in three cohorts with recorded sampling times. **C**, Paired fasting-to-hyperinsulinaemia change in
adipose tissue (GSE26637; 10 people); each dot one person, bars medians. **D**, Per-person change in
PER1 under hydrocortisone against placebo (GSE148640) and under clamp against fasting (GSE26637).
**E**, Hydrocortisone minus placebo in adipose tissue, both biopsies at 07:00 (GSE148640; 10
people); numbers above, participants changing upwards. **F**, Four days of dexamethasone in skeletal
muscle (GSE17992; 10 men; two-colour arrays, log ratios). **G**, Cushing syndrome adipose tissue
(GSE66446), centred on the control median. **H**, The clock-output genes in insulin-resistant adipose
tissue, Cushing syndrome, hydrocortisone and three months of thiazolidinedione; of the nine genes
of the set, eight were measured in all four contrasts. **I**, The same eight genes in
insulin-resistant adipose tissue against their change after three months of thiazolidinedione
(GSE13070); dashed line, least-squares fit, whose slope is given above the panel.

![](figs/F4.png)

## Figure 4 | Transcriptional context of the receptor axis in adipose tissue

The receptor axis is the genome-wide ranking of genes by their partial Spearman correlation with
NR3C1 in AAGMEx adipose tissue (*n* = 232), adjusted for sex, age, LEP and CD68.
**A**, Transcription-factor activity inferred from DoRothEA A–C regulons, correlated with the
clock-output score across individuals; the eight highest and eight lowest of 134 factors; marker
size, targets in the regulon. Clock genes are removed from every regulon. **B**, Hallmark
collections in the receptor axis (*z* against size-matched random sets; dot size, genes in the set).
**C**, Human Phenotype Ontology terms in the receptor axis of each tissue (AAGMEx adipose and
muscle, GSE159984 islets); the islet ranking is not adjusted for donor covariates. **D**, Spearman correlation of each gene with LEP and with the
clock-output score. **E**, Distribution of the receptor axis across 34,020 genes, with three genes
marked. **F**, **G**, Reactome and GO biological process in the receptor axis.

![](figs/F5.png)

## Figure 5 | A forced-oscillator model of the glucocorticoid signal

Poincaré amplitude–phase oscillators forced by a 24-h zeitgeber; effective forcing F·K is the
product of the amplitude of the cortisol rhythm and receptor availability. Populations of 40–150
cells with dispersed intrinsic periods; no parameter is fitted to the data.
**A**, Ensemble synchrony *R* as a function of effective forcing and intrinsic tissue period. White
line, *R* = 0.6; markers, the two regimes referred to in the text. **B**, Ensemble amplitude (blue,
left axis) and synchrony (orange, right axis) against effective forcing. **C**, Change in clock gene
expression predicted when coupling falls from 1.4 to 0.5, against the change observed in
insulin-resistant adipose tissue (GSE13070). **D**, Difference in level and in amplitude between low
and high coupling, for three widths of the biopsy sampling window. **E**, Spearman correlation
between true coupling and the amplitude inferred by our estimator from model-generated cohorts, at
three noise levels; shaded band, the range observed with NR3C1. **F**, Summary of the model.
