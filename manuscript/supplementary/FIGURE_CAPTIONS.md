# Figure captions

Panels carry only their letter; the description lives here.

---

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
individuals of the same cohorts; bars, 95% CI; marker size, *n*. Dashed lines in **F**, **G** are least-squares fits. Expression is centred per sample throughout.

## Figure 2 | The clock-output score tracks glucocorticoid receptor level across tissues

**A**, Random-effects meta-analysis (DerSimonian–Laird on Fisher *z*) of the Spearman correlation
between NR3C1 and the clock-output score, across 13 cohort–tissue estimates grouped by tissue.
Circles, individual estimates with 95% CI, size proportional to *n*; diamonds, pooled estimates;
dotted line, overall estimate. **B**, Clock-output score by NR3C1 quintile in pancreatic islets
(GSE159984) and GTEx visceral adipose tissue; each dot one person, bars medians. **C**, The same
meta-analysis for the glucocorticoid signalling score, stratified by living donors (purple) and
post-mortem tissue (grey). **D**, Receptor–clock correlation before and after regressing out each
donor's mean NR3C1 across their tissues, in the five cohort–tissue pairs with several tissues per
person. **E**, Cultured islet preparations (GSE159984; 45 samples); dashed line, least-squares fit.

## Figure 3 | Hydrocortisone, but not acute insulin, shifts clock gene expression

**A**, Sampling scheme of the two perturbations. **B**, Change in clock gene expression after three
hours of hyperinsulinaemic clamp (GSE13070 muscle) compared with the change over a morning interval
in three cohorts with recorded sampling times. **C**, Paired fasting-to-hyperinsulinaemia change in
adipose tissue (GSE26637; 10 people); each dot one person, bars medians. **D**, Per-person change in
PER1 under hydrocortisone against placebo (GSE148640) and under clamp against fasting (GSE26637).
**E**, Hydrocortisone minus placebo in adipose tissue, both biopsies at 07:00 (GSE148640; 10
people); numbers above, participants changing upwards. **F**, Four days of dexamethasone in skeletal
muscle (GSE17992; 10 men; two-colour arrays, log ratios). **G**, Cushing syndrome adipose tissue
(GSE66446), centred on the control median. **H**, The nine clock genes in insulin-resistant adipose
tissue, Cushing syndrome, hydrocortisone and three months of thiazolidinedione. **I**, The same
clock genes in insulin-resistant adipose tissue against their change after three months of
thiazolidinedione (GSE13070); dashed line, least-squares fit.

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

## Figure 5 | A forced oscillator model of glucocorticoid entrainment

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

---

# Extended Data captions

## Extended Data Fig. 1 | The receptor–clock association in GTEx under technical covariates

**A**, Spearman correlation between NR3C1 and the clock-output score in three GTEx v11 tissues,
adjusted for sex, age and Hardy scale, and additionally for ischaemic time and RIN. **B**,
Correlation of ischaemic time and RIN with the clock-output score in the same tissues. **C**,
Visceral adipose tissue (*n* = 587); dashed line, least-squares fit. Expression is centred per
sample.

## Extended Data Fig. 2 | Inferred clock amplitude in two muscle cohorts

**A–C**, Inferred amplitude per person in GSE13070 (51 resistant, 18 sensitive) and GSE182120 (25
T2D, 24 NGT) under three estimators: projection onto an 800-gene template, clock genes with anchored
ordering, and clock genes with free ordering. Grey, sensitive or NGT; orange, resistant or T2D;
bars, medians; *P*, Mann–Whitney. **D**, Difference in medians for each estimator and cohort.
**E**, The two estimators compared across all samples. **F**, Ratio of the 5–95 percentile range
between groups; *P* from 3,000 permutations, with the Levene *P* (median-centred) alongside.
Coloured markers, *P* < 0.05.

## Extended Data Fig. 3 | Preprocessing decisions

**A**, Receptor–clock correlation against the median correlation of a random gene with the
clock-output score, in ten cohort–tissue matrices, before and after centring each sample; lines
join the same matrix. **B**, Genes matched in MSigDB v2026.2 (H + C2 + C5) using the platform
symbol and after resolving Entrez identifiers to current symbols; numbers, genes gained. **C**,
Probes carrying a retired symbol, by platform.

## Extended Data Fig. 4 | Regulatory architecture of the glucocorticoid axis in adipose tissue

**A**, Chord diagram of the five transcription factors whose inferred activity (DoRothEA A–C
regulons) correlates most strongly with the glucocorticoid signalling score, and their curated
targets with |ρ| > 0.35 against that score. Upper arcs, factors, with the correlation of their
inferred activity given outside; lower arcs, targets, coloured by the sign of their correlation;
chords, curated factor-to-target edges. **B**, Correlation of each of the 38 FOXO3 targets with
the glucocorticoid signalling score; TXNIP and FBXO32, referred to in the text, are in bold. AAGMEx adipose tissue, *n* = 232, sample-centred expression. The
ranking of factor activity against clock output is Fig. 4A and is not repeated here. Full regulon
in Supplementary Table S18.

## Extended Data Fig. 5 | Within-person response to the clamp and to thiazolidinedione

**A–D**, Median per-person change of each gene-set score in the four contrasts of GSE13070, with
the Wilcoxon *P* and the percentage of people changing downwards. Coloured bars, *P* < 0.05;
*n*, pairs.

## Extended Data Fig. 6 | An empirically derived glucocorticoid signature

**A**, Per-person response to hydrocortisone against placebo across 12,641 genes (GSE148640; 10
people); highlighted, the 68 genes selected for the empirical signature. **B**, Genes by the number
of participants changing in the same direction. **C**, Score difference for the empirical signature
and for the pre-specified nine-gene set in four cohorts not used to build either.

## Extended Data Fig. 7 | Model robustness

**A**, Ensemble amplitude against effective forcing for four integration steps. **B**, Level and
amplitude differences between low and high coupling, for a morning window and for uniform sampling
hours; mean ± s.d. of five seeds. **C**, Ensemble synchrony against effective forcing for three
values of the dispersion of intrinsic cellular periods.
