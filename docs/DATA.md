# Data

## Public downloads expected under `data/raw/`

Eighteen GEO series and GTEx v11. `bash download_public.sh` fetches all of them, with the platform
tables and annotation files the pipeline reads, and `python/analyses/00_prepare_data.py` lists
each one and stops the run if any is missing.

GSE13070, GSE17992, GSE25462, GSE26637, GSE40589, GSE41168, GSE66446, GSE78721, GSE87761, GSE95674,
GSE95675, GSE107537, GSE108539, GSE129843, GSE148640, GSE159984, GSE182117, GSE182120.

**Two series listed in Supplementary Table S1 are not used.** GSE104674 (adipose tissue at four
circadian times) and GSE11868 (a 3 h clamp in six lean subjects) were downloaded and listed by the
releases up to v2.9.1, and no step read them. They are no longer downloaded. The series the
analyses use are therefore seventeen of the nineteen in that table, plus GSE40589.

Seventeen of them are listed, with tissue, platform, sample number and role, in
Supplementary Table S1 of the manuscript. The eighteenth, GSE40589, was added in v2.2.0 as a
time-matched control for the clamp (last section); its place in the manuscript is not decided.

GSE107537 is the external control for the estimator: a controlled 10-h delay of the habitual sleep
period with recorded sampling times, in blood mononuclear cells. Its platform annotation (GPL24324)
is needed to resolve transcript cluster identifiers to gene symbols.

GSE95674 and GSE95675 are the two tissues of one study (AAGMEx). GSE182117 and GSE182120 are
separate studies despite adjacent accessions.

For GTEx, the gene-reads GCT of each tissue plus the subject-phenotype and sample-attribute
annotation files. Individual-level controlled-access data require dbGaP authorisation; the
analyses here use the open-access expression matrices and annotations.

The DoRothEA A–C regulons (`dorothea_ABC.tsv`, exported from OmniPath), the MSigDB collections and
the NCBI gene annotation used to resolve identifiers to current symbols are fetched by the same
script.

## From the downloads to the matrices

Releases up to v1.3.1 did not include the step that turns the public downloads into the matrices
under `results/new/*.pkl`. It was reconstructed by Hugo Tovar and is part of `main` since v2.0.0:

| script | what it produces |
|---|---|
| `download_public.sh` | every public input, from GEO FTP, the GTEx bucket, Broad, OmniPath and recount3; `T2D_ARCHS4_DL=1` also fetches the ARCHS4 HDF5 |
| `python/analyses/00b_build_matrices.py` | the expression matrices and their `*_meta.pkl` tables: annotation by Entrez for the Affymetrix platforms, by platform symbol plus NCBI synonyms for Illumina HT-12, Ensembl to current symbol for RNA-seq and GTEx; log2 CPM with the library size taken before gene filtering; one donor per column for GTEx |
| `python/analyses/00d_reference_tables.py` | regenerates the tables of `data/reference/` and compares each with the versioned copy (`results/new/reference_tables_check.tsv`) |
| `python/analyses/01c_islets_deseq2.py` | the PyDESeq2 contrast of the 85 islet donors of GSE159984 |
| `python/analyses/00c_archs4_counts.py` | optional: counts of GSE159984 from ARCHS4 when NCBI does not serve them |

Decisions taken in that reconstruction that a reader should know:

- `Human_GRCh38_p13_annot_tsv.gz` is derived from `Homo_sapiens.gene_info.gz` (NCBI Gene FTP); the
  scripts read only `GeneID` and `Symbol`.
- Platform tables are the `!platform_table` block of the family SOFT file, read in streaming mode.
- GSE108539 is built from recount3 (SRP127581): coverage sums divided by the average read length
  of `sra.recount_qc`.
- GSE78721 has 130 samples from 95 people (35 contribute two depots and GEO gives no person
  identifier). One sample per person is used, with the depot as a covariate.
- GSE129843 is used from the authors' CPM table, as no raw counts are deposited.
- `gene_info` and the OmniPath export have no version; their date and sha256 are appended to
  `data/raw/CHECKSUMS.txt`.

### If NCBI answers with a CAPTCHA

The NCBI-generated counts are served by `www.ncbi.nlm.nih.gov`, which intermittently replies with
an HTML challenge page instead of the file. `download_public.sh` retries five times with
increasing waits and deletes any HTML it receives. GSE159984 is the one series the pipeline cannot
do without: if NCBI keeps refusing, `T2D_ARCHS4_DL=1 bash download_public.sh` followed by
`python python/analyses/00c_archs4_counts.py` builds the counts from ARCHS4 (a 62 GB download).

### Software that is cloned rather than installed

CHIRAL (step 06) is not on PyPI. `download_public.sh` clones its Python implementation,
github.com/AndreaSalati/pyCHIRAL, into `data/raw/pychiral` at commit `e1ea76e`, the revision used
in the reproduction of November 2026, and records it in `CHECKSUMS.txt`. `PYCHIRAL_COMMIT` selects
another revision; `PYCHIRAL_SRC` or `SCRITMO_SRC` point the pipeline at another copy. It also
needs `tqdm` and `anndata`, both in `requirements.txt`. Without it step 06 runs with our own
estimator only, marks the CHIRAL rows as not run and says so; `T2D_REQUIRE_CHIRAL=1` makes that
an error.

### The islet counts, and a fallback that is not equivalent

The islet estimates use the counts NCBI generates for GSE159984 (HISAT2 and featureCounts against
RefSeq), which `download_public.sh` fetches. `python/analyses/00c_archs4_counts.py` is a fallback
for machines that cannot reach that server: it extracts the same samples from ARCHS4, quantified
with kallisto against GENCODE, and writes them under the name the pipeline reads. The counts are
not identical, so islet numbers obtained that way are not the ones reported; nothing in the
outputs records which source was used, so whoever uses the fallback has to say so.

### R packages for the comparison with other methods

Step 02f compares the phase template with ZeitZeiger and tauFisher. Both are R packages that are
not on CRAN:

```bash
# R 4.3 with the dependencies, for instance from conda-forge and bioconda
micromamba create -p renv -c conda-forge -c bioconda r-base=4.3 r-abind r-bbmle r-data.table r-foreach \
    bioconductor-limma r-pma r-statmod bioconductor-sva r-metacycle r-scales r-dplyr r-stringr r-nnet r-magrittr r-fda
git clone --depth 1 https://github.com/hugheylab/zeitzeiger.git && renv/bin/R CMD INSTALL zeitzeiger
git clone --depth 1 https://github.com/micnngo/tauFisher.git   && renv/bin/R CMD INSTALL tauFisher
T2D_RSCRIPT=$PWD/renv/bin/Rscript bash run_all.sh
```

The versions used for the values in `CHANGELOG.md` were zeitzeiger 2.1.3 and tauFisher 1.0.0. Without
R the step runs our own estimator only and says so.

## Versioned inputs under `data/reference/`, and how they compare with their regeneration

These tables were shipped with the earlier releases, which had no code to produce most of them.
Step 02d now regenerates every one from the public data and compares it with the versioned copy.
Every figure is drawn from tables computed in the run. The copies stay here as the reference for
the comparison below, and for one more purpose, stated so that nobody has to find it in the code:
step 16 builds the signature table twice. `signatures_published_reference.tsv` starts from the
versioned copy of `signatures_reference.tsv` and replaces three of its twelve columns by their
recomputation; it is the table behind the published values. `signatures_fully_regenerated.tsv` has
the twelve columns recomputed. The same statistics are computed on both, the difference is written
gene by gene, and the figures read the second. `T2D_SIGNATURES=published` makes them read the
first. In the same way, Fig. 1f and 1g plot the leave-one-out validation of the run;
`T2D_LOO_COLUMN=amp_inf` plots the estimator as it was published.

The comparison of a complete run (11 November 2026, public data as downloaded that week) is below.
The correlation is between the versioned and the regenerated values over the rows they share. It
is the answer to which results of the earlier releases the present code reproduces.

| table | agreement | note |
|---|---|---|
| `GSE13070_IR_vs_sensitive_log2.tsv` | r = 1.000 in both tissues | the unbalanced contrast (35 against 6 in adipose tissue) is now regenerated exactly |
| `GSE26637_clamp_adipose_per_person.tsv` | r = 1.000 | |
| `observed_signature_for_model.tsv` | r = 1.000 | |
| `FOXO3_regulon_ED4.tsv` | r = 1.000 | |
| `AAGMEx_GCreceptor_vs_adiposity.tsv` | r = 1.000 | |
| `AAGMEx_NR3C1_enrichment*.tsv` | r = 0.996 and 0.998 | |
| `GRAS_signature.tsv` | r = 0.996 | |
| `AAGMEx_rank_vs_NR3C1.tsv` | r = 0.995 | |
| `time_calibration_vs_clamp.tsv` | r = 0.979 to 1.000 | |
| `AUDIT_global_axis.tsv` | r = 0.990 (uncentred), 0.984 (centred) | the median of the centred null agrees less (r = 0.49) |
| `GSE17992_dex_log2ratio_muscle.tsv` | r = 0.980 to 0.989 | |
| `GSE148640_HC_response_stats.tsv` | r = 0.955 to 0.983 | |
| `signatures_reference.tsv` | r = 1.000 in the five GSE13070 columns; 0.98 hydrocortisone and Cushing; 0.96 and 0.93 dexamethasone; 0.98 and 0.85 the two morning contrasts; 0.82 myotubes | all twelve columns are regenerated |
| `within_donor_decomposition.tsv` | r = 0.925 and 0.993 | |
| `ncbi_synonym_map.tsv` | 99.1% of shared entries identical | the NCBI file changes daily |
| `GSE129843_amplitude_per_person.tsv` | r = 0.70 to 0.95 | |
| `template_2harm.tsv` | r = 0.93 and 0.84 for the first harmonic, 0.84 and 0.64 for the second; 0.55 for r²; 0.49 and 0.08 for the two amplitudes | **reproduces only in part.** The regenerated template comes from recount3 counts (13,430 genes) and the versioned one from the NCBI counts of the original analysis (14,024). Every phase and amplitude estimate downstream uses the regenerated template |
| `amplitude_phase_distributions_two_cohorts.tsv` | amplitude: ρ = 0.71 to 0.78. Phase, compared as a circular quantity and allowing a global offset and either direction: concentration of the difference 0.06 to 0.22, median error 4.6 to 6.3 h once the offset is removed, where chance gives 6 h | **the phases do not reproduce, and it is not a matter of convention.** In these cohorts every biopsy was taken at about the same hour, so the true phases differ little and what the estimator returns for one sample is mostly noise: two builds of the template give unrelated phases to the same sample. This table feeds the dispersion analysis, whose phase results should not be relied on |
| `GSE182117_inferred_free_vs_measured.tsv` | r = −0.22 for inferred amplitude | **does not reproduce.** No figure reads it since v2.0.0, when Fig. 1f moved to the leave-one-out table |
| `msigdb_symbol_coverage.tsv` | not comparable | a five-row coverage summary computed on a different gene universe |
| `symbol_audit_renamed_genes.tsv` | not comparable | the versioned table is a subset whose selection rule was not recorded (427 rows against 10,926) |

What this means for the manuscript is recorded in the reproduction report kept with the manuscript
files, statement by statement. It is not in this repository.

## Convention

Scores computed across individuals use per-sample centred expression. Differences between groups,
including the signature table, are computed without centring. See the README.

## The time-matched control for the clamp (GSE40589)

Subcutaneous adipose biopsies taken 40 min before and 240 min after the start of an infusion that
begins at 08:00, with saline alone, with a euglycaemic hyperinsulinaemic clamp or with a
hyperglycaemic clamp (Tuvia et al., Diabetes 2021;70:1985). `download_public.sh` fetches the series
matrix and the platform table (`GPL6848_annot.txt`, taken from the SOFT file because GEO has no
curated annotation for that platform).

The series does not declare the participant. `28_clamp_vs_saline.py` pairs the two biopsies of an
experiment by the consecutive numbering of the sample titles and checks each pair with genes that
identify an individual (the Genetic and HLA sets of `python/lib/gene_sets.py`). A pair is excluded
from the primary analysis when its identity correlation falls below the 90th percentile of the
correlations between different pairs. The script also reports the analysis with every pair and
one that does not use the pairing at all.
