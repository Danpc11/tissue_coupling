# Data

## Public downloads expected under `data/raw/`

Nineteen GEO series and GTEx v11. Accessions, tissue, platform, sample number and the role each
series plays are listed in Supplementary Table S1 of the manuscript:

GSE11868, GSE13070, GSE17992, GSE25462, GSE26637, GSE41168, GSE66446, GSE78721, GSE87761,
GSE95674, GSE95675, GSE104674, GSE107537, GSE108539, GSE129843, GSE148640, GSE159984, GSE182117,
GSE182120.

GSE107537 is the external control for the estimator: a controlled 10-h delay of the habitual sleep
period with recorded sampling times, whose published result is dampened amplitude at a preserved
phase. Its platform annotation (GPL24324) is needed to resolve transcript cluster identifiers to
gene symbols.

GSE95674 and GSE95675 are the two tissues of one study (AAGMEx). GSE182117 and GSE182120 are
separate studies despite adjacent accessions.

For GTEx, the gene-reads GCT of each tissue plus the subject-phenotype and sample-attribute
annotation files. Individual-level controlled-access data require dbGaP authorisation; the
analyses here use the open-access expression matrices and annotations.

The DoRothEA A–C regulons (`dorothea_ABC.tsv`, exported from OmniPath) are also expected under
`data/raw/`. Platform annotation files for the microarray series, and the NCBI gene annotation used to resolve
Entrez identifiers to current symbols, are also expected here. `00_prepare_data.py` lists each file
it needs and fails with the missing name rather than silently skipping it.

## Versioned inputs under `data/reference/`

These are inputs that `run_all.sh` does not regenerate. Each is here for a stated reason, not for
convenience.

| file | why it is versioned |
|---|---|
| `GSE13070_IR_vs_sensitive_log2.tsv` | the insulin-resistant against insulin-sensitive contrast of GSE13070 is unbalanced (35 against 6 in adipose tissue) and was not reproduced consistently; it is used only to define a profile of change compared with other profiles, and is excluded from the differential expression of Fig. 1 |
| `signatures_reference.tsv` | the twelve-contrast signature table. `23_build_signatures_table.py` regenerates three of its columns from public data and **checks each against this reference**, aborting if any disagrees; the rest are carried over |
| `template_2harm.tsv` | the first- and second-harmonic template, derived from cohorts with recorded sampling times |
| the remaining files | per-cohort intermediates that a single step produces and several figures read; versioned so a figure can be redrawn without rerunning the whole pipeline |

## Convention

Scores computed across individuals use per-sample centred expression. Differences between groups,
including the signature table, are computed without centring. See the README.

## Reconstruction of the preparation step (branch `repro-hugo`, November 2026)

The published package did not include the step that turns the public downloads into the matrices under `results/new/*.pkl`. This branch adds it:

| script | what it produces |
|---|---|
| `download_public.sh` | every public input listed above, from GEO FTP, the GTEx bucket, Broad, OmniPath and recount3; `T2D_ARCHS4=1` also fetches the ARCHS4 HDF5 |
| `python/analyses/00b_build_matrices.py` | the 26 `*.pkl` matrices (annotation by Entrez for microarrays, log2 CPM for RNA-seq, one donor per column for GTEx), the `*_meta.pkl` tables, `GSE108539_counts.tsv.gz` from recount3, and the copy of `data/reference/` into `results/new/` |
| `python/analyses/00d_reference_tables.py` | regenerates the 23 tables of `data/reference/` that had no producer, writes them to `results/new/` and compares each with the versioned copy (`results/new/reference_tables_check.tsv`) |
| `python/analyses/01c_islets_deseq2.py` | `GSE159984_islets_T2D_vs_ND.tsv`, the PyDESeq2 contrast of the 85 islet donors (27 T2D against 58 non-diabetic) that Fig. 1b-c reads |
| `python/analyses/00c_archs4_counts.py` | `GSE159984_raw_counts.tsv.gz` from ARCHS4, read by HTTP range if the HDF5 is not on disk |

Substitutes for files that are only served by `www.ncbi.nlm.nih.gov`, which answers non-interactive clients with a CAPTCHA:

- `Human_GRCh38_p13_annot_tsv.gz` is derived from `Homo_sapiens.gene_info.gz` (NCBI Gene FTP); the scripts read only `GeneID` and `Symbol`.
- The NCBI-generated counts (`GSE108539_raw_counts_GRCh38_p13_NCBI.tsv.gz`, `GSE159984_raw_counts_GRCh38.p13_NCBI.tsv.gz`) download with a plain `curl` when NCBI is not rate-limiting the client; when it answers with a CAPTCHA page instead, `00b_build_matrices.py` falls back to recount3 for GSE108539 (SRP127581) and `00c_archs4_counts.py` can supply the 44 cultured preparations of GSE159984 from ARCHS4 (the 86 donor islets, GPL9115, are not in ARCHS4 nor recount3).
- Platform tables (GPL570, GPL15207, GPL17586, GPL24324) are the `!platform_table` block of the family SOFT file, read in streaming mode.
- CHIRAL (step 06) needs the Python port by Andrea Salati, which is not on PyPI: `git clone --depth 1 https://github.com/AndreaSalati/pyCHIRAL.git` and `PYCHIRAL_SRC=<clone>/src` (or a clone of scRITMO and `SCRITMO_SRC`); it also needs `tqdm`.
