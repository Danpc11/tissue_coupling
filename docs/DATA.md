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
