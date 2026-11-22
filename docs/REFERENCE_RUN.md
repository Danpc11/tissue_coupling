# The run the reported numbers come from

The numbers of the manuscript, its figures and its Supplementary Tables are those of one complete
run of the pipeline, from the public downloads to the last check, in the reference configuration
(profile `full`, default settings).

This file is rewritten by `tests/write_reference_run.py` from the `run_manifest.json` of that run:
release, commit, run identifier, number and hash of the inputs and of the outputs, duration, and
which external benchmarks were computed. Until it has been rewritten for the release being
archived, what can be said is this:

| | |
|---|---|
| Supplementary Tables and figures of the manuscript | written from a complete run of release 2.13.7 |
| earlier complete run on record | release 2.13.2, run `20261117T221711Z-b6a76416-16f4e1`, commit `b6a76416b066dc437457c1d7eb71d850b177029d` (52 steps, 76 inputs, 176 outputs, 49 minutes) |
| difference between the two | the pooled estimate is 0.415 (two-level) and 0.456 (three-level) in the run of 2.13.7 and was 0.416 and 0.458 in that of 2.13.2; two cohort estimates differ in the third decimal (GSE78721, 0.542 against 0.545; GSE25462, 0.553 against 0.550). The cause has not been identified |
| release 2.13.8 | changes how two supplementary panels are drawn, adds one table (`GRAS_validation.tsv`) and corrects labels; it changes no analysis whose numbers are in the text |

To record the run of a release:

```bash
bash run_all.sh                              # one complete run, default settings
python tests/write_reference_run.py          # rewrites this file from results/new/run_manifest.json
```

The script refuses a run that did not pass the verification, was not in the reference
configuration or had uncommitted changes in versioned files. `run_manifest.json`,
`input_files.tsv` and `output_files.tsv` of that run are deposited with the archived release, so
that a later reproduction can be compared file by file.

The history of this repository was rewritten on 17 November 2026 (see the entry of v2.13.2 in
`CHANGELOG.md`); the commit above is an identifier after the rewriting.
