# The run the reported numbers come from

This file is written by `tests/write_reference_run.py` from the `run_manifest.json` of the run; it is not
edited by hand. The numbers of the manuscript, its figures and its Supplementary Tables are those of
this one complete run of the pipeline, from the public downloads to the last check.

| | |
|---|---|
| release | 2.13.9 |
| commit | `1301046a4281362c432cdadb5ba9bdeae585379c` |
| run identifier | `20261122T202227Z-1301046a-2fae26` |
| profile | full (default settings) |
| started | 2026-11-22T20:22:27+00:00 |
| finished | 2026-11-22T21:09:27+00:00 |
| duration | 47.0 minutes |
| input files | 75 (sha256 of their hashes `0e34a95a713832f3eb20f178cc54de54d1cf54a8cea47e13284c61c3ae710a79`) |
| output files | 181 (sha256 of their hashes `a2f61d910a2a670621ea846bd5ac466bf7bcd29c07e17c52bc2ba180cfe42267`) |
| islet counts | NCBI |
| uncommitted changes in versioned files | none |
| CHIRAL (step 06) | computed |
| ZeitZeiger and tauFisher (step 02f) | computed |
| final check | `verify_package.py` clean; the run reproduces every result of the release, the external benchmarks included |

The record of the run is its `run_manifest.json`, with `input_files.tsv` and `output_files.tsv`: the
sha256 of every input and of every output. They are deposited with the archived release, so that a
later reproduction can be compared file by file.

The history of this repository was rewritten on 17 November 2026 (see the entry of v2.13.2 in
`CHANGELOG.md`); commits are identified as they are after the rewriting.
