# The run the reported numbers come from

The numbers of the manuscript are those of one complete run of the pipeline on the authors'
cluster, from the public downloads to the last check.

| | |
|---|---|
| release | 2.13.2 |
| commit | `b6a76416b066dc437457c1d7eb71d850b177029d` |
| run identifier | `20261117T221711Z-b6a76416-16f4e1` |
| profile | full (default settings) |
| steps | 52 |
| input files | 76 |
| output files | 176 |
| duration | 49 minutes, with 32 processes in the comparison of mechanisms |
| islet counts | NCBI |
| final check | `verify_package.py` clean; "This is a complete reproduction of release 2.13.2" |
| not run | the comparison with ZeitZeiger and tauFisher (step 02f), for lack of the two R packages on that machine; its numbers come from a run on another machine |

Releases 2.13.3 and 2.13.4 change checks, messages and documentation, not the analyses; the run was
not repeated with them.

The record of that run is its `run_manifest.json`, with `input_files.tsv` and `output_files.tsv`:
the sha256 of every input and of every output. They are to be deposited with the archived release,
so that a later reproduction can be compared file by file.

The history of this repository was rewritten on 17 November 2026 (see the entry of v2.13.2 in
`CHANGELOG.md`); the commit above is an identifier after the rewriting.
