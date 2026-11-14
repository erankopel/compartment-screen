# compartment-screen

Pre-registration and code for a compartment screen on fly connectomes. ASCII, no em dashes.
The repository name is a placeholder; rename it before the first push if you prefer another.

## Status

| stage | state |
|---|---|
| Registration (this release, v0.1.0-prereg) | prepared 14 Nov 2026; DOI and publication time to be entered here |
| Confirmatory run on BANC release 888 | NOT RUN. The BANC edge list has not been downloaded. |
| Results release | after the run |

## What is registered

`prereg/PREREG.md` states the data, the graph, the sets (by names only), the statistics, five criteria
(C1 to C4 for the compartment hypothesis, C5 for the coupling of mushroom body and central complex),
the decision rule, what is reported without being scored, and the protocol. The three code files in
`prereg/` are the frozen code. `prereg/MANIFEST.sha256` fixes all four files:

    cd prereg && sha256sum -c MANIFEST.sha256

## Layout

- `prereg/` the registration: PREREG.md, screen_lib.py, transfer_screen.py, crosswalk_types.json, MANIFEST.sha256.
- `exploratory/` everything that led to it and is NOT confirmatory: scripts t0 to t7, their results and
  logs (male CNS v1.0, null models, FlyWire v783 rehearsal), with its own README.
- `results/` empty until the run.
- `DEVIATIONS.md` empty until the run.

## The run (once, after the registration is public)

    B=https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888
    curl -O $B/banc_888_meta.feather
    curl -O $B/banc_888_edgelist_simple_v3.feather
    cd prereg
    python3 transfer_screen.py banc_888_id ../banc_888_meta.feather ../banc_888_edgelist_simple_v3.feather results_banc.json 1 5
    cp results_banc.json ../results/

Environment of the rehearsal: Python 3.13, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, pyarrow 26.0.1, two
BLAS threads. The script writes the md5 and sha256 of both inputs and the sha256 of the code into the
result file; the md5 values must equal those in PREREG.md, Section 2.

## Related work

- maleCNS-depth (arXiv:2609.33054), the study this continues.
- Bates et al. (2026), the BANC connectome, doi:10.1038/s41586-026-10735-w.
- Klein and Barron (2025), phenomenal interface theory, doi:10.1098/rstb.2024.0301.

## Licence

MIT for the code. The connectome data belong to their producers and are not redistributed here.
