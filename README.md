# compartment-screen

Pre-registration and code for a compartment screen on fly connectomes. ASCII, no em dashes.

## Status

| stage | state |
|---|---|
| Registration (release v0.1.0-prereg) | doi:10.5281/zenodo.24223802 (concept doi:10.5281/zenodo.24223801), record created 2026-11-14T15:59:33Z |
| Confirmatory run on BANC release 888 | one run, 2026-11-14, 16:03:09 to 16:07:21 UTC, after the registration |
| Results (release v0.2.0) | H-EB confirmed (C1 to C4 hold), H-PIT confirmed (C5 holds), primary analysis; see RESULTS.md, including the values at 5 synapses, where C2 and C3 would not hold |

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
- `RESULTS.md` verdicts and every reported number.
- `results/` the result file, the run log and the wrapper used to start the run.
- `DEVIATIONS.md` none from the registered analysis; four notes for the record.

## The run (made once, on 14 Nov 2026; the commands are kept for reference and must not be repeated as a confirmatory run)

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
