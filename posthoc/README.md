# posthoc/

POST HOC, EXPLORATORY. Written on 14 Nov 2026 after the confirmatory run. Nothing in this folder is
confirmatory and nothing in it was registered. Report: `POSTHOC.md`. ASCII, no em dashes.

Question: why do the registered criteria C2 and C3 hold in BANC with all edges and not at 5 synapses?

| script | what it does | output |
|---|---|---|
| load3.py | loader for BANC 888, male CNS v1.0 and FlyWire v783 with the registered names-only sets; imports the frozen `prereg/screen_lib.py` unchanged | |
| p1_sweep.py | EB statistics at thresholds 1 to 20 in each data set | results/p1_sweep_<data>.json, .log |
| p2_scale.py | synapse counts per neuron, per edge and per family block; BANC proofread flags | results/p2_scale.json, .log |
| p3_anatomy.py | slow mode, classes, unmixed starts, diagnostic subsets at given thresholds | results/p3_anatomy_<data>.json, .log |
| p4_thinning.py | thinning of the male or FlyWire data to the BANC density, uniformly and with the BANC block pattern; three seeds | results/p4_thinning_<source>.json, .log |
| p5_attribution.py | the block pattern applied to one group of connections at a time; one seed | results/p5_attribution_<source>.json, .log |
| p6_summary.py | prints the tables of POSTHOC.md from the result files | results/p6_summary.txt |

Run (two BLAS threads, about 7 GB of memory; run the male jobs one at a time):

    export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
    export BANC_DIR=<folder with banc_888_meta.feather and banc_888_edgelist_simple_v3.feather>
    export FLYWIRE_DIR=<folder with fafb_783_meta.feather and fafb_783_simple_edgelist.feather>
    export MALE_DIR=<folder with neuron_edges.npz and body-annotations.feather>
    for d in banc male flywire; do python3 p1_sweep.py $d; done
    python3 p2_scale.py
    python3 p3_anatomy.py banc 1 2 3 5; python3 p3_anatomy.py male 1 5 10; python3 p3_anatomy.py flywire 1 5
    python3 p4_thinning.py male; python3 p4_thinning.py flywire        # needs results/p2_scale.json
    python3 p5_attribution.py male; python3 p5_attribution.py flywire  # needs results/p2_scale.json
    python3 p6_summary.py > results/p6_summary.txt

Notes on the runs of 14 Nov 2026 (Python 3.13.18, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, pyarrow 26.0.1):
- The first male sweep stopped at threshold 12 with a division by zero in the frozen `set_stats` (the two
  leading eigenvalues of the mushroom-body set coincide there). The traceback is kept in
  results/p1_sweep_male.log. p1_sweep.py now records such a set as an error and goes on; thresholds 12, 15
  and 20 of the male are in results/p1_sweep_male_hi.json.
- p5_attribution.py was first written with per-edge string arrays and stopped by hand for its memory use
  before it produced results; the version here uses integer codes.
- The confirmatory files in `results/` of the repository were not touched and the registered analysis was not
  rerun. Hemibrain was not opened.
