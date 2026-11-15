# posthoc/

POST HOC, EXPLORATORY. Written on 14 Nov 2026 after the confirmatory run. Nothing in this folder is
confirmatory and nothing in it was registered. Reports: `POSTHOC.md` (part 1) and `POSTHOC_2.md` (part 2).
ASCII, no em dashes.

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

Part 2 (measurement or biology; report: `POSTHOC_2.md`):

| script | what it does | output |
|---|---|---|
| p1_sweep.py with BANC_EDGE_FILE | the same sweep on the BANC v2 edge list (the table of the BANC paper) | results/p1_sweep_banc_v2.json, .log |
| p7_synapse_table.py | one pass over the BANC v3 synapse table: per neuron, synapses by kind of partner (identified or not); checks the counts against the registered edge list | results/p7_synapse_table.json, .log, p7_neuron_counts.npz |
| p8_capture_positions.py | share of synapses that end on identified neurons, by set and family; synapses of the set in the documented problem regions; profile along the sectioning axis | results/p8_capture_positions.json, .log |
| p12_family_tables.py | per family and per extrinsic type: synapses per neuron in the three data sets (run before p9) | results/p8_family_totals.json, p11_extrinsic_types.json, p12_family_tables.log |
| p9_types.py | 3,876 cell types present in the three data sets: synapses made and received per neuron, by transmitter; two-way fit of the family-block ratios | results/p9_types.json, .log |
| p10_neuropils.py | the same ratios by the neuropil that holds most of a type's output, beside the documented recall of the detector per cutout; by cell class | results/p10_neuropils.json, .log |
| p11_regional.py | synapses inside the EB set by neuropil in BANC and FlyWire; regional thinning of FlyWire | results/p11_regional.json, .log |
| p13_v2_extrinsic.py with BANC_EDGE_FILE | the extrinsic-type table of p12 on the BANC v2 edge list | results/p13_v2_extrinsic.json, .log |

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

Part 2 needs further files (none of them is hemibrain):

    # BANC bucket, https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/
    compiled_data/banc_888/banc_888_edgelist_simple_v2.feather          (md5 OUQG+Km98JPIlfla/09sSQ==)
    neuron_connectivity/v888/synapses_v3_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet   (5.6 GB)
    neuron_connectivity/v888/banc_888_edgelist_neuropil_v3.feather      (5.1 GB)
    compiled_data/fafb_783/fafb_783_synapses.parquet                    (2.0 GB)
    # BANC Dataverse deposit, doi:10.7910/DVN/7WTH1N
    251013_synapse_evaluation.csv, banc_problem_regions.csv
    # male CNS flat connectome: body-neurotransmitters (saved as body-nt.feather in MALE_DIR)

    BANC_EDGE_FILE=<...>/banc_888_edgelist_simple_v2.feather python3 p1_sweep.py banc _v2 1 2 3 4 5 6 8 10 12 15
    python3 p7_synapse_table.py <synapses_v3 parquet>
    python3 p8_capture_positions.py <banc_problem_regions.csv>
    python3 p12_family_tables.py
    python3 p9_types.py
    python3 p10_neuropils.py <edgelist_neuropil_v3 feather> <251013_synapse_evaluation.csv>
    python3 p11_regional.py <edgelist_neuropil_v3 feather> <fafb_783_synapses.parquet>

Notes on the runs of 14 Nov 2026 (Python 3.13.18, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, pyarrow 26.0.1):
- The first male sweep stopped at threshold 12 with a division by zero in the frozen `set_stats` (the two
  leading eigenvalues of the mushroom-body set coincide there). The traceback is kept in
  results/p1_sweep_male.log. p1_sweep.py now records such a set as an error and goes on; thresholds 12, 15
  and 20 of the male are in results/p1_sweep_male_hi.json.
- p5_attribution.py was first written with per-edge string arrays and stopped by hand for its memory use
  before it produced results; the version here uses integer codes.
- Part 2: the first version of p7_synapse_table.py ran out of memory and was rewritten (sorted merge instead of
  repeated sorting); its second run is the one reported. Its check passed: the counts between graph neurons equal
  the registered edge list. The file results/p7_eb_synapses.npz (10 MB, every synapse with an end in the EB set)
  is not in the repository; p7 writes it again.
- The confirmatory files in `results/` of the repository were not touched and the registered analysis was not
  rerun. Hemibrain was not opened.
