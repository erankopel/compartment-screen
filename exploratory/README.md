# Tests for the outline "A necessary-condition screen for self-measuring systems"

Run on 14 Nov 2026 (second batch the same day: T1, T1c, T2c). ASCII, no em dashes. Nothing here is
pre-registered: T0, T1 (parts a and b), T2 (parts 1 to 3), T2c and T4a are arithmetic, derivations or
reproductions; T1c, T2 part 4, T3 and T4b are EXPLORATORY.

| script | what it does | data | status |
|---|---|---|---|
| t0_floor.py | Margolus-Levitin and Landauer floors at 310 K, Planck capacity ceiling, product rule (screen bits x round rate <= P / kT ln 2) | none | arithmetic |
| t1_carrier.py | Tier 1 carrier screen: runs verify.py of the reference implementation (Zenodo 10.5281/zenodo.21837368), recomputes Table 3 and Table 2 of BioSystems 270, 106004, and derives the screen in rounds (phase alphabet after n cycles; cycles with at least one bit = pi Q / 16) | none (needs the archive) | reproduction + derivation |
| t1c_transmitters.py | flow-weighted predicted transmitter of the signals into, inside and out of the 727-neuron set, with comparison sets | male CNS v1.0 | EXPLORATORY |
| t2_record_rounds.py | three-qubit loop of arXiv:2608.13764 built from its circuit: closed forms, reference point, record bound on 20,000 random loops, index, index against record | none | reproduction + exploratory frontier |
| t2_core.py | the function definitions of t2_record_rounds.py, importable | none | helper |
| t2b_slice_variants.py | the 81 x 81 slice of arXiv:2609.27908 Sec. 6 under four guessed readings of the grid (all wrong, kept for the record) | none | superseded by t2c |
| t2c_slice_published_grid.py | the same slice on the grid of the archived figure script (Zenodo 10.5281/zenodo.22830716): theta/pi in [0.01, 0.25], phi/pi in [0.01, 0.49]; reproduces every published count | none | reproduction (exact) |
| t4a_frame_mismatch.py | boundary bit prepared in one frame and measured in another: p = sin^2(theta/2), eta_TV = abs(cos theta), eta_KL = cos^2 theta, r rounds, wide screens | none | exact + numeric check |
| t3_subsystem.py | the 727-neuron cut of E1 (maleCNS-depth) as a compartment: escape, residence, internal relaxation, interfaces, conditioned-walk Dobrushin coefficient, comparison sets, inside/outside TV | male CNS v1.0 | EXPLORATORY |
| t4b_exit_channel.py | exit channel of the cut (first exit time and neuron), 240 sampled starts, neuron and cell-type level | male CNS v1.0 | EXPLORATORY |
| t4b_exit_channel_allpairs.py | the same for all 727 starts, neuron level, t <= 96 | male CNS v1.0 | EXPLORATORY |
| screen_lib.py | the screen as a library: one function per statistic (closure, stays / mixes, conditioned-walk Dobrushin coefficient, exit channel, labelled flow shares, reach), used unchanged for every data set | any | candidate for the frozen confirmatory code |
| transfer_screen.py | applies screen_lib.py to a connectome in the compiled layout of the BANC data tutorial, with sets defined by names only (crosswalk_types.json); records SHA-256 of inputs and code | compiled meta + edge list | candidate for the frozen confirmatory code; rehearsed on FlyWire v783 |
| crosswalk_types.json | family rules for the EB set, class names for CX, MB and AL, and the amine type list | written by t6 | candidate for the freeze |
| t6_male_typed_nulls.py, t6b_male_nulls.py | the screen on the male CNS with name-defined sets at >= 1 and >= 5 synapses, and under the null models N1, N2, N3 of the repository (t6b reruns the nulls after the first run hit the memory limit) | male CNS v1.0 | EXPLORATORY |
| t5_pit_union.py | central complex and mushroom body: is the union a compartment, and how fast does a walker leaving one reach the other (against a stationary walker and an antennal-lobe control) | male CNS v1.0 | EXPLORATORY |
| t4b_exit_comparison.py | the exit-channel statistic for comparison sets (central complex, mushroom body, Kenyon cells, antennal lobe), 240 starts each | male CNS v1.0 | EXPLORATORY |

## To rerun

    pip install numpy scipy pandas pyarrow
    python3 t0_floor.py; python3 t4a_frame_mismatch.py; python3 t2_record_rounds.py; python3 t2c_slice_published_grid.py
    python3 t1_carrier.py /path/to/paper3-criteria   # unpacked Kopel_2026_Criteria_Code_v1.0.0.zip

    git clone https://github.com/erankopel/maleCNS-depth
    export MALECNS_REPO=$PWD/maleCNS-depth MALECNS_DATA=$PWD/data
    python3 $MALECNS_REPO/code/fetch_data.py      # 1.1 GB, md5 checked against the pre-registration
    python3 $MALECNS_REPO/code/build_graph.py
    python3 t3_subsystem.py                        # about 10 min, 2 threads
    python3 t4b_exit_channel.py                    # about 2 min
    python3 t4b_exit_channel_allpairs.py           # about 15 min
    python3 t4b_exit_comparison.py                 # about 6 min
    python3 t1c_transmitters.py                    # about 3 min
    python3 t5_pit_union.py                        # about 3 min
    python3 t6_male_typed_nulls.py                 # about 35 min; needs 5 GB free memory
    # rehearsal on FlyWire v783 (files from gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/fafb_783/)
    python3 transfer_screen.py fafb_783_id fafb_783_meta.feather fafb_783_simple_edgelist.feather results/t7_flywire_rehearsal.json 1 5

The connectome scripts import code/confirmatory/common.py of the repository unchanged and take the cut
from exploratory/results/e1_slow_mode.json (key w1, cut_body_ids). Versions used here: Python 3,
numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, pyarrow 26.0.1 (the confirmatory runs of the repository used
numpy 2.4.4, scipy 1.17.1, pandas 3.0.2, pyarrow 25.0.1 and one BLAS thread; the runs here used two).
Seeds: 20261114 everywhere.

## Freeze package

prereg_screen_banc.zip holds PREREG.md, the three code files (screen_lib.py, transfer_screen.py,
crosswalk_types.json) and MANIFEST.sha256. The copies there are the ones to deposit; the final files
reproduce the FlyWire rehearsal exactly (results/t7_flywire_rehearsal.json records their hashes).
The BANC metadata table was read for names only; the BANC edge list has not been downloaded.

## Known gaps

- Nothing is frozen and no confirmatory test was run. BANC and hemibrain connectivity have not been
  downloaded or opened. FlyWire v783 was used for the rehearsal and so cannot serve as confirmatory data.
- The amine statistic does not carry across data sets (transmitter predictions and type names differ);
  it is reported, not scored.
- T1: criterion C5 for chemical labels is not evaluated (needs measured signal lifetimes per modulator).
