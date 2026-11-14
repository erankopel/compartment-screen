# Pre-registration: a compartment screen on a second fly connectome (BANC)

Version of 14 Nov 2026, prepared for deposit. ASCII, no em dashes.
Status: thresholds and primary analysis approved by E. Kopel on 14 Nov 2026. NOT YET DEPOSITED.
The BANC edge list must not be downloaded or opened before this document and the three code files
are deposited with a public time stamp. The run starts only after that deposit.

## 1. Background (all of it exploratory)

- Male CNS v1.0 (the data of the maleCNS-depth study, arXiv:2609.33054). The smallest-conductance sweep
  cut of the synapse-flow chain at >= 1 synapse is a set of 727 neurons of the ellipsoid body and its
  anterior visual input. A walker stays in this set about twice as long as it takes to mix inside it,
  and the law of its first exit depends only weakly on where it started.
- Null models on the male data (generators of E3 of the maleCNS-depth repository): a null that keeps the
  synapse counts between every pair of cell types (N3) reproduces every statistic below; nulls that keep
  only superclass-to-superclass counts (N2) or only degrees (N1) do not. The statistics are therefore a
  property of the wiring between cell types.
- FlyWire v783 (female brain) was used for a full rehearsal of the code in this package with sets defined
  by names only. All criteria below hold there. FlyWire also informed the thresholds, so it is not
  confirmatory.
- Biological hypothesis: phenomenal interface theory (Klein and Barron 2025, doi:10.1098/rstb.2024.0301;
  Barron and Klein 2016, doi:10.1073/pnas.1520084113) names the central complex and the mushroom bodies,
  with their aminergic networks, as the insect interface. This registration tests structural statements
  only. It does not test any statement about experience.

## 2. Data

BANC, release 888 (Bates et al. 2026, doi:10.1038/s41586-026-10735-w), female brain and nerve cord,
compiled flat files of the BANC data tutorial:
  https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/

| file | bytes | md5 (base64, bucket header) | crc32c | generation | last modified |
|---|---|---|---|---|---|
| banc_888_meta.feather | 57550730 | YjyrACb8v8r2J/2q+en7HQ== | +decfQ== | 1794372660597213 | 2026-11-11 04:51:00 GMT |
| banc_888_edgelist_simple_v3.feather | 359161658 | CFQrB3HbdBjtR0vmDcmIbA== | zvJNpA== | 1786578377086929 | 2026-08-12 23:46:17 GMT |

sha256 of banc_888_meta.feather as downloaded on 14 Nov 2026:
60a7c2252841c9f05aef0803c84bc3be7802b211f6bd2a40cef123ee3ad6265f
The run script records md5 (base64) and sha256 of both files. If either md5 differs from the table, the
run stops and the difference is logged as a deviation before anything else is done.

### What was inspected before the freeze

- File headers of both files (size, hashes, dates). The edge list was not downloaded.
- From the metadata table, names only: the column names; the columns banc_888_id, super_class,
  cell_class, cell_type; and value counts of region, status and proofread. No transmitter column and
  no connectivity.
- Result of that inspection: 188,508 rows. One change was made to the graph rule because of it: the
  super class not_a_neuron (195 rows) is excluded besides glia and trachea.

## 3. Graph

- Nodes: rows of the metadata table whose super_class is not null and is not glia, trachea or
  not_a_neuron (146,769 rows; 139,458 of them flagged proofread; 118,670 with a cell type). The
  proofread flag is not used.
- Edges: rows of the edge list with both ends among the nodes, weight = count.
- Chain: random walk on the giant strongly connected component, P = row-normalised counts (rows are
  presynaptic). Stationary distribution by power iteration from the uniform vector, L1 tolerance 1e-12,
  at most 5000 iterations; the last change is reported.
- Primary analysis: all edges (count >= 1). The same numbers at count >= 5 are reported, not scored.

## 4. Sets (names only, file crosswalk_types.json)

- EB: a neuron belongs to EB if any comma-separated part of its cell_type matches one of
  ^ER[1-6]  ^EPG$  ^PEG$  ^PEN_  ^EL$  ^ExR[12346](?![0-9])  ^TuBu  ^MeTu2  ^LPsP$  ^GLNO$  ^AOTU046$  ^TuTuB
  (male: 691 neurons, identical to the 53 male types with at least half of their neurons in the sweep
  cut; FlyWire: 683; BANC by the names-only check: 595 neurons in 54 types).
- CX: cell_class in central_complex_intrinsic_neuron, central_complex_input_neuron,
  central_complex_output_neuron, plus EB (BANC: 3,164).
- MB: cell_class in kenyon_cell, mushroom_body_output_neuron, mushroom_body_dopaminergic_neuron (BANC: 4,957).
- AL: cell_class in antennal_lobe_projection_neuron, antennal_lobe_local_neuron (BANC: 1,128).
- Sets are restricted to the giant strongly connected component. No neuron is added or removed by hand.

## 5. Statistics (file screen_lib.py)

For a set M with stationary mass pi(M):
- loss per step: the share of pi(M) that leaves M in one step.
- stays: 1 / (1 - lam1), lam1 the leading eigenvalue of P restricted to M (the walk stopped on leaving).
- mixes: 1 / (1 - |lam2| / lam1), lam2 the second eigenvalue by modulus. Dense eigenvalues for sets of up
  to 1600 neurons, ARPACK (k = 6, ncv = 60, tol 1e-9, seeded) above.
- Dobrushin coefficient of the walk conditioned to stay in M (Doob transform with the right Perron
  vector; neurons with a zero entry are dropped and counted): the largest total variation distance
  between rows of the r-step kernel, over all pairs, at r = 8.
- exit channel: for 240 starting neurons drawn with numpy default_rng(20261114), the total variation
  distance between the laws of (first exit step, exit neuron) for every pair; steps are summed until
  the largest remaining mass is below 1e-4 or 200 steps. Reported: mean, median, largest, share of
  pairs below 0.5.
- reach ratio from A to B at one step: the probability that a walker leaving A (stationary exit flow)
  lands in B, divided by the probability that a stationary walker outside B enters B in one step.

## 6. Hypotheses and criteria (primary analysis)

H-EB, the ellipsoid-body set is a compartment. Confirmed only if C1 to C4 all hold.
- C1 closure: loss per step of EB < 0.20, and smaller than that of MB and smaller than that of AL.
- C2 mixing before leaving: stays / mixes > 1.3 for EB, and < 1 for MB and for AL.
- C3 forgetting inside: Dobrushin coefficient of the conditioned walk in EB at r = 8 < 0.6.
- C4 exit channel: mean distance for EB < 0.37; share of pairs below 0.5 > 0.70; mean for EB at least
  0.10 below the mean for MB.

H-PIT, central complex and mushroom body are preferentially coupled (secondary; judged on its own).
- C5: the one-step reach ratio from MB to CX is at least 3 times the one-step reach ratio from AL to CX.

Exploratory values behind the thresholds (male >= 1, male >= 5, FlyWire >= 1, FlyWire >= 5):
- EB loss per step 0.116, 0.102, 0.094, 0.074 (MB 0.367, 0.612, 0.413, 0.661; N2 0.99).
- EB stays / mixes 1.92, 1.70, 2.80, 2.49 (MB 0.21, 0.44, 0.40, 0.64; AL 0.20, 0.19, 0.14, 0.12; N2 0.03 to 0.08).
- EB Dobrushin coefficient at r = 8: 0.255, 0.364, 0.210, 0.381 (N2 1.00).
- EB exit mean 0.267, 0.235, 0.231, 0.197; pairs below 0.5: 0.948, 0.977, 0.966, 0.991 (MB mean 0.499, 0.647, 0.539, 0.705).
- C5 contrast 9.2, 10.8, 13.4, 17.2 (N2 1.2; N1 1.1).
Rule used: geometric midpoint between the weakest of the four exploratory values and the nearest value
that marks no compartment (a comparison set, N2, or the natural boundary 1), rounded; the 0.10 gap in C4
is half the smallest gap seen.

## 7. Reported, not scored

Amine shares of the flow into, inside and out of EB and MB (by the type-name list of the crosswalk and by
the data set's own transmitter prediction); the loss per step of the union of CX and MB against CX alone;
the largest exit distance; the Dobrushin coefficient at r = 4 and 16; the reverse reach from CX to MB;
the exit channels of MB and CX; every number at count >= 5.
The amine statistic is not scored because transmitter predictions do not carry across data sets: the
classifier labels nearly all Kenyon cells as dopamine in FlyWire (5,172 of 5,177) and in the raw male
predictions (4,058 of 4,064), and only the male consensus corrects this to acetylcholine.

## 8. Protocol

- Environment: Python 3.13, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, pyarrow 26.0.1; two BLAS threads.
- Command, run once:
  python3 transfer_screen.py banc_888_id banc_888_meta.feather banc_888_edgelist_simple_v3.feather results_banc.json 1 5
- The three code files are fixed by the sha256 values in MANIFEST.sha256. The script writes their
  hashes, and the hashes of both inputs, into the result file.
- No rerun after a result exists. If the script stops with a technical error before any statistic is
  printed, only the technical fault is fixed, the change is logged as a deviation with a diff, and the
  script is run again. Definitions, sets, thresholds and seeds are never changed.
- Everything is reported whatever the outcome. A failed criterion is reported as failed; later
  analyses of why are exploratory and labelled as such.

## 9. Limitations known before the run

- In BANC the visual input types of the EB set are less completely typed (MeTu2: 27 neurons against 101
  in the male and in FlyWire), and 28,099 nodes have no cell type. The rule is applied as it stands.
- Synapse detection and proofreading differ from the male data set, so counts are not comparable edge by
  edge; this is why the primary analysis uses all edges and the thresholds leave a margin.
- Because N3 reproduces the statistics, this run tests whether a type-level architecture recurs in
  another animal, not the arrangement of individual neurons.
- Thresholds were partly calibrated on FlyWire, a female brain; BANC is a female brain and nerve cord.
- A random walk on synapse counts is a structural model. A step is not a unit of time, and nothing here
  concerns neural activity or experience.
