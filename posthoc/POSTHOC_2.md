# Post hoc exploratory analysis, part 2: measurement or biology?

ASCII, no em dashes. Written on 14 Nov 2026, after the confirmatory run and after part 1 (`POSTHOC.md`).

**Status: exploratory. Nothing here is confirmatory, nothing here was registered, and no criterion is scored.**
The registered verdict is not changed. The confirmatory files in `results/` were not touched, the registered
analysis was not rerun, and hemibrain was not opened.

## Question

Part 1 found that the BANC result weakens with the synapse threshold because the connections between the ring
neurons and the rest of the EB set carry few synapses in BANC. Is that a property of the animal or of the
measurement?

## Short answer

- **The low synapse density of the BANC edge list is measurement.** BANC detects more synapses per neuron than
  the other two data sets, but only 22 percent of them end on an identified neuron.
- **Several technical causes are excluded** for the ring-neuron pattern: the detector version and size cutoff,
  missing sections, damaged regions, and unattached presynaptic arbors.
- **Part of the pattern is measurement with high likelihood:** the connections of the aminergic neurons of the set
  (EL, ExR2, ExR3). Aminergic neurons have few detected output synapses in BANC across the whole brain.
- **The core of the pattern is not decided:** BANC finds unusually few synapses made by ring neurons (lowest 5 to 9
  percent of 3,876 cell types against both other animals). It is not a uniform regional effect. No ground truth
  exists for the detector in the ellipsoid body. Under-detection of this class of synapse and a real difference
  in this animal both remain possible.
- What would decide it: a dense ground-truth cutout in the ellipsoid body of BANC.

## Sources read

- BANC data bucket, `README.txt`, `CHANGELOG.md` and the notes in `documentation/` (edge lists v2 and v3, synapse
  tables, dense evaluation, reviewed samples, problem regions, backbone proofreading, metrics). Base:
  https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/
- BANC Dataverse deposit, doi:10.7910/DVN/7WTH1N (version 3.0): `251013_synapse_evaluation.csv`,
  `2026-05-14_aelysia_synapse_sample_complete_v3.csv`, `2024-09-20_aelysia_synapse_sample_complete_v2.csv`,
  `banc_problem_regions.csv` (md5 checked against the deposit).
- BANC paper: Bates, Phelps et al., Distributed control circuits across a brain-and-cord connectome, Nature 656
  (27 August 2026), doi:10.1038/s41586-026-10735-w. Checked on 15 Nov 2026 in the publisher's PDF: the passages on
  synapse detection, attachment, thresholds and the central complex were searched and read; the paper was not read
  end to end.
- BANC preprint, Methods, same tool, https://pmc.ncbi.nlm.nih.gov/articles/PMC12324551.

What they say:
- The registered edge list (`edgelist_simple_v3`) is built from the v3 detector with synapses of at least 10
  voxels. A second list (`v2`, at least 5 voxels, an earlier detector) is the one the BANC paper uses.
- Dense evaluation of the detector (v2) in 16 cutouts: pooled precision 0.87, recall 0.78. Recall per cutout
  runs from 0.53 (mushroom-body medial lobe) to 0.91 (SMP). The one central-complex cutout (fan-shaped body)
  has recall 0.72. **There is no cutout in the ellipsoid body.**
- Reviewed samples of single detections in the ellipsoid body: 31 of 32 true (v2); 33 of 40 true (v3, a sample
  stratified by detector score, so not a population precision). These measure false positives, not misses.
- BANC paper (2026): precision 0.87, recall 0.78, F-score 0.83. 74 percent of presynaptic ends and 23 percent
  of postsynaptic ends are connected to a proofread neuron (preprint: 65 and 22). 18 percent of the synaptic
  links have an identified neuron on both sides, against 41.9 percent in FAFB, 35.3 in the hemibrain, 44.4 in
  MANC and 40.1 within the neuropils of maleCNS (the paper's own comparison).
- The analyses of the paper use the v2 synapse table with postsynapse size at least 5. On thresholds for
  connections the Methods say that no synapse threshold was applied in the connectivity analyses, except in the
  influence analysis, where connections of fewer than five synapses were removed. Betweenness uses edges of at
  least five synapses, the links between networks and the counts of direct partners of Fig. 6 use at least 5,
  and the circuit diagrams show connections of at least 10.
- The authors write that, because the detector relies on postsynaptic profiles, they may under-detect some
  classes of synaptic connection with less distinct postsynaptic sites, and that Kenyon cells have 36.6
  outgoing links on average against around 200 in FAFB. The preprint calls the mushroom body a known problem
  area for the detection method.
- The paper divides the nervous system into 13 networks by spectral clustering. Of these, the central complex
  network and the olfactory networks have low influence on effector neurons, weak input from other networks and
  few ascending or descending neurons; the authors call them supervisory. This is their own analysis, made with
  the thresholds above, and it is not a test of the criteria used here.
- Neurons are proofread to the backbone standard, which the notes call a lower bar than full proofreading.

## Findings from the data

### 1. BANC detects more and attaches less (p7, p8, p12)

One pass over the BANC v3 synapse table (259,393,451 rows; after dropping autapses 255,774,906 synapses). Check:
the counts between graph neurons equal the registered edge list exactly (12,305,034 pairs, 39,297,815 synapses).

| | BANC | male | FlyWire |
|---|---|---|---|
| synapses made per neuron, identified partners only (the edge lists) | 268 | 745 | 495 |
| synapses made per neuron, all partners | 1,205 | | |

- Of all BANC synapses, 77.7 percent have their postsynaptic end on a segment nobody has identified, and 28.0
  percent their presynaptic end. Of the synapses a graph neuron makes, 22 percent end on a graph neuron (median
  per neuron 0.22). The BANC paper gives 23 percent (preprint 22). Both ends lie on graph neurons for 15.4
  percent of the synapses counted here (39,297,815 of 255,774,906); the paper gives 18 percent for identified
  neurons on both sides, and about 40 percent for the other data sets.
- Over 3,876 cell types present in all three data sets, the synapses made per neuron in BANC (all partners) are
  1.79 times the male (quartiles 1.26 to 2.34) and 2.45 times FlyWire (1.86 to 3.30). The edge-list counts are
  0.32 times the male.
- So the factor 0.36 between the BANC and the male edge lists is attachment of postsynaptic twigs, not detection.

### 2. The EB set is attached better than average (p8, p10)

- Share of the synapses made by EB neurons that end on a graph neuron: 0.39 (ring neurons 0.45, EPG 0.40, EL
  0.50, ExR 0.44). For the ring, columnar and extrinsic families the medians rank between the 88th and the
  100th percentile of all graph neurons. The visual input families are lower (TuBu 0.12, MeTu2 0.16, TuTuB 0.19,
  AOTU046 0.30).
- In the ellipsoid-body neuropil 90 percent of presynaptic ends and 51 percent of postsynaptic ends are on graph
  neurons; over the 24 neuropils compared, the postsynaptic share runs from 0.17 to 0.53.
- Inputs to EB neurons come from graph neurons in 78 percent of cases (whole graph 72). The presynaptic arbors
  of the set are therefore not missing.

### 3. Excluded technical causes

- **Detector version and size cutoff (p1 on the v2 list).** With the v2 edge list (the table of the BANC paper)
  stays / mixes is 1.42, 1.25, 1.08, 0.81, 0.99 at thresholds 1 to 5 and the Dobrushin coefficient at r = 8 is
  0.45, 0.67, 0.84, 0.93, 0.98. The registered list gives 1.40, 1.27, 1.10, 0.87, 0.84 and 0.41, 0.57, 0.69,
  0.84, 0.90. Same behaviour. With all edges the v2 list gives loss 0.128, exit mean 0.30 and 87 percent of
  pairs below 0.5, as the registered list does. On the v2 list ring neurons make 285 synapses per EL neuron and
  694 per ExR1 neuron (registered list 392 and 925; FlyWire 1,562 and 4,930; male 2,130 and 7,128) (p13).
  (Detections below 10 voxels of the v3 detector are only in a 15 GB raw file that was not used.)
- **Missing sections (p8).** Along the sectioning axis, 3 of the 1,418 sections through the set fall below 20
  percent of the local median; filling them would add 0.01 percent of the synapses.
- **Damaged regions (p8).** 5.5 percent of the synapses inside the EB set lie in the 22 documented problem
  regions (dorsal esophageal crush 12,323; dorsal wavy patch 1,742; of 257,977).

### 4. Aminergic synapses depend on the detector across the whole brain (p9)

Synapses made per neuron, ratio over cell types grouped by the transmitter of the type (male consensus):

| transmitter (types) | BANC : male | BANC : FlyWire | FlyWire : male |
|---|---|---|---|
| all (3,876) | 1.79 | 2.45 | 0.71 |
| acetylcholine (2,174) | 1.85 | 2.67 | 0.67 |
| GABA (659) | 1.90 | 2.26 | 0.84 |
| glutamate (878) | 1.67 | 2.26 | 0.74 |
| dopamine (35) | 0.31 | 0.86 | 0.36 |
| serotonin (7) | 0.48 | 1.93 | 0.27 |
| octopamine (13) | 0.49 | 1.34 | 0.31 |

- The three detectors disagree most on aminergic neurons, in the order male, FlyWire, BANC. This is a
  property of the detectors, not of the ellipsoid body.
- EL (octopamine in the BANC annotations) has the 7th lowest ratio of the 3,876 types against the male (0.12)
  and the 14th lowest against FlyWire (0.43). ExR2 is dopaminergic and ExR3 serotonergic.
- GABA types as a group are ordinary (1.90), so the ring neurons are not low because they are GABAergic.

### 5. The central complex is the lowest region, the ellipsoid body lowest of all (p10)

Cell types grouped by the neuropil that holds most of their output synapses (24 neuropils with at least 8 types):

| neuropil (types) | BANC : male (rank among all types) | BANC : FlyWire (rank) |
|---|---|---|
| ellipsoid body (13) | 0.76 (0.07) | 1.32 (0.10) |
| fan-shaped body (114) | 0.94 (0.13) | 1.43 (0.12) |
| antennal lobe (127) | 0.95 (0.13) | 1.89 (0.26) |
| AOTU (22) | 1.08 (0.18) | 1.70 (0.20) |
| the other 20 | 1.18 to 2.61 | 1.49 to 11.8 |

- Groups by name: ring neurons 0.72 (rank 0.06) and 1.31 (0.09); fan-shaped body tangential neurons 0.94 and
  1.37; hDelta and vDelta 0.90 and 1.13; Kenyon cells 0.40 (0.01) and 1.35; dopaminergic PAM and PPL 0.28 and
  0.70 (0.01).
- Recall of the detector in the densely evaluated cutouts against the ratio of the types of the same neuropil:
  rank correlation 0.60 (male) and 0.68 (FlyWire) over 9 cutouts. Suggestive only. The fan-shaped body has the
  second lowest recall of the central brain (0.72) and the lowest ratio of the nine.

### 6. But it is not a uniform regional effect (p11)

Synapses inside the EB set by neuropil: 72 percent lie in the ellipsoid-body neuropil, where BANC has 0.56 of the
FlyWire count. Thinning FlyWire regionally (every synapse of the set kept with the ratio of its neuropil, no
family information) leaves the set closed:

| threshold | BANC | FlyWire, regional thinning | FlyWire, family-block pattern (part 1) |
|---|---|---|---|
| stays / mixes, 1 | 1.40 | 3.08 | 1.89 |
| stays / mixes, 2 | 1.27 | 3.12 | 1.65 |
| stays / mixes, 3 | 1.10 | 3.01 | 1.26 |
| stays / mixes, 5 | 0.84 | 3.06 | 0.58 |
| Dobrushin r = 8, 5 | 0.90 | 0.45 | 0.87 |

So a lower yield in the region does not produce the weakening. The loss has to differ between connection types
inside the same neuropil.

### 7. Which connections inside the neuropil: those of the ring neurons (p9, p12)

- The BANC : source ratio of a family block is largely a factor of the presynaptic family: it explains 67
  percent (male) and 57 percent (FlyWire) of the variance of the log ratios; presynaptic and postsynaptic
  factors together explain 78 and 81 percent. The largest remainders are ExR onto ring neurons and ring neurons
  onto ExR.
- For the same postsynaptic cells, the synapses from ring neurons are reduced far more than those from the
  columnar neurons (EPG, PEN, PEG). Ratio BANC : male and BANC : FlyWire:

| onto | from ring neurons | from columnar neurons |
|---|---|---|
| EL | 0.18, 0.25 | 0.39, 0.42 |
| ExR1 | 0.13, 0.19 | 0.64, 1.28 |
| ExR2 | 0.21, 0.61 | 0.43, 0.73 |
| ExR4 | 0.19, 0.61 | 0.50, 1.19 |
| ExR6 | 0.27, 0.58 | 0.62, 0.86 |

- Counting all partners, ring neurons make 1,123 synapses per neuron in BANC against 1,715 in the male and 914
  in FlyWire: ratios 0.66 and 1.23, where the typical cell type has 1.79 and 2.45.
- The glutamatergic ExR4 and ExR6 connect to the columnar neurons at ordinary ratios (0.50 to 0.62 against the
  male, 0.82 to 1.82 against FlyWire).

## Reading

| candidate cause | status |
|---|---|
| lower synapse density of the BANC edge list | measurement (attachment of postsynaptic twigs); established |
| detector version, size cutoff | excluded (the v2 list behaves the same) |
| missing sections, damaged regions | excluded as main cause |
| unattached presynaptic arbors of the set | excluded |
| GABA as such | excluded (GABA types are ordinary) |
| a uniformly lower yield in the ellipsoid-body neuropil | present (0.56 of FlyWire), but does not produce the weakening |
| connections of the aminergic neurons (EL, ExR2, ExR3) | measurement likely (brain-wide detector dependence) |
| few detected synapses of ring neurons | **not decided**: under-detection of this synapse class, or the animal |

For the undecided row, what points to measurement: the two other animals, measured with two other detectors,
agree with each other (ring neurons onto EL: 2,130 and 1,562 synapses per EL neuron; BANC 392); the shortfall
covers the whole central complex, fan-shaped body included; the only central-complex cutout has below-average
recall; the authors expect under-detection for some synapse classes. What keeps it open: no ground truth in the
ellipsoid body, and the pattern is specific to one cell class inside one neuropil.

## Consequences

- For the paper: the confirmation in BANC holds with all edges (and with either synapse set of BANC). Robustness
  to thresholds cannot be assessed in BANC, because the ring-neuron connections rest on few detected synapses
  there. That should be stated as a limit of the data set, with the facts above, and not as biology.
- The BANC paper applies no synapse threshold in its connectivity analyses in general. Its influence analysis,
  betweenness, links between networks and counts of direct partners use connections of at least 5 synapses, and
  its circuit diagrams connections of at least 10. That is the range in which criteria C2 and C3 do not hold in
  BANC, so a reader of that paper may ask for the thresholded result, and the limit belongs in the main text.
  (Corrected on 15 Nov 2026. The earlier wording said that the paper itself works with connections of at least
  5 synapses, which overstated it: the published Methods apply the threshold to some analyses, not to all.)
- What would decide the open row: a densely annotated cutout in the ellipsoid body of BANC (recall of the
  detector for ring-neuron synapses). The BANC team may have one or could say; the project takes questions at
  https://github.com/jasper-tms/the-BANC-fly-connectome/issues.
- For a later registration on any connectome: fix the synapse set and its cutoffs by name, and state thresholds
  in units that carry across data sets.

## Limits of this analysis

- Post hoc, many choices, no correction for them.
- Ratios between data sets mix detector, proofreading and animal. Only differences between cell types or
  regions within the same pair of data sets are interpreted.
- Cell types are matched by identical type names; transmitters are the male consensus of the type.
- The FlyWire per-synapse table holds 79 percent of the synapses of the FlyWire edge list, and its neuropil
  labels outside the ellipsoid body do not match those of BANC (noduli, AOTU), so only the ellipsoid-body ratio
  of p11 is reliable. The male files used here have no neuropil labels.
- Statements of the BANC paper were checked in the publisher's PDF (the passages named under Sources, not the
  whole paper); the preprint was read through an extraction tool. The attachment percentages agree with the
  counts made here.
