# Post hoc exploratory analysis: why the BANC result weakens with the synapse threshold

ASCII, no em dashes. Written on 14 Nov 2026, AFTER the confirmatory run (registration
doi:10.5281/zenodo.24223802, results doi:10.5281/zenodo.24224168).

**Status: exploratory. Nothing here is confirmatory, nothing here was registered, and no criterion is scored.**
The registered verdict (primary analysis, all edges) is not changed by this analysis. The confirmatory
result files in `results/` were not touched and the registered analysis was not rerun; the BANC rows at
thresholds 1 and 5 below are recomputed by the sweep and equal the registered numbers.

## Question

In BANC the registered criteria C2 (stays / mixes > 1.3) and C3 (Dobrushin coefficient at r = 8 < 0.6) hold
with all edges and would not hold at 5 synapses. In the male CNS and in FlyWire they hold at both. Why?

## Data and code

- BANC release 888, male CNS v1.0, FlyWire v783: the same files as in the registration and in `exploratory/`.
- Statistics: the frozen `prereg/screen_lib.py`, imported unchanged (sha256 e062d134...). Sets: the registered
  names-only definitions (`prereg/crosswalk_types.json`).
- Scripts `p1` to `p6` and their result files and logs are in this folder. `p6_summary.py` prints every number
  quoted below from the result files (`results/p6_summary.txt`).
- Check of the pipeline: the sweep reproduces the earlier values at thresholds 1 and 5 in all three data sets
  (BANC 1.396 and 0.842; male 1.919 and 1.697; FlyWire 2.801 and 2.490).

## Findings

### 1. The weakening is not specific to 5 synapses: in BANC it starts at 2 (p1)

EB set, three data sets (BANC / male / FlyWire):

| threshold | stays / mixes | Dobrushin, r = 8 | exit mean | giant component, share of nodes |
|---|---|---|---|---|
| 1 | 1.40 / 1.92 / 2.80 | 0.41 / 0.26 / 0.21 | 0.31 / 0.27 / 0.23 | 0.98 / 0.99 / 0.97 |
| 2 | 1.27 / 1.89 / 2.82 | 0.57 / 0.28 / 0.24 | 0.30 / 0.25 / 0.22 | 0.96 / 0.98 / 0.94 |
| 3 | 1.10 / 1.82 / 2.74 | 0.69 / 0.32 / 0.28 | 0.32 / 0.24 / 0.19 | 0.91 / 0.97 / 0.92 |
| 4 | 0.87 / 1.76 / 2.59 | 0.84 / 0.34 / 0.34 | 0.34 / 0.24 / 0.20 | 0.85 / 0.96 / 0.91 |
| 5 | 0.84 / 1.70 / 2.49 | 0.90 / 0.36 / 0.38 | 0.36 / 0.24 / 0.20 | 0.78 / 0.95 / 0.89 |
| 8 | 0.33 / 1.57 / 2.47 | 1.00 / 0.50 / 0.49 | 0.45 / 0.24 / 0.20 | 0.54 / 0.91 / 0.83 |
| 10 | 0.65 / 1.49 / 2.23 | 1.00 / 0.57 / 0.61 | 0.49 / 0.23 / 0.21 | 0.39 / 0.88 / 0.77 |
| 12 | (degenerate) / 1.42 / 1.51 | 1.00 / 0.64 / 0.76 | 0.38 / 0.24 / 0.21 | 0.28 / 0.83 / 0.71 |
| 15 | (degenerate) / 1.05 / 0.85 | 0.98 / 0.74 / 0.91 | 0.37 / 0.25 / 0.23 | 0.19 / 0.74 / 0.60 |

- First threshold at which stays / mixes is at or below 1.3: BANC 2, male 15, FlyWire 15.
- First threshold at which the Dobrushin coefficient at r = 8 is at or above 0.6: BANC 3, male 12, FlyWire 10.
- So the same statistics weaken with the threshold in all three animals; BANC reaches that point at a
  threshold three to seven times lower. **In BANC, C2 holds only when edges of a single synapse are counted.**
- From threshold 12 the BANC graph itself has fallen apart (giant component below 30 percent of the nodes)
  and the numbers are not interpretable.
- The loss per step (C1) and the exit channel (C4) stay inside the registered limits up to threshold 5 in
  BANC (loss 0.13 to 0.12; exit mean 0.31 to 0.36 against the limit 0.37; pairs below 0.5: 0.87 to 0.78).

### 2. BANC has far fewer synapses, in the whole graph and inside EB (p2)

| | BANC | male | FlyWire |
|---|---|---|---|
| synapses per neuron, whole graph | 268 | 745 | 495 |
| synapses inside EB per EB neuron | 434 | 1,213 | 669 |
| synapses per edge inside EB, mean | 6.4 | 12.3 | 9.7 |
| synapses per edge inside EB, median | 3 | 5 | 4 |
| share of edges inside EB with 1 synapse | 0.33 | 0.20 | 0.24 |
| share of edges inside EB below 5 synapses | 0.66 | 0.46 | 0.54 |
| share of EB synapses carried by edges below 5 | 0.19 | 0.08 | 0.11 |
| EB neurons whose strongest input inside EB is below 5 | 57 | 3 | 11 |

BANC over male: 0.36 (synapses per neuron, whole graph and inside EB). BANC over FlyWire: 0.54 (whole graph),
0.65 (inside EB). 579 of the 595 EB neurons of BANC carry the proofread flag.

### 3. The lower density alone does not explain it (p4, uniform thinning)

Every synapse of the male (or FlyWire) data kept independently with the probability that brings the data to
the BANC density, then the same statistics. Mean over three seeds; the ranges are in `results/p6_summary.txt`.

| threshold | BANC | male, thinned uniformly | FlyWire, thinned uniformly |
|---|---|---|---|
| stays / mixes, 1 | 1.40 | 1.92 | 3.29 |
| stays / mixes, 2 | 1.27 | 1.81 | 3.34 |
| stays / mixes, 3 | 1.10 | 1.67 | 3.28 |
| stays / mixes, 5 | 0.84 | 1.51 | 3.34 |
| Dobrushin r = 8, 1 | 0.41 | 0.27 | 0.20 |
| Dobrushin r = 8, 5 | 0.90 | 0.59 | 0.42 |

At the BANC density the male would still pass C2 at 5 synapses (1.51) and sit at the limit of C3 (0.59; one
seed of three at 0.60). FlyWire would pass both with room. Density moves the male towards the thresholds but
does not reproduce BANC, and it changes nothing at threshold 1.

### 4. The pattern of the loss explains it (p4, block-pattern thinning)

Inside EB the BANC : source ratio of synapses per possible pair differs strongly between family-to-family
blocks (0.05 to 0.72 against the male, among the 42 blocks with at least 2,000 synapses in the male).
Thinning each block by its own ratio:

| threshold | BANC | male, BANC pattern | FlyWire, BANC pattern |
|---|---|---|---|
| stays / mixes, 1 | 1.40 | 1.22 | 1.89 |
| stays / mixes, 2 | 1.27 | 0.88 | 1.65 |
| stays / mixes, 3 | 1.10 | 0.70 | 1.26 |
| stays / mixes, 5 | 0.84 | 0.30 | 0.58 |
| Dobrushin r = 8, 1 | 0.41 | 0.47 | 0.38 |
| Dobrushin r = 8, 2 | 0.57 | 0.62 | 0.50 |
| Dobrushin r = 8, 3 | 0.69 | 0.77 | 0.65 |
| Dobrushin r = 8, 5 | 0.90 | 0.98 | 0.87 |

With the BANC pattern both other animals show the BANC behaviour: a monotone fall with the threshold, and
failure of both criteria by 5 synapses. The pattern is sufficient, and in the male more than sufficient (it
overshoots BANC), so the family level is a coarse description of what differs.

### 5. Which connections: ring neurons to and from the rest of the set (p5)

Edges inside EB in five groups; share of the BANC : source synapse ratio (uniform would be 0.36 for the male,
0.65 for FlyWire):

| group | against male | against FlyWire |
|---|---|---|
| columnar and extrinsic neurons onto ring neurons (EL, ExR, EPG, PEN, PEG onto ER) | 0.18 | 0.47 |
| ring neurons onto columnar and extrinsic neurons (ER onto EL, ExR, EPG, ...) | 0.22 | 0.43 |
| ring neurons onto ring neurons | 0.37 | 0.67 |
| among columnar and extrinsic neurons | 0.45 | 0.78 |
| visual input tail (TuBu, MeTu2, TuTuB, AOTU046) | 0.27 | 0.31 |

Examples, synapses per connected pair (BANC / male / FlyWire): ER onto EL 2.7 / 8.4 / 7.6; ER onto ExR
4.2 / 18.6 / 16.7; ER onto EPG 3.7 / 10.7 / 5.9; EL onto ER 1.4 / 5.2 / 2.8; ExR onto ER 2.0 / 12.1 / 6.2.

Attribution (one seed; BANC pattern in one group only, uniform thinning elsewhere; stays / mixes and Dobrushin
at r = 8 at thresholds 1, 2, 5):

| variant | male source | FlyWire source |
|---|---|---|
| uniform | 1.93, 1.82, 1.53; 0.26, 0.33, 0.58 | 3.27, 3.32, 3.39; 0.20, 0.25, 0.42 |
| only ring neurons onto the rest | 1.39, 1.17, 1.32; 0.46, 0.59, 0.88 | 2.33, 2.20, 1.71; 0.35, 0.45, 0.75 |
| only the rest onto ring neurons | 1.40, 1.18, 0.64; 0.26, 0.35, 0.66 | 2.91, 2.85, 2.61; 0.19, 0.24, 0.40 |
| only ring onto ring | 1.87, 1.74, 1.48; 0.27, 0.35, 0.62 | 3.27, 3.34, 3.38; 0.22, 0.26, 0.45 |
| only columnar | 2.30, 2.27, 2.17; 0.26, 0.32, 0.58 | 3.01, 3.10, 3.29; 0.23, 0.28, 0.45 |
| only the visual tail | 1.91, 1.78, 1.50; 0.26, 0.33, 0.59 | 3.33, 3.42, 3.39; 0.19, 0.24, 0.43 |
| full pattern | 1.18, 0.81, 0.34; 0.48, 0.64, 0.98 | 1.88, 1.59, 0.48; 0.39, 0.52, 0.89 |
| full pattern except ring neurons onto the rest | 1.82, 1.65, 1.17; 0.24, 0.30, 0.62 | 2.78, 2.76, 2.40; 0.21, 0.24, 0.43 |

- In both source animals the output of the ring neurons onto the rest of the set is the group that matters:
  weakening it alone raises the Dobrushin coefficient to or near the BANC levels, and restoring it alone
  removes most of the effect of the full pattern.
- In the male the input to the ring neurons from the rest also matters (it alone brings stays / mixes to
  1.40, 1.18, 0.64). In FlyWire that group has a small effect.
- The visual tail, where BANC is least completely typed (MeTu2: 27 neurons against 101), is not the cause.
  Neither is ring onto ring, nor the columnar group. A diagnostic subset without the tail gives the BANC
  numbers almost unchanged (stays / mixes 1.36, 1.26, 1.10, 0.85 at thresholds 1, 2, 3, 5).

### 6. The slow mode is the same in all three animals (p3)

Second eigenvalue of the walk conditioned to stay in EB, and what it separates:

| | threshold | second eigenvalue | stationary mass on ER3a + ER3d + ER3m |
|---|---|---|---|
| BANC | 1 | 0.854 | 0.42 |
| BANC | 2 | 0.884 | 0.58 |
| BANC | 3 | 0.917 | 0.73 |
| BANC | 5 | 0.957 | 0.84 (ER3d alone 0.68) |
| male | 1 | 0.793 | 0.21 |
| male | 5 | 0.834 | 0.25 |
| male | 10 | 0.882 | 0.39 |
| FlyWire | 1 | 0.766 | 0.17 |
| FlyWire | 5 | 0.829 | 0.14 |

In every data set the slow mode has the ring-neuron types ER3a, ER3d and ER3m (with their inputs TuBu02 and
TuBu03) at one end and the compass neurons (EPG, PEN, PEG) with ER2, ER4d and ER3w at the other. BANC with
all edges (0.854) sits where the male sits at 7 synapses (0.855, from the mixing time of the sweep). With rising threshold the BANC walk is caught
in the ER3 group: 84 percent of the stationary mass at 5 synapses.

### 7. The failure of C3 at 5 synapses is global, not a few neurons (p3)

At 5 synapses in BANC, 515 of 550 starting neurons are farther than 0.3 (total variation) from the stationary
law after 8 steps, and 33 percent of all pairs are farther than 0.6 from each other. Removing the starts one
by one, worst first, the coefficient stays above 0.6 until more than 500 of the 550 are removed. With all
edges no start is farther than 0.3 and no pair above 0.6.

## Reading

- The registered result stands as registered. Its limit is sharper than RESULTS.md says: in BANC, C2 needs
  the single-synapse edges, the least reliable edges of any connectome.
- The cause is not that BANC has fewer synapses overall. It is that the two-way coupling between the ring
  neurons and the rest of the set is carried by few synapses per pair in BANC (1 to 4 in the blocks listed
  above, against 4 to 19 in the male), so any threshold removes it first, and the set separates along a slow mode that all three
  animals share.
- **Not decided by these data:** whether that weak coupling in BANC is biology (an individual difference)
  or measurement (synapse detection with a recall that differs between connection types). Two of the three
  animals agree with each other on the ring-neuron output weights (ER onto EL and ER onto ExR are close in
  the male and in FlyWire, and three to four times lower in BANC), which points to something specific to
  BANC, but does not say which. The synapse-detection validation of BANC has not been read for this note.
- Positive result of the look: the slow mode of the compartment (ER3a, ER3d, ER3m against the compass side)
  recurs in three animals. It was not registered and would need its own test.

## Limits of this analysis

- Post hoc, after the outcome was known; many analysis choices; no correction for them. Hypothesis
  generating only.
- Thinning is a model of lower synapse recall (independent loss of synapses). Real detection errors need
  not be independent.
- The block pattern is defined at the level of twelve families; it overshoots BANC in the male.
- The attribution uses one seed. The uniform and full-pattern variants use three seeds (ranges in the
  summary file).
- In p4 and p5 the stationary distribution is iterated 300 and 60 times; it enters only the loss per step,
  which is therefore not quoted from p5.
- Hemibrain was not opened.

## Suggestions for a later registration (not acted on)

- State thresholds in units that carry across data sets (for example a multiple of the median synapse
  count of the edges inside the set), not as absolute counts.
- Register the threshold at which a criterion first fails as a robustness index.
- If the slow mode is to be tested, register the two groups of types before opening a fourth data set.
