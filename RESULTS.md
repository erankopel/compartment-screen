# Confirmatory results, 14 Nov 2026

ASCII, no em dashes.

- Registration: release `v0.1.0-prereg`, doi:10.5281/zenodo.24223802 (concept doi:10.5281/zenodo.24223801),
  record created 2026-11-14T15:59:33Z.
- Run: one run, 16:03:09 to 16:07:21 UTC on 14 Nov 2026, after the registration. Completed without error
  (log line saved results_banc.json; the exit code line in the log is not informative, see
  DEVIATIONS.md), no rerun. Log: `results/run.log`. Result file: `results/results_banc.json`
  (sha256 5ffee1df970754889e7f7309b56cf91e9a699d5b4b2fc8a01ec2d0e59210f364).
- Code: the three files of `prereg/` taken from the Zenodo archive of the registration; their sha256
  equal `prereg/MANIFEST.sha256` and are written into the result file.
- Inputs: BANC release 888, `banc_888_meta.feather` (md5 YjyrACb8v8r2J/2q+en7HQ==) and
  `banc_888_edgelist_simple_v3.feather` (md5 CFQrB3HbdBjtR0vmDcmIbA==, sha256
  8c296e946f3c69a8c7222f30ad75fa8a98eeb189124fec6df829c9125f4be64b). Both md5 values equal the
  registration and were checked before the run. The edge list was downloaded at 16:03 UTC, after the
  registration.
- Environment: Python 3.13.18, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, pyarrow 26.0.1, two BLAS threads.
- Graph: 146,769 nodes; 12,305,034 edges and 39,297,815 synapses inside it; giant strongly connected
  component 144,122 nodes at count >= 1 (stationary distribution in 477 iterations, last change 9.5e-13).
  Sets in the component: EB 589 of 595, CX 3,065, MB 4,898, AL 1,128.

## Verdicts against the registered criteria (primary analysis, all edges)

| | criterion (PREREG Sec. 6) | value | met |
|---|---|---|---|
| C1 | loss per step of EB < 0.20 | 0.1333 | yes |
| C1 | and smaller than MB and than AL | MB 0.3111, AL 0.4225 | yes |
| C2 | stays / mixes of EB > 1.3 | 1.3958 (stays 9.57, mixes 6.86) | yes |
| C2 | and < 1 for MB and for AL | MB 0.3547, AL 0.1903 | yes |
| C3 | Dobrushin coefficient of the conditioned walk in EB at r = 8 < 0.6 | 0.4133 (4 neurons dropped) | yes |
| C4 | exit channel of EB, mean < 0.37 | 0.3067 | yes |
| C4 | share of pairs below 0.5 > 0.70 | 0.8667 | yes |
| C4 | mean at least 0.10 below MB | MB 0.5816, gap 0.2749 | yes |
| **H-EB** | C1 to C4 all hold | | **CONFIRMED** |
| C5 | one-step reach ratio MB to CX at least 3 times that of AL to CX | 11.45 against 0.95, contrast 12.06 | yes |
| **H-PIT** | C5 (secondary, judged on its own) | | **CONFIRMED** |

## Reported, not scored (PREREG Sec. 7)

At count >= 1:
- Loss per step of the union of CX and MB: 0.1903, against 0.1660 for CX alone.
- CX as a whole: stays / mixes 0.7856.
- Dobrushin coefficient in EB at r = 4: 0.7532; at r = 16: 0.1215.
- Exit channel of EB: median 0.2559, largest 1.0000 (89 steps, mass left below 1e-4).
  MB: mean 0.5816, pairs below 0.5: 0.3331. CX: mean 0.4439, pairs below 0.5: 0.6657.
- Reach MB to CX within 1 to 4 steps: 0.0396, 0.0719, 0.0989, 0.1219 (stationary walker 0.0035, 0.0067,
  0.0099, 0.0130); ratio at 4 steps 9.38. AL to CX: ratio 0.95 at 1 step, 2.41 at 4 steps.
  Reverse, CX to MB: ratio 6.74 at 1 step, 6.38 at 4 steps.
- Amine share of the flow, type-name list: whole chain 0.0021; EB into 0.0060, inside 0.0788, out 0.1439;
  MB inside 0.1376. By the data set's own transmitter prediction (11,428 neurons labelled aminergic):
  whole chain 0.0472; EB into 0.2165, inside 0.0898, out 0.1554.

At count >= 5 (giant component 113,971 nodes; EB 555, CX 2,602, MB 2,147, AL 1,068 in it):
- EB: loss per step 0.1246; stays / mixes 0.8415 (stays 19.36, mixes 23.01); Dobrushin coefficient at
  r = 8 0.8954 (r = 4: 0.9788, r = 16: 0.7257; 5 neurons dropped); exit channel mean 0.3605, pairs
  below 0.5: 0.7822, largest 1.0000.
- MB: loss 0.7094, stays / mixes 2.8602, exit mean 0.7188. AL: loss 0.2774, stays / mixes 0.3994.
  CX: loss 0.1923, stays / mixes 0.3625. Union of CX and MB: loss 0.2278.
- Reach ratio MB to CX 8.11 at 1 step, AL to CX 0.98, contrast 8.26. CX to MB 1.58 at 1 step.
- Had the criteria been scored at this threshold, C1, C4 and C5 would hold and C2 and C3 would not.

## Plain reading (not a new hypothesis)

- Both registered hypotheses hold in the primary analysis. C2 holds by a small margin (1.40 against
  1.3); the other criteria hold with room.
- The result does not survive the 5-synapse threshold in this data set: there EB no longer mixes before
  it leaves (0.84) and the conditioned walk has not forgotten its start after 8 steps (0.90). In BANC
  an edge carries 3.2 synapses on average, the threshold keeps 79 percent of the giant component and
  44 percent of the mushroom-body set, so the threshold is harsher here than in the male data. The
  registration made all edges the primary analysis for this reason, but the dependence on weak
  connections is a limit of the result and has to be stated with it.
- CX as a whole is not a compartment by the stays / mixes measure in BANC (0.79), and adding MB to CX
  again makes the set leakier, not tighter. The closed part is the ellipsoid-body set.
- Some pairs of starting neurons in EB are fully distinguishable from outside (largest exit distance
  1.00), unlike in the exploratory data sets (0.59 to 0.75 over sampled pairs).
- By the null models run before the registration, these statistics are a property of the wiring
  between cell types. The run therefore shows that this type-level architecture recurs in a third
  animal. It says nothing about neural activity or experience.
