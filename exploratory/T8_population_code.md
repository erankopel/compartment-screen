# T8. A population-code row for Table 5, evaluated for the heading bump

ASCII, no em dashes. Written on 15 Nov 2026, after the confirmatory run; revised the same day after the papers
were read in full text (section Revision at the end). Exploratory arithmetic from published numbers. It reads no
connectome edges and scores no registered criterion.

Scripts and numbers:
- `t8_population_code.py` -> `results/t8_population_code.json`, `.log` (the row and its arithmetic)
- `t8b_standing_bouts.py` -> `results/t8b_standing_bouts.json`, `.log`, `results/t8b_fig4_points.csv`
  (the single standing bouts of Seelig and Jayaraman 2015, Fig. 4h and 4i, read from the vector graphics of the
  published PDF; rerun with the path of the PDF as argument, needs pdftocairo)

## Why

Tier 1 of the screen (Kopel, BioSystems 270, 106004) tests whether a physical carrier can hold a label. The
candidate compartment found on the connectome is the heading system of the ellipsoid body, whose label is not a
frequency or a phase but the position of a bump of activity on a ring of neurons. Table 5 of the article has no
row for that scheme. This note adds one and evaluates it.

## The row

| Scheme | Alphabet bound, in place of C1 | Readout floor, in place of C2 | C3 to C7 |
|---|---|---|---|
| Population code on a ring (activity bump) | K <= pi sqrt(tau_b / 2t) after a hold of t | sigma_0 = w / sqrt(N_s): bump width over the root of the spike count in the read window | unchanged, except that C5 splits: the hold is set by tau_b, the rewrite by a separate drive |

- Without input the bump position diffuses, var(theta) = 2 D t = 2 t / tau_b, the same law as the phase of an
  oscillator in Section 4.2 of the article, with the bump's coherence time tau_b = 1 / D in place of tau_coh.
  Bins of width twice the accumulated error stay distinguishable, as in the phase row.
- There is no carrier frequency and so no Q. The alphabet depends on the hold time alone. One bit survives
  (pi^2 / 8) tau_b = 1.23 tau_b, the same factor as for a phase label.
- Readout: each spike of the bump is a draw from a profile of standard deviation w (FWHM / 2.355), so N_s
  spikes locate the bump to w / sqrt(N_s). The two errors add in quadrature.
- C5 is two-sided for a passive carrier, tau_read <= tau_coh <= tau_update, because one time constant sets both
  how long a label lasts and how soon it can be replaced. For a driven bump these are two mechanisms, so a long
  tau_b does not make the label hard to reassign. The window becomes two conditions: the hold time for the
  needed alphabet, and the speed of the drive.

## Published numbers

| quantity | value | source |
|---|---|---|
| wedges of the ellipsoid body | 16, of 22.5 deg | Seelig and Jayaraman 2015, Methods |
| bump width (FWHM) | 82.3 +/- 11.5 deg with one stripe; 90.9 +/- 11.2 deg walking in darkness | same, figure legends |
| standing bouts in darkness, at restart | 6.7 +/- 5.1 s; change of the population-vector angle 0.017 +/- 0.76 rad; r = 0.7, slope 0.96 +/- 0.17; n = 499 printed, 11 flies | same, Fig. 4h, 4i |
| the same, three more panels in darkness | 6.6 +/- 5.1 s, 0.09 +/- 1 rad (n = 449); 6.2 +/- 4.5 s, 0.03 +/- 0.8 (n = 374); 6.1 +/- 4.47 s, 0.04 +/- 0.9 (n = 374) | same, Extended Data Fig. 10b, d, f |
| the same with a visible stripe (closed loop) | 4.85 +/- 3.0 s, 0.04 +/- 0.74 rad (n = 388) | same, Extended Data Fig. 10h |
| standing bout | at least 20 frames (2.35 s); estimates averaged over 5 frames (0.59 s) before the stop and after the restart | same, Methods |
| persistence | sometimes more than 30 s | same, main text |
| indicator | GCaMP6f, decay about 0.38 s, imaged at 8.5 Hz | same, Methods |
| bump width | about 100 deg; fixed at 90 deg in the models of Kim et al. | Turner-Evans et al. 2017; Kim et al. 2017 |
| turns | bumps analysed up to the bin 150 to 180 deg / s; model set to follow at least 200 deg / s; P-EN bandwidth 145 +/- 82 deg / s | Turner-Evans et al. 2017 |
| firing rate | P-EN neurons 3.9 +/- 2.6 Hz while the fly stands (loose patch). No E-PG firing rate was found | Turner-Evans et al. 2017 |
| E-PG neurons | 45 (BANC), 46 (male), 47 (FlyWire), own counts; 46 in the hemibrain; 54 in the model of Turner-Evans et al. | Hulse et al. 2021, Table 4 |
| artificial bump | written by light on a region 22.5 deg wide; the previous bump declines and disappears; afterwards the new bump "either stayed in the induced location for several seconds (fig. S2F) or slowly drifted away" | Kim et al. 2017 |

Seelig JD, Jayaraman V. Nature 521, 186-191 (2015), doi:10.1038/nature14446.
Kim SS, Rouault H, Druckmann S, Jayaraman V. Science 356, 849-853 (2017), doi:10.1126/science.aal4835.
Turner-Evans D et al. eLife 6, e23496 (2017), doi:10.7554/eLife.23496.
Hulse BK et al. eLife 10, e66039 (2021), doi:10.7554/eLife.66039.
Fisher YE et al. Nature 576, 121-125 (2019), doi:10.1038/s41586-019-1772-4 (searched for E-PG firing rates; it
reports E-PG membrane potentials and spike rates of ring neurons).

## Evaluation

### Hold time: the pooled bound

For a diffusion observed over bouts of random length the pooled variance of the change is 2 D E[t], so
D <= 0.76^2 / (2 x 6.7) = 0.043 rad^2 / s and tau_b >= 23 s. This is an upper bound on D: the 0.76 rad also holds
the measurement error of two population-vector estimates and the bouts after which the position is unrelated.
The other panels in darkness give 0.052 to 0.076 rad^2 / s (tau_b at least 13 to 19 s).

The panel with a visible stripe gives the same spread (0.74 rad over 4.85 s) as the panels in darkness measured
the same way (0.8 and 1.0 rad). With a stripe in view and the fly standing the heading has a visual anchor, so
most of the pooled spread is not loss of a memory in darkness.

### Hold time: the single bouts (T8b)

Fig. 4 is stored as vector graphics in the published PDF, so the bouts can be read one by one.

- Reading checked twice. The 445 markers of panel i give 6.71 +/- 5.14 s and 0.020 +/- 0.762 rad, and panel h
  gives r = 0.703: the printed statistics. The markers of h and i come in the same order, and the difference
  computed from h equals the value plotted in i for every bout (largest disagreement 0.004 rad).
- The legend prints n = 499; each panel holds 445 markers. Extended Data Fig. 10a prints n = 449 for the
  6-mm ball. The reason for the difference is not known.
- The spread grows with the duration of the bout, but slowly: rank correlation of |d| with t 0.15 (p = 0.002).

| bout duration (s) | n | SD (rad) | robust SD, 1.4826 MAD | within one wedge | beyond pi / 2 |
|---|---|---|---|---|---|
| up to 3.2 | 102 | 0.71 | 0.48 | 0.64 | 0.05 |
| 3.2 to 4.2 | 90 | 0.62 | 0.39 | 0.69 | 0.02 |
| 4.2 to 5.5 | 62 | 0.80 | 0.47 | 0.61 | 0.05 |
| 5.5 to 7.5 | 63 | 0.75 | 0.44 | 0.62 | 0.06 |
| 7.5 to 11 | 59 | 0.83 | 0.63 | 0.46 | 0.05 |
| 11 to 36 | 69 | 0.91 | 0.71 | 0.42 | 0.09 |

- Least squares d^2 = a + 2 D t: a = 0.43 rad^2, D = 0.011 rad^2 / s (bootstrap 95 percent -0.0001 to 0.025).
  Most of the variance is there at the shortest bouts.
- Maximum likelihood for d ~ (1 - p) Normal(0, s0^2 + 2 D t) + p Uniform(-pi, pi):

| | estimate | bootstrap 95 percent |
|---|---|---|
| p, share of bouts after which the position is unrelated | 0.105 | 0.060 to 0.154 |
| s0, error that does not depend on duration | 0.32 rad (0.82 wedge) | 0.10 to 0.44 |
| D | 0.012 rad^2 / s | 0.0049 to 0.020 |
| tau_b = 1 / D | 83 s | 49 to 204 |
| one bit kept, 1.23 tau_b | 103 s | 61 to 251 |

- Without diffusion the fit is worse (likelihood ratio p = 6e-5); without the unrelated share far worse
  (p = 6e-30); letting the unrelated share grow with duration does not help (p = 0.78; time constant 534 s).
- So three things are in the 0.76 rad: an error of about 0.3 rad that is there at any duration (it includes
  the measurement error and the 0.59 s after the restart), about one bout in ten after which the position is
  unrelated to the one before, at any duration, and a slow diffusion, 3.6 times slower than the pooled bound.
- The longest bout is 36 s. One bit for 103 s is an extrapolation.

### Alphabet and hold

| hold time | 0.1 s | 0.3 s | 1 s | 3 s | 6.7 s | 10 s | 30 s |
|---|---|---|---|---|---|---|---|
| alphabet, pooled bound (tau_b = 23 s) | 34 | 19 | 11 | 6.2 | 4.1 | 3.4 | 2.0 |
| alphabet, single bouts, diffusion only (tau_b = 83 s) | 64 | 37 | 20 | 12 | 7.8 | 6.4 | 3.7 |
| alphabet, single bouts, with s0 | 9.6 | 9.4 | 8.8 | 7.5 | 6.1 | 5.4 | 3.5 |

| alphabet kept | 2 (one bit) | 4 | 8 | 16 (one wedge) |
|---|---|---|---|---|
| hold time, pooled bound | 29 s | 7.2 s | 1.8 s | 0.45 s |
| hold time, single bouts, diffusion only | 103 s | 26 s | 6.4 s | 1.6 s |

- After the mean standing bout of 6.7 s the pooled bound leaves about four distinguishable headings (two
  bits), the single bouts six to eight; after 30 s, the duration the authors report as sometimes exceeded,
  two (one bit) and three to four.
- The single-bout rows hold for the nine bouts in ten that keep the position.
- The anatomy caps the alphabet near 16 (wedges) for a reader that takes one column at a time; 45 to 47 E-PG
  neurons carry it.

### Readout floor

- Resolving 8 positions takes 2.4 to 3.6 spikes from the bump in the read window (widths 82 to 100 deg), 16
  positions 10 to 14, 32 positions 39 to 57.
- With about a quarter of the E-PG neurons inside the half-maximum width (10 to 13 neurons) and a read window
  of 0.1 s, that asks 2.3 to 2.8 Hz of each of them for 8 positions, 9 to 11 Hz for 16, 37 to 45 Hz for 32.
- No measured E-PG firing rate was found (Fisher et al. 2019 report membrane potentials). The only patched
  columnar type, P-EN, fires 3.9 +/- 2.6 Hz while the fly stands. If E-PG neurons fire at rates of that order,
  one window of 0.1 s resolves about 8 positions and 16 need a longer window or higher rates. This stays a
  requirement, not a verdict.

### Rewrite and writing

- Rewrite: at 150 to 200 deg / s the bump crosses one wedge in 0.11 to 0.15 s.
- Writing (C6): Kim et al. 2017 wrote a bump at a chosen place with light on a region 22.5 deg wide. The
  previous bump declined and disappeared, and after the light ended the new bump "either stayed in the induced
  location for several seconds (fig. S2F) or slowly drifted away". In darkness "the bump drifted gradually
  around the EB", not to random distant places. One label at a time, writable from outside, kept by the circuit.
- Accuracy of the rewrite (model): in the model of Turner-Evans et al. 2017, which has no intrinsic noise, the
  error between the bump and the integrated turning grows as sigma0^2 + 2 D t with D = 1.82e-3 rad^2 / s under
  a simulated turning input: 0.19 rad after 10 s, 0.47 rad after 60 s. This is not a hold time.

### Against the oscillatory carriers of the article

A phase label on cortical gamma has five bins within one cycle and keeps one bit for 0.185 s (T1). The bump
keeps one bit at least 155 times longer by the pooled bound (29 s, 290 read windows of 0.1 s), and about 560
times longer by the single bouts (103 s, extrapolated).

## Reading

- The candidate's own code passes Tier 1 in its population-code form: a small alphabet (two to three bits)
  held for seconds, one bit for half a minute or more, rewritable in about a tenth of a second per wedge, and
  writable from outside.
- This is the scheme that meets the hold requirement where oscillatory phase cannot: among the neural rhythms
  of Table 3 of the article the longest-lived phase bit is that of theta, 1.2 s (T1).
- Passing is necessary, not sufficient, as for every row of the screen.

## Limits

- Diffusion plus an unrelated share is a model chosen after looking at the points. The unrelated share may be
  jumps, fading and return of the bump, or turning during the 0.59 s after the restart; the data do not say.
- The points were read from a figure. The reading reproduces the printed statistics and agrees between the two
  panels, but the figure holds 445 markers where the legend prints 499.
- All hold estimates are for flies standing on a ball in darkness between walking bouts, bouts of 2.6 to 36 s.
  In tethered flight in darkness the bump drifts within seconds and seldom follows the fly's turns (Kim et al.
  2017, Fig. 3; figure reading).
- Calcium imaging at 8.5 Hz with an indicator decay of 0.38 s; spikes were not measured. The readout floor is
  stated as a requirement.
- The Extended Data panels are images in the PDF and were used through their printed statistics only.

## Revision, 15 Nov 2026 (after reading the papers)

- Seelig and Jayaraman 2015 and Kim et al. 2017 were read in the publishers' PDFs. The numbers of the first
  version stand. Before, the statistics of Fig. 4 rested on two readings through an extraction tool, and Kim et
  al. on the abstract.
- New: the single bouts (T8b), the other panels of Extended Data Fig. 10, the writing experiment of Kim et al.,
  the firing rates asked by the readout floor, the hemibrain count of E-PG neurons.
- Corrected: the first version used the diffusion coefficient of the model of Turner-Evans et al. 2017 as a
  second estimate of the hold (tau_b = 550 s; alphabets of 165 to 9.5). That coefficient is the growth of the
  integration error under a turning input in a model without intrinsic noise. It says nothing about a standing
  bump. It is now reported as accuracy of the rewrite, and the second hold estimate comes from the single bouts.
  Its value and convention (1.82e-3 rad^2 / s, 2 D t) were confirmed in the text.
- Hulse et al. 2021 is the hemibrain paper. It was searched for cell counts and for statements on physiology.
  The hemibrain connectivity files stay unopened.
