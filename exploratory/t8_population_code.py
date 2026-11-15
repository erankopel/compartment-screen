"""T8. Tier 1 for the candidate itself: a population-code row for Table 5 of the criteria article
(Kopel, BioSystems 270, 106004), evaluated for the heading bump of the fly ellipsoid body.

Written on 15 Nov 2026, AFTER the confirmatory run; revised the same day after the papers were read in full text
(see T8_population_code.md, section Revision). Exploratory arithmetic from published numbers; it reads no
connectome edges and scores no registered criterion. ASCII, no em dashes.

The row. A bump of activity on a ring of neurons labels an angle. Without input its position diffuses:
    var(theta) = 2 t / tau_b            (tau_b = 1 / D, the bump's coherence time; D in rad^2 / s)
which is the phase-diffusion law of the phase-of-firing row with tau_b in place of the oscillator's tau_coh.
Bins of width twice the accumulated error stay distinguishable, so
    alphabet, in place of C1:   K <= pi / sigma_theta = pi * sqrt(tau_b / (2 t))
    one bit (K = 2) survives    t = (pi^2 / 8) tau_b = 1.23 tau_b    (the same 1.23 as for a phase label)
There is no carrier frequency and therefore no Q: the alphabet is set by the hold time alone.
    readout floor, in place of C2:   sigma_0 = w / sqrt(N_s)
for N_s spikes from the bump in the read window, each spike a draw from a profile of standard deviation w
(w = FWHM / 2.355). Both errors add in quadrature:  K <= pi / sqrt(w^2 / N_s + 2 t / tau_b).
C5 splits in two: the hold is set by tau_b, the rewrite by a separate drive (the turn signal), so a long tau_b does
not make the label hard to reassign, unlike a passive oscillator.

Two estimates of D are used:
  - the pooled bound: D <= SD^2 / (2 mean t) from the statistics printed in the figure legends (an upper bound,
    because the spread also holds measurement error and bouts after which the position is unrelated);
  - the fit to the single bouts of Fig. 4i (t8b_standing_bouts.py; read from its result file if present).

Published numbers used (sources in the output file), all checked in the publishers' PDFs unless noted:
  Seelig and Jayaraman 2015, Nature 521, 186, doi:10.1038/nature14446
     16 wedges of 22.5 deg; bump FWHM 82.3 +/- 11.5 deg (one stripe), 90.9 +/- 11.2 deg (darkness);
     standing bouts in darkness, Fig. 4i: duration 6.7 +/- 5.1 s, change of the population-vector angle
     0.017 +/- 0.76 rad; Extended Data Fig. 10b, d, f, h: three more panels in darkness and one with a visible stripe;
     persistence sometimes more than 30 s; standing bouts of at least 20 frames (2.35 s); indicator decay about 0.38 s
  Kim et al. 2017, Science 356, 849, doi:10.1126/science.aal4835
     bump width fixed at 90 deg in the models to match the measurements; stimulated regions 22.5 deg wide
  Turner-Evans et al. 2017, eLife 6, e23496, doi:10.7554/eLife.23496 (PubMed Central text and one page reading)
     both bumps about 100 deg wide; model with 54 E-PG neurons; P-EN neurons fire 3.9 +/- 2.6 Hz while the fly
     stands (loose patch); circuit built to follow at least 200 deg / s. The model's diffusion coefficient
     (1.82e-3 rad^2 / s, fit of sigma0^2 + 2 D t) is the growth of the ERROR between bump and integrated turning
     under a simulated turning input in a model without intrinsic noise. It is an accuracy of the rewrite, not a
     drift of a standing bump, and is kept apart from the hold estimates.
  Hulse et al. 2021, eLife 10, e66039, doi:10.7554/eLife.66039: 46 EPG neurons in the hemibrain (their Table 4)
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = {'sd_change_rad': 0.76, 'mean_bout_s': 6.7, 'sd_bout_s': 5.1, 'n_bouts_printed': 499, 'n_markers_in_figure': 445, 'n_flies': 11,
       'wedges': 16, 'min_standing_s': 2.35,
       'fwhm_deg': {'one stripe': 82.3, 'darkness': 90.9, 'Turner-Evans 2017': 100.0},
       'other_panels': {'ED Fig. 10b, darkness, 6-mm ball, before restart': [6.6, 1.0, 449],
                        'ED Fig. 10d, darkness, 10-mm ball, before restart': [6.2, 0.8, 374],
                        'ED Fig. 10f, darkness, 10-mm ball, at restart': [6.1, 0.9, 374],
                        'ED Fig. 10h, visible stripe, closed loop, before restart': [4.85, 0.74, 388]},
       'model_integration_error': {'D_rad2_per_s': 1.82e-3, 'sigma0_sq_rad2': 1.27e-3, 'input': 'Ornstein-Uhlenbeck turning, tau 120 ms, sigma 50 deg/s'},
       'turn_deg_per_s': [150, 180, 200], 'rig_max_mean_turn_deg_per_s': 35, 'persistence_s': 30, 'indicator_decay_s': 0.38,
       'epg_counts': {'BANC': 45, 'male': 46, 'FlyWire': 47, 'hemibrain (Hulse et al. 2021, Table 4)': 46, 'model of Turner-Evans 2017': 54},
       'pen_rate_standing_hz': [3.9, 2.6]}
res = {'note': 'EXPLORATORY, written after the confirmatory run; published numbers only; revised 15 Nov 2026', 'sources': SRC}

# 1. bump diffusion. For a diffusion with bouts of random length the pooled variance of the change is 2 D E[t]
#    exactly, so the mean bout length is the right divisor. Measurement error of the two population-vector estimates
#    and bouts after which the position is unrelated are inside the 0.76 rad, so D is an upper bound, tau_b a lower.
D_meas = SRC['sd_change_rad'] ** 2 / (2 * SRC['mean_bout_s'])
est = {'pooled bound, Fig. 4i (upper bound on D)': D_meas}
fit_file = os.path.join(HERE, 'results', 't8b_standing_bouts.json')
if os.path.exists(fit_file):
    FIT = json.load(open(fit_file))
    est['single bouts of Fig. 4i (t8b fit)'] = FIT['mixture']['D_rad2_per_s']
    res['t8b'] = {k: FIT['mixture'][k] for k in ('p_unrelated', 's0_rad', 'D_rad2_per_s', 'ci95')}
res['diffusion'] = {}
print('bump diffusion and coherence time (var = 2 D t = 2 t / tau_b)')
for nm, D in est.items():
    tau = 1 / D
    one_bit = math.pi ** 2 / 8 * tau
    res['diffusion'][nm] = {'D_rad2_per_s': D, 'tau_b_s': tau, 'one_bit_s': one_bit}
    print(f'  {nm:<42} D = {D:.4g} rad^2/s   tau_b = {tau:6.1f} s   one bit kept {one_bit:6.1f} s')
print('pooled bounds from the other panels (mean t, SD, n as printed in the legend)')
res['pooled_bounds_other_panels'] = {}
for nm, (mt, sd, n) in SRC['other_panels'].items():
    D = sd ** 2 / (2 * mt)
    res['pooled_bounds_other_panels'][nm] = {'D_rad2_per_s': D, 'tau_b_s': 1 / D}
    print(f'  {nm:<58} D <= {D:.3f} rad^2/s   tau_b >= {1 / D:5.1f} s')

# 2. alphabet against hold time, and hold time against alphabet
print('\nalphabet K <= pi sqrt(tau_b / 2t) after a hold of t (diffusion only)')
res['alphabet_vs_hold'] = {}
for nm, D in est.items():
    row = {t: math.pi / math.sqrt(2 * D * t) for t in (0.1, 0.3, 1.0, 3.0, 6.7, 10.0, 30.0)}
    res['alphabet_vs_hold'][nm] = row
    print(f'  {nm:<42} ' + '  '.join(f't={t:g}s: {k:5.1f}' for t, k in row.items()))
print('hold time for which an alphabet K survives, t = pi^2 tau_b / (2 K^2)')
res['hold_vs_alphabet'] = {}
for nm, D in est.items():
    row = {K: math.pi ** 2 / (2 * D * K * K) for K in (2, 4, 8, 16)}
    res['hold_vs_alphabet'][nm] = row
    print(f'  {nm:<42} ' + '  '.join(f'K={K}: {t:7.2f} s' for K, t in row.items()))

# 3. readout floor: spikes needed in the read window to resolve K positions, N_s = (w K / pi)^2, and the firing rate
#    this asks of each neuron inside the half-maximum width when the window is 0.1 s
READ = 0.1
print('\nreadout floor: spikes from the bump needed in one read window, N_s = (w K / pi)^2, w = FWHM / 2.355')
res['spikes_needed'] = {}
for nm, f in SRC['fwhm_deg'].items():
    w = math.radians(f) / 2.355
    row = {K: (w * K / math.pi) ** 2 for K in (4, 8, 16, 32)}
    inside = f / 360 * SRC['epg_counts']['male']
    rate = {K: v / (inside * READ) for K, v in row.items()}
    res['spikes_needed'][nm] = {'w_rad': w, 'w_deg': f / 2.355, 'N_s': row, 'epg_inside_fwhm': inside, 'rate_hz_per_neuron_at_0p1_s': rate}
    print(f'  FWHM {f:5.1f} deg ({nm:<17}) w = {f / 2.355:4.1f} deg   ' + '  '.join(f'K={K}: {n:6.1f}' for K, n in row.items()) +
          f'   | {inside:4.1f} E-PG inside the FWHM; rate per neuron for a 0.1 s window: ' + ', '.join(f'K={K}: {v:.1f} Hz' for K, v in rate.items()))

# 4. both errors together, for the pooled bound
print('\nalphabet with both errors, K <= pi / sqrt(w^2 / N_s + 2 D t), pooled bound, FWHM 90.9 deg (darkness)')
w = math.radians(90.9) / 2.355
res['combined'] = {}
for Ns in (10, 30, 100, 300):
    row = {t: math.pi / math.sqrt(w * w / Ns + 2 * D_meas * t) for t in (0.0, 0.1, 1.0, 6.7, 30.0)}
    res['combined'][Ns] = row
    print(f'  N_s = {Ns:>3}: ' + '  '.join(f't={t:g}s: {k:5.1f}' for t, k in row.items()))

# 5. the rewrite, the accuracy of the rewrite in the model, and the comparison with a phase label on gamma (T1)
wedge = 360 / SRC['wedges']
res['rewrite'] = {f'{v} deg/s': wedge / v for v in SRC['turn_deg_per_s'] + [SRC['rig_max_mean_turn_deg_per_s']]}
print('\nrewrite: time for the bump to cross one wedge of 22.5 deg at a turn speed: ' +
      ', '.join(f'{k}: {v:.2f} s' for k, v in res['rewrite'].items()))
M = SRC['model_integration_error']
res['integration_error_model'] = {f'{t} s': math.sqrt(M['sigma0_sq_rad2'] + 2 * M['D_rad2_per_s'] * t) for t in (1, 10, 60, 300)}
print('accuracy of the rewrite in the model of Turner-Evans 2017 (error between bump and integrated turning, not a hold): SD after ' +
      ', '.join(f'{k}: {v:.2f} rad' for k, v in res['integration_error_model'].items()))
gamma_one_bit = 0.185          # s, T1: a gamma phase label keeps one bit for 7.4 cycles
res['comparison_gamma'] = {'gamma_one_bit_s': gamma_one_bit, **{nm: v['one_bit_s'] / gamma_one_bit for nm, v in res['diffusion'].items()}}
res['rounds_one_bit'] = {nm: v['one_bit_s'] / READ for nm, v in res['diffusion'].items()}
print('one bit is kept longer than by a gamma phase label (0.185 s) by a factor: ' +
      ', '.join(f'{nm}: {res["comparison_gamma"][nm]:.0f}' for nm in res['diffusion']) + '; in read windows of 0.1 s: ' +
      ', '.join(f'{nm}: {v:.0f}' for nm, v in res['rounds_one_bit'].items()))
res['checks'] = {'observed_persistence_s': SRC['persistence_s'],
                 'alphabet_after_30_s_pooled_bound': math.pi / math.sqrt(2 * D_meas * 30),
                 'alphabet_after_mean_bout_pooled_bound': math.pi / math.sqrt(2 * D_meas * SRC['mean_bout_s']),
                 'hold_over_indicator_decay': SRC['mean_bout_s'] / SRC['indicator_decay_s']}
print(f'check: after 30 s the pooled bound gives K = {res["checks"]["alphabet_after_30_s_pooled_bound"]:.2f} (one bit), the duration the authors '
      f'report as sometimes exceeded; after the mean bout of 6.7 s, K = {res["checks"]["alphabet_after_mean_bout_pooled_bound"]:.2f}')
json.dump(res, open(os.path.join(HERE, 'results', 't8_population_code.json'), 'w'), indent=1)
print('saved')
