"""T8b. The standing bouts of Seelig and Jayaraman 2015 (Nature 521, 186, doi:10.1038/nature14446), Fig. 4h and 4i,
read point by point, and what they say about the hold time of the heading bump.

Written on 15 Nov 2026, AFTER the confirmatory run. Exploratory: it reads no connectome edges and scores no
registered criterion. The models were chosen after looking at the points. ASCII, no em dashes.

Page 5 of the published PDF stores Fig. 4 as vector graphics, so every marker of panels h and i has its plotted
coordinates. The script converts the page with pdftocairo, reads the markers, calibrates the axes from the tick
marks, and checks the reading in two ways: against the statistics printed in the legend, and panel against
panel (the markers of h and i come in the same order, so the difference computed from h must equal the value
plotted in i, bout by bout).

Then, for the difference d of the population-vector angle over a standing bout of duration t:
  1. the spread of d in bins of t;
  2. least squares  d^2 = a + 2 D t;
  3. maximum likelihood for   d ~ (1 - p) Normal(0, s0^2 + 2 D t) + p Uniform(-pi, pi):
     a share p of bouts after which the position is unrelated to the one before, and for the others a
     duration-independent error s0 (it contains the measurement error of the two estimates and whatever happens in
     the 0.59 s after the restart over which the second estimate is averaged) plus diffusion;
     compared with the same model without diffusion, without the unrelated share, and with a share that grows
     with t.
Usage: python3 t8b_standing_bouts.py <nature14446.pdf>      (needs pdftocairo, poppler)
       python3 t8b_standing_bouts.py                         (analysis only, from results/t8b_fig4_points.csv)
"""
import csv, json, math, os, re, subprocess, sys, tempfile
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')
CSV = os.path.join(R, 't8b_fig4_points.csv')
LEGEND = {'t_mean_s': 6.7, 't_sd_s': 5.1, 'd_mean_rad': 0.017, 'd_sd_rad': 0.76, 'r': 0.7, 'n_printed': 499}
SEED = 20261115
LOG = []


def log(s=''):
    print(s, flush=True)
    LOG.append(s)


def wrap(a):
    return np.angle(np.exp(1j * a))


def extract(pdf):
    """markers of Fig. 4h and 4i from page 5 of the PDF, in data units"""
    with tempfile.TemporaryDirectory() as tmp:
        svg = os.path.join(tmp, 'p5.svg')
        subprocess.run(['pdftocairo', '-svg', '-f', '5', '-l', '5', pdf, svg], check=True)
        txt = open(svg).read()
    num = r'-?\d+\.?\d*(?:e-?\d+)?'
    marks, lines = [], []
    for e in re.findall(r'<path[^>]*?/>', txt):
        m = re.search(r'transform="matrix\(([^)]*)\)"', e)
        if not m:
            continue
        a, b, c, dd, tx, ty = [float(v) for v in m.group(1).split(',')]
        d = re.search(r' d="([^"]*)"', e).group(1)
        shape = re.sub(num, 'N', d).strip()
        xy = np.array([float(v) for v in re.findall(num, d)]).reshape(-1, 2)
        X, Y = a * xy[:, 0] + c * xy[:, 1] + tx, b * xy[:, 0] + dd * xy[:, 1] + ty
        if shape.startswith('M N N C N N N N N N C N N N N N N C N N N N N N C') and 1.4 < X.max() - X.min() < 1.7:
            marks.append(((X.min() + X.max()) / 2, (Y.min() + Y.max()) / 2))
        elif shape == 'M N N L N N':
            lines.append((X[0], Y[0], X[1], Y[1]))
    M, L = np.array(marks), np.array(lines)
    ys = np.sort(M[:, 1])
    cut = (ys[np.argmax(np.diff(ys))] + ys[np.argmax(np.diff(ys)) + 1]) / 2     # largest gap separates h (above) from i
    out = {}
    for nm, sel in (('h', M[:, 1] < cut), ('i', M[:, 1] > cut)):
        P = M[sel]
        x0, x1, y0, y1 = P[:, 0].min() - 3, P[:, 0].max() + 3, P[:, 1].min() - 3, P[:, 1].max() + 3
        near = L[(L[:, 0] > x0) & (L[:, 0] < x1) & (L[:, 1] > y0) & (L[:, 1] < y1)]
        ln = np.hypot(near[:, 2] - near[:, 0], near[:, 3] - near[:, 1])
        vert = near[(np.abs(near[:, 2] - near[:, 0]) < 1e-6) & (ln < 0.8)]          # ticks of the x axis
        hori = near[(np.abs(near[:, 3] - near[:, 1]) < 1e-6) & (ln < 0.8)]          # ticks of the y axis
        out[nm] = (P, np.sort(vert[:, 0]), np.sort(hori[:, 1]))
    Ph, xt, yt = out['h']
    assert len(xt) == 3 and len(yt) == 3, (xt, yt)
    before = np.polyval(np.polyfit(xt, [-math.pi, 0, math.pi], 1), Ph[:, 0])
    restart = np.polyval(np.polyfit(yt, [math.pi, 0, -math.pi], 1), Ph[:, 1])       # page y runs downward
    Pi, xt, yt = out['i']
    assert len(xt) == 3 and len(yt) == 5, (xt, yt)
    t = np.polyval(np.polyfit(xt, [10, 20, 30], 1), Pi[:, 0])
    d = np.polyval(np.polyfit(yt, [math.pi, math.pi / 2, 0, -math.pi / 2, -math.pi], 1), Pi[:, 1])
    assert len(t) == len(before)
    return t, d, before, restart


if len(sys.argv) > 1:
    t, d, before, restart = extract(sys.argv[1])
    os.makedirs(R, exist_ok=True)
    with open(CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['bout', 'duration_s', 'difference_rad_panel_i', 'before_stop_rad_panel_h', 'at_restart_rad_panel_h'])
        for k in range(len(t)):
            w.writerow([k + 1, f'{t[k]:.4f}', f'{d[k]:.4f}', f'{before[k]:.4f}', f'{restart[k]:.4f}'])
rows = list(csv.DictReader(open(CSV)))
t = np.array([float(r['duration_s']) for r in rows]); d = np.array([float(r['difference_rad_panel_i']) for r in rows])
before = np.array([float(r['before_stop_rad_panel_h']) for r in rows]); restart = np.array([float(r['at_restart_rad_panel_h']) for r in rows])
res = {'note': 'EXPLORATORY, written after the confirmatory run; points read from the vector graphics of Fig. 4h, 4i',
       'source': 'Seelig and Jayaraman 2015, Nature 521, 186, doi:10.1038/nature14446, Fig. 4h, 4i', 'legend': LEGEND}

# checks of the reading
n = len(t)
dh = wrap(restart - before)
res['check'] = {'markers_per_panel': n, 't_mean_s': float(t.mean()), 't_sd_s': float(t.std(ddof=1)), 't_min_s': float(t.min()), 't_max_s': float(t.max()),
                'd_mean_rad': float(d.mean()), 'd_sd_rad': float(d.std(ddof=1)), 'r_panel_h': float(np.corrcoef(before, restart)[0, 1]),
                'panel_h_vs_i_max_abs_difference_rad': float(np.abs(wrap(dh - d)).max()),
                'sign': 'the value plotted in panel i equals (at restart) minus (before stop)'}
c = res['check']
log(f'markers per panel {n} (the legend prints n = {LEGEND["n_printed"]})')
log(f'duration {c["t_mean_s"]:.2f} +/- {c["t_sd_s"]:.2f} s (legend 6.7 +/- 5.1), range {c["t_min_s"]:.2f} to {c["t_max_s"]:.2f} s')
log(f'difference {c["d_mean_rad"]:.3f} +/- {c["d_sd_rad"]:.3f} rad (legend 0.017 +/- 0.76); r of panel h {c["r_panel_h"]:.3f} (legend 0.7)')
log(f'panel h against panel i, bout by bout: largest disagreement {c["panel_h_vs_i_max_abs_difference_rad"]:.4f} rad')

# 1. spread in bins of duration
WEDGE = math.pi / 8
log('\n1. spread of the difference in bins of bout duration')
log('   bin (s)     n  mean t    SD  1.4826 MAD  median |d|  within one wedge  beyond pi/2')
res['bins'] = []
for a, b in zip([0, 3.2, 4.2, 5.5, 7.5, 11], [3.2, 4.2, 5.5, 7.5, 11, 40]):
    m = (t >= a) & (t < b); x = d[m]
    row = {'from_s': a, 'to_s': b, 'n': int(m.sum()), 'mean_t_s': float(t[m].mean()), 'sd_rad': float(x.std(ddof=1)),
           'mad_sd_rad': float(1.4826 * np.median(np.abs(x - np.median(x)))), 'median_abs_rad': float(np.median(np.abs(x))),
           'share_within_one_wedge': float((np.abs(x) <= WEDGE).mean()), 'share_beyond_half_pi': float((np.abs(x) > math.pi / 2).mean())}
    res['bins'].append(row)
    log(f'   {a:>4} to {b:<4} {row["n"]:>4} {row["mean_t_s"]:6.2f} {row["sd_rad"]:6.3f} {row["mad_sd_rad"]:9.3f} {row["median_abs_rad"]:10.3f} '
        f'{row["share_within_one_wedge"]:14.2f} {row["share_beyond_half_pi"]:12.2f}')
sp = spearmanr(np.abs(d), t)
res['all'] = {'sd_rad': float(d.std(ddof=1)), 'mad_sd_rad': float(1.4826 * np.median(np.abs(d - np.median(d)))), 'median_abs_rad': float(np.median(np.abs(d))),
              'share_within_one_wedge': float((np.abs(d) <= WEDGE).mean()), 'share_beyond_half_pi': float((np.abs(d) > math.pi / 2).mean()),
              'spearman_abs_d_vs_t': float(sp.statistic), 'spearman_p': float(sp.pvalue)}
log(f'   all: SD {res["all"]["sd_rad"]:.3f}, 1.4826 MAD {res["all"]["mad_sd_rad"]:.3f}, within one wedge {res["all"]["share_within_one_wedge"]:.2f}, '
    f'beyond pi/2 {res["all"]["share_beyond_half_pi"]:.3f}; rank correlation of |d| with t {sp.statistic:.3f} (p = {sp.pvalue:.4f})')

# 2. least squares d^2 = a + 2 D t
rng = np.random.default_rng(SEED)
A = np.c_[np.ones(n), t]
co = np.linalg.lstsq(A, d ** 2, rcond=None)[0]
bs = np.array([np.linalg.lstsq(A[k], d[k] ** 2, rcond=None)[0] for k in (rng.integers(0, n, n) for _ in range(4000))])
lo, hi = np.percentile(bs[:, 1], [2.5, 97.5]) / 2
res['least_squares'] = {'a_rad2': float(co[0]), 'D_rad2_per_s': float(co[1] / 2), 'D_ci95': [float(lo), float(hi)], 'bootstrap': 4000}
log(f'\n2. least squares d^2 = a + 2 D t: a = {co[0]:.3f} rad^2, D = {co[1] / 2:.4f} rad^2/s (bootstrap 95 percent {lo:.4f} to {hi:.4f})')
log(f'   pooled bound of T8, D <= SD^2 / (2 mean t) = {d.var(ddof=1) / (2 * t.mean()):.4f} rad^2/s')


# 3. mixture
def nll(th, tt, dd, kind='const'):
    p0, s0, D = 1 / (1 + math.exp(-th[0])), math.exp(th[1]), math.exp(th[2])
    p = p0 if kind == 'const' else 1 - (1 - p0) * np.exp(-tt / math.exp(th[3]))
    v = s0 ** 2 + 2 * D * tt
    return -np.log((1 - p) * np.exp(-dd ** 2 / (2 * v)) / np.sqrt(2 * math.pi * v) + p / (2 * math.pi)).sum()


def fit(f, starts):
    best = None
    for s in starts:
        r = minimize(f, s, method='Nelder-Mead', options={'xatol': 1e-7, 'fatol': 1e-9, 'maxiter': 40000})
        if best is None or r.fun < best.fun:
            best = r
    return best


full = fit(lambda th: nll(th, t, d), ([-2, -1.5, -5], [-1, -1, -3], [-2, -2, -7], [-3, -1, -4]))
noD = fit(lambda th: nll([th[0], th[1], -40], t, d), ([-2, -1.5], [-1, -0.5]))
noP = fit(lambda th: nll([-40, th[0], th[1]], t, d), ([-1, -3], [-0.5, -5]))
grow = fit(lambda th: nll(th, t, d, 'grow'), ([-2, -1.5, -5, 3], [-3, -1.2, -6, 4], [-2, -1.5, -4, 5]))
p, s0, D = 1 / (1 + math.exp(-full.x[0])), math.exp(full.x[1]), math.exp(full.x[2])
bb = []
for _ in range(1000):
    k = rng.integers(0, n, n)
    r = minimize(lambda th: nll(th, t[k], d[k]), full.x, method='Nelder-Mead', options={'maxiter': 6000})
    bb.append([1 / (1 + math.exp(-r.x[0])), math.exp(r.x[1]), math.exp(r.x[2])])
bb = np.array(bb)
ci = {k: [float(v) for v in np.percentile(bb[:, j], [2.5, 97.5])] for j, k in enumerate(('p', 's0_rad', 'D_rad2_per_s'))}
lr = lambda small, big, df: float(chi2.sf(2 * (small - big), df))
res['mixture'] = {'p_unrelated': p, 's0_rad': s0, 'D_rad2_per_s': D, 'ci95': ci, 'bootstrap': 1000, 'nll': float(full.fun),
                  'without_diffusion': {'p': 1 / (1 + math.exp(-noD.x[0])), 's0_rad': math.exp(noD.x[0 + 1]), 'nll': float(noD.fun), 'p_value_against_full': lr(noD.fun, full.fun, 1)},
                  'without_unrelated_share': {'s0_rad': math.exp(noP.x[0]), 'D_rad2_per_s': math.exp(noP.x[1]), 'nll': float(noP.fun), 'p_value_against_full': lr(noP.fun, full.fun, 1)},
                  'unrelated_share_growing_with_t': {'p_at_0': 1 / (1 + math.exp(-grow.x[0])), 'tau_loss_s': math.exp(grow.x[3]), 'nll': float(grow.fun),
                                                     'p_value_full_against_this': lr(full.fun, grow.fun, 1)}}
m = res['mixture']
log('\n3. d ~ (1 - p) Normal(0, s0^2 + 2 D t) + p Uniform(-pi, pi), maximum likelihood')
log(f'   p = {p:.3f} ({ci["p"][0]:.3f} to {ci["p"][1]:.3f}); s0 = {s0:.3f} rad ({ci["s0_rad"][0]:.3f} to {ci["s0_rad"][1]:.3f}); '
    f'D = {D:.4f} rad^2/s ({ci["D_rad2_per_s"][0]:.4f} to {ci["D_rad2_per_s"][1]:.4f}); -log L = {full.fun:.2f}')
log(f'   without diffusion: -log L = {noD.fun:.2f} (likelihood ratio p = {m["without_diffusion"]["p_value_against_full"]:.1e})')
log(f'   without the unrelated share: -log L = {noP.fun:.2f} (p = {m["without_unrelated_share"]["p_value_against_full"]:.1e})')
log(f'   unrelated share growing with t: -log L = {grow.fun:.2f}, time constant {m["unrelated_share_growing_with_t"]["tau_loss_s"]:.0f} s '
    f'(no better: p = {m["unrelated_share_growing_with_t"]["p_value_full_against_this"]:.2f})')

# 4. consequences for the population-code row
tau = 1 / D
tau_ci = [1 / ci['D_rad2_per_s'][1], 1 / ci['D_rad2_per_s'][0]]
one_bit = math.pi ** 2 / 8 * tau
res['row'] = {'tau_b_s': tau, 'tau_b_ci95_s': tau_ci, 'one_bit_s': one_bit, 'one_bit_ci95_s': [math.pi ** 2 / 8 * v for v in tau_ci],
              'retained_share': 1 - p,
              'alphabet_diffusion_only': {str(x): math.pi / math.sqrt(2 * D * x) for x in (1.0, 3.0, 6.7, 10.0, 30.0)},
              'alphabet_with_s0': {str(x): math.pi / math.sqrt(s0 ** 2 + 2 * D * x) for x in (0.0, 1.0, 3.0, 6.7, 10.0, 30.0)},
              'hold_for_alphabet_diffusion_only_s': {str(K): math.pi ** 2 / (2 * D * K * K) for K in (2, 4, 8, 16)},
              's0_in_wedges': s0 / WEDGE}
r_ = res['row']
log(f'\n4. for the row: tau_b = 1 / D = {tau:.0f} s ({tau_ci[0]:.0f} to {tau_ci[1]:.0f}); one bit kept {one_bit:.0f} s '
    f'({r_["one_bit_ci95_s"][0]:.0f} to {r_["one_bit_ci95_s"][1]:.0f}) in the {1 - p:.0%} of bouts that keep the position')
log('   alphabet after a hold of t, diffusion only: ' + ', '.join(f'{k} s: {v:.1f}' for k, v in r_['alphabet_diffusion_only'].items()))
log('   alphabet after a hold of t, with s0:        ' + ', '.join(f'{k} s: {v:.1f}' for k, v in r_['alphabet_with_s0'].items()))
log('   hold time for an alphabet K, diffusion only: ' + ', '.join(f'K = {k}: {v:.1f} s' for k, v in r_['hold_for_alphabet_diffusion_only_s'].items()))
log(f'   s0 = {s0 / WEDGE:.2f} wedges of 22.5 deg')
json.dump(res, open(os.path.join(R, 't8b_standing_bouts.json'), 'w'), indent=1)
open(os.path.join(R, 't8b_standing_bouts.log'), 'w').write('\n'.join(LOG) + '\n')
log('saved')
