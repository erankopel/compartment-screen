"""T3 (and T4b). EXPLORATORY, not pre-registered. The 727-neuron anterior-visual-pathway and
central-complex set found by E1 of the maleCNS-depth repository (smallest-conductance sweep cut of the
weight >= 1 synapse-flow chain), looked at as a candidate compartment:

 A. does the walk mix inside the set faster than it leaves it (metastability)?
 B. where are its entry and exit interfaces, and do they overlap?
 C. the same quantities for other anatomically defined sets (annotation class), for comparison;
 D. T4b: after r steps, how much of the difference between two starting neurons of the set is still
    visible inside the set and how much outside it.

Uses the helpers of code/confirmatory/common.py of the repository unchanged. Two BLAS threads are
allowed here (the confirmatory runs used one). ASCII, no em dashes.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components

REPO = os.environ['MALECNS_REPO']
sys.path.insert(0, os.path.join(REPO, 'code', 'confirmatory'))
from common import load_edges, giant_scc, row_stochastic, stationary, annotations, log, _json_default  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', 't3_subsystem.json')
t00 = time.time()

pre, post, w, body_ids = load_edges()
N = len(body_ids)
ann = annotations(body_ids, cols=('type', 'class', 'superclass', 'somaSide'))
e1 = json.load(open(os.path.join(REPO, 'exploratory', 'results', 'e1_slow_mode.json')))
cut_ids = np.array(e1['w1']['cut_body_ids'], dtype=np.int64)
in_cut = np.isin(body_ids, cut_ids)
log(f'G = {N:,}; cut ids found in G: {int(in_cut.sum())} of {len(cut_ids)}')


def top_counts(values, weights, top=10):
    u, inv = np.unique(values.astype(str), return_inverse=True)
    s = np.bincount(inv, weights=weights, minlength=len(u))
    o = np.argsort(-s)[:top]
    tot = s.sum()
    return [[str(u[i]), float(s[i] / tot)] for i in o]


def sub_spectrum(PSS):
    n = PSS.shape[0]
    if n <= 1600:
        ev = np.linalg.eigvals(PSS.toarray())
    else:
        ev = sla.eigs(PSS, k=6, which='LM', tol=1e-9, maxiter=20000, ncv=60,
                      return_eigenvectors=False, rng=np.random.default_rng(11))
    ev = ev[np.argsort(-np.abs(ev))]
    return ev


def subsystem(P, pi, S, name):
    """Substochastic kernel of the walk killed on leaving S."""
    idx = np.flatnonzero(S)
    PSS = P[idx][:, idx].tocsr()
    stay = np.asarray(PSS.sum(axis=1)).ravel()
    mass = float(pi[idx].sum())
    Fout = float((pi[idx] * (1 - stay)).sum())
    ev = sub_spectrum(PSS)
    lam1 = float(np.real(ev[0]))
    lam2 = float(np.abs(ev[1]))
    ncomp, lab = connected_components(PSS, directed=True, connection='strong')
    big = int(np.bincount(lab).max())
    r = {'name': name, 'size': int(len(idx)), 'stationary_mass': mass, 'escape_per_step': Fout / mass,
         'kac_sojourn_steps': mass / Fout, 'lambda1': lam1, 'lambda2_abs': lam2,
         'lambda2_is_complex': bool(abs(np.imag(ev[1])) > 1e-9),
         'residence_steps_qsd': 1 / (1 - lam1), 'relaxation_steps': 1 / (1 - lam2 / lam1),
         'ratio_residence_over_relaxation': (1 / (1 - lam1)) / (1 / (1 - lam2 / lam1)),
         'largest_internal_scc': big}
    log(f'{name:<34s} n {len(idx):>6,} mass {mass:.4f} escape {Fout / mass:.4f} lam1 {lam1:.4f} '
        f'|lam2| {lam2:.4f} residence {r["residence_steps_qsd"]:.2f} relaxation {r["relaxation_steps"]:.2f} '
        f'ratio {r["ratio_residence_over_relaxation"]:.2f} internal SCC {big}')
    return r, idx, PSS


def chain(thresh):
    keep, A = giant_scc(pre, post, w, N, thresh=thresh)
    P = row_stochastic(A)
    pi, it, dl = stationary(P, tol=1e-13, maxit=20000)
    log(f'\n=== weight >= {thresh}: SCC n = {len(keep):,}, nnz = {P.nnz:,}, stationary in {it} iterations '
        f'(delta {dl:.1e})')
    return keep, P, pi


res = {'note': 'EXPLORATORY, not pre-registered; male CNS v1.0; cut ids from exploratory/results/e1_slow_mode.json'}

# ----------------------------------------------------------------------------- weight >= 1
keep, P, pi = chain(1)
n = len(keep)
S = in_cut[keep]
typ, cls, sup = ann['type'][keep], ann['class'][keep], ann['superclass'][keep]
r0, idx, PSS = subsystem(P, pi, S, 'AVP/CX cut (727)')
r0['check_against_repo'] = {'cut_mass_repo': e1['w1']['cut_mass'], 'escape_repo': e1['w1']['escape_small']}
res['w1'] = {'n': n, 'cut': r0}

# layers inside the cut: anterior visual pathway (MeTu, TuBu) versus the rest (ring, columnar neurons)
is_avp = np.array([t.startswith('MeTu') or t.startswith('TuBu') for t in typ])
layer = np.where(~S, 0, np.where(is_avp, 1, 2))          # 0 outside, 1 AVP, 2 core
rows = np.repeat(np.arange(n), np.diff(P.indptr)); cols = P.indices
f = pi[rows] * P.data
flow = np.zeros((3, 3))
np.add.at(flow, (layer[rows], layer[cols]), f)
res['w1']['layers'] = {'sizes': {'outside': int((layer == 0).sum()), 'AVP (MeTu, TuBu)': int((layer == 1).sum()),
                                 'core (rest of the cut)': int((layer == 2).sum())},
                       'mass': [float(pi[layer == k].sum()) for k in range(3)],
                       'flow_per_step_rows_from_cols_to': flow.tolist()}
log('layers: sizes ' + str(res['w1']['layers']['sizes']))
log('flow per step (rows = from outside, AVP, core; columns = to outside, AVP, core):')
for k in range(3):
    log('   ' + '  '.join(f'{x:.3e}' for x in flow[k]))
for nm, m in (('AVP layer alone', layer == 1), ('core alone', layer == 2)):
    rr, _, _ = subsystem(P, pi, m, nm)
    res['w1'][nm] = rr

# entry and exit interfaces of the cut
Sr, Sc = S[rows], S[cols]
fin, fout = (~Sr) & Sc, Sr & (~Sc)
inflow = np.bincount(cols[fin], weights=f[fin], minlength=n)      # per neuron of S
outflow = np.bincount(rows[fout], weights=f[fout], minlength=n)


def carriers(x, share):
    o = np.argsort(-x)
    c = np.cumsum(x[o]) / x.sum()
    return o[:int(np.searchsorted(c, share) + 1)]


iface = {}
for nm, x in (('entry', inflow), ('exit', outflow)):
    c50, c90 = carriers(x, 0.5), carriers(x, 0.9)
    iface[nm] = {'total_flow': float(x.sum()), 'neurons_for_50_percent': int(len(c50)),
                 'neurons_for_90_percent': int(len(c90)), 'types_by_flow': top_counts(typ, x, 10)}
e90, x90 = set(carriers(inflow, 0.9).tolist()), set(carriers(outflow, 0.9).tolist())
iface['jaccard_of_90_percent_sets'] = len(e90 & x90) / len(e90 | x90)
iface['neurons_in_neither_90_percent_set'] = int(S.sum() - len(e90 | x90))
iface['external_sources_by_type'] = top_counts(typ[rows[fin]], f[fin], 10)
iface['external_sources_by_superclass'] = top_counts(sup[rows[fin]], f[fin], 6)
iface['external_targets_by_type'] = top_counts(typ[cols[fout]], f[fout], 10)
iface['external_targets_by_superclass'] = top_counts(sup[cols[fout]], f[fout], 6)
res['w1']['interfaces'] = iface
log(f'entry: 50% of inflow on {iface["entry"]["neurons_for_50_percent"]} neurons, 90% on '
    f'{iface["entry"]["neurons_for_90_percent"]}; types {iface["entry"]["types_by_flow"][:6]}')
log(f'exit : 50% of outflow on {iface["exit"]["neurons_for_50_percent"]} neurons, 90% on '
    f'{iface["exit"]["neurons_for_90_percent"]}; types {iface["exit"]["types_by_flow"][:6]}')
log(f'Jaccard of the two 90% sets {iface["jaccard_of_90_percent_sets"]:.3f}; neurons in neither '
    f'{iface["neurons_in_neither_90_percent_set"]}')
log(f'external sources {iface["external_sources_by_type"][:6]}')
log(f'external targets {iface["external_targets_by_type"][:6]}')

# sojourn of a walker that enters with the stationary entry distribution
D = PSS.toarray()
a = inflow[idx] / inflow[idx].sum()
v, surv = a.copy(), []
for t in range(0, 129):
    surv.append(float(v.sum()))
    v = v @ D
res['w1']['sojourn'] = {'survival': {str(t): surv[t] for t in (0, 1, 2, 4, 8, 16, 32, 64, 128)},
                        'mean_steps': float(sum(surv)), 'kac': r0['kac_sojourn_steps']}
log('survival inside after t steps (entry distribution): ' +
    ', '.join(f't={t}: {surv[t]:.3f}' for t in (1, 2, 4, 8, 16, 32, 64)) + f'; mean {sum(surv):.2f} steps')

# the walk conditioned to stay (Doob transform) and its Dobrushin coefficient
evals, V = np.linalg.eig(D)
k1 = int(np.argmax(np.real(evals)))
lam1 = float(np.real(evals[k1]))
hvec = np.abs(np.real(V[:, k1]))
Q = D * hvec[None, :] / (lam1 * hvec[:, None])
log(f'Doob transform: row sums in [{Q.sum(1).min():.6f}, {Q.sum(1).max():.6f}], min h {hvec.min():.2e}')
dob, Qr, done = {}, np.eye(len(idx)), 0
for r in (1, 2, 4, 8, 16, 32, 64):
    for _ in range(r - done):
        Qr = Qr @ Q
    done = r
    m = 0.0
    for i in range(len(idx) - 1):
        m = max(m, 0.5 * float(np.abs(Qr[i][None, :] - Qr[i + 1:]).sum(axis=1).max()))
    dob[str(r)] = m
res['w1']['dobrushin_conditioned_walk'] = dob
log('Dobrushin coefficient of the conditioned walk, all pairs: ' + ', '.join(f'r={k}: {x:.4f}' for k, x in dob.items()))
json.dump(res, open(OUT, 'w'), indent=1, default=_json_default)

# comparison sets
comp = []
sets = {'class CX (all)': cls == 'CX',
        'Kenyon cells': cls == 'Kenyon_Cell',
        'mushroom body (KC, MBON, DAN)': np.isin(cls, ['Kenyon_Cell', 'MBON', 'DAN']),
        'antennal lobe (ALPN, ALLN)': np.isin(cls, ['ALPN', 'ALLN']),
        'antennal lobe with olfactory sensory': np.isin(cls, ['ALPN', 'ALLN', 'olfactory']),
        'class CX outside the cut': (cls == 'CX') & ~S,
        'cut plus the rest of class CX': (cls == 'CX') | S}
for nm, m in sets.items():
    if m.sum() < 50:
        continue
    rr, _, _ = subsystem(P, pi, m, nm)
    comp.append(rr)
rng = np.random.default_rng(20261114)
rnd = []
for _ in range(20):
    m = np.zeros(n, bool); m[rng.choice(n, int(S.sum()), replace=False)] = True
    idr = np.flatnonzero(m)
    stay = np.asarray(P[idr][:, idr].sum(axis=1)).ravel()
    rnd.append(float((pi[idr] * (1 - stay)).sum() / pi[idr].sum()))
res['w1']['comparison_sets'] = comp
res['w1']['random_sets_of_same_size_escape'] = {'min': min(rnd), 'max': max(rnd), 'mean': float(np.mean(rnd))}
log(f'20 random sets of 727 neurons: escape per step in [{min(rnd):.4f}, {max(rnd):.4f}]')
json.dump(res, open(OUT, 'w'), indent=1, default=_json_default)

# ----------------------------------------------------------------------------- T4b, weight >= 1
K = 200
src = np.sort(rng.choice(idx, K, replace=False))
X = np.zeros((K, n)); X[np.arange(K), src] = 1.0
PT = P.T.tocsr()
t4b, done, t0 = [], 0, time.time()
for r in (1, 2, 4, 8, 16, 32):
    for _ in range(r - done):
        X = (PT @ X.T).T
    done = r
    Xin = np.ascontiguousarray(X[:, idx])
    din_max = dtot_max = dout_max = 0.0
    din_sum = dtot_sum = 0.0
    for i in range(K - 1):
        dt = 0.5 * np.abs(X[i][None, :] - X[i + 1:]).sum(axis=1)
        di = 0.5 * np.abs(Xin[i][None, :] - Xin[i + 1:]).sum(axis=1)
        dtot_max, din_max, dout_max = max(dtot_max, dt.max()), max(din_max, di.max()), max(dout_max, (dt - di).max())
        dtot_sum += dt.sum(); din_sum += di.sum()
    npairs = K * (K - 1) / 2
    row = {'r': r, 'mass_still_inside_mean': float(Xin.sum(axis=1).mean()),
           'tv_total_max': float(dtot_max), 'tv_inside_max': float(din_max), 'tv_outside_max': float(dout_max),
           'tv_total_mean': float(dtot_sum / npairs), 'tv_inside_mean': float(din_sum / npairs),
           'tv_outside_mean': float((dtot_sum - din_sum) / npairs)}
    t4b.append(row)
    log(f'T4b r={r:>2d}: inside mass {row["mass_still_inside_mean"]:.3f}; TV max total {dtot_max:.4f} inside '
        f'{din_max:.4f} outside {dout_max:.4f}; mean total {row["tv_total_mean"]:.4f} inside '
        f'{row["tv_inside_mean"]:.4f} outside {row["tv_outside_mean"]:.4f}  ({time.time() - t0:.0f} s)')
res['w1']['t4b_inside_outside'] = {'starts': K, 'seed': 20261114, 'rows': t4b}
del X, PT
json.dump(res, open(OUT, 'w'), indent=1, default=_json_default)

# ----------------------------------------------------------------------------- weight >= 5
keep5, P5, pi5 = chain(5)
S5 = in_cut[keep5]
r5, _, _ = subsystem(P5, pi5, S5, 'AVP/CX cut in the w >= 5 chain')
typ5 = ann['type'][keep5]
avp5 = np.array([t.startswith('MeTu') or t.startswith('TuBu') for t in typ5])
r5core, _, _ = subsystem(P5, pi5, S5 & ~avp5, 'core alone, w >= 5')
res['w5'] = {'n': int(len(keep5)), 'cut': r5, 'core': r5core}
res['wall_s'] = round(time.time() - t00, 1)
json.dump(res, open(OUT, 'w'), indent=1, default=_json_default)
log(f'saved {OUT} ({res["wall_s"]} s)')
