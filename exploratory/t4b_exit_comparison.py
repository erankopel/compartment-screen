"""T4b (third part). EXPLORATORY. The exit-channel statistic of t4b_exit_channel.py (first exit step
and neuron) for comparison sets, 240 sampled starts each, so that the value for the ellipsoid-body cut
can be read against other anatomically defined sets. ASCII, no em dashes."""
import json, os, sys, time
import numpy as np
REPO = os.environ['MALECNS_REPO']
sys.path.insert(0, os.path.join(REPO, 'code', 'confirmatory'))
from common import load_edges, giant_scc, row_stochastic, annotations, log, _json_default  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
pre, post, w, body_ids = load_edges()
N = len(body_ids)
e1 = json.load(open(os.path.join(REPO, 'exploratory', 'results', 'e1_slow_mode.json')))
in_cut = np.isin(body_ids, np.array(e1['w1']['cut_body_ids'], dtype=np.int64))
keep, A = giant_scc(pre, post, w, N, thresh=1)
P = row_stochastic(A)
ann = annotations(body_ids, cols=('type', 'class'))
cls, typ = ann['class'][keep], ann['type'][keep]
S0 = in_cut[keep]
sets = {'ellipsoid-body cut (727)': S0,
        'class CX (all)': cls == 'CX',
        'class CX outside the cut': (cls == 'CX') & ~S0,
        'mushroom body (KC, MBON, DAN)': np.isin(cls, ['Kenyon_Cell', 'MBON', 'DAN']),
        'Kenyon cells': cls == 'Kenyon_Cell',
        'antennal lobe (ALPN, ALLN)': np.isin(cls, ['ALPN', 'ALLN'])}
K, TOL, TCAP = 240, 1e-4, 400
out = {'note': 'EXPLORATORY; time-and-neuron exit record; K sampled starts per set', 'K': K, 'seed': 20261114}
for name, S in sets.items():
    rng = np.random.default_rng(20261114)
    idx, oidx = np.flatnonzero(S), np.flatnonzero(~S)
    D = P[idx][:, idx].toarray()
    E = P[idx][:, oidx].tocsc()
    nzc = np.flatnonzero(np.diff(E.indptr) > 0)
    E = E[:, nzc].tocsr(); ET = E.T.tocsr()
    sel = np.sort(rng.choice(len(idx), min(K, len(idx)), replace=False))
    k = len(sel)
    V = np.zeros((k, len(idx))); V[np.arange(k), sel] = 1.0
    d = np.zeros((k, k)); t0 = time.time(); t = 0
    while True:
        t += 1
        R = np.asarray((ET @ V.T).T)
        for i in range(k - 1):
            d[i, i + 1:] += np.abs(R[i][None, :] - R[i + 1:]).sum(axis=1)
        V = V @ D
        left = float(V.sum(axis=1).max())
        if left < TOL or t >= TCAP:
            break
    x = 0.5 * d[np.triu_indices(k, 1)]
    ts = typ[idx][sel]; iu = np.triu_indices(k, 1); same = ts[iu[0]] == ts[iu[1]]
    r = {'size': int(len(idx)), 'starts': int(k), 'steps': t, 'mass_left_bound': left, 'targets': int(len(nzc)),
         'max': float(x.max()), 'mean': float(x.mean()), 'median': float(np.median(x)),
         'share_below_0.5': float((x < 0.5).mean()), 'share_above_0.9': float((x > 0.9).mean()),
         'mean_same_type': float(x[same].mean()) if same.any() else None,
         'mean_different_type': float(x[~same].mean())}
    out[name] = r
    log(f'{name:<32s} n {len(idx):>5,} steps {t:>3d} max {r["max"]:.3f} mean {r["mean"]:.3f} median {r["median"]:.3f} '
        f'<0.5: {r["share_below_0.5"]:.3f} >0.9: {r["share_above_0.9"]:.3f} same type {r["mean_same_type"]} '
        f'diff type {r["mean_different_type"]:.3f} ({time.time() - t0:.0f} s)')
    json.dump(out, open(os.path.join(HERE, 'results', 't4b_exit_comparison.json'), 'w'), indent=1, default=_json_default)
log('saved')
