"""T4b (second part). EXPLORATORY, not pre-registered. The exit channel of the 727-neuron set.

A walker starts on neuron a of the set S. An outside observer sees nothing until the walker first
leaves S, at step t onto the outside neuron j. By the strong Markov property the pair (t, j) is a
sufficient statistic of the start for everything that happens outside afterwards, so

    d_out(a, b) = TV( law of (t, j) from a , law of (t, j) from b )

is the whole difference between two internal starts that can ever reach the outside (for an observer
who knows the entry time; d_out without t is the same for one who does not). The Dobrushin
coefficient of this exit channel is max over pairs. ASCII, no em dashes.
"""
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
S = in_cut[keep]
typ = annotations(body_ids, cols=('type',))['type'][keep]
idx, out_idx = np.flatnonzero(S), np.flatnonzero(~S)
D = P[idx][:, idx].toarray()
E = P[idx][:, out_idx].tocsc()
nzc = np.flatnonzero(np.diff(E.indptr) > 0)
E = E[:, nzc].tocsr()
m = E.shape[1]
log(f'set {len(idx)}; outside neurons reached directly from the set: {m:,}; exit edges {E.nnz:,}')
tt = typ[out_idx[nzc]]
ut, inv = np.unique(tt.astype(str), return_inverse=True)
import scipy.sparse as sp
Agg = sp.csr_matrix((np.ones(m), (np.arange(m), inv)), shape=(m, len(ut)))     # neurons -> types
ET = (E @ Agg).tocsr()
log(f'exit targets grouped into {len(ut)} cell types')

is_avp = np.array([t.startswith('MeTu') or t.startswith('TuBu') for t in typ[idx]])
rng = np.random.default_rng(20261114)
K = 240
sel = np.sort(rng.choice(len(idx), K, replace=False))
core = ~is_avp[sel]
V = np.zeros((K, len(idx))); V[np.arange(K), sel] = 1.0
iu = np.triu_indices(K, 1)


def pair_l1(R):
    out = np.zeros((K, K))
    for i in range(K - 1):
        out[i, i + 1:] = np.abs(R[i][None, :] - R[i + 1:]).sum(axis=1)
    return out


TMAX = 160
d_time_neuron = np.zeros((K, K)); d_time_type = np.zeros((K, K))
Rsum = np.zeros((K, m)); RTsum = np.zeros((K, len(ut)))
surv, t0 = [], time.time()
for t in range(1, TMAX + 1):
    R = np.asarray((E.T @ V.T).T)            # K x m: probability of first exit at step t onto neuron j
    RT = np.asarray((ET.T @ V.T).T)
    d_time_neuron += pair_l1(R); d_time_type += pair_l1(RT)
    Rsum += R; RTsum += RT
    V = V @ D
    surv.append(float(V.sum(axis=1).max()))
    if t in (1, 2, 4, 8, 16, 32, 64, 128, TMAX):
        log(f't = {t:>3d}: largest mass still inside {surv[-1]:.2e} ({time.time() - t0:.0f} s)')
resid = surv[-1]
d_time_neuron *= 0.5; d_time_type *= 0.5
d_neuron = 0.5 * pair_l1(Rsum); d_type = 0.5 * pair_l1(RTsum)


def summary(Dm, mask=None):
    x = Dm[iu]
    if mask is not None:
        x = x[mask]
    return {'pairs': int(len(x)), 'max': float(x.max()), 'mean': float(x.mean()), 'median': float(np.median(x)),
            'p10': float(np.percentile(x, 10)), 'share_below_0.5': float((x < 0.5).mean()),
            'share_below_0.25': float((x < 0.25).mean()), 'share_above_0.9': float((x > 0.9).mean())}


cc = core[iu[0]] & core[iu[1]]
res = {'note': 'EXPLORATORY; exit channel of the E1 cut in the weight >= 1 chain; lower bounds on the maximum '
               'over all pairs because only K starts are used',
       'K': K, 'seed': 20261114, 'core_starts': int(core.sum()), 'TMAX': TMAX,
       'truncation_error_bound_on_each_TV': resid, 'outside_neurons_reached': int(m), 'target_types': int(len(ut))}
for nm, Dm in (('time_and_neuron', d_time_neuron), ('neuron_only', d_neuron),
               ('time_and_type', d_time_type), ('type_only', d_type)):
    res[nm] = {'all_pairs': summary(Dm), 'core_core_pairs': summary(Dm, cc)}
    a, c = res[nm]['all_pairs'], res[nm]['core_core_pairs']
    log(f'{nm:<16s} all pairs: max {a["max"]:.3f} mean {a["mean"]:.3f} median {a["median"]:.3f} '
        f'<0.5: {a["share_below_0.5"]:.3f} >0.9: {a["share_above_0.9"]:.3f} | core-core: max {c["max"]:.3f} '
        f'mean {c["mean"]:.3f} median {c["median"]:.3f} <0.5: {c["share_below_0.5"]:.3f}')

# same cell type versus different cell type (core-core pairs), time and neuron
ts = typ[idx][sel]
same = (ts[iu[0]] == ts[iu[1]])
x = d_time_neuron[iu]
res['time_and_neuron']['same_type_pairs'] = {'pairs': int(same.sum()), 'mean': float(x[same].mean()),
                                             'median': float(np.median(x[same]))}
res['time_and_neuron']['different_type_pairs'] = {'pairs': int((~same).sum()), 'mean': float(x[~same].mean()),
                                                  'median': float(np.median(x[~same]))}
log(f'same type pairs: {int(same.sum())}, mean {x[same].mean():.3f}; different type: mean {x[~same].mean():.3f}')
json.dump(res, open(os.path.join(HERE, 'results', 't4b_exit_channel.json'), 'w'), indent=1, default=_json_default)
log('saved')
