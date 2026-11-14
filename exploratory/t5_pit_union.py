"""T5. EXPLORATORY. A first structural look at the anatomical claim of phenomenal interface theory
(Klein and Barron 2025, Phil. Trans. R. Soc. B 380, 20240301): that the central complex and the
mushroom bodies, with their aminergic inputs, work together as one system. Asked here in the terms
of T3: is the union a compartment, and how directly are the two parts coupled?
Male CNS v1.0, weight >= 1 chain, helpers of the maleCNS-depth repository unchanged.
ASCII, no em dashes."""
import json, os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
REPO = os.environ['MALECNS_REPO']
sys.path.insert(0, os.path.join(REPO, 'code', 'confirmatory'))
from common import load_edges, giant_scc, row_stochastic, stationary, annotations, log, _json_default  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
pre, post, w, body_ids = load_edges()
N = len(body_ids)
e1 = json.load(open(os.path.join(REPO, 'exploratory', 'results', 'e1_slow_mode.json')))
in_cut = np.isin(body_ids, np.array(e1['w1']['cut_body_ids'], dtype=np.int64))
keep, A = giant_scc(pre, post, w, N, thresh=1)
P = row_stochastic(A)
pi, it, dl = stationary(P, tol=1e-13, maxit=20000)
n = len(keep)
cls = annotations(body_ids, cols=('class',))['class'][keep]
S0 = in_cut[keep]
CX = (cls == 'CX') | S0
MB = np.isin(cls, ['Kenyon_Cell', 'MBON', 'DAN'])
AL = np.isin(cls, ['ALPN', 'ALLN'])
res = {'note': 'EXPLORATORY, not pre-registered'}


def subsystem(M, name):
    idx = np.flatnonzero(M)
    PSS = P[idx][:, idx].tocsr()
    stay = np.asarray(PSS.sum(axis=1)).ravel()
    mass = float(pi[idx].sum()); F = float((pi[idx] * (1 - stay)).sum())
    ev = sla.eigs(PSS, k=6, which='LM', tol=1e-9, maxiter=20000, ncv=60, return_eigenvectors=False,
                  rng=np.random.default_rng(11))
    ev = ev[np.argsort(-np.abs(ev))]
    l1, l2 = float(np.real(ev[0])), float(np.abs(ev[1]))
    r = {'size': int(len(idx)), 'mass': mass, 'escape': F / mass, 'stays': 1 / (1 - l1), 'mixes': 1 / (1 - l2 / l1),
         'ratio': (1 / (1 - l1)) / (1 / (1 - l2 / l1))}
    log(f'{name:<40s} n {len(idx):>5,} mass {mass:.4f} escape {F / mass:.4f} stays {r["stays"]:.2f} mixes {r["mixes"]:.2f} ratio {r["ratio"]:.2f}')
    return r


res['sets'] = {nm: subsystem(M, nm) for nm, M in (('central complex (class CX and the cut)', CX), ('mushroom body', MB),
                                                 ('union of central complex and mushroom body', CX | MB),
                                                 ('union of central complex and antennal lobe', CX | AL))}

# direct flow between the parts
lab = np.zeros(n, int); lab[CX] = 1; lab[MB] = 2
rows = np.repeat(np.arange(n), np.diff(P.indptr)); cols = P.indices
f = pi[rows] * P.data
flow = np.zeros((3, 3)); np.add.at(flow, (lab[rows], lab[cols]), f)
out_cx, out_mb = flow[1, 0] + flow[1, 2], flow[2, 0] + flow[2, 1]
res['direct_flow'] = {'matrix_rows_from_cols_to_[outside,CX,MB]': flow.tolist(),
                      'share_of_CX_outflow_to_MB': float(flow[1, 2] / out_cx),
                      'share_of_MB_outflow_to_CX': float(flow[2, 1] / out_mb),
                      'mass_CX': float(pi[CX].sum()), 'mass_MB': float(pi[MB].sum())}
log(f'direct: share of CX outflow that lands in MB {flow[1, 2] / out_cx:.4f} (MB mass {pi[MB].sum():.4f}); '
    f'share of MB outflow that lands in CX {flow[2, 1] / out_mb:.4f} (CX mass {pi[CX].sum():.4f})')


# indirect coupling: a walker leaves the source set; probability that it reaches the target within r steps
def reach(source, target, R=16):
    Sr, Sc = source[rows], source[cols]
    m = Sr & ~Sc
    d = np.bincount(cols[m], weights=f[m], minlength=n); d /= d.sum()      # where walkers land on leaving
    PT = P.T.tocsr()
    hit, cum = [], 0.0
    for r in range(1, R + 1):
        a = float(d[target].sum()); cum += a; hit.append(cum)
        d = d.copy(); d[target] = 0.0
        d = PT @ d
    return hit


def reach_stationary(target, R=16):
    d = pi.copy(); d[target] = 0; d /= d.sum()
    PT = P.T.tocsr(); hit, cum = [], 0.0
    for r in range(1, R + 1):
        d = PT @ d
        a = float(d[target].sum()); cum += a; hit.append(cum)
        d[target] = 0.0
    return hit


h_mb_cx, h_cx_mb = reach(MB, CX), reach(CX, MB)
b_cx, b_mb = reach_stationary(CX), reach_stationary(MB)
h_al_cx = reach(AL, CX)
res['reach_within_r_steps'] = {'MB_to_CX': h_mb_cx, 'CX_to_MB': h_cx_mb, 'AL_to_CX': h_al_cx,
                               'stationary_walker_to_CX': b_cx, 'stationary_walker_to_MB': b_mb}
for nm, h, b in (('MB -> CX', h_mb_cx, b_cx), ('AL -> CX', h_al_cx, b_cx), ('CX -> MB', h_cx_mb, b_mb)):
    log(f'{nm}: reached within 1, 2, 4, 8, 16 steps: ' + ', '.join(f'{h[k - 1]:.4f}' for k in (1, 2, 4, 8, 16)) +
        ' | stationary walker: ' + ', '.join(f'{b[k - 1]:.4f}' for k in (1, 2, 4, 8, 16)))
json.dump(res, open(os.path.join(HERE, 'results', 't5_pit_union.json'), 'w'), indent=1, default=_json_default)
log('saved')
