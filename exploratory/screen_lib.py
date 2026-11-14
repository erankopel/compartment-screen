"""Compartment screen for a directed weighted connectome: one function per statistic, used unchanged for
every data set (male CNS, its null models, FlyWire, and later the confirmatory data set).

Chain: random walk on the giant strongly connected component, P = row-normalised synapse counts
(rows = presynaptic), pi = stationary distribution. ASCII, no em dashes.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components

SEED = 20261114


def build_chain(pre, post, w, n, tol=1e-12, maxit=5000):
    """Power iteration from the uniform vector, L1 stopping rule; if the tolerance is not met after maxit
    iterations the last vector is used and the last change (delta) is reported."""
    A = sp.csr_matrix((np.asarray(w, np.float64), (pre, post)), shape=(n, n))
    A.sum_duplicates()
    _, lab = connected_components(A, directed=True, connection='strong')
    keep = np.flatnonzero(lab == np.bincount(lab).argmax())
    A = A[keep][:, keep].tocsr()
    rs = np.asarray(A.sum(axis=1)).ravel()
    P = (sp.diags(1.0 / rs) @ A).tocsr()
    PT = P.T.tocsr()
    v = np.full(len(keep), 1.0 / len(keep))
    for it in range(1, maxit + 1):
        vn = PT @ v
        vn /= vn.sum()
        d = np.abs(vn - v).sum()
        v = vn
        if d < tol:
            break
    rows = np.repeat(np.arange(len(keep)), np.diff(P.indptr))
    f = v[rows] * P.data
    return {'keep': keep, 'P': P, 'PT': PT, 'pi': v, 'iters': it, 'delta': float(d), 'rows': rows, 'cols': P.indices, 'f': f,
            'n': len(keep)}


def set_stats(ch, M):
    """Escape per step, residence 1/(1 - lam1) of the walk stopped on leaving M, internal relaxation
    1/(1 - |lam2|/lam1), and their ratio."""
    P, pi = ch['P'], ch['pi']
    idx = np.flatnonzero(M)
    if len(idx) < 3:
        return {'size': int(len(idx))}
    PSS = P[idx][:, idx].tocsr()
    stay = np.asarray(PSS.sum(axis=1)).ravel()
    mass = float(pi[idx].sum())
    F = float((pi[idx] * (1 - stay)).sum())
    if len(idx) <= 1600:
        ev = np.linalg.eigvals(PSS.toarray())
    else:
        ev = sla.eigs(PSS, k=6, which='LM', tol=1e-9, maxiter=20000, ncv=60, return_eigenvectors=False,
                      rng=np.random.default_rng(11))
    ev = ev[np.argsort(-np.abs(ev))]
    l1, l2 = float(np.real(ev[0])), float(np.abs(ev[1]))
    return {'size': int(len(idx)), 'mass': mass, 'escape': F / mass, 'lam1': l1, 'lam2': l2,
            'stays': 1 / (1 - l1), 'mixes': 1 / (1 - l2 / l1), 'ratio': (1 - l2 / l1) / (1 - l1)}


def escape_only(ch, M):
    """Share of the stationary mass of M that leaves M in one step."""
    rows, cols, f, pi = ch['rows'], ch['cols'], ch['f'], ch['pi']
    return float(f[M[rows] & ~M[cols]].sum() / pi[M].sum())


def dobrushin_conditioned(ch, M, rs=(4, 8, 16)):
    """Dobrushin coefficient, over all pairs, of the walk conditioned to stay in M (Doob transform with
    the right Perron vector h). Neurons with h = 0 (no way to remain) are dropped and counted."""
    idx = np.flatnonzero(M)
    D = ch['P'][idx][:, idx].toarray()
    ev, V = np.linalg.eig(D)
    k = int(np.argmax(ev.real))
    l1, h = float(ev[k].real), np.abs(V[:, k].real)
    ok = h > 1e-10 * h.max()
    D, h = D[ok][:, ok], h[ok]
    Q = D * h[None, :] / (l1 * h[:, None])
    Q /= Q.sum(axis=1, keepdims=True)
    out, Qr, done = {'dropped': int((~ok).sum())}, np.eye(len(h)), 0
    for r in rs:
        Qr = Qr @ np.linalg.matrix_power(Q, r - done)
        done = r
        m = 0.0
        for i in range(len(h) - 1):
            m = max(m, 0.5 * float(np.abs(Qr[i][None, :] - Qr[i + 1:]).sum(axis=1).max()))
        out[str(r)] = m
    return out


def exit_channel(ch, M, K=240, tol=1e-4, tcap=200):
    """Total variation distance between the laws of (first exit step, exit neuron) for pairs of starting
    neurons of M; K sampled starts."""
    P = ch['P']
    rng = np.random.default_rng(SEED)
    idx, oidx = np.flatnonzero(M), np.flatnonzero(~M)
    D = P[idx][:, idx].toarray()
    E = P[idx][:, oidx].tocsc()
    ET = E[:, np.flatnonzero(np.diff(E.indptr) > 0)].T.tocsr()
    k = min(K, len(idx))
    sel = np.sort(rng.choice(len(idx), k, replace=False))
    V = np.zeros((k, len(idx)))
    V[np.arange(k), sel] = 1.0
    d, t = np.zeros((k, k)), 0
    while True:
        t += 1
        R = np.asarray((ET @ V.T).T)
        for i in range(k - 1):
            d[i, i + 1:] += np.abs(R[i][None, :] - R[i + 1:]).sum(axis=1)
        V = V @ D
        left = float(V.sum(axis=1).max())
        if left < tol or t >= tcap:
            break
    x = 0.5 * d[np.triu_indices(k, 1)]
    return {'starts': int(k), 'steps': t, 'mass_left_bound': left, 'max': float(x.max()), 'mean': float(x.mean()),
            'median': float(np.median(x)), 'share_below_half': float((x < 0.5).mean())}


def label_shares(ch, M, flag):
    """Share of the flow carried by presynaptic neurons with flag (e.g. aminergic), for the whole chain
    and for the flow into, inside and out of M."""
    rows, cols, f = ch['rows'], ch['cols'], ch['f']
    a, Mr, Mc = flag[rows], M[rows], M[cols]

    def sh(m):
        return float(f[m & a].sum() / f[m].sum())
    return {'whole': float(f[a].sum() / f.sum()), 'into': sh(~Mr & Mc), 'inside': sh(Mr & Mc), 'out': sh(Mr & ~Mc)}


def reach(ch, source, target, R=4):
    """Probability that a walker leaving `source` (stationary exit flow) has entered `target` within
    1..R steps, and the same for a stationary walker outside `target`."""
    rows, cols, f, PT, pi, n = ch['rows'], ch['cols'], ch['f'], ch['PT'], ch['pi'], ch['n']
    m = source[rows] & ~source[cols]
    d = np.bincount(cols[m], weights=f[m], minlength=n)
    d /= d.sum()
    hit, cum = [], 0.0
    for _ in range(R):
        cum += float(d[target].sum()); hit.append(cum)
        d = d.copy(); d[target] = 0.0
        d = PT @ d
    b = pi.copy(); b[target] = 0.0; b /= b.sum()
    base, cum = [], 0.0
    for _ in range(R):
        b = PT @ b
        cum += float(b[target].sum()); base.append(cum)
        b[target] = 0.0
    return {'hit': hit, 'baseline': base, 'ratio_1': hit[0] / base[0], 'ratio_R': hit[-1] / base[-1]}


def screen(ch, sets, amine, full=True, log=print):
    """sets: dict with boolean masks over ch['keep'] for 'EB', 'CX', 'MB', 'AL'. amine: boolean mask."""
    out = {'n': ch['n'], 'iters': ch['iters'], 'delta': ch['delta'], 'sets': {}}
    for nm in ('EB', 'CX', 'MB', 'AL'):
        out['sets'][nm] = set_stats(ch, sets[nm])
        s = out['sets'][nm]
        if 'escape' in s:
            log(f'  {nm}: n {s["size"]:>5} mass {s["mass"]:.4f} escape {s["escape"]:.4f} stays {s["stays"]:.2f} '
                f'mixes {s["mixes"]:.2f} ratio {s["ratio"]:.2f}')
    out['escape_union_CX_MB'] = escape_only(ch, sets['CX'] | sets['MB'])
    log(f'  escape of the union of CX and MB {out["escape_union_CX_MB"]:.4f}')
    out['dobrushin_EB'] = dobrushin_conditioned(ch, sets['EB'])
    out['exit_EB'] = exit_channel(ch, sets['EB'])
    log(f'  Dobrushin EB {out["dobrushin_EB"]}; exit EB mean {out["exit_EB"]["mean"]:.3f} max {out["exit_EB"]["max"]:.3f}')
    if full:
        out['exit_MB'] = exit_channel(ch, sets['MB'])
        out['exit_CX'] = exit_channel(ch, sets['CX'])
        log(f'  exit MB mean {out["exit_MB"]["mean"]:.3f} max {out["exit_MB"]["max"]:.3f}; exit CX mean '
            f'{out["exit_CX"]["mean"]:.3f} max {out["exit_CX"]["max"]:.3f}')
    out['amine'] = {nm: label_shares(ch, sets[nm], amine) for nm in ('EB', 'CX', 'MB')}
    a = out['amine']['EB']
    log(f'  amine share: whole {a["whole"]:.4f}; EB into {a["into"]:.3f} inside {a["inside"]:.3f} out {a["out"]:.3f}; '
        f'MB inside {out["amine"]["MB"]["inside"]:.3f}')
    out['reach'] = {'MB_to_CX': reach(ch, sets['MB'], sets['CX']), 'CX_to_MB': reach(ch, sets['CX'], sets['MB']),
                    'AL_to_CX': reach(ch, sets['AL'], sets['CX'])}
    r = out['reach']
    log(f'  reach MB->CX 1 step {r["MB_to_CX"]["hit"][0]:.4f} (x{r["MB_to_CX"]["ratio_1"]:.1f}), 4 steps '
        f'{r["MB_to_CX"]["hit"][3]:.4f} (x{r["MB_to_CX"]["ratio_R"]:.1f}); CX->MB x{r["CX_to_MB"]["ratio_1"]:.1f}, '
        f'x{r["CX_to_MB"]["ratio_R"]:.1f}; AL->CX x{r["AL_to_CX"]["ratio_1"]:.1f}, x{r["AL_to_CX"]["ratio_R"]:.1f}')
    return out
