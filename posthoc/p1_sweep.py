"""P1. POST HOC, EXPLORATORY (after the confirmatory run of 14 Nov 2026; not confirmatory, no thresholds scored).

Threshold sweep: the EB statistics of the registered screen at synapse-count thresholds 1 to 20, in the three
data sets, with the share of edges and synapses kept and the strongly connected classes of the EB-restricted
graph. Statistics from the frozen screen_lib.py.
Usage: python3 p1_sweep.py <banc | male | flywire> [<output tag> <threshold> ...]
The frozen set_stats divides by zero when the two leading eigenvalues of a set coincide (a set that has fallen
apart at a high threshold); such a set is recorded as an error and the sweep goes on. ASCII, no em dashes.
"""
import json, os, sys, time
import numpy as np
from scipy.sparse.csgraph import connected_components
import load3
from load3 import sl

name = sys.argv[1]
THS = [int(x) for x in sys.argv[3:]] or [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]
TAG = sys.argv[2] if len(sys.argv) > 2 else ''
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', f'p1_sweep_{name}{TAG}.json')


def log(s):
    print(s, flush=True)


d = load3.load(name)
pre, post, w, N, SETS, fam = d['pre'], d['post'], d['w'], d['N'], d['sets'], d['fam']
EB = SETS['EB']
inEB = EB[pre] & EB[post]
tot_e, tot_s, eb_e, eb_s = len(w), int(w.sum()), int(inEB.sum()), int(w[inEB].sum())
log(f'{name}: N {N:,}; edges {tot_e:,}; synapses {tot_s:,}; inside EB: edges {eb_e:,}, synapses {eb_s:,}')
res = {'note': 'POST HOC, EXPLORATORY', 'data': name, 'N': N, 'edges': tot_e, 'synapses': tot_s,
       'EB_edges': eb_e, 'EB_synapses': eb_s, 'runs': []}


def classes(ch, M):
    """Strongly connected classes of the graph restricted to M, with their Perron roots (classes of one
    neuron have root 0 unless they carry a self loop)."""
    idx = np.flatnonzero(M)
    PSS = ch['P'][idx][:, idx].tocsr()
    nc, lab = connected_components(PSS, directed=True, connection='strong')
    size = np.bincount(lab)
    out = []
    for c in np.argsort(-size)[:12]:
        mem = np.flatnonzero(lab == c)
        if len(mem) == 1:
            root = float(PSS[mem[0], mem[0]])
        else:
            root = float(np.max(np.linalg.eigvals(PSS[mem][:, mem].toarray()).real))
        f, n = np.unique(fam[ch['keep']][idx[mem]], return_counts=True)
        out.append({'size': int(len(mem)), 'root': root, 'families': dict(zip(f.tolist(), n.tolist()))})
    out.sort(key=lambda x: -x['root'])
    ev = np.linalg.eigvals(PSS.toarray())
    ev = ev[np.argsort(-np.abs(ev))][:6]
    return {'n_classes': int(nc), 'largest': int(size.max()), 'singletons': int((size == 1).sum()),
            'top_by_root': out[:5], 'top_eigenvalues': [[float(x.real), float(x.imag)] for x in ev]}


for th in THS:
    t0 = time.time()
    m = w >= th
    maxit = 5000 if (name != 'flywire' or th in (1, 5)) else 1500
    ch = sl.build_chain(pre[m], post[m], w[m], N, maxit=maxit)
    k = ch['keep']
    r = {'threshold': th, 'scc': ch['n'], 'iters': ch['iters'], 'delta': ch['delta'],
         'kept_edges': float(m.mean()), 'kept_synapses': float(w[m].sum() / tot_s),
         'kept_EB_edges': float((m & inEB).sum() / eb_e), 'kept_EB_synapses': float(w[m & inEB].sum() / eb_s),
         'in_scc': {a: int(b[k].sum()) for a, b in SETS.items()}, 'sets': {}}
    for nm in ('EB', 'CX', 'MB', 'AL'):
        try:
            r['sets'][nm] = sl.set_stats(ch, SETS[nm][k])
        except ZeroDivisionError:
            r['sets'][nm] = {'error': 'leading eigenvalues coincide', 'ratio': float('nan'), 'escape': float('nan'),
                             'stays': float('nan'), 'mixes': float('nan')}
    r['dobrushin_EB'] = sl.dobrushin_conditioned(ch, EB[k])
    r['exit_EB'] = sl.exit_channel(ch, EB[k])
    r['classes_EB'] = classes(ch, EB[k])
    r['wall_s'] = round(time.time() - t0)
    s, db, ex, cl = r['sets']['EB'], r['dobrushin_EB'], r['exit_EB'], r['classes_EB']
    log(f'th {th:>2}: scc {ch["n"]:,} EB {r["in_scc"]["EB"]} kept syn {r["kept_synapses"]:.3f} EBsyn {r["kept_EB_synapses"]:.3f} | '
        f'loss {s["escape"]:.4f} stays {s["stays"]:.2f} mixes {s["mixes"]:.2f} ratio {s["ratio"]:.3f} | '
        f'Dob4 {db["4"]:.3f} Dob8 {db["8"]:.3f} Dob16 {db["16"]:.3f} drop {db["dropped"]} | exit mean {ex["mean"]:.3f} '
        f'below {ex["share_below_half"]:.3f} | classes {cl["n_classes"]} largest {cl["largest"]} | '
        f'CX {r["sets"]["CX"]["ratio"]:.2f} MB {r["sets"]["MB"]["ratio"]:.2f} AL {r["sets"]["AL"]["ratio"]:.2f} | {r["wall_s"]} s')
    res['runs'].append(r)
    json.dump(res, open(OUT, 'w'), indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x))
log('saved ' + OUT)
