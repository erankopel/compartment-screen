"""P3. POST HOC, EXPLORATORY (after the confirmatory run of 14 Nov 2026; not confirmatory).

Anatomy of the EB set at a given threshold: which neurons carry the slow mode (stays / mixes) and which
starting neurons have not mixed after 8 steps (Dobrushin coefficient). Also the same statistics for two
diagnostic subsets: the core (ER, EPG, PEG, PEN, EL, ExR, LPsP, GLNO) and the visual input tail (TuBu, MeTu2,
TuTuB, AOTU046). The subsets are a diagnostic, not a new definition of the set.
Usage: python3 p3_anatomy.py <banc | male | flywire> <threshold> [<threshold> ...]
ASCII, no em dashes.
"""
import json, os, sys
import numpy as np
from scipy.sparse.csgraph import connected_components
import load3
from load3 import sl, conditioned_Q

name = sys.argv[1]
THS = [int(x) for x in sys.argv[2:]]
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', f'p3_anatomy_{name}.json')
CORE = ['ER', 'EPG', 'PEG', 'PEN', 'EL', 'ExR', 'LPsP', 'GLNO']
TAIL = ['TuBu', 'MeTu2', 'TuTuB', 'AOTU046']


def log(s):
    print(s, flush=True)


def comp(labels):
    u, c = np.unique(labels, return_counts=True)
    return dict(sorted(zip(u.tolist(), c.tolist()), key=lambda t: -t[1]))


d = load3.load(name)
pre, post, w, N, SETS, fam, side, ctype = d['pre'], d['post'], d['w'], d['N'], d['sets'], d['fam'], d['side'], d['ctype']
side = np.array([s[:1].upper() for s in side])
EB = SETS['EB']
res = {'note': 'POST HOC, EXPLORATORY', 'data': name, 'runs': []}

for th in THS:
    m = w >= th
    ch = sl.build_chain(pre[m], post[m], w[m], N, maxit=600)
    k = ch['keep']
    M = EB[k]
    r = {'threshold': th, 'EB_in_scc': int(M.sum())}
    log(f'=== {name}, threshold {th}: EB neurons in the component {int(M.sum())}')

    # 1. slow mode of the conditioned walk
    Q, kept, l1 = conditioned_Q(ch, M)
    f_, s_, t_ = fam[k][kept], side[k][kept], ctype[k][kept]
    n = len(kept)
    ev, V = np.linalg.eig(Q)
    o = np.argsort(-np.abs(ev))
    ev, V = ev[o], V[:, o]
    evT, U = np.linalg.eig(Q.T)
    nu = np.abs(U[:, np.argmax(evT.real)].real); nu /= nu.sum()      # stationary law of the conditioned walk
    r['lam1'] = l1
    r['eigen_Q'] = [[float(x.real), float(x.imag)] for x in ev[:6]]
    v2 = V[:, 1].real
    v2 = v2 / np.abs(v2).max()
    pos, neg = v2 > 0.1, v2 < -0.1
    if nu[pos].sum() > nu[neg].sum():
        pos, neg = neg, pos                          # 'minor' = the side with less stationary mass
    r['slow_mode'] = {'eigenvalue': [float(ev[1].real), float(ev[1].imag)],
                      'minor_side': {'neurons': int(pos.sum()), 'mass': float(nu[pos].sum()), 'families': comp(f_[pos]),
                                     'sides': comp(s_[pos]), 'types': dict(list(comp(t_[pos]).items())[:12])},
                      'major_side': {'neurons': int(neg.sum()), 'mass': float(nu[neg].sum()), 'families': comp(f_[neg]),
                                     'sides': comp(s_[neg])},
                      'near_zero': int((~pos & ~neg).sum())}
    sm = r['slow_mode']
    sgn = -1.0 if (v2[pos].mean() if pos.any() else 0) < 0 else 1.0      # orient: the minor side is positive
    tc = {}
    for t in np.unique(t_):
        mk = t_ == t
        tc[t] = [int(mk.sum()), float(sgn * v2[mk].mean()), float(nu[mk].sum())]
    sm['type_coordinate_n_mean_mass'] = dict(sorted(tc.items(), key=lambda kv: -kv[1][1]))
    sm['major_side']['types'] = dict(list(comp(t_[neg]).items())[:12])
    log(f'  eigenvalues of the conditioned walk: ' + ', '.join(f'{x.real:+.4f}{x.imag:+.4f}i' for x in ev[:5]))
    log(f'  slow mode, minor side: {sm["minor_side"]["neurons"]} neurons, mass {sm["minor_side"]["mass"]:.4f}, families {sm["minor_side"]["families"]}, sides {sm["minor_side"]["sides"]}')
    log(f'     types {sm["minor_side"]["types"]}')
    log(f'  slow mode, major side: {sm["major_side"]["neurons"]} neurons, mass {sm["major_side"]["mass"]:.4f}, families {sm["major_side"]["families"]}')

    # 2. strongly connected classes inside EB
    idx = np.flatnonzero(M)
    PSS = ch['P'][idx][:, idx].tocsr()
    nc, lab = connected_components(PSS, directed=True, connection='strong')
    size = np.bincount(lab)
    cls = []
    for c in np.argsort(-size)[:15]:
        mem = np.flatnonzero(lab == c)
        root = float(PSS[mem[0], mem[0]]) if len(mem) == 1 else float(np.max(np.linalg.eigvals(PSS[mem][:, mem].toarray()).real))
        cls.append({'size': int(len(mem)), 'root': root, 'families': comp(fam[k][idx[mem]]), 'sides': comp(side[k][idx[mem]])})
    cls.sort(key=lambda x: -x['root'])
    single = size[lab] == 1
    r['classes'] = {'n': int(nc), 'largest': int(size.max()), 'singletons': int((size == 1).sum()),
                    'singleton_families': comp(fam[k][idx[single]]), 'top_by_root': cls[:6]}
    log(f'  classes inside EB: {nc}; largest {int(size.max())}; single neurons {int((size == 1).sum())} {r["classes"]["singleton_families"]}')
    for c in cls[:4]:
        log(f'     class size {c["size"]:>4} root {c["root"]:.4f} families {c["families"]} sides {c["sides"]}')

    # 3. who has not mixed after 8 steps
    Q8 = np.linalg.matrix_power(Q, 8)
    to_stat = 0.5 * np.abs(Q8 - nu[None, :]).sum(axis=1)
    ecc = np.zeros(n); arg = np.zeros(n, int)
    allpairs = []
    for i in range(n):
        dd = 0.5 * np.abs(Q8[i][None, :] - Q8).sum(axis=1)
        ecc[i], arg[i] = dd.max(), int(dd.argmax())
        allpairs.append(dd[i + 1:])
    allpairs = np.concatenate(allpairs)
    i0 = int(ecc.argmax()); j0 = int(arg[i0])
    # greedy: drop the start with the largest distance to the stationary law until the coefficient is below 0.6
    order = np.argsort(-to_stat)
    removed, cur = 0, float(ecc.max())
    keepmask = np.ones(n, bool)
    steps = []
    while cur >= 0.6 and removed < n - 2:
        keepmask[order[removed]] = False; removed += 1
        if removed in (1, 2, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300) or removed % 25 == 0:
            R = Q8[keepmask]
            cur = max(float((0.5 * np.abs(R[i][None, :] - R[i + 1:]).sum(axis=1)).max()) for i in range(len(R) - 1))
            steps.append([removed, cur])
    fam_tab = {}
    for f in np.unique(f_):
        mk = f_ == f
        fam_tab[f] = {'neurons': int(mk.sum()), 'median_to_stationary': float(np.median(to_stat[mk])),
                      'max_to_stationary': float(to_stat[mk].max()), 'above_0.3': int((to_stat[mk] > 0.3).sum())}
    r['unmixed_r8'] = {'dobrushin': float(ecc.max()), 'pair': [t_[i0], s_[i0], t_[j0], s_[j0]],
                       'pairwise_percentiles_50_90_99': [float(x) for x in np.percentile(allpairs, [50, 90, 99])],
                       'share_pairs_above_0.6': float((allpairs > 0.6).mean()),
                       'starts_above_0.3_to_stationary': int((to_stat > 0.3).sum()),
                       'by_family': fam_tab, 'greedy_removal_steps': steps,
                       'removed_families_first_50': comp(f_[order[:50]])}
    u = r['unmixed_r8']
    log(f'  r = 8: coefficient {u["dobrushin"]:.4f}, attained by {u["pair"]}; pairwise distances median {u["pairwise_percentiles_50_90_99"][0]:.3f}, '
        f'90% {u["pairwise_percentiles_50_90_99"][1]:.3f}, 99% {u["pairwise_percentiles_50_90_99"][2]:.3f}; pairs above 0.6: {u["share_pairs_above_0.6"]:.4f}')
    log(f'  starts farther than 0.3 from the stationary law: {u["starts_above_0.3_to_stationary"]} of {n}; greedy removal (starts removed, coefficient): {steps}')
    log('  by family (neurons, median and max distance to the stationary law, starts above 0.3): ' +
        '; '.join(f'{f} {v["neurons"]} {v["median_to_stationary"]:.2f} {v["max_to_stationary"]:.2f} {v["above_0.3"]}' for f, v in fam_tab.items()))

    # 4. diagnostic subsets
    r['subsets'] = {}
    for nm, fams in (('core', CORE), ('tail', TAIL), ('core_left', CORE), ('core_right', CORE)):
        S = EB[k] & np.isin(fam[k], fams)
        if nm.endswith('left'):
            S &= side[k] == 'L'
        if nm.endswith('right'):
            S &= side[k] == 'R'
        st = sl.set_stats(ch, S)
        db = sl.dobrushin_conditioned(ch, S) if st.get('size', 0) >= 3 else {}
        ex = sl.exit_channel(ch, S) if nm in ('core', 'tail') else {}
        r['subsets'][nm] = {'stats': st, 'dobrushin': db, 'exit': ex}
        log(f'  subset {nm:<10}: n {st["size"]:>4} loss {st["escape"]:.4f} stays {st["stays"]:.2f} mixes {st["mixes"]:.2f} ratio {st["ratio"]:.3f} | '
            f'Dobrushin r=8 {db.get("8", float("nan")):.3f} (dropped {db.get("dropped", 0)})' + (f' | exit mean {ex["mean"]:.3f} below {ex["share_below_half"]:.3f}' if ex else ''))
    res['runs'].append(r)
    json.dump(res, open(OUT, 'w'), indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x))
log('saved ' + OUT)
