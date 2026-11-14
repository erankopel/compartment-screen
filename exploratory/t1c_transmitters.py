"""T1c. EXPLORATORY. The chemical-identity row of Table 5 (BioSystems 270, 106004) looked at on the
727-neuron ellipsoid-body set of T3: which predicted transmitters carry the flow into the set, inside
it and out of it, against the whole-CNS baseline. The chemical identity of a signal is taken as the
consensus transmitter prediction of the presynaptic neuron (body-nt.feather of male CNS v1.0).
Neuropeptides and co-transmission are not in these predictions. ASCII, no em dashes."""
import json, os, sys
import numpy as np
REPO = os.environ['MALECNS_REPO']
sys.path.insert(0, os.path.join(REPO, 'code', 'confirmatory'))
from common import load_edges, giant_scc, row_stochastic, stationary, annotations, consensus_nt, log, _json_default  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
pre, post, w, body_ids = load_edges()
N = len(body_ids)
e1 = json.load(open(os.path.join(REPO, 'exploratory', 'results', 'e1_slow_mode.json')))
in_cut = np.isin(body_ids, np.array(e1['w1']['cut_body_ids'], dtype=np.int64))
keep, A = giant_scc(pre, post, w, N, thresh=1)
P = row_stochastic(A)
pi, it, dl = stationary(P, tol=1e-13, maxit=20000)
S = in_cut[keep]
nt = consensus_nt(body_ids)[keep].astype(str)
typ = annotations(body_ids, cols=('type',))['type'][keep]
n = len(keep)
rows = np.repeat(np.arange(n), np.diff(P.indptr)); cols = P.indices
f = pi[rows] * P.data
Sr, Sc = S[rows], S[cols]
AMINES = ['dopamine', 'serotonin', 'octopamine', 'histamine']


nt_names, nt_code = np.unique(nt, return_inverse=True)       # integer codes: strings per edge would not fit in memory
typ_names, typ_code = np.unique(typ.astype(str), return_inverse=True)
is_amine_node = np.isin(nt, AMINES)


def comp(mask):
    s = np.bincount(nt_code[rows[mask]], weights=f[mask], minlength=len(nt_names)); s = s / s.sum()
    d = {str(nt_names[i]): float(s[i]) for i in np.argsort(-s) if s[i] > 0}
    p = s[s > 0]
    return {'shares': d, 'entropy_bits': float(-(p * np.log2(p)).sum()),
            'classes_above_1_percent': int((s > 0.01).sum()),
            'aminergic_share': float(sum(d.get(a, 0.0) for a in AMINES))}


res = {'note': 'EXPLORATORY; flow-weighted consensus transmitter of the presynaptic neuron; weight >= 1 chain',
       'whole_CNS': comp(np.ones(len(f), bool)),
       'into_the_set': comp(~Sr & Sc),
       'inside_the_set': comp(Sr & Sc),
       'out_of_the_set': comp(Sr & ~Sc)}
# neuron counts by transmitter inside the set, and the aminergic types of the set and of its sources
u, cnt = np.unique(nt[S], return_counts=True)
res['neurons_in_set_by_transmitter'] = {str(a): int(b) for a, b in sorted(zip(u, cnt), key=lambda x: -x[1])}
am_in = S & np.isin(nt, AMINES)
u, cnt = np.unique(np.char.add(np.char.add(typ[am_in].astype(str), ' / '), nt[am_in]), return_counts=True)
res['aminergic_types_in_set'] = {str(a): int(b) for a, b in sorted(zip(u, cnt), key=lambda x: -x[1])}
m = ~Sr & Sc & is_amine_node[rows]
src = rows[m]
s_node = np.bincount(src, weights=f[m], minlength=n) / f[~Sr & Sc].sum()      # share of inflow per source neuron
lab = np.char.add(np.char.add(typ.astype(str), ' / '), nt)
nz = np.flatnonzero(s_node > 0)
u, inv = np.unique(lab[nz], return_inverse=True)
s = np.bincount(inv, weights=s_node[nz])
res['aminergic_sources_share_of_inflow'] = {str(u[i]): float(s[i]) for i in np.argsort(-s)[:10]}
json.dump(res, open(os.path.join(HERE, 'results', 't1c_transmitters.json'), 'w'), indent=1, default=_json_default)
for k in ('whole_CNS', 'into_the_set', 'inside_the_set', 'out_of_the_set'):
    r = res[k]
    log(f'{k:<16s} entropy {r["entropy_bits"]:.2f} bits, classes > 1%: {r["classes_above_1_percent"]}, aminergic '
        f'{r["aminergic_share"]:.3f}; ' + ', '.join(f'{a} {b:.3f}' for a, b in list(r['shares'].items())[:7]))
log('neurons in set by transmitter: ' + str(res['neurons_in_set_by_transmitter']))
log('aminergic types in set: ' + str(res['aminergic_types_in_set']))
log('aminergic sources (share of inflow): ' + str({k: round(v, 4) for k, v in res['aminergic_sources_share_of_inflow'].items()}))

# comparison sets: aminergic share of the flow into, inside and out of each
cls = annotations(body_ids, cols=('class',))['class'][keep]
cmp_sets = {'ellipsoid-body set': S, 'class CX (all)': cls == 'CX', 'class CX outside the set': (cls == 'CX') & ~S,
            'mushroom body (KC, MBON, DAN)': np.isin(cls, ['Kenyon_Cell', 'MBON', 'DAN']),
            'antennal lobe (ALPN, ALLN)': np.isin(cls, ['ALPN', 'ALLN'])}
res['comparison_sets'] = {}
for nm, M in cmp_sets.items():
    Mr, Mc = M[rows], M[cols]
    r = {'into': comp(~Mr & Mc), 'inside': comp(Mr & Mc), 'out': comp(Mr & ~Mc)}
    res['comparison_sets'][nm] = {k: {'aminergic_share': v['aminergic_share'], 'entropy_bits': v['entropy_bits'],
                                      'top': list(v['shares'].items())[:4]} for k, v in r.items()}
    log(f'{nm:<32s} aminergic share into {r["into"]["aminergic_share"]:.3f} inside {r["inside"]["aminergic_share"]:.3f} '
        f'out {r["out"]["aminergic_share"]:.3f}; entropy bits {r["into"]["entropy_bits"]:.2f} / {r["inside"]["entropy_bits"]:.2f} / {r["out"]["entropy_bits"]:.2f}')
json.dump(res, open(os.path.join(HERE, 'results', 't1c_transmitters.json'), 'w'), indent=1, default=_json_default)
