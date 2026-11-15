"""P12. POST HOC, EXPLORATORY (step 2: measurement or biology; not confirmatory).

Two small tables used by P9 and by the report:
  a. per EB family: synapses made and received per neuron in the three data sets; for BANC both over identified
     neurons (the edge list) and over all partners, identified or not (from P7)  -> results/p8_family_totals.json
  b. per extrinsic type (EL, ExR1, ExR2, ExR3, ExR4, ExR6): synapses per neuron with the ring neurons and with the
     columnar neurons (EPG, PEN, PEG), in both directions, in the three data sets -> results/p11_extrinsic_types.json
Run after p7_synapse_table.py and before p9_types.py. ASCII, no em dashes.
"""
import json, os, re
import numpy as np
import load3
from load3 import FAM_NAMES

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')
z = np.load(os.path.join(R, 'p7_neuron_counts.npz'))
out, inc = z['out'][:, 2, :].astype(float), z['inc'][:, 2, :].astype(float)
XT = ['EL', 'ExR1', 'ExR2', 'ExR3', 'ExR4', 'ExR6']
NT = {'EL': 'octopamine', 'ExR1': 'acetylcholine', 'ExR2': 'dopamine', 'ExR3': 'serotonin', 'ExR4': 'glutamate', 'ExR6': 'glutamate'}   # BANC neurotransmitter_verified


def grp(t):
    for x in XT:
        if any(re.match('^' + x + '(?![0-9])', p.strip()) for p in t.split(',')):
            return x
    return ''


tab, ext = {}, {}
for name in ('banc', 'male', 'flywire'):
    d = load3.load(name)
    N, fam, pre, post, w = d['N'], d['fam'], d['pre'], d['post'], d['w']
    o = np.bincount(pre, weights=w, minlength=N); i = np.bincount(post, weights=w, minlength=N)
    for f in FAM_NAMES:
        mk = fam == f
        tab.setdefault(f, {})[name] = [int(mk.sum()), float(o[mk].mean()), float(i[mk].mean())]
        if name == 'banc':
            tab[f]['banc_all'] = [float(out[mk].sum(axis=1).mean()), float(inc[mk].sum(axis=1).mean())]
    tab.setdefault('_whole', {})[name] = [N, float(o.mean()), float(i.mean())]
    if name == 'banc':
        tab['_whole']['banc_all'] = [float(out.sum(axis=1).mean()), float(inc.sum(axis=1).mean())]
    g = np.array([grp(t) for t in d['ctype']]); g[fam == ''] = ''
    ER, COL = fam == 'ER', np.isin(fam, ['EPG', 'PEN', 'PEG'])
    ext[name] = {}
    for x in XT:
        X = g == x
        n = int(X.sum())
        if n:
            f_ = lambda a, b: float(w[a[pre] & b[post]].sum()) / n
            ext[name][x] = {'n': n, 'ER_to_X': f_(ER, X), 'X_to_ER': f_(X, ER), 'COL_to_X': f_(COL, X), 'X_to_COL': f_(X, COL)}
    del d
json.dump(tab, open(os.path.join(R, 'p8_family_totals.json'), 'w'), indent=1)
json.dump({'transmitter': NT, 'data': ext}, open(os.path.join(R, 'p11_extrinsic_types.json'), 'w'), indent=1)
print('a. per neuron means: made (BANC identified partners | BANC all partners | male | FlyWire), ratios BANC all : male, BANC all : FlyWire; then received, same order')
for f in ['_whole'] + FAM_NAMES:
    b, ba, m, fw = tab[f]['banc'], tab[f]['banc_all'], tab[f]['male'], tab[f]['flywire']
    print(f'  {f:<8} made {b[1]:7.0f} {ba[0]:7.0f} {m[1]:7.0f} {fw[1]:7.0f}  {ba[0] / m[1]:.2f} {ba[0] / fw[1]:.2f} | received {b[2]:7.0f} {ba[1]:7.0f} {m[2]:7.0f} {fw[2]:7.0f}  {ba[1] / m[2]:.2f} {ba[1] / fw[2]:.2f}')
print('b. synapses per neuron of the extrinsic type X (BANC / male / FlyWire), ratios BANC : male and BANC : FlyWire')
for x in XT:
    print(f'  {x:<5} {NT[x]:<13} n {ext["banc"][x]["n"]}/{ext["male"][x]["n"]}/{ext["flywire"][x]["n"]}')
    for k, lab in (('ER_to_X', 'ring neurons onto X'), ('X_to_ER', 'X onto ring neurons'), ('COL_to_X', 'EPG, PEN, PEG onto X'), ('X_to_COL', 'X onto EPG, PEN, PEG')):
        b, m, f = (ext[n][x][k] for n in ('banc', 'male', 'flywire'))
        print(f'       {lab:<22} {b:8.0f} {m:8.0f} {f:8.0f}   {b / m:.2f}   {b / f:.2f}')
