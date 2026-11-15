"""P13. POST HOC, EXPLORATORY (step 2; not confirmatory).

The table of p12 (synapses per neuron of the extrinsic types EL and ExR1 with the ring neurons and with the
columnar neurons) on the BANC v2 edge list, the synapse set used by the BANC paper.
Usage: BANC_EDGE_FILE=<banc_888_edgelist_simple_v2.feather> python3 p13_v2_extrinsic.py     ASCII, no em dashes.
"""
import json, os, re, sys
import numpy as np
import load3

if 'BANC_EDGE_FILE' not in os.environ:
    sys.exit('set BANC_EDGE_FILE to the v2 edge list')
HERE = os.path.dirname(os.path.abspath(__file__))
d = load3.load('banc')
fam, ct, pre, post, w = d['fam'], d['ctype'], d['pre'], d['post'], d['w']
EB = d['sets']['EB']
ER, COL = fam == 'ER', np.isin(fam, ['EPG', 'PEN', 'PEG'])
res = {'note': 'POST HOC, EXPLORATORY', 'edge_list': os.path.basename(os.environ['BANC_EDGE_FILE']),
       'synapses_in_graph': int(w.sum()), 'synapses_inside_EB': int(w[EB[pre] & EB[post]].sum()), 'types': {}}
print(f'{res["edge_list"]}: synapses in the graph {res["synapses_in_graph"]:,}; inside the EB set {res["synapses_inside_EB"]:,}')
for x in ('EL', 'ExR1', 'ExR2', 'ExR3', 'ExR4', 'ExR6'):
    X = np.array([any(re.match('^' + x + '(?![0-9])', p.strip()) for p in t.split(',')) for t in ct]) & (fam != '')
    n = int(X.sum())
    f = lambda a, b: float(w[a[pre] & b[post]].sum()) / n
    res['types'][x] = {'n': n, 'ER_to_X': f(ER, X), 'X_to_ER': f(X, ER), 'COL_to_X': f(COL, X), 'X_to_COL': f(X, COL)}
    r = res['types'][x]
    print(f'  {x:<5} n {n}: ring neurons onto X {r["ER_to_X"]:.0f}; X onto ring neurons {r["X_to_ER"]:.0f}; EPG, PEN, PEG onto X {r["COL_to_X"]:.0f}; X onto EPG, PEN, PEG {r["X_to_COL"]:.0f}')
json.dump(res, open(os.path.join(HERE, 'results', 'p13_v2_extrinsic.json'), 'w'), indent=1)
