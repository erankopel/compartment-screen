"""T2c. The slice of arXiv:2609.27908 Sec. 6 on the grid of the archived figure script
(Zenodo 10.5281/zenodo.22830716, figures/make_fig1_index_landscape.py): theta/pi in [0.01, 0.25],
phi/pi in [0.01, 0.49], 81 x 81, cutoff 200, criterion lambda_min >= 0. Computed with the independent
implementation of t2_core.py (circuit built here, nothing imported from the archive).
ASCII, no em dashes."""
import importlib.util, json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('t2', os.path.join(HERE, 't2_core.py'))
t2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t2)
pi = math.pi
k0, b0 = 0.43230980 * pi, 0.23903823 * pi
ths = np.linspace(.01, .25, 81) * pi
phs = np.linspace(.01, .49, 81) * pi
out = {}
for tol in (0.0, 1e-12):
    Ns, prods = [], []
    for ph in phs:
        for th in ths:
            A, c, env, off = t2.channel_data(t2.unitary(th, k0, ph, b0))
            n1 = t2.index(A, c, cap=200, tol=tol)
            Ns.append(n1)
            if 1 < n1 <= 200:
                prods.append(n1 * (1 - float(np.abs(np.linalg.eigvals(A)).max())))
    Ns, prods = np.array(Ns), np.array(prods)
    r = {'share_index_1': float((Ns == 1).mean()), 'distinct': int(len(np.unique(Ns[Ns <= 200]))),
         'cutoff': int((Ns > 200).sum()), 'n': int(len(prods)), 'max_resolved': int(Ns[Ns <= 200].max()),
         'product_min': float(prods.min()), 'product_max': float(prods.max()), 'product_median': float(np.median(prods))}
    out[f'tol {tol:g}'] = r
    print(f'tol {tol:<6g}', r)
out['published'] = {'share_index_1': 0.388, 'distinct': 93, 'cutoff': 18, 'n': 3997, 'max_resolved': 199,
                    'product_min': 0.60, 'product_max': 1.50, 'product_median': 1.047}
print('published ', out['published'])
json.dump(out, open(os.path.join(HERE, 'results', 't2c_slice_published_grid.json'), 'w'), indent=1)
