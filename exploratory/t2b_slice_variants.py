"""T2b. The 81 x 81 slice of arXiv:2609.27908 Sec. 6 under different readings of the grid and of the
numerical tolerance, to see which reading reproduces the published counts. ASCII, no em dashes."""
import json, math, os, sys
import numpy as np
sys.argv = [sys.argv[0]]
import importlib.util
spec = importlib.util.spec_from_file_location('t2', os.path.join(os.path.dirname(os.path.abspath(__file__)), 't2_core.py'))
t2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t2)
pi = math.pi
k0, b0 = 0.43230980 * pi, 0.23903823 * pi
variants = {'inclusive [0, pi/2]': np.linspace(0, pi / 2, 81),
            'endpoint excluded': np.linspace(0, pi / 2, 81, endpoint=False),
            'cell centred': (np.arange(81) + 0.5) * (pi / 2) / 81,
            'open interval (83 points, ends dropped)': np.linspace(0, pi / 2, 83)[1:-1]}
out = {}
for name, grid in variants.items():
    for tol in (1e-12, 1e-9, 0.0):
        Ns, prods = [], []
        for th in grid:
            for ph in grid:
                A, c, env, off = t2.channel_data(t2.unitary(th, k0, ph, b0))
                n1 = t2.index(A, c, cap=200, tol=tol)
                Ns.append(n1)
                if 1 < n1 <= 200:
                    prods.append(n1 * (1 - float(np.abs(np.linalg.eigvals(A)).max())))
        Ns, prods = np.array(Ns), np.array(prods)
        r = {'share_index_1': float((Ns == 1).mean()), 'distinct': int(len(np.unique(Ns[Ns <= 200]))),
             'cutoff': int((Ns > 200).sum()), 'n': int(len(prods)), 'min': float(prods.min()),
             'max': float(prods.max()), 'median': float(np.median(prods))}
        out[f'{name}; tol {tol:g}'] = r
        print(f'{name:<40s} tol {tol:<6g} N=1 {r["share_index_1"]:.3f} distinct {r["distinct"]:>3d} cutoff '
              f'{r["cutoff"]:>3d} n {r["n"]} product [{r["min"]:.2f}, {r["max"]:.2f}] median {r["median"]:.3f}', flush=True)
print('paper: N=1 0.388, distinct 93, cutoff 18, n 3997, product [0.60, 1.50], median 1.047')
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 't2b_slice_variants.json'), 'w'), indent=1)
