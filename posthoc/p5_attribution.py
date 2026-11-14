"""P5. POST HOC, EXPLORATORY (after the confirmatory run of 14 Nov 2026; not confirmatory).

Attribution of the thinning test. The edges inside EB are split into five groups by the families of their ends:
  into_ER    columnar and extrinsic neurons (EPG, PEG, PEN, EL, ExR, LPsP, GLNO) onto ring neurons (ER)
  ER_ER      ring neurons onto ring neurons
  out_of_ER  ring neurons onto the columnar and extrinsic neurons
  columnar   among the columnar and extrinsic neurons
  tail       any edge with an end in the visual input tail (TuBu, MeTu2, TuTuB, AOTU046)
'only g': the BANC block pattern of P4 is applied in group g, uniform thinning (p_in) in the other groups.
'all but g': the block pattern everywhere inside EB except in group g, which is thinned uniformly.
Outside EB every edge is thinned with p_out. One seed; thresholds 1, 2 and 5. The stationary distribution is
iterated 60 times only, so the loss per step is not reported from this script.
Usage: python3 p5_attribution.py <male | flywire>      ASCII, no em dashes.
"""
import json, os, sys, time
import numpy as np
import load3
from load3 import sl, FAM_NAMES

src = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', f'p5_attribution_{src}.json')
SC = json.load(open(os.path.join(HERE, 'results', 'p2_scale.json')))['data']
B, S = SC['banc'], SC[src]
p_in = B['EB']['synapses_inside_per_neuron'] / S['EB']['synapses_inside_per_neuron']
p_out = B['synapses_per_neuron'] / S['synapses_per_neuron']
THS = [1, 2, 5]
SEED = 20261115
COL = ['EPG', 'PEG', 'PEN', 'EL', 'ExR', 'LPsP', 'GLNO']
TAIL = ['TuBu', 'MeTu2', 'TuTuB', 'AOTU046']


def log(s):
    print(s, flush=True)


d = load3.load(src)
pre, post, w, N, fam = d['pre'], d['post'], d['w'], d['N'], d['fam']
EB = d['sets']['EB']
inEB = EB[pre] & EB[post]
fi = np.array([FAM_NAMES.index(f) if f else -1 for f in fam])
nF = len(FAM_NAMES)
pb = np.zeros((nF, nF))
for i, fa in enumerate(FAM_NAMES):
    for j, fb in enumerate(FAM_NAMES):
        s_, b_ = S['blocks'].get(f'{fa}>{fb}'), B['blocks'].get(f'{fa}>{fb}')
        if s_ and b_:
            pb[i, j] = min(1.0, (b_['synapses'] / b_['possible']) / (s_['synapses'] / s_['possible']))
ia, ib = fi[pre].astype(np.int8), fi[post].astype(np.int8)          # integer family codes (string arrays per edge cost memory)
COLI, TAILI, ERI = [FAM_NAMES.index(x) for x in COL], [FAM_NAMES.index(x) for x in TAIL], FAM_NAMES.index('ER')
col_a, col_b, er_a, er_b = np.isin(ia, COLI), np.isin(ib, COLI), ia == ERI, ib == ERI
tail = inEB & (np.isin(ia, TAILI) | np.isin(ib, TAILI))
del ia, ib
GROUPS = {'into_ER': inEB & col_a & er_b, 'ER_ER': inEB & er_a & er_b, 'out_of_ER': inEB & er_a & col_b,
          'columnar': inEB & col_a & col_b, 'tail': tail}
cover = np.zeros(len(w), int)
for g in GROUPS.values():
    cover += g
assert (cover[inEB] == 1).all(), 'groups must partition the edges inside EB'
p_block_in = np.zeros(len(w)); p_block_in[inEB] = pb[fi[pre[inEB]], fi[post[inEB]]]
res = {'note': 'POST HOC, EXPLORATORY', 'source': src, 'p_inside_EB': p_in, 'p_elsewhere': p_out, 'groups': {}, 'runs': []}
for g, mk in GROUPS.items():
    res['groups'][g] = {'edges': int(mk.sum()), 'synapses': int(w[mk].sum()), 'share_of_EB_synapses': float(w[mk].sum() / w[inEB].sum()),
                        'banc_over_source': float((p_block_in[mk] * w[mk]).sum() / w[mk].sum())}
    log(f'group {g:<10}: edges {int(mk.sum()):>6,}, synapses {int(w[mk].sum()):>8,} ({res["groups"][g]["share_of_EB_synapses"]:.3f} of EB), '
        f'synapses kept under the BANC pattern {res["groups"][g]["banc_over_source"]:.3f} (uniform {p_in:.3f})')


def run(label, pv):
    rng = np.random.default_rng(SEED)
    wt = rng.binomial(w, pv).astype(np.int32)
    keep = wt > 0
    for th in THS:
        t0 = time.time()
        m = keep & (wt >= th)
        ch = sl.build_chain(pre[m], post[m], wt[m], N, maxit=60)
        k = ch['keep']
        st = sl.set_stats(ch, EB[k]); db = sl.dobrushin_conditioned(ch, EB[k])
        res['runs'].append({'variant': label, 'threshold': th, 'EB_in_scc': int(EB[k].sum()),
                            'stays': st['stays'], 'mixes': st['mixes'], 'ratio': st['ratio'], 'dobrushin': db})
        log(f'  {label:<22} th {th}: EB {int(EB[k].sum())} | stays {st["stays"]:.2f} mixes {st["mixes"]:.2f} ratio {st["ratio"]:.3f} | '
            f'Dob8 {db["8"]:.3f} | {time.time() - t0:.0f} s')
        json.dump(res, open(OUT, 'w'), indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x))


base = np.where(inEB, p_in, p_out)
full = np.where(inEB, p_block_in, p_out)
run('uniform', base)
run('block pattern', full)
for g, mk in GROUPS.items():
    pv = base.copy(); pv[mk] = p_block_in[mk]
    run(f'only {g}', pv)
for g, mk in GROUPS.items():
    pv = full.copy(); pv[mk] = p_in
    run(f'all but {g}', pv)
log('saved ' + OUT)
