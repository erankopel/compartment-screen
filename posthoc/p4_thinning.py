"""P4. POST HOC, EXPLORATORY (after the confirmatory run of 14 Nov 2026; not confirmatory).

Thinning test. Question: if the male (or FlyWire) connectome had the synapse density of BANC, would its EB set
behave like the EB set of BANC? Every synapse of the source data set is kept independently with a probability p
(binomial thinning of the edge counts; edges that lose all synapses disappear):
  U  uniform: p_in for edges inside EB, p_out for all other edges, where p_in and p_out are the ratios
     BANC : source of synapses per neuron inside EB and in the whole graph;
  B  block pattern: inside EB, p of an edge is the ratio BANC : source of synapses per possible pair of its
     family block (from p2_scale.json; capped at 1; 0 where BANC has no edge in the block); p_out elsewhere.
Then the frozen statistics of the EB set at thresholds 1, 2, 3 and 5. Three seeds per variant.
The stationary distribution is iterated 300 times only; it enters the loss per step and nothing else.
Usage: python3 p4_thinning.py <male | flywire>      ASCII, no em dashes.
"""
import json, os, sys, time
import numpy as np
import load3
from load3 import sl, FAM_NAMES

src = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', f'p4_thinning_{src}.json')
SC = json.load(open(os.path.join(HERE, 'results', 'p2_scale.json')))['data']
B, S = SC['banc'], SC[src]
p_in = B['EB']['synapses_inside_per_neuron'] / S['EB']['synapses_inside_per_neuron']
p_out = B['synapses_per_neuron'] / S['synapses_per_neuron']
THS = [1, 2, 3, 5]
SEEDS = [20261115, 20261116, 20261117]


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
p_block = np.full(len(w), p_out)
p_block[inEB] = pb[fi[pre[inEB]], fi[post[inEB]]]
p_unif = np.where(inEB, p_in, p_out)
log(f'{src}: p inside EB {p_in:.4f}, elsewhere {p_out:.4f}; block probabilities inside EB: synapse-weighted mean '
    f'{float((p_block[inEB] * w[inEB]).sum() / w[inEB].sum()):.4f}, range {pb[pb > 0].min():.3f} to {pb.max():.3f}')
res = {'note': 'POST HOC, EXPLORATORY', 'source': src, 'p_inside_EB': p_in, 'p_elsewhere': p_out, 'runs': []}


def stats(pp, qq, ww, label, seed):
    for th in THS:
        t0 = time.time()
        m = ww >= th
        ch = sl.build_chain(pp[m], qq[m], ww[m], N, maxit=300)
        k = ch['keep']
        st = sl.set_stats(ch, EB[k])
        db = sl.dobrushin_conditioned(ch, EB[k])
        ex = sl.exit_channel(ch, EB[k])
        r = {'variant': label, 'seed': seed, 'threshold': th, 'scc': ch['n'], 'EB_in_scc': int(EB[k].sum()), 'delta': ch['delta'],
             'stats': st, 'dobrushin': db, 'exit': ex}
        res['runs'].append(r)
        log(f'  {label:<8} seed {seed} th {th}: scc {ch["n"]:,} EB {r["EB_in_scc"]} | loss {st["escape"]:.4f} stays {st["stays"]:.2f} mixes {st["mixes"]:.2f} '
            f'ratio {st["ratio"]:.3f} | Dob8 {db["8"]:.3f} (drop {db["dropped"]}) | exit mean {ex["mean"]:.3f} below {ex["share_below_half"]:.3f} | {time.time() - t0:.0f} s')
        json.dump(res, open(OUT, 'w'), indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x))


for label, pv in (('uniform', p_unif), ('block', p_block)):
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        wt = rng.binomial(w, pv).astype(np.int32)
        keep = wt > 0
        ie = inEB & keep
        log(f'{label}, seed {seed}: edges kept {int(keep.sum()):,} of {len(w):,}; synapses {int(wt.sum()):,}; inside EB: edges {int(ie.sum()):,}, '
            f'synapses {int(wt[ie].sum()):,} (per neuron {wt[ie].sum() / EB.sum():.0f}; BANC {B["EB"]["synapses_inside_per_neuron"]:.0f}); '
            f'mean per edge {wt[keep].mean():.2f} (BANC {B["mean_per_edge"]:.2f}), inside EB {wt[ie].mean():.2f} (BANC {B["EB"]["mean_per_edge"]:.2f})')
        stats(pre[keep], post[keep], wt[keep], label, seed)
log('saved ' + OUT)
