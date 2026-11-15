"""P9. POST HOC, EXPLORATORY (step 2: measurement or biology; not confirmatory).

a. Brain-wide, by cell type present in all three data sets: synapses made per neuron (outputs) and received per
   neuron (inputs). For BANC the outputs are counted over ALL partners, identified or not (from P7), so the count
   does not depend on how many postsynaptic twigs were attached to neurons. Ratios BANC : male, BANC : FlyWire and
   FlyWire : male, grouped by the transmitter of the type (male consensus), and the rank of the EB families.
b. The family blocks inside EB (from P2): is the BANC : source ratio of a block the product of a factor of the
   presynaptic family and a factor of the postsynaptic family (weighted two-way fit of the log ratios), and do
   those factors agree with the per-family totals measured in a.
Needs results/p7_neuron_counts.npz and results/p2_scale.json. ASCII, no em dashes.
"""
import json, os
import numpy as np
import pyarrow.feather as feather
import load3
from load3 import FAM_NAMES, family

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')
res = {'note': 'POST HOC, EXPLORATORY'}


def by_type(types, *arrays, min_n=2):
    u, inv = np.unique(types, return_inverse=True)
    n = np.bincount(inv)
    out = {}
    means = [np.bincount(inv, weights=a) / n for a in arrays]
    for i, t in enumerate(u):
        if t != '?' and n[i] >= min_n:
            out[t] = (int(n[i]),) + tuple(float(m[i]) for m in means)
    return out


# BANC: outputs and inputs over all partners (P7) and over graph partners
d = load3.load('banc')
z = np.load(os.path.join(R, 'p7_neuron_counts.npz'))
o_all, i_all = z['out'][:, 2, :].sum(axis=1).astype(float), z['inc'][:, 2, :].sum(axis=1).astype(float)
o_gr, i_gr = z['out'][:, 2, 0].astype(float), z['inc'][:, 2, 0].astype(float)
B = by_type(d['ctype'], o_all, i_all, o_gr, i_gr)
del d
d = load3.load('flywire')
N = d['N']
F = by_type(d['ctype'], np.bincount(d['pre'], weights=d['w'], minlength=N), np.bincount(d['post'], weights=d['w'], minlength=N))
del d
d = load3.load('male')
N = d['N']
mo, mi = np.bincount(d['pre'], weights=d['w'], minlength=N), np.bincount(d['post'], weights=d['w'], minlength=N)
M = by_type(d['ctype'], mo, mi)
ids = np.load(os.path.join(load3.MALE_DIR, 'neuron_edges.npz'))['body_ids']
nt = feather.read_table(os.path.join(load3.MALE_DIR, 'body-nt.feather'), columns=['body', 'consensus_nt']).to_pandas()
nt = nt.drop_duplicates('body').set_index('body').reindex(ids)['consensus_nt'].astype(object).fillna('?').astype(str).to_numpy()
tnt = {}
u, inv = np.unique(d['ctype'], return_inverse=True)
for i, t in enumerate(u):
    v, c = np.unique(nt[inv == i], return_counts=True)
    tnt[t] = v[np.argmax(c)]
del d

common = sorted(t for t in B if t in M and t in F and M[t][1] >= 100 and F[t][1] >= 50)
print(f'cell types (at least 2 neurons each) in BANC {len(B)}, male {len(M)}, FlyWire {len(F)}; in all three with enough synapses: {len(common)}')
rows = {}
for t in common:
    rows[t] = {'nt': tnt.get(t, '?'), 'n': [B[t][0], M[t][0], F[t][0]],
               'out_BM': B[t][1] / M[t][1], 'out_BF': B[t][1] / F[t][1], 'out_FM': F[t][1] / M[t][1],
               'in_BM': B[t][2] / M[t][2], 'in_BF': B[t][2] / F[t][2], 'in_FM': F[t][2] / M[t][2],
               'outgraph_BM': B[t][3] / M[t][1], 'ingraph_BM': B[t][4] / M[t][2]}
keys = ['out_BM', 'out_BF', 'out_FM', 'in_BM', 'in_BF', 'in_FM', 'outgraph_BM', 'ingraph_BM']
A = {k: np.array([rows[t][k] for t in common]) for k in keys}
NT = np.array([rows[t]['nt'] for t in common])
print('\na. ratios per cell type, median (25th to 75th percentile). out = synapses made per neuron, in = received per neuron.')
print('   B = BANC counted over all partners, M = male, F = FlyWire; outgraph, ingraph = BANC counted over identified neurons only (the edge list)')
tab = {}
for lab, mk in [('all types', np.ones(len(common), bool))] + [(x, NT == x) for x in ('acetylcholine', 'gaba', 'glutamate', 'dopamine', 'serotonin', 'octopamine', 'histamine')]:
    if mk.sum() < 5:
        continue
    tab[lab] = {'types': int(mk.sum())}
    line = f'  {lab:<14} {int(mk.sum()):>5} types | '
    for k in keys:
        q = np.percentile(A[k][mk], [25, 50, 75])
        tab[lab][k] = q.tolist()
        line += f'{k} {q[1]:.2f} ({q[0]:.2f} to {q[2]:.2f}) | '
    print(line)
res['by_transmitter'] = tab


def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


res['spearman'] = {'out_BM_vs_out_FM': spearman(A['out_BM'], A['out_FM']), 'in_BM_vs_in_FM': spearman(A['in_BM'], A['in_FM']),
                   'out_BM_vs_in_BM': spearman(A['out_BM'], A['in_BM'])}
print('   rank correlation across types: out B:M with out F:M %.2f; in B:M with in F:M %.2f; out B:M with in B:M %.2f' %
      tuple(res['spearman'].values()))
print('   spread across types (ratio of the 90th to the 10th percentile): ' +
      ', '.join(f'{k} {np.percentile(A[k], 90) / np.percentile(A[k], 10):.1f}' for k in ('out_BM', 'out_BF', 'out_FM', 'in_BM', 'in_BF', 'in_FM')))
res['spread_p90_over_p10'] = {k: float(np.percentile(A[k], 90) / np.percentile(A[k], 10)) for k in keys}

# EB families: their types, with the rank of each ratio among all matched types
print('\n   EB families (types in all three data sets): ratio and rank among all matched types')
fam_rows = {}
for f in FAM_NAMES:
    ts = [t for t in common if family(t) >= 0 and FAM_NAMES[family(t)] == f]
    if not ts:
        continue
    w = np.array([rows[t]['n'][0] for t in ts], float)
    fr = {'types': ts, 'nt': sorted({rows[t]['nt'] for t in ts})}
    for k in ('out_BM', 'out_BF', 'out_FM', 'in_BM', 'in_BF', 'in_FM'):
        v = float(np.average([rows[t][k] for t in ts], weights=w))
        fr[k] = v
        fr[k + '_rank'] = float((A[k] < v).mean())
    fam_rows[f] = fr
    print(f'   {f:<8} {len(ts):>2} types {"/".join(fr["nt"]):<24} out B:M {fr["out_BM"]:.2f} (rank {fr["out_BM_rank"]:.2f}) B:F {fr["out_BF"]:.2f} ({fr["out_BF_rank"]:.2f}) F:M {fr["out_FM"]:.2f} ({fr["out_FM_rank"]:.2f}) | '
          f'in B:M {fr["in_BM"]:.2f} ({fr["in_BM_rank"]:.2f}) B:F {fr["in_BF"]:.2f} ({fr["in_BF_rank"]:.2f}) F:M {fr["in_FM"]:.2f} ({fr["in_FM_rank"]:.2f})')
res['eb_families'] = fam_rows
res['types'] = rows

# b. two-way fit of the block ratios
sc = json.load(open(os.path.join(R, 'p2_scale.json')))['data']
tot = json.load(open(os.path.join(R, 'p8_family_totals.json')))
print('\nb. family blocks inside EB: log ratio = constant + presynaptic family + postsynaptic family (weights: square root of the source synapses)')
res['two_way'] = {}
for src in ('male', 'flywire'):
    S, Bk = sc[src]['blocks'], sc['banc']['blocks']
    blocks = [(k, np.log((Bk[k]['synapses'] / Bk[k]['possible']) / (S[k]['synapses'] / S[k]['possible'])), np.sqrt(S[k]['synapses']))
              for k in S if S[k]['synapses'] >= 2000 and k in Bk]
    pre = sorted({k.split('>')[0] for k, _, _ in blocks}); post = sorted({k.split('>')[1] for k, _, _ in blocks})
    X = np.zeros((len(blocks), 1 + len(pre) + len(post))); y = np.zeros(len(blocks)); w = np.zeros(len(blocks))
    for r, (k, v, ww) in enumerate(blocks):
        a, b = k.split('>')
        X[r, 0] = 1; X[r, 1 + pre.index(a)] = 1; X[r, 1 + len(pre) + post.index(b)] = 1; y[r] = v; w[r] = ww
    sw = np.sqrt(w)
    beta = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    fit = X @ beta
    ybar = np.average(y, weights=w)
    r2 = 1 - (w * (y - fit) ** 2).sum() / (w * (y - ybar) ** 2).sum()
    # one-way fits for comparison
    def oneway(cols):
        Xc = X[:, [0] + cols]
        b1 = np.linalg.lstsq(Xc * sw[:, None], y * sw, rcond=None)[0]
        return 1 - (w * (y - Xc @ b1) ** 2).sum() / (w * (y - ybar) ** 2).sum()
    r2_pre, r2_post = oneway(list(range(1, 1 + len(pre)))), oneway(list(range(1 + len(pre), X.shape[1])))
    resid = y - fit
    worst = sorted(zip(np.abs(resid), [k for k, _, _ in blocks], resid), reverse=True)[:6]
    # marginal factors measured independently of the blocks: per-neuron totals of the family (P8)
    key = 'male' if src == 'male' else 'flywire'
    marg_out = {f: tot[f]['banc_all'][0] / tot[f][key][1] for f in pre}       # BANC outputs to all partners : source outputs
    marg_in = {f: tot[f]['banc_all'][1] / tot[f][key][2] for f in post}        # BANC inputs from all partners : source inputs
    # typical block ratio of a family as presynaptic and as postsynaptic partner (weighted mean of the log ratios)
    row_eff = {a: float(np.exp(np.average([v for k, v, _ in blocks if k.split('>')[0] == a], weights=[ww for k, _, ww in blocks if k.split('>')[0] == a]))) for a in pre}
    col_eff = {b: float(np.exp(np.average([v for k, v, _ in blocks if k.split('>')[1] == b], weights=[ww for k, _, ww in blocks if k.split('>')[1] == b]))) for b in post}
    c_out = spearman(np.array([row_eff[a] for a in pre]), np.array([marg_out[a] for a in pre]))
    c_in = spearman(np.array([col_eff[b] for b in post]), np.array([marg_in[b] for b in post]))
    # prediction of each block from the two marginals alone (no fitting to the blocks except one constant)
    pred = np.array([np.log(marg_out[k.split('>')[0]]) + np.log(marg_in[k.split('>')[1]]) for k, _, _ in blocks])
    const = np.average(y - pred, weights=w)
    r2_marg = 1 - (w * (y - pred - const) ** 2).sum() / (w * (y - ybar) ** 2).sum()
    res['two_way'][src] = {'blocks': len(blocks), 'r2_two_way': float(r2), 'r2_pre_only': float(r2_pre), 'r2_post_only': float(r2_post),
                           'r2_from_marginals': float(r2_marg), 'row_effect': row_eff, 'col_effect': col_eff, 'marginal_out': marg_out,
                           'marginal_in': marg_in, 'spearman_row_vs_marginal_out': c_out, 'spearman_col_vs_marginal_in': c_in,
                           'largest_residuals': [[k, float(r)] for _, k, r in worst]}
    print(f'  against {src}: {len(blocks)} blocks; variance of the log ratios explained: presynaptic family alone {r2_pre:.2f}, postsynaptic family alone {r2_post:.2f}, '
          f'both {r2:.2f}; predicted from the per-family totals alone (one constant) {r2_marg:.2f}')
    print('     as presynaptic partner (typical block ratio; BANC outputs to all partners : source outputs): ' + '; '.join(f'{a} {row_eff[a]:.2f}, {marg_out[a]:.2f}' for a in pre))
    print('     as postsynaptic partner (typical block ratio; BANC inputs from all partners : source inputs): ' + '; '.join(f'{b} {col_eff[b]:.2f}, {marg_in[b]:.2f}' for b in post))
    print(f'     rank correlation of the two columns: presynaptic {c_out:.2f}, postsynaptic {c_in:.2f}; largest residuals of the two-way fit (log): ' +
          ', '.join(f'{k} {r:+.2f}' for _, k, r in worst))
json.dump(res, open(os.path.join(R, 'p9_types.json'), 'w'), indent=1)
print('saved')
