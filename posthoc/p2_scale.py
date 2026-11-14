"""P2. POST HOC, EXPLORATORY (after the confirmatory run of 14 Nov 2026; not confirmatory).

Synapse-count scale and data quality of the EB set in the three data sets: synapses per neuron and per edge,
the weight distribution of edges inside EB, family-by-family blocks (connection density, mean weight), and the
proofreading flags of BANC. ASCII, no em dashes.
"""
import json, os
import numpy as np
import load3
from load3 import FAM_NAMES

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', 'p2_scale.json')
res = {'note': 'POST HOC, EXPLORATORY', 'data': {}}


def q(x, ps=(10, 25, 50, 75, 90)):
    return [float(v) for v in np.percentile(x, ps)] if len(x) else []


for name in ('banc', 'male', 'flywire'):
    d = load3.load(name)
    pre, post, w, N, fam = d['pre'], d['post'], d['w'], d['N'], d['fam']
    EB = d['sets']['EB']
    nEB = int(EB.sum())
    inEB = EB[pre] & EB[post]
    we = w[inEB]
    out_all = np.bincount(pre, weights=w, minlength=N)
    in_all = np.bincount(post, weights=w, minlength=N)
    out_eb = np.bincount(pre[inEB], weights=we, minlength=N)
    in_eb = np.bincount(post[inEB], weights=we, minlength=N)
    r = {'N': N, 'edges': int(len(w)), 'synapses': int(w.sum()), 'synapses_per_neuron': float(w.sum() / N),
         'mean_per_edge': float(w.mean()), 'median_per_edge': float(np.median(w)),
         'share_edges_w1': float((w == 1).mean()), 'share_edges_below5': float((w < 5).mean()),
         'share_synapses_in_edges_below5': float(w[w < 5].sum() / w.sum()),
         'EB': {'neurons': nEB, 'edges': int(inEB.sum()), 'synapses': int(we.sum()),
                'synapses_inside_per_neuron': float(we.sum() / nEB), 'mean_per_edge': float(we.mean()),
                'weight_percentiles_10_25_50_75_90': q(we),
                'share_edges_w1': float((we == 1).mean()), 'share_edges_below5': float((we < 5).mean()),
                'share_synapses_in_edges_below5': float(we[we < 5].sum() / we.sum()),
                'out_all_per_neuron_median': float(np.median(out_all[EB])),
                'in_all_per_neuron_median': float(np.median(in_all[EB])),
                'out_inside_per_neuron_median': float(np.median(out_eb[EB])),
                'in_inside_per_neuron_median': float(np.median(in_eb[EB])),
                'side': dict(zip(*[x.tolist() for x in np.unique(d['side'][EB], return_counts=True)]))},
         'families': {}, 'blocks': {}}
    fi = np.array([FAM_NAMES.index(f) if f else -1 for f in fam])
    for i, f in enumerate(FAM_NAMES):
        mk = fi == i
        r['families'][f] = {'neurons': int(mk.sum()),
                            'out_all_median': float(np.median(out_all[mk])) if mk.any() else None,
                            'in_all_median': float(np.median(in_all[mk])) if mk.any() else None,
                            'out_inside_median': float(np.median(out_eb[mk])) if mk.any() else None,
                            'in_inside_median': float(np.median(in_eb[mk])) if mk.any() else None,
                            'in_inside_min': float(in_eb[mk].min()) if mk.any() else None,
                            'neurons_with_inside_input_below5_max_edge': None}
    # strongest incoming and outgoing edge inside EB for every EB neuron (what a threshold removes)
    mx_in = np.zeros(N); np.maximum.at(mx_in, post[inEB], we)
    mx_out = np.zeros(N); np.maximum.at(mx_out, pre[inEB], we)
    for i, f in enumerate(FAM_NAMES):
        mk = fi == i
        r['families'][f]['neurons_with_inside_input_below5_max_edge'] = int((mx_in[mk] < 5).sum())
        r['families'][f]['neurons_with_inside_output_below5_max_edge'] = int((mx_out[mk] < 5).sum())
    r['EB']['neurons_strongest_inside_input_below5'] = int((mx_in[EB] < 5).sum())
    r['EB']['neurons_strongest_inside_output_below5'] = int((mx_out[EB] < 5).sum())
    a, b = fi[pre[inEB]], fi[post[inEB]]
    nf = np.bincount(fi[fi >= 0], minlength=len(FAM_NAMES))
    key = a * len(FAM_NAMES) + b
    cnt = np.bincount(key, minlength=len(FAM_NAMES) ** 2)
    syn = np.bincount(key, weights=we, minlength=len(FAM_NAMES) ** 2)
    syn5 = np.bincount(key[we >= 5], weights=we[we >= 5], minlength=len(FAM_NAMES) ** 2)
    for i, fa in enumerate(FAM_NAMES):
        for j, fb in enumerate(FAM_NAMES):
            kk = i * len(FAM_NAMES) + j
            poss = int(nf[i]) * int(nf[j]) - (int(nf[i]) if i == j else 0)
            if cnt[kk]:
                r['blocks'][f'{fa}>{fb}'] = {'possible': poss, 'edges': int(cnt[kk]), 'synapses': int(syn[kk]),
                                             'density': float(cnt[kk] / poss), 'mean_weight': float(syn[kk] / cnt[kk]),
                                             'synapses_in_edges_ge5': int(syn5[kk])}
    if name == 'banc':
        ex = d['extra']
        pr = np.array([str(x) for x in ex['proofread']])
        st = np.array([str(x) for x in ex['status']])
        r['banc_flags'] = {'proofread_EB': dict(zip(*[x.tolist() for x in np.unique(pr[EB], return_counts=True)])),
                           'proofread_all': dict(zip(*[x.tolist() for x in np.unique(pr, return_counts=True)])),
                           'status_EB': dict(zip(*[x.tolist() for x in np.unique(st[EB], return_counts=True)]))}
    res['data'][name] = r
    e = r['EB']
    print(f'{name}: syn/neuron {r["synapses_per_neuron"]:.0f}; per edge mean {r["mean_per_edge"]:.2f} median {r["median_per_edge"]:.0f}; '
          f'edges w=1 {r["share_edges_w1"]:.3f}, <5 {r["share_edges_below5"]:.3f} (synapses in them {r["share_synapses_in_edges_below5"]:.3f})')
    print(f'   EB: {e["neurons"]} neurons, {e["edges"]:,} edges, {e["synapses"]:,} synapses; inside per neuron {e["synapses_inside_per_neuron"]:.0f}; '
          f'per edge mean {e["mean_per_edge"]:.2f}, pct {e["weight_percentiles_10_25_50_75_90"]}; w=1 {e["share_edges_w1"]:.3f}, <5 {e["share_edges_below5"]:.3f} '
          f'(synapses in them {e["share_synapses_in_edges_below5"]:.3f}); strongest inside input <5: {e["neurons_strongest_inside_input_below5"]}, '
          f'output <5: {e["neurons_strongest_inside_output_below5"]}', flush=True)
    del d

json.dump(res, open(OUT, 'w'), indent=1)
B, M, F = (res['data'][k] for k in ('banc', 'male', 'flywire'))
print('\nfamily: neurons (banc male flywire) | median inside input per neuron | neurons whose strongest inside input < 5')
for f in FAM_NAMES:
    x = [res['data'][k]['families'][f] for k in ('banc', 'male', 'flywire')]
    print(f'  {f:<8} {x[0]["neurons"]:>4} {x[1]["neurons"]:>4} {x[2]["neurons"]:>4} | {x[0]["in_inside_median"]:>7.0f} {x[1]["in_inside_median"]:>7.0f} '
          f'{x[2]["in_inside_median"]:>7.0f} | {x[0]["neurons_with_inside_input_below5_max_edge"]:>3} {x[1]["neurons_with_inside_input_below5_max_edge"]:>3} '
          f'{x[2]["neurons_with_inside_input_below5_max_edge"]:>3}')
print('\nblocks with at least 2000 synapses in the male: synapses per possible pair, density, mean weight (banc / male / flywire) and ratio banc:male')
rows = []
for kk, m in M['blocks'].items():
    if m['synapses'] >= 2000:
        b, f = B['blocks'].get(kk), F['blocks'].get(kk)
        rows.append((m['synapses'], kk, b, m, f))
ratios_w, ratios_d, ratios_s, wts = [], [], [], []
for _, kk, b, m, f in sorted(rows, key=lambda t: -t[0]):
    def g(x, key):
        return x[key] if x else float('nan')
    sp = [g(x, 'synapses') / x['possible'] if x else float('nan') for x in (b, m, f)]
    print(f'  {kk:<14} syn/pair {sp[0]:6.2f} {sp[1]:6.2f} {sp[2]:6.2f} | density {g(b, "density"):.3f} {g(m, "density"):.3f} {g(f, "density"):.3f} | '
          f'mean w {g(b, "mean_weight"):6.2f} {g(m, "mean_weight"):6.2f} {g(f, "mean_weight"):6.2f} | banc:male syn/pair {sp[0] / sp[1]:.2f}, mean w {g(b, "mean_weight") / m["mean_weight"]:.2f}')
    if b:
        ratios_s.append(sp[0] / sp[1]); ratios_w.append(b['mean_weight'] / m['mean_weight']); ratios_d.append(b['density'] / m['density']); wts.append(m['synapses'])
ratios_s, ratios_w, ratios_d, wts = map(np.array, (ratios_s, ratios_w, ratios_d, wts))
summ = {'blocks_compared': int(len(wts)), 'median_ratio_syn_per_pair': float(np.median(ratios_s)),
        'median_ratio_mean_weight': float(np.median(ratios_w)), 'median_ratio_density': float(np.median(ratios_d)),
        'weighted_mean_ratio_syn_per_pair': float((ratios_s * wts).sum() / wts.sum()),
        'weighted_mean_ratio_mean_weight': float((ratios_w * wts).sum() / wts.sum()),
        'ratio_inside_synapses_per_neuron_banc_male': B['EB']['synapses_inside_per_neuron'] / M['EB']['synapses_inside_per_neuron'],
        'ratio_inside_synapses_per_neuron_banc_flywire': B['EB']['synapses_inside_per_neuron'] / F['EB']['synapses_inside_per_neuron'],
        'ratio_synapses_per_neuron_whole_banc_male': B['synapses_per_neuron'] / M['synapses_per_neuron'],
        'ratio_synapses_per_neuron_whole_banc_flywire': B['synapses_per_neuron'] / F['synapses_per_neuron']}
res['summary_banc_over_male'] = summ
json.dump(res, open(OUT, 'w'), indent=1)
print('\nsummary', json.dumps(summ, indent=1))
if 'banc_flags' in B:
    print('BANC flags', B['banc_flags'])
