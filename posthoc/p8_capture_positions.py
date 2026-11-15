"""P8. POST HOC, EXPLORATORY (step 2: measurement or biology; not confirmatory).

From the outputs of P7 (BANC v3 synapse table, size >= 10 voxels, autapses dropped):
  a. capture: of the synapses made by (or onto) a neuron, the share whose partner is a graph neuron, another row of
     the metadata table, no segment at all, or a segment nobody has identified (a fragment);
  b. positions: the share of the synapses of the EB set inside the 22 problem regions of the BANC volume
     (banc_problem_regions.csv of the BANC Dataverse deposit), and the profile of synapse counts along the
     sectioning axis (z, 45 nm sections), where missing or damaged sections show as gaps.
Usage: python3 p8_capture_positions.py <banc_problem_regions.csv>      ASCII, no em dashes.
"""
import csv, json, os, sys
import numpy as np
import load3
from load3 import FAM_NAMES

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')
KIND = ['graph neuron', 'other table row', 'no segment', 'unidentified fragment']
d = load3.load('banc')
N, EB, fam, SETS = d['N'], d['sets']['EB'], d['fam'], d['sets']
z = np.load(os.path.join(R, 'p7_neuron_counts.npz'))
out, inc, tot = z['out'][:, 2, :].astype(float), z['inc'][:, 2, :].astype(float), z['total'][2].astype(float)
res = {'note': 'POST HOC, EXPLORATORY', 'kinds': KIND}

print('a. CAPTURE (size >= 10)')
print(f'all synapses: {int(tot.sum()):,}; pre and post both graph neurons {tot[0, 0] / tot.sum():.3f}; post a fragment {tot[:, 3].sum() / tot.sum():.3f}; '
      f'post no segment {tot[:, 2].sum() / tot.sum():.3f}; pre a fragment {tot[3, :].sum() / tot.sum():.3f}; pre no segment {tot[2, :].sum() / tot.sum():.3f}')
res['all_synapses'] = {'n': int(tot.sum()), 'by_prekind_postkind': tot.tolist()}


def share(a, mask):
    s = a[mask].sum(axis=0)
    return (s / s.sum()).tolist(), float(s.sum())


rows = {}
groups = [('whole graph', np.ones(N, bool)), ('EB set', EB), ('CX set', SETS['CX']), ('MB set', SETS['MB']), ('AL set', SETS['AL'])] + \
         [(f, fam == f) for f in FAM_NAMES]
print(f'{"group":<12} {"n":>6} | outputs: per neuron, share onto graph neurons / other rows / no segment / fragments | inputs: per neuron, share from graph neurons / other / none / fragments')
for nm, mk in groups:
    so, no = share(out, mk); si, ni = share(inc, mk)
    rows[nm] = {'neurons': int(mk.sum()), 'outputs_per_neuron': no / mk.sum(), 'output_shares': so, 'inputs_per_neuron': ni / mk.sum(), 'input_shares': si}
    print(f'{nm:<12} {int(mk.sum()):>6} | {no / mk.sum():8.0f}  ' + ' '.join(f'{x:.3f}' for x in so) + f' | {ni / mk.sum():8.0f}  ' + ' '.join(f'{x:.3f}' for x in si))
res['capture'] = rows
# per neuron capture of outputs: distribution over graph neurons with at least 100 outputs, and where the EB families fall
tot_out = out.sum(axis=1)
big = tot_out >= 100
cap = np.where(tot_out > 0, out[:, 0] / np.maximum(tot_out, 1), np.nan)
pct = np.percentile(cap[big], [5, 25, 50, 75, 95])
print('output capture per neuron (share onto graph neurons), neurons with at least 100 outputs: percentiles 5, 25, 50, 75, 95: ' + ' '.join(f'{x:.3f}' for x in pct))
res['output_capture_percentiles_5_25_50_75_95'] = pct.tolist()
tot_in = inc.sum(axis=1)
capi = np.where(tot_in > 0, inc[:, 0] / np.maximum(tot_in, 1), np.nan)
pcti = np.percentile(capi[tot_in >= 100], [5, 25, 50, 75, 95])
print('input capture per neuron (share from graph neurons), neurons with at least 100 inputs: percentiles 5, 25, 50, 75, 95: ' + ' '.join(f'{x:.3f}' for x in pcti))
res['input_capture_percentiles_5_25_50_75_95'] = pcti.tolist()
for f in FAM_NAMES:
    mk = (fam == f) & big
    if mk.any():
        med = float(np.median(cap[mk]))
        rows[f]['output_capture_median'] = med
        rows[f]['output_capture_percentile_in_graph'] = float((cap[big] < med).mean())
print('median output capture by family and its rank among all graph neurons: ' +
      '; '.join(f'{f} {rows[f]["output_capture_median"]:.3f} (rank {rows[f]["output_capture_percentile_in_graph"]:.2f})' for f in FAM_NAMES if 'output_capture_median' in rows[f]))

print('\nb. POSITIONS')
E = np.load(os.path.join(R, 'p7_eb_synapses.npz'))['rows']
gp, gq, kq, x4, y4, zz = E[:, 0], E[:, 1], E[:, 3], E[:, 5] * 4, E[:, 6] * 4, E[:, 7]
pre_eb = (gp >= 0) & EB[np.clip(gp, 0, N - 1)]
post_eb = (gq >= 0) & EB[np.clip(gq, 0, N - 1)]
inside = pre_eb & post_eb
regs = list(csv.DictReader(open(sys.argv[1])))
hit_any = np.zeros(len(E), bool)
per = {}
for r in regs:
    h = ((x4 >= int(r['min_x'])) & (x4 <= int(r['max_x'])) & (y4 >= int(r['min_y'])) & (y4 <= int(r['max_y'])) &
         (zz >= int(r['min_z'])) & (zz <= int(r['max_z'])))
    hit_any |= h
    per[r['issue']] = per.get(r['issue'], 0) + int((h & inside).sum())
print(f'synapses with an end in EB: {len(E):,}; both ends in EB: {int(inside.sum()):,}')
print(f'both ends in EB, inside a problem region: {int((hit_any & inside).sum()):,} ({(hit_any & inside).sum() / inside.sum():.4f}); by issue: {per}')
print(f'presynaptic neuron in EB (any partner), inside a problem region: {(hit_any & pre_eb).sum() / pre_eb.sum():.4f}')
lo = [np.percentile(v[inside], [0.5, 50, 99.5]).tolist() for v in (x4, y4, zz)]
print(f'extent of the synapses inside EB (0.5, 50, 99.5 percentiles; 4 nm voxels in x and y, sections in z): x {lo[0]}, y {lo[1]}, z {lo[2]}')
res['positions'] = {'eb_involved': int(len(E)), 'inside_EB': int(inside.sum()), 'inside_EB_in_problem_region': int((hit_any & inside).sum()),
                    'by_issue': per, 'pre_in_EB_in_problem_region_share': float((hit_any & pre_eb).sum() / pre_eb.sum()), 'extent': lo}
# profile along z for the synapses inside EB
z0, z1 = np.percentile(zz[inside], [2.5, 97.5]).astype(int)
cnt = np.bincount(zz[inside] - zz[inside].min(), minlength=int(zz[inside].max() - zz[inside].min() + 1))
zs = np.arange(zz[inside].min(), zz[inside].max() + 1)
core = (zs >= z0) & (zs <= z1)
c = cnt[core].astype(float)
k = 25
pad = np.pad(c, k, mode='edge')
med = np.array([np.median(pad[i:i + 2 * k + 1]) for i in range(len(c))])
low = c < 0.2 * med
runs, i = [], 0
while i < len(low):
    if low[i]:
        j = i
        while j < len(low) and low[j]:
            j += 1
        runs.append((int(zs[core][i]), int(j - i)))
        i = j
    else:
        i += 1
print(f'sections {z0} to {z1} (central 95 percent of the synapses inside EB): {len(c)} sections, median {np.median(c):.0f} synapses per section; '
      f'sections below 20 percent of the local median: {int(low.sum())} ({low.mean():.3f}); empty sections: {int((c == 0).sum())}; '
      f'longest low run: {max((r[1] for r in runs), default=0)} sections')
print('  low runs (first section, length): ' + str(runs[:30]))
res['z_profile'] = {'z_from': int(z0), 'z_to': int(z1), 'sections': int(len(c)), 'median_per_section': float(np.median(c)),
                    'low_sections': int(low.sum()), 'low_share': float(low.mean()), 'empty_sections': int((c == 0).sum()), 'low_runs': runs,
                    'expected_loss_if_low_sections_were_filled': float((med[low] - c[low]).sum() / (c.sum() + (med[low] - c[low]).sum()))}
print(f'  synapses that would be added if every low section had the local median: {res["z_profile"]["expected_loss_if_low_sections_were_filled"]:.4f} of the total')
json.dump(res, open(os.path.join(R, 'p8_capture_positions.json'), 'w'), indent=1)
print('saved')
