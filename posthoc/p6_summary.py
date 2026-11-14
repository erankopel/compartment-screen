"""P6. POST HOC, EXPLORATORY. Reads the result files of P1 to P5 and prints the tables quoted in POSTHOC.md.
No computation on connectome data. ASCII, no em dashes."""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')


def J(name):
    return json.load(open(os.path.join(R, name)))


sw = {n: {r['threshold']: r for r in J(f'p1_sweep_{n}.json')['runs']} for n in ('banc', 'male', 'flywire')}
sw['male'].update({r['threshold']: r for r in J('p1_sweep_male_hi.json')['runs']})
NS = {n: J(f'p1_sweep_{n}.json')['N'] for n in sw}
print('TABLE 1. Threshold sweep, EB set (banc / male / flywire). Criteria of the registration for reference: stays/mixes > 1.3, Dobrushin r=8 < 0.6,')
print('         exit mean < 0.37, pairs below 0.5 > 0.70, loss < 0.20. Not scored here.')
print('  th | stays/mixes          | Dobrushin r=8        | exit mean            | pairs below 0.5      | loss per step        | giant component share | EB synapses kept')
for th in (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20):
    cells = []
    for key in ('ratio', 'dob', 'exit', 'below', 'loss', 'scc', 'kept'):
        c = []
        for n in ('banc', 'male', 'flywire'):
            r = sw[n][th]
            v = {'ratio': r['sets']['EB']['ratio'], 'dob': r['dobrushin_EB']['8'], 'exit': r['exit_EB']['mean'],
                 'below': r['exit_EB']['share_below_half'], 'loss': r['sets']['EB']['escape'], 'scc': r['scc'] / NS[n],
                 'kept': r['kept_EB_synapses']}[key]
            c.append(f'{v:6.3f}' if v < 100 else f'{v:6.0f}')
        cells.append(' '.join(c))
    print(f'  {th:>2} | ' + ' | '.join(cells))
print('  first threshold at which stays/mixes <= 1.3: ' + ', '.join(
    f'{n} {next((t for t in sorted(sw[n]) if sw[n][t]["sets"]["EB"]["ratio"] <= 1.3), None)}' for n in sw))
print('  first threshold at which Dobrushin r=8 >= 0.6: ' + ', '.join(
    f'{n} {next((t for t in sorted(sw[n]) if sw[n][t]["dobrushin_EB"]["8"] >= 0.6), None)}' for n in sw))

sc = J('p2_scale.json')
print('\nTABLE 2. Synapse-count scale (banc / male / flywire)')
D = sc['data']
rows = [('synapses per neuron, whole graph', lambda d: d['synapses_per_neuron'], '{:8.0f}'),
        ('synapses per edge, whole graph (mean)', lambda d: d['mean_per_edge'], '{:8.2f}'),
        ('EB neurons', lambda d: d['EB']['neurons'], '{:8.0f}'),
        ('synapses inside EB per EB neuron', lambda d: d['EB']['synapses_inside_per_neuron'], '{:8.0f}'),
        ('synapses per edge inside EB (mean)', lambda d: d['EB']['mean_per_edge'], '{:8.2f}'),
        ('synapses per edge inside EB (median)', lambda d: d['EB']['weight_percentiles_10_25_50_75_90'][2], '{:8.0f}'),
        ('share of edges inside EB with 1 synapse', lambda d: d['EB']['share_edges_w1'], '{:8.3f}'),
        ('share of edges inside EB below 5 synapses', lambda d: d['EB']['share_edges_below5'], '{:8.3f}'),
        ('share of EB synapses in edges below 5', lambda d: d['EB']['share_synapses_in_edges_below5'], '{:8.3f}'),
        ('EB neurons whose strongest input inside EB is below 5', lambda d: d['EB']['neurons_strongest_inside_input_below5'], '{:8.0f}')]
for lab, f, fmt in rows:
    print(f'  {lab:<55}' + ' '.join(fmt.format(f(D[n])) for n in ('banc', 'male', 'flywire')))
print('  summary banc over male:', {k: round(v, 3) for k, v in sc['summary_banc_over_male'].items()})
print('  BANC proofread flag in EB:', D['banc']['banc_flags']['proofread_EB'])
print('  blocks between ring neurons and the rest, synapses per connected pair (banc / male / flywire):')
for kk in ('EL>ER', 'ExR>ER', 'EPG>ER', 'PEN>ER', 'PEG>ER', 'ER>EL', 'ER>ExR', 'ER>EPG', 'ER>PEN', 'ER>PEG', 'ER>ER', 'PEN>EPG', 'EPG>PEN'):
    v = [D[n]['blocks'].get(kk) for n in ('banc', 'male', 'flywire')]
    print(f'    {kk:<8} mean weight ' + ' '.join(f'{x["mean_weight"]:6.2f}' for x in v) + ' | density ' + ' '.join(f'{x["density"]:.3f}' for x in v) +
          ' | synapses per possible pair ' + ' '.join(f'{x["synapses"] / x["possible"]:6.2f}' for x in v))

print('\nTABLE 3. Thinning to the BANC density (mean and range over three seeds). BANC itself: ' +
      '; '.join(f'th {t}: {sw["banc"][t]["sets"]["EB"]["ratio"]:.2f}, {sw["banc"][t]["dobrushin_EB"]["8"]:.2f}' for t in (1, 2, 3, 5)))
for src in ('male', 'flywire'):
    try:
        t4 = J(f'p4_thinning_{src}.json')
    except FileNotFoundError:
        print(f'  {src}: not available'); continue
    print(f'  source {src} (p inside EB {t4["p_inside_EB"]:.3f}, elsewhere {t4["p_elsewhere"]:.3f}); unthinned: ' +
          '; '.join(f'th {t}: {sw[src][t]["sets"]["EB"]["ratio"]:.2f}, {sw[src][t]["dobrushin_EB"]["8"]:.2f}' for t in (1, 2, 3, 5)))
    for var in ('uniform', 'block'):
        for th in (1, 2, 3, 5):
            rr = [r for r in t4['runs'] if r['variant'] == var and r['threshold'] == th]
            if not rr:
                continue
            a = np.array([r['stats']['ratio'] for r in rr]); b = np.array([r['dobrushin']['8'] for r in rr])
            e = np.array([r['exit']['mean'] for r in rr]); h = np.array([r['exit']['share_below_half'] for r in rr])
            print(f'    {var:<8} th {th}: seeds {len(rr)} | stays/mixes {a.mean():.2f} ({a.min():.2f} to {a.max():.2f}) | Dobrushin r=8 {b.mean():.2f} ({b.min():.2f} to {b.max():.2f}) | '
                  f'exit mean {e.mean():.3f} | pairs below 0.5 {h.mean():.3f}')

print('\nTABLE 4. Attribution (one seed): stays/mixes and Dobrushin r=8 at thresholds 1, 2, 5')
for src in ('male', 'flywire'):
    try:
        t5 = J(f'p5_attribution_{src}.json')
    except FileNotFoundError:
        print(f'  {src}: not available'); continue
    print(f'  source {src}; groups (share of EB synapses in the source; synapses kept under the BANC pattern; uniform keeps {t5["p_inside_EB"]:.3f}):')
    for g, v in t5['groups'].items():
        print(f'    {g:<10} {v["share_of_EB_synapses"]:.3f}  {v["banc_over_source"]:.3f}')
    seen = []
    for r in t5['runs']:
        if r['variant'] not in seen:
            seen.append(r['variant'])
    for var in seen:
        rr = {r['threshold']: r for r in t5['runs'] if r['variant'] == var}
        print(f'    {var:<22} ' + ' | '.join(f'th {t}: {rr[t]["ratio"]:5.2f}, {rr[t]["dobrushin"]["8"]:.2f}' for t in (1, 2, 5)))

print('\nTABLE 5. Slow mode of the walk conditioned to stay in EB')
for n in ('banc', 'male', 'flywire'):
    a = J(f'p3_anatomy_{n}.json')
    for r in a['runs']:
        sm, u, s = r['slow_mode'], r['unmixed_r8'], r['subsets']
        tc = sm['type_coordinate_n_mean_mass']
        er3 = sum(v[2] for t, v in tc.items() if t.startswith(('ER3a', 'ER3d', 'ER3m')))
        er3d = sum(v[2] for t, v in tc.items() if t.startswith('ER3d'))
        top = list(tc.items())[:4]
        print(f'  {n:<8} th {r["threshold"]:>2}: second eigenvalue {sm["eigenvalue"][0]:.3f}; stationary mass on ER3a + ER3d + ER3m {er3:.3f} (ER3d {er3d:.3f}); '
              f'types at the positive end: {", ".join(f"{t} ({v[1]:+.2f})" for t, v in top)}')
        print(f'             r = 8: coefficient {u["dobrushin"]:.3f}; pairs above 0.6: {u["share_pairs_above_0.6"]:.3f}; starts farther than 0.3 from the stationary law: '
              f'{u["starts_above_0.3_to_stationary"]}; core subset: stays/mixes {s["core"]["stats"]["ratio"]:.3f}, Dobrushin r=8 {s["core"]["dobrushin"]["8"]:.3f}; '
              f'classes inside EB {r["classes"]["n"]} (largest {r["classes"]["largest"]})')
