"""T6. EXPLORATORY (male CNS v1.0). The screen of screen_lib.py for sets defined by cell-type names,
at >= 1 and >= 5 synapses, and under the three null models of the maleCNS-depth repository (E3: targets
permuted within one group, within superclass pairs, within cell-type pairs). Also writes the crosswalk
of type names used to carry the sets to other data sets. ASCII, no em dashes."""
import json, os, sys, time
import numpy as np
REPO = os.environ['MALECNS_REPO']
sys.path.insert(0, os.path.join(REPO, 'code', 'confirmatory'))
from common import load_edges, annotations, consensus_nt, log, _json_default  # noqa: E402
import screen_lib as sl  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', 't6_male_typed_nulls.json')
pre, post, w, body_ids = load_edges()
N = len(body_ids)
ann = annotations(body_ids, cols=('type', 'class', 'superclass', 'flywireType'))
typ, cls, fw = ann['type'].astype(str), ann['class'].astype(str), ann['flywireType'].astype(str)
nt = consensus_nt(body_ids).astype(str)
e1 = json.load(open(os.path.join(REPO, 'exploratory', 'results', 'e1_slow_mode.json')))
in_cut = np.isin(body_ids, np.array(e1['w1']['cut_body_ids'], dtype=np.int64))

# sets by male cell type: types with at least half of their neurons in the E1 cut
ut, inv = np.unique(typ, return_inverse=True)
tot, inc = np.bincount(inv), np.bincount(inv, weights=in_cut)
eb_types = [t for t, a, b in zip(ut, inc, tot) if t != '?' and a / b >= 0.5]
EB = np.isin(typ, eb_types)
AMINES = ['dopamine', 'serotonin', 'octopamine']
amine = np.isin(nt, AMINES)
SETS = {'EB': EB, 'CX': (cls == 'CX') | EB, 'MB': np.isin(cls, ['Kenyon_Cell', 'MBON', 'DAN']), 'AL': np.isin(cls, ['ALPN', 'ALLN'])}
log(f'EB set by type: {len(eb_types)} male types, {int(EB.sum())} neurons, {int((EB & in_cut).sum())} of them in the E1 cut (727)')

# crosswalk to FlyWire-style names (the flywireType column of the male annotations; composite names are also split)
def names(mask):
    s = set()
    for x in np.unique(fw[mask]):
        if x != '?':
            s.add(x); s.update(p.strip() for p in x.split(','))
    return sorted(s)

uf, invf = np.unique(fw, return_inverse=True)
totf, amf = np.bincount(invf), np.bincount(invf, weights=amine)
amine_fw = sorted({p.strip() for t, a, b in zip(uf, amf, totf) if t != '?' and a / b >= 0.5 for p in [t] + t.split(',')})
cross = {'EB_family_regex': ['^ER[1-6]', '^EPG$', '^PEG$', '^PEN_', '^EL$', '^ExR[12346](?![0-9])', '^TuBu', '^MeTu2', '^LPsP$', '^GLNO$', '^AOTU046$', '^TuTuB'],
         'EB_rule': 'a neuron is in EB if any comma-separated part of its cell type name matches one of EB_family_regex',
         'EB_types_male': eb_types, 'EB_types_flywire_names': names(EB), 'amine_types_flywire_names': amine_fw,
         'amine_definition': 'male consensus transmitter in dopamine, serotonin, octopamine for at least half of the neurons of the type',
         'class_map': {'CX': ['central_complex_intrinsic_neuron', 'central_complex_input_neuron', 'central_complex_output_neuron'],
                       'MB': ['kenyon_cell', 'mushroom_body_output_neuron', 'mushroom_body_dopaminergic_neuron'],
                       'AL': ['antennal_lobe_projection_neuron', 'antennal_lobe_local_neuron']}}
json.dump(cross, open(os.path.join(HERE, 'crosswalk_types.json'), 'w'), indent=1)
log(f'crosswalk: {len(cross["EB_types_flywire_names"])} EB names, {len(amine_fw)} amine type names')
amine_by_type = np.isin(fw, amine_fw) | False
log(f'amine neurons: consensus {int(amine.sum())}; by type-name list {int(amine_by_type.sum())}; agree on {int((amine & amine_by_type).sum())}')


def run(p, q, ww, label, full):
    t0 = time.time()
    ch = sl.build_chain(p, q, ww, N)
    k = ch['keep']
    log(f'--- {label}: SCC {ch["n"]:,}, stationary in {ch["iters"]} iterations')
    r = sl.screen(ch, {a: b[k] for a, b in SETS.items()}, amine[k], full=full, log=log)
    r['amine_by_type_list'] = {nm: sl.label_shares(ch, SETS[nm][k], amine_by_type[k]) for nm in ('EB', 'MB')}
    r['label'] = label; r['wall_s'] = round(time.time() - t0)
    return r, ch


def null_edges(m, groups, seed):
    rng = np.random.default_rng(seed)
    p, q, ww = pre[m], post[m].copy(), w[m]
    if groups is None:
        q = q[rng.permutation(len(q))]
    else:
        g = groups[m]
        o0 = np.argsort(g, kind='stable'); o1 = np.lexsort((rng.random(len(g)), g))
        q2 = q.copy(); q2[o0] = q[o1]; q = q2
    key = p.astype(np.int64) * N + q
    u, inv2 = np.unique(key, return_inverse=True)
    return (u // N).astype(np.int64), (u % N).astype(np.int64), np.bincount(inv2, weights=ww)


res = {'note': 'EXPLORATORY', 'runs': []}
m1 = w >= 1
r, ch = run(pre[m1], post[m1], w[m1], 'observed, >= 1 synapse', True)
k = ch['keep']
r['original_cut'] = {'stats': sl.set_stats(ch, in_cut[k]), 'exit': sl.exit_channel(ch, in_cut[k])}
res['runs'].append(r); json.dump(res, open(OUT, 'w'), indent=1, default=_json_default)
m5 = w >= 5
r, _ = run(pre[m5], post[m5], w[m5], 'observed, >= 5 synapses', True)
res['runs'].append(r); json.dump(res, open(OUT, 'w'), indent=1, default=_json_default)

sc_code = np.unique(ann['superclass'].astype(str), return_inverse=True)[1]
typ_g = np.where(typ == '?', np.char.add('body', body_ids.astype(str)), typ)
ty_code = np.unique(typ_g, return_inverse=True)[1]
groups = {'N1 configuration': None,
          'N2 superclass blocks': sc_code[pre].astype(np.int64) * (sc_code.max() + 1) + sc_code[post],
          'N3 cell-type blocks': ty_code[pre].astype(np.int64) * (ty_code.max() + 1) + ty_code[post]}
for name, reps in (('N3 cell-type blocks', 3), ('N2 superclass blocks', 2), ('N1 configuration', 1)):
    for s in range(reps):
        r, _ = run(*null_edges(m1, groups[name], 1000 + 10 * s + 1), f'{name}, seed {1000 + 10 * s + 1}', False)
        res['runs'].append(r); json.dump(res, open(OUT, 'w'), indent=1, default=_json_default)
log('saved ' + OUT)
