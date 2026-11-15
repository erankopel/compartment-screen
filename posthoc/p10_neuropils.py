"""P10. POST HOC, EXPLORATORY (step 2: measurement or biology; not confirmatory).

Is the shortfall of detected synapses regional? One pass over the per-neuropil edge list of BANC
(neuron_connectivity/v888/banc_888_edgelist_neuropil_v3.feather: every pre and post segment, identified or not,
with the neuropil of the synapses). Produces:
  a. per neuropil: synapses, and the share whose presynaptic and whose postsynaptic end is a graph neuron;
  b. per cell type present in the three data sets (P9): the neuropil holding most of its output synapses; then,
     per neuropil, the median over its types of the ratio BANC : male and BANC : FlyWire of synapses made per
     neuron (BANC counted over all partners);
  c. the same ratios beside the recall of the BANC detector in the densely evaluated cutouts
     (251013_synapse_evaluation.csv of the BANC Dataverse deposit), neuropil by neuropil;
  d. the ratios by BANC cell class and for groups of central-complex and mushroom-body types.
Usage: python3 p10_neuropils.py <edgelist_neuropil feather> <251013_synapse_evaluation.csv>   ASCII, no em dashes.
"""
import csv, json, os, re, sys, time
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.feather as feather
import load3

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')
path, evalcsv = sys.argv[1], sys.argv[2]


def log(s):
    print(s, flush=True)


d = load3.load('banc')
N, ctype = d['N'], d['ctype']
del d
meta = feather.read_table(load3.PATHS['banc'][1], columns=['banc_888_id', 'super_class', 'cell_class', 'cell_type']).to_pandas()
sc = meta['super_class'].astype(object)
ok = sc.notna() & ~sc.isin(load3.CROSS['graph_exclude_super_class'])
ids = np.sort(meta[ok.to_numpy()]['banc_888_id'].astype(str).astype(np.int64).to_numpy())
assert len(ids) == N

rd = pa.ipc.open_file(pa.memory_map(path))
names, code = [], {}
out_np = None
tot = {}
t0 = time.time()
gg = 0
for i in range(rd.num_record_batches):
    b = rd.get_batch(i)
    p = pc.cast(b.column('pre_root_888'), pa.int64()).to_numpy(zero_copy_only=False)
    q = pc.cast(b.column('post_root_888'), pa.int64()).to_numpy(zero_copy_only=False)
    c = b.column('count').to_numpy(zero_copy_only=False).astype(np.int64)
    npl = b.column('neuropil').to_numpy(zero_copy_only=False).astype(object)
    keep = p != q
    p, q, c, npl = p[keep], q[keep], c[keep], npl[keep]
    u, inv = np.unique(np.array(['NA' if x is None else x for x in npl]), return_inverse=True)
    for x in u:
        if x not in code:
            code[x] = len(names); names.append(x)
    k = np.array([code[x] for x in u])[inv]
    ip = np.clip(np.searchsorted(ids, p), 0, N - 1); gp = ids[ip] == p
    iq = np.clip(np.searchsorted(ids, q), 0, N - 1); gq = ids[iq] == q
    if out_np is None:
        out_np = np.zeros((N, 256), np.int64)
    np.add.at(out_np, (ip[gp], k[gp]), c[gp])
    for kk in np.unique(k):
        m = k == kk
        t = tot.setdefault(int(kk), np.zeros(4, np.int64))
        t += np.array([c[m].sum(), c[m & gp].sum(), c[m & gq].sum(), c[m & gp & gq].sum()])
    gg += int(c[gp & gq].sum())
    if i % 400 == 0:
        log(f'  batch {i} of {rd.num_record_batches}, {time.time() - t0:.0f} s')
out_np = out_np[:, :len(names)]
log(f'neuropils {len(names)}; synapses between graph neurons {gg:,} (registered edge list 39,297,815); {time.time() - t0:.0f} s')
res = {'note': 'POST HOC, EXPLORATORY', 'graph_graph_synapses': gg, 'neuropils': {}}
for nm in names:
    t = tot[code[nm]]
    res['neuropils'][nm] = {'synapses': int(t[0]), 'pre_graph_share': float(t[1] / t[0]), 'post_graph_share': float(t[2] / t[0]),
                            'both_graph_share': float(t[3] / t[0])}

# dominant output neuropil of each matched type
T = json.load(open(os.path.join(R, 'p9_types.json')))['types']
u, inv = np.unique(ctype, return_inverse=True)
tsum = np.zeros((len(u), len(names)))
np.add.at(tsum, inv, out_np)
dom = {}
for i, t in enumerate(u):
    if t in T and tsum[i].sum() > 0:
        j = int(np.argmax(tsum[i]))
        dom[t] = (names[j], float(tsum[i, j] / tsum[i].sum()))
allr = {k: np.array([T[t][k] for t in T]) for k in ('out_BM', 'out_BF', 'out_FM', 'in_BM', 'in_BF')}


def rank(k, v):
    return float((allr[k] < v).mean())


log('\na, b. neuropils with at least 8 matched types that place at least half of their output there:')
log('   neuropil: synapses, share with identified pre, share with identified post | types | median out B:M (rank), out B:F (rank), F:M')
tab = {}
for nm in names:
    ts = [t for t, (n_, s_) in dom.items() if n_ == nm and s_ >= 0.5]
    if len(ts) >= 8:
        v = {k: float(np.median([T[t][k] for t in ts])) for k in ('out_BM', 'out_BF', 'out_FM', 'in_BM', 'in_BF')}
        tab[nm] = {'types': len(ts), **v, 'out_BM_rank': rank('out_BM', v['out_BM']), 'out_BF_rank': rank('out_BF', v['out_BF']), **res['neuropils'][nm]}
for nm, v in sorted(tab.items(), key=lambda kv: kv[1]['out_BM']):
    log(f'   {nm:<10} {v["synapses"]:>11,} {v["pre_graph_share"]:.2f} {v["post_graph_share"]:.2f} | {v["types"]:>4} | {v["out_BM"]:.2f} ({v["out_BM_rank"]:.2f})  {v["out_BF"]:.2f} ({v["out_BF_rank"]:.2f})  {v["out_FM"]:.2f}')
res['by_dominant_neuropil'] = tab

log('\nc. densely evaluated cutouts: recall of the detector beside the type ratios of the same neuropil')
ev = list(csv.DictReader(open(evalcsv)))
MAP = {'ME': ['ME'], 'LO': ['LO'], 'LAL': ['LAL'], 'SPS': ['SPS'], 'GNG': ['GNG'], 'FB': ['FB'], 'SMP': ['SMP'], 'AVLP': ['AVLP'],
       'MB_ML': ['MB_ML', 'MB_VL', 'MB_PED', 'MB'], 'MB_CA': ['MB_CA']}
cut = []
for r in ev:
    key = re.sub(r'_(left|right)$', '', r['neuropil'])
    cand = [nm for nm in names if any(nm == x or nm.startswith(x + '_') for x in MAP.get(key, [key]))]
    ts = [t for t, (n_, s_) in dom.items() if n_ in cand and s_ >= 0.5]
    row = {'cutout': r['cutouts'], 'neuropil': r['neuropil'], 'region': r['region'], 'recall': float(r['recall']), 'precision': float(r['precision']),
           'matched_neuropil_names': cand, 'types': len(ts)}
    if len(ts) >= 5:
        row['out_BM'] = float(np.median([T[t]['out_BM'] for t in ts])); row['out_BF'] = float(np.median([T[t]['out_BF'] for t in ts]))
    cut.append(row)
    log(f'   {r["neuropil"]:<22} recall {row["recall"]:.2f} precision {row["precision"]:.2f} | types {len(ts):>4} | ' +
        (f'out B:M {row["out_BM"]:.2f}  out B:F {row["out_BF"]:.2f}' if 'out_BM' in row else 'too few matched types') + f' | {cand[:4]}')
res['cutouts'] = cut
ok_rows = [r for r in cut if 'out_BM' in r]
if len(ok_rows) >= 4:
    def spearman(a, b):
        ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
        return float(np.corrcoef(ra, rb)[0, 1])
    rc = np.array([r['recall'] for r in ok_rows])
    res['cutout_rank_correlation'] = {'n': len(ok_rows), 'recall_vs_out_BM': spearman(rc, np.array([r['out_BM'] for r in ok_rows])),
                                      'recall_vs_out_BF': spearman(rc, np.array([r['out_BF'] for r in ok_rows]))}
    log(f'   rank correlation over {len(ok_rows)} cutouts: recall with out B:M {res["cutout_rank_correlation"]["recall_vs_out_BM"]:.2f}, '
        f'with out B:F {res["cutout_rank_correlation"]["recall_vs_out_BF"]:.2f}')

log('\nd. by BANC cell class (at least 8 matched types), the 10 lowest by out B:M; and groups by type name')
g = meta[ok.to_numpy()]
cls_of = g.groupby(g['cell_type'].astype(str))['cell_class'].agg(lambda s: s.astype(object).dropna().value_counts().index[0] if s.notna().any() else '?').to_dict()
cl = np.array([cls_of.get(t, '?') for t in T])
rows = []
for c in sorted(set(cl)):
    mk = cl == c
    if mk.sum() >= 8:
        rows.append((float(np.median(allr['out_BM'][mk])), c, int(mk.sum()), float(np.median(allr['out_BF'][mk]))))
rows.sort()
res['by_class'] = {c: {'types': n, 'out_BM': m, 'out_BM_rank': rank('out_BM', m), 'out_BF': bf, 'out_BF_rank': rank('out_BF', bf)} for m, c, n, bf in rows}
for m, c, n, bf in rows[:10]:
    log(f'   {c:<42} {n:>4} types  out B:M {m:.2f} (rank {rank("out_BM", m):.2f})  out B:F {bf:.2f} (rank {rank("out_BF", bf):.2f})')
GROUPS = {'ring neurons (ER)': r'^ER[1-6]', 'ExR': r'^ExR', 'EPG, PEG, PEN, EL': r'^(EPG|PEG|PEN_|EL$)', 'FB tangential': r'^FB[1-9]',
          'hDelta, vDelta': r'^(hDelta|vDelta)', 'FC, FR, FS': r'^(FC|FR|FS)[1-9]', 'PFN, PFL, PFG, PFR': r'^(PFN|PFL|PFG|PFR)',
          'LNO, GLNO, LCNO': r'^(LNO|GLNO|LCNO)', 'TuBu': r'^TuBu', 'MeTu': r'^MeTu', 'Kenyon cells': r'^KC', 'MBON': r'^MBON',
          'dopaminergic PAM, PPL': r'^(PAM|PPL)'}
names_T = list(T)
res['groups'] = {}
for nm, pat in GROUPS.items():
    mk = np.array([bool(re.match(pat, t)) for t in names_T])
    if mk.sum() >= 2:
        m, bf, fm = (float(np.median(allr[k][mk])) for k in ('out_BM', 'out_BF', 'out_FM'))
        res['groups'][nm] = {'types': int(mk.sum()), 'out_BM': m, 'out_BM_rank': rank('out_BM', m), 'out_BF': bf, 'out_BF_rank': rank('out_BF', bf), 'out_FM': fm}
        log(f'   {nm:<24} {int(mk.sum()):>4} types  out B:M {m:.2f} (rank {rank("out_BM", m):.2f})  out B:F {bf:.2f} (rank {rank("out_BF", bf):.2f})  F:M {fm:.2f}')
json.dump(res, open(os.path.join(R, 'p10_neuropils.json'), 'w'), indent=1)
log('saved')
