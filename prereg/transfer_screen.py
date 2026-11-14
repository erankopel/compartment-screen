"""Transfer of the screen to another connectome in the compiled layout of the BANC data tutorial
(<name>_meta.feather with columns <name>_id, super_class, cell_class, cell_type, ...; an edge list with
columns pre, post, count). Sets are defined by names only, from crosswalk_types.json:
  EB  neurons whose cell_type matches a family rule (EB_family_regex; any comma-separated part of the name)
  CX  the three central-complex classes, plus EB
  MB  Kenyon cells, mushroom-body output neurons, mushroom-body dopaminergic neurons
  AL  antennal-lobe projection and local neurons
  amine  neurons whose cell_type is in amine_types_flywire_names
Graph: every row of the metadata table whose super_class is not null and not excluded in the crosswalk
(glia, trachea, not_a_neuron); edges with both ends in it.
Usage: python3 transfer_screen.py <id column> <meta.feather> <edgelist.feather | NONE> <out.json> [min count ...]
With NONE the script reads names only (no transmitter column, no edges), prints the sets and stops.
Used on FlyWire v783 as an exploratory rehearsal. ASCII, no em dashes."""
import hashlib, json, os, re, sys, time
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.feather as feather
import screen_lib as sl

HERE = os.path.dirname(os.path.abspath(__file__))
idcol, meta_path, edge_path, out_path = sys.argv[1:5]
thresholds = [int(x) for x in sys.argv[5:]] or [1, 5]
cross = json.load(open(os.path.join(HERE, 'crosswalk_types.json')))


def log(s):
    print(s, flush=True)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def md5_b64(path, chunk=1 << 22):
    """MD5 in base64, the form in which the storage bucket reports it (x-goog-hash)."""
    import base64
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(chunk), b''):
            h.update(b)
    return base64.b64encode(h.digest()).decode()


META_ONLY = edge_path == 'NONE'
meta = feather.read_table(meta_path, columns=[idcol, 'super_class', 'cell_class', 'cell_type'] +
                          ([] if META_ONLY else ['neurotransmitter_predicted'])).to_pandas()
sc = meta['super_class'].astype(object)
ok = sc.notna() & ~sc.isin(cross['graph_exclude_super_class'])
meta = meta[ok.to_numpy()].copy()
ids = meta[idcol].astype(str).astype(np.int64).to_numpy()
order = np.argsort(ids); ids = ids[order]; meta = meta.iloc[order]
N = len(ids)
ctype = meta['cell_type'].astype(object).fillna('?').astype(str).to_numpy()
cclass = meta['cell_class'].astype(object).fillna('?').astype(str).to_numpy()
ntp = (np.full(N, '?') if META_ONLY else
       meta['neurotransmitter_predicted'].astype(object).fillna('?').astype(str).to_numpy())
_rx = [re.compile(r) for r in cross['EB_family_regex']]
EB = np.array([any(r.match(p.strip()) for p in t.split(',') for r in _rx) for t in ctype])
SETS = {'EB': EB, 'CX': np.isin(cclass, cross['class_map']['CX']) | EB,
        'MB': np.isin(cclass, cross['class_map']['MB']), 'AL': np.isin(cclass, cross['class_map']['AL'])}
amine = np.isin(ctype, cross['amine_types_flywire_names'])
amine_pred = np.isin(ntp, ['dopamine', 'serotonin', 'octopamine'])
log(f'graph nodes {N:,}; sets: ' + ', '.join(f'{k} {int(v.sum())}' for k, v in SETS.items()) +
    f'; amine by type list {int(amine.sum())}, by this data set\'s own prediction {int(amine_pred.sum())}')
present = sorted(set(ctype[EB])); missing = [r for r in cross['EB_family_regex'] if not any(re.match(r, p.strip()) for t in present for p in t.split(','))]
log(f'EB types selected: {len(present)}: {present}; family rules with no match: {missing}')
for cl in sum(cross['class_map'].values(), []):
    log(f'  class {cl}: {int((cclass == cl).sum())}')
if META_ONLY:
    log('names only: stopping before any edge is read')
    sys.exit(0)

rd = pa.ipc.open_file(pa.memory_map(edge_path))
P_, Q_, W_ = [], [], []
for i in range(rd.num_record_batches):
    b = rd.get_batch(i)
    p = pc.cast(b.column('pre'), pa.int64()).to_numpy(zero_copy_only=False)
    q = pc.cast(b.column('post'), pa.int64()).to_numpy(zero_copy_only=False)
    c = b.column('count').to_numpy(zero_copy_only=False)
    ip = np.clip(np.searchsorted(ids, p), 0, N - 1); iq = np.clip(np.searchsorted(ids, q), 0, N - 1)
    m = (ids[ip] == p) & (ids[iq] == q)
    P_.append(ip[m].astype(np.int32)); Q_.append(iq[m].astype(np.int32)); W_.append(c[m].astype(np.int32))
pre, post, w = np.concatenate(P_), np.concatenate(Q_), np.concatenate(W_)
log(f'edges inside the graph: {len(pre):,}; synapses {int(w.sum()):,}')

res = {'meta': meta_path, 'edges': edge_path, 'md5_base64': {'meta': md5_b64(meta_path), 'edges': md5_b64(edge_path)},
       'sha256': {'meta': sha256(meta_path), 'edges': sha256(edge_path),
       'screen_lib': sha256(os.path.join(HERE, 'screen_lib.py')), 'crosswalk': sha256(os.path.join(HERE, 'crosswalk_types.json')),
       'transfer_screen': sha256(os.path.abspath(__file__))}, 'nodes': int(N), 'set_sizes': {k: int(v.sum()) for k, v in SETS.items()},
       'EB_names_missing': missing, 'runs': []}
for th in thresholds:
    t0 = time.time()
    m = w >= th
    ch = sl.build_chain(pre[m], post[m], w[m], N)
    k = ch['keep']
    log(f'--- >= {th} synapses: SCC {ch["n"]:,} of {N:,}; stationary in {ch["iters"]} iterations (last change {ch["delta"]:.1e}); '
        f'EB neurons in the SCC {int(EB[k].sum())}')
    r = sl.screen(ch, {a: b[k] for a, b in SETS.items()}, amine[k], full=True, log=log)
    r['amine_own_prediction'] = {nm: sl.label_shares(ch, SETS[nm][k], amine_pred[k]) for nm in ('EB', 'MB')}
    r['threshold'] = th; r['wall_s'] = round(time.time() - t0)
    res['runs'].append(r)
    json.dump(res, open(out_path, 'w'), indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x))
log('saved ' + out_path)
