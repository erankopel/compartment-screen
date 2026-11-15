"""POST HOC, EXPLORATORY (written after the confirmatory run of 14 Nov 2026; nothing here is confirmatory).

Loader for the three data sets, with the set definitions of the registration (names only):
  banc     BANC release 888          compiled layout (meta feather + edge list)
  flywire  FlyWire v783              compiled layout
  male     male CNS v1.0             neuron_edges.npz + body-annotations.feather
Returns a dict: pre, post, w (int arrays over node indices 0..N-1), N, sets (EB, CX, MB, AL boolean masks),
fam (EB family of every node, '' outside EB), ctype, side, extra (data set specific columns).
The statistics come from the frozen screen_lib.py, imported unchanged. ASCII, no em dashes.
"""
import json, os, re, sys
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.feather as feather

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('SCREEN_DATA', os.path.join(HERE, '..'))
# frozen code: ../prereg in the repository; SCREEN_FROZEN overrides
FROZEN = os.environ.get('SCREEN_FROZEN') or next(
    (d for d in (os.path.join(HERE, '..', 'prereg'), os.path.join(BASE, 'prereg_screen_banc')) if os.path.isdir(d)), '')
# data folders (override with BANC_DIR, FLYWIRE_DIR, MALE_DIR); file names as distributed
BANC_DIR = os.environ.get('BANC_DIR', os.path.join(BASE, 'confirm', 'run'))
FLYWIRE_DIR = os.environ.get('FLYWIRE_DIR', os.path.join(BASE, 'data_fafb'))
MALE_DIR = os.environ.get('MALE_DIR', os.path.join(BASE, 'data'))
sys.path.insert(0, FROZEN)
import screen_lib as sl  # noqa: E402  (frozen, sha256 e062d134...)

CROSS = json.load(open(os.path.join(FROZEN, 'crosswalk_types.json')))
RX = [re.compile(r) for r in CROSS['EB_family_regex']]
FAM_NAMES = ['ER', 'EPG', 'PEG', 'PEN', 'EL', 'ExR', 'TuBu', 'MeTu2', 'LPsP', 'GLNO', 'AOTU046', 'TuTuB']
# BANC_EDGE_FILE: another edge list of the same release (e.g. the v2 synapse set); default is the registered v3 list
BANC_EDGES = os.environ.get('BANC_EDGE_FILE', os.path.join(BANC_DIR, 'banc_888_edgelist_simple_v3.feather'))
PATHS = {'banc': ('banc_888_id', os.path.join(BANC_DIR, 'banc_888_meta.feather'), BANC_EDGES),
         'flywire': ('fafb_783_id', os.path.join(FLYWIRE_DIR, 'fafb_783_meta.feather'), os.path.join(FLYWIRE_DIR, 'fafb_783_simple_edgelist.feather'))}


def family(t):
    """Index of the first family rule matched by any comma-separated part of the type name, or -1."""
    for p in t.split(','):
        for i, r in enumerate(RX):
            if r.match(p.strip()):
                return i
    return -1


def _fam_array(ctype):
    idx = np.array([family(t) for t in ctype])
    return idx >= 0, np.array([FAM_NAMES[i] if i >= 0 else '' for i in idx])


def load(name):
    if name == 'male':
        d = np.load(os.path.join(MALE_DIR, 'neuron_edges.npz'))
        ids = d['body_ids']
        N = len(ids)
        ann = feather.read_table(os.path.join(MALE_DIR, 'body-annotations.feather'),
                                 columns=['bodyId', 'type', 'class', 'somaSide', 'superclass']).to_pandas()
        ann = ann.drop_duplicates('bodyId').set_index('bodyId').reindex(ids)
        ctype = ann['type'].astype(object).fillna('?').astype(str).to_numpy()
        cls = ann['class'].astype(object).fillna('?').astype(str).to_numpy()
        side = ann['somaSide'].astype(object).fillna('?').astype(str).to_numpy()
        EB, fam = _fam_array(ctype)
        sets = {'EB': EB, 'CX': (cls == 'CX') | EB, 'MB': np.isin(cls, ['Kenyon_Cell', 'MBON', 'DAN']),
                'AL': np.isin(cls, ['ALPN', 'ALLN'])}
        return {'name': name, 'pre': d['pre'].astype(np.int32), 'post': d['post'].astype(np.int32),
                'w': d['w'].astype(np.int32), 'N': N, 'sets': sets, 'fam': fam, 'ctype': ctype, 'side': side,
                'extra': {}}
    idcol, mp, ep = PATHS[name]
    cols = [idcol, 'super_class', 'cell_class', 'cell_type', 'side']
    if name == 'banc':
        cols += ['proofread', 'status', 'input_connections', 'output_connections']
    meta = feather.read_table(mp, columns=cols).to_pandas()
    sc = meta['super_class'].astype(object)
    ok = sc.notna() & ~sc.isin(CROSS['graph_exclude_super_class'])
    meta = meta[ok.to_numpy()].copy()
    ids = meta[idcol].astype(str).astype(np.int64).to_numpy()
    order = np.argsort(ids); ids = ids[order]; meta = meta.iloc[order]
    N = len(ids)
    ctype = meta['cell_type'].astype(object).fillna('?').astype(str).to_numpy()
    cclass = meta['cell_class'].astype(object).fillna('?').astype(str).to_numpy()
    side = meta['side'].astype(object).fillna('?').astype(str).to_numpy()
    EB, fam = _fam_array(ctype)
    sets = {'EB': EB, 'CX': np.isin(cclass, CROSS['class_map']['CX']) | EB,
            'MB': np.isin(cclass, CROSS['class_map']['MB']), 'AL': np.isin(cclass, CROSS['class_map']['AL'])}
    extra = {}
    if name == 'banc':
        for c in ('proofread', 'status', 'input_connections', 'output_connections'):
            extra[c] = meta[c].to_numpy()
    rd = pa.ipc.open_file(pa.memory_map(ep))
    P_, Q_, W_ = [], [], []
    for i in range(rd.num_record_batches):
        b = rd.get_batch(i)
        p = pc.cast(b.column('pre'), pa.int64()).to_numpy(zero_copy_only=False)
        q = pc.cast(b.column('post'), pa.int64()).to_numpy(zero_copy_only=False)
        c = b.column('count').to_numpy(zero_copy_only=False)
        ip = np.clip(np.searchsorted(ids, p), 0, N - 1); iq = np.clip(np.searchsorted(ids, q), 0, N - 1)
        m = (ids[ip] == p) & (ids[iq] == q)
        P_.append(ip[m].astype(np.int32)); Q_.append(iq[m].astype(np.int32)); W_.append(c[m].astype(np.int32))
    return {'name': name, 'pre': np.concatenate(P_), 'post': np.concatenate(Q_), 'w': np.concatenate(W_), 'N': N,
            'sets': sets, 'fam': fam, 'ctype': ctype, 'side': side, 'extra': extra}


def conditioned_Q(ch, M):
    """Transition matrix of the walk conditioned to stay in M (the Doob transform used by the frozen
    dobrushin_conditioned), with the indices (into ch['keep']) of the neurons kept."""
    idx = np.flatnonzero(M)
    D = ch['P'][idx][:, idx].toarray()
    ev, V = np.linalg.eig(D)
    k = int(np.argmax(ev.real))
    l1, h = float(ev[k].real), np.abs(V[:, k].real)
    ok = h > 1e-10 * h.max()
    D, h = D[ok][:, ok], h[ok]
    Q = D * h[None, :] / (l1 * h[:, None])
    Q /= Q.sum(axis=1, keepdims=True)
    return Q, idx[ok], l1
