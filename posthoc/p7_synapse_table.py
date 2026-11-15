"""P7. POST HOC, EXPLORATORY (step 2: measurement or biology; not confirmatory).

One pass over the BANC v3 synapse table (259 million rows; file
synapses_v3_human_readable_id_size_prerootid_postrootid_prex_prey_prez.parquet of the BANC bucket). This parquet
holds only synapses of size >= 10 voxels (the size bins below 10 stay empty; the smaller detections are only in the
raw csv.gz, which was not used). Produces:
  - per graph neuron, the number of outgoing and incoming synapses by the kind of partner: 0 a graph neuron, 1 a
    row of the metadata table outside the graph, 2 no segment (root 0), 3 a segment that is not in the metadata
    table (a fragment nobody has identified)                                   -> results/p7_neuron_counts.npz
  - the neuron-to-neuron counts inside the graph (check against the registered edge list)
  - every synapse with an end in the EB set (ends, partner kinds, size, position) -> results/p7_eb_synapses.npz
Autapses are dropped, as in the distributed edge list. Check: the counts must equal the registered edge list.
Usage: python3 p7_synapse_table.py <synapse parquet>      ASCII, no em dashes.
"""
import json, os, sys, time
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.feather as feather
import load3

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results')
path = sys.argv[1]


def log(s):
    print(s, flush=True)


d = load3.load('banc')
N, EB = d['N'], d['sets']['EB']
meta = feather.read_table(load3.PATHS['banc'][1], columns=['banc_888_id', 'super_class']).to_pandas()
all_ids = np.sort(meta['banc_888_id'].astype(str).astype(np.int64).to_numpy())
sc = meta['super_class'].astype(object)
ok = sc.notna() & ~sc.isin(load3.CROSS['graph_exclude_super_class'])
ids = np.sort(meta[ok.to_numpy()]['banc_888_id'].astype(str).astype(np.int64).to_numpy())
assert len(ids) == N
del d                                                        # the registered edge list is loaded again for the check


def lookup(r):
    """graph index (or -1) and partner kind 0..3 of root ids r"""
    i = np.clip(np.searchsorted(ids, r), 0, N - 1)
    g = np.where(ids[i] == r, i, -1)
    j = np.clip(np.searchsorted(all_ids, r), 0, len(all_ids) - 1)
    inmeta = all_ids[j] == r
    kind = np.where(g >= 0, 0, np.where(inmeta, 1, np.where(r == 0, 2, 3))).astype(np.int8)
    return g, kind


out_c = np.zeros(N * 12, np.int64)      # [neuron, size bin, partner kind]
in_c = np.zeros(N * 12, np.int64)
tot = np.zeros((3, 4, 4), np.int64)     # [size bin, kind of pre, kind of post], all synapses
K, Cn = np.zeros(0, np.int64), np.zeros(0, np.int32)       # sorted keys (pair * 3 + size bin) and their counts
ebrows = []


def merge(K, Cn, u, c):
    """add the sorted unique keys u with counts c into the sorted arrays K, Cn"""
    pos = np.searchsorted(K, u)
    hit = (pos < len(K)) & (K[np.minimum(pos, max(len(K) - 1, 0))] == u) if len(K) else np.zeros(len(u), bool)
    Cn[pos[hit]] += c[hit]
    if (~hit).any():
        K = np.insert(K, pos[~hit], u[~hit]); Cn = np.insert(Cn, pos[~hit], c[~hit])
    return K, Cn


pf = pq.ParquetFile(path)
t0, nrow, naut = time.time(), 0, 0
for bi, b in enumerate(pf.iter_batches(batch_size=2_000_000, columns=['size', 'pre_root_id', 'post_root_id', 'pre_x', 'pre_y', 'pre_z'])):
    sz = b.column('size').to_numpy(zero_copy_only=False)
    pr = b.column('pre_root_id').to_numpy(zero_copy_only=False)
    po = b.column('post_root_id').to_numpy(zero_copy_only=False)
    nrow += len(sz)
    keep = pr != po
    naut += int((~keep).sum())
    sz, pr, po = sz[keep], pr[keep], po[keep]
    x, y, z = (b.column(c).to_numpy(zero_copy_only=False)[keep] for c in ('pre_x', 'pre_y', 'pre_z'))
    gp, kp = lookup(pr)
    gq, kq = lookup(po)
    sb = np.where(sz >= 10, 2, np.where(sz >= 5, 1, 0)).astype(np.int8)
    np.add.at(tot, (sb, kp, kq), 1)
    m = gp >= 0
    out_c += np.bincount(gp[m].astype(np.int64) * 12 + sb[m] * 4 + kq[m], minlength=N * 12)
    m = gq >= 0
    in_c += np.bincount(gq[m].astype(np.int64) * 12 + sb[m] * 4 + kp[m], minlength=N * 12)
    m = (gp >= 0) & (gq >= 0)
    k3 = (gp[m].astype(np.int64) * N + gq[m]) * 3 + sb[m]
    u, c = np.unique(k3, return_counts=True)
    K, Cn = merge(K, Cn, u, c.astype(np.int32))
    del u, c, k3
    e = ((gp >= 0) & EB[np.clip(gp, 0, N - 1)]) | ((gq >= 0) & EB[np.clip(gq, 0, N - 1)])
    if e.any():
        ebrows.append(np.stack([gp[e], gq[e], kp[e], kq[e], sz[e], x[e], y[e], z[e]], axis=1).astype(np.int64))
    del b, sz, pr, po, x, y, z, gp, gq, kp, kq, sb, m, e
    if bi % 20 == 0:
        pa.default_memory_pool().release_unused()
        log(f'  batch {bi}: {nrow:,} rows, {len(K):,} keys, {time.time() - t0:.0f} s')
pair, sbin = K // 3, K % 3
up, inv = np.unique(pair, return_inverse=True)
C = np.zeros((len(up), 3), np.int32)
C[inv, sbin] = Cn
del K, Cn, pair, sbin, inv
reg = load3.load('banc')
reg_pre, reg_post, reg_w = reg['pre'], reg['post'], reg['w']
del reg
e_pre, e_post = (up // N).astype(np.int32), (up % N).astype(np.int32)
log(f'rows {nrow:,}; autapses dropped {naut:,}; graph pairs {len(up):,}; {time.time() - t0:.0f} s')

# check against the registered edge list (size >= 10)
m10 = C[:, 2] > 0
reg_key = reg_pre.astype(np.int64) * N + reg_post
o = np.argsort(reg_key)
same_pairs = bool(m10.sum() == len(reg_key) and np.array_equal(up[m10], reg_key[o]))
same_counts = bool(same_pairs and np.array_equal(C[m10, 2], reg_w[o]))
log(f'check against the registered edge list: pairs at size >= 10: {int(m10.sum()):,} (registered {len(reg_key):,}); '
    f'synapses {int(C[:, 2].sum()):,} (registered {int(reg_w.sum()):,}); identical pairs {same_pairs}; identical counts {same_counts}')

np.savez_compressed(os.path.join(OUT, 'p7_neuron_counts.npz'), out=out_c.reshape(N, 3, 4), inc=in_c.reshape(N, 3, 4), total=tot)
E = np.concatenate(ebrows) if ebrows else np.zeros((0, 8), np.int64)
np.savez_compressed(os.path.join(OUT, 'p7_eb_synapses.npz'), rows=E,
                    columns=np.array(['pre_index', 'post_index', 'pre_kind', 'post_kind', 'size', 'x16', 'y16', 'z45']))
summ = {'note': 'POST HOC, EXPLORATORY', 'rows': int(nrow), 'autapses_dropped': int(naut), 'graph_pairs_any_size': int(len(up)),
        'check_registered': {'pairs_ge10': int(m10.sum()), 'registered_pairs': int(len(reg_key)), 'synapses_ge10': int(C[:, 2].sum()),
                             'registered_synapses': int(reg_w.sum()), 'identical_pairs': same_pairs, 'identical_counts': same_counts},
        'graph_synapses_by_size_bin': {'lt5': int(C[:, 0].sum()), '5to9': int(C[:, 1].sum()), 'ge10': int(C[:, 2].sum())},
        'all_synapses_by_size_bin_prekind_postkind': tot.tolist(), 'eb_rows': int(len(E))}
json.dump(summ, open(os.path.join(OUT, 'p7_synapse_table.json'), 'w'), indent=1)
log(json.dumps({k: v for k, v in summ.items() if k != 'all_synapses_by_size_bin_prekind_postkind'}))
log('saved')
