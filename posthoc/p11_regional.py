"""P11. POST HOC, EXPLORATORY (step 2: measurement or biology; not confirmatory).

Is the block pattern of P4 a regional pattern? The synapses between neurons of the EB set are counted neuropil by
neuropil in BANC (per-neuropil edge list) and in FlyWire (per-synapse table fafb_783_synapses.parquet of the BANC
bucket, which holds 79 percent of the synapses of the FlyWire edge list; its counts are scaled to the edge list).
Then FlyWire is thinned REGIONALLY: a synapse inside the EB set is kept with the BANC : FlyWire ratio of its
neuropil, every other synapse with p_out as in P4. No family information is used. Three seeds; thresholds 1, 2, 3, 5.
The male data have no neuropil labels in the files used here, so the test is made on FlyWire only.
Usage: python3 p11_regional.py <banc edgelist_neuropil feather> <fafb_783_synapses.parquet>   ASCII, no em dashes.
"""
import json, os, sys, time
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyarrow.feather as feather
import load3
from load3 import sl

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')
banc_np, fw_syn = sys.argv[1], sys.argv[2]
THS, SEEDS = [1, 2, 3, 5], [20261115, 20261116, 20261117]


def log(s):
    print(s, flush=True)


def graph_ids(name):
    idcol, mp, _ = load3.PATHS[name]
    meta = feather.read_table(mp, columns=[idcol, 'super_class']).to_pandas()
    sc = meta['super_class'].astype(object)
    ok = sc.notna() & ~sc.isin(load3.CROSS['graph_exclude_super_class'])
    return np.sort(meta[ok.to_numpy()][idcol].astype(str).astype(np.int64).to_numpy())


# 1. BANC: synapses inside the EB set by neuropil
d = load3.load('banc')
N, EB = d['N'], d['sets']['EB']
banc_inside = int(d['w'][EB[d['pre']] & EB[d['post']]].sum())
del d
ids = graph_ids('banc')
rd = pa.ipc.open_file(pa.memory_map(banc_np))
B = {}
for i in range(rd.num_record_batches):
    b = rd.get_batch(i)
    p = pc.cast(b.column('pre_root_888'), pa.int64()).to_numpy(zero_copy_only=False)
    q = pc.cast(b.column('post_root_888'), pa.int64()).to_numpy(zero_copy_only=False)
    ip = np.clip(np.searchsorted(ids, p), 0, N - 1); iq = np.clip(np.searchsorted(ids, q), 0, N - 1)
    e = (ids[ip] == p) & (ids[iq] == q) & EB[ip] & EB[iq] & (p != q)
    if e.any():
        c = b.column('count').to_numpy(zero_copy_only=False)[e]
        nl = b.column('neuropil').to_numpy(zero_copy_only=False)[e]
        for x, cc in zip(nl, c):
            x = 'NA' if x is None else x
            B[x] = B.get(x, 0) + int(cc)
log(f'BANC: synapses inside EB by neuropil sum to {sum(B.values()):,} (edge list {banc_inside:,})')

# 2. FlyWire: per-pair neuropil composition inside the EB set
d = load3.load('flywire')
N, EB, pre, post, w = d['N'], d['sets']['EB'], d['pre'], d['post'], d['w']
inEB = EB[pre] & EB[post]
fw_inside = int(w[inEB].sum())
ids = graph_ids('flywire')
pf = pq.ParquetFile(fw_syn)
P_, Q_, L_ = [], [], []
nrows = 0
for b in pf.iter_batches(batch_size=2_000_000, columns=['pre', 'post', 'neuropil']):
    p = pc.cast(b.column('pre'), pa.int64()).to_numpy(zero_copy_only=False)
    q = pc.cast(b.column('post'), pa.int64()).to_numpy(zero_copy_only=False)
    nrows += len(p)
    ip = np.clip(np.searchsorted(ids, p), 0, N - 1); iq = np.clip(np.searchsorted(ids, q), 0, N - 1)
    e = (ids[ip] == p) & (ids[iq] == q) & EB[ip] & EB[iq] & (p != q)
    if e.any():
        P_.append(ip[e]); Q_.append(iq[e])
        L_.append(np.array(['NA' if x is None else x for x in b.column('neuropil').to_numpy(zero_copy_only=False)[e]]))
fp, fq, fl = np.concatenate(P_), np.concatenate(Q_), np.concatenate(L_)
scale = fw_inside / len(fp)
log(f'FlyWire: synapse table rows {nrows:,}; inside EB {len(fp):,} rows against {fw_inside:,} synapses in the edge list (scale {scale:.3f})')
u, c = np.unique(fl, return_counts=True)
F = dict(zip(u.tolist(), (c * scale).tolist()))

# 3. ratios by neuropil
MAJOR = [x for x in F if F[x] >= 3000 and x in B and x not in ('NA', 'meshes', 'outside_neuropil', 'outside')]
oth_b = sum(v for k, v in B.items() if k not in MAJOR); oth_f = sum(v for k, v in F.items() if k not in MAJOR)
ratio = {x: B[x] / F[x] for x in MAJOR}
ratio['other'] = oth_b / oth_f
p_in = sum(B.values()) / fw_inside
log('synapses inside the EB set by neuropil: BANC, FlyWire (scaled), ratio')
for x in sorted(MAJOR, key=lambda k: -F[k]):
    log(f'   {x:<8} {B[x]:>9,} {F[x]:>11,.0f}  {ratio[x]:.3f}')
log(f'   other    {oth_b:>9,} {oth_f:>11,.0f}  {ratio["other"]:.3f}')
log(f'   all      {sum(B.values()):>9,} {fw_inside:>11,}  {p_in:.3f}')
res = {'note': 'POST HOC, EXPLORATORY', 'banc_by_neuropil': B, 'flywire_by_neuropil_scaled': F, 'ratio': ratio, 'p_inside_overall': p_in,
       'flywire_table_coverage': len(fp) / fw_inside, 'runs': []}

# 4. regional thinning of FlyWire
SC = json.load(open(os.path.join(R, 'p2_scale.json')))['data']
p_out = SC['banc']['synapses_per_neuron'] / SC['flywire']['synapses_per_neuron']
pr = np.array([min(1.0, ratio.get(x, ratio['other'])) for x in fl])
key_tab = fp.astype(np.int64) * N + fq
o = np.argsort(key_tab, kind='stable'); key_tab, pr = key_tab[o], pr[o]
st = np.flatnonzero(np.r_[True, key_tab[1:] != key_tab[:-1]])
pair_keys = key_tab[st]
pair_p = np.add.reduceat(pr, st) / np.diff(np.r_[st, len(pr)])      # mean keep probability of the pair's synapses
ekey = pre.astype(np.int64) * N + post
pv = np.full(len(w), p_out)
j = np.clip(np.searchsorted(pair_keys, ekey[inEB]), 0, len(pair_keys) - 1)
found = pair_keys[j] == ekey[inEB]
pin = np.where(found, pair_p[j], min(1.0, ratio['other']))
pv[inEB] = pin
log(f'edges inside EB: {int(inEB.sum()):,}; with neuropil information {int(found.sum()):,}; synapse-weighted mean keep probability inside EB '
    f'{float((pin * w[inEB]).sum() / w[inEB].sum()):.3f} (uniform variant of P4: {SC["banc"]["EB"]["synapses_inside_per_neuron"] / SC["flywire"]["EB"]["synapses_inside_per_neuron"]:.3f})')
for seed in SEEDS:
    rng = np.random.default_rng(seed)
    wt = rng.binomial(w, pv).astype(np.int32)
    for th in THS:
        t0 = time.time()
        m = wt >= th
        ch = sl.build_chain(pre[m], post[m], wt[m], N, maxit=300)
        k = ch['keep']
        st_ = sl.set_stats(ch, EB[k]); db = sl.dobrushin_conditioned(ch, EB[k]); ex = sl.exit_channel(ch, EB[k])
        res['runs'].append({'seed': seed, 'threshold': th, 'EB_in_scc': int(EB[k].sum()), 'stats': st_, 'dobrushin': db, 'exit': ex})
        log(f'  regional seed {seed} th {th}: EB {int(EB[k].sum())} | loss {st_["escape"]:.4f} stays {st_["stays"]:.2f} mixes {st_["mixes"]:.2f} ratio {st_["ratio"]:.3f} | '
            f'Dob8 {db["8"]:.3f} | exit mean {ex["mean"]:.3f} below {ex["share_below_half"]:.3f} | {time.time() - t0:.0f} s')
        json.dump(res, open(os.path.join(R, 'p11_regional.json'), 'w'), indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x))
log('saved')
