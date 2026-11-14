"""T4a. The boundary channel with a frame mismatch: exact contraction coefficients and a numeric check.

Model (Fields, Glazebrook and Levin: a bit is prepared on the boundary with one reference frame and
measured with another). One writer prepares a qubit in its z basis, one reader measures in a basis
rotated by theta. Claim to check:
  flip probability p = sin^2(theta/2), a binary symmetric channel;
  Dobrushin coefficient eta_TV = |cos theta|; KL contraction coefficient eta_KL = cos^2 theta;
  r rounds: eta_TV = |cos theta|^r, eta_KL = cos^(2r) theta;
  on a screen of N independent bits eta_KL stays cos^2 theta, eta_TV does not.
ASCII, no em dashes.
"""
import json
import math
import os

import numpy as np

rng = np.random.default_rng(20261114)
res = {}


def bsc(p):
    return np.array([[1 - p, p], [p, 1 - p]])


def kl(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = a > 0
    return float((a[m] * np.log(a[m] / b[m])).sum())


def tv(a, b):
    return 0.5 * float(np.abs(np.asarray(a) - np.asarray(b)).sum())


# 1. the quantum step: prepare |0> or |1>, measure in the basis rotated by theta about y
def flip_probability(theta):
    Ry = np.array([[math.cos(theta / 2), -math.sin(theta / 2)], [math.sin(theta / 2), math.cos(theta / 2)]])
    basis = Ry                                  # columns are the reader's basis vectors
    P = np.abs(basis.T @ np.eye(2)) ** 2        # P[outcome, prepared]
    return float(P[1, 0]), float(P[0, 1])


worst = 0.0
for th in np.linspace(0, math.pi, 181):
    p01, p10 = flip_probability(th)
    worst = max(worst, abs(p01 - math.sin(th / 2) ** 2), abs(p10 - math.sin(th / 2) ** 2))
res['max_error_flip_probability'] = worst

# 2. contraction coefficients of one round, sampled
rows = []
for deg in (5, 10, 20, 30, 45, 60, 80):
    th = math.radians(deg)
    p = math.sin(th / 2) ** 2
    W = bsc(p)
    # Dobrushin: exact, the two rows
    eta_tv = tv(W[0], W[1])
    # TV ratio over random pairs of inputs
    a = rng.random(20000); b = rng.random(20000)
    pa = np.stack([a, 1 - a], 1); pb = np.stack([b, 1 - b], 1)
    qa, qb = pa @ W, pb @ W
    tv_ratio = (np.abs(qa - qb).sum(1) / np.abs(pa - pb).sum(1))
    # KL ratio: the supremum is approached for nearby inputs around the uniform distribution
    kl_ratio_random = max(kl(qa[i], qb[i]) / kl(pa[i], pb[i]) for i in range(2000))
    eps = 1e-4
    kl_ratio_local = kl(np.array([0.5 + eps, 0.5 - eps]) @ W, np.array([0.5, 0.5]) @ W) / \
        kl([0.5 + eps, 0.5 - eps], [0.5, 0.5])
    rows.append({'theta_deg': deg, 'p': p, 'eta_TV_exact': eta_tv, 'abs_cos': abs(math.cos(th)),
                 'TV_ratio_min': float(tv_ratio.min()), 'TV_ratio_max': float(tv_ratio.max()),
                 'eta_KL_claim': math.cos(th) ** 2, 'KL_ratio_random_max': kl_ratio_random,
                 'KL_ratio_local': kl_ratio_local})
res['one_round'] = rows

# 3. r rounds: powers of the channel
comp = []
for deg in (10, 30, 60):
    th = math.radians(deg)
    W = bsc(math.sin(th / 2) ** 2)
    Wr = np.eye(2)
    for r in range(1, 65):
        Wr = Wr @ W
        if r in (1, 2, 4, 8, 16, 32, 64):
            comp.append({'theta_deg': deg, 'r': r, 'eta_TV': tv(Wr[0], Wr[1]),
                         'abs_cos_pow_r': abs(math.cos(th)) ** r,
                         'info_bits_uniform_input': 1 + sum(x * math.log2(x) for x in Wr[0] if x > 0)})
res['r_rounds'] = comp


# 4. rounds budget: after how many rounds is the KL coefficient below a level
def rounds_to(level, th):
    return math.log(level) / (2 * math.log(abs(math.cos(th))))


res['rounds_budget'] = [{'theta_deg': d, 'eta_KL': math.cos(math.radians(d)) ** 2,
                         'rounds_to_half': rounds_to(0.5, math.radians(d)),
                         'rounds_to_1_percent': rounds_to(0.01, math.radians(d))}
                        for d in (1, 2, 5, 10, 20, 30, 45, 60)]


# 5. a screen of N independent bits: Dobrushin coefficient of the product channel
def dobrushin_product(p, N):
    # worst pair of inputs is a pair of complementary strings; the output Hamming weights relative
    # to one of them are Binomial(N, p) and Binomial(N, 1 - p), a sufficient statistic
    k = np.arange(N + 1)
    logc = np.array([math.lgamma(N + 1) - math.lgamma(i + 1) - math.lgamma(N - i + 1) for i in k])
    a = np.exp(logc + k * math.log(p) + (N - k) * math.log(1 - p))
    b = a[::-1]
    return tv(a, b)


def dobrushin_bruteforce(p, N):
    W = bsc(p)
    WN = np.array([[1.0]])
    for _ in range(N):
        WN = np.kron(WN, W)
    return max(tv(WN[i], WN[j]) for i in range(2 ** N) for j in range(i))


wide = []
th = math.radians(60)
p = math.sin(th / 2) ** 2
for N in (1, 2, 4, 6, 8, 16, 64, 256):
    row = {'N': N, 'eta_TV_product': dobrushin_product(p, N), 'eta_KL_product': math.cos(th) ** 2}
    if N <= 6:
        row['eta_TV_bruteforce'] = dobrushin_bruteforce(p, N)
    wide.append(row)
res['wide_screen_theta_60deg'] = wide

# 6. the two-channel comparison of Tier 4 in this model
gap = []
for din, dout in ((5, 30), (10, 30), (10, 60), (20, 20), (30, 10)):
    ci, co = math.cos(math.radians(din)) ** 2, math.cos(math.radians(dout)) ** 2
    gap.append({'theta_in_deg': din, 'theta_out_deg': dout, 'eta_in': ci, 'eta_out': co,
                'gap_certified_one_round': ci > co,
                'ratio_out_over_in_after_8_rounds': (co / ci) ** 8})
res['two_channel_gap'] = gap

here = os.path.dirname(os.path.abspath(__file__))
json.dump(res, open(os.path.join(here, 'results', 't4a_frame_mismatch.json'), 'w'), indent=1)

print(f'max error of p = sin^2(theta/2) over 181 angles: {worst:.1e}')
print('theta  eta_TV   |cos|    TV ratio [min, max]   cos^2    KL ratio: random max, local')
for r in rows:
    print(f'{r["theta_deg"]:>4d}  {r["eta_TV_exact"]:.5f}  {r["abs_cos"]:.5f}  [{r["TV_ratio_min"]:.5f}, '
          f'{r["TV_ratio_max"]:.5f}]   {r["eta_KL_claim"]:.5f}  {r["KL_ratio_random_max"]:.5f}  '
          f'{r["KL_ratio_local"]:.5f}')
print('r rounds (theta, r, eta_TV, |cos|^r, bits kept for a uniform input):')
for x in comp:
    print(f'  {x["theta_deg"]:>3d} {x["r"]:>3d}  {x["eta_TV"]:.6f}  {x["abs_cos_pow_r"]:.6f}  '
          f'{x["info_bits_uniform_input"]:.6f}')
print('rounds budget (theta, eta_KL, rounds to 1/2, rounds to 1 percent):')
for x in res['rounds_budget']:
    print(f'  {x["theta_deg"]:>3d}  {x["eta_KL"]:.5f}  {x["rounds_to_half"]:8.1f}  {x["rounds_to_1_percent"]:8.1f}')
print('wide screen at theta = 60 deg (N, eta_TV of the product, brute force, eta_KL):')
for x in wide:
    print(f'  {x["N"]:>4d}  {x["eta_TV_product"]:.6f}  {x.get("eta_TV_bruteforce", float("nan")):.6f}  '
          f'{x["eta_KL_product"]:.4f}')
print('two-channel gap:')
for x in gap:
    print(f'  in {x["theta_in_deg"]:>2d}, out {x["theta_out_deg"]:>2d}: eta_in {x["eta_in"]:.4f}, eta_out '
          f'{x["eta_out"]:.4f}, certified {x["gap_certified_one_round"]}, ratio after 8 rounds '
          f'{x["ratio_out_over_in_after_8_rounds"]:.3g}')
