"""T2. Record bound and rounds budget on the three-qubit loop of arXiv:2608.13764, built here from
the circuit itself (an independent implementation: nothing is taken from the closed forms except to
compare with them).

 1. Reproduction: Bloch data, Q of Theorem 9 and the reference-point values of the paper.
 2. Record bound on random loops: mu <= F/2 and (2 mu)^2 + T^2 <= 1, with F and T computed from the
    states of the discarded registers.
 3. Rounds: the entanglement-breaking index N on the 81 x 81 slice of arXiv:2609.27908 Sec. 6, and
    on random loops; N (1 - rho(A)).
 4. EXPLORATORY: the largest index found among loops whose discarded registers hold a record of the
    message with trace distance at least T (per round).
ASCII, no em dashes.
"""
import json
import math
import os

import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1.0, -1.0]).astype(complex)
PAULI = [X, Y, Z]
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)


def Ry(a):
    return np.array([[math.cos(a / 2), -math.sin(a / 2)], [math.sin(a / 2), math.cos(a / 2)]], complex)


def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


P0, P1 = np.diag([1.0, 0]).astype(complex), np.diag([0, 1.0]).astype(complex)
# SWAP of M and L with F in between (register order M, F, L)
SWAP_ML = np.zeros((8, 8), complex)
for m in range(2):
    for f in range(2):
        for l in range(2):
            SWAP_ML[(l << 2) | (f << 1) | m, (m << 2) | (f << 1) | l] = 1


def unitary(theta, kappa, phi, beta):
    UW = kron(P0, Ry(math.pi - 2 * theta), I2) + kron(P1, Ry(2 * theta), I2)
    UK = math.cos(kappa / 2) * np.eye(8) - 1j * math.sin(kappa / 2) * kron(I2, Z, Y)
    Ufb = math.cos(phi) * np.eye(8) - 1j * math.sin(phi) * SWAP_ML
    return kron(Ry(beta), I2, I2) @ Ufb @ UK @ UW


def channel_data(U):
    """Bloch pair (A, c) of the message channel and the states of the discarded registers."""
    V = U[:, [0, 4]]                    # isometry: |m> -> U |m, 0, 0>, columns m = 0, 1
    T4 = V.reshape(2, 4, 2)             # indices: M', (F L), m

    def phi_of(rho):
        return np.einsum('aek,kl,bel->ab', T4, rho, T4.conj())

    A = np.zeros((3, 3)); c = np.zeros(3)
    out_I = phi_of(I2 / 2)
    for i in range(3):
        c[i] = np.real(np.trace(PAULI[i] @ out_I))
        for j in range(3):
            A[i, j] = 0.5 * np.real(np.trace(PAULI[i] @ phi_of(PAULI[j])))
    env = [np.einsum('ae,af->ef', T4[:, :, m], T4[:, :, m].conj()) for m in (0, 1)]
    off = phi_of(np.array([[0, 0], [1, 0]], complex))        # Phi(|1><0|)
    return A, c, env, off


def choi_pt(A, c):
    """Partial transpose (on the reference) of the normalised Choi state, from the Bloch pair."""
    J = np.kron(I2, I2 + sum(c[i] * PAULI[i] for i in range(3)))
    for j in range(3):
        J = J + np.kron(PAULI[j], sum(A[i, j] * PAULI[i] for i in range(3)))
    return J / 4


def mu_of(A, c):
    return -float(np.linalg.eigvalsh(choi_pt(A, c))[0])


def index(A, c, cap=400, tol=1e-12):
    An, cn = A.copy(), c.copy()
    for n in range(1, cap + 1):
        if mu_of(An, cn) <= tol:
            return n
        cn = A @ cn + c
        An = A @ An
    return cap + 1                      # not entanglement breaking up to the cap


def psd_sqrt(R):
    w, v = np.linalg.eigh(R)
    return (v * np.sqrt(np.clip(w, 0, None))) @ v.conj().T


def fidelity(R0, R1):
    return float(np.linalg.svd(psd_sqrt(R0) @ psd_sqrt(R1), compute_uv=False).sum())


def trace_distance(R0, R1):
    return 0.5 * float(np.abs(np.linalg.eigvalsh(R0 - R1)).sum())


def closed_form(theta, kappa, phi, beta):
    s, nu = math.sin(2 * theta), math.sin(phi) * math.cos(phi) * math.cos(kappa)
    wv = math.sin(kappa) * math.sin(phi) ** 2 * math.cos(2 * theta)
    cp2 = math.cos(phi) ** 2
    At = np.array([[s * cp2, -s * nu, -wv], [s * nu, s * cp2, 0], [0, 0, cp2]])
    ct = np.array([0, math.sin(kappa) * math.sin(phi) * math.cos(phi) * math.cos(2 * theta),
                   math.sin(phi) ** 2 * math.cos(kappa)])
    Rb = np.array([[math.cos(beta), 0, math.sin(beta)], [0, 1, 0], [-math.sin(beta), 0, math.cos(beta)]])
    return Rb @ At, Rb @ ct


def Qform(theta, kappa, phi):
    C, S, u = math.cos(2 * theta) ** 2, math.sin(phi) ** 2, math.sin(kappa / 2) ** 2
    q = u * (1 - u)
    return ((1 - C) * (1 - S) * (1 - 4 * q * S) * (1 - u * S) * (1 - (1 - u) * S)
            - q * S ** 2 * (1 - u * S - (1 - u) * C) * (1 - (1 - u) * S - u * C))


res = {}
rng = np.random.default_rng(20261114)
pi = math.pi

# ---------------------------------------------------------------- 1. reproduction
ref = (0.16345853 * pi, 0.43230980 * pi, 0.20061939 * pi, 0.23903823 * pi)
A, c, env, off = channel_data(unitary(*ref))
vstar = np.linalg.solve(np.eye(3) - A, c)
res['reference_point'] = {
    'norm_A': float(np.linalg.norm(A, 2)), 'paper_norm_A': 0.71331,
    'v_star': vstar.tolist(), 'paper_v_star': [0.09282, 0.56347, 0.12724],
    'mu': mu_of(A, c), 'paper_mu': 0.20528,
    'Q_from_det': float(-16 * np.linalg.det(choi_pt(A, c)).real), 'Q_formula': Qform(*ref[:3]), 'paper_Q': 0.20381,
    'index': index(A, c), 'paper_index': 3}
dA = dc = dQ = dF = 0.0
for _ in range(2000):
    p = (rng.uniform(0, pi / 2), rng.uniform(0, pi), rng.uniform(0, pi / 2), rng.uniform(0, 2 * pi))
    A, c, env, off = channel_data(unitary(*p))
    Ac, cc = closed_form(*p)
    dA, dc = max(dA, np.abs(A - Ac).max()), max(dc, np.abs(c - cc).max())
    dQ = max(dQ, abs(-16 * np.linalg.det(choi_pt(A, c)).real - Qform(*p[:3])))
    gam = abs(math.sin(2 * p[0]) * math.cos(p[2])) * math.sqrt(1 - math.sin(p[2]) ** 2 * math.sin(p[1]) ** 2)
    dF = max(dF, abs(fidelity(env[0], env[1]) - gam))
res['closed_forms_vs_circuit_2000_points'] = {'max_abs_dev_A': float(dA), 'max_abs_dev_c': float(dc),
                                              'max_abs_dev_Q': float(dQ), 'max_abs_dev_F_vs_gamma': float(dF)}

# ---------------------------------------------------------------- 2 and 4. random loops
Npts = 20000
rec = np.zeros((Npts, 7))
for k in range(Npts):
    p = (rng.uniform(0, pi / 2), rng.uniform(0, pi), rng.uniform(0, pi / 2), rng.uniform(0, 2 * pi))
    A, c, env, off = channel_data(unitary(*p))
    mu = mu_of(A, c)
    F, Tr = fidelity(env[0], env[1]), trace_distance(env[0], env[1])
    rho = float(np.abs(np.linalg.eigvals(A)).max())
    rec[k] = (mu, F, Tr, index(A, c), rho, float(np.linalg.norm(A, 2)), p[2])
mu, F, Tr, Nidx, rho = rec[:, 0], rec[:, 1], rec[:, 2], rec[:, 3], rec[:, 4]
mup = np.clip(mu, 0, None)
res['record_bound_20000_random_loops'] = {
    'max_of_mu_minus_half_F': float((mu - F / 2).max()),
    'max_of_(2mu)^2_plus_T^2': float(((2 * mup) ** 2 + Tr ** 2).max()),
    'violations_mu_le_half_F': int((mu > F / 2 + 1e-10).sum()),
    'violations_duality': int(((2 * mup) ** 2 + Tr ** 2 > 1 + 1e-10).sum()),
    'share_not_EB_at_one_round': float((mu > 1e-12).mean())}
ok = (Nidx > 1) & (Nidx <= 200)
prod = Nidx[ok] * (1 - rho[ok])
res['one_clock_random_loops'] = {'n': int(ok.sum()), 'min': float(prod.min()), 'max': float(prod.max()),
                                 'median': float(np.median(prod)),
                                 'share_in_0.60_1.50': float(((prod >= 0.6) & (prod <= 1.5)).mean())}
front = []
for lo in (0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999):
    m = Tr >= lo
    if m.sum() == 0:
        continue
    nn = Nidx[m]
    front.append({'T_at_least': lo, 'loops': int(m.sum()), 'bound_on_mu': 0.5 * math.sqrt(1 - lo ** 2),
                  'max_mu_found': float(mup[m].max()), 'share_EB_at_one_round': float((nn == 1).mean()),
                  'median_index': float(np.median(nn)), 'p99_index': float(np.percentile(nn, 99)),
                  'max_index': int(nn.max()), 'reached_cap_400': int((nn > 400).sum())})
res['frontier_index_vs_record'] = front

# weak-feedback family: how the budget scales when the record per round is made strong
fam = []
for phi in (0.02, 0.05, 0.1, 0.2):
    for T_target in (0.5, 0.9, 0.99, 0.999):
        best = None
        for th in np.linspace(0.001, pi / 4, 600):
            A, c, env, off = channel_data(unitary(th, pi / 2, phi * pi, 0.0))
            Tr1 = trace_distance(env[0], env[1])
            if Tr1 >= T_target:
                n1 = index(A, c, cap=5000)
                if best is None or n1 > best[0]:
                    best = (n1, th, Tr1, mu_of(A, c))
        if best:
            fam.append({'phi_over_pi': phi, 'T_at_least': T_target, 'max_index': int(best[0]),
                        'theta_over_pi': best[1] / pi, 'T': best[2], 'mu_one_round': best[3]})
res['weak_feedback_family_kappa_half_pi_beta_0'] = fam

# ---------------------------------------------------------------- 3. the 81 x 81 slice of Paper II
k0, b0 = ref[1], ref[3]
grid = np.linspace(0, pi / 2, 81)
Ns, prods = [], []
for th in grid:
    for ph in grid:
        A, c, env, off = channel_data(unitary(th, k0, ph, b0))
        n1 = index(A, c, cap=200)
        Ns.append(n1)
        if 1 < n1 <= 200:
            prods.append(n1 * (1 - float(np.abs(np.linalg.eigvals(A)).max())))
Ns, prods = np.array(Ns), np.array(prods)
res['paperII_slice_81x81'] = {
    'share_index_1': float((Ns == 1).mean()), 'paper_share_index_1': 0.388,
    'distinct_values_up_to_200': int(len(np.unique(Ns[Ns <= 200]))), 'paper_distinct_values': 93,
    'reach_cutoff_200': int((Ns > 200).sum()), 'paper_reach_cutoff': 18,
    'n_with_1_lt_N_le_200': int(len(prods)), 'paper_n': 3997,
    'product_min': float(prods.min()), 'product_max': float(prods.max()), 'product_median': float(np.median(prods)),
    'paper_product': [0.60, 1.50, 1.047]}

here = os.path.dirname(os.path.abspath(__file__))
json.dump(res, open(os.path.join(here, 'results', 't2_record_rounds.json'), 'w'), indent=1)
print(json.dumps(res, indent=1))
