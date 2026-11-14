"""T2 core (functions only, extracted from t2_record_rounds.py). Record bound and rounds budget on the three-qubit loop of arXiv:2608.13764, built here from
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


