"""T0. Tier 0 arithmetic: what the universal floors do and do not constrain.

Exact arithmetic from CODATA constants; the biological inputs are order-of-magnitude textbook values
and are labelled as such. No data are read. ASCII, no em dashes.
"""
import json
import math
import os

h = 6.62607015e-34          # J s (exact, SI 2019)
hbar = h / (2 * math.pi)
kB = 1.380649e-23           # J / K (exact)
c = 299792458.0             # m / s (exact)
G = 6.67430e-11             # m^3 kg^-1 s^-2 (CODATA 2018)
eV = 1.602176634e-19
T = 310.0
LN2 = math.log(2)

kT = kB * T
landauer = kT * LN2
out = {'T_K': T, 'kT_J': kT, 'kT_meV': kT / eV * 1e3, 'landauer_J_per_bit': landauer}

# Margolus-Levitin: an orthogonalisation (one distinguishable step) needs t >= h / (4 E)
ml_kT = h / (4 * kT)
ml_landauer = h / (4 * landauer)
out['margolus_levitin_s'] = {'E = kT': ml_kT, 'E = kT ln2': ml_landauer}
out['FGL2021_empirical_time_limit_s'] = 200e-15
out['ratio_200fs_over_ML_at_kT'] = 200e-15 / ml_kT

# Planck units and the holographic capacity ceiling S = A / (4 l_P^2) nats
lP = math.sqrt(hbar * G / c ** 3)
tP = lP / c
out['planck_length_m'] = lP
out['planck_time_s'] = tP


def sphere_area_from_volume(V):
    r = (3 * V / (4 * math.pi)) ** (1 / 3)
    return 4 * math.pi * r * r, r


def holo_bits(A):
    return A / (4 * lP * lP) / LN2


caps = {}
for name, A in (('human brain, 1.3 L sphere', sphere_area_from_volume(1.3e-3)[0]),
                ('neuron soma, r = 10 um', 4 * math.pi * (10e-6) ** 2),
                ('protein site, r = 1.5 nm (FGL 2021 lower limit 1 to 2 nm)', 4 * math.pi * (1.5e-9) ** 2)):
    caps[name] = {'area_m2': A, 'holographic_ceiling_bits': holo_bits(A)}
out['holographic_ceiling'] = caps

# Energy budget: classical bits that can be rewritten per second at the Landauer price
P_brain = 20.0              # W, textbook value for the adult human brain
n_neurons = 8.6e10          # textbook count
n_syn_lo, n_syn_hi = 1e14, 1e15
rate = P_brain / landauer
out['brain_power_W'] = P_brain
out['max_bit_rewrites_per_s_brain'] = rate
out['max_bit_rewrites_per_s_per_neuron'] = rate / n_neurons
dG_ATP = 50e3 / 6.02214076e23      # J per ATP at about 50 kJ/mol in vivo (order of magnitude)
out['ATP_per_s_at_20W'] = P_brain / dG_ATP
out['landauer_units_per_ATP'] = dG_ATP / landauer

rounds = {'200 fs (FGL 2021 time limit)': 200e-15, '1 ps': 1e-12, '1 ms': 1e-3,
          '25 ms (40 Hz gamma cycle)': 25e-3, '100 ms (alpha cycle)': 0.1, '300 ms': 0.3}
tab = {}
for name, dt in rounds.items():
    bits = rate * dt
    tab[name] = {'round_s': dt,
                 'max_screen_bits_brain': bits,
                 'max_screen_bits_per_neuron': bits / n_neurons,
                 'max_bits_per_synapse_lo_hi': [bits / n_syn_hi, bits / n_syn_lo],
                 'ML_ticks_per_round_at_kT': dt / ml_kT,
                 'planck_times_per_round': dt / tP}
out['screen_capacity_by_round'] = tab

# The product bound: (screen bits) x (round frequency) <= P / (kT ln 2)
out['product_bound_bits_per_s'] = rate
out['max_round_frequency_Hz'] = {'screen = 1 bit per neuron': rate / n_neurons,
                                 'screen = 1 bit per synapse (1e14)': rate / n_syn_lo,
                                 'screen = 1 bit per synapse (1e15)': rate / n_syn_hi}

here = os.path.dirname(os.path.abspath(__file__))
json.dump(out, open(os.path.join(here, 'results', 't0_floor.json'), 'w'), indent=1)


def f(x):
    return f'{x:.3g}'


print(f'kT at 310 K            : {f(kT)} J = {kT / eV * 1e3:.2f} meV')
print(f'Landauer kT ln2        : {f(landauer)} J per bit')
print(f'Margolus-Levitin h/4E  : {f(ml_kT)} s at E = kT; {f(ml_landauer)} s at E = kT ln2')
print(f'200 fs / ML(kT)        : {200e-15 / ml_kT:.1f}')
print(f'Planck length, time    : {f(lP)} m, {f(tP)} s')
for k, v in caps.items():
    print(f'holographic ceiling    : {k}: {f(v["holographic_ceiling_bits"])} bits')
print(f'bit rewrites per s     : brain {f(rate)}, per neuron {f(rate / n_neurons)}')
print(f'ATP per s at 20 W      : {f(P_brain / dG_ATP)}; Landauer units per ATP {dG_ATP / landauer:.1f}')
for k, v in tab.items():
    print(f'round {k:<30s}: brain {f(v["max_screen_bits_brain"])} bits, per neuron '
          f'{f(v["max_screen_bits_per_neuron"])}, ML ticks {f(v["ML_ticks_per_round_at_kT"])}, '
          f'Planck times {f(v["planck_times_per_round"])}')
for k, v in out['max_round_frequency_Hz'].items():
    print(f'max round frequency    : {k}: {f(v)} Hz')
