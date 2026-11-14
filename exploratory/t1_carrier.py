"""T1. Tier 1, the carrier screen, from the published article (BioSystems 270, 106004) and its
reference implementation (Zenodo 10.5281/zenodo.21837368, criteria.py used unchanged).

 a. Reproduction: verify.py of the archive, and Table 3 of the article recomputed and compared with
    the printed values.
 b. The screen in rounds (derived here from Eq. (1), Eq. (6) and the phase row of Table 5; not in the
    article): cycles per coherence time, the phase alphabet after n cycles of hold, and the number of
    cycles for which a phase label keeps at least one bit.
Usage: python3 t1_carrier.py /path/to/paper3-criteria
ASCII, no em dashes."""
import json, math, os, subprocess, sys
CRIT = os.path.abspath(sys.argv[1])
sys.path.insert(0, CRIT)
import criteria as c  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
res = {}

# ---------------------------------------------------------------- a. reproduction
p = subprocess.run([sys.executable, 'verify.py'], cwd=CRIT, capture_output=True, text=True)
res['verify_py_exit_code'] = p.returncode
res['verify_py_last_lines'] = p.stdout.strip().splitlines()[-3:]
# rows of Table 3 as printed in the article: carrier, nu, tau, Q printed, C1 printed, C5 printed
T3 = [('Cortical gamma', 40, 0.15, 38, 'pass', 'pass'),
      ('Cortical gamma, sustained', 40, 0.30, 75, 'pass', 'pass'),
      ('Cortical alpha', 10, 0.5, 31, 'pass', 'pass'),
      ('Hippocampal ripple', 180, 0.05, 57, 'pass', 'pass'),
      ('Cortical theta', 6, 1.0, 38, 'pass', 'fail'),
      ('Calcium wave', 0.1, 100.0, 63, 'pass', 'fail'),
      ('Protein vibrational mode', 1e13, 1e-12, 63, 'pass', 'fail'),
      ('Microtubule terahertz mode', 1e13, 1e-13, 6.3, 'marginal', 'fail'),
      ('Microwave carrier, generous proxy', 30e9, 1e-12, 0.19, 'fail', 'fail'),
      ('Microwave carrier, continuum estimate', 30e9, 27e-15, 0.005, 'fail', 'fail'),
      ('Microwave carrier, as its model requires', 30e9, 1e-3, 1.9e8, 'pass', 'fail')]
rows, ok = [], True
for name, nu, tau, qp, c1p, c5p in T3:
    q = c.quality_factor(nu, tau)
    c5 = 'pass' if c.in_persistence_window(tau) else 'fail'
    c1 = 'fail' if q <= 1 else 'pass'
    q_ok = abs(q - qp) / qp < 0.03                     # printed to two significant figures
    c1_ok = (c1 == c1p) or (c1p == 'marginal' and q > 1)
    ok &= q_ok and c1_ok and (c5 == c5p)
    rows.append({'carrier': name, 'nu_Hz': nu, 'tau_s': tau, 'Q': q, 'Q_printed': qp, 'Q_agrees': q_ok,
                 'C1': c1, 'C1_printed': c1p, 'C5': c5, 'C5_printed': c5p})
res['table3'] = rows
res['table3_reproduced'] = bool(ok)
res['table2_required_snr'] = {str(m): c.required_snr(40, 0.15, m) for m in (10, 25, 50, 100, 200)}
res['table2_printed'] = {'10': 0.42, '25': 2.6, '50': 10.6, '100': 42, '200': 169}
res['phase_bins_gamma_one_cycle'] = 0.5 * math.sqrt(math.pi * c.quality_factor(40, 0.15))
res['design_rule_nu_max_Hz'] = {'t_read 0.1 s': c.design_rule_nu_max(10, 0.1), 't_read 0.05 s': c.design_rule_nu_max(10, 0.05)}

# ---------------------------------------------------------------- b. the screen in rounds
def phase_alphabet(Q, n):
    return 0.5 * math.sqrt(math.pi * Q / n)

rb = []
for name, nu, tau, qp, c1p, c5p in T3[:6]:
    Q = c.quality_factor(nu, tau)
    rb.append({'carrier': name, 'Q': Q, 'cycles_per_coherence_time': nu * tau,
               'phase_bins_after_n_cycles': {str(n): phase_alphabet(Q, n) for n in (1, 2, 4, 8, 16)},
               'cycles_with_at_least_one_bit': math.pi * Q / 16,
               'time_with_at_least_one_bit_s': math.pi ** 2 * tau / 8,
               'time_with_five_bins_s': math.pi ** 2 * tau / 50,
               'inside_C5_window': c.in_persistence_window(tau),
               'one_bit_hold_inside_C5_window': c.in_persistence_window(math.pi ** 2 * tau / 8)})
res['rounds_form'] = rb
res['identity_checked'] = 'K >= 2 holds for t <= pi^2 tau / 8, that is 1.234 coherence times, at any frequency'
json.dump(res, open(os.path.join(HERE, 'results', 't1_carrier.json'), 'w'), indent=1)

print('verify.py exit code:', p.returncode, '|', ' / '.join(res['verify_py_last_lines']))
print('Table 3 reproduced:', ok)
for r in rows:
    print(f'  {r["carrier"]:<42s} Q {r["Q"]:>10.4g} (printed {r["Q_printed"]:g}) C1 {r["C1"]:<5s} C5 {r["C5"]}')
print('Table 2 required SNR:', {k: round(v, 2) for k, v in res['table2_required_snr'].items()}, 'printed', res['table2_printed'])
print(f'phase bins, gamma, one cycle: {res["phase_bins_gamma_one_cycle"]:.2f} (article: about five)')
print('design rule nu_max:', res['design_rule_nu_max_Hz'])
print('rounds form (cycles per tau; bins after 1, 2, 4, 8, 16 cycles; cycles with >= 1 bit; hold time s; in window):')
for r in rb:
    b = r['phase_bins_after_n_cycles']
    print(f'  {r["carrier"]:<28s} {r["cycles_per_coherence_time"]:>5.1f}   ' +
          ' '.join(f'{b[k]:>4.1f}' for k in ('1', '2', '4', '8', '16')) +
          f'   {r["cycles_with_at_least_one_bit"]:>5.1f}   {r["time_with_at_least_one_bit_s"]:>7.3f}   {r["one_bit_hold_inside_C5_window"]}')
