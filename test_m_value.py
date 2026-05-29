"""Regression test: empirical and rigorous M literals in m_value.py match fresh computations.

Run as:
    .venv/bin/python test_m_value.py

This asserts that the empirical literal stored in m_value.py agrees with the
value computed by prime_digit_sums_constraint_checker_with_search.compute_M at
the paper's working parameters (eta, nu) = (0.0545, 0.2859), and that the rigorous
closed-form certificate from Y_* <= 9.10e31 matches its stored exact integer.
Exits 0 on success, non-zero on mismatch.
"""

import sys

import mpmath as mp

from m_value import M_LITERAL, M_RIGOROUS_CERTIFIED_INT, M_THEOREM_HEADLINE_INT
from prime_digit_sums_constraint_checker_with_search import Inputs, compute_M

mp.mp.dps = 90

inp = Inputs()  # defaults to eta=0.0545, nu=0.2859 — paper's working parameters
result = compute_M(inp)
M_computed = result["M"]
M_from_literal = mp.mpf(M_LITERAL)

rel_diff = abs(M_computed - M_from_literal) / M_from_literal
tol = mp.mpf("1e-39")

if rel_diff >= tol:
    print(f"FAIL: M_LITERAL does not match fresh computation.")
    print(f"  M_LITERAL  = {M_LITERAL}")
    print(f"  M_computed = {mp.nstr(M_computed, 40)}")
    print(f"  rel diff   = {mp.nstr(rel_diff, 5)}  (tolerance {mp.nstr(tol, 2)})")
    sys.exit(1)

print(f"OK: M_LITERAL matches fresh computation (rel diff < 1e-39).")
print(f"  M_LITERAL  = {M_LITERAL}")
print(f"  M_computed = {mp.nstr(M_computed, 40)}")

M_rigorous = int(mp.ceil((mp.mpf(9) / (2 * mp.log(10))) * mp.mpf("9.10e31") + mp.mpf(9) / 2))
if M_rigorous != M_RIGOROUS_CERTIFIED_INT:
    print("FAIL: rigorous closed-form M certificate drifted.")
    print(f"  stored   = {M_RIGOROUS_CERTIFIED_INT}")
    print(f"  computed = {M_rigorous}")
    sys.exit(1)
if not (M_rigorous < M_THEOREM_HEADLINE_INT):
    print("FAIL: rigorous closed-form M certificate is not below the theorem headline.")
    sys.exit(1)

print("OK: rigorous closed-form M certificate matches and is below the theorem headline.")
