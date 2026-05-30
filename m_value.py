"""Threshold values used by verify_q0.py and find_terminal_seeds.py.

This module intentionally separates the empirical bisection diagnostic from the
rigorous theorem threshold.  The empirical value is produced by the constraint
checker at the paper's working parameters (eta, nu) = (0.0545, 0.2859), using
the script's bisection for the characteristic-function threshold.  The theorem
itself uses the rigorous closed-form certificate Y_* <= 9.10e31 and states
M < 1.78e32.
The empirical chain reads:
- the bisected Y_43* is ~6.27e31;
- Y47 ~3.61e31 after the c43* safety adjustment;
- Y_major ~7.44e31;
- Y_minor ~7.59e31 is the empirical binding term;
- the integral digit-length specialization uses
  M = ceil((9/2)*Y_required/log 10 + 9/2).

For the seed primes in verify_q0.py (q ~ 10^{135}), q >= M is satisfied
with a ~10^{100} margin.

To re-verify, run the constraint checker at the paper's working parameters
(the defaults):

    .venv/bin/python prime_digit_sums_constraint_checker_with_search.py

and compare with M_EMPIRICAL_LITERAL below.  test_m_value.py automates this
check.
"""

# M = ceil((9/2) * Y_required / log 10 + 9/2), derived in the paper's
# conversion from the log-scale threshold Y_required to the integer digit-sum
# threshold.
M_EMPIRICAL_LITERAL = "148288679844722033445044981196451.0e0"

# Float form for coarse diagnostics only; exact threshold comparisons must use
# M_EMPIRICAL_INT or M_THEOREM_HEADLINE_INT.
M_EMPIRICAL_FLOAT = float(M_EMPIRICAL_LITERAL)

# Exact-integer form for bignum threshold comparisons q >= M where q is an
# exact integer.  This is the integer reading of M_EMPIRICAL_LITERAL (the
# digits before the decimal point); test_m_value.py guards against drift.
M_EMPIRICAL_INT = 148288679844722033445044981196451
assert M_EMPIRICAL_INT == int(M_EMPIRICAL_LITERAL.split(".")[0])

# Rigorous closed-form paper certificate from Y_* <= 9.10e31.  This is not
# used by the seed search, but documents the theorem-level bound alongside the
# sharper empirical literal above.
M_RIGOROUS_CERTIFIED_INT = 177843590339381623423137292296355
M_THEOREM_HEADLINE_INT = 178000000000000000000000000000000
assert M_EMPIRICAL_INT < M_RIGOROUS_CERTIFIED_INT < M_THEOREM_HEADLINE_INT
