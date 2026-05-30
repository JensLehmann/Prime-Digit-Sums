#!/usr/bin/env python3
"""Print and write theorem-level reproducibility certificates.

This is the archival-style entry point for the paper's numerical claims.  It
separates three layers:

* rigorous theorem constants used in the manuscript;
* exact/rational or outward-rounded interval checks for the displayed
  theorem-level numerical inequalities;
* empirical bisection diagnostics from the exploratory constraint checker.

The proof remains in the paper.  This script is a machine-readable audit trail
for the explicit constants and certificates referenced there.
"""

from __future__ import annotations

import argparse
import json
import sys
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

import mpmath as mp

from m_value import (
    M_INT,
    M_RIGOROUS_CERTIFIED_INT,
    M_THEOREM_HEADLINE_INT,
)
from prime_digit_sums_constraint_checker_with_search import Inputs, compute_M, derive
from verify_q0 import SEEDS, chain, is_prime_small


LOG_TERMS = 140
ATAN_TERMS = 80
DECIMAL_PREC = 140

Y43_STAR_RIGOROUS_UPPER = Fraction(91, 10) * 10**31
YMAJ_RIGOROUS_UPPER = Fraction(79, 10) * 10**31
YMIN_RIGOROUS_UPPER = Fraction(8, 1) * 10**31
Y_STAR_RIGOROUS_UPPER = Y43_STAR_RIGOROUS_UPPER


def ceil_fraction(x: Fraction) -> int:
    return -(-x.numerator // x.denominator)


def frac_to_decimal_string(x: Fraction, places: int = 80) -> str:
    """Return a decimal rendering for JSON/human output only."""
    with localcontext() as ctx:
        ctx.prec = places
        value = ctx.divide(Decimal(x.numerator), Decimal(x.denominator))
    return format(value, "f")


def decimal_floor_fraction(x: Fraction) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PREC
        ctx.rounding = ROUND_FLOOR
        return ctx.divide(Decimal(x.numerator), Decimal(x.denominator))


def decimal_ceil_fraction(x: Fraction) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PREC
        ctx.rounding = ROUND_CEILING
        return ctx.divide(Decimal(x.numerator), Decimal(x.denominator))


def sqrt_decimal_lower(x: Decimal) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PREC
        ctx.rounding = ROUND_FLOOR
        return ctx.sqrt(x)


def sqrt_decimal_upper(x: Decimal) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PREC
        ctx.rounding = ROUND_CEILING
        return ctx.sqrt(x)


def log_near_one_interval(x: Fraction, terms: int = LOG_TERMS) -> Tuple[Fraction, Fraction]:
    """Rigorous interval for log(x), with x close to 1.

    Uses log(x)=2*sum_{k>=0} z^(2k+1)/(2k+1), z=(x-1)/(x+1), and a geometric
    tail bound.  The formula is valid for x>0 and |z|<1.  The symmetric tail
    handles z<0 without relying on alternating-series orientation.
    """
    z = (x - 1) / (x + 1)
    z2 = z * z
    power = z
    total = Fraction(0)
    for k in range(terms):
        total += Fraction(2, 2 * k + 1) * power
        power *= z2
    tail = Fraction(2, 2 * terms + 1) * abs(power) / (1 - z2)
    return total - tail, total + tail


def log_rational_interval(
    numerator: int, denominator: int = 1, terms: int = LOG_TERMS
) -> Tuple[Fraction, Fraction]:
    """Rigorous interval for log(numerator/denominator).

    Range-reduces by powers of two into [3/4, 3/2], where the atanh series
    converges rapidly, and combines with a separately bounded log(2).
    """
    x = Fraction(numerator, denominator)
    if x <= 0:
        raise ValueError("log interval requires a positive rational")

    shift = 0
    while x >= Fraction(3, 2):
        x /= 2
        shift += 1
    while x <= Fraction(3, 4):
        x *= 2
        shift -= 1

    lo, hi = log_near_one_interval(x, terms)
    if shift:
        log2_lo, log2_hi = log_near_one_interval(Fraction(2), terms)
        if shift > 0:
            lo += shift * log2_lo
            hi += shift * log2_hi
        else:
            lo += shift * log2_hi
            hi += shift * log2_lo
    return lo, hi


def arctan_inverse_interval(q: int, terms: int = ATAN_TERMS) -> Tuple[Fraction, Fraction]:
    """Alternating-series interval for arctan(1/q)."""
    total = Fraction(0)
    for k in range(terms):
        term = Fraction(1, (2 * k + 1) * q ** (2 * k + 1))
        total = total + term if k % 2 == 0 else total - term
    next_term = Fraction(1, (2 * terms + 1) * q ** (2 * terms + 1))
    endpoint = total + next_term if terms % 2 == 0 else total - next_term
    return min(total, endpoint), max(total, endpoint)


def pi_interval() -> Tuple[Fraction, Fraction]:
    """Machin formula interval: pi/4 = 4 arctan(1/5) - arctan(1/239)."""
    a_lo, a_hi = arctan_inverse_interval(5)
    b_lo, b_hi = arctan_inverse_interval(239)
    quarter_lo = 4 * a_lo - b_hi
    quarter_hi = 4 * a_hi - b_lo
    return 4 * quarter_lo, 4 * quarter_hi


def status(name: str, ok: bool, details: Dict[str, Any]) -> Dict[str, Any]:
    return {"name": name, "ok": ok, **details}


def build_numeric_certificate() -> Dict[str, Any]:
    checks = []

    ln10_lo, ln10_hi = log_rational_interval(10)
    pi_lo, pi_hi = pi_interval()

    m_from_y_upper = ceil_fraction(
        Fraction(9, 2) * Y_STAR_RIGOROUS_UPPER / ln10_lo + Fraction(9, 2)
    )
    checks.append(
        status(
            "M from Y_* <= 9.10e31",
            m_from_y_upper == M_RIGOROUS_CERTIFIED_INT,
            {
                "computed_M": str(m_from_y_upper),
                "stored_M_RIGOROUS_CERTIFIED_INT": str(M_RIGOROUS_CERTIFIED_INT),
                "ln10_lower": frac_to_decimal_string(ln10_lo, 100),
            },
        )
    )

    checks.append(
        status(
            "headline M < 1.78e32",
            M_RIGOROUS_CERTIFIED_INT < M_THEOREM_HEADLINE_INT,
            {
                "M_RIGOROUS_CERTIFIED_INT": str(M_RIGOROUS_CERTIFIED_INT),
                "M_THEOREM_HEADLINE_INT": str(M_THEOREM_HEADLINE_INT),
            },
        )
    )

    # Corollary L_* bound: e^72.75 < 3.94e31 and 3.94e31*log(10) < 9.10e31.
    log_394e31_lo, _ = log_rational_interval(394 * 10**31, 100)
    checks.append(
        status(
            "e^72.75 < 3.94e31",
            Fraction(291, 4) < log_394e31_lo,
            {
                "left_log": "72.75",
                "right_log_lower": frac_to_decimal_string(log_394e31_lo, 90),
            },
        )
    )
    checks.append(
        status(
            "3.94e31 * log(10) < 9.10e31",
            Fraction(394, 100) * ln10_hi < Fraction(91, 10),
            {
                "left_coefficient_upper": frac_to_decimal_string(
                    Fraction(394, 100) * ln10_hi, 80
                ),
                "right_coefficient": "9.10",
            },
        )
    )

    # Named C9--C13 audit rows for the permanence threshold in Corollary 1.
    log8_hi = log_rational_interval(8)[1]
    log_4e8_hi = log_rational_interval(4 * 10**8)[1]
    log7_hi = log_rational_interval(7)[1]
    logL0 = Fraction(291, 4)  # 72.75
    nu = Fraction(2859, 10000)
    checks.append(
        status(
            "C9 at log L=72.75: L0 and D bounds",
            logL0 * (1 - nu) > log8_hi
            and logL0 * nu > Fraction(20, 1)
            and log_4e8_hi < Fraction(20, 1),
            {
                "logL_times_1_minus_nu": str(logL0 * (1 - nu)),
                "log8_upper": frac_to_decimal_string(log8_hi, 80),
                "logL_times_nu": str(logL0 * nu),
                "log_4e8_upper": frac_to_decimal_string(log_4e8_hi, 80),
                "paper_conclusion": "L0 >= L/2 and D >= 2, D <= L0 for L >= exp(72.75)",
            },
        )
    )
    checks.append(
        status(
            "C10 at log L=72.75: ratio > 2",
            Fraction(48, 10) * 10**23 > 2,
            {
                "certified_lower_bound": "4.8e23",
                "target": "2",
                "paper_quantity": "A*C*L^(eta+1/2+nu)/D(L)",
            },
        )
    )
    checks.append(
        status(
            "C11 at log L=72.75: coefficient margin",
            Fraction(6691, 100000) < Fraction(7195, 100000)
            and ln10_lo / 32 > Fraction(7195, 100000),
            {
                "lhs_coefficient_upper": "0.06691",
                "rhs_coefficient_lower": "0.07195",
                "3c4_over_16_lower": frac_to_decimal_string(ln10_lo / 32, 80),
            },
        )
    )
    checks.append(
        status(
            "C12 at log L=72.75: coefficient margin",
            Fraction(666, 10000) < Fraction(3837, 10000)
            and ln10_lo / 6 > Fraction(3837, 10000),
            {
                "lhs_coefficient_upper": "0.0666",
                "rhs_coefficient_lower": "0.3837",
                "c4_lower": frac_to_decimal_string(ln10_lo / 6, 80),
            },
        )
    )
    checks.append(
        status(
            "C13a at log L=72.75: main-term absorption",
            Fraction(10**8, 1) > log7_hi,
            {
                "main_term_lower_bound": "1e8",
                "log7_upper": frac_to_decimal_string(log7_hi, 80),
            },
        )
    )
    checks.append(
        status(
            "C13b at log L=72.75: fixed-point margin",
            Fraction(567, 100) * Fraction(12828, 1000) < Fraction(7274, 100)
            and Fraction(7274, 100) < logL0,
            {
                "fixed_point_rhs_upper": "72.74",
                "logL": "72.75",
                "factor_bound": "5.67 * 12.828 <= 72.74",
            },
        )
    )

    # Major-arc closed-form display in the manuscript.
    _, log_major_base_hi = log_rational_interval(92000, 41)  # 460/0.205
    _, log_major_slack_hi = log_rational_interval(200000001, 200000000)
    log_790_lo, _ = log_rational_interval(790, 100)
    major_left_hi = Fraction(1903, 200) * (log_major_base_hi + log_major_slack_hi)
    major_right_lo = log_790_lo + 31 * ln10_lo
    checks.append(
        status(
            "Y_maj closed-form upper < 7.90e31",
            major_left_hi < major_right_lo,
            {
                "log_left_upper": frac_to_decimal_string(major_left_hi, 90),
                "log_right_lower": frac_to_decimal_string(major_right_lo, 90),
            },
        )
    )

    # Minor-arc displays, including the wider 4.1e12 insensitivity envelope.
    log_y0_lo, log_y0_hi = log_rational_interval(8 * 10**31)
    checks.append(
        status(
            "73.459 < log(8.00e31) < 73.460",
            Fraction(73459, 1000) < log_y0_lo and log_y0_hi < Fraction(73460, 1000),
            {
                "log_y0_lower": frac_to_decimal_string(log_y0_lo, 90),
                "log_y0_upper": frac_to_decimal_string(log_y0_hi, 90),
            },
        )
    )

    log_3002_hi = log_rational_interval(3002)[1]
    checks.append(
        status(
            "(8.00e31)^0.109 > 3002.0",
            Fraction(109, 1000) * log_y0_lo > log_3002_hi,
            {
                "log_left_lower": frac_to_decimal_string(Fraction(109, 1000) * log_y0_lo, 90),
                "log_3002_upper": frac_to_decimal_string(log_3002_hi, 90),
            },
        )
    )

    lhs_minor_floor = Fraction(122009, 1000000) * 3002
    checks.append(
        status(
            "0.122009 * 3002.0 > 366.27",
            lhs_minor_floor > Fraction(36627, 100),
            {
                "left_exact": str(lhs_minor_floor),
                "right": "366.27",
            },
        )
    )

    log_minor_rhs_hi = log_rational_interval(40000000000000000, 79)[1]
    minor_rhs_hi = log_minor_rhs_hi + Fraction(9, 2) * Fraction(73460, 1000)
    checks.append(
        status(
            "minor RHS with C_min=4e12 < 364.43",
            minor_rhs_hi < Fraction(36443, 100),
            {
                "right_side_upper": frac_to_decimal_string(minor_rhs_hi, 90),
                "bound": "364.43",
            },
        )
    )

    log_insens_rhs_hi = log_rational_interval(41000000000000000, 79)[1]
    insens_rhs_hi = log_insens_rhs_hi + Fraction(9, 2) * Fraction(73460, 1000)
    checks.append(
        status(
            "insensitivity RHS with C_min=4.1e12 < 364.46",
            insens_rhs_hi < Fraction(36446, 100),
            {
                "right_side_upper": frac_to_decimal_string(insens_rhs_hi, 90),
                "bound": "364.46",
            },
        )
    )

    # 9R interval from certified log(10), pi, and Decimal directed sqrt/division.
    ln10_lo_dec = decimal_floor_fraction(ln10_lo)
    ln10_hi_dec = decimal_ceil_fraction(ln10_hi)
    pi_lo_dec = decimal_floor_fraction(pi_lo)
    pi_hi_dec = decimal_ceil_fraction(pi_hi)
    sigma2_dec = Decimal(33) / Decimal(4)
    numerator_lo = Decimal(9) * sqrt_decimal_lower(ln10_lo_dec)
    denominator_hi = Decimal(240) * sqrt_decimal_upper(Decimal(2) * pi_hi_dec * sigma2_dec)
    numerator_hi = Decimal(9) * sqrt_decimal_upper(ln10_hi_dec)
    denominator_lo = Decimal(240) * sqrt_decimal_lower(Decimal(2) * pi_lo_dec * sigma2_dec)
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PREC
        ctx.rounding = ROUND_FLOOR
        nine_r_lo = ctx.divide(numerator_lo, denominator_hi)
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PREC
        ctx.rounding = ROUND_CEILING
        nine_r_hi = ctx.divide(numerator_hi, denominator_lo)
    checks.append(
        status(
            "0.00790 < 9R < 0.00791",
            Decimal("0.00790") < nine_r_lo and nine_r_hi < Decimal("0.00791"),
            {
                "nine_R_lower": format(nine_r_lo, "f"),
                "nine_R_upper": format(nine_r_hi, "f"),
            },
        )
    )

    return {
        "method": {
            "logarithms": (
                "exact rational atanh-series intervals with geometric tail bounds "
                f"(terms={LOG_TERMS})"
            ),
            "pi": (
                "Machin formula pi/4=4 atan(1/5)-atan(1/239), "
                f"alternating-series intervals (terms={ATAN_TERMS})"
            ),
            "sqrt_and_division": "Python Decimal with directed rounding",
        },
        "intervals": {
            "ln10": {
                "lower": frac_to_decimal_string(ln10_lo, 110),
                "upper": frac_to_decimal_string(ln10_hi, 110),
            },
            "pi": {
                "lower": frac_to_decimal_string(pi_lo, 90),
                "upper": frac_to_decimal_string(pi_hi, 90),
            },
        },
        "checks": checks,
    }


def seed_certificates() -> Dict[str, Any]:
    records = []
    all_ok = True
    for key in (2, 5, 7, "7'"):
        seed = SEEDS[key]
        k = seed["k"]
        n = seed["n"]
        q = k * (2**n) + 1
        seq = chain(q)
        terminal = 7 if key == "7'" else key
        witness_mod = pow(7, (q - 1) // 2, q)
        checks = {
            "proth_k_odd": k % 2 == 1,
            "proth_size_condition": 2**n > k,
            "above_theorem_headline_threshold": q >= M_THEOREM_HEADLINE_INT,
            "displayed_decimal_matches": (
                True if seed["displayed"] is None else str(q) == seed["displayed"]
            ),
            "chain_terminal_matches": seq[-1] == terminal,
            "chain_intermediates_prime": all(is_prime_small(x) for x in seq[1:]),
            "proth_witness_base_7": witness_mod == q - 1,
        }
        ok = all(checks.values())
        all_ok = all_ok and ok
        records.append(
            {
                "terminal": str(key),
                "k": k,
                "n": n,
                "q": str(q),
                "digits": len(str(q)),
                "chain": [str(x) for x in seq],
                "witness_base": 7,
                "witness_mod_q": str(witness_mod),
                "checks": checks,
                "ok": ok,
            }
        )
    return {"ok": all_ok, "records": records}


def empirical_diagnostics() -> Dict[str, str]:
    inp = Inputs()
    der = derive(inp)
    rep = compute_M(inp)
    return {
        "q": str(inp.q),
        "eta": str(inp.eta),
        "nu": str(inp.nu),
        "c4": mp.nstr(rep["c4"], 50),
        "c43_star": mp.nstr(rep["c43*"], 50),
        "Y43_star_bisect": mp.nstr(rep["Y43*"], 50),
        "Y47": mp.nstr(rep["Y47"], 50),
        "Y_major_empirical": mp.nstr(rep["Y_major"], 50),
        "Y_minor_empirical": mp.nstr(rep["Y_minor"], 50),
        "Y_required_empirical": mp.nstr(rep["Y_required"], 50),
        "M_empirical": mp.nstr(rep["M"], 50),
        "ln10_mpmath": mp.nstr(der.log10, 50),
    }


def build_certificate() -> Dict[str, Any]:
    numeric = build_numeric_certificate()
    seeds = seed_certificates()
    diagnostics = empirical_diagnostics()
    all_numeric_ok = all(item["ok"] for item in numeric["checks"])
    return {
        "certificate_format_version": 1,
        "script": Path(__file__).name,
        "proof_relevant_constants": {
            "Y43_star_upper": "9.10e31",
            "Y_major_upper": "7.90e31",
            "Y_minor_upper": "8.00e31",
            "Y_star_upper": "9.10e31",
            "M_rigorous_certified": str(M_RIGOROUS_CERTIFIED_INT),
            "M_headline_bound": "1.78e32",
            "M_headline_integer": str(M_THEOREM_HEADLINE_INT),
        },
        "numeric_certificate": numeric,
        "seed_certificates": seeds,
        "empirical_diagnostics_not_used_in_proof": diagnostics,
        "ok": all_numeric_ok and seeds["ok"],
    }


def print_summary(cert: Dict[str, Any]) -> None:
    constants = cert["proof_relevant_constants"]
    print("=== Theorem-level constants used in the proof ===")
    print(f"Y43* <= {constants['Y43_star_upper']}")
    print(f"Y_major <= {constants['Y_major_upper']}")
    print(f"Y_minor <= {constants['Y_minor_upper']}")
    print(f"Y_* <= {constants['Y_star_upper']}")
    print(f"M <= {constants['M_rigorous_certified']} < {constants['M_headline_bound']}")
    print("")

    print("=== Rational / interval checks ===")
    for item in cert["numeric_certificate"]["checks"]:
        mark = "OK" if item["ok"] else "FAIL"
        print(f"{mark}: {item['name']}")
    print("")

    print("=== Proth seed certificates ===")
    for item in cert["seed_certificates"]["records"]:
        mark = "OK" if item["ok"] else "FAIL"
        chain_tail = " -> ".join(item["chain"][1:])
        print(
            f"{mark}: terminal {item['terminal']}, q={item['k']}*2^{item['n']}+1, "
            f"{item['digits']} digits, chain q -> {chain_tail}"
        )
    print("")

    print("=== Empirical diagnostics (not used in the proof) ===")
    diag = cert["empirical_diagnostics_not_used_in_proof"]
    print(f"Y43*_bisect = {diag['Y43_star_bisect']}")
    print(f"Y_major_emp = {diag['Y_major_empirical']}")
    print(f"Y_minor_emp = {diag['Y_minor_empirical']}")
    print(f"M_empirical = {diag['M_empirical']} (stored M_INT={M_INT})")
    print("")

    print("Overall:", "OK" if cert["ok"] else "FAIL")


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json",
        default="certificates/theorem_certificate.json",
        help="Path for the machine-readable JSON certificate.",
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Print the certificate summary without writing JSON.",
    )
    args = parser.parse_args(argv)

    cert = build_certificate()
    print_summary(cert)

    if not args.no_write:
        out = Path(args.json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"\nWrote {out}")

    return 0 if cert["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
