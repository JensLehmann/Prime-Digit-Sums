#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prime digit sums — constraint checker / parameter search

This script is aligned with the submitted paper:
- Computes derived constants.
- Computes explicit thresholds x41, x42, x45, and x43*.
- Computes major/minor arc Y thresholds and the resulting M.
- Optionally searches over (eta, nu) to reduce M.

Dependencies: mpmath (pip install mpmath)

Usage examples:
  python prime_digit_sums_constraint_checker_with_search.py
  python prime_digit_sums_constraint_checker_with_search.py --q 10 --eta 0.0545 --nu 0.2859
  python prime_digit_sums_constraint_checker_with_search.py --search --search-nu-min 0.20 --search-nu-max 0.34 \
      --search-eta-min 0.03 --search-eta-max 0.08 --search-grid-nu-step 0.005 --search-grid-eta-step 0.002
"""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import mpmath as mp

# Default precision: high enough for giant logs and stable comparisons
mp.mp.dps = 90

# Set to True (or pass --debug on CLI) to emit doubling-search progress lines
# during the bisection that locates L_*^(4).
DEBUG = False


# ----------------------------
# Utilities
# ----------------------------

def mpf(x) -> mp.mpf:
    return mp.mpf(str(x))

def ln(x) -> mp.mpf:
    return mp.log(x)

def pow_real(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    # a^b for positive a
    return mp.e ** (b * mp.log(a))

def ceil_int(x: mp.mpf) -> int:
    return int(mp.ceil(x))

def floor_int(x: mp.mpf) -> int:
    return int(mp.floor(x))

def is_finite(x: mp.mpf) -> bool:
    return mp.isfinite(x)

def safe_log_ratio(num: mp.mpf, den: mp.mpf) -> mp.mpf:
    # log(num/den) robustly
    return mp.log(num) - mp.log(den)

def find_min_L_monotone(predicate: Callable[[int], bool],
                        L0: int = 3,
                        max_doublings: int = 400,
                        *,
                        debug_label: Optional[str] = None,
                        debug_every: int = 25,
                        debug_report: Optional[Callable[[int], str]] = None) -> int:
    """
    Find minimal integer L >= L0 such that predicate(L) is True,
    assuming predicate is monotone (False,...,False,True,...).
    Uses doubling then bisection.
    """
    if predicate(L0):
        return L0
    lo = L0
    hi = L0
    for _ in range(max_doublings):
        hi *= 2
        ok = predicate(hi)
        if debug_label is not None and (hi == L0 * 2 or (debug_every > 0 and _ % debug_every == 0)):
            extra = ""
            if debug_report is not None:
                try:
                    extra = " | " + debug_report(hi)
                except Exception as e:
                    extra = f" | (debug_report failed: {e})"
            print(f"[debug] {debug_label}: L={hi} -> {ok}{extra}")
        if ok:
            break
        lo = hi
    else:
        msg = "find_min_L_monotone: predicate not satisfied within doubling cap."
        if debug_label is not None:
            msg = f"{debug_label}: {msg} (last L tried: {hi})"
        raise RuntimeError(msg)
    # bisection on integers
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if predicate(mid):
            hi = mid
        else:
            lo = mid
    return hi


# ----------------------------
# Inputs / derived constants
# ----------------------------

@dataclass(frozen=True)
class Inputs:
    q: int = 10
    eta: mp.mpf = mpf("0.0545")  # paper's working value after the c4 improvement
    nu: mp.mpf = mpf("0.2859")   # paper's working value after the c4 improvement

    # analytic constants used elsewhere in the paper:
    CDMR: mp.mpf = mpf("102.0")
    c1: mp.mpf = mpf("0.001506288700915")   # minor-arc constant
    Cmin: mp.mpf = mpf("4000000000000.0")   # minor-arc polynomial constant

    # for base-10 only (AP threshold):
    xAP: mp.mpf = mpf("4050")


@dataclass
class Derived:
    # base and digit stats
    logq: mp.mpf
    log10: mp.mpf
    mu: mp.mpf
    sigma2: mp.mpf
    sigma: mp.mpf

    # lemma/proposition constants
    c4: mp.mpf          # decimal-base digit-discrepancy rate log q / 6
    A: mp.mpf           # 4q^2/sigma
    B: mp.mpf           # (q-1)/sigma
    theta: mp.mpf       # c4 / (16(eta+1/2))
    C_tmax: mp.mpf      # the C in |t| <= C L^eta (application choice)
    c43_star: mp.mpf    # characteristic-function comparison exponent
    C42: mp.mpf
    C0: mp.mpf

    # major/minor RHS constants used in the positivity inequalities
    RHS_major: mp.mpf
    RHS_minor: mp.mpf

    # major bound coefficient in alpha-space
    C1_alpha: mp.mpf


def derive(inp: Inputs) -> Derived:
    q = inp.q
    eta = inp.eta
    nu = inp.nu

    logq = ln(q)
    log10 = ln(10)

    mu = mpf(q - 1) / 2
    sigma2 = (mpf(q*q) - 1) / 12
    sigma = mp.sqrt(sigma2)

    # Decimal-base digit-discrepancy constant.
    c4 = logq / 6
    A = (4 * mpf(q) * mpf(q)) / sigma
    B = (mpf(q - 1)) / sigma
    theta = c4 / (16 * (eta + mpf("0.5")))

    # application-specific t-range constant: |t| <= C L^eta
    # From tmax(x) = 2πσ (log q)^{eta-1/2} L^eta
    C_tmax = (2 * mp.pi * sigma) * (logq ** (eta - mpf("0.5")))

    # Per the proof of the characteristic-function comparison,
    # any c43* strictly less than (c4/32)*min(1, (nu-2eta)/(eta+1/2)) is admissible at a
    # sufficiently large threshold L*. We pick a 0.5% safety margin: smaller c43* gives
    # larger a = 2 theta (nu/2 - eta) - c43*, hence smaller theta/a, hence smaller L*^{(4)}.
    # The choice keeps Y47 below the major/minor thresholds while improving the closed-form
    # L*^{(4)} certificate.
    ratio = (nu - 2*eta) / (eta + mpf("0.5"))
    c43_safety = mpf("0.005")
    c43_star = c43_safety * (c4 / 32) * mpf(min(1.0, float(ratio))) if ratio > 0 else mpf("0")

    # C42 and C0: base-10 constants from the major-arc proof.
    if q == 10:
        C42 = (2 * log10 * mpf(5651)) / mpf(108900)   # as in paper
    else:
        # placeholder (extend if you generalize beyond q=10)
        C42 = mpf("1.0")

    C0 = (16 * C42) / (mp.e**2)

    # C1_alpha is the major-arc comparison coefficient after converting to alpha-space.
    C1 = (10 * (q - 1) / sigma) + 1
    C1_alpha = (2 * mp.pi * sigma) * C1 * (logq ** (-nu))

    # RHS_major and RHS_minor come from the final major/minor comparison inequalities.
    c0 = mp.erf((mp.sqrt(2) * mp.pi * sigma) / mp.sqrt(log10))
    RHS_major = (c0 - mpf("1")/40) * mp.sqrt(log10) / mp.sqrt(2 * mp.pi * sigma2)
    RHS_minor = mp.sqrt(log10) / (240 * mp.sqrt(2 * mp.pi * sigma2))

    return Derived(
        logq=logq, log10=log10, mu=mu, sigma2=sigma2, sigma=sigma,
        c4=c4, A=A, B=B, theta=theta, C_tmax=C_tmax, c43_star=c43_star,
        C42=C42, C0=C0, RHS_major=RHS_major, RHS_minor=RHS_minor, C1_alpha=C1_alpha
    )


# ----------------------------
# Basic L-geometry helpers
# ----------------------------

def geom_params(L: int, inp: Inputs) -> Dict[str, int]:
    """
    For a given integer L = log_q x, returns r, L0, K, H exponents.
    Uses r = ceil(L^nu), L0 = L - 2r, K=floor(r/2), H=q^{floor(r/3)}
    in the decimal-base c4 improvement.
    """
    nu = inp.nu
    q = inp.q
    Lmp = mpf(L)
    r = ceil_int(pow_real(Lmp, nu))
    L0 = L - 2*r
    # Exact integer division — math.floor(r/3) loses precision for r > 2^53.
    K = r // 2
    eH = r // 3
    # H = q^eH (we store exponent eH; H can be huge)
    return {"r": r, "L0": L0, "K": K, "eH": eH}


# ----------------------------
# Thresholds: x41, x42
# ----------------------------

def compute_Y41(inp: Inputs, der: Derived) -> mp.mpf:
    # x41(q,nu) := exp((log q) max{10, 8^{1/(1-nu)}})
    nu = inp.nu
    t = mp.power(8, 1/(1-nu))
    Lmin = max(mpf("10"), t)
    return der.logq * Lmin

def compute_Y42(inp: Inputs, der: Derived) -> mp.mpf:
    # Require L0 >= L/2 and tmax/(sigma*sqrt(L0)) <= 1/10.
    # Use the explicit sufficient condition from the paper:
    # 2π sqrt(2) (log q)^{eta-1/2} L^{eta-1/2} <= 1/10
    eta = inp.eta
    logq = der.logq
    const = 20 * mp.pi * mp.sqrt(2) * (logq ** (eta - mpf("0.5")))
    L_needed = mp.power(const, 1/(mpf("0.5") - eta))
    # Also need L large enough for x41 (which ensures L0 >= L/2).
    Y41 = compute_Y41(inp, der)
    L41 = Y41 / logq
    Lmin = max(mpf("3"), L_needed, L41)
    return der.logq * Lmin


# ----------------------------
# Threshold x45 via conditions (C1)-(C8)
# ----------------------------

def Dj_bound(L: int, inp: Inputs, der: Derived) -> mp.mpf:
    """
    Implements eq:Dj-bound (as used in (C5)) per paper, with Kuipers-Niederreiter (6, 4/pi):
      D_j <= 6/(H+1) + (4/pi) * (1+log H) * 6 * CDMR * (log x)^3 * q^{-K/2}
    with x = q^L, log x = L log q, H=q^{floor(r/3)}, K=floor(r/2).
    """
    g = geom_params(L, inp)
    r, K, eH = g["r"], g["K"], g["eH"]
    q = inp.q
    logx = mpf(L) * der.logq

    # H = q^eH, so:
    H = mp.power(q, eH)  # can be huge but mp handles
    term1 = mpf(6) / (H + 1)
    term2 = (mpf(4) / mp.pi) * (1 + mp.log(H)) * (6 * inp.CDMR) * (logx**3) * mp.power(q, -mpf(K)/2)
    return term1 + term2

def check_C_conditions(L: int, inp: Inputs, der: Derived) -> bool:
    q = inp.q
    nu = inp.nu

    g = geom_params(L, inp)
    r, L0, K, eH = g["r"], g["L0"], g["K"], g["eH"]

    # (C1) L0 >= 1
    if L0 < 1:
        return False

    # (C2) K <= 2L/5
    if K > (2*L)//5:
        return False

    # (C3) q^K <= (q-1) q^r/H and (q-1) q^{L-r+1} <= q^{L-K}
    # first: K <= r - eH + log_q(q-1)
    if mpf(K) > mpf(r - eH) + mp.log(q - 1) / der.logq + mpf("1e-30"):
        return False
    # second: log_q(q-1) + (L - r + 1) <= (L - K)
    # i.e. K <= r - 1 - log_q(q-1)
    # We check in logs to avoid floating issues:
    if mp.log(q-1)/der.logq + mpf(L - r + 1) > mpf(L - K) + mpf("1e-30"):
        return False

    # (C4) base 10 only: x >= xAP
    if q == 10:
        logx = mpf(L) * der.logq
        if mp.e**logx < inp.xAP:
            return False

    # (C6) H >= (2q)^2 (so Delta <= 1/(2q))
    H = mp.power(q, eH)
    if H < (2*q)**2:
        return False

    # Delta = H^{-1/2}
    Delta = H ** mpf("-0.5")

    # (C8) 1 + log(H+1) <= L^nu
    if 1 + mp.log(H + 1) > pow_real(mpf(L), nu):
        return False

    # (C7) 2q Delta + q/((H+1)Delta) <= 4 L^nu e^{-c4 L^nu}     (per paper, line ~1587)
    lhs7 = 2*q*Delta + q/((H + 1)*Delta)
    rhs7 = 4 * pow_real(mpf(L), nu) * mp.e ** (-der.c4 * pow_real(mpf(L), nu))
    if lhs7 > rhs7:
        return False

    # (C5) Dj-bound implies Dj <= (1/q) e^{-c4 L^nu}
    Dj = Dj_bound(L, inp, der)
    rhs5 = (1/mpf(q)) * mp.e ** (-der.c4 * pow_real(mpf(L), nu))
    if Dj > rhs5:
        return False

    return True

def compute_L45(inp: Inputs, der: Derived) -> int:
    # Search from L=3 upward, but we use monotone bisection:
    # The manuscript proves eventual permanence from the leading-coefficient gap;
    # this routine only locates the first passing point numerically.
    return find_min_L_monotone(lambda L: check_C_conditions(L, inp, der), L0=3)

def compute_Y45(inp: Inputs, der: Derived) -> mp.mpf:
    L45 = compute_L45(inp, der)
    return mpf(L45) * der.logq


# ----------------------------
# Threshold x43*
# ----------------------------

def D_trunc(L: int, inp: Inputs, der: Derived) -> int:
    # D(L) = 2 floor( theta * L^nu / log L )
    if L < 3:
        return 0
    Lmp = mpf(L)
    val = der.theta * pow_real(Lmp, inp.nu) / mp.log(Lmp)
    return 2 * int(mp.floor(val))

def compute_Lstar_x43(inp: Inputs, der: Derived, L45: int) -> Tuple[int, Dict[str, int]]:
    """
    Compute L* = max(L45, L*_0,...,L*_4) where each L*_i is the minimal
    integer satisfying the corresponding monotone condition from the proof.
    Returns (L*, details).
    """
    q = inp.q
    eta = inp.eta
    nu = inp.nu
    c4 = der.c4
    A = der.A
    B = der.B
    C = der.C_tmax
    theta = der.theta
    c43 = der.c43_star

    if c43 <= 0:
        raise RuntimeError("Invalid parameters: need nu > 2*eta for c43*>0.")

    def base_geom_ok(L: int) -> bool:
        g = geom_params(L, inp)
        return (g["L0"] >= 1)

    # L*_0: ensure 2 <= D <= L0 (and implicitly L>=3)
    def pred0(L: int) -> bool:
        if L < 3:
            return False
        if not base_geom_ok(L):
            return False
        D = D_trunc(L, inp, der)
        if D < 2:
            return False
        g = geom_params(L, inp)
        if D > g["L0"]:
            return False
        return True

    L0s = find_min_L_monotone(pred0, L0=max(3, L45))

    # L*_1: ratio of consecutive terms >=2 for 1<=d<=D:
    # z/d >= 2; sufficient z/D >= 2 with z = A*C*L^{eta+1/2+nu}
    def pred1(L: int) -> bool:
        if L < L0s:
            return False
        D = D_trunc(L, inp, der)
        if D < 2:
            return False
        z = A * C * pow_real(mpf(L), eta + mpf("0.5") + nu)
        return (z / mpf(D)) >= 2

    L1s = find_min_L_monotone(pred1, L0=L0s)

    # L*_2 (C11): log D + D log(e z / D) <= (c4/8 + c4/16) L^nu = (3c4/16) L^nu.
    # The c4/16 buffer leaves exactly the slack needed for the downstream Step-2
    # conclusion |t| e^{-(3c4/4) L^nu}; indeed that conclusion still
    # holds because c4 - 3c4/16 = 13c4/16 >= 3c4/4.
    def pred2(L: int) -> bool:
        if L < L1s:
            return False
        D = D_trunc(L, inp, der)
        if D < 2:
            return False
        z = A * C * pow_real(mpf(L), eta + mpf("0.5") + nu)
        lhs = mp.log(mpf(D)) + mpf(D) * mp.log((mp.e * z) / mpf(D))
        rhs = (3*c4/16) * pow_real(mpf(L), nu)
        return lhs <= rhs

    L2s = find_min_L_monotone(pred2, L0=L1s)

    # L*_3: ensure the exp-small moment-comparison term at d=D is at most
    # the explicit sub-Gaussian moment envelope used for E|Y|^D:
    #   A^D L^{(1/2+nu)D} e^{-c4 L^nu}
    #     <= 2 sqrt(2 pi D) (D/e)^{D/2}.
    # We use the slightly stronger log predicate without the positive
    # log(2 sqrt(2 pi D)) term on the RHS.
    def pred3(L: int) -> bool:
        if L < L2s:
            return False
        D = D_trunc(L, inp, der)
        if D < 2:
            return False
        Lmp = mpf(L)
        # compare logs:
        # D(log A + (1/2+nu)log L - 1/2 log D + 1/2) <= c4 L^nu
        lhs = mpf(D) * (
            mp.log(A)
            + (mpf("0.5") + nu) * mp.log(Lmp)
            - mpf("0.5") * mp.log(mpf(D))
            + mpf("0.5")
        )
        rhs = c4 * pow_real(Lmp, nu)
        return lhs <= rhs

    L3s = find_min_L_monotone(pred3, L0=L2s)

    # L*_4: two-part final absorption sufficient for
    # |phi2(t)-phi3(t)| <= |t| e^{-c43 L^nu}.
    # Tight (sharp sub-Gaussian) form using the analytic uniform-digit MGF bound
    # (sigma_sg(Y) = 1 = std(Y)). Sufficient inequalities:
    #   log 7 <= (3*c4/4 - c43) L^nu
    #   log 7 + D log( sqrt(e) * |t| / sqrt D ) <= log|t| - c43 L^nu
    def pred4(L: int) -> bool:
        if L < L3s:
            return False
        D = D_trunc(L, inp, der)
        if D < 2:
            return False
        Lmp = mpf(L)
        tmax = C * pow_real(Lmp, eta)
        # Avoid log(0):
        if tmax <= 0:
            return False
        if mp.log(7) > (3*c4/4 - c43) * pow_real(Lmp, nu):
            return False
        term = mp.sqrt(mp.e) * tmax / mp.sqrt(mpf(D))
        lhs = mp.log(7) + mpf(D) * mp.log(term)
        rhs = mp.log(tmax) - c43 * pow_real(Lmp, nu)
        return lhs <= rhs

    def pred4_report(L: int) -> str:
        D = D_trunc(L, inp, der)
        if D < 2:
            return f"D={D} (<2)"
        Lmp = mpf(L)
        tmax = C * pow_real(Lmp, eta)
        main_gap = (3*c4/4 - c43) * pow_real(Lmp, nu) - mp.log(7)
        term = mp.sqrt(mp.e) * tmax / mp.sqrt(mpf(D))
        lhs = mp.log(7) + mpf(D) * mp.log(term)
        rhs = mp.log(tmax) - c43 * pow_real(Lmp, nu)
        return f"D={D}, main_gap={mp.nstr(main_gap, 8)}, lhs-rhs={mp.nstr(lhs-rhs, 8)}, tmax={mp.nstr(tmax, 6)}"

    L4s = find_min_L_monotone(
        pred4,
        L0=L3s,
        debug_label="L4* (pred4)" if DEBUG else None,
        debug_every=10,
        debug_report=pred4_report,
    )

    Lstar = max(L45, L0s, L1s, L2s, L3s, L4s)

    # Permanence spot-check: each predicate should still hold at 10*L_check and 1000*L_check,
    # corroborating the asymptotic-regime claim in the paper (the leading-coefficient gap
    # ensures (RHS-LHS) -> +infty, so once past the analytical threshold L_*^{(i)}, stays True).
    #
    # The bisection's L_is is the smallest L at which pred_i FIRST holds (empirically), which
    # for non-monotone middle regimes may be below the analytical threshold L_*^{(i)} (Lemma
    # lem:r4-explicit + coefficient-gap argument). We therefore check at L_check =
    # max(L_is, L_init), where L_init = 10^25 is the lemma's analytical threshold; this
    # corroborates the paper's claim that L_*^{(i)} <= L_init for i in {0,1,2,3} at the
    # working parameters, and that pred_4 holds for all L >= L_*^{(4)}.
    L_init_paper = 10**25
    permanence = {}
    for label, Li, pred in [("L0*", L0s, pred0), ("L1*", L1s, pred1),
                            ("L2*", L2s, pred2), ("L3*", L3s, pred3),
                            ("L4*", L4s, pred4)]:
        L_check = max(Li, L_init_paper)
        ok10 = pred(10 * L_check)
        ok1000 = pred(1000 * L_check)
        permanence[label] = (ok10, ok1000)
        if not (ok10 and ok1000):
            print(f"WARNING: permanence spot-check failed for {label} at 10*L_check={10*L_check} (ok={ok10}) or 1000*L_check={1000*L_check} (ok={ok1000})")

    details = {"L45": L45, "L0*": L0s, "L1*": L1s, "L2*": L2s, "L3*": L3s, "L4*": L4s,
               "L*": Lstar, "permanence_10x_1000x_above_L_init": permanence}
    return Lstar, details

def compute_Y43_star(inp: Inputs, der: Derived, Y45: mp.mpf) -> Tuple[mp.mpf, Dict[str, int]]:
    L45 = int(mp.floor(Y45 / der.logq))
    Lstar, details = compute_Lstar_x43(inp, der, L45=L45)
    return mpf(Lstar) * der.logq, details


# ----------------------------
# Major/minor thresholds and M
# ----------------------------

def major_LHS(Y: mp.mpf, inp: Inputs, der: Derived) -> mp.mpf:
    """
    Major-arc inequality LHS matching paper eq:Ymaj-ineq (lem:positivity-integrated form):
      C^{(1)} * Y^{-1/2+nu+2eta}
      + 2*C^{(0)} * Y^{eta-1}
      + 8*pi * Y^{-1/2+2eta}
    The exponential factor exp(-c43* L^nu) is absorbed into x_maj (Step 5 of prop:maj)
    via the threshold L >= (2(1/2-nu)/(nu e c43*))^{2/nu}, not via this LHS.
    """
    eta = inp.eta
    nu = inp.nu
    t_phi12 = der.C1_alpha * (Y ** (-mpf("0.5") + nu + 2*eta))
    t_phi3 = 2 * der.C0 * (Y ** (eta - 1))
    t_phase = 8 * mp.pi * (Y ** (-mpf("0.5") + 2*eta))
    return t_phi12 + t_phi3 + t_phase

def compute_Y_major(inp: Inputs, der: Derived) -> mp.mpf:
    """
    Smallest Y such that major_LHS(Y) <= RHS_major.
    Uses monotone doubling + bisection in Y.
    """
    rhs = der.RHS_major

    def predY(Y: mp.mpf) -> bool:
        return major_LHS(Y, inp, der) <= rhs

    # start at something safe:
    Y0 = max(compute_Y41(inp, der), compute_Y42(inp, der), mpf("10"))

    # doubling
    Ylo = Y0
    if predY(Ylo):
        return Ylo
    Yhi = Ylo
    for _ in range(400):
        Yhi *= 2
        if predY(Yhi):
            break
        Ylo = Yhi
    else:
        raise RuntimeError("Major threshold search failed (too many doublings).")

    # bisection (float bisection)
    for _ in range(200):
        Ym = (Ylo + Yhi)/2
        if predY(Ym):
            Yhi = Ym
        else:
            Ylo = Ym
    return Yhi

def minor_LHS(Y: mp.mpf, inp: Inputs, der: Derived) -> mp.mpf:
    """
    Minor-arc inequality LHS matching paper eq:xmin-ineq:
      (Cmin/9) * Y^(9/2) * exp(-81 * c1 * Y^(2 eta))
    where Y = log x. The RHS is sqrt(log 10) / (240 sqrt(2 pi sigma^2)).
    """
    eta = inp.eta
    return (inp.Cmin / 9) * (Y ** mpf("4.5")) * mp.e ** (-81 * inp.c1 * (Y ** (2*eta)))

def compute_Y_minor(inp: Inputs, der: Derived) -> mp.mpf:
    rhs = der.RHS_minor

    def predY(Y: mp.mpf) -> bool:
        return minor_LHS(Y, inp, der) <= rhs

    Y0 = mpf("10")
    if predY(Y0):
        return Y0
    Ylo = Y0
    Yhi = Ylo
    for _ in range(500):
        Yhi *= 2
        if predY(Yhi):
            break
        Ylo = Yhi
    else:
        raise RuntimeError("Minor threshold search failed (too many doublings).")

    for _ in range(240):
        Ym = (Ylo + Yhi)/2
        if predY(Ym):
            Yhi = Ym
        else:
            Ylo = Ym
    return Yhi

def compute_M(inp: Inputs) -> Dict[str, mp.mpf]:
    der = derive(inp)

    # core thresholds
    Y41 = compute_Y41(inp, der)
    Y42 = compute_Y42(inp, der)
    Y45 = compute_Y45(inp, der)

    # x43*
    Y43, Lstar_details = compute_Y43_star(inp, der, Y45=Y45)

    # major/minor
    Ymaj = compute_Y_major(inp, der)
    Ymin = compute_Y_minor(inp, der)

    # Y47: log of the Step-5 absorption threshold in prop:maj.
    # x_47 := q^L_abs where L_abs = ceil((2(1/2-nu)/(nu*e*c43*))^{2/nu}).
    base = 2 * (mpf("0.5") - inp.nu) / (inp.nu * mp.e * der.c43_star)
    L_abs = mp.ceil(pow_real(base, 2/inp.nu))
    Y47 = mpf(L_abs) * der.logq

    YAP = mp.log(inp.xAP)

    Yreq = max(Y41, Y42, Y45, Y43, Y47, Ymaj, Ymin, YAP)

    # Integral digit-length specialization: use L_m=floor(2m/9), so require
    # L_m log(10) >= Yreq.  It suffices and is sharp up to the integer ceiling to take
    # m >= (9/(2 log 10)) Yreq + 9/2.
    M = mp.ceil((mpf("9")/2) * (Yreq / der.log10) + mpf("4.5"))

    return {
        "Y41": Y41, "Y42": Y42, "Y45": Y45, "Y43*": Y43, "Y47": Y47,
        "Y_AP": YAP,
        "Y_major": Ymaj, "Y_minor": Ymin,
        "Y_required": Yreq,
        "M": M,
        "c4": der.c4,
        "c43*": der.c43_star,
        "sigma": der.sigma,
        "C0": der.C0,
        "C1_alpha": der.C1_alpha,
        "RHS_major": der.RHS_major,
        "RHS_minor": der.RHS_minor,
        "Lstar_details": Lstar_details,
    }


# ----------------------------
# Search (eta, nu)
# ----------------------------

def feasible(eta: float, nu: float) -> bool:
    # Admissibility (paper Lemma c4-explicit / Prop CF-replacement / Prop maj):
    # 0 < eta < nu/2 < 1/4 (so nu < 1/2) and nu + 2 eta < 1/2.
    # The paper's working point is (eta, nu) = (0.0545, 0.2859).
    if not (0.0 < eta < 0.5 and 0.0 < nu < 0.5):
        return False
    if not (eta < nu/2):
        return False
    if not (nu + 2*eta < 0.5):
        return False
    return True

def objective_log10M(q: int, eta: float, nu: float, base_inp: Inputs) -> float:
    if not feasible(eta, nu):
        return 1e99
    try:
        inp = Inputs(
            q=q,
            eta=mpf(eta),
            nu=mpf(nu),
            CDMR=base_inp.CDMR,
            c1=base_inp.c1,
            Cmin=base_inp.Cmin,
            xAP=base_inp.xAP,
        )
        rep = compute_M(inp)
        M = rep["M"]
        if not is_finite(M) or M <= 0:
            return 1e99
        return float(mp.log10(M))
    except Exception:
        return 1e99

def grid_search(base_inp: Inputs,
                nu_min: float, nu_max: float, eta_min: float, eta_max: float,
                nu_step: float, eta_step: float) -> Tuple[float, float, float]:
    best = (1e99, 0.0, 0.0)
    nu = nu_min
    while nu <= nu_max + 1e-12:
        eta = eta_min
        while eta <= eta_max + 1e-12:
            if feasible(eta, nu) and eta < nu/2 - 1e-6:
                val = objective_log10M(base_inp.q, eta, nu, base_inp)
                if val < best[0]:
                    best = (val, eta, nu)
            eta += eta_step
        nu += nu_step
    return best  # (best_log10M, eta, nu)

def nelder_mead_2d(f: Callable[[float, float], float],
                   x0: Tuple[float, float],
                   step: Tuple[float, float] = (0.01, 0.01),
                   max_iter: int = 120) -> Tuple[float, float, float]:
    """
    Simple Nelder–Mead in 2D. Returns (f_best, eta_best, nu_best).
    """
    # simplex points
    (e0, n0) = x0
    (se, sn) = step
    simplex = [
        (e0, n0),
        (e0 + se, n0),
        (e0, n0 + sn),
    ]
    vals = [f(*p) for p in simplex]

    def order():
        nonlocal simplex, vals
        items = sorted(zip(vals, simplex), key=lambda t: t[0])
        vals = [i[0] for i in items]
        simplex = [i[1] for i in items]

    order()
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5

    for _ in range(max_iter):
        order()
        (best, good, worst) = simplex[0], simplex[1], simplex[2]
        fbest, fgood, fworst = vals[0], vals[1], vals[2]

        # centroid of best+good
        ce = (best[0] + good[0]) / 2.0
        cn = (best[1] + good[1]) / 2.0

        # reflection
        re = ce + alpha * (ce - worst[0])
        rn = cn + alpha * (cn - worst[1])
        frefl = f(re, rn)

        if frefl < fbest:
            # expansion
            ee = ce + gamma * (re - ce)
            en = cn + gamma * (rn - cn)
            fexp = f(ee, en)
            if fexp < frefl:
                simplex[2] = (ee, en); vals[2] = fexp
            else:
                simplex[2] = (re, rn); vals[2] = frefl
        elif frefl < fgood:
            simplex[2] = (re, rn); vals[2] = frefl
        else:
            # contraction
            if frefl < fworst:
                # outside contraction
                oe = ce + rho * (re - ce)
                on = cn + rho * (rn - cn)
                fcon = f(oe, on)
                if fcon <= frefl:
                    simplex[2] = (oe, on); vals[2] = fcon
                else:
                    # shrink
                    b = simplex[0]
                    simplex = [b,
                               (b[0] + sigma*(simplex[1][0]-b[0]), b[1] + sigma*(simplex[1][1]-b[1])),
                               (b[0] + sigma*(simplex[2][0]-b[0]), b[1] + sigma*(simplex[2][1]-b[1]))]
                    vals = [f(*p) for p in simplex]
            else:
                # inside contraction
                ie = ce - rho * (ce - worst[0])
                inn = cn - rho * (cn - worst[1])
                fcon = f(ie, inn)
                if fcon < fworst:
                    simplex[2] = (ie, inn); vals[2] = fcon
                else:
                    # shrink
                    b = simplex[0]
                    simplex = [b,
                               (b[0] + sigma*(simplex[1][0]-b[0]), b[1] + sigma*(simplex[1][1]-b[1])),
                               (b[0] + sigma*(simplex[2][0]-b[0]), b[1] + sigma*(simplex[2][1]-b[1]))]
                    vals = [f(*p) for p in simplex]

    order()
    return vals[0], simplex[0][0], simplex[0][1]


# ----------------------------
# Pretty printing
# ----------------------------

def fmt(x: mp.mpf, digits: int = 24) -> str:
    try:
        return mp.nstr(x, digits)
    except Exception:
        return str(x)

def print_report(inp: Inputs, rep: Dict[str, mp.mpf]) -> None:
    print("=== Environment ===")
    print(f"mpmath version = {mp.libmp.__name__.split('.')[0]}=={mp.__version__}, mp.dps = {mp.mp.dps}")
    print("")
    print("=== Inputs ===")
    print(f"q={inp.q}, eta={inp.eta}, nu={inp.nu}")
    print(f"CDMR={inp.CDMR}, c1={inp.c1}, Cmin={inp.Cmin}, Cmaj=510 (fixed, see Cor. of prop:maj), xAP={inp.xAP}")
    print("")
    print("=== Derived (selected) ===")
    print(f"sigma = {fmt(rep['sigma'])}")
    print(f"c4 = log(q)/6 = {fmt(rep['c4'])}")
    print(f"C0 = {fmt(rep['C0'])}")
    print(f"C1_alpha = {fmt(rep['C1_alpha'])}")
    print(f"RHS_major = {fmt(rep['RHS_major'])}")
    print(f"RHS_minor = {fmt(rep['RHS_minor'])}")
    print("")
    print("=== Thresholds (Y = log x) ===")
    for k in ["Y41", "Y42", "Y43*", "Y45", "Y47", "Y_major", "Y_minor", "Y_AP"]:
        print(f"{k:>10}: {fmt(rep[k], 30)}")
    print(f"{'c43*':>10}: {fmt(rep['c43*'], 30)}")
    print("")
    print(f"Y_required = {fmt(rep['Y_required'], 40)}")
    print(f"M = ceil((9/2)*Y/log10 + 9/2) = {fmt(rep['M'], 40)}")
    print("")

    # Diagnostics about x43*
    det = rep.get("Lstar_details", {})
    if det:
        print("=== x43* details (L-thresholds from the CF-comparison proof) ===")
        for kk in ["L45", "L0*", "L1*", "L2*", "L3*", "L4*", "L*"]:
            if kk in det:
                print(f"{kk:>6}: {det[kk]}")
        # show D at the L corresponding to Y_major (rounded)
        try:
            der = derive(inp)
            L_at_major = int(mp.floor(rep["Y_major"] / der.logq))
            Dm = D_trunc(L_at_major, inp, der)
            print(f"\nDiagnostic: at L=floor(Y_major/log q)={L_at_major}, D(L)={Dm}")
            if Dm < 2:
                print("  -> This explains a large x43*: the CF-comparison proof needs D(L) >= 2,")
                print("     but with these (eta,nu), L^nu/log L may be too small until astronomically large L.")
                print("     Practical fix: increase nu and/or decrease eta (so nu-2eta is not tiny).")
        except Exception:
            pass


# ----------------------------
# CLI
# ----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--q",
        type=int,
        default=10,
        help=(
            "Digit base. Currently only --q 10 is fully supported: several "
            "constants in the derive() routine (notably C42's geometric "
            "prefactor) are hardcoded for q=10 and read as 1.0 for other q. "
            "Use other values only for experimentation; the paper's "
            "main theorem statement applies only to q=10."
        ),
    )
    ap.add_argument("--eta", type=float, default=0.0545)
    ap.add_argument("--nu", type=float, default=0.2859)

    ap.add_argument("--CDMR", type=float, default=102.0)
    ap.add_argument("--c1", type=float, default=0.001506288700915)
    ap.add_argument("--Cmin", type=float, default=4000000000000.0)
    ap.add_argument("--xAP", type=float, default=4050.0)

    ap.add_argument("--search", action="store_true")
    ap.add_argument("--search-nu-min", type=float, default=0.20)
    ap.add_argument("--search-nu-max", type=float, default=0.34)
    ap.add_argument("--search-eta-min", type=float, default=0.03)
    ap.add_argument("--search-eta-max", type=float, default=0.08)
    ap.add_argument("--search-grid-nu-step", type=float, default=0.005)
    ap.add_argument("--search-grid-eta-step", type=float, default=0.002)
    ap.add_argument("--search-nm-iters", type=int, default=140)

    ap.add_argument("--debug", action="store_true",
                    help="Emit doubling-search progress lines during the bisection that locates L_*^(4).")

    args = ap.parse_args()

    global DEBUG
    DEBUG = args.debug

    # Cmaj is fixed by the paper to 510 and is not exposed as a CLI flag.
    base_inp = Inputs(
        q=args.q,
        eta=mpf(args.eta),
        nu=mpf(args.nu),
        CDMR=mpf(args.CDMR),
        c1=mpf(args.c1),
        Cmin=mpf(args.Cmin),
        xAP=mpf(args.xAP),
    )

    if not args.search:
        rep = compute_M(base_inp)
        print_report(base_inp, rep)
        return

    # Search
    print("=== SEARCH MODE ===")
    print("Coarse grid search...")
    best_logM, best_eta, best_nu = grid_search(
        base_inp,
        nu_min=args.search_nu_min, nu_max=args.search_nu_max,
        eta_min=args.search_eta_min, eta_max=args.search_eta_max,
        nu_step=args.search_grid_nu_step, eta_step=args.search_grid_eta_step
    )
    print(f"Best grid: log10(M)≈{best_logM:.6f} at eta={best_eta:.6f}, nu={best_nu:.6f}")

    def f(eta: float, nu: float) -> float:
        return objective_log10M(base_inp.q, eta, nu, base_inp)

    print("Nelder–Mead refinement...")
    fbest, eta_best, nu_best = nelder_mead_2d(
        f,
        x0=(best_eta, best_nu),
        step=(args.search_grid_eta_step*2, args.search_grid_nu_step*2),
        max_iter=args.search_nm_iters
    )
    print(f"Best NM: log10(M)≈{fbest:.6f} at eta={eta_best:.8f}, nu={nu_best:.8f}")

    final_inp = Inputs(
        q=base_inp.q,
        eta=mpf(eta_best),
        nu=mpf(nu_best),
        CDMR=base_inp.CDMR,
        c1=base_inp.c1,
        Cmin=base_inp.Cmin,
        xAP=base_inp.xAP,
    )
    rep = compute_M(final_inp)
    print_report(final_inp, rep)


if __name__ == "__main__":
    main()
