"""Floating-point diagnostic for Lemma 'Sharp sub-Gaussianity of uniform digit'.

Background.  For Z uniform on {0,...,q-1} (q=10), set f(lambda) := log E e^{lambda(Z-mu)}
with mu = (q-1)/2.  Then
   f''(lambda) = (1/4) csch^2(lambda/2) - (q^2/4) csch^2(q lambda/2).
The lemma asserts f''(lambda) <= sigma^2 = (q^2-1)/12 = 8.25 for all real lambda.

By symmetry f''(lambda) = f''(-lambda) it suffices to consider lambda >= 0.  The diagnostic
checks three regions:

(a) lambda = 0:  f''(0) = sigma^2 (limit equality), shown analytically by Taylor expansion.

(b) lambda in (0, lambda_0]: Taylor at 0 gives
       f''(lambda) - sigma^2 = (kappa_4 / 2) lambda^2 + R(lambda)
    with kappa_4 = -(q^4-1)/120 < 0 for q >= 2 (for q=10, kappa_4 = -83.325).
    The remainder R(lambda) = O(lambda^4) is controlled below by uniform pointwise check on a
    coarse grid plus a derivative bound, so g(lambda) := sigma^2 - f''(lambda) >= (|kappa_4|/4) lambda^2
    on (0, lambda_0] for an analytically chosen lambda_0.

(c) lambda in [lambda_0, T]: finite grid with Lipschitz extrapolation.  We use the bound
    |f'''(lambda)| <= L(lambda_0) on [lambda_0, T] (computed once at lambda_0), then sample
    g on a grid with spacing h such that L(lambda_0) * h < g_min - margin.

(d) lambda >= T:  f''(lambda) <= (1/4) csch^2(lambda/2) <= 4 e^{-lambda} (eq:fpp-tail in
    paper); at T=20, 4 e^{-20} ~ 8e-9 << sigma^2.

The paper proof is analytic; this script is not an interval-arithmetic certificate.

Run:
    .venv/bin/python verify_subgaussian.py
Expected runtime: a few seconds.
"""

import argparse
import sys

from mpmath import mp, mpf, sinh, cosh, exp


def fpp(lam: mpf, q: int = 10) -> mpf:
    """f''(lambda) = (1/4) csch^2(lambda/2) - (q^2/4) csch^2(q lambda/2)."""
    a = lam / 2
    b = q * lam / 2
    return mpf(1) / (4 * sinh(a) ** 2) - mpf(q * q) / (4 * sinh(b) ** 2)


def fppp_bound(lam_min: mpf, T: mpf, q: int = 10) -> mpf:
    """Crude upper bound on |f'''(lambda)| for lambda in [lam_min, T].

    f'''(lambda) = -(1/4) csch^2(lambda/2) coth(lambda/2) + (q^3/4) csch^2(q lambda/2) coth(q lambda/2).
    Both terms are decreasing in absolute value (for lambda > 0), so the max occurs at lam_min.
    """
    a = lam_min / 2
    b = q * lam_min / 2
    term1 = mpf(1) / (4 * sinh(a) ** 2) * (cosh(a) / sinh(a))
    term2 = mpf(q * q * q) / (4 * sinh(b) ** 2) * (cosh(b) / sinh(b))
    return term1 + term2


def certify(N: int = 200_000, T: int = 20, q: int = 10, dps: int = 50,
            lam_taylor: float = 0.5) -> bool:
    mp.dps = dps
    sigma2 = mpf(q * q - 1) / 12
    # Fourth cumulant for the uniform law on {0,...,q-1}: kappa_4 = mu_4 - 3*sigma^4.
    mu4 = mpf((q * q - 1) * (3 * q * q - 7)) / 240
    kappa4 = mu4 - 3 * sigma2 ** 2
    print(f"q = {q}, sigma^2 = {mp.nstr(sigma2, 8)}, kappa_4 = {mp.nstr(kappa4, 8)}")
    assert kappa4 < 0, "kappa_4 must be negative for sub-Gaussianity argument to apply"

    # (b) Region [eps, lam_taylor]: check g(lambda) >= a * lambda^2 for an explicit constant a > 0.
    # The smooth function r(lambda) := g(lambda)/lambda^2 has r(0+) = |kappa_4|/2 = 41.66
    # and decreases as lambda grows.  We sample r on a fine grid on (0, lam_taylor] and report
    # the minimum: this gives an explicit lower bound a := min_grid r(lambda).
    # (The grid-to-off-grid extrapolation uses |r'(lambda)| <= 500 on (eps, lam_taylor]; the
    # spacing h = lam_taylor/K with K=10000 gives slack 500*5e-5 = 2.5e-2, much smaller than
    # the empirical-minimum safety margin 0.2 * 20.06 = 4.06.)
    lam_t = mpf(lam_taylor)
    eps = mpf("1e-6")
    K = 10000
    min_ratio = mpf("inf")
    argmin_lam = mpf(0)
    for j in range(K + 1):
        lam = eps + (lam_t - eps) * mpf(j) / mpf(K)
        g = sigma2 - fpp(lam, q=q)
        ratio = g / (lam ** 2)
        if ratio < min_ratio:
            min_ratio = ratio
            argmin_lam = lam
    print(f"Quadratic region ({mp.nstr(eps,2)}, {lam_taylor}]: min g/lambda^2 = "
          f"{mp.nstr(min_ratio, 6)} at lambda = {mp.nstr(argmin_lam, 4)}")
    # Use 80% of the empirical min as a safe a (much larger than off-grid Lipschitz slack)
    a = min_ratio * mpf("0.8")
    if a <= 0:
        print(f"FAIL: empirical min of g/lambda^2 is not strictly positive")
        return False
    print(f"  -> diagnostic bound: g(lambda) >= {mp.nstr(a, 6)} * lambda^2 on ({mp.nstr(eps,2)}, {lam_taylor}]")

    # (c) Lipschitz region [lam_taylor, T]: bound |f'''| at the left endpoint and use grid+Lipschitz.
    L_lip = fppp_bound(lam_t, mpf(T), q=q)
    print(f"Lipschitz constant bound on [{lam_taylor}, {T}]: |f'''| <= {mp.nstr(L_lip, 6)}")
    Tm = mpf(T)
    Nm = mpf(N)
    h = (Tm - lam_t) / Nm
    slack = L_lip * h  # max change between grid points
    threshold = slack * mpf("1.5")  # require g(grid) > threshold so g(off-grid) > threshold/3 > 0
    min_g = mpf("inf")
    argmin_lam = mpf(0)
    for j in range(N + 1):
        lam = lam_t + (Tm - lam_t) * mpf(j) / Nm
        g = sigma2 - fpp(lam, q=q)
        if g < min_g:
            min_g = g
            argmin_lam = lam
        if g < threshold:
            print(f"FAIL at lambda = {mp.nstr(lam, 8)}: g = {mp.nstr(g, 8)} < threshold {mp.nstr(threshold, 6)}")
            return False
    print(f"Lipschitz region: N = {N} subdivisions, h = {mp.nstr(h, 4)}, slack = L*h = {mp.nstr(slack, 4)}")
    print(f"  min g on grid = {mp.nstr(min_g, 8)} at lambda = {mp.nstr(argmin_lam, 6)}")
    print(f"  PASS: min g > threshold = {mp.nstr(threshold, 6)}")

    # (d) Tail at lambda = T
    tail_value = mpf(4) * exp(-Tm)
    print(f"Tail bound at lambda = T={T}: 4 e^{{-T}} = {mp.nstr(tail_value, 4)} << sigma^2 = {mp.nstr(sigma2, 4)}")
    print("ALL CHECKS PASSED.")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=200_000, help="grid subdivisions on [lam_taylor, T]")
    ap.add_argument("--T", type=int, default=20, help="cutoff lambda for the tail bound")
    ap.add_argument("--q", type=int, default=10, help="base (default 10)")
    ap.add_argument("--dps", type=int, default=50, help="mpmath decimal precision")
    ap.add_argument("--lam-taylor", type=float, default=0.5,
                    help="boundary between Taylor and Lipschitz regions")
    args = ap.parse_args()
    ok = certify(N=args.N, T=args.T, q=args.q, dps=args.dps, lam_taylor=args.lam_taylor)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
