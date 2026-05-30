"""Search for Proth seed primes q_0 = k * 2^n + 1 >= M such that the digit-sum
chain s(q_0), s(s(q_0)), ... terminates at a specified one-digit prime t in
{2, 5, 7} with every intermediate value prime.

Usage:
    .venv/bin/python find_terminal_seeds.py [--targets 2,5] [--k-max 5000] [--n-range 425,445] [--max-candidates 10000000]
"""

from __future__ import annotations

import argparse
import sys
import time
from typing import Optional

from m_value import M_EMPIRICAL_INT, M_RIGOROUS_CERTIFIED_INT, M_THEOREM_HEADLINE_INT


THRESHOLDS = {
    "empirical": M_EMPIRICAL_INT,
    "rigorous": M_RIGOROUS_CERTIFIED_INT,
    "headline": M_THEOREM_HEADLINE_INT,
}


def digit_sum(n: int) -> int:
    return sum(int(d) for d in str(n))


def is_prime_small(n: int) -> bool:
    """Trial division for small n (used on intermediate digit-sum values < 10000)."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    p = 3
    while p * p <= n:
        if n % p == 0:
            return False
        p += 2
    return True


def chain_terminal(n: int) -> Optional[int]:
    """Compute the iterated digit-sum chain starting at n.
    Returns terminal one-digit value if the entire chain (including n if multi-digit)
    consists of primes; otherwise returns None.

    The starting value n itself is NOT required to be prime (this is checked separately
    by the Proth witness for the candidate q_0).
    """
    cur = n
    while cur >= 10:
        cur = digit_sum(cur)
        if not is_prime_small(cur):
            return None
    return cur if cur in (2, 3, 5, 7) else None


def proth_witness(q: int, witness: int = 7) -> bool:
    """Proth's theorem: q = k * 2^n + 1 with k odd, 2^n > k is prime iff
    there exists a witness `a` with a^((q-1)/2) ≡ -1 mod q.

    Returns True if witness `a` confirms primality. False does NOT prove compositeness
    (need other witnesses), but any confirming supplied witness proves q prime.
    """
    return pow(witness, (q - 1) // 2, q) == q - 1


def proth_is_prime(q: int, witnesses: tuple = (7,)) -> Optional[int]:
    for a in witnesses:
        if proth_witness(q, a):
            return a
    return None


def search_seeds(targets: tuple, k_max: int, n_min: int, n_max: int, max_candidates: int,
                 threshold_name: str, witnesses: tuple):
    threshold = THRESHOLDS[threshold_name]
    print(f"# Searching for Proth seeds q_0 = k * 2^n + 1 >= {threshold_name} threshold ({threshold})")
    print(f"# k odd, 1 <= k <= {k_max}; n in [{n_min}, {n_max}]; targets = {targets}")
    print(f"# Proth witnesses tried, in order: {witnesses}")
    print(f"# Will examine up to {max_candidates} (k, n) pairs.")
    print()

    t0 = time.time()
    examined = 0
    seeds_found = {t: None for t in targets}

    for n in range(n_min, n_max + 1):
        for k in range(1, k_max + 1, 2):  # k odd
            if 2 ** n <= k:
                continue
            q0 = k * (2 ** n) + 1
            if q0 < threshold:
                continue
            examined += 1
            if examined > max_candidates:
                print(f"# Examined {max_candidates} candidates; stopping.", flush=True)
                return seeds_found

            ds = digit_sum(q0)
            # quick reject: ds must be a prime <= 10000
            if ds < 2 or ds > 9999:
                continue
            if not is_prime_small(ds):
                continue
            t = chain_terminal(q0)
            if t not in targets:
                continue
            if seeds_found.get(t) is not None:
                continue  # already found one for this terminal
            # Now do the expensive Proth primality test
            witness = proth_is_prime(q0, witnesses=witnesses)
            if witness is not None:
                seeds_found[t] = (k, n, q0, ds, witness)
                print(f"# FOUND terminal-{t} seed: k={k}, n={n}, s(q_0)={ds}, witness={witness}", flush=True)
                print(f"#   q_0 = {q0}")
                # Check if all targets found
                if all(seeds_found[x] is not None for x in targets):
                    print(f"# All targets found in {time.time()-t0:.1f}s after examining {examined} pairs.")
                    return seeds_found

        if n % 5 == 0:
            elapsed = time.time() - t0
            print(f"# ... progress: n={n}, examined={examined} pairs, elapsed={elapsed:.1f}s", flush=True)

    print(f"# Search ended. Examined {examined} pairs in {time.time()-t0:.1f}s.")
    return seeds_found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=str, default="2,5", help="Comma-separated terminal targets")
    ap.add_argument("--k-max", type=int, default=5000)
    ap.add_argument("--n-range", type=str, default="425,445")
    ap.add_argument("--max-candidates", type=int, default=10**7)
    ap.add_argument("--threshold", choices=sorted(THRESHOLDS), default="headline",
                    help="Minimum q_0 threshold to enforce; headline is the theorem-level 1.78e32 bound")
    ap.add_argument("--witnesses", type=str, default="7",
                    help="Comma-separated Proth witnesses to try, in order")
    args = ap.parse_args()

    targets = tuple(int(x) for x in args.targets.split(","))
    n_min, n_max = (int(x) for x in args.n_range.split(","))
    witnesses = tuple(int(x) for x in args.witnesses.split(","))

    seeds = search_seeds(targets, args.k_max, n_min, n_max, args.max_candidates,
                         args.threshold, witnesses)
    print()
    print("# === Results ===")
    for t, seed in seeds.items():
        if seed is None:
            print(f"# Terminal {t}: NO seed found")
        else:
            k, n, q0, ds, witness = seed
            print(f"# Terminal {t}: k={k}, n={n}, s(q_0)={ds}, witness={witness}")
            print(f"#   q_0 = {q0}")


if __name__ == "__main__":
    main()
