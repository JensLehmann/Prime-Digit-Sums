"""Verify the Proth seed primes used in the paper's applications.

For each terminal t in {2, 5, 7}, the script checks the seed q_t = k * 2^n + 1 used
in the paper:
- Proth form: k odd, 2^n > k.
- 135-digit (resp. 134-, 132-digit) decimal expansion as displayed in the paper.
- Iterated digit-sum chain q_t -> ... -> t with every intermediate value prime.
- Proth primality witness: 7^((q_t - 1)/2) ≡ -1 (mod q_t).
- q_t >= 1.78e32, the theorem-level headline threshold from m_value.py.  This is stronger than
  comparing against the empirical M ≈ 1.48e32.

The seeds in SEEDS are sized at ~10^135, comfortably above M, so the q_t >= M
comparison holds with a ~10^100 margin.
"""

import sys

from m_value import M_THEOREM_HEADLINE_INT


SEEDS = {
    # Seed for A070027 infinitude; chain terminates at 7.
    7: dict(k=43917, n=430, displayed="121767334956703826188105215725011807003818007039279736531761951741852323974"
                                               "899660487046555992875497851303129811791960664013058089156609"),
    # Second, larger terminal-7 seed used as a robustness backup.
    "7'": dict(k=44853, n=450, displayed=None),
    # Seed for terminal 2.
    2: dict(k=3557, n=431, displayed=None),
    # Seed for terminal 5.
    5: dict(k=4027, n=426, displayed=None),
}

def is_prime_small(n: int) -> bool:
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


def digit_sum(n: int) -> int:
    return sum(int(d) for d in str(n))


def chain(n: int):
    seq = [n]
    while seq[-1] >= 10:
        seq.append(digit_sum(seq[-1]))
    return seq


def main() -> int:
    failures = 0
    for t in (2, 5, 7, "7'"):
        seed = SEEDS[t]
        k, n = seed["k"], seed["n"]
        q = k * (2 ** n) + 1
        digits = len(str(q))

        label = "7 (second seed)" if t == "7'" else str(t)
        print(f"--- Terminal t = {label}: q = {k} * 2^{n} + 1 ---")
        ok = True

        # Proth-form preconditions
        if k % 2 != 1:
            print(f"  FAIL: k = {k} is not odd"); ok = False
        if 2 ** n <= k:
            print(f"  FAIL: 2^{n} <= k = {k}"); ok = False

        # Size
        if q < M_THEOREM_HEADLINE_INT:
            print(f"  FAIL: q < theorem threshold"); ok = False
        else:
            ratio = q // M_THEOREM_HEADLINE_INT
            print(f"  OK: q has {digits} digits and q >= 1.78e32 (integer ratio floor {ratio})")

        # Displayed decimal expansion, included in the paper for the terminal-7 seed.
        if seed["displayed"] is not None:
            if str(q) == seed["displayed"]:
                print(f"  OK: displayed decimal expansion (length {digits}) matches str(q)")
            else:
                print(f"  FAIL: displayed expansion does not match str(q)")
                ok = False

        # Chain
        target_terminal = 7 if t == "7'" else t
        seq = chain(q)
        all_prime = all(is_prime_small(x) for x in seq[1:])
        terminal = seq[-1]
        print(f"  chain: q -> " + " -> ".join(str(x) for x in seq[1:]))
        if terminal != target_terminal:
            print(f"  FAIL: terminal {terminal} != target {target_terminal}"); ok = False
        elif not all_prime:
            print(f"  FAIL: not all intermediate values prime"); ok = False
        else:
            print(f"  OK: chain has all-prime intermediates terminating at {target_terminal}")

        # Proth witness
        witness = pow(7, (q - 1) // 2, q)
        if witness == q - 1:
            print(f"  OK: Proth witness 7^((q-1)/2) ≡ -1 (mod q) — q is prime")
        else:
            print(f"  FAIL: Proth witness check failed (got {witness}, expected q-1)")
            ok = False

        if not ok:
            failures += 1
        print()

    if failures == 0:
        print(f"All {len(SEEDS)} seeds verified.")
        return 0
    else:
        print(f"{failures} of {len(SEEDS)} seeds FAILED verification.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
