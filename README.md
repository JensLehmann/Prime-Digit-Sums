# Prime-Digit-Sums — verification code

Companion code for the paper *An Explicit Surjectivity Threshold for Digit Sums of Primes*. These scripts verify the explicit constants stated in the paper:

- **the seed primes $q_t$** used in the proof of infinitude of OEIS A070027 and in the terminal-stratification application, and
- **the explicit threshold $M$** in the main surjectivity theorem, via the constraint chains (C1)–(C8) for the moment-comparison threshold $x_{45}$ and (C9)–(C13) for the characteristic-function threshold $x_{43}^\ast$, together with the major/minor-arc thresholds $Y_{\mathrm{maj}}$, $Y_{\mathrm{min}}$, the Step-5 absorption threshold $Y_{47}$, and the AP threshold $x_{\mathrm{AP}}$ as defined in the derivation of $M$.

## Setup

Tested with Python 3.9.6 on macOS (Darwin 25.3). Only dependency is `mpmath` (arbitrary-precision arithmetic).

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

## Reproducing the paper's numerical claims

### One-command smoke test

Runs the theorem certificate, Proth seed verification, and stored-$M$ regression in sequence:

```sh
.venv/bin/python run_all.py
```

Expected runtime: < 5 seconds.

### Theorem-level certificate

This is the primary reproducibility entry point for the theorem-level constants. It prints the
rigorous constants used in the manuscript, checks the displayed numerical inequalities by exact
rational or outward-rounded interval arithmetic, verifies the Proth seed certificates, and writes a
machine-readable JSON certificate:

```sh
.venv/bin/python certify_theorem.py
```

By default this writes:

```text
certificates/theorem_certificate.json
```

The script distinguishes:
- **proof-relevant constants:** the closed-form bounds used in the theorem, including
  $Y_{43}^{\ast}\le 9.10\times 10^{31}$ and
  $M\le 177843590339381623423137292296355<1.78\times10^{32}$;
- **machine-checkable certificates:** rational/interval checks for the rounded inequalities and
  exact Proth witness checks for the seed primes;
- **empirical diagnostics:** the sharper bisection values, including
  $M_{\mathrm{emp}}\approx 1.48\times10^{32}$, which are not used in the proof.

The JSON also exposes named records for the threshold components
($Y_{41},Y_{42},Y_{43}^{\ast},Y_{45},Y_{47},Y_{\mathrm{maj}},Y_{\mathrm{min}},x_{\mathrm{AP}}$)
so the numerical layer can be compared directly with the manuscript's threshold table.

Expected runtime: < 2 seconds.

The code verifies numerical inequalities, integer threshold conversions, deterministic Proth
certificates, and reproducibility of the exploratory bisection diagnostics. It does not formalize
the analytic reductions from DMR/MR10/IK; those reductions are proved in the manuscript, with source
equation locators recorded in the audit appendix.

### Seed primes for each terminal $t \in \{2, 5, 7\}$

Verifies the four seed primes used in the paper: one per terminal one-digit prime in
the iterated digit-sum chain, plus a robustness backup for $t=7$. The seed for $t=7$
is the one in the proof of the A070027 theorem; the backup for $t=7$ appears in the
following remark; the seeds for $t=2,5$ appear in the terminal-stratification application.

```sh
.venv/bin/python verify_q0.py
```

For each seed $q_t = k \cdot 2^n + 1$ the script checks:
- Proth form: $k$ odd, $2^n > k$.
- Decimal expansion (for $t = 7$, additionally cross-checks the displayed 135-digit string).
- Iterated digit-sum chain: $q_t \xrightarrow{s} \cdots \xrightarrow{s} t$ with every
  intermediate value prime.
- **Proth primality witness:** $7^{(q_t - 1)/2} \equiv -1 \pmod{q_t}$. This is a *deterministic* primality certificate via Proth's theorem (1878), not a probabilistic Miller–Rabin test: once verified, $q_t$ is unconditionally prime.
- $q_t \ge M$ (the threshold from the main surjectivity theorem).

The four seeds:
| Terminal $t$ | $k$ | $n$ | digits | chain |
|---|---|---|---|---|
| $2$ | $3557$ | $431$ | $134$ | $q_2 \to 641 \to 11 \to 2$ |
| $5$ | $4027$ | $426$ | $132$ | $q_5 \to 599 \to 23 \to 5$ |
| $7$ | $43917$ | $430$ | $135$ | $q_7 \to 601 \to 7$ |
| $7$ (backup) | $44853$ | $450$ | $141$ | $q_7' \to 601 \to 7$ |

Expected runtime: < 1 second.

### Searching for new seed primes

The seeds for $t = 2, 5$ were found by `find_terminal_seeds.py`, which scans
Proth-form candidates $k\cdot 2^n + 1$:

```sh
.venv/bin/python find_terminal_seeds.py --targets 2,5 --k-max 5000 --n-range 425,445
```

Found both seeds in under one second on a 2024-era laptop.  The reported seeds are
the *first* hits returned by the search (loop order $n$ first, then odd $k$, over the box
$(n,k)\in[425,445]\times[1,5000]$); they are not canonical in any other sense.
Different box parameters generally yield different seeds.

### Main theorem — constraint checker and threshold $M$

Computes the explicit constants and empirical/bisected thresholds in the derivation of $M$ using `mpmath` at 90 decimal digits of precision. This script is useful for exploring the constraint chain and reproducing the sharper diagnostic values printed in the paper. The theorem-level certificate is `certify_theorem.py` above. The defaults are the paper's working parameters:

```sh
.venv/bin/python prime_digit_sums_constraint_checker_with_search.py
```

(equivalent to `--eta 0.0545 --nu 0.2859`, the paper's working parameters.) With these, the script computes:

```
Y41         ≈ 4.24e+01      (major-arc proposition, Step 2)
Y42         ≈ 2.37e+04      (major-arc proposition, Step 4)
Y43*        ≈ 6.27e+31      (CF-replacement threshold, empirical bisection)
Y45         ≈ 3.54e+09      (moment-comparison threshold)
Y47         ≈ 3.61e+31      (Step 5 absorption threshold)
Y_major     ≈ 7.44e+31      (solving the major inequality)
Y_minor     ≈ 7.59e+31      (solving the minor inequality) ← binding
Y_AP        ≈ 8.31e+00      (AP lower-bound threshold log 4050)
c43*        ≈ 1.91e-05      (1/200 of the maximum admissible value)
M           ≈ 1.48e+32      (= ⌈(9/2) · Y_required / log 10 + 9/2⌉)
```

Expected runtime: < 1 second.

The near-binding thresholds are $Y_{43}^{\ast}$, $Y_{\mathrm{min}}$, and $Y_{\mathrm{maj}}$, with $Y_{\mathrm{min}}$ empirically binding. The headline bound in the main theorem is $M < 1.78\times10^{32}$ (rigorous, with the closed-form certificate giving $M \le 177843590339381623423137292296355$ from $Y_{\ast}\le 9.10\times10^{31}$); the empirical $M \approx 1.48 \times 10^{32}$ is reported only as a sharper numerical observation.

### Optional: search over $(\eta, \nu)$

To explore alternative $(\eta,\nu)$ in the admissible region:

```sh
.venv/bin/python prime_digit_sums_constraint_checker_with_search.py --search \
    --search-nu-min 0.20 --search-nu-max 0.34 \
    --search-eta-min 0.03 --search-eta-max 0.08 \
    --search-grid-nu-step 0.005 --search-grid-eta-step 0.002
```

Expected runtime: ~minutes.

## Implementation notes

### Map from paper constants to code locations

| Paper symbol | Definition in paper | Code location |
|---|---|---|
| $c_4$ | $\log q / 6$, decimal-base constant in `lem:c4-explicit` | `derive()`, `c4 = logq / 6` |
| $\theta$ | $c_4 / (16(\eta + 1/2))$, in the CF-replacement proposition | `derive()`, `theta = c4 / (16 * (eta + 0.5))` |
| $c_{43}^{\ast}$ | $\tfrac{1}{200}\cdot(c_4/32)\min\{1,(\nu-2\eta)/(\eta+1/2)\}$ | `derive()`, `c43_safety = mpf("0.005")` and downstream |
| $C_{\mathrm{DMR}}$ | $102$, explicit DMR Lemma 4.3 input (`lem:DMR43`) | `Inputs.CDMR: mp.mpf = mpf("102.0")` |
| $C_{\tau}$ | $10^3$, explicit $X^{1/3}$ short-interval divisor-sum input for carry propagation | used only in the paper's Type-II constant bookkeeping |
| $C_{\mathrm{II}}$ | $<64000$, Type-II minor-arc constant (`prop:typeII-explicit`) | still dominated by the Type-I constant in `lem:Cmin`, so not a separate script input |
| $C_{\min}$ | $4\cdot 10^{12}$, explicit minor-arc prime-sum bound (`lem:Cmin`) | `Inputs.Cmin: mp.mpf = mpf("4000000000000.0")` |
| $C_{\mathrm{maj}}$ | $510$ rounded envelope; decomposed proof constants $C^{(1)}<460$, $C^{(0)}<0.52$, numerical major-arc corollary | `derive()`, `C1_alpha = pi*(180+sqrt(33)) * logq**(-nu)` (the decomposed constants are what enter $Y_{\mathrm{maj}}$ and the sharper theorem headline; the rounded $C_{\mathrm{maj}}=510$ remains a proposition-level envelope) |
| $Y_{43}^{\ast}$ | $L_{\ast}^{(4)} \log q$, analytical threshold corollary | `compute_Y43_star()` |
| $Y_{\mathrm{maj}}$ | explicit major-arc positivity threshold | `compute_Y_major()` |
| $Y_{\mathrm{min}}$ | explicit minor-arc threshold | `compute_Y_minor()` |
| $M$ | main theorem, with integral digit length $L_m=\lfloor 2m/9\rfloor$ | `certify_theorem.py` prints and certifies the rigorous theorem-level value; `compute_M()` returns the empirical bisected value; `m_value.py` stores both literals |

### Safety margin on $c_{43}^{\ast}$

The characteristic-function comparison replacing DMR Prop. 4.1 proves its bound for **any** $c_{43}^{\ast} < \frac{c_4}{32}\min\{1,\ (\nu-2\eta)/(\eta+1/2)\}$. The script uses a $0.5\%$ safety margin (i.e., $c_{43}^{\ast} = \frac{1}{200}$ of the max admissible) so the strict inequality is realized as a checkable predicate at finite $L$ while keeping the absorption threshold $Y_{47}$ below the major/minor scale:
```python
c43_safety = mpf("0.005")
c43_star = c43_safety * (c4/32) * min(1, (nu-2*eta)/(eta+0.5))
```
The choice balances the closed-form $L_{\ast}^{(4)}$ certificate against the $Y_{47}$ absorption threshold.

### Erdős–Turán constants

The script implements the interval-discrepancy bound used in the moment-comparison step with the Kuipers–Niederreiter interval-discrepancy constants $(6, 4/\pi)$.

### Major-arc treatment

The major-arc step uses the paper's inline-keep form of Step 5 of the major-arc proposition: the exponential factor $e^{-c_{43}^\ast L^\nu}$ is folded into the overall threshold $x_{\mathrm{maj}}$ rather than being absorbed into a separate $x_{47}^\ast$. The value $C_{\mathrm{maj}}=510$ is fixed in the numerical corollary of that proposition and is not exposed as a CLI flag.

### Strict sub-Gaussianity of the uniform digit

The paper's sharp sub-Gaussian lemma uses the MGF bound $\mathbb{E}\,e^{\lambda(Z-\mu)} \le e^{\sigma^2 \lambda^2 / 2}$ for the uniform digit distribution (variance proxy = true variance), then derives the even-moment bound by an explicit tail integral.

The paper now proves this lemma analytically. The script `verify_subgaussian.py` is retained only as a floating-point diagnostic for the shape of the same inequality at $q = 10$:

- a quadratic-bound region $\lambda \in (0, 0.5]$: $g(\lambda) := \sigma^2 - f''(\lambda) \ge 16\,\lambda^2$, checked by $10^4$ grid samples with Lipschitz extrapolation;
- a Lipschitz region $\lambda \in [0.5, 20]$: $g(\lambda) \ge 5.0$ on a grid of $2 \times 10^5$ samples with $|f'''| \le 23$ Lipschitz slack.

The tail bound for $\lambda \ge 20$ is closed-form ($4 e^{-20} \approx 8 \times 10^{-9} \ll \sigma^2 = 8.25$). Run:

```sh
.venv/bin/python verify_subgaussian.py
```

Expected runtime: < 10 seconds.

## Files

| File | Purpose |
|---|---|
| `run_all.py` | One-command smoke test for the theorem certificate, seed certificates, and $M$ regression. |
| `verify_q0.py` | Seed-prime verification for terminals $t\in\{2,5,7\}$; Proth witness and digit chain checks. |
| `certify_theorem.py` | Single theorem-level reproducibility script; writes `certificates/theorem_certificate.json`. |
| `find_terminal_seeds.py` | Search over Proth-form seed primes leading to a prescribed terminal one-digit prime. |
| `prime_digit_sums_constraint_checker_with_search.py` | Main-theorem constraint checker and $(\eta,\nu)$ grid search; computes $M$. |
| `verify_subgaussian.py` | Floating-point diagnostic for the sharp sub-Gaussian lemma (the paper proof is analytic). |
| `m_value.py` | Single source of truth for $M$ (canonical literal used by the seed scripts). |
| `test_m_value.py` | Regression test verifying `M_LITERAL` matches a fresh `compute_M(Inputs())` run; run directly (no pytest needed). |
| `requirements.txt` | Python dep `mpmath` (`>=1.3,<2`; $M$ reproduces under 1.3.x and 1.4.x). |
| `certificates/theorem_certificate.json` | Machine-readable certificate generated by `certify_theorem.py`. |

## Reproducibility

The theorem-level certificate is generated by `certify_theorem.py` and recorded in
`certificates/theorem_certificate.json`.  The exploratory constraint checker can be rerun directly
at the paper's working parameters $(\eta, \nu) = (0.0545, 0.2859)$ to reproduce the empirical
constants reported in the derivation of $M$. The seed verification (`verify_q0.py`) and the
regression test (`test_m_value.py`) both run to completion under one second with deterministic
output.
