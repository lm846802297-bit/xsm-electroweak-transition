# xSM: replace the scalar zero mode once

This increment extends the [positive-mass, unresummed scalar audit](../positive-mass-scalar-tadpole-audit/README.md) without redoing its mass scan. It tests a **specified** positive scalar screening mass in the static, diagonal $H=0$ sector. The fixture self-energies are not computed xSM self-energies. It does not modify the EFT, include gauge-field rings, or estimate a phase transition or gravitational waves. Starting point: [xSM Draft #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1), exact head 96553f8bb36f07799f1cfc3fb230d58e8792a4e3, main 0e421748932d8b0a0ec635bc909a2f9be06f62c6.

## Splitting and derivative

For each real scalar mode at fixed $T$, write $X=M^2(s)$, $P=\Pi(s,T)$ and $Y=X+P$, with $X,Y>0$. Let $D=\partial_s X$ and $E=\partial_s P$. Four doublet modes share one $X$ at $H=0$; the singlet adds one mode. Use the unresummed free thermal potential $V_{\rm th}(X)$ **without** the 4D vacuum term. The DR-regulated 3D zero-mode contribution in 4D free-energy units is $Z(X)=-TX^{3/2}/(12\pi)$ per real degree of freedom.

The matched-accounting identity for this restricted one-loop potential is
$$
V_{\rm test} = \sum_i n_i \big[V_{\rm th}(X_i)-Z(X_i)+Z(Y_i)\big].
$$
In the same renormalization convention as the earlier audit, $V_{\rm th}-Z(X)$ corresponds to hard modes **minus** the zero-temperature vacuum potential, up to field-independent constants. It is not the complete hard effective potential. Its argument stays $X$, while the retained zero mode uses $Y$. This is the scalar zero-mode substitution in [Arnold and Espinosa, §2.1, Eqs. (2.4)–(2.7)](https://arxiv.org/pdf/hep-ph/9212235#page=6), compared to the explicit 3D soft determinant in [Schicho, Tenkanen and Österman, Eqs. (3.57), (B.2)–(B.3)](https://arxiv.org/pdf/2102.11145#page=50). The massless hard scalar contribution retained by [Liu and Xie, Eqs. (2.8)–(2.9)](https://arxiv.org/pdf/2101.10469#page=4) is not to be added a second time.

Differentiating gives
$$
\partial_s V_{\rm test}
=\sum_i n_i\left[\frac{D_i T^2}{2} F(\sqrt{X_i}/T)
-\frac{T}{8\pi}\big((D_i+E_i)\sqrt{Y_i}-D_i\sqrt{X_i}\big)\right],
$$
using the previous audit's thermal Bose kernel $F(x)$. The ring change is not correctly found by evaluating the original thermal kernel at $Y$ for every Matsubara mode. Keeping $E_i=0$ is correct only if the *used* screening shift is field-independent to the chosen accuracy. A field-dependent prescribed shift differentiates with the potential; dropping its derivative omits $-n_i T\sqrt{Y_i}E_i/(8\pi)$. This exercise does not by itself justify promoting a field-dependent screening prescription into xSM EFT matching.

## Double counting and first insertion

* Adding $Z(Y)$ to $V_{\rm th}(X)$ without subtracting $Z(X)$ duplicates the bare scalar zero mode. At the tadpole level the difference from correct replacement is **exactly** $\sum_i n_i\,\partial_s Z(X_i)$, even for $P=0$.
* The one-loop baseline here has **no** explicit self-energy insertion. Its replacement adds all zero-mode insertions starting with the term $P\,Z'(X)$. If a fixed-order baseline already contains precisely that one-insertion term in the same static conventions, subtract this overlap from the extra ring term as well. Its tadpole is
$$
n_i\partial_s[P_i Z'(X_i)]
=-\frac{n_iT}{8\pi}
 \left(E_i\sqrt{X_i}+\frac{P_iD_i}{2\sqrt{X_i}}\right).
$$
This additional subtraction is conditional on what the fixed-order baseline actually contains, not a universal second subtraction. At fixed positive $X_i$, replacing $(P_i,E_i)$ by $\epsilon(P_i,E_i)$ makes the subtracted difference quadratic in $\epsilon$. The numeric control checks the coefficient at zero and the factor-of-four decrease on halving a small $\epsilon$.
* When $M\sim gT$ and $P\sim g^2T^2$, $P/X$ is generally $O(1)$, so the small-$\epsilon$ overlap check is not a justification to truncate the physical soft mass replacement at one insertion. With $D/T\sim g^2$, the leading hard tadpole is $O(g^2T^3)$ and the soft replacement from $D$ is $O(g^3T^3)$. The $E$ term counts according to the order actually assigned to $E/T$. For example $E/T\sim g^4$ would put its soft contribution at $O(g^5T^3)$, not imply this audit computes all terms at that order. The standard leading thermal masses proportional to $T^2$ use $E=0$ in the usual leading xSM scalar counting. No statement is made near a mass zero, tachyonic background, or a critical region in which this loop expansion ceases to control long-distance physics.

## Numerical checks

The controlled fixtures use $T=120$ GeV, masses $(M_H,M_S)=(60,180)$ GeV, and $(D_H,D_S)=(22.5,-90)$ GeV. Screening raises the masses at the evaluation point to $(75,190)$ GeV. The field-dependent test takes $(E_H,E_S)=(7,-8)$ GeV; the constant shift test instead sets both to zero. These are synthetic numbers. The output preserves every species' signed terms before summation.

| Prescription | Total scalar derivative (GeV³) |
|---|---:|
| Unresummed | 19,893.85 |
| Replace the zero mode once, constant shifts | 17,745.26 |
| Replace the zero mode once, specified nonzero derivatives of the shifts | 14,975.96 |
| Incorrect, add the screened zero mode without subtraction (same nonzero derivatives) | 66,542.16 |
| Incorrect, omit $E_i$ in the field-dependent fixture | 17,745.26 |

At 65-digit working precision, the closed formula agrees with (i) a spatial radial integral with explicit DR massless-propagator subtraction for each zero mode, and (ii) an independent five-point field derivative of the directly integrated thermal logarithm plus the zero-mode potential replacement. Their maximum scaled differences over the supplied cases and species are $5.43\times10^{-65}$ and $1.22\times10^{-49}$, using $\max(1,|A|,|B|)$ in fixed GeV units. These are numerical consistency measures, not theoretical error bars.

## Provenance and replay

[DRalgo's latest stable release checked in this wake, v1.6.0](https://github.com/DR-algo/DRalgo/releases/tag/v1.6.0), resolves to cb27ed3fc48be05c056b3298d0dd930969d15e8a. Its Python 2.0 is described as outlook, not a stable release found in the published branches/releases. No supported local Wolfram kernel was found, and the project's Actions history showed ordinary tests, not DRalgo. No package installation, DRalgo/Mathics calculation, or CI dispatch was performed. This is a standalone numerical check using the cited analytic results.

The only numerical inputs are the bundled fixtures. Primary environment: GPT-managed Linux, Python 3.12.14, mpmath 1.3.0. Replay:

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_screened
    python screened.py /tmp/xsm-screened-replay.json
    cmp native-results.json /tmp/xsm-screened-replay.json

The byte comparison assumes consistent runtime versions.
