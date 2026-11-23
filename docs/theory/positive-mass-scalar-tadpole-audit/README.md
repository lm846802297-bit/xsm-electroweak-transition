# xSM: positive-mass correction to the scalar thermal tadpole

## Scope

[Liu and Xie, Eqs. (2.8)–(2.9)](https://arxiv.org/pdf/2101.10469#page=4) retain the gauge-invariant $T^2$ terms and include a singlet tadpole $(a_1+b_3)T^2 S/12$. This is a check of the finite positive-mass scalar one-loop contribution underlying that term, limited to the symmetric doublet background $H=0$. It is not the [temporal gauge threshold](../soft-temporal-threshold-audit/README.md), a new full 3D parameter match, a scan of paper benchmarks, or a phase-transition result.

Starting repository: [draft PR #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1) at 7d6bdaac820470a77425b064f11f6b83f6caa667, main 0e421748932d8b0a0ec635bc909a2f9be06f62c6. The [open question](../../open-questions.md#q2-finite-mass-thermal-tadpole) concerned matching this correction with separate hard and soft modes. Here we evaluate the free scalar determinant contribution at prescribed **nonnegative** masses. In particular, do not analytically continue these real numerical kernels through a tachyonic background or identify the fixture masses with the physical masses of the false vacuum.

## Action normalization

For each real scalar, the thermal part of the free determinant obeys $2\,\partial V_{\rm th}/\partial M^2=T^2 F(M/T)$. At $H=0$, there are four real degenerate doublet modes and one singlet, without Higgs–singlet mass mixing. With the 4D action convention of Liu and Xie,
$$
\begin{aligned}
M_H^2(s)&=-\mu^2+\tfrac12a_1s+\tfrac12a_2s^2,&
M_S^2(s)&=b_2+2b_3s+3b_4s^2,\\
\alpha_H&=2\partial_sM_H^2=a_1+2a_2s,&
\alpha_S&=\tfrac12\partial_sM_S^2=b_3+3b_4s,\\
\partial_s V_{\rm scalar,th}(s,T)&=
T^2[\alpha_H F(x_H)+\alpha_S F(x_S)],&
x_i&=M_i(s)/T .
\end{aligned}
$$
Under a constant field translation applied to the action and evaluation point together, the masses and slopes at that physical point are invariant. The formula is a local derivative; away from $s=0$ it is not a globally linear term. $F(0)=1/12$ gives exactly $m_1 T^2$ at $s=0$, and the local linear-in-$s$ term in $\partial_s V_{\rm scalar,th}$ reproduces the scalar contribution $c_ST^2 s$ in (2.9). The bosonic kernel is
$$
F(x)=\frac1{2\pi^2}\int_0^\infty
\frac{q^2\,dq}{\sqrt{q^2+x^2}(\exp[\sqrt{q^2+x^2}]-1)}
=\frac{x}{2\pi^2}\sum_{n=1}^\infty\frac{K_1(nx)}{n}.
$$
The integrand is positive and strictly decreases with $x>0$. If $\alpha_H=-\alpha_S\ne0$, the finite-mass thermal scalar derivative cancels only at equal masses. An exact cancellation in the leading approximation therefore does not bound the absolute correction.

## Hard/soft identity and high-temperature error bound

For a common positive input mass and $\nu=\bar\mu/T$ define $F=I_{\rm th}/T^2$ and split the full one-loop Matsubara tadpole into its hard contribution, soft zero mode, and subtracted zero-temperature vacuum term:
$$
F(x)=F_{\rm hard,\overline{MS}}(x;\nu)-\frac{x}{4\pi}
-F_{\rm vac,\overline{MS}}(x;\nu).
$$
Using $\int_p(p^2+M^2)^{-1}=-M/(4\pi)$ in 3D dimensional regularization, the linear term is solely the scalar zero mode. The vacuum subtraction is $x^2[2\ln(x/\nu)-1]/(16\pi^2)$. The hard term through order $x^2$ is $1/12-x^2[\ln(\nu/(4\pi))+\gamma_E]/(8\pi^2)$. See [Schicho, Tenkanen and Österman, Eqs. (B.2), (B.3), (B.9), (B.21)](https://arxiv.org/pdf/2102.11145#page=50), for the primed sum, scalar zero-mode determinant, and hard mass expansion. Cancellation of $\nu$ gives
$$
F(x) = F_{\mathrm{HT2}}(x)+R(x),\qquad
F_{\mathrm{HT2}}(x)=\frac1{12}-\frac{x}{4\pi}
+\frac{x^2}{8\pi^2}[\ln(4\pi/x)-\gamma_E+1/2].
$$
After subtracting the first two terms from the nonzero-mode propagator, its exact remainder integrand is $M^4/[P^4(P^2+M^2)]>0$, with $P^2\ge (2\pi T)^2$. Thus
$$
\frac{\zeta(3)x^4}{128\pi^4[1+(x/2\pi)^2]}\le R(x)
\le\frac{\zeta(3)x^4}{128\pi^4}.
$$
The finite identity and bound hold at any positive $x$; the geometric **hard** series in $M^2$ is guaranteed convergent for $|x|<2\pi$. Neither fact makes $F_{\mathrm{HT2}}$ useful at large $x$ where the true thermal result is Boltzmann suppressed. The script also integrates each hard pair at fixed $n$ analytically in spatial momentum, given by $[w_n+x^2/(2w_n)-\sqrt{w_n^2+x^2}]/(2\pi)$ with $w_n=2\pi n$. The code rationalizes the subtraction, sums 64 explicit pairs, and applies the same inequality using $n\ge65$ to enclose the Hurwitz-zeta tail.

**Double counting:** In a dimensionally reduced calculation, include nonzero modes in hard matching and count the scalar zero mode once in the soft/ultrasoft theory at its consistently matched mass. Do not add the full $F(x)$ to the existing $T^2/12$ tadpole or to a calculation already containing that soft scalar determinant. Likewise, the $\overline{\rm MS}$ vacuum term subtracted here must not be interpreted as a thermal mass correction. The identity uses the *same* mass everywhere; replacing only the zero-mode mass with a screened mass requires a consistent matching or subtract-and-add ring prescription. Plugging a thermal screening mass into all Matsubara frequencies, or holding only its derivatives fixed, is not a complete higher-order resummation. No such correction is made to this repo's EFT or published phase-transition analysis.

## Numerical checks and implications

Independent quadrature, a truncated Bessel sum with analytic positive tail bound, and an interval reconstruction from the nonzero Matsubara pairs agree on the supplied grid $0\le x\le4$. Maximum radial/Bessel discrepancy is $2.78\times10^{-17}$ in $F$. Tests also compare three points with a 55-decimal-digit quadrature and check the split for three values of $\nu$. Analytic tail bounds are for the exact arithmetic truncation; the floating-point checks separately allow a numerical discrepancy. No Mathematica/Wolfram kernel was found on this Linux machine; DRalgo was not executed.

| $M/T$ | $F(x)/F(0)$ |
|---:|---:|
| 0.10 | 0.9117360 |
| 0.25 | 0.7977487 |
| 0.50 | 0.6421773 |
| 1.00 | 0.4191490 |
| 1.50 | 0.2737726 |
| 2.00 | 0.1782425 |

Ten percent suppression relative to the massless scalar loop occurs at $M/T=0.11433755$; this is not a 10% error threshold for an electroweak observable. At $M/T=1$, the complete combination through $x^2$ differs from the exact scalar kernel by 0.273% of the exact value; the relative truncation error grows to 9.96% at 2 and exceeds 100% at 3. The linear zero-mode term alone is not a reliable extrapolation. As a cancellation stress test at $T=120$ GeV with local slopes $+45$ and $-45$ GeV, equal masses give zero, while $x_H=0.5,\ x_S=1.5$ give $+19{,}893.85$ GeV$^3$ despite zero at massless order. These are synthetic mass/slope fixtures, not predictions about the physical frequency or size of this cancellation. Neglecting the mass dependence, or using a relative ratio to a near-zero leading sum, can hide these absolute shifts.

## Replay and provenance

Numerical sources are the definitions above; the scalar input grid and fixtures are in `inputs.json`. The native result includes each sector separately, the quadrature estimates, analytic bounds, split pieces, and source/checksum metadata. Primary run and tests used GPT-managed Linux, Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and mpmath 1.3.0. GitHub Actions was not run.

    sha256sum -c SHA256SUMS.txt
    OPENBLAS_NUM_THREADS=1 python -m unittest -v test_tadpole
    OPENBLAS_NUM_THREADS=1 python tadpole.py /tmp/xsm-tadpole-replay.json
    cmp native-results.json /tmp/xsm-tadpole-replay.json

Exact byte replay assumes matching runtime versions.
