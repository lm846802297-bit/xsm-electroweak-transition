# Second-stage temporal EFT: factor-of-two audit

## Status

**Open candidate discrepancy.** The algebraic mismatch is established for the displayed operator basis and published LO matching relations. A final claim about an error in the published matching awaits a full independent 1PI calculation and an executed DRalgo comparison.

## Published temporal action and matching

The temporal part of the 3D action in Ramsey-Musolf, Tenkanen and Tran, JHEP 09 (2026) 260, Eq. (B.1), contains the operators

\[
\mathcal L_{3,\mathrm{temporal}}\supset
x_3 S_3 A_0^a A_0^a+y_3S_3^2A_0^aA_0^a
+x'_3S_3B_0^2+y'_3S_3^2B_0^2.
\]

The LO matching relations in Eqs. (B.18)--(B.21) are

\[
\begin{aligned}
x_3&=\frac{\sqrt T}{(4\pi)^2}g^2a_1,&
 y_3&=\frac{T}{(4\pi)^2}\frac12g^2a_2,\\
x'_3&=\frac{\sqrt T}{(4\pi)^2}g'^2a_1,&
 y'_3&=\frac{T}{(4\pi)^2}\frac12g'^2a_2.
\end{aligned}
\]

The 4D field translation S4=S4'+delta induces

\[
a'_1=a_1+2a_2\delta,\qquad \delta_3=\frac{\delta}{\sqrt T}.
\]

## Direct transformation of the 3D action

For the displayed operators, substituting S3=S3'+delta3 gives

\[
x'_3=x_3+2y_3\delta_3,\qquad y'_3=y_3.
\]

However, substituting the shifted 4D portal parameter into Eq. (B.18) gives

\[
\begin{aligned}
\Delta x_3^{\mathrm{paper}}
&=\frac{\sqrt T}{(4\pi)^2}g^2(2a_2\delta)\\
&=\frac{2T}{(4\pi)^2}g^2a_2\delta_3
=4y_3\delta_3.
\end{aligned}
\]

The action transformation requires 2y3 delta3. The mismatch is therefore

\[
\boxed{\Delta x_3^{\mathrm{paper}}-2y_3\delta_3=2y_3\delta_3.}
\]

The same result follows for the primed U(1) coefficient x3'. If Eqs. (B.18) and (B.19) were smaller by a factor of two, the displayed temporal action would transform covariantly while Eqs. (B.20) and (B.21) remain unchanged.

## Independent hard-mode determinant check

At leading order in the Higgs mass expansion, the singlet-dependent Higgs mass insertion is

\[
m_H^2(S)=m_{H,0}^2+\frac{a_1}{2}S+\frac{a_2}{2}S^2.
\]

The one-loop hard contribution to the coefficient of A0^2 is a common loop function multiplying this mass insertion. Consequently, the coefficients of S A0^2 and S^2 A0^2 must inherit the same common loop factor multiplying a1/2 and a2/2, respectively. This relative normalization supports the factor-one-half candidate for x3 when y3 is taken from Eq. (B.20). The full gauge-covariant 1PI calculation, including the seagull and two-linear-A0 diagrams, is still required before calling this a confirmed correction.

## DRalgo source-level cross-check

At DRalgo commit e00c28e3d01c83f2f2f205ba508b21709eef2d41:

- tests/smRealSinglet.m (blob 2f1ac03cb04effe1b548311325fffdfa320884f3) contains the expected temporal-coupling tensor components
  \[
  \texttt{lamVVSL[1]}=\frac{g^2\sqrt T\,\mu_m}{16\pi^2},\qquad
  \texttt{lamVL[5]}=\frac{g^2T\,\lambda_m}{8\pi^2}.
  \]
- The DRalgo manual identifies lamVVSL as a V^2 S tensor and lamVL as a V^2 S^2 tensor.
- The xSM example constructs scalar cubic and quartic tensors with GradCubic and GradQuartic. In the symmetric tensor convention, the corresponding action coefficients are the cubic tensor divided by 2! for two identical vector fields, and the quartic tensor divided by 2!2!.

Thus the corresponding action coefficients are

\[
x_3^{\mathrm{DRalgo}}=\frac{g^2\sqrt T\,a_1}{2(4\pi)^2},\qquad
y_3^{\mathrm{DRalgo}}=\frac{g^2T\,a_2}{2(4\pi)^2},
\]

which obey the finite-shift identity and agree with the factor-one-half candidate for x3. This is a source-level comparison, not an executed Mathematica result; the tensor-to-action normalization must be checked in a runnable DRalgo session.

## Required closure tests

1. Derive Gamma(A0 A0 S) and Gamma(A0 A0 S S) directly from the same hard-mode diagrams, retaining all symmetry factors.
2. Reproduce lamVVSL and lamVL numerically with DRalgo and convert its tensor convention to the coefficients in Eq. (B.1).
3. Compare the SU(2) and U(1) sectors independently.
4. Recompute the second-stage tadpole bbar(1,3) with the verified coefficients.
5. Only then decide whether the issue is a typo in Eqs. (B.18)--(B.19), a convention mismatch, or an additional matching contribution.

Primary sources: arXiv:2409.17554, Appendix B.1 and Eqs. (B.18)--(B.21); DRalgo examples/15_smRealSinglet.m, tests/smRealSinglet.m, and manual/main.tex at the commit above.
