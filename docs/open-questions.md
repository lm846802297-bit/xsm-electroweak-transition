# Open theoretical audit questions

## Q1. Second-stage dimensional reduction and singlet-basis covariance

**Status: candidate issue identified; not yet resolved. Do not cite as a confirmed error.**

The first-stage hard-to-soft scalar matching (Appendix B, Eqs. B.5--B.10 of Ramsey-Musolf, Tenkanen and Tran) passes the exact singlet-translation covariance tests in this repository.

The next step integrates out temporal gauge fields. In the temporal 3D action, the retained operators include

\[
\mathcal L_{\rm temporal}\supset x_3 S_3 A_0^aA_0^a+y_3 S_3^2 A_0^aA_0^a
+x'_3S_3B_0^2+y'_3S_3^2B_0^2.
\]

The quoted hard-mode matching relations are

\[
x_3=\frac{\sqrt T}{(4\pi)^2}g^2a_1,\qquad
 y_3=\frac{T}{(4\pi)^2}\frac12g^2a_2,
\]

and similarly with g^2 -> g'^2 for the hypercharge coefficients. The 4D-to-3D field translation is delta_3=delta/sqrt(T), while the 4D portal map is a1' = a1 + 2 a2 delta.

If only the displayed temporal operators are retained, a direct shift of the temporal action requires

\[
x_3' = x_3+2y_3\delta_3.
\]

However, substituting the shifted 4D portal parameter into the quoted hard-mode matching relation gives

\[
x_3(a_1')-x_3(a_1)=\frac{2T}{(4\pi)^2}g^2a_2\delta_3
=4y_3\delta_3,
\]

which differs by a factor of two from the shift of the truncated temporal action.

This is an **operator-truncation covariance question**, not yet evidence of an error in the published matching. Possible resolutions to test are:

1. normalization and combinatorial conventions of temporal operators and their 1PI vertices;
2. higher-dimensional operators such as S_3^3 A_0^2 and their mixing under finite field translations;
3. power counting when the singlet background shift is not parametrically small;
4. whether the second-stage matching is intended only in the zero-VEV basis;
5. an independent DRalgo reproduction of temporal-sector vertices.

A subsequent hard-mode determinant and DRalgo source-level audit sharpened this question: the same hard Higgs loop multiplies the a1/2 S and a2/2 S^2 mass insertions, and the DRalgo xSM tensor test supports a factor-one-half candidate for the action coefficient x3. The algebraic mismatch is established; a complete gauge-covariant 1PI calculation and executed DRalgo run are still required before declaring a published error. See docs/theory/temporal-factor-two-audit.md.

The same issue can propagate into

\[
\bar b_{1,3}=b_{1,3}-\frac{1}{4\pi}(3m_Dx_3+m_D'x'_3),
\]

and must be resolved before using the second-stage EFT in a basis-independent scan.

## Q2. Finite-mass thermal tadpole

The LO term (a1+b3) T^2/12 assumes the massless high-temperature limit. A mass-dependent thermal function is needed when the effective singlet mass along the false-vacuum and bounce configurations is not small compared with T. Matching must separate hard and soft modes to avoid double counting the zero mode.

## Q3. Zero-temperature parameter conversion

The 2026 analysis includes a higher-order relation between Lagrangian parameters and physical observables. The complete reproduction must use the same one-loop zero-temperature renormalization prescription before comparing numerical benchmarks.

## Q4. Coordinate dependence of the path fraction

Liu--Xie classify a transition using a condition of the form |v_s^i/v_s^f| < 0.01. This ratio is not invariant under S -> S + delta, although Delta v_s = v_s^f-v_s^i is invariant. The quoted percentage is therefore coordinate-dependent. The reanalysis will compare free-energy differences, bounce action, nucleation/percolation temperatures, Delta v_s, alpha and beta/H.
