# Static hard temporal portal audit in the xSM

## Scope and inspected baseline

Pinned open [draft PR #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1), head f365f6021e5a469c189aaf834ad7acf419dc6ab6, on research/eft-matching-audit; main was 0e421748932d8b0a0ec635bc909a2f9be06f62c6. This checkout already contains LO/NLO first-stage basis covariance, symbolic and 60-digit checks, and the 100-point native ensemble. This increment is a separate evaluation of the open static 1PI temporal couplings. The prior [archive-covariance run](https://github.com/lm846802297-bit/xsm-electroweak-transition/actions/runs/45052620212) succeeded on 706a5d995d8801887c3e65c99a771d037de6c752; f365f602 is one later commit, changing only the provenance, manifest and archived pytest output. The two scalar xSM coefficients below are not a substitute for the full transition and GW calculation.

Primary sources: Liu and Xie [Eq. (2.8)–(2.9)](https://arxiv.org/pdf/2101.10469#page=4); Schicho, Tenkanen and Österman [Eqs. (4.8), (4.24)–(4.27), (4.36)–(4.39)](https://arxiv.org/pdf/2102.11145); Ramsey-Musolf, Tenkanen and Tran [Eqs. (B.1), (B.18)–(B.21)](https://arxiv.org/pdf/2409.17554#page=26). The first paper's Eq. (2.8) is the independent 4D high-temperature tadpole issue, not the temporal portal matching examined here.

## Derivation: both static 1PI contributions

Use a spatially uniform singlet S and temporal gauge field, zero background electroweak doublet, dimensional regularization with d=3−2ε spatial dimensions, and nonzero bosonic Matsubara frequencies only. Let Q = g A₀ᵃ τᵃ/2 + g′ B₀ I/2. The Higgs doublet comprises two complex fields, so its one-loop contribution is Tr log, not ½ Tr log. Its scalar mass in this background is

    M²(S) = cH + (a1/2) S + (a2/2) S².

For a Matsubara term the operator is D + 2ω Q + Q², where D=p²+ω²+M². The coefficient through second order in the temporal background is

    tr(Q²) Σhard [D⁻¹ − 2ω² D⁻²].

The first term is the seagull; the second is the two-linear-vertex term, including its relative sign and log-expansion symmetry factor. Differentiating in M² creates the static hard Higgs mass insertion. At M²=0 the derivative is

    −I₂ + 4 Σhard ω² / (p²+ω²)³ = (3−d) I₂.

Here I_n=Σhard (p²+ω²)^(-n). Spatial integration by parts gives Σ p²/P⁶ = (d/4) I₂. In the matched massless limit, I₂ = 1/[(4π)² ε]+O(1). The finite 2/(4π)² survives via the product of an evanescent difference and the pole. Setting d=3 in the two contributions before the Matsubara sum incorrectly erases it. The result concerns the local, zero-external-momentum static 1PI terms. At H=0 the background singlet couples directly only to Higgs fluctuations at this order; pure gauge, ghost and fermion determinant terms do not carry S. This does not address higher-order gauge dependence or field rescalings multiplying loop-generated terms.

For Pauli generators and Y=I/2, both tr(τᵃτᵇ/4)=δᵃᵇ/2 and tr(Y²)=1/2; the mixed trace vanishes. Convert 4D zero modes with field_4=√T field_3 and L_3=L_4/T. In the action convention x3 S3 A0² + y3 S3² A0² with no factorials this gives

    x3 = g² √T a1 / [2 (4π)²],       x3′ = x3 with g² → g′²,
    y3 = g² T a2 / [2 (4π)²],        y3′ = y3 with g² → g′².

The printed 2026 linear couplings are twice these values, while their quadratic couplings agree. Schicho et al. (4.36)–(4.39) display the same linear normalization. That paper also retains a contribution proportional to a1² in y3. The second background mass derivative of the temporal kernel is (2d−10) I₃; I₃=2 ζ(3)/[(4π)^4 T²] at d=3. Thus the independently derived a1² contribution to y3 is −g²T ζ(3) a1²/[2(4π)^4 T²], agreeing with (4.38), with the same U(1) replacement in (4.39). Retain the same mass/power counting when comparing these truncations. The leading coefficients obey x3(a1+2a2 δ)=x3(a1)+2y3 δ/√T. With other matched quantities fixed, inserting the verified x3 candidate in the *displayed* leading second-stage tadpole formula halves the temporal-gauge part of that threshold, not the complete tadpole or any observable.

## Independent finite-mass check

For a unit-charge complex scalar, use the finite-temperature determinant at imaginary chemical potential. Its temporal A0² coefficient is T⁻¹ ∫ nB(E)[1+nB(E)] d³p/(2π)³. Differentiating it and subtracting the bosonic zero-mode derivative −T/(8π M) yields the hard derivative. The radial quadrature and differentiated Bessel-K series are separately implemented. For M/T=0.02 to 4 on the supplied grid, the methods agree within 2.3e−16 absolute; with M/T=0.02 the hard derivative is 0.0126649937036292, approaching 1/(8π²)=0.0126651479552922. The first finite-mass coefficient is −8 ζ(3)/(4π)^4 for (M/T)². These test masses are positive; do not substitute tachyonic background masses in the holonomy formula.

The analytic result and numerical checks strongly support a factor-one-half candidate for the displayed linear action coupling at the stated order. No Mathematica/DRalgo execution, contact with authors, full finite-mass scalar tadpole implementation, transition scan, or change to the repository is claimed.

## Replay

With Python 3.12, SciPy 1.17.0, SymPy 1.14.0 and mpmath 1.3.0:

    python -m unittest -v test_temporal_static
    python temporal_static.py replay.json
    cmp replay.json static-temporal-results.json
    sha256sum -c SHA256SUMS.txt

The byte comparison assumes the same runtime versions; use a numerical tolerance when version metadata differs.
