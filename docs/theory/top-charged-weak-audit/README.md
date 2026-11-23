# Top pole matching: the charged W and Goldstone block

This audit evaluates only the charged W and charged Goldstone, massless-bottom 1PI terms of the top-quark pole projection in [Niemi, Schicho and Tenkanen (NST)](https://doi.org/10.1103/PhysRevD.103.115035), [arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3#page=22), Appendix A, (A20)–(A23), (A26), (A30), (A36)–(A37), (A39), (A42) and (A61). The [published erratum](https://doi.org/10.1103/PhysRevD.109.039902) does not alter the selected terms. We work in NST's MS-bar convention with a positive self-energy in the inverse propagator, in Feynman–'t Hooft gauge. The singlet parameters do not enter this block.

The [pre-execution brief](task-brief.md) and [standards task record](../../../research/standards-provenance/tasks/20261123-top-charged-weak.json) bind the actual lock on xSM main 5d354669. The same numerical inputs were checked against the [historical Draft #1 massless-vector input file](https://github.com/lm846802297-bit/xsm-electroweak-transition/blob/b1ba61e6a27224696822db20367da2e364e9e33b/docs/theory/top-massless-vector-audit/inputs.json) by exact SHA-256. The old result is not rerun or added here. The present local checkout uses the lock-enabled baseline, not the older scientific branch which does not yet carry that lock. No changes to the lock or any earlier archive.

## Independent projections and physical cut

Use a massless bottom, unit charged-current CKM weight and equal W and charged Goldstone masses. The two terms from (A61) are:

$$
(\Sigma_s+\Sigma_v)_{\rm ch}
=-\frac{g_0^2}{4}(D-2) B_1(k;0,M_W)
 -\frac{y_{t,0}^2}{2} B_1(k;0,M_W),\qquad
D=4-2\epsilon,\quad k^2=M_t^2.
$$

The left-handed W vertices contract the shifted bottom momentum with a Dirac factor of magnitude D minus 2. At the charged Goldstone vertices only the top Yukawa survives. On external on-shell spinors each chiral momentum projector contributes half the vector numerator, with no bottom chirality-flip scalar term. Let the Feynman parameter weight the massive W or Goldstone propagator; the selected numerator after shifting is proportional to that parameter. At the pole define:

$$
b=\frac{M_W^2}{M_t^2},\quad z=1-b,\quad L=\ln\frac{\bar\mu^2}{M_t^2},\quad
\frac{\Delta}{M_t^2}=x(x-z),\qquad
I=\int_0^1 dx\,x\left[L-\ln(x(x-z)-i0)\right].
$$

The shifted spinor projection integrates this expression directly, splitting at the interior branch point. In units of the common loop factor it gives the finite W and Goldstone projections below; the constant in the W term is the product of its evanescent numerator and the divergent moment:

$$
16\pi^2\Sigma_W=\frac{g_0^2}{4}(2I-1),\qquad
16\pi^2\Sigma_G=\frac{y_{t,0}^2}{2}I.
$$

The separate PV route evaluates the scalar moment analytically and reduces the vector moment using (A42); these formulas strip the loop factor:

$$
B_{0,\mathrm{fin}}^{\rm red}=L+2-z\ln|z|-b\ln b+i\pi\max(z,0),
\qquad
B_{1,\mathrm{fin}}^{\rm red}=-\frac{zB_{0,\mathrm{fin}}^{\rm red}+b(1+L-\ln b)}2 .
$$

The vector pole coefficient is minus one-half. The evanescent W correction and the exact threshold and branch controls are:

$$
\Sigma_{W,\rm ev}=-\frac{g_0^2}{64\pi^2},\qquad
\left.B_{1,\mathrm{fin}}^{\rm red}\right|_{b=1}=-\frac{L+1}{2},\qquad
\operatorname{Im} B_{1,\mathrm{fin}}^{\rm red}=-\frac{\pi}{2}\max(z,0)^2 .
$$

The imaginary part is continuous at threshold. No mass is taken to zero except the bottom's.

The absorptive part of the **sum**, not either gauge-fixed term separately, reproduces the width obtained by contracting the tree decay with the physical W polarization sum. The two-body phase space gives:

$$
\Gamma_{bW}=2M_t\,\mathrm{Im}(\Sigma_W+\Sigma_G)
=\frac{g_0^2 M_t}{64\pi}(1-b)^2(2+1/b).
$$

This is a tree-level normalization test, not a precision top-width prediction. It relies on the one-loop pole convention and uses the physical two-body phase-space spin sum independently of the PV route.

## Benchmark and scope

| Scale (GeV) | Relative top-Yukawa-squared shift from this block |
|---:|---:|
| 91.1876 | +0.00421722046158799 |
| 80.379 | +0.00308996637239025 |
| 172.76 | +0.00992634930441539 |

At the Z scale the additive squared-Yukawa shift is +0.00415238170934985, and the absorptive check is 1.488515292713186 GeV. Setting D to four too early raises the relative Z-scale result to +0.00556696022736549. Only the real on-shell top self-energy term in (A30) is projected. The other top self-energies, all 1PR tadpoles, and the W and coupling terms in (A30) are not included here. No full EW conversion, thermal matching, bounce or gravitational-wave prediction is made.

## Replay and record

Run from this directory with Python 3.12.14 and mpmath 1.3.0 (65 working digits):

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_charged_top
    python charged_top.py replay.json
    cmp native-results.json replay.json

The optional CLI comparison uses the immutable input JSON when available; the main numerical replay needs only the local input JSON. The [run manifest](run_manifest.json) records the exact source/input hashes, environment, and signed native output. The run was performed in a managed Linux shell; no CI run was requested or executed.
