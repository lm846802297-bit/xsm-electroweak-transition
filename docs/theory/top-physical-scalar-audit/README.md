# Top pole matching: physical h1 and h2 exchange

This audit evaluates the two **physical CP-even scalar 1PI diagrams** in the top self-energy in [Niemi, Schicho and Tenkanen (NST), Phys. Rev. D 103, 115035](https://doi.org/10.1103/PhysRevD.103.115035), [arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3#page=22), (A20)–(A23), (A26), (A30), (A36)–(A39), (A42), and (A61). The [published erratum](https://doi.org/10.1103/PhysRevD.109.039902) does not change these terms. This does not include Goldstones, 1PR tadpoles, the other top diagrams, or the W and coupling contributions to the top-Yukawa conversion.

The benchmark uses the published [Ramsey-Musolf, Tenkanen and Tran, JHEP 09 (2026) 260, §4.1](https://link.springer.com/content/pdf/10.1007/JHEP09(2026)260.pdf#page=15), with heavy mass 350 GeV and positive sine of the mixing angle 0.1. Its previously audited [input bytes](https://github.com/lm846802297-bit/xsm-electroweak-transition/blob/6ff2c3e5af9550359692cb90b77bb84a50b1c15f/docs/theory/paper-benchmark-input-audit/inputs.json) are hash-pinned and were checked. The light mass 125.10 GeV and the electroweak pole inputs follow NST (A20); this is not a re-extraction of the journal authors' running parameters. Only the squared doublet overlaps enter, so changing the sign convention of the mixing angle leaves this particular sub-block unchanged. The singlet-potential coefficients do not enter these two one-loop diagrams.

The local checkout starts at the lock-enabled [exact predecessor](https://github.com/lm846802297-bit/xsm-electroweak-transition/commit/b22fbc1c221720982529cb260f54385ea95a27fa). See the [pre-execution brief](task-brief.md) and [standards task record](../../../research/standards-provenance/tasks/20261123-top-physical-scalar.json).

## Two projections

Use the NST inverse fermion propagator with positive self-energy, MS-bar subtraction, the Feynman prescription and real top pole momentum. For the physical scalar exchange the normalized scalar Yukawa factors and the squared-mixing relation are:

$$
y_{t,0}^2=\frac{g_0^2 M_t^2}{2M_W^2},\quad
g_0^2=\frac{8M_W^2G_\mu}{\sqrt2},\qquad
w_1=\cos^2\theta=0.99,\quad w_2=\sin^2\theta=0.01 .
$$

The internal top is the first propagator, and x weights the scalar. After shifting the internal top momentum, its numerator between external spinors contains the mass term plus x times the external slash; the physical scalar Yukawa vertices make their signs add. The direct finite quadrature for either scalar is therefore:

$$
b=\frac{M_h^2}{M_t^2},\quad L=\ln\frac{\bar\mu^2}{M_t^2},\qquad
K_h=\int_0^1dx\,(1+x)\left[L-\ln\!\big((1-x)^2+bx-i0\big)\right].
$$

The scalar and vector loop moments are instead computed in the second route with the **exact hash-pinned analytic quadratic PV primitive and A42 vector reduction** from the preceding neutral audit, without using its spinor trace or prior numerical outputs. With B1 defined by (A39), B1 carries minus x, giving B0 minus B1. There is no evanescent numerator in this scalar block; its stripped pole coefficient is the integral of 1+x, equal to 3/2 for either mass. The relative shift inserted in (A30) is:

$$
\left.\frac{\delta y_t^2}{y_{t,0}^2}\right|_{h_1,h_2}
=2\operatorname{Re}\sum_i\Sigma_{h_i}
=\frac{y_{t,0}^2}{16\pi^2}\sum_i w_i\operatorname{Re}K_{h_i}.
$$

The stripped signed kernel routes agree within 2.38 × 10⁻⁶⁶ at all three input scales. Both physical loops are real because the external top pole is below the top-plus-scalar cut. An auxiliary equal-internal-mass control checks the exact threshold and the causal imaginary part above it with the on-shell numerator held fixed; it is not a prediction of an off-shell self-energy.

## Benchmark result

| Scale (GeV) | h1 unweighted | h2 unweighted | Relative squared-Yukawa shift, weighted 1PI sum |
|---:|---:|---:|---:|
| 91.1876 | −0.006433074767 | −0.020084939812 | −0.006569593417 |
| 80.379 | −0.008793089634 | −0.022444954680 | −0.008929608285 |
| 172.76 | +0.005519535631 | −0.008132329414 | +0.005383016981 |

At the Z scale, the additive shift of the squared top Yukawa is −0.006468587496. The change from zero mixing in this sub-block is −0.000136518650 relative to the tree squared Yukawa at all three tested scales; the common scale logarithm cancels out of that difference. At zero mixing only the light-state contribution remains; when the two scalar masses coincide the result is independent of their mixing angle. The code checks these limits separately.

## Heavy mass is not by itself a decoupling test

For fixed light mass, weak inputs and renormalization scale, expand only the heavy-state kernel. The remainder can be bounded without relying on a numerical fit:

$$
A(b)=\frac32(L-\ln b)+\frac54,\qquad
0\le A(b)-K_h\le 2\left[\left(1+\frac1b\right)\ln\left(1+\frac1b\right)+\frac{\ln b}{b}\right].
$$

Indeed the difference is the integral of the positive log of one plus the ratio of the squared remainder in the Feynman denominator to bx. Bounding the weight by two and that squared remainder by one gives the stated integral bound. With fixed nonzero mixing the difference from the light-only result grows logarithmically in magnitude. This artificial fixed-mixing family at arbitrarily large mass is not a claim about perturbative scalar-potential trajectories. If the mixing **also** decreases, for example as the inverse heavy mass, the squared mixing suppresses the log and the contribution tends to the light-only result. In the recorded scaling diagnostic the heavy mass of 35 TeV gives −9.70765 × 10⁻⁸ instead of the fixed-mixing −9.70765 × 10⁻⁴. This is decoupling of the selected 1PI difference along the stated trajectory, not a proof of full-model low-energy decoupling.

## Reproduce

From this directory, Python 3.12.14 with mpmath 1.3.0:

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_physical_top
    python physical_top.py replay.json
    cmp native-results.json replay.json

The script checks the sibling PV source/input hashes before reuse. The optional historical fixture check requires only the pinned input JSON from the archived benchmark. The [run manifest](run_manifest.json) records the actual inputs, source, environment, signed native output and validation. No CI, full gauge or Yukawa conversion, dimensional reduction, phase transition or gravitational-wave prediction was performed.
