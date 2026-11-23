# Top pole matching: Z and neutral Goldstone exchange

This audit computes the Z and neutral Goldstone top self-energy diagrams in [Niemi, Schicho and Tenkanen (NST)](https://doi.org/10.1103/PhysRevD.103.115035), [arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3#page=22), especially (A20)–(A23), (A26), (A30), (A36)–(A39), (A42) and (A61). The [published erratum](https://doi.org/10.1103/PhysRevD.109.039902) does not change these two terms. The calculation is at zero temperature in NST's Feynman–'t Hooft gauge and MS-bar scheme, with the self-energy entering the inverse propagator with a plus sign. This gauge-fixed subset is not asserted to be independently gauge invariant.

The input file is pinned to the actual [charged weak branch baseline](https://github.com/lm846802297-bit/xsm-electroweak-transition/commit/8f368bf32a9fc7ebbd43f64adfb9b576baf363e1). The prior numerical calculation is not rerun or included. The [pre-execution task brief](task-brief.md) and [task record](../../../research/standards-provenance/tasks/20261123-top-neutral-weak.json) identify the exact repository standards lock and gate criteria.

## Independent signed routes

Write the Z vertex with tree-level chiral couplings evaluated from the Fermi constant and W and Z pole masses:

$$
g_0^2=\frac{8M_W^2G_\mu}{\sqrt2},\quad s_0^2=1-\frac{M_W^2}{M_Z^2},\quad
g_{Z,0}^2=\frac{g_0^2}{1-s_0^2},\quad
g_L=g_{Z,0}\!\left(\frac12-\frac23s_0^2\right),\quad
g_R=-g_{Z,0}\frac23s_0^2 .
$$

Use the propagator order massive top then massive vector. Let the Feynman parameter x weight the latter. At the real on-shell momentum, the shifted fermion numerator projects onto x times the external momentum plus the mass term. The gamma contraction separates the mass-flip product of left and right couplings from the chirality-preserving term. For the neutral Goldstone, the relative phase of the pseudoscalar vertices and the anticommuting gamma-five reverse the mass numerator while preserving the slash term. With the NST overall signs, the two stripped, projected numerators are:

$$
N_Z(D,x)=-D g_Lg_R+(D-2)\frac{g_L^2+g_R^2}{2}x,\qquad
N_{G^0}(x)=\frac{y_{t,0}^2}{2}(x-1).
$$

For squared invariant ratio r, squared mass ratio b and logarithm L, the integrals use:

$$
r=\frac{k^2}{M_t^2},\qquad b=\frac{M_Z^2}{M_t^2},\qquad
L=\ln\frac{\bar\mu^2}{M_t^2},\qquad
\frac{\Delta}{M_t^2}=1+(b-1-r)x+rx^2,\qquad
K(x)=L-\ln(\Delta/M_t^2-i0).
$$

The direct route integrates the four-dimensional numerators against the causal log and adds the evanescent numerator times the MS-bar pole **before** subtraction. The finite Z term generated this way is:

$$
\Sigma_{Z,\rm ev}=\frac{2g_Lg_R-(g_L^2+g_R^2)/2}{16\pi^2}.
$$

The second route evaluates the scalar logarithm analytically using a quadratic primitive, then reduces the vector moment by (A42). It uses the two **printed** electroweak polynomials in (A61) as coefficients, independently confirming that they are the combinations below. In this route the pole coefficients are +1 and −1/2 for the scalar and vector moments, in units of the common loop factor:

$$
\Sigma_Z=-Dg_Lg_R B_0-(D-2)\frac{g_L^2+g_R^2}{2}B_1,\qquad
\Sigma_{G^0}=-\frac{y_{t,0}^2}{2}(B_0+B_1).
$$

In particular the minus sign is outside **both** neutral-Goldstone PV terms. The two paths agree at all three scales within 2.90 × 10⁻⁶⁹ for each dimensionless self-energy. The benchmark pole has r=1 and lies below the top-plus-Z cut, so it is real. A separate diagnostic at equal internal masses checks the vanishing imaginary part at the exact threshold, the positive imaginary scalar moment above threshold, and a pseudothreshold with its repeated root outside the integration domain. This diagnostic varies the scalar integral invariant with the on-shell numerator held fixed; it is not an off-shell self-energy or a prediction of a decay width.

## Result and limits

| Renormalization scale (GeV) | Relative squared-Yukawa shift, neutral block only |
|---:|---:|
| 91.1876 | +0.00016375310025433245 |
| 80.379 | +0.00045650594181043973 |
| 172.76 | −0.00131893275341788782 |

At the Z scale the individual real self-energies are −0.000971336056312390 (Z) and +0.001053212606439556 (neutral Goldstone). Taking twice the real sum gives the relative correction in the table, or +0.000161235435647436 additively to the squared Yukawa. The small net value reflects a cancellation at this scale. Setting the vector numerator dimension to four too early would instead yield +0.001395274139279575 relative. This subset excludes the charged and massless-vector top terms, physical Higgs loops, tadpoles and the other (A30) W/coupling corrections. It does not determine the full, gauge-complete top Yukawa or any thermal or transition result.

## Reproduction

Python 3.12.14, mpmath 1.3.0; run from this directory. The script checks the prior input fingerprint in the sibling audit before using the matching numerical inputs:

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_neutral_top
    python neutral_top.py replay.json
    cmp native-results.json replay.json

The [run manifest](run_manifest.json) links the source and input hashes, environment, signed native results and local validation. No CI execution was requested or run.
