# xSM vacuum inputs: closed fermion loops in the W/Z gauge matching

This audit completes the **closed-fermion 1PI contribution** to the W and Z self-energy terms of the one-loop gauge-coupling conversion in Appendix A of [Niemi, Schicho and Tenkanen (NST)](https://doi.org/10.1103/PhysRevD.103.115035). The new calculation is the five light charged-current families. The massive top–bottom W loop and all neutral-current Z fermion loops are read from immutable earlier audits, not recalculated. A separate algebraic check tracks the shared 1PR Higgs tadpoles, which are not part of the numerical fermion bubble sum.

Base: [xSM draft #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1), exact head 6726ea60b15301e542ff2e18fc6588a8e0d2de73. The formulas are locked to [arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3#page=17), (A20)–(A21), (A26)–(A30), (A36)–(A44), (A58)–(A59), and the [published erratum](https://doi.org/10.1103/PhysRevD.109.039902). The erratum records a sign bug in the old numerical implementation of (A28), not a change to the printed formula. Same MS-bar pole convention, same pole masses and same tree Fermi-constant input as the two archived audits.

## The light W block

NST (A58) takes the top as the only massive internal fermion and sums over six doublet families. After removing the massive top–bottom entry, the multiplicity is two quark families of three colors plus three lepton families. Flavor mixing does not change the inclusive multiplicities in this mass approximation, by CKM unitarity and the vanishing down-type masses. The calculation treats the neutrinos as massless.

$$
n_{\rm light}=2(3)+3=9,\qquad s=k^2,\qquad
g_0^2=\frac{8M_W^2G_\mu}{\sqrt2}.
$$

In one route, shift the Feynman parameter integral for the two-left-current trace in dimensional regularization. The effective mass squared is negative for nonzero timelike momentum, with the indicated branch prescription. At zero the integrals are scaleless. For the timelike case, the finite result is:

$$
u=x(1-x),\qquad
\Pi_{W,\rm light}^{\rm fin}(s)
=\frac{n_{\rm light}g_0^2}{16\pi^2}
 \int_0^1dx\;2us\ln\frac{\bar\mu^2}{-us-i0}
=\frac{n_{\rm light}g_0^2s}{48\pi^2}
 \left[\ln\frac{\bar\mu^2}{s}+\frac53+i\pi\right].
$$

For the second route use the (A58) massless bracket, (A41)–(A43), and the exact shift identity for B1:

$$
B_1=-\frac12 B_0,\qquad 2B_{00}=-\frac{s}{2(D-1)}B_0,\qquad
K=s B_1,\qquad
2B_{00}-K=\frac{s}{2}\!\left(1-\frac1{D-1}\right)B_0.
$$

Subtract the pole only after multiplication. With the common loop denominator stripped from B0, the product is

$$
s\left(\frac13-\frac{\epsilon}{9}+\cdots\right)
 \left[\frac1\epsilon+\ln\frac{\bar\mu^2}{s}+2+i\pi\right].
$$

The two numerical paths agree within 3.1 × 10⁻⁶⁴ GeV². The absorptive part divided by the W pole mass gives 2.0450021317 GeV, coinciding with nine massless tree widths calculated from a single chiral decay vertex. This checks the sign and normalization of the cut, not the physical W total width. Replacing the tensor denominator by its four-dimensional value too early overestimates the real pole contribution by 17.4407531967 GeV².

## Combining the W and Z terms

The dependency hashes in the inputs are for the exact native JSON files from the [W top–bottom audit](../ew-w-top-bottom-audit/) at [29ba0e97](https://github.com/lm846802297-bit/xsm-electroweak-transition/commit/29ba0e97b845d9c51913663b217639112ed95086) and the [Z fermion audit](../ew-z-fermion-audit/) at [6726ea60](https://github.com/lm846802297-bit/xsm-electroweak-transition/commit/6726ea60b15301e542ff2e18fc6588a8e0d2de73). At the reference scale 91.1876 GeV:

| Isolated 1PI quantity | Value |
|---|---:|
| Light W real part at the pole (GeV²) | +100.4063437623 |
| Light W imaginary part at the pole (GeV²) | +164.3752263473 |
| All closed fermion W real part at the pole (GeV²) | +179.4569304622 |
| All closed fermion W real part at zero (GeV²) | +94.0199242500 |
| All closed fermion Z real part at the pole (GeV²) | +299.3029432823 |
| A28 relative squared SU(2) coupling shift | +0.01322393861 |
| A29 relative squared hypercharge coupling shift | +0.05007569924 |
| Additive shift in squared SU(2) coupling | +0.00563716298 |
| Additive shift in squared hypercharge coupling | +0.00612694217 |
| Relative shift in squared VEV from this W input only | +0.01455240254 |

To isolate this class, use its defined transverse self-energies in the self-energy parts of (A28) and (A29) and project the real pole parts. This does not include the explicit non-self-energy term in (A28).

$$
d_g=\frac{\operatorname{Re}\Pi_W(M_W^2)-\Pi_W(0)}{M_W^2},
\qquad
d_{g'}=d_g+
 \frac{\operatorname{Re}\Pi_Z(M_Z^2)-(M_Z^2/M_W^2)\operatorname{Re}\Pi_W(M_W^2)}
 {M_Z^2-M_W^2}.
$$

As an independent completeness check in the full theory, the logarithmic slopes with respect to the log of the squared renormalization scale, and their pole coefficients, reproduce the fermion pieces of the squared-coupling running. With the normalization Q=T3+Y, the three generations give the Weyl-charge sums below (the universal two-thirds is per Weyl fermion).

$$
b_2^{\rm f}=\frac23\,3\frac{3+1}{2}=4,\qquad
b_Y^{\rm f}=\frac23\,3\!\left[
6\!\left(\frac16\right)^2+3\!\left(\frac23\right)^2
+3\!\left(\frac13\right)^2+2\!\left(\frac12\right)^2+1\right]=\frac{20}{3}.
$$

$$
\frac{\partial d_g}{\partial\ln\bar\mu^2}
=\frac{b_2^{\rm f}g_0^2}{16\pi^2},\qquad
\frac{\partial d_{g'}}{\partial\ln\bar\mu^2}
=\frac{b_Y^{\rm f}g_0'^2}{16\pi^2}.
$$

The test obtains these slopes both from changing the scale and directly from the composed ultraviolet pole coefficients. Keeping the massive top in the full-theory loop is consistent with this test; no top threshold decoupling to a lower-flavor effective theory is being asserted.

## Shared 1PR tadpoles

The first term of (A58) and of (A59) contains the same Higgs tadpole combination. Write it as

$$
T_h=\frac{\cos\theta}{m_{h_1}^2}t_{h_1}+
    \frac{\sin\theta}{m_{h_2}^2}t_{h_2},\qquad
\Pi_W^{\rm tad}=-\frac12 vg^2T_h,\qquad
\Pi_Z^{\rm tad}=-\frac12 v(g^2+g'^2)T_h.
$$

At one loop the coupling and mass factors of this already one-loop piece can be evaluated at tree level. The two insertions are momentum independent and share the same fractional mass factor:

$$
\frac{\Pi_W^{\rm tad}}{M_W^2}=\frac{\Pi_Z^{\rm tad}}{M_Z^2}
=-\frac{2T_h}{v}.
$$

Momentum independence cancels this term in (A28); the equal fractional factors cancel it in (A29). Neither identity removes it from the VEV definition. Using NST's positive self-energy pole convention and expanding the defining relation gives

$$
\frac{\delta v^2}{v_0^2}
 =\frac{\operatorname{Re}\Pi_W(M_W^2)}{M_W^2}-d_g ,
 \qquad
\left.\frac{\delta v^2}{v_0^2}\right|_{\rm tad}
=-\frac{2T_h}{v}.
$$

The 1PR shift is in general nonzero in this tree-minimum tadpole convention. It was **not** evaluated numerically or included in the VEV row above. Bosonic (including scalar) W/Z self-energies and the full non-self-energy matching are beyond this 1PI calculation; the independent top self-energy needed for the Yukawa conversion is still missing. This is not a full electroweak conversion, dimensional-reduction run, phase-transition or gravitational-wave result.

## Replay and dependencies

Run from this directory in the repository checkout. The two exact native JSON files are already archived in sibling directories. Their SHA-256 is checked before the calculation; there is no dependency on downloading papers or regenerating previous runs. The CLI also accepts explicit paths to the two pinned JSONs.

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_fermion_gauge
    python fermion_gauge.py replay.json
    cmp native-results.json replay.json

Runs use Python 3.12.14, mpmath 1.3.0, 65 working digits. Exact-byte replay assumes matching software versions. No independent GitHub Actions execution is claimed.
