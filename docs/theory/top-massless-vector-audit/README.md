# xSM vacuum input: massless gluon and photon terms of the top pole

This audit evaluates the two massless-vector exchange diagrams in the top self-energy in Appendix A of [Niemi, Schicho and Tenkanen (NST)](https://doi.org/10.1103/PhysRevD.103.115035). It follows the [closed-fermion W/Z input audit](../ew-fermion-gauge-matching-audit/) on [Draft #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1), exact head 13223296900022a9c04b4dad977b14888a97b0bf. The archived W/Z input file is hash-pinned and was checked byte-for-byte against the repository; no prior loop is rerun or added to the new numerical block.

Use [NST, arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3#page=17), (A20)–(A23), (A26), (A30), (A36)–(A37), (A39), and the last two terms of (A61), with the positive self-energy sign in the inverse fermion propagator. Appendix A.2 explicitly fixes the squared QCD coupling to 1.48409 and neglects its running. The [published erratum](https://doi.org/10.1103/PhysRevD.109.039902) does not change this prescription or the two selected terms. The fixed value is a reproduction input, not a claim about the independently evolved QCD coupling at each scale. The electric charge instead comes from the same tree Fermi-constant and pole inputs used previously, not a separately chosen fine-structure constant:

$$
g_0^2=\frac{8M_W^2G_\mu}{\sqrt2},\qquad
e_0^2=g_0^2\left(1-\frac{M_W^2}{M_Z^2}\right),\qquad
y_{t,0}^2=\frac{g_0^2 M_t^2}{2M_W^2},\qquad Q_t=\frac23,\quad C_F=\frac43 .
$$

## One-loop pole projection

The selected scalar-plus-vector combination, as printed in (A61), is

$$
(\Sigma_s+\Sigma_v)_{g+\gamma}
=-\left(C_Fg_s^2+Q_t^2 e_0^2\right)
   [D\,B_0+(D-2)B_1]_{(k^2=M_t^2;\,M_t,0)} .
$$

Let x weight the massless propagator in the Feynman representation of (A37). At the pole the shifted invariant is: 

$$
\Delta(x)=(1-x)^2M_t^2.
$$

The endpoint logarithm is integrable: the pole is at this cut's endpoint, and there is no absorptive term in this particular pole projection. The code also verifies a real approach to the endpoint, varying only this denominator, not the on-shell projector. No claim is made about the on-shell wave-function residue or the weak decay width.

Let the squared-scale logarithm be defined by

$$
L=\ln\frac{\bar\mu^2}{M_t^2},\qquad D=4-2\epsilon .
$$

The PV route integrates the scalar and vector moments separately, using the vector decomposition in (A39). With the overall loop factor stripped, their values before subtraction are

$$
16\pi^2 B_0=\frac1\epsilon+L+2,\qquad
16\pi^2 B_1=-\frac1{2\epsilon}-\frac L2-\frac32 .
$$

The second route uses the external on-shell spinor projection of the shifted numerator. External spinors satisfy the massive Dirac equation, and the term linear in the shifted momentum integrates to zero. The fermion momentum and the dimension-dependent numerator are

$$
q_f=xk-\ell,\qquad D-(D-2)x=(4-2x)-2\epsilon(1-x).
$$

This gives the independent numerical integral

$$
\int_0^1dx\left[
 (4-2x)\big(L-2\ln(1-x)\big)-2(1-x)
 \right]=3L+4 .
$$

Both routes give the following combined pole and finite part:

$$
16\pi^2[D B_0+(D-2)B_1]=\frac3\epsilon+3L+4 .
$$

Setting the numerator dimension to four too early gives a constant 5 instead of 4. The spurious constant is absent when the dimension-dependent numerator multiplies the ultraviolet poles before subtraction. Direct quadrature, PV and the closed expression agree to better than 3.0 × 10⁻⁶⁷ for the stripped finite bracket at the three input scales. This pole projection of the closed massless-vector class also has a longitudinal gauge check: between on-shell spinors each longitudinal-momentum contraction may be replaced, with the corresponding sign, by the internal inverse fermion propagator; after cancellation the remainder is odd in the exchanged momentum against its squared massless denominator and vanishes in dimensional regularization. No massive weak-gauge subset is claimed separately gauge independent by this argument.

## Coefficient and benchmark

At the pole scale the gluon result reproduces the one-loop term of the standard [quark pole-to-MS mass relation](https://arxiv.org/pdf/hep-ph/9911434#page=36) (Chetyrkin and Steinhauser), with the photon term obtained by the Abelian charge replacement:

$$
\left.\frac{\delta m_t}{M_t}\right|_{\bar\mu=M_t,g}
 =-C_F\frac{\alpha_s}{\pi},
\qquad
\left.\frac{\delta m_t}{M_t}\right|_{\bar\mu=M_t,\gamma}
 =-Q_t^2\frac{\alpha_0}{\pi},
\qquad
\alpha_s=\frac{g_s^2}{4\pi},\quad \alpha_0=\frac{e_0^2}{4\pi}.
$$

Only the selected top self-energy terms of (A30) are applied here, at first order in the loop corrections:

$$
\left.\frac{\delta y_t^2}{y_{t,0}^2}\right|_{g+\gamma}
 =2\operatorname{Re}(\Sigma_s+\Sigma_v)_{g+\gamma}
 =-\frac{C_Fg_s^2+Q_t^2e_0^2}{8\pi^2}\,(3L+4).
$$

| Renormalization scale (GeV) | Gluon relative term | Photon relative term | Sum |
|---:|---:|---:|---:|
| 91.1876 | −0.004162571884 | −0.000088881349 | −0.004251453233 |
| 80.379 | +0.014808974904 | +0.000316208754 | +0.015125183658 |
| 172.76 | −0.100246503621 | −0.002140514265 | −0.102387017886 |

The result near the Z scale is small because the explicit logarithm nearly cancels the constant; the selected one-loop sum crosses zero at about 88.69794 GeV. At the Z scale the additive squared-Yukawa shift in this subset is −0.004186088160, whereas the prematurely four-dimensional relative result would be −0.029848207705. These isolated one-loop terms neither predict the full running top mass or Yukawa nor estimate the total error. The other top self-energy terms (massive weak bosons, scalar/Goldstone, tadpoles) and the (A30) W and coupling terms are not part of this audit. No full electroweak match, dimensional reduction, thermal prediction, bounce or gravitational-wave calculation is performed.

## Replay

From the repository's audit directory, Python 3.12.14 and mpmath 1.3.0 (65 working digits):

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_top_vector
    python top_vector.py replay.json
    cmp native-results.json replay.json

The script hash-checks the adjacent archived W/Z input JSON. A stand-alone harness can supply the identical file via the command-line option; no paper PDF is required for replay or redistributed. Exact bytes assume matching software versions. No GitHub Actions execution is claimed.
