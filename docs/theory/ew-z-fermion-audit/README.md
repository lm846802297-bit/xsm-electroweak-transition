# xSM vacuum input: closed fermion loops in the Z pole self-energy

This audit computes the one-loop closed-fermion 1PI block of the transverse Z self-energy in the scheme of Niemi, Schicho and Tenkanen (NST). It advances the Z input after the [charged-current W top–bottom block](../ew-w-top-bottom-audit/README.md), without recalculating that block or the [vacuum scalar block](../vacuum-scalar-pole-audit/README.md). The pinned [xSM draft #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1) head is **29ba0e97b845d9c51913663b217639112ed95086**. Each fermion group is reported separately; NST's exact massless approximation for every fermion except top is used.

## Source and scheme

Use [NST, Phys. Rev. D 103, 115035 (2021)](https://doi.org/10.1103/PhysRevD.103.115035), specifically [arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3#page=17), (A20), (A21), (A26), (A29), (A36)–(A44), (A59) and (A60), and the [published erratum](https://doi.org/10.1103/PhysRevD.109.039902). Both exact PDF hashes are in the inputs. The erratum notes an old numerical implementation error in (A28) and does not modify (A29), (A59) or (A60). We use NST's inverse-propagator convention with positive self-energy, dimensional regularization and the finite part after subtracting the displayed MS-bar pole. Within the one-loop block all couplings are tree-level expressions from the same pole-mass and Fermi-constant inputs:

$$
g_0^2={8 M_W^2 G_\mu\over\sqrt2},\qquad
s_{W,0}^2=1-{M_W^2\over M_Z^2},\qquad
g_{Z,0}^2={g_0^2\over1-s_{W,0}^2},\qquad
g_0'^2={g_0^2(M_Z^2-M_W^2)\over M_W^2}.
$$

## Independent loop routes

Define the vertex as follows, where the color and flavor multiplicities are applied once outside the trace:

$$
\frac{g_{Z,0}}2\gamma^\mu(v_f-a_f^{\rm ax}\gamma_5),\qquad
v_f=T_3^f-2Q_f s_{W,0}^2,\qquad a_f^{\rm ax}=T_3^f .
$$

The NST coefficient in (A60) is distinct from the axial charge:

$$
a_f^{\rm NST}=4Q_f s_{W,0}^2(2Q_f s_{W,0}^2-2T_3^f),\qquad
v_f^2+(a_f^{\rm ax})^2={1+a_f^{\rm NST}\over2},\qquad
v_f^2-(a_f^{\rm ax})^2={a_f^{\rm NST}\over2}.
$$

For the first route, take the even part of the transverse tensor trace. The terms with one axial vertex do not survive the symmetric two-point integration. The two-axial trace can use anticommuting gamma-five in this non-anomalous two-point block; in the momentum numerator it gives the same sign as the two-vector trace and flips the fermion mass numerator. After the Feynman shift set:

$$
s=k^2,\qquad u=x(1-x),\qquad \Delta=m_f^2-us,\qquad D=4-2\epsilon .
$$

Dropping external momentum tensors, the trace coefficient after stripping the overall couplings is proportional to:

$$
[v_f^2+(a_f^{\rm ax})^2]\left[(2/D-1)\ell^2+us\right]
+[v_f^2-(a_f^{\rm ax})^2]m_f^2 .
$$

The regulated tadpole identity under the shifted integral is:

$$
\int\frac{\ell^2}{(\ell^2-\Delta+i0)^2}=\frac D2 A(\sqrt{\Delta}),\qquad
(2/D-1)\frac D2=1-\frac D2.
$$

Combining the pole and finite parts before subtraction therefore gives:

$$
\Pi_{Z,f}^{\rm fin}(s)=\frac{N_c^f g_{Z,0}^2}{16\pi^2}
\int_0^1\!dx\;2\left[v_f^2us-(a_f^{\rm ax})^2\Delta\right]
\ln\!\left(\frac{\bar\mu^2}{\Delta-i0}\right).
$$

The coefficient of the MS-bar pole, in the same normalization, is the integral of the polynomial multiplying the log:

$$
\frac{N_c^f g_{Z,0}^2}{16\pi^2}\left[
 (v_f^2+(a_f^{\rm ax})^2)\frac{s}{3}-2(a_f^{\rm ax})^2m_f^2\right].
$$

The second route uses (A59)/(A60) and the (A41)–(A43) tensor reduction with the equal-mass identities:

$$
B_1=-\frac12 B_0,\qquad
2B_{00}=\frac{(2m_f^2-s/2)B_0+A(m_f)}{D-1},\qquad
K=A(m_f)+\left[m_f^2-\frac{s}{2}-\frac{a_f^{\rm NST}}{1+a_f^{\rm NST}}m_f^2\right] B_0 .
$$

The code evaluates this expression without dividing by the chiral weight and uses the analytic equal-mass finite B0 rather than reusing the direct parameter integral. The finite term from the pole numerator times the order-epsilon term in the tensor denominator is retained. The two paths agree within **9.2 × 10⁻⁶⁴ GeV²** over both momenta, all three scales and all fermion groups. Checks include the vector-current Ward condition (vanishing vector part at zero momentum), the axial zero-momentum term via (A44), the log-scale derivative, and the analytic massless limit. The below-threshold top loop is evaluated at zero and at the physical Z pole. The massless loop is scaleless at zero and has an absorptive part at the pole.

## Pole values and A29

At the scale 91.1876 GeV the finite block is:

| Contribution | Re at zero (GeV²) | Re at Z pole (GeV²) | Im at Z pole (GeV²) |
|---|---:|---:|---:|
| Top | +198.7761754661 | +181.0812540104 | 0 |
| All others, massless | 0 | +118.2216892719 | +222.8426343070 |
| All closed fermion loops | +198.7761754661 | +299.3029432823 | +222.8426343070 |

For the massless groups, analytic continuation gives a positive imaginary part in the NST sign convention. Their imaginary parts divided by the Z pole mass reproduce each tree-level partial width calculated from the respective chiral vertex; their sum is **2.4437822062 GeV**. This is a normalization and branch-cut check, not a prediction of the measured total width including electroweak, mass and QCD corrections. Setting the tensor denominator to its four-dimensional value before subtraction would instead give **+265.3865793572 GeV²** for the combined real part, missing **+33.9163639252 GeV²**.

Only the real part enters the isolated Z-input term of (A29):

$$
\left.\frac{\delta g'^2}{g_0'^2}\right|_{Z,f}=
\frac{\operatorname{Re}\Pi_{Z,f}(M_Z^2)}{M_Z^2-M_W^2}.
$$

This is **+0.09764978775** for top alone (additive change to the squared coupling **+0.01194780326**) and **+0.16140195761** for the full closed fermion Z block (additive **+0.01974811088**). Neither is the physical total shift. The W and coupling-shift terms in (A29) are not included here, nor are the other Z 1PI loops or the 1PR Higgs tadpoles (even if a fermion generates their subloop). The masses and weak gauge inputs are shared with the 350 GeV xSM benchmark, but this pure gauge–fermion block does not depend on its portal or singlet parameters. No full electroweak match, DRalgo comparison, thermal prediction, bounce or gravitational-wave calculation is performed.

## Replay

Run from this directory with Python 3.12.14 and mpmath 1.3.0 (65 working digits):

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_zferm
    python zferm.py replay.json
    cmp native-results.json replay.json

The input JSON is the only numerical input; the paper files are not redistributed. Exact bytes are expected with the versions specified. No GitHub Actions execution was used.
