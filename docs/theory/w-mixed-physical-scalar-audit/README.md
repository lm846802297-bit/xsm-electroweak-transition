# W self-energy: mixed physical Higgs and charged boson bubbles

This audit evaluates only the momentum-dependent physical Higgs–W and physical Higgs–charged-Goldstone bubbles singled out in (A58) of [Niemi, Schicho and Tenkanen (NST)](https://doi.org/10.1103/PhysRevD.103.115035), [arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3#page=21). Other relevant equations are (A20), (A26)–(A30), and (A36)–(A44). The [published erratum](https://doi.org/10.1103/PhysRevD.109.039902) reports a sign error in the historical numerical implementation of (A28), not a change to its printed expression.

Use the [published 350 GeV benchmark](https://link.springer.com/content/pdf/10.1007/JHEP09(2026)260.pdf#page=15), §4.1, with the same NST gauge inputs and squared doublet overlaps as the preceding audited [fixture](../top-shared-1pr-tadpole-audit/). The previous [benchmark input audit](https://github.com/lm846802297-bit/xsm-electroweak-transition/blob/6ff2c3e5af9550359692cb90b77bb84a50b1c15f/docs/theory/paper-benchmark-input-audit/inputs.json) is pinned. The exact baseline is ca8b2ec6654a533939ab775c5f3eae1070cc0c2f. The [pre-execution brief](task-brief.md) and [task record](../../../research/standards-provenance/tasks/20261123-w-mixed-physical-scalar.json) state the scope and actual standards lock.

## Diagram class and signs

In Feynman–'t Hooft gauge the internal W and the charged Goldstone have the same mass. Let q_i be the signed doublet overlap of the physical scalar. Expanding the doublet kinetic term gives the WWh_i vertex with magnitude g_0 M_W q_i and a Lorentz metric, and the Wh_iG vertex with coefficient g_0 q_i/2 multiplying the difference of incoming scalar momenta. Its conjugate supplies the matching charged vertex. The transverse projection removes the external-momentum part of the derivative numerator and leaves four times the loop-momentum tensor. Together with the bosonic propagator signs in the NST inverse-propagator convention, the metric vector bubble enters positively, and the derivative Goldstone bubble enters negatively. These factors reproduce the two selected coefficients in (A58):

$$
\Pi_{W,i}^{\rm mix}(s)=g_0^2 q_i^2\big[ M_W^2 B_0(s;M_{h_i},M_W)-B_{00}(s;M_{h_i},M_W)\big],
\quad g_0^2=\frac{8M_W^2G_\mu}{\sqrt2},\qquad
(q_1^2,q_2^2)=(0.99,0.01).
$$

The scalar seagulls, 1PR tadpoles, and all other W/Z loops are out of scope. In particular a momentum-independent seagull would also affect the W-input part of (A30) even though it cancels from (A28). This selected pair of bubbles is not asserted separately gauge invariant.

## Independent numerical routes

Use the squared masses a and b and invariant s:

$$
a=M_{h_i}^2,\qquad b=M_W^2,\qquad
\Delta=a(1-x)+bx-sx(1-x),\qquad D=4-2\epsilon .
$$

For the direct route the scalar integral is the causal logarithm. For the tensor, D-dimensional angular integration followed by the regulated total derivative gives the factor one-half inside the shifted tadpole integral. Only then is its UV pole subtracted:

$$
\int d^D\ell\,\frac{\ell^2}{D(\ell^2-\Delta+i0)^2}
=\frac12\int d^D\ell\,\frac{1}{\ell^2-\Delta+i0},
$$

$$
16\pi^2 B_0^{\rm fin}=\int_0^1 dx\,\ln\frac{\bar\mu^2}{\Delta-i0},
\qquad
16\pi^2 B_{00}^{\rm fin}=\frac12\int_0^1dx\,\Delta\left[1+\ln\frac{\bar\mu^2}{\Delta-i0}\right].
$$

The second route uses the exact SHA-pinned analytic scalar-PV backend from the [neutral audit](../top-neutral-weak-audit/) and A42 for nonzero momentum, and closed linear-log integrals for both moments at zero. It then applies (A41), with a check using (A44) at zero:

$$
B_{00}=\frac{2a B_0+A(M_W)-(b-a-s)B_1}{2(D-1)},\qquad
B_{00}(0)=\frac{A(M_W)+aB_0(0)}{D}.
$$

For the (A41) denominator retain the finite term from its UV numerator; in stripped units it is:

$$
16\pi^2\Delta B_{00}^{\rm ev}
=\frac{1}{9}\left[\frac32(a+b)-\frac{s}{2}\right].
$$

Direct and reduced bubbles agree within 1.24 × 10⁻⁶⁴ GeV². Both physical W evaluations are real: at the pole the Feynman squared mass is a(1−x)+b x², strictly positive, and the physical cut is at the sum of the internal masses. A diagnostic with equal internal squared masses m² verifies the causal terms below, at and above the cut, in stripped units:

$$
\operatorname{Im}(16\pi^2 B_0)=\pi\beta,\qquad
\operatorname{Im}(16\pi^2 B_{00})=-\frac{\pi s}{12}\beta^3,\qquad
\beta=\sqrt{1-4m^2/s}\quad (s>4m^2).
$$

This auxiliary test is not a W-width calculation.

## W-input projection

Only the selected pair is inserted into the W-input terms:

$$
d_g^{\rm mix}
=\frac{\operatorname{Re}\Pi_W^{\rm mix}(M_W^2)-\Pi_W^{\rm mix}(0)}{M_W^2},
\qquad
\left.\frac{\delta y_t^2}{y_{t,0}^2}\right|_{\rm mix,W\,input}
=d_g^{\rm mix}-\frac{\operatorname{Re}\Pi_W^{\rm mix}(M_W^2)}{M_W^2}
=-\frac{\Pi_W^{\rm mix}(0)}{M_W^2}.
$$

| Scale (GeV) | Pi_W(0) (GeV²) | Re Pi_W(pole) (GeV²) | A28 relative | A30 W-input only |
|---:|---:|---:|---:|---:|
| 91.1876 | −13.805290979 | −12.262052894 | +0.000238862369 | +0.002136782741 |
| 80.379 | −14.258919752 | −13.082419939 | +0.000182098624 | +0.002206995396 |
| 172.76 | −11.507827488 | −8.107194462 | +0.000526349932 | +0.001781181375 |

Differentiating the scale log analytically gives g_0²/(12 × 16π²) for A28. The code checks that as well as the predicted UV/log coefficient separately for each momentum and scalar. The deliberately wrong early-D=4 result uses the same (A41) representation at both momenta: at the Z scale it gives +0.000088891284 instead of +0.000238862369 for A28. Their difference is g_0²/(18 × 16π²), not an attempt to reproduce the historical sign bug.

No top self-energy, common tadpole or prior closed-fermion result is part of the sum; do not use it as a full W input or full Yukawa conversion. The remaining bosonic W/Z sectors and full pole matching are not evaluated. No thermal, bounce or gravitational-wave calculation is performed.

## Replay

From this directory with Python 3.12.14 and mpmath 1.3.0, 65 working digits:

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_mixed_w
    python mixed_w.py replay.json
    cmp native-results.json replay.json

The script validates exact sibling source and input fingerprints. The [run manifest](run_manifest.json) archives typed provenance, signed numerical results and validation. No CI execution was requested or run.
