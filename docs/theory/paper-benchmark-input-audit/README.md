# xSM: source-version and scalar-tensor input audit at the 350 GeV benchmark

This is a **tree input and normalization audit**, not a thermal or collider prediction. It uses the published [Ramsey-Musolf, Tenkanen and Tran, JHEP 09 (2026) 260, §4.1](https://link.springer.com/content/pdf/10.1007/JHEP09(2026)260.pdf#page=15) point: $m_{h_2}=350$ GeV, $\sin\theta=+0.1$, $a_2=3$, $b_3=40$ GeV, $b_4=0.3$. For the tree map use $m_{h_1}=125.1$ GeV and $v_h=246$ GeV (the paper's approximate electroweak VEV, taken exactly as our reference). These labels do not provide the authors' one-loop renormalized 4D inputs. The paper explicitly upgrades the tree relations at zero temperature before its thermodynamics calculation.

Baseline: [xSM Draft #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1) at 6621e5f0ced37458fd5918936bfa1781596cf44f, main 0e421748932d8b0a0ec635bc909a2f9be06f62c6. This does not repeat the [positive thermal scalar mass audit](../positive-mass-scalar-tadpole-audit/README.md), [screened zero mode audit](../screened-scalar-zero-mode-audit/README.md), or temporal-gauge thresholds.

## Journal vs. preprint: a signed and ordered matrix check

Use the canonical real doublet coordinates $r_i$ in $H^\dagger H=\sum_{i=1}^4r_i^2/2$, and let $(r_4,S)=(v_h,0)$ in the zero-temperature vacuum. Write $h=r_4-v_h$. In Eq. (2.1), $c_H=\mu^2$ is the signed coefficient of $H^\dagger H$. Stationarity requires $c_H=-\lambda v_h^2$ and $b_1=-a_1v_h^2/4$. The neutral Hessian in the $(h,S)$ basis is
$$
M_0^2=\begin{pmatrix}
2\lambda v_h^2 & a_1v_h/2\\
a_1v_h/2 & b_2+a_2v_h^2/2
\end{pmatrix}.
$$
The displayed Eq. (2.3) rotation is the same in both sources: $(h_1,h_2)^t=R(h,S)^t$ with first row $(c_\theta,s_\theta)$ and second row $(-s_\theta,c_\theta)$. Independent symbolic differentiation and rotation gives $M_0^2=R^t\,\mathrm{diag}(m_{h_1}^2,m_{h_2}^2)\,R$, agreeing with the [journal Eqs. (2.4)–(2.8)](https://link.springer.com/content/pdf/10.1007/JHEP09(2026)260.pdf#page=6).

The [literal arXiv:2409.17554v1 Eqs. (2.4), (2.5), (2.7)](https://arxiv.org/pdf/2409.17554v1#page=6) instead interchange the two diagonal mass assignments and reverse the signed off-diagonal entry, compared with the journal. With **the same** Eq. (2.3) rotation the resulting Hessian is diagonal with the masses in reverse order. Merely sorting eigenvalues would still produce $125.1,350$ GeV and would miss the error: the physically light state would have 1% rather than 99% doublet content. This is not fixed by just reflecting $S$; that leaves the two original diagonal entries unchanged and would require transforming the other odd-singlet couplings. We do not infer which input was used in the authors' calculation.

| Tree input | Journal, consistent with Eq. (2.3) | Literal arXiv v1 formulas |
|---|---:|---:|
| $c_H$ (GeV²) | −8,359.25495 | −60,715.75005 |
| $\lambda$ | 0.13813297 | 1.00330078 |
| $a_1$ (GeV) | −86.43446968 | +86.43446968 |
| $b_1$ (GeV³) | +1,307,667.09177 | −1,307,667.09177 |
| $b_2$ (GeV²) | +30,657.50010 | −74,055.49010 |

As an additional signed check the canonical potential's third derivative $\partial^3V/(\partial h_1^2\partial h_2)$ equals +84.28686 GeV at this tree point in the specified eigenstate phase convention, and reproduces journal Eq. (2.13). We do not compute a width or experimental limit.

## Scalar tensor contract

The [pinned DRalgo real-singlet example](https://github.com/DR-algo/DRalgo/blob/e00c28e3d01c83f2f2f205ba508b21709eef2d41/examples/15_smRealSinglet.m#L55-L109) constructs the mass, cubic, quartic and tadpole potentials before calling GradMass, GradCubic, GradQuartic and GradTadpole. Its source was read back at Git blob 49a6150ac7d5a13807668e27b09ab9267b788cad and SHA-256 9a5029e3594a60cd559a063419285f35e74bb445f5961360eaee47acec10d3e2. In the example's potential coefficient names, the tree inputs are $m1=c_H$, $\mu_\Sigma=b_2$, $\lambda_{1H}=\lambda$, $\lambda_\Sigma=b_4$, $\lambda_m=a_2$, $\mu_m=a_1$, $c3=b_3$, $\mu_1=b_1$, when the generated quadratic invariants equal $H^\dagger H$ and $S^2$, respectively. Do not replace $\mu_m$ by $a_1/2$ again at the potential coefficient input: the example already inserts the factor in the cubic potential.

The separately derived symmetric derivatives at the unshifted field origin are:

| Canonical index pattern | Derivative |
|---|---:|
| $S$ | $b_1$ |
| $ii$, $SS$ | $c_H$, $b_2$ |
| $iiS$, $SSS$ | $a_1/2$, $2b_3$ |
| $iiii$, $iijj$ ($i\ne j$ in the doublet) | $6\lambda$, $2\lambda$ |
| $iiSS$, $SSSS$ | $a_2$, $6b_4$ |
| $ijS$ with $i\ne j$ | 0 |

These are positive-potential derivatives, not Minkowski Feynman-rule factors; the signs of their coefficients can be negative. Index patterns are given in a canonical real basis and do not assume a specific GroupMath component order. Once a Wolfram runtime is available, verify the generated invariants and component ordering explicitly and compare the `Grad*` outputs to this table **before** matching. No Wolfram, GroupMath or DRalgo code was run.

For the engineering **tree-only** 3D comparison, $F_4=\sqrt T F_3$, $V_3(F_3)=V_4(\sqrt T F_3)/T$, so the mass-squared coefficients are unchanged, quartics multiply by $T$, cubics by $\sqrt T$, and the tadpole divides by $\sqrt T$. At the paper's $82.12$ GeV temperature label, the engineering $a_{1,3}$ is −783.26991 GeV$^{3/2}$, while the corresponding $iiS$ derivative is −391.63495 GeV$^{3/2}$. The JSON gives both $82.12$ and $64.75$ GeV as normalization points only. These are **not** matched 3D couplings or thermal curvatures, and carry no hard-mode thermal shifts or finite renormalization. The needed next input step is the one-loop pole-to-running conversion, including the 1PR tadpole convention, of [Niemi, Schicho and Tenkanen, App. A.2](https://arxiv.org/pdf/2103.07467#page=17), checking the [journal erratum](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.109.039902) and fixing the common running scale and SM inputs before comparing actual matching.

## Bare positive-mass domain

At the symmetric doublet background $H=0$, the reconstructed **bare 4D tree** doublet curvature is $c_H+a_1s/2+a_2s^2/2$. It is negative at $s=0$ and strictly positive only for $s<-61.62295$ GeV or $s>90.43444$ GeV. The bare singlet curvature has a positive minimum $28,879.72$ GeV² for all real $s$. The controlled positive-mass thermal kernels from prior audits therefore cannot be blindly evaluated using the bare $H=0,s=0$ inputs of this real point. This does not establish where the authors' thermal minima are, nor what sign their screened or matched curvatures have.

## Local replay

Only the bundled scalar tree inputs are needed. Python 3.12.14, SymPy 1.14.0, GPT-managed Linux.

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_benchmark
    python benchmark.py /tmp/xsm-tree-input-replay.json
    cmp native-results.json /tmp/xsm-tree-input-replay.json

Exact byte replay assumes the same environment; no CI was dispatched for this audit.
