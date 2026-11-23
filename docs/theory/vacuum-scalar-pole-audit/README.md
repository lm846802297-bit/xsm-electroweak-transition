# xSM: massive-scalar part of the one-loop pole-to-running input map

**Scope.** On the exact [canonical 350 GeV tree point](../paper-benchmark-input-audit/README.md), this is the zero-temperature, finite $\overline{\rm MS}$ contribution with only the two **physical massive scalars** running in the internal loops of the scalar two-point functions and their 1PR tadpoles. It does not calculate Goldstone, gauge, ghost, fermion or full electroweak input corrections. It does not produce a complete running input tuple, thermal DR matching, or a transition prediction. Starting xSM draft PR #1 head: 6ff2c3e5af9550359692cb90b77bb84a50b1c15f, main 0e421748932d8b0a0ec635bc909a2f9be06f62c6.

## Source pins and conventions

- [Niemi, Schicho and Tenkanen, arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3), Eqs. (4)–(9), Appendix A.2 and A.3, in particular (A21), (A26)–(A38), (A45)–(A51), and (A54)–(A57). The exact downloaded v3 PDF SHA-256 is in `inputs.json`.
- [APS Erratum, Phys. Rev. D 109, 039902](https://doi.org/10.1103/PhysRevD.109.039902), PDF retrieved from the [publisher's full-text endpoint](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevD.109.039902/fulltext), also hashed in `inputs.json`. Its two distinct fixes: the $\sqrt T$ powers in thermal Eqs. (38)/(39), and a sign error in the old **numerical implementation** of (A28). The erratum explicitly says the equation as printed is correct. The corrected temperature power does not independently settle the factor-one-half question in the thermal action. Do not use the old benchmark numerics as validation of (A28).
- The NST Eq. (4) rotation is inverse in angle relative to [RMT journal Eq. (2.3)](https://link.springer.com/content/pdf/10.1007/JHEP09(2026)260.pdf#page=6): to retain the same mass-state phases and $b_3=+40$ GeV use $\sin\theta_{\rm NST}=-0.1$ for the RMT point $\sin\theta_{\rm RMT}=+0.1$. Canonical scalar inputs are $m_1=125.1$, $m_2=350$, $v=246$ GeV, $a_2=3$, $b_3=40$ GeV and $b_4=0.3$. The retained $v=246$ is the earlier tree fixture, not the $G_\mu$-inferred exact tree value. No mixing-angle renormalization condition is introduced beyond the fixed-angle input of NST A.2.

## Finite scalar two-point contribution and 1PR term

Let $g_{ijk}$ and $q_{ijkl}$ be ordinary derivatives of the scalar potential evaluated at the tree vacuum in the ordered mass basis. NST uses signed vertices $c=-g,-q$ in its formulas. The code independently differentiates the reconstructed tree potential and verifies **all seven** relevant massive-scalar couplings against the literal NST Eqs. (A45)–(A51) after converting the angle. With the inverse propagator $p^2-m_i^2(\bar\mu)+\Pi_i(p^2)$, the one-loop relation is $m_i^2(\bar\mu)=M_i^2+\Re\Pi_i(M_i^2)$ to this order.

The scalar bubbles and four-leg seagulls are
$$
\Pi_i^{1\rm PI}(p^2)=\frac12\sum_{ab}g_{iab}^2 B_0^{\rm fin}(p^2;m_a^2,m_b^2)
+\frac12\sum_a q_{iiaa}A^{\rm fin}(m_a^2).
$$
Here $A^{\rm fin}(a)=a[1+\ln(\bar\mu^2/a)]/(16\pi^2)$ and $B_0^{\rm fin}$ is minus the standard Feynman-parameter log of $[ua+(1-u)b-u(1-u)p^2-i0]/\bar\mu^2$, divided by $16\pi^2$. Both explicitly subtract the $1/\epsilon$ pole in NST's convention. Ordered $(a,b)$ runs independently over 1,2 and accounts for both mixed orderings.

Let $t_i^{\rm CW}$ be the positive-potential Coleman–Weinberg derivative (including its algebraic sign), $t_i^{\rm CW}=-\tfrac12\sum_a g_{iaa}A^{\rm fin}(m_a^2)$. At a tree-stationary but not one-loop-stationary expansion point, the induced local vacuum translation is $\delta u_i=-t_i^{\rm CW}/M_i^2$. Thus
$$
\Pi_i^{1\rm PR}=-\sum_k g_{iik}\,\delta u_k
=\sum_k \frac{g_{iik} t_k^{\rm CW}}{M_k^2}.
$$
In NST's notation the scalar part of its signed tadpole is $t_i^{\rm NST}=-t_i^{\rm CW}$, so $c_{iik} t_k^{\rm NST}/M_k^2$ gives the same term. Do not both shift the expansion vacuum by this local translation **and** add the same 1PR term to the inverse propagator. This is not the electroweak-input shift of $v$ extracted from $G_\mu,M_W$. The two definitions must not be substituted for each other. Off-diagonal self-energies do not affect nondegenerate pole *eigenvalues* at first order; they would matter for a different mixing-angle definition or near degeneracy.

## Results and independently tested limits

At NST's reference scale $\bar\mu=M_Z=91.1876$ GeV, the signed pieces (GeV²) are:

| External state at its pole | $\Re\Pi^{1\rm PI}$ | $\Pi^{1\rm PR}$ | $\Re\Pi^{\rm subset}$ |
|---|---:|---:|---:|
| $h_1$ | −6510.24268 | +6314.38490 | −195.85778 |
| $h_2$ | −6574.83540 | +22296.59061 | +15721.75521 |

This highlights the cancellation in the lighter-state contribution and why discarding the signed tadpole is not a safe approximation even in this partial block. For $p^2=M_2^2$ the $h_1h_1$ internal channel is open and its absorptive part is carried in the JSON; we use only the real part in the pole projection and do not calculate a width. Three supplied $\bar\mu$ values illustrate the scale dependence of the **incomplete** subset and do not purport to test complete RG cancellation.

Two evaluations of each $B_0$ (integrable Feynman-parameter quadrature partitioned at real roots, analytic roots/primitive) agree to $4.3\times 10^{-73}$ absolute in these runs. A separate five-point derivative of the matrix Coleman–Weinberg potential verifies both $t^{\rm CW}$ and $-\Pi^{1\rm PI}(0)$ with scaled residuals below $1.4\times10^{-55}$ and $3.7\times10^{-41}$, respectively. These are numerical consistency checks, not theoretical error estimates.

### Mapping just the computed block

For $P_i=\Re\Pi_i(M_i^2)$ of this selected subset, $P_H=c_\theta^2P_1+s_\theta^2P_2$, $P_S=s_\theta^2P_1+c_\theta^2P_2$, and $A=2c_\theta s_\theta(P_1-P_2)/v$ in the RMT angle convention. The partial response from the **explicit scalar mass poles only** in NST (A31)–(A35) is $(\delta c_H,\delta b_2,\delta a_1,\delta b_1,\delta\lambda)=(-P_H/2,P_S,A,-v^2 A/4,P_H/(2v^2))$. The extra coefficient multiplying $\delta v_{\rm EW}/v$ is $(0,-a_2v^2,-a_1,b_1,-2\lambda)$. The JSON records the coefficient but deliberately does **not** fill in $\delta v_{\rm EW}/v$ or add the partial response to claim a full running tuple.

For the erratum-sensitive gauge input, using NST (A20)/(A26) and $\bar\mu=M_Z$, the explicitly known analytic bracket in (A28) contributes **−0.00959357134** to $\delta g^2/g_0^2$, *prior to* its separately required term $[\Re\Pi_W(M_W^2)-\Pi_W(0)]/M_W^2$. This is a sign sentinel, not the full one-loop matching of $g$. Momentum-independent 1PR pieces cancel in that difference, but the separate $\Pi_W(M_W^2)$ used for the VEV relation may not discard them. The next gate is the complete finite $W,Z,t,h_i$ two-point and electroweak input block in the same scheme, including gauge/Goldstone/ghost and fermion terms with the prescribed treatment of poles, and then an independent common-scale check. No Wolfram or DRalgo run was made.

## Replay

GPT-managed Linux; Python 3.12.14, SymPy 1.14.0, mpmath 1.3.0. The input JSON records source hashes, and source documents themselves are not redistributed. Primary replay:

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_pole
    python pole.py /tmp/xsm-vacuum-pole-replay.json
    cmp native-results.json /tmp/xsm-vacuum-pole-replay.json

Exact byte replay assumes the same runtime versions. No CI was dispatched.
