# xSM zero-temperature inputs: the top–bottom part of the W two-point function

This audit advances only the closed third-generation quark-doublet contribution in the one-loop W self-energy. It complements the [massive-scalar pole block](../vacuum-scalar-pole-audit/README.md). At this order the selected W term depends only on the shared electroweak tree inputs, not the 350 GeV scalar mass or portal parameters. It does not re-evaluate either scalar self-energy. No full electroweak conversion, top-quark self-energy, matching to a thermal EFT, or transition result is implied. Pinned starting [xSM Draft #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1) head 895adc43b8a9a80720452722c04d7f2e538acc1f; main 0e421748932d8b0a0ec635bc909a2f9be06f62c6.

## Source and scheme

Use [Niemi, Schicho and Tenkanen, arXiv:2103.07467v3, Appendix A](https://arxiv.org/pdf/2103.07467v3#page=17), specifically (A20), (A26)–(A30), (A36)–(A44) and the first fermion term of (A58). This is the source's inverse-propagator sign with $+\Pi_W$ and its $\overline{\rm MS}$ loop normalization. At this order use $g_0^2=8M_W^2G_\mu/\sqrt2$ in the loop and the source's $m_b=0$, $V_{tb}=1$, $N_c=3$. $M_t$ and $M_W$ are pole mass inputs, not running masses inserted in the one-loop integrals.

The [APS erratum](https://doi.org/10.1103/PhysRevD.109.039902) states that (A28) is correct as printed; the old **numerical implementation** had an electroweak sign error. Source PDF hashes are in the input file. The prior scalar audit evaluated the explicit non-self-energy bracket of (A28). Here we compute only one new, signed self-energy component, and do not call their sum the full answer.

## Two routes with the same normalization

Let $x$ be the Feynman parameter weighting the massive internal top, $s=k^2$, and $Q(x)=x[M_t^2-(1-x)s]$. The two external momenta used are $s=0$ and $s=M_W^2<M_t^2$, below the $t b$ cut. Then NST's third-generation term is
$$
\Pi_W^{tb}(s)=g_0^2 N_c\,\big[2B_{00}(s;0,M_t)-K(s;0,M_t;1,1,0)\big],\qquad
K=A(M_t)+sB_1 .
$$
Here $A(M)$ denotes NST's one-point integral (the argument in the code is $M^2$). Shifting the loop momentum gives $B_1^{\rm fin}=-\int_0^1 dx\,x\ln(\bar\mu^2/Q)/(16\pi^2)$. Under the regulated integral, Lorentz symmetry and integration by parts give $2B_{00}=\int_0^1 dx\, A(\sqrt{Q(x)})$. This parameter expression is checked against (A41)/(A42), where
$$
2B_{00}=\frac{A(M_t)-(M_t^2-s)B_1}{D-1},\qquad D=4-2\epsilon .
$$
For the latter route keep the finite term from $1/(D-1)=1/3+2\epsilon/9+\cdots$ multiplying the pole in the numerator. Equation (A44) gives a separate check at $s=0$, $2B_{00}(0)=2A(M_t)/D$, also expanded **before** subtraction. The finite result is independent of which exact relation is used. The intentionally incorrect early-$D=4$ comparison uses the *same* (A41) form at both momenta, rather than comparing two differently truncated representations.

The coefficient of $1/\epsilon$ in $[2B_{00}-K]$ is $[-M_t^2/2+s/3]/(16\pi^2)$ prior to subtraction. The code retains it explicitly but does not treat this isolated term as the full electroweak counterterm. For the finite part, integrating at $s=0$ gives $-M_t^2[1/4+\ln(\bar\mu^2/M_t^2)/2]/(16\pi^2)$, and its first $s$ derivative is $[\ln(\bar\mu^2/M_t^2)/3+1/9]/(16\pi^2)$. The scalar $B_0$ is also checked against its closed log form for one massless propagator. The direct parameter integration and the PV reduction agree within $6.1\times10^{-64}$ GeV² on the bracket across the supplied scales.

## Consequences for the input equations, not a full match

The selected term in (A28) is $d_g^{tb}=[\Pi_W^{tb}(M_W^2)-\Pi_W^{tb}(0)]/M_W^2$. In the VEV relation before (A29), $\delta v/v=[\Pi_W(M_W^2)/M_W^2-d_g]/2$, so the difference of this **particular W block** reduces to $\Pi_W^{tb}(0)/(2M_W^2)$. For the W-input part of (A30) it reduces to $-\Pi_W^{tb}(0)/M_W^2$. That is not the full Yukawa shift: among the missing terms is the top self-energy $\Sigma_s+\Sigma_v$.

At $\bar\mu=M_Z=91.1876$ GeV:

| Selected top–bottom contribution | Value |
|---|---:|
| $\Pi_W^{tb}(0)$, finite (GeV²) | +94.01992425 |
| $\Pi_W^{tb}(M_W^2)$, finite (GeV²) | +79.05058670 |
| A28 difference $d_g^{tb}$ | −0.00231695385 |
| W-input contribution to $\delta v/v$ | +0.00727620127 |
| W-input-only contribution to $\delta y_t^2/y_{t,0}^2$ | −0.01455240254 |

Setting $D=4$ before expanding the tensor denominator, despite keeping the scalar finite parts, would instead give $d_g^{tb}=-0.00141712734$: the missing finite piece is $-g_0^2 N_c/(9\,16\pi^2)=-0.00089982651$. This demonstrates a distinct check from the erratum's original implementation sign issue. Values at $M_W$ and $M_t$ are diagnostics of the scale dependence of an isolated component. No scale cancellation, effective $v$, or complete $g$ or $y_t$ has been computed.

## Replay

The eight-file audit is self-contained. Only the input JSON is numerical input; paper PDFs are not redistributed. GPT-managed Linux, Python 3.12.14, mpmath 1.3.0; 65 working digits. No Wolfram, DRalgo, or independent GitHub Actions run was used.

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_wtb
    python wtb.py /tmp/xsm-wtb-replay.json
    cmp native-results.json /tmp/xsm-wtb-replay.json

Exact byte replay assumes matching software versions.
