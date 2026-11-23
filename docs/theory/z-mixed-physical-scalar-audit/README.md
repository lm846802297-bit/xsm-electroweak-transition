# Mixed physical-Higgs bubbles in the vacuum Z self-energy

This audit isolates the two physical-Higgs–Z and two physical-Higgs–neutral-Goldstone 1PI diagrams in NST Eq. (A59). It gives their signed contribution to (A29) with the W input fixed, at the existing 350 GeV benchmark. Scalar seagulls, common 1PR tadpoles, fermions and other W/Z contributions are not included.

## Source and normalization

[Niemi, Schicho and Tenkanen, Phys. Rev. D 103 115035](https://doi.org/10.1103/PhysRevD.103.115035), [arXiv v3](https://arxiv.org/abs/2103.07467v3), §II, (A26)–(A29), (A36)–(A44), (A59), and [APS erratum](https://doi.org/10.1103/PhysRevD.109.039902). The erratum does not modify the selected A59 terms. The frozen input is [Ramsey-Musolf, Tenkanen and Tran, JHEP 09 (2026) 260](https://doi.org/10.1007/JHEP09(2026)260), §4.1, as previously audited: heavy scalar 350 GeV, light scalar 125.10 GeV, journal mixing sine 0.1, pole W 80.379 GeV, pole Z 91.1876 GeV, Fermi constant 1.1663787e−5 GeV−2. In the NST zero-temperature basis the singlet VEV at tree level is zero. Work in Feynman gauge and MS-bar with the NST inverse-propagator-plus-self-energy sign. Exact physics, source PDF and prior algorithm hashes are pinned via [inputs.json](inputs.json).

In the neutral doublet component $\phi^0=(v+h_d+iG^0)/\sqrt2$, the unbroken neutral rotation gives zero photon charge and $|q_Z|=\sqrt{g_0^2+g_0'^2}/2=g_{Z0}/2$. The scalar kinetic term contains the $Z_\mu Z^\mu(v+h_d)^2$ mass term and the current $q_Z Z^\mu(h_d\partial_\mu G^0-G^0\partial_\mu h_d)$ up to the overall sign of the neutral gauge field. Writing $h_d=\sum_i d_i h_i$, $d_i^2=(0.99,0.01)$, differentiating the identical Z legs gives the $ZZh_i$ coefficient $g_{Z0} M_Z d_i$. In the derivative $Zh_iG^0$ vertex the coefficient $q_Z d_i$ multiplies the momentum difference; after transverse projection this difference is twice the loop momentum, giving $g_{Z0}d_i$. Neither mixed two-line loop has an identical internal particle factor. The vector and derivative-scalar loops in NST's transverse sign convention give, respectively, the metric $+M_Z^2B_0$ and $-B_{00}$ terms in
$$
  \Pi_{Z,\mathrm{mixed}}(s)=g_{Z0}^2\sum_i d_i^2
  \left[M_Z^2\, B_0(s;m_i,M_Z)-B_{00}(s;m_i,M_Z)\right],
  \qquad g_{Z0}^2=4\sqrt{2}G_\mu M_Z^2.
$$
The latter uses $m_{G^0}=M_Z$ at $\xi=1$. Equal Goldstone/vector masses are part of the fixed prescription, not a blind W-to-Z replacement or a claim that the separate two diagrams are physical.

## Signed checks

Write $a=m_i^2$, $b=M_Z^2$, $u=x(1-x)$, and $\Delta=(1-x)a+xb-su$. The direct D-dimensional shifted integration with transverse tensor IBP has MS-bar finite integrand $bL-\tfrac12\Delta(1+L)$, with $L=\ln[\bar\mu^2/(\Delta-i0)]$; its prefactor is $g_{Z0}^2d_i^2/(16\pi^2)$. Independent algebraic (A41)/(A42) reduction of the tensor integral retains the finite contribution from $1/[2(D-1)]$ times the pole: $(3(a+b)-s)/18$ in $16\pi^2B_{00}$. Omitting it by putting $D=4$ early gives the wrong signed shift. The direct and reduced values agree within $10^{-48}$ GeV² at both $s=0$ and $s=M_Z^2$, individually for both Higgs masses and each scale. The analytic slope per state in $\ln\bar\mu^2$ is $g_{Z0}^2d_i^2[b-(a+b)/4+s/12]/(16\pi^2)$.

As a separate route, the on-shell physical $h_i Z$ cut of the transverse projection uses the massive-vector polarization sum, with the angular average of its spatial longitudinal part $|\mathbf p|^2/(3M_Z^2)$. This leaves the spectral density
$$
 \rho_i(t)=\frac{g_{Z0}^2 d_i^2 M_Z^2}{16\pi}\,
  \beta(t)\left[1+\frac{t\beta(t)^2}{12M_Z^2}\right],\qquad
 \beta(t)=\frac{\sqrt{[t-(m_i+M_Z)^2][t-(m_i-M_Z)^2]}}{t}, \qquad t>(m_i+M_Z)^2 .
$$
This follows from distinguishable two-body phase space and the optical-theorem factor $1/2$. It independently reproduces the signed imaginary part of the combined loops at $1.25$ times each squared threshold: $+44.3568933219$ and $+0.467641121846$ GeV², respectively. The physical Z pole is below both cuts. Applying two subtractions removes the UV polynomial; the curvature is
$$
 \Pi_{Z,i}''(s)=\frac{2}{\pi}\int_{(m_i+M_Z)^2}^\infty
       \frac{\rho_i(t)}{(t-s)^3}\,dt,
$$
which at the pole matches differentiating the shifted integral twice, replacing its integrand by $u^2[b/\Delta^2+1/(2\Delta)]$. It gives $1.53839393943\times10^{-8}$ and $1.44842465007\times10^{-11}$ GeV−2 for the two states, with maximum replay discrepancy $1.8\times10^{-77}$ GeV−2 at 65 working digits. This spectral test does not fix either subtraction constant; those are checked by the signed PV reduction in the specified scheme. The vanishing physical absorptive piece and the positive cut test are not claims about other bosonic channels.

## Benchmark result at the Z scale

| Internal Higgs | $\Pi_Z(0)$ (GeV²) | $\mathrm{Re}\,\Pi_Z(M_Z^2)$ (GeV²) |
|---|---:|---:|
| $h_1$ | −22.7433397692428 | −19.7752177485705 |
| $h_2$ | +0.821457998172492 | +0.782193421456794 |
| Sum | −21.9218817710703 | −18.9930243271137 |

Holding all W terms fixed, (A29) gives $\delta g'^2_Z=(g_0^2/M_W^2)\operatorname{Re}\Pi_Z(M_Z^2)=-0.00125316626077627$, or $\delta g'^2_Z/g_0'^2=-0.0102421689334729$. The wrong early-$D=4$ relative value is −0.00329070278235969. The first column is a validation point, not an input to the A29 projection. No part of the earlier W audits or common 1PR calculation was added; neither this signed projection nor its scale slope is the full running hypercharge coupling.

## Replay

Python 3.12 and mpmath 1.3.0, from this directory:

    PYTHONDONTWRITEBYTECODE=1 python -m unittest -v test_mixed_z.py
    output=$(mktemp)
    PYTHONDONTWRITEBYTECODE=1 python mixed_z.py "$output"
    cmp native-results.json "$output"
    rm "$output"
    sha256sum -c SHA256SUMS.txt

[Signed native results](native-results.json), [run manifest](run_manifest.json), [local test log](validation.log) and [exact standards proof](provenance-validation.log). No CI run is claimed.
