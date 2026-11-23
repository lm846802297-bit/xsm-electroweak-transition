# xSM soft-to-ultrasoft temporal threshold at one loop

**Scope:** This extends the [hard static temporal portal audit](../static-temporal-portal-audit/README.md), not the 4D high-temperature thermal tadpole of Liu–Xie [Eq. (2.8)](https://arxiv.org/pdf/2101.10469#page=4). The checked starting point was [draft PR #1](https://github.com/lm846802297-bit/xsm-electroweak-transition/pull/1) at 5db03592fe3719234cc5264d770065e1908eda47; main was 0e421748932d8b0a0ec635bc909a2f9be06f62c6. Here we integrate only the three SU(2), one U(1), and optionally eight color temporal fields in the Gaussian **static** soft EFT. This is neither the full tadpole nor a thermal transition calculation. All fields, masses, and couplings below are 3D unless marked 4D.

## Action normalization and derivation

Write the canonically normalized background $\rho=H^\dagger H$ with $H=(0,\sqrt{\rho})^T$. Match the polynomial coefficients in
$$
V_{\mathrm{local}}= b_1S+\tfrac12 b_2 S^2+c_H\rho+\lambda\rho^2+
\tfrac12 a_1\rho S+\tfrac12 a_2\rho S^2+\cdots.
$$
A temporal action term $x_i S X_i^2+y_i S^2 X_i^2$ contributes $2x_iS+2y_iS^2$ to the squared mass of a real mode. Put $d_i(S)=m_{D,i}^2+2x_iS+2y_iS^2$, $u_i=x_i+2y_iS$ and $m_i=\sqrt{d_i(S)}$. For $i=A,B,C$, let $h_i=(h_3,h_3',\omega_3)$, with $x_C=y_C=0$. From the primary 2026 action Eq. (B.1), $A^{1,2}$ each have squared mass $d_A+2h_A\rho$ and the neutral block is
$$
{\cal M}^2_{A^3,B} =
\begin{pmatrix}d_A+2h_A\rho&-k\rho\\-k\rho&d_B+2h_B\rho\end{pmatrix},
\qquad k=h_3''.
$$
The eight color modes each have squared mass $d_C+2h_C\rho$. The off-diagonal term is $-k\rho$ because the action has $k H^\dagger \tau^aH A_0^a B_0$ and uses Pauli matrices, not $\tau^a/2$ here. In 3D dimensional regularization the real-scalar static contribution is $J(m)=-m^3/(12\pi)$, so the threshold potential is $\Delta V=-\mathrm{Tr}[({\cal M}^2)^{3/2}]/(12\pi)$.

At any constant $S=s$ and $\rho=0$, writing $n=(3,1,8)$, the expansion in $\delta S=S-s$ gives:
$$
\begin{aligned}
\Delta b_1 &= -\frac1{4\pi}\sum_i n_i m_i u_i,\\
\Delta b_2 &= -\frac1{2\pi}\sum_i n_i m_i y_i-\frac1{4\pi}\sum_i n_i u_i^2/m_i,\\
\Delta c_H &= -\frac1{4\pi}\sum_i n_i m_i h_i,\\
\Delta a_1 &= -\frac1{2\pi}\sum_i n_i h_i u_i/m_i,\\
\Delta a_2 &= -\frac1{2\pi}\sum_i n_i h_i y_i/m_i+\frac1{4\pi}\sum_i n_i h_i u_i^2/m_i^3,\\
\Delta\lambda &= -\frac1{8\pi}\left(\sum_i n_i h_i^2/m_i+\frac{k^2}{m_A+m_B}\right).
\end{aligned}
$$
The last term stays finite in the degenerate limit $m_A=m_B=m>0$. Neither the first $\rho$ derivative nor its $S$ derivatives can receive a single off-diagonal mass insertion. Exact recentering $S=S'+d$ gives $m_{D,i}^{\prime 2}=d_i(d)$, $x_i'=x_i+2y_i d$ and $y_i'=y_i$, with $s'=s-d$. Both determinant and derivatives are identical at the same physical point. This is action-level covariance; substituting just a shifted portal into a hard matching formula without the corresponding temporal mass change is not the same test.

## Independent momentum check and conventions

Differentiate $\tfrac12\,\mathrm{Tr}\log(p^2+{\cal M}^2)$ under the integral, rather than differentiating the sum of eigenvalue powers. Single-propagator terms use $\int_p[(p^2+m^2)^{-1}-p^{-2}]=-m/(4\pi)$, equivalent to setting the scaleless integral to zero in 3D dimensional regularization. The two-neutral-propagator term uses $\int_p[(p^2+m_A^2)(p^2+m_B^2)]^{-1}=1/[4\pi(m_A+m_B)]$. The replay instead evaluates these convergent radial integrals directly at 65 decimal digits. The action derivatives at zero external momenta are not a determination of a real-time propagator pole or of a screening pole at nonzero spatial momentum. No 3D logarithmic UV pole is subtracted in this one-loop potential. The earlier **hard** UV-pole-times-evanescent calculation is separate.

At $s=0$ the singlet derivatives reproduce [Schicho, Tenkanen and Österman, JHEP06(2021)130, Eqs. (4.60), (4.61), (4.64), (4.65)](https://inspirehep.net/files/2e6c3509f604fc4d17e096fd09f8f988) and [Ramsey-Musolf, Tenkanen and Tran, JHEP09(2026)260, Eqs. (B.26), (B.27)](https://link.springer.com/content/pdf/10.1007/JHEP09(2026)260.pdf#page=26), as functions of $x_i,y_i,m_i$. The explicit Higgs-operator shifts $\Delta c_H$, $\Delta\lambda$ give the displayed one-loop bracket of journal (B.30) and EW bracket of (B.28), respectively. In the [2026 arXiv PDF](https://arxiv.org/pdf/2409.17554#page=27) those numbers are (B.27), (B.28), (B.31) and (B.29), respectively. The printed (B.2) includes $1/2$ and $1/4$ on several doublet terms: if read literally with $\rho=\phi^\dagger\phi$ of (B.1), then $c_H=\mu^2_3/2$ and $\lambda=\lambda_3/4$ in our definition, and shifts of those literal symbols require multiplying our coefficient shifts by 2 and 4. We compare action coefficients and brackets, not equate potentially mismatched symbol conventions. The optional color contribution in $\Delta\lambda$ is not in the EW bracket and is higher order for $\omega_3\sim g^4T$.

## 2021 literature-attribution erratum

The previous hard audit said that Schicho et al. (4.36), (4.37) **display** the same linear normalization as the independent half candidate. They do not. The [2021 version of record](https://inspirehep.net/files/2e6c3509f604fc4d17e096fd09f8f988), printed p. 30, displays $x_3=Tg^2\mu_m/(4\pi)^2$ and $x_3'=Tg'^2\mu_m/(4\pi)^2$, namely $T$ (not $\sqrt T$) and no $1/2$. The [arXiv v1 PDF](https://arxiv.org/pdf/2102.11145v1#page=33) displays the same. In their four-dimensional potential (4.1) the term is $\mu_m\sigma\phi^\dagger\phi/2$, so $\mu_m$ has mass dimension 1 there. Their 3D quadratic temporal masses (4.58), (4.59) use $2x_i s_3$, implying $[x_i]=3/2$ since $[s_3]=1/2$. With those literal definitions, the printed $T\mu_m$ has dimension 2 and cannot be directly compared as an already normalized 3D action coefficient. We do not infer what the authors used numerically. The independent hard determinant (and the separate, source-only DRalgo tensor-to-action comparison recorded on this branch) is unaffected. The 2021 soft formulas are meaningful as functions of consistently normalized 3D $x_i$, independently of this attribution correction.

## Matched comparison and limits

Keep $y_i$, $m_{D,i}$ and $h_i$ fixed and compare only the 2026 printed linear $x_i$ with the hard-static candidate $x_i/2$, at $S=0$. The temporal pieces $\Delta b_1$ and $\Delta a_1$ halve, the $x_i^2$ contributions to $\Delta b_2$ and $\Delta a_2$ quarter, and the $y_i$ terms stay unchanged. At nonzero $S$ the masses and $u_i$ also depend on these couplings, so no universal rescaling follows. Under the illustrative strict counting $a_{1,4}/T,a_{2,4}\sim g^2$, $m_D\sim gT$, and $h_{A,B},k\sim g^2T$, one has $x\sim g^4T^{3/2}$ and $y\sim g^4T$; the temporal $\Delta b_1$ is $O(g^5T^{5/2})$, the $my$ contribution to $\Delta b_2$ is $O(g^5T^2)$, and the $x^2/m$ term is $O(g^7T^2)$. These displayed diagnostic pieces are not the complete prediction at either higher order and do not change the formal $O(g^4)$ analysis. For the two supplied fixtures we use leading Debye masses and leading electroweak $h$ values, and the displayed leading nonzero $\omega_3$. The two configurations are synthetic checks, not sampled parameter points used in a phase-transition result.

At the $T=120$ GeV test point, the temporal $\Delta b_1$ changes from $-35.1366012$ to $-17.5683006$ GeV$^{5/2}$. The temporal $\Delta b_2$ changes only from $-4.70851469$ to $-4.70539867$ GeV$^2$: its $y$ piece is $-4.70435999$ GeV$^2$ in both cases. At the equal-Debye fixture, the quartic neutral-mixing contribution is $-0.07934200793$ GeV and agrees with directly differentiating the eigenvalue trace and integrating the two propagators. Across both runs, the maximum scaled error between analytic derivatives, radial integrals and the shift identity is $2.17\times10^{-66}$ (scale $\max(1,|u|,|v|)$ in fixed GeV-based units). Positivity of temporal squared masses is required; no extrapolation through a vanishing soft mass is tested.

No local Wolfram kernel was found; DRalgo was not run. No 4D thermal tadpole recalculation, full two-loop soft matching, transition observable, or gravitational-wave prediction is claimed.

## Replay

GPT-managed Linux, Python 3.12.14, mpmath 1.3.0. The 2021 journal PDF inspected by layout-preserving text extraction had SHA-256 ee791c67600d5dc240f8fe30355bdff026891cdb1ba11c244722e61642717443; no PDF is required at replay. The only numerical inputs are the inline fixtures. GitHub Actions was not run.

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_threshold
    python threshold.py /tmp/xsm-threshold-replay.json
    cmp threshold-results.json /tmp/xsm-threshold-replay.json

The exact byte comparison assumes the same recorded runtime; with a different runtime compare parsed values within tolerance.
