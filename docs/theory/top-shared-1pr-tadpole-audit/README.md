# Common one-loop tadpoles in the xSM top and VEV inputs

The same 1PR Higgs tadpole enters the W, Z, and top two-point functions in the tree-minimum scheme of [Niemi, Schicho and Tenkanen (NST)](https://doi.org/10.1103/PhysRevD.103.115035). This audit computes the signed one-point functions and verifies their **paired cancellation in (A30)**. It does not calculate other two-point terms or a full electroweak match.

Sources: [NST arXiv:2103.07467v3](https://arxiv.org/pdf/2103.07467v3#page=20), §II, (A20)–(A23), (A26)–(A30), (A36), (A45)–(A48), (A54)–(A55), (A58)–(A59) and (A61); their [published erratum](https://doi.org/10.1103/PhysRevD.109.039902); [Martin and Patel, arXiv:1808.07615](https://arxiv.org/pdf/1808.07615#page=12), (3.1)–(3.3), (3.69)–(3.70) for the independent potential route; the [published Ramsey-Musolf–Tenkanen–Tran benchmark](https://link.springer.com/content/pdf/10.1007/JHEP09(2026)260.pdf#page=15), §4.1.

The scalar masses are 125.10 and 350 GeV. The journal mixing sine is +0.1; in the NST ordered basis the sine is −0.1. In particular the heavy eigenstate's doublet component in this signed convention is negative. Use the Fermi-input tree VEV, 246.2196507941 GeV, rather than the earlier rounded engineering fixture of 246 GeV; the scalar-potential benchmark values are a2=3, b3=40 GeV and b4=0.3. The prior [physical-scalar top audit](../top-physical-scalar-audit/) pins the exact common EW input. The [task brief](task-brief.md) and [task record](../../../research/standards-provenance/tasks/20261123-top-shared-1pr-tadpole.json) record the gate and actual lock at base 7dc6adc2.

## Two independent signed evaluations

The loop convention is four-dimensional zero-temperature MS-bar, positive self-energy in the inverse propagator, and Feynman–'t Hooft gauge. Only top is retained among fermions in (A54)–(A55). The diagram route evaluates their four signed classes (top, vector, Goldstone and physical scalar), using the printed (A45)–(A48) vertex coefficients. The vector numerator must be multiplied by the UV tadpole before subtraction:

$$
A_{\mathrm{fin}}(m)=\frac{m^2}{16\pi^2}\left[1+\ln\frac{\bar\mu^2}{m^2}\right],\qquad
\big[(D-1)A(m)\big]_{\mathrm{fin}}=3A_{\mathrm{fin}}(m)-\frac{2m^2}{16\pi^2}.
$$

For the other route construct the field-dependent CP-even Hessian directly from the interaction-basis tree potential. At each displaced point diagonalize it and evaluate the actual one-loop determinant functions; no printed trilinear formula is used in this route. Let p=v+h be the displaced doublet background, S the singlet, and G the tree-potential Goldstone squared mass. Standard MS-bar functions are:

$$
f_s(q)=\frac{q^2}{64\pi^2}\left[\ln\frac{q}{\bar\mu^2}-\frac32\right],\qquad
f_v(q)=\frac{q^2}{64\pi^2}\left[\ln\frac{q}{\bar\mu^2}-\frac56\right],\qquad
G=\frac{1}{p}\frac{\partial V_{\mathrm{tree}}}{\partial p}.
$$

Here q is a squared mass. In the background R-xi representation the vector, Goldstone and ghost determinants contribute, per massive vector at xi=1:

$$
3f_v(M_V^2(p)) + f_s(M_V^2(p)+G)-f_s(M_V^2(p)).
$$

The last subtraction combines the unphysical longitudinal pole and ghost pair. At the **tree stationary point**, G=0, their common gauge-fixing mass derivatives cancel, leaving the Goldstone-potential derivative; no equality of the off-shell potentials of different gauge-fixing implementations is required. The top has multiplicity −12, and each massive physical scalar has +1. The potential derivative defines the NST tadpole with the minus sign and gives the linearized displacement independently through the inverse interaction-basis Hessian:

$$
t_i=-\frac{\partial V_1}{\partial h_i},\qquad
\begin{pmatrix}\delta h\\\delta S\end{pmatrix}_{\mathrm{local}}
=-H_0^{-1}\nabla_{(h,S)}V_1 .
$$

With 65 working digits two centered five-point stencil steps (10⁻⁹ and 10⁻¹¹ GeV) give maximum scaled diagram-versus-potential differences of 8.26 × 10⁻⁴⁴ and 8.13 × 10⁻⁵². The signs compare separately for each state and each loop class. The diagonalized tree Hessian reproduces both squared pole inputs; the numerical exercise is first-order response, not minimization of the full loop-corrected potential.

## Pole projection and no double count

Define the common doublet displacement with the signed NST overlaps:

$$
T_h=\frac{\cos\theta}{M_{h_1}^2}t_1+\frac{\sin\theta}{M_{h_2}^2}t_2.
$$

Equations (A58), (A59) and (A61) put the same value in the two-point functions as:

$$
\frac{\Pi_W^{\rm tad}}{M_W^2}=\frac{\Pi_Z^{\rm tad}}{M_Z^2}=-\frac{2T_h}{v},\qquad
(\Sigma_s+\Sigma_v)^{\rm tad}=-\frac{T_h}{v}.
$$

These are real and independent of external momentum. The W difference in (A28) is zero, and the equal fractional W/Z shifts cancel in (A29). The key additional cancellation in (A30) is:

$$
\left.\frac{\delta y_t^2}{y_{t,0}^2}\right|_{\rm shared\,tad}
=\left.\frac{\delta g^2}{g_0^2}\right|_{\rm shared\,tad}
-\frac{\Pi_W^{\rm tad}(M_W^2)}{M_W^2}
+2(\Sigma_s+\Sigma_v)^{\rm tad}
=0+\frac{2T_h}{v}-\frac{2T_h}{v}=0 .
$$

The VEV matching is different. Its pole-defined first-order term is minus twice the fractional local displacement:

$$
\left.\frac{\delta v_{\rm pole}^2}{v_0^2}\right|_{\rm tad}=-\frac{2T_h}{v},
\qquad
\left.\frac{\delta v}{v_0}\right|_{\rm local}=\frac{T_h}{v}.
$$

Do **not** insert the displacement as an extra VEV correction. In (A30) include both tadpole pieces once, or cancel them analytically and exclude both. Do not also add scalar external-pole tadpoles, which probe different self-energies, or add any previously computed 1PI contribution twice. This cancellation does not drop tadpoles from all other matching equations.

## Reference numbers at the Z scale

| Signed term | h1 tadpole (GeV³) | h2 tadpole (GeV³) |
|---|---:|---:|
| Top | +76037.060048 | −7642.012071 |
| Vector | +5522.953793 | −555.077741 |
| Goldstone | +4905.434372 | −3859.058319 |
| Massive physical scalars | −475590.806576 | +94668.109728 |
| Total | −389125.358363 | +82611.961597 |

The shared projection gives T_h = −24.8070291114 GeV. The W-input part of (A30) is −0.201503243396 and its top part is +0.201503243396: their sum vanishes, not a physical 20% shift. The isolated pole VEV-squared term is +0.201503243396; the distinct local fractional displacement is −0.100751621698. The maximum numerical residual of (A30) at all three scales is 2.97 × 10⁻⁶⁷. The large individual gauge-fixed terms motivate treating them consistently, not treating their size as a prediction. Bosonic momentum-dependent W/Z terms, the full electroweak/VEV match, and any thermal, bounce or GW output are outside scope.

## Replay

With Python 3.12.14 and mpmath 1.3.0, run from this directory:

    sha256sum -c SHA256SUMS.txt
    python -m unittest -v test_shared_tadpole
    python shared_tadpole.py replay.json
    cmp native-results.json replay.json

The input check also supports the exact optional historic benchmark fixture. The sibling top physical-scalar input is hash-checked at execution. [Run manifest](run_manifest.json) records the actual environment, signed results and validation. No CI run was requested or used.
