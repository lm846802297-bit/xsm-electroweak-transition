# Physical-Higgs seagulls in the vacuum W self-energy

This audit closes only the explicitly physical-Higgs-containing **1PI bosonic W diagrams** in Eq. (A58) at one published xSM parameter point: the two new W–W–h–h seagulls plus the previously verified mixed Higgs–W and Higgs–charged-Goldstone bubbles. The common one-point insertions are separately 1PR, and are not counted here. This is an algebraic gauge-fixed diagram subset, not an independently gauge-invariant observable, the complete W/Z self-energy, or a complete pole-to-running input.

## Primary sources and fixed inputs

[Niemi, Schicho and Tenkanen (NST), Phys. Rev. D 103 115035](https://doi.org/10.1103/PhysRevD.103.115035), [arXiv v3](https://arxiv.org/abs/2103.07467v3), Eqs. (A20), (A26)–(A30), (A36), (A58), with [APS erratum](https://doi.org/10.1103/PhysRevD.109.039902). The displayed seagull term is in (A58), and is not amended by the erratum. The parameter point is [Ramsey-Musolf, Tenkanen and Tran, JHEP 09 (2026) 260](https://doi.org/10.1007/JHEP09(2026)260), §4.1: 350 GeV heavy scalar and positive journal sin(theta)=0.1. We use NST's tree-minimum basis with zero singlet VEV; the sign of the light/heavy doublet components drops out of diagonal contractions. Light-scalar mass 125.10 GeV, pole W 80.379 GeV, and Gmu = 1.1663787e−5 GeV−2. The mu-bar values MZ, MW and Mt match the two pinned bubble inputs; no conversion to a new physical scheme is inferred. The pin for the journal input and both NST PDFs is inherited through the hashed [prior input](../w-mixed-physical-scalar-audit/inputs.json). The already included [1PR tadpole audit](../top-shared-1pr-tadpole-audit/) and [mixed bubbles](../w-mixed-physical-scalar-audit/) are untouched.

## Quartic coefficient and sign

Let the neutral doublet excitation be $h_d = d_1 h_1+d_2 h_2$, with $d_1^2=0.99$, $d_2^2=0.01$. The charged-vector squared mass in the scalar background, derived from the doublet kinetic term, is
$$
 m_W^2(h)={g_0^2\over4}(v+d_1h_1+d_2h_2)^2,\qquad
 g_0^2=4\sqrt2G_\mu M_W^2 .
$$
The coincident contraction gives $\tfrac12\sum_{ij}(\partial_i\partial_j m_W^2)\langle h_i h_j\rangle$. In the tree mass basis $\langle h_i h_j\rangle=\delta_{ij}G_i$; the Hessian has diagonal $g_0^2 d_i^2/2$. Equivalently, the $W^+W^-h_ih_i$ vertex carries this same $g_0^2 d_i^2/2$ and the identical-scalar loop has symmetry factor $1/2$. Thus the coefficient of each A36 tadpole integral is $g_0^2d_i^2/4$. The positive Euclidean mass shift from the coincident propagator enters the NST inverse-propagator-plus-$\Pi$ convention with opposite sign. This is also the Wick sign of their $1/[i(2\pi)^D(p^2-m^2+i0)]$ integral. The resulting finite NST seagull is
$$
 \Pi_{W,\mathrm{sea}}(k^2)= {g_0^2\over4} \sum_i d_i^2 A_{\rm fin}(M_{h_i}),
 \qquad A_{\rm fin}(m)={m^2\over16\pi^2}\left(1+\ln{\bar\mu^2\over m^2}\right).
$$
There is no transverse dimension-dependent numerator, so no extra $D\times1/\epsilon$ term here, unlike in the separately audited tensor bubble. The subtraction is the A36 MS-bar pole. The term has no external momentum or absorptive cut.

## Independent regulated check

With Euclidean proper-time cutoff $\tau\geq \Lambda^{-2}$, set $z=m^2/\Lambda^2$. The heat-kernel momentum integral gives
$$
 16\pi^2 I_\Lambda=\int_{\Lambda^{-2}}^\infty {e^{-\tau m^2}\over \tau^2}\,d\tau
 =\Lambda^2 e^{-z}-m^2E_1(z).
$$
The dimensionful matching subtraction is $\Lambda^2+m^2[\ln(\bar\mu^2/\Lambda^2)+\gamma_E]$ in the same units. After subtracting and reversing the Euclidean sign, $A_{\Lambda\to\overline{\mathrm{MS}}}$ tends to $A_{\rm fin}$. The remainder before reversal, in units $m^2/(16\pi^2)$, is $f(z)=(e^{-z}-1)/z-E_1(z)-\ln z-\gamma_E+1$; hence $f(0)=0$ and $f'(z)=(1-e^{-z}-z)/z^2\in[-1/2,0]$ for positive $z$. It follows that the signed error after reversal is nonnegative and bounded by $m^4/(32\pi^2\Lambda^2)$ before the quartic coefficient. This bound was checked individually at $\Lambda/m=10^3,10^{12},10^{26}$, with precision guards. At the last ratio the largest permitted one-diagram discrepancy is $5.229\times10^{-52}$ GeV². Omitting the finite Euler constant conversion offsets the total by $-6.51261976616$ GeV², not a cutoff effect. Independently the derivative with respect to $\ln\bar\mu^2$ is $+11.2828188182874$ GeV². The exact sum of mixing weights preserves the degenerate-mass result and zero mixing projects only onto the doublet.

## Signed result at $\bar\mu=M_Z$

| 1PI piece | $\Pi_W(0)$, GeV² | $\mathrm{Re}\,\Pi_W(M_W^2)$, GeV² | A28 relative | A30 relative, W input only |
|---|---:|---:|---:|---:|
| $h_1$ seagull | +3.84377867176811 | +3.84377867176811 | 0 | |
| $h_2$ seagull | −1.39717290595516 | −1.39717290595516 | 0 | |
| Seagulls | +2.44660576581295 | +2.44660576581295 | 0 | −0.000378685605610880 |
| Seagulls plus pinned mixed bubbles | −11.3586852127162 | −9.81544712847906 | +0.000238862368713873 | +0.00175809713556018 |

The linearized W-self-energy part of (A28) is $(\mathrm{Re}\,\Pi_W(M_W^2)-\Pi_W(0))/M_W^2$; therefore it vanishes identically for the seagull. Substituting this into (A30), its W-input-only piece becomes $-\Pi_W(0)/M_W^2$. This projection does not modify or recompute the separate top self-energy. The mixed-bubble values are taken from their immutable prior signed file, not fitted (their archived precision controls the total). The sign of the individual diagram block changes with scale, as expected for a finite subset of running-parameter contributions. No longitudinal/gauge-independence claim, Z matching, other W sectors, complete pole match, thermal analysis or transition predictions follow.

## Replay

Using Python 3.12 and mpmath 1.3.0, from this directory:

    python -m unittest -v test_seagull.py
    output=$(mktemp)
    python seagull.py "$output"
    cmp native-results.json "$output"
    sha256sum -c SHA256SUMS.txt

See [run manifest](run_manifest.json), [native signed results](native-results.json), [test log](validation.log), and [standards validators](provenance-validation.log). The guarded cutoff route and the main A36 path are independent under the explicitly matched scheme; this does not constitute independent validation of all of (A58).
