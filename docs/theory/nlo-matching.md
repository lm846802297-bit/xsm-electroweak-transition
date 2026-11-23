# NLO dimensional-reduction matching audit

Reference: Ramsey-Musolf, Tenkanen and Tran, JHEP 09 (2026) 260, arXiv:2409.17554, Appendix B, Eqs. (B.5)--(B.10). The NLO hard-mode tadpole structure is independently present in Schicho, Tenkanen and Osterman, JHEP 06 (2021) 130, Eq. (4.42).

## LO matching

[
\begin{aligned}
a_{1,3}&=\sqrt T a_1,&a_{2,3}&=T a_2,\\
b_{1,3}&=\frac{1}{\sqrt T}\left[b_1+\frac{T^2}{12}(a_1+b_3)\right],\\
b_{2,3}&=b_2+T^2\left(\frac{a_2}{6}+\frac{b_4}{4}\right),\\
b_{3,3}&=\sqrt T b_3,&b_{4,3}&=T b_4.
\end{aligned}
]

Define
[
L_b=2\ln\frac{\Lambda}{T}-2[\ln(4\pi)-\gamma_E],\quad L_f=L_b+4\ln2,
]
[
c=\frac12\left[\ln\frac{8\pi}{9}+\frac{\zeta'(2)}{\zeta(2)}-2\gamma_E\right].
]

The NLO tadpole is
[
\begin{aligned}
b_{1,3}={}&\frac{1}{\sqrt T}\Bigg\{b_1+\frac{T^2}{12}(a_1+b_3)-\frac{L_b}{(4\pi)^2}(a_1c_H+b_3b_2)\\
&+\frac{T^2}{(4\pi)^2}\left[\frac{2+3L_b}{48}(3g^2+g'^2)a_1-\frac{L_b}{2}\left((\lambda+\frac{7a_2}{12})a_1+(\frac{a_2}{3}+\frac{3b_4}{2})b_3\right)-\frac18a_1y_t^2(3L_b-L_f)\right]\Bigg\}\\
&-\frac{1}{(4\pi)^2}\left[2b_{3,3}b_{4,3}-\frac12a_{1,3}(3g_3^2+g_3'^2-2a_{2,3})\right]\left(c+\ln\frac{3T}{\Lambda_3}\right).
\end{aligned}
]

The remaining NLO equations for portal couplings, singlet mass, cubic and quartic couplings are implemented in `src/xsm_ewpt/matching.py`. The parent-theory Higgs quadratic parameter is `cH`, opposite in sign to positive `mu2_LX` in Liu--Xie.

The zero-mode normalization gives `delta_3=delta/sqrt(T)`. The code tests matching/translation covariance symbolically and on a deterministic 100-point ensemble. This milestone does not include temporal gauge-field integration, bounce, percolation or GW spectra.
