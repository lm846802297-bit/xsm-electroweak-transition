# Exact field-basis map

Write the general renormalizable potential as
[
V=c_H X+\lambda X^2+\frac{a_1}{2}XS+\frac{a_2}{2}XS^2
+b_1S+\frac{b_2}{2}S^2+\frac{b_3}{3}S^3+\frac{b_4}{4}S^4,
\qquad X=H^\dagger H.
]
The coefficient `c_H` is stored directly. Liu--Xie write `c_H=-mu2_LX`; the later EFT convention uses `c_H=+mu2_RMT`.

For `S_old=S_new+delta`, exact re-expansion gives
[
\begin{aligned}
c_H'&=c_H+\frac12a_1\delta+\frac12a_2\delta^2,\\
a_1'&=a_1+2a_2\delta,\quad a_2'=a_2,\\
b_1'&=b_1+b_2\delta+b_3\delta^2+b_4\delta^3,\\
b_2'&=b_2+2b_3\delta+3b_4\delta^2,\\
b_3'&=b_3+3b_4\delta,\quad b_4'=b_4.
\end{aligned}
]
The shift changes the Higgs quadratic coefficient; omitting it is an error.

Set `b1_LX=0` and `S_LX=S_0+v_s`. Then
[
\begin{aligned}
\mu^2_{RMT}&=-\mu^2_{LX}+\frac12a_{1,LX}v_s+\frac12a_2v_s^2,\\
a_{1,RMT}&=a_{1,LX}+2a_2v_s,\\
b_{1,RMT}&=b_{2,LX}v_s+b_{3,LX}v_s^2+b_4v_s^3,\\
b_{2,RMT}&=b_{2,LX}+2b_{3,LX}v_s+3b_4v_s^2,\\
b_{3,RMT}&=b_{3,LX}+3b_4v_s,\quad b_{4,RMT}=b_4.
\end{aligned}
]
The zero-temperature stationarity conditions are `mu2_RMT=-lambda v^2` and `b1_RMT=-a1_RMT v^2/4`.

The 3D zero-mode normalization is `S_4=sqrt(T) S_3`; hence `delta_3=delta/sqrt(T)`.
