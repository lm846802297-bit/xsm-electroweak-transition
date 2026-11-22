"""Exact NLO covariance check, retaining terms through one loop factor."""

import sympy as sp


def test_nlo_scalar_matching_commutes_with_field_translation_symbolically():
    T, d, a1, a2, b1, b2, b3, b4, lam, cH = sp.symbols(
        "T d a1 a2 b1 b2 b3 b4 lam cH", nonzero=True
    )
    g2, gp2, yt2, Lb, Lf, E, L3 = sp.symbols("g2 gp2 yt2 Lb Lf E L3")
    rT = sp.sqrt(T)

    def match(p):
        cH_, lam_, a1_, a2_, b1_, b2_, b3_, b4_ = p
        A1 = rT * a1_
        A2 = T * a2_
        B1 = (b1_ + T**2 * (a1_ + b3_) / 12) / rT
        B2 = b2_ + T**2 * (a2_ / 6 + b4_ / 4)
        B3 = rT * b3_
        B4 = T * b4_
        G3, GP3 = T * g2, T * gp2

        A1n = A1 + rT * E * (
            Lb * (sp.Rational(3, 4) * (3 * g2 + gp2) * a1_ - 2 * b3_ * a2_ - (6 * lam_ + 2 * a2_) * a1_)
            - 3 * Lf * yt2 * a1_
        )
        A2n = A2 + T * E * (
            Lb * (sp.Rational(3, 4) * (3 * g2 + gp2) - 6 * lam_ - 2 * a2_ - 3 * b4_) * a2_
            - 3 * Lf * yt2 * a2_
        )
        B1n = B1 + E / rT * (
            -Lb * (a1_ * cH_ + b3_ * b2_)
            + T**2 * (
                (2 + 3 * Lb) / 48 * (3 * g2 + gp2) * a1_
                - Lb / 2 * ((lam_ + sp.Rational(7, 12) * a2_) * a1_ + (a2_ / 3 + sp.Rational(3, 2) * b4_) * b3_)
                - sp.Rational(1, 8) * a1_ * yt2 * (3 * Lb - Lf)
            )
        ) - E * (2 * B3 * B4 - sp.Rational(1, 2) * A1 * (3 * G3 + GP3 - 2 * A2)) * L3
        B2n = B2 - Lb * E * (2 * b3_**2 + sp.Rational(1, 2) * a1_**2 + 2 * a2_ * cH_ + 3 * b4_ * b2_)
        B2n += T**2 * E * (
            (2 + 3 * Lb) / 24 * (3 * g2 + gp2) * a2_
            - Lb * ((lam_ + sp.Rational(7, 12) * a2_ + sp.Rational(1, 2) * b4_) * a2_ + sp.Rational(9, 4) * b4_**2)
            - sp.Rational(1, 4) * a2_ * yt2 * (3 * Lb - Lf)
        )
        B2n += E * ((3 * G3 + GP3) * A2 - 2 * A2**2 - 6 * B4**2) * L3
        B3n = B3 - rT * E * 3 * Lb * (sp.Rational(1, 2) * a1_ * a2_ + 3 * b4_ * b3_)
        B4n = B4 - T * E * Lb * (a2_**2 + 9 * b4_**2)
        return [A1n, A2n, B1n, B2n, B3n, B4n]

    p = [cH, lam, a1, a2, b1, b2, b3, b4]
    shifted = [
        cH + a1 * d / 2 + a2 * d**2 / 2,
        lam,
        a1 + 2 * a2 * d,
        a2,
        b1 + b2 * d + b3 * d**2 + b4 * d**3,
        b2 + 2 * b3 * d + 3 * b4 * d**2,
        b3 + 3 * b4 * d,
        b4,
    ]

    matched = match(p)
    matched_shifted = [
        matched[0] + 2 * matched[1] * d / rT,
        matched[1],
        matched[2] + matched[3] * d / rT + matched[4] * d**2 / T + matched[5] * d**3 / T**sp.Rational(3, 2),
        matched[3] + 2 * matched[4] * d / rT + 3 * matched[5] * d**2 / T,
        matched[4] + 3 * matched[5] * d / rT,
        matched[5],
    ]
    shifted_matched = match(shifted)

    for lhs, rhs in zip(matched_shifted, shifted_matched):
        residual = sp.series(sp.expand(lhs - rhs), E, 0, 2).removeO().expand()
        assert sp.simplify(residual) == 0
