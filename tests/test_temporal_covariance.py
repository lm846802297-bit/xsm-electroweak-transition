"""Symbolic audit of the published LO temporal-coupling field-shift relation.

This is a diagnostic test: it records the factor-two residual in the published
x3 matching coefficient and verifies the coefficient required by the action's
finite field translation. It does not declare the paper corrected.
"""

import sympy as sp


def test_published_temporal_x3_has_factor_two_shift_residual():
    T, g2, a1, a2, delta = sp.symbols("T g2 a1 a2 delta", positive=True)
    delta3 = delta / sp.sqrt(T)
    loop = (4 * sp.pi) ** 2

    x3_published = sp.sqrt(T) * g2 * a1 / loop
    y3_published = T * g2 * a2 / (2 * loop)
    a1_shifted = a1 + 2 * a2 * delta

    delta_x_from_matching = sp.simplify(
        (sp.sqrt(T) * g2 * a1_shifted / loop) - x3_published
    )
    delta_x_from_action = sp.simplify(2 * y3_published * delta3)

    assert sp.simplify(delta_x_from_matching - 2 * delta_x_from_action) == 0
    assert sp.simplify(delta_x_from_matching - delta_x_from_action) != 0


def test_factor_half_candidate_for_x3_restores_covariance():
    T, g2, a1, a2, delta = sp.symbols("T g2 a1 a2 delta", positive=True)
    delta3 = delta / sp.sqrt(T)
    loop = (4 * sp.pi) ** 2

    x3_candidate = sp.sqrt(T) * g2 * a1 / (2 * loop)
    y3 = T * g2 * a2 / (2 * loop)
    a1_shifted = a1 + 2 * a2 * delta
    delta_x = sp.simplify(
        sp.sqrt(T) * g2 * a1_shifted / (2 * loop) - x3_candidate
    )
    assert sp.simplify(delta_x - 2 * y3 * delta3) == 0
