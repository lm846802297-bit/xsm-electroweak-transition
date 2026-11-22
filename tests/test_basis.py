import math

from xsm_ewpt.basis import (
    XSMParameters,
    higgs_stationarity_residual,
    map_lx_to_zero_vev_basis,
    shift_singlet_basis,
    stationarity_residual,
)


def test_lx_to_zero_vev_map_satisfies_zero_vev_tadpole_and_higgs_stationarity():
    v = 246.0
    vs = 115.0
    lam = 0.129
    a1 = -210.0
    a2 = 1.1
    b3 = 45.0
    b4 = 0.35
    b2 = -(b3 * vs**2 + b4 * vs**3 + 0.25 * a1 * v**2 + 0.5 * a2 * v**2 * vs) / vs
    mu2_lx = lam * v**2 + 0.5 * vs * (a1 + a2 * vs)
    p_lx = XSMParameters(-mu2_lx, lam, a1, a2, 0.0, b2, b3, b4)

    assert abs(stationarity_residual(p_lx, v, vs)) < 1e-8
    assert abs(higgs_stationarity_residual(p_lx, v, vs)) < 1e-8

    p_zero = map_lx_to_zero_vev_basis(p_lx, vs)
    assert math.isclose(p_zero.b1, -0.25 * p_zero.a1 * v**2, rel_tol=1e-12, abs_tol=1e-8)
    assert math.isclose(p_zero.cH, -lam * v**2, rel_tol=1e-12, abs_tol=1e-8)
    assert abs(stationarity_residual(p_zero, v, 0.0)) < 1e-8
    assert abs(higgs_stationarity_residual(p_zero, v, 0.0)) < 1e-8


def test_shift_round_trip():
    p = XSMParameters(-7800.0, 0.13, 210.0, 0.7, 1000.0, 8000.0, 35.0, 0.25)
    p2 = shift_singlet_basis(shift_singlet_basis(p, 83.0), -83.0)
    for name in p.__dataclass_fields__:
        assert math.isclose(getattr(p, name), getattr(p2, name), rel_tol=2e-13, abs_tol=2e-10)


def test_potential_identity_symbolically():
    import sympy as sp

    X, S, d = sp.symbols("X S d")
    cH, lam, a1, a2, b1, b2, b3, b4 = sp.symbols(
        "cH lam a1 a2 b1 b2 b3 b4"
    )
    old = (
        cH * X
        + lam * X**2
        + a1 * X * (S + d) / 2
        + a2 * X * (S + d) ** 2 / 2
        + b1 * (S + d)
        + b2 * (S + d) ** 2 / 2
        + b3 * (S + d) ** 3 / 3
        + b4 * (S + d) ** 4 / 4
    )
    cp = cH + a1 * d / 2 + a2 * d**2 / 2
    a1p = a1 + 2 * a2 * d
    b1p = b1 + b2 * d + b3 * d**2 + b4 * d**3
    b2p = b2 + 2 * b3 * d + 3 * b4 * d**2
    b3p = b3 + 3 * b4 * d
    new = (
        cp * X
        + lam * X**2
        + a1p * X * S / 2
        + a2 * X * S**2 / 2
        + b1p * S
        + b2p * S**2 / 2
        + b3p * S**3 / 3
        + b4 * S**4 / 4
    )
    offset = sp.expand(old - new)
    assert sp.diff(offset, X) == 0
    assert sp.diff(offset, S) == 0
