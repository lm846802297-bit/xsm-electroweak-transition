import math

from xsm_ewpt.basis import XSMParameters, shift_singlet_basis
from xsm_ewpt.matching import C_THERMAL, C_THERMAL_DECIMAL, MatchingScales, match_scalar_sector_lo, match_scalar_sector_nlo, shift_3d_basis


def sample_parameters():
    # Coefficients are in GeV powers appropriate to the 4D potential.
    return XSMParameters(
        cH=-7800.0,
        lam=0.129,
        a1=210.0,
        a2=0.72,
        b1=0.0,
        b2=8200.0,
        b3=38.0,
        b4=0.28,
    )


def sample_scales():
    return MatchingScales(
        temperature=95.0,
        parent_scale=173.0,
        eft_scale=95.0,
        g2=0.425,
        gp2=0.127,
        yt2=0.94,
    )


def assert_close_fields(x, y, rtol=2e-8, atol=2e-8):
    for name in x.__dataclass_fields__:
        xv, yv = getattr(x, name), getattr(y, name)
        assert math.isclose(xv, yv, rel_tol=rtol, abs_tol=atol), (name, xv, yv, xv - yv)


def test_lo_tadpole_matches_lx_high_temperature_formula():
    p = sample_parameters()
    s = sample_scales()
    q = match_scalar_sector_lo(p, s)
    expected_lo = (p.b1 + s.temperature**2 * (p.a1 + p.b3) / 12.0) / math.sqrt(s.temperature)
    assert math.isclose(q.b1_3, expected_lo, rel_tol=1e-14, abs_tol=1e-12)
    assert math.isclose(q.b2_3, p.b2 + s.temperature**2 * (p.a2 / 6.0 + p.b4 / 4.0), rel_tol=1e-14)


def test_lo_matching_is_covariant_under_singlet_translation():
    p = sample_parameters()
    s = sample_scales()
    delta = -41.0
    path_a = shift_3d_basis(match_scalar_sector_lo(p, s), delta / math.sqrt(s.temperature))
    path_b = match_scalar_sector_lo(shift_singlet_basis(p, delta), s)
    assert_close_fields(path_a, path_b, rtol=2e-13, atol=2e-11)


def test_nlo_matching_is_covariant_under_singlet_translation():
    p = sample_parameters()
    s = sample_scales()
    delta = 73.0
    p_shift = shift_singlet_basis(p, delta)

    q_then_shift = shift_3d_basis(match_scalar_sector_nlo(p, s), delta / math.sqrt(s.temperature))
    q_shift_then_match = match_scalar_sector_nlo(p_shift, s)

    assert_close_fields(q_then_shift, q_shift_then_match, rtol=4e-7, atol=4e-7)


def test_nlo_matching_covariance_over_deterministic_parameter_ensemble():
    import random

    rng = random.Random(20261123)
    for _ in range(100):
        p = XSMParameters(
            cH=rng.uniform(-2.0e4, -2.0e3),
            lam=rng.uniform(0.08, 0.25),
            a1=rng.uniform(-300.0, 300.0),
            a2=rng.uniform(-0.5, 1.5),
            b1=rng.uniform(-1.0e5, 1.0e5),
            b2=rng.uniform(-2.0e4, 2.0e4),
            b3=rng.uniform(-100.0, 100.0),
            b4=rng.uniform(0.05, 0.8),
        )
        temperature = rng.uniform(50.0, 150.0)
        scales = MatchingScales(
            temperature=temperature,
            parent_scale=rng.uniform(100.0, 300.0),
            eft_scale=rng.uniform(50.0, 200.0),
            g2=0.42,
            gp2=0.13,
            yt2=0.94,
        )
        delta = rng.uniform(-150.0, 150.0)
        matched_then_shifted = shift_3d_basis(
            match_scalar_sector_nlo(p, scales), delta / math.sqrt(temperature)
        )
        shifted_then_matched = match_scalar_sector_nlo(
            shift_singlet_basis(p, delta), scales
        )
        assert_close_fields(
            matched_then_shifted,
            shifted_then_matched,
            rtol=2e-11,
            atol=2e-10,
        )


def test_thermal_constant_matches_zeta_definition_at_high_precision():
    import mpmath as mp

    mp.mp.dps = 60
    expected = mp.mpf("0.5") * (
        mp.log(8 * mp.pi / 9)
        + mp.diff(mp.zeta, 2) / mp.zeta(2)
        - 2 * mp.euler
    )
    assert abs(mp.mpf(C_THERMAL_DECIMAL) - expected) < mp.mpf("1e-50")
    assert abs(mp.mpf(str(C_THERMAL)) - expected) < mp.mpf("1e-16")
