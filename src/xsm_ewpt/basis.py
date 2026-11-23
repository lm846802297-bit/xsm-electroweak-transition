"""Exact singlet-field translations for the general renormalizable xSM potential.

We write the Higgs-quadratic term as cH * (H^dagger H), so cH is the
coefficient itself rather than a convention-dependent positive mu^2.
The remaining potential is

 V = cH X + lam X^2 + (a1/2) X S + (a2/2) X S^2
     + b1 S + (b2/2) S^2 + (b3/3) S^3 + (b4/4) S^4,
 X = H^dagger H.

A translation S_old = S_new + delta induces the parameter map below.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class XSMParameters:
    cH: float
    lam: float
    a1: float
    a2: float
    b1: float
    b2: float
    b3: float
    b4: float


def shift_singlet_basis(p: XSMParameters, delta: float) -> XSMParameters:
    """Return parameters after substituting S_old = S_new + delta."""
    return XSMParameters(
        cH=p.cH + 0.5 * p.a1 * delta + 0.5 * p.a2 * delta**2,
        lam=p.lam,
        a1=p.a1 + 2.0 * p.a2 * delta,
        a2=p.a2,
        b1=p.b1 + p.b2 * delta + p.b3 * delta**2 + p.b4 * delta**3,
        b2=p.b2 + 2.0 * p.b3 * delta + 3.0 * p.b4 * delta**2,
        b3=p.b3 + 3.0 * p.b4 * delta,
        b4=p.b4,
    )


def map_lx_to_zero_vev_basis(p_lx: XSMParameters, v_s: float) -> XSMParameters:
    """Map Liu--Xie (b1=0, <S>=v_s) to the basis with <S>=0.

    Liu--Xie use V ⊃ -mu2_LX H^dagger H, while the 2026 paper writes
    V ⊃ +mu2_RMT H^dagger H. Therefore cH=-mu2_LX=mu2_RMT after
    the field translation. The caller supplies p_lx.cH=-mu2_LX.
    """
    if abs(p_lx.b1) > 1e-10:
        raise ValueError("The Liu--Xie input basis must satisfy b1=0")
    return shift_singlet_basis(p_lx, v_s)


def stationarity_residual(p: XSMParameters, v: float, v_s: float) -> float:
    """dV/dS at H^dagger H=v^2/2 and S=v_s."""
    return (
        p.b1
        + p.b2 * v_s
        + p.b3 * v_s**2
        + p.b4 * v_s**3
        + 0.25 * p.a1 * v**2
        + 0.5 * p.a2 * v**2 * v_s
    )


def higgs_stationarity_residual(p: XSMParameters, v: float, v_s: float) -> float:
    """dV/d(H^dagger H) at H^dagger H=v^2/2 and S=v_s."""
    return p.cH + p.lam * v**2 + 0.5 * p.a1 * v_s + 0.5 * p.a2 * v_s**2
