"""LO and NLO high-temperature dimensional-reduction matching for xSM.

The NLO scalar matching relations are transcribed from Appendix B,
Eqs. (B.5)--(B.10), of Ramsey-Musolf, Tenkanen and Tran,
JHEP 09 (2026) 260, arXiv:2409.17554.

The 3D scalar potential convention is
 V3 = 1/2 mu3_sq H3^dag H3 + 1/4 lam3 (H3^dag H3)^2
    + 1/4 a1_3 H3^dag H3 S3 + 1/4 a2_3 H3^dag H3 S3^2
    + b1_3 S3 + 1/2 b2_3 S3^2 + 1/3 b3_3 S3^3 + 1/4 b4_3 S3^4.

For the current covariance audit, mu3_sq and lam3 are not yet matched;
the scalar couplings entering the singlet-tadpole sector are implemented.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log, pi, sqrt

from .basis import XSMParameters


C_THERMAL_DECIMAL = "-0.348722736352290904004127649150998928954493426168099873381728682"
C_THERMAL = float(C_THERMAL_DECIMAL)
EULER_GAMMA = 0.5772156649015328606065120900824024


@dataclass(frozen=True)
class MatchingScales:
    temperature: float
    parent_scale: float
    eft_scale: float
    g2: float
    gp2: float
    yt2: float

    def __post_init__(self) -> None:
        if min(self.temperature, self.parent_scale, self.eft_scale) <= 0:
            raise ValueError("Temperature and renormalization scales must be positive")
        if min(self.g2, self.gp2, self.yt2) < 0:
            raise ValueError("Squared couplings must be non-negative")

    @property
    def Lb(self) -> float:
        return 2.0 * log(self.parent_scale / self.temperature) - 2.0 * (
            log(4.0 * pi) - EULER_GAMMA
        )

    @property
    def Lf(self) -> float:
        return self.Lb + 4.0 * log(2.0)

    @property
    def c(self) -> float:
        return C_THERMAL

    @property
    def eft_log(self) -> float:
        return self.c + log(3.0 * self.temperature / self.eft_scale)


@dataclass(frozen=True)
class ThreeDParameters:
    a1_3: float
    a2_3: float
    b1_3: float
    b2_3: float
    b3_3: float
    b4_3: float
    g3_sq: float
    gp3_sq: float


def match_scalar_sector_lo(p: XSMParameters, s: MatchingScales) -> ThreeDParameters:
    """Return leading-order hard-scale matching for the scalar sector."""
    T = s.temperature
    rootT = sqrt(T)
    return ThreeDParameters(
        a1_3=rootT * p.a1,
        a2_3=T * p.a2,
        b1_3=(p.b1 + T**2 * (p.a1 + p.b3) / 12.0) / rootT,
        b2_3=p.b2 + T**2 * (p.a2 / 6.0 + p.b4 / 4.0),
        b3_3=rootT * p.b3,
        b4_3=T * p.b4,
        g3_sq=T * s.g2,
        gp3_sq=T * s.gp2,
    )


def match_scalar_sector_nlo(p: XSMParameters, s: MatchingScales) -> ThreeDParameters:
    """Return NLO matched scalar-sector parameters (B.5)--(B.10).

    The formulas include the 3D logarithmic terms in b1_3 and b2_3.
    The 3D couplings appearing inside those explicitly two-loop terms are
    evaluated at LO, as required at this perturbative order.
    """
    T = s.temperature
    Lb = s.Lb
    Lf = s.Lf
    g2 = s.g2
    gp2 = s.gp2
    yt2 = s.yt2
    loop = (4.0 * pi) ** 2
    rootT = sqrt(T)

    # LO 3D couplings used in the explicitly two-loop 3D logarithms.
    a1_3_lo = rootT * p.a1
    a2_3_lo = T * p.a2
    b3_3_lo = rootT * p.b3
    b4_3_lo = T * p.b4
    g3_sq_lo = T * g2
    gp3_sq_lo = T * gp2

    a1_3 = rootT * p.a1 + rootT / loop * (
        Lb
        * (
            0.75 * (3.0 * g2 + gp2) * p.a1
            - 2.0 * p.b3 * p.a2
            - (6.0 * p.lam + 2.0 * p.a2) * p.a1
        )
        - 3.0 * Lf * yt2 * p.a1
    )

    a2_3 = T * p.a2 + T / loop * (
        Lb
        * (
            0.75 * (3.0 * g2 + gp2)
            - 6.0 * p.lam
            - 2.0 * p.a2
            - 3.0 * p.b4
        )
        * p.a2
        - 3.0 * Lf * yt2 * p.a2
    )

    b1_3 = (
        p.b1
        + T**2 * (p.b3 + p.a1) / 12.0
        - Lb / loop * (p.a1 * p.cH + p.b3 * p.b2)
        + T**2
        / loop
        * (
            ((2.0 + 3.0 * Lb) / 48.0) * (3.0 * g2 + gp2) * p.a1
            - Lb
            / 2.0
            * (
                (p.lam + 7.0 * p.a2 / 12.0) * p.a1
                + (p.a2 / 3.0 + 1.5 * p.b4) * p.b3
            )
            - 0.125 * p.a1 * yt2 * (3.0 * Lb - Lf)
        )
    ) / rootT - 1.0 / loop * (
        2.0 * b3_3_lo * b4_3_lo
        - 0.5
        * a1_3_lo
        * (3.0 * g3_sq_lo + gp3_sq_lo - 2.0 * a2_3_lo)
    ) * s.eft_log

    b2_3 = (
        p.b2
        + T**2 * (p.a2 / 6.0 + p.b4 / 4.0)
        - Lb / loop * (2.0 * p.b3**2 + 0.5 * p.a1**2 + 2.0 * p.a2 * p.cH + 3.0 * p.b4 * p.b2)
        + T**2
        / loop
        * (
            ((2.0 + 3.0 * Lb) / 24.0) * (3.0 * g2 + gp2) * p.a2
            - Lb
            * (
                (p.lam + 7.0 * p.a2 / 12.0 + 0.5 * p.b4) * p.a2
                + 2.25 * p.b4**2
            )
            - 0.25 * p.a2 * yt2 * (3.0 * Lb - Lf)
        )
    ) + 1.0 / loop * (
        (3.0 * g3_sq_lo + gp3_sq_lo) * a2_3_lo
        - 2.0 * a2_3_lo**2
        - 6.0 * b4_3_lo**2
    ) * s.eft_log

    b3_3 = rootT * p.b3 - rootT / loop * 3.0 * Lb * (
        0.5 * p.a1 * p.a2 + 3.0 * p.b4 * p.b3
    )

    b4_3 = T * p.b4 - T / loop * Lb * (p.a2**2 + 9.0 * p.b4**2)

    return ThreeDParameters(
        a1_3=a1_3,
        a2_3=a2_3,
        b1_3=b1_3,
        b2_3=b2_3,
        b3_3=b3_3,
        b4_3=b4_3,
        g3_sq=g3_sq_lo,
        gp3_sq=gp3_sq_lo,
    )


def shift_3d_basis(q: ThreeDParameters, delta_3: float) -> ThreeDParameters:
    """Apply S3_old = S3_new + delta_3 to the matched 3D potential."""
    return ThreeDParameters(
        a1_3=q.a1_3 + 2.0 * q.a2_3 * delta_3,
        a2_3=q.a2_3,
        b1_3=q.b1_3 + q.b2_3 * delta_3 + q.b3_3 * delta_3**2 + q.b4_3 * delta_3**3,
        b2_3=q.b2_3 + 2.0 * q.b3_3 * delta_3 + 3.0 * q.b4_3 * delta_3**2,
        b3_3=q.b3_3 + 3.0 * q.b4_3 * delta_3,
        b4_3=q.b4_3,
        g3_sq=q.g3_sq,
        gp3_sq=q.gp3_sq,
    )
