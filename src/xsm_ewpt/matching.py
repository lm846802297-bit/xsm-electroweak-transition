"""LO/NLO high-temperature xSM scalar dimensional-reduction matching."""
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
    def __post_init__(self):
        if min(self.temperature,self.parent_scale,self.eft_scale)<=0: raise ValueError("Scales must be positive")
        if min(self.g2,self.gp2,self.yt2)<0: raise ValueError("Squared couplings must be non-negative")
    @property
    def Lb(self): return 2*log(self.parent_scale/self.temperature)-2*(log(4*pi)-EULER_GAMMA)
    @property
    def Lf(self): return self.Lb+4*log(2)
    @property
    def c(self): return C_THERMAL
    @property
    def eft_log(self): return self.c+log(3*self.temperature/self.eft_scale)

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
    T=s.temperature; r=sqrt(T)
    return ThreeDParameters(r*p.a1,T*p.a2,(p.b1+T*T*(p.a1+p.b3)/12)/r,
      p.b2+T*T*(p.a2/6+p.b4/4),r*p.b3,T*p.b4,T*s.g2,T*s.gp2)

def match_scalar_sector_nlo(p: XSMParameters, s: MatchingScales) -> ThreeDParameters:
    T=s.temperature; Lb=s.Lb; Lf=s.Lf; g2=s.g2; gp2=s.gp2; yt2=s.yt2
    loop=(4*pi)**2; r=sqrt(T)
    a13lo=r*p.a1; a23lo=T*p.a2; b33lo=r*p.b3; b43lo=T*p.b4
    g3lo=T*g2; gp3lo=T*gp2
    a13=r*p.a1+r/loop*(Lb*(.75*(3*g2+gp2)*p.a1-2*p.b3*p.a2-(6*p.lam+2*p.a2)*p.a1)-3*Lf*yt2*p.a1)
    a23=T*p.a2+T/loop*(Lb*(.75*(3*g2+gp2)-6*p.lam-2*p.a2-3*p.b4)*p.a2-3*Lf*yt2*p.a2)
    b13=(p.b1+T*T*(p.b3+p.a1)/12-Lb/loop*(p.a1*p.cH+p.b3*p.b2)+T*T/loop*(
      (2+3*Lb)/48*(3*g2+gp2)*p.a1-Lb/2*((p.lam+7*p.a2/12)*p.a1+(p.a2/3+1.5*p.b4)*p.b3)
      -.125*p.a1*yt2*(3*Lb-Lf)))/r-1/loop*(2*b33lo*b43lo-.5*a13lo*(3*g3lo+gp3lo-2*a23lo))*s.eft_log
    b23=(p.b2+T*T*(p.a2/6+p.b4/4)-Lb/loop*(2*p.b3**2+.5*p.a1**2+2*p.a2*p.cH+3*p.b4*p.b2)
      +T*T/loop*((2+3*Lb)/24*(3*g2+gp2)*p.a2-Lb*((p.lam+7*p.a2/12+.5*p.b4)*p.a2+2.25*p.b4**2)
      -.25*p.a2*yt2*(3*Lb-Lf))+1/loop*((3*g3lo+gp3lo)*a23lo-2*a23lo**2-6*b43lo**2)*s.eft_log)
    b33=r*p.b3-r/loop*3*Lb*(.5*p.a1*p.a2+3*p.b4*p.b3)
    b43=T*p.b4-T/loop*Lb*(p.a2**2+9*p.b4**2)
    return ThreeDParameters(a13,a23,b13,b23,b33,b43,g3lo,gp3lo)

def shift_3d_basis(q: ThreeDParameters, delta_3: float) -> ThreeDParameters:
    return ThreeDParameters(q.a1_3+2*q.a2_3*delta_3,q.a2_3,
      q.b1_3+q.b2_3*delta_3+q.b3_3*delta_3**2+q.b4_3*delta_3**3,
      q.b2_3+2*q.b3_3*delta_3+3*q.b4_3*delta_3**2,
      q.b3_3+3*q.b4_3*delta_3,q.b4_3,q.g3_sq,q.gp3_sq)
