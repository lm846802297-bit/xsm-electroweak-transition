from xsm_ewpt.basis import XSMParameters
from xsm_ewpt.matching import MatchingScales, match_scalar_sector_lo, match_scalar_sector_nlo, shift_3d_basis

def test_lo_tadpole():
    p=XSMParameters(-10000,.13,100,.4,0,1000,20,.2)
    s=MatchingScales(100,200,100,.42,.13,.94)
    q=match_scalar_sector_lo(p,s)
    assert abs(q.b1_3-(10000/12)/10) < 1e-12

def test_nlo_finite():
    p=XSMParameters(-10000,.13,100,.4,0,1000,20,.2)
    s=MatchingScales(100,200,100,.42,.13,.94)
    q=match_scalar_sector_nlo(p,s)
    assert all(abs(getattr(q,k))<1e10 for k in q.__dataclass_fields__)
    shifted=shift_3d_basis(q,0.3)
    assert shifted.a1_3 == q.a1_3+2*q.a2_3*.3
