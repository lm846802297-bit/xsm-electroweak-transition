from xsm_ewpt.basis import XSMParameters, shift_singlet_basis, stationarity_residual, higgs_stationarity_residual

def test_shift_roundtrip():
    p=XSMParameters(-12000,.13,120,.4,0,4000,30,.2)
    q=shift_singlet_basis(shift_singlet_basis(p,75),-75)
    for name in p.__dataclass_fields__:
        assert abs(getattr(p,name)-getattr(q,name)) < 1e-9

def test_stationarity_after_map():
    v=246.; vs=35.; a2=.4; b3=20.; b4=.2; lam=.13
    a1=-2*a2*vs
    b2=-(.25*a1*v*v+b3*vs*vs+b4*vs**3)/vs
    b1=0.
    cH=-lam*v*v-.5*a1*vs-.5*a2*vs*vs
    p=XSMParameters(cH,lam,a1,a2,b1,b2,b3,b4)
    q=shift_singlet_basis(p,vs)
    assert abs(stationarity_residual(q,v,0)) < 1e-8
    assert abs(higgs_stationarity_residual(q,v,0)) < 1e-8
