import sympy as sp

def test_thermal_tadpole_coefficient():
    a1,b3,T=sp.symbols("a1 b3 T")
    assert sp.expand((a1+b3)*T**2/12).coeff(a1) == T**2/12
