"""Static one-loop hard temporal xSM portal: DR + independent holonomy check."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import platform
import mpmath as mp
import scipy
from scipy.integrate import quad
from scipy.special import k1, zeta
import sympy as sp

REPO_HEAD='f365f6021e5a469c189aaf834ad7acf419dc6ab6'
MASS_POINTS=(.02,.05,.1,.2,.5,1.,2.,4.)

def symbolic():
    d,w,Q=sp.symbols('d w Q')
    # d-dimensional spatial IBP: sum/spatial integral over p²/(P²)^3
    # is d*I2/4. Thus omega²/(P²)^3 contributes (1-d/4)*I2.
    mass_factor=sp.expand(-1+4*(1-d/4))
    second_factor=sp.expand(2-12*(1-d/6))
    sx=sp.Matrix([[0,1],[1,0]])/2
    sy=sp.Matrix([[0,-sp.I],[sp.I,0]])/2
    sz=sp.diag(sp.Rational(1,2),sp.Rational(-1,2))
    y=sp.eye(2)/2
    ts=(sx,sy,sz)
    assert all(sp.trace(a*b)==(sp.Rational(1,2) if a==b else 0) for a in ts for b in ts)
    assert all(sp.trace(a*y)==0 for a in ts)
    rtT,T,a1,a2,dd,g=sp.symbols('rtT T a1 a2 delta_3 g',positive=True)
    zz=sp.symbols('zeta3')
    I3_T=2*zz/((4*sp.pi)**4*T**2)
    result={'first_mass_derivative_over_I2':str(mass_factor),'second_mass_derivative_over_I3':str(second_factor),
            'SU2_U1_mixed_trace':str(sp.trace(sz*y))}
    for name,G in [('SU2',sz),('U1',y)]:
        C=g**2*sp.trace(G*G)
        baseline=2/(4*sp.pi)**2
        x=C*baseline*(a1/2)*rtT
        yy=C*baseline*(a2/2)*T
        yy_a1_squared=C*sp.Rational(1,2)*(a1/2)**2*(-4)*I3_T*T
        rmtx=g**2*rtT*a1/(4*sp.pi)**2
        rmty=g**2*T*a2/(2*(4*sp.pi)**2)
        sct=g**2*T*(-zz*a1**2/(2*(4*sp.pi)**4*T**2))
        shift=(x.subs(a1,a1+2*a2*dd*rtT)-x-2*yy*dd).subs(rtT**2,T)
        result[name]={'generator_trace':str(sp.trace(G*G)),'x3':str(sp.simplify(x)),
                      'y3':str(sp.simplify(yy)),'y3_a1_squared':str(sp.simplify(yy_a1_squared)),
                      'RMT_x_over_derived_x':str(sp.simplify(rmtx/x)),
                      'RMT_y_over_derived_y':str(sp.simplify(rmty/yy)),
                      'Schicho_a1_squared_over_derived':str(sp.simplify(sct/yy_a1_squared)),
                      'shift_candidate_difference':str(sp.simplify(shift))}
    return result

def epsilon_limit(e,temperature=1.,mu=1.):
    # Sum-integral measure (e^gamma mu²/(4π))^ε ∫ d^dp/(2π)^d.
    t=mp.mpf(temperature); mu=mp.mpf(mu)
    d=3-2*e
    i2=2*t*(mp.e**mp.euler*mu**2/(4*mp.pi))**e
    i2*=mp.gamma(2-d/2)/(4*mp.pi)**(d/2)*(2*mp.pi*t)**(d-4)*mp.zeta(4-d)
    return 2*e*i2

def leading():return 1/(8*math.pi**2)
def next_mass_slope():return -8*float(zeta(3))/(4*math.pi)**4

def holonomy_quad(r):
    r=float(r)
    if not (r>0 and math.isfinite(r)):raise ValueError("M/T must be real and positive")
    def integrand(u):
        E=math.hypot(u,r)
        n=math.exp(-E)/(-math.expm1(-E))
        return u*u/E*n*(1+n)*(1+2*n)
    thermal=quad(integrand,0,math.inf,epsabs=1e-12,epsrel=5e-14,limit=500)[0]
    return 1/(8*math.pi*r)-thermal/(4*math.pi**2)

def holonomy_bessel(r):
    r=float(r)
    if not (r>0 and math.isfinite(r)):raise ValueError("M/T must be real and positive")
    limit=max(50,math.ceil(52/r))
    thermal=math.fsum(n*float(k1(n*r)) for n in range(1,limit+1))
    return 1/(8*math.pi*r)-r*thermal/(4*math.pi**2)

def main():
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args()
    with mp.workdps(70):
        eps=['1e-2','1e-4','1e-8','1e-16','1e-25']
        lim=mp.mpf(1)/(8*mp.pi**2)
        epsrows=[{'epsilon':e,'value':mp.nstr(v:=epsilon_limit(mp.mpf(e)),58),'difference':mp.nstr(v-lim,48)} for e in eps]
    rows=[]
    for r in MASS_POINTS:
        q,b=holonomy_quad(r),holonomy_bessel(r)
        rows.append({'M_over_T':r,'holonomy_quad':q,'holonomy_bessel':b,'abs_route_difference':abs(q-b),
                   'relative_to_massless':q/leading(),'slope_from_massless':(q-leading())/r**2,
                   'massless_plus_first_mass_correction':leading()+next_mass_slope()*r*r})
    out={'baseline_draft_head':REPO_HEAD,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'static one-loop hard doublet determinant in a constant singlet and temporal gauge background at H=0; no DRalgo runtime',
        'references':{'Liu_Xie':'https://arxiv.org/pdf/2101.10469#page=4',
                      'Schicho_Tenkanen_Osterman':'https://arxiv.org/pdf/2102.11145#page=33',
                      'Ramsey_Musolf_Tenkanen_Tran':'https://arxiv.org/pdf/2409.17554#page=26'},
        'symbolic':symbolic(),'epsilon_convergence':epsrows,'massless_hard_derivative':leading(),
        'leading_M_over_T_squared_coefficient':next_mass_slope(),
        'holonomy':rows,'max_abs_holonomy_route_difference':max(x['abs_route_difference'] for x in rows),
        'runtime':{'python':platform.python_version(),'scipy':scipy.__version__,'sympy':sp.__version__,'mpmath':mp.__version__}}
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(out['symbolic'])
    print('max independent holonomy route absolute difference',out['max_abs_holonomy_route_difference'])
if __name__=='__main__':main()
