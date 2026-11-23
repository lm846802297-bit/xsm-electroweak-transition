"""Independent positive-mass scalar thermal tadpole audit; not a resummation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import k1, zeta

ROOT=Path(__file__).resolve().parent
pi=math.pi
gamma=float(np.euler_gamma)
F0=1/12
ZETA3=float(zeta(3,1))


def load():
    return json.loads((ROOT/"inputs.json").read_text())


def checked(x):
    x=float(x)
    if not math.isfinite(x) or x<0:
        raise ValueError("m/T must be finite and nonnegative; no tachyonic continuation")
    return x


def radial(x):
    x=checked(x)
    if x==0: return F0,0.
    def f(q):
        E=math.hypot(q,x)
        return q*q/E*math.exp(-E)/(-math.expm1(-E))/(2*pi*pi)
    return quad(f,0,math.inf,epsabs=1e-13,epsrel=1e-13,limit=200)


def bessel(x,cut=44):
    x=checked(x)
    if x==0: return F0,0.,0
    N=max(1,math.ceil(cut/x))
    n=np.arange(1,N+1,dtype=float)
    terms=k1(n*x)/n
    total=x/(2*pi*pi)*float(np.sum(terms))
    # K1((N+j)x)<=exp(-jx)K1(Nx), from its positive cosh integral.
    tail=x/(2*pi*pi)*float(terms[-1])*math.exp(-x)/(-math.expm1(-x))
    return total,tail,N


def highT2(x):
    x=checked(x)
    if x==0: return F0
    return F0-x/(4*pi)+x*x/(8*pi*pi)*(math.log(4*pi/x)-gamma+.5)


def bound_R(x):
    x=checked(x)
    u=x**4*ZETA3/(128*pi**4)
    return u/(1+(x/(2*pi))**2),u


def paired_R(x,N=64):
    x=checked(x)
    if N<1: raise ValueError("include at least one explicit hard pair")
    w=2*pi*np.arange(1,N+1,dtype=float)
    # Rationalized (w + x²/(2w) - hypot(w,x))/(2*pi).
    partial=float(np.sum((x**4/(4*w*w))/(w+x*x/(2*w)+np.hypot(w,x))))/(2*pi)
    tail_hi=x**4*float(zeta(3,N+1))/(128*pi**4)
    tail_lo=tail_hi/(1+(x/(2*pi*(N+1)))**2)
    return partial+tail_lo,partial+tail_hi


def split(x,mu_over_T,N=64):
    x=checked(x)
    mu=float(mu_over_T)
    if not math.isfinite(mu) or mu<=0: raise ValueError("mu/T must be positive")
    coef=x*x/(16*pi*pi)
    hard2=F0-2*coef*(gamma+math.log(mu/(4*pi)))
    vac=coef*(2*math.log(x/mu)-1) if x else 0.
    zero=-x/(4*pi)
    Rlo,Rhi=paired_R(x,N)
    return {"mu_over_T":mu,"vacuum_ms":vac,"zero_matsubara_DR":zero,
            "hard_ms_interval":[hard2+Rlo,hard2+Rhi],
            "thermal_reconstructed_interval":[(hard2-vac+zero)+Rlo,(hard2-vac+zero)+Rhi]}


def example(row):
    T,a,b,x,y=map(float,[row[k] for k in ("T","alpha_H","alpha_S","x_H","x_S")])
    H,S=radial(x)[0],radial(y)[0]
    return dict(row,Higgs=a*T*T*H,singlet=b*T*T*S,
                thermal_scalar_total=T*T*(a*H+b*S),massless_eq_2_8_local=T*T*(a+b)*F0)


def main():
    cli=argparse.ArgumentParser();cli.add_argument("out",type=Path);a=cli.parse_args()
    cfg=load();rows=[]
    for x in cfg["x_scan"]:
        exact,e=radial(x); b,bb,N=bessel(x,cfg["bessel_cutoff_xn"])
        t=highT2(x);R0,R1=bound_R(x);r0,r1=paired_R(x,cfg["hard_explicit_pairs"])
        splits=[split(x,u,cfg["hard_explicit_pairs"]) for u in cfg["mu_over_T"]]
        lim=2e-13
        assert abs(exact-b)<=bb+lim
        assert R0-lim<=r0<=r1<=R1+lim
        for d in splits:
            l,h=d["thermal_reconstructed_interval"]
            assert l-lim<=exact<=h+lim
        rows.append(dict(x=x,F_thermal=exact,quad_reported_error=e,F_bessel_summed=b,
                         bessel_tail_bound=bb,bessel_terms=N,LO_ratio=exact/F0,
                         linear_soft_only=F0-x/(4*pi),HT_through_x2_log_and_constant=t,
                         rigorous_R_global=[R0,R1],rigorous_R_paired=[r0,r1],observed_R=exact-t,
                         split=splits,absolute_radial_bessel_difference=abs(exact-b)))
    root=brentq(lambda x:radial(x)[0]/F0-.9,*cfg["x_root_10pct_bracket"],xtol=1e-13)
    max_diff=max(r["absolute_radial_bessel_difference"] for r in rows)
    report=dict(scope="positive scalar mass H=0 static 4D thermal determinant; not full EFT or resummation",
                runtime={"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
                rows=rows,x_threshold_LO_fraction_0p9=root,examples=[example(e) for e in cfg["examples"]],
                max_absolute_radial_bessel_difference=max_diff,
                sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ("tadpole.py","inputs.json","test_tadpole.py")})
    a.out.write_text(json.dumps(report,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("x = m/T where exact/massless = 0.9:",root)
    print("max radial vs Bessel:",max_diff)
    for r in rows:
        print(f"x={r['x']:<5g} F/F0={r['LO_ratio']:.12f} R/F = {r['observed_R']/r['F_thermal']:.7g}")
    for v in report["examples"]: print(v["label"],v["thermal_scalar_total"],v["massless_eq_2_8_local"])


if __name__=="__main__": main()
