"""NST Eq. (A58) massive top / massless bottom W fermion subset."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import mpmath as mp

ROOT=Path(__file__).resolve().parent


def load():return json.loads((ROOT/"inputs.json").read_text())


def sigma(p,M2,mu):
    d=mp.mpf(p)
    if M2<=0 or d>=M2: raise ValueError("require a massive top and p² below its threshold")
    return lambda x:x*(M2-(1-x)*d)


def loopfac():return 1/(16*mp.pi**2)


def a0(M2,mu):return M2*(1+mp.log(mu**2/M2))*loopfac() if M2 else mp.mpf(0)


def params(p,M2,mu):
    d=sigma(p,M2,mu)
    L=lambda x:mp.log(mu**2/d(x)) if x else mp.mpf(0)
    A=a0(M2,mu)
    # Feynman parameter x weights the second propagator (massive top).
    pts=[0,mp.mpf(".125"),1]
    B0=mp.quad(L,pts)*loopfac()
    B1=-mp.quad(lambda x:x*L(x),pts)*loopfac()
    twiceB00=mp.quad(lambda x:a0(d(x),mu),pts)
    k=A+p*B1  # NST A43: a=1, b=1, c=0, m1=0
    return {"A":A,"B0":B0,"B1":B1,"2B00":twiceB00,"K":k,
            "shape":twiceB00-k,"shape_UV_pole_coefficient":(-M2/2+p/3)*loopfac()}


def reduction(p,M2,mu):
    q=params(p,M2,mu)
    A=q["A"]; poleA=M2*loopfac(); poleB1=-loopfac()/2
    if p==0:
        B1PV=-(mp.log(mu**2/M2)/2+mp.mpf(1)/4)*loopfac()
    else:
        # A42 reconstructs B1 independently of the parameter-weighted route.
        B1PV=((M2-p)*q["B0"]-A)/(2*p)
    num=A-(M2-p)*B1PV
    res=poleA-(M2-p)*poleB1
    # Keep A41 as the common algebraic form at both momenta before the
    # intentionally incorrect comparison. A44 is a separate exact check.
    twiceB00=num/3+2*res/9
    n4=num/3
    k=A+p*B1PV
    if p==0:
        assert abs(twiceB00-(A/2+poleA/4))<mp.mpf("1e-59")
    return {"shape":twiceB00-k,"shape_if_D_set_to_4_too_early":n4-k,
            "B1_by_A42_when_nonzero":B1PV,"finite_epsilon_times_pole":twiceB00-n4}


def analytic_zero(M2,mu):
    return -M2*(mp.mpf(1)/4+mp.log(mu**2/M2)/2)*loopfac()


def analytic_zero_slope(M2,mu):
    return (mp.log(mu**2/M2)/3+mp.mpf(1)/9)*loopfac()


def analytic_b0(p,M2,mu):
    u=p/M2
    L=mp.log(mu**2/M2)
    return (1+L if p==0 else 2+L+(1-u)/u*mp.log1p(-u))*loopfac()


def calc(cfg,mu):
    if mp.mpf(cfg["mb_GeV"])!=0: raise ValueError("this A58 reduction takes mb=0")
    gf,W,t=(mp.mpf(cfg[a]) for a in ("GF_GeV_minus2","MW_GeV","Mt_GeV"))
    g0sq=8*W**2*gf/mp.sqrt(2)
    pref=g0sq*cfg["Nc"]*mp.mpf(cfg["Vtb"])**2
    m2=t*t;s=W*W
    P0,Ps=params(0,m2,mu),params(s,m2,mu)
    R0,Rs=reduction(0,m2,mu),reduction(s,m2,mu)
    W0,Ws=[pref*x["shape"] for x in (P0,Ps)]
    dg=(Ws-W0)/s
    dv=(Ws/s-dg)/2
    dy=dg-Ws/s
    # Secant derivative at zero is compared with the independent integrated
    # small-s formula; first neglected term here is O(h²).
    h=m2*mp.mpf("1e-22")
    num_sl=(params(h,m2,mu)["shape"]-params(-h,m2,mu)["shape"])/(2*h)
    maxerr=max(abs(P0["shape"]-R0["shape"]),abs(Ps["shape"]-Rs["shape"]))
    return {"mu_bar":mu,"g0sq_GF_tree":g0sq,
      "PiW_tb_MSfinite_at_zero_GeV2":W0,"PiW_tb_MSfinite_at_pole_GeV2":Ws,
      "PiW_tb_UV_pole_coeff_at_zero_GeV2":pref*P0["shape_UV_pole_coefficient"],
      "PiW_tb_UV_pole_coeff_at_pole_GeV2":pref*Ps["shape_UV_pole_coefficient"],
      "A28_tb_momentum_difference_only":dg,"dv_over_v_from_tb_W_input_only":dv,
      "A30_relative_y_t_squared_from_tb_W_input_only_no_top_self_energy":dy,
      "A28_tb_wrong_if_D4_before_UV_subtraction":pref*(Rs["shape_if_D_set_to_4_too_early"]-R0["shape_if_D_set_to_4_too_early"])/s,
      "zero_analytic_check_error":abs(P0["shape"]-analytic_zero(m2,mu)),
      "B0_closed_form_max_error":max(abs(P0["B0"]-analytic_b0(0,m2,mu)),abs(Ps["B0"]-analytic_b0(s,m2,mu))),
      "zero_slope_analytic_check_error":abs(num_sl-analytic_zero_slope(m2,mu)),
      "max_two_methods_abs_error":maxerr,"parametric":{"zero":P0,"pole":Ps},"PV_reduction":{"zero":R0,"pole":Rs}}


def output():
    cfg=load(); vals=[calc(cfg,mp.mpf(m)) for m in cfg["mu_bar_GeV"]]
    f=lambda n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest()
    return {"inputs":cfg,"rows":vals,
            "runtime":{"python":platform.python_version(),"mpmath":mp.__version__,"digits":mp.mp.dps},
            "sha256":{n:f(n) for n in ("inputs.json","wtb.py","test_wtb.py")}}


def pack(x):
    if isinstance(x,dict):return {k:pack(v) for k,v in x.items()}
    if isinstance(x,list):return [pack(v) for v in x]
    return mp.nstr(x,59) if isinstance(x,mp.mpf) else x


def main():
    p=argparse.ArgumentParser();p.add_argument("output",type=Path);a=p.parse_args()
    mp.mp.dps=65
    q=output();a.output.write_text(json.dumps(pack(q),indent=2,sort_keys=True,allow_nan=False)+"\n")
    for x in q["rows"]:
        print([mp.nstr(x[k],13) for k in ("mu_bar","PiW_tb_MSfinite_at_zero_GeV2","PiW_tb_MSfinite_at_pole_GeV2",
                    "A28_tb_momentum_difference_only","dv_over_v_from_tb_W_input_only","A30_relative_y_t_squared_from_tb_W_input_only_no_top_self_energy",
                    "A28_tb_wrong_if_D4_before_UV_subtraction","max_two_methods_abs_error","zero_slope_analytic_check_error")])


if __name__=="__main__":main()
