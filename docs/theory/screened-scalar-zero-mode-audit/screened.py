"""Static scalar zero-mode replacement, not complete xSM matching."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import mpmath as mp

ROOT=Path(__file__).resolve().parent


def load():
    return json.loads((ROOT/"fixtures.json").read_text())


def evaluate(mode,s,kind):
    X0,X1,X2,P0,P1,P2=(mp.mpf(mode[k]) for k in ("X0","X1","X2","P0","P1","P2"))
    if kind=="none": P0=P1=P2=mp.mpf(0)
    elif kind=="constant": P1=P2=mp.mpf(0)
    elif kind!="field": raise ValueError("unknown shift")
    X=X0+s*(X1+s*X2)
    P=P0+s*(P1+s*P2)
    if min(X,X+P)<=0: raise ValueError("both scalar squared masses must stay positive")
    return X,X1+2*s*X2,P,P1+2*s*P2


def th(T,X):
    v=X/T**2
    return T**4/(2*mp.pi**2)*mp.quad(
        lambda q:q*q*mp.log1p(-mp.exp(-mp.sqrt(q*q+v))),[0,1,5,mp.inf])


def It(T,X):
    v=X/T**2
    def integrand(q):
        e=mp.sqrt(q*q+v)
        return q*q/e/mp.expm1(e)
    return T**2/(2*mp.pi**2)*mp.quad(integrand,[0,1,5,mp.inf])


def z(T,X):
    if X<=0: raise ValueError("strict positive squared mass required")
    return -T*X*mp.sqrt(X)/(12*mp.pi)


def static_derivative_radial(T,X,D):
    # radial DR subtraction of a massless propagator at integrand level
    return -T*X*D/(4*mp.pi**2)*mp.quad(lambda p:1/(p*p+X),[0,mp.sqrt(X),mp.inf])


def direct_potential(T,mode,s,kind):
    X,_,P,_=evaluate(mode,s,kind)
    return mode["n"]*(th(T,X)-z(T,X)+z(T,X+P))


def fivepoint(fun,s,h):
    return (fun(s-2*h)-8*fun(s-h)+8*fun(s+h)-fun(s+2*h))/(12*h)


def record_mode(T,mode,s,kind):
    X,D,P,E=evaluate(mode,s,kind)
    pi=mp.pi
    n=mode["n"]; m=mp.sqrt(X); ms=mp.sqrt(X+P)
    dth=n*D*It(T,X)/2
    dZ0=-n*T*m*D/(8*pi)
    dZs=-n*T*ms*(D+E)/(8*pi)
    hard_minus_vac=dth-dZ0
    equivalent=hard_minus_vac+dZs
    radial=dth+n*(static_derivative_radial(T,X+P,D+E)-static_derivative_radial(T,X,D))
    numeric=fivepoint(lambda y:direct_potential(T,mode,y,kind),s,mp.mpf("1e-13"))
    double_count=dth+dZs
    frozen_slope=hard_minus_vac-n*T*ms*D/(8*pi)
    # Shifting the entire thermal function is shown only as a diagnostic
    # of the different mode content; it is not a full alternate scheme.
    all_mode_shift=n*(D+E)*It(T,X+P)/2
    first_insertion=-n*T*(m*E+P*D/(2*m))/(8*pi)
    return {"multiplicity":n,"M2":X,"M2_prime":D,"Pi":P,"Pi_prime":E,"Mscreened":ms,
            "thermal_unresummed":dth,"bare_static":dZ0,"screened_static":dZs,
            "thermal_nonzero_minus_vacuum":hard_minus_vac,
            "ring_only":dZs-dZ0,"screened_once":equivalent,
            "radial_DR_derivative":radial,"fivepoint_thermal_plus_ring_derivative":numeric,
            "wrong_add_without_subtract":double_count,
            "wrong_omit_Pi_prime":frozen_slope,
            "all_mode_thermal_only_diagnostic":all_mode_shift,
            "first_insertion_tadpole_if_already_in_fixed_order":first_insertion,
            "ring_remainder_after_first_insertion":dZs-dZ0-first_insertion}


def overlap_control(T,mode,s,eps_list):
    X,D,P,E=evaluate(mode,s,"field"); n=mode["n"]
    first=-n*T*(mp.sqrt(X)*E+P*D/(2*mp.sqrt(X)))/(8*mp.pi)
    direct=lambda e:-n*T*((D+e*E)*mp.sqrt(X+e*P)-D*mp.sqrt(X))/(8*mp.pi)
    numeric=mp.diff(direct,mp.mpf(0))
    rows=[]
    for e in map(mp.mpf,eps_list):
        q=direct(e); r=q-e*first
        rows.append({"epsilon":e,"ring_tadpole":q,"single_insertion":e*first,"subtract_insertion_remainder":r})
    return dict(coefficient=first,derivative_at_zero=numeric,rows=rows)


def scaled(a,b):
    return abs(a-b)/max(1,abs(a),abs(b))


def result():
    cfg=load();T,s=mp.mpf(cfg["T"]),mp.mpf(cfg["s"]);records=[]
    sumkeys=["thermal_unresummed","bare_static","screened_static","ring_only","screened_once",
             "wrong_add_without_subtract","wrong_omit_Pi_prime","all_mode_thermal_only_diagnostic"]
    for c in cfg["cases"]:
        modes={m["label"]:record_mode(T,m,s,c["Pi"]) for m in cfg["mode_data"]}
        maxrad=max(scaled(r["screened_once"],r["radial_DR_derivative"]) for r in modes.values())
        maxdiff=max(scaled(r["screened_once"],r["fivepoint_thermal_plus_ring_derivative"]) for r in modes.values())
        assert maxrad < mp.mpf("1e-59") and maxdiff < mp.mpf("1e-43")
        totals={k:sum(r[k] for r in modes.values()) for k in sumkeys}
        records.append(dict(case=c,modes=modes,totals=totals,scaled_radial=maxrad,scaled_fivepoint=maxdiff))
    oc={m["label"]:overlap_control(T,m,s,cfg["insertion_eps"]) for m in cfg["mode_data"]}
    f=lambda n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest()
    return dict(scope="positive masses, H=0, static one-loop scalar zero-mode replacement only",
                runtime={"python":platform.python_version(),"mpmath":mp.__version__,"digits":mp.mp.dps},
                sha256={k:f(k) for k in ("fixtures.json","screened.py","test_screened.py")},
                records=records,overlap_control=oc)


def pack(x):
    if isinstance(x,dict):return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [pack(v) for v in x]
    return mp.nstr(x,61) if isinstance(x,mp.mpf) else x


def main():
    p=argparse.ArgumentParser();p.add_argument("out",type=Path);a=p.parse_args()
    mp.mp.dps=65
    r=result()
    for row in r["records"]:
        print(row["case"]["name"],"radial/finite-difference scaled",*[mp.nstr(row[k],8) for k in ("scaled_radial","scaled_fivepoint")])
        print({k:mp.nstr(v,13) for k,v in row["totals"].items()})
    a.out.write_text(json.dumps(pack(r),indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":main()
