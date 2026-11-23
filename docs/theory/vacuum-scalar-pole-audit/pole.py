"""NST finite MS subset with only massive physical scalars inside each loop."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parent


def load(): return json.loads((ROOT/"inputs.json").read_text())


def tree(cfg):
    L,M,s,v,a,b,q=(sp.Rational(cfg[k]) for k in ("m1","m2","sin_journal","v_tree","a2","b3","b4"))
    c=sp.sqrt(1-s*s)
    w,z=sp.symbols("h1 h2",real=True)
    h=c*w-s*z
    S=s*w+c*z
    X=(v+h)**2/2
    lam=(c*c*L*L+s*s*M*M)/(2*v*v)
    a1=2*c*s*(L*L-M*M)/v
    b2=s*s*L*L+c*c*M*M-a*v*v/2
    P=-lam*v*v*X+lam*X*X+(a1*S+a*S*S)*X/2-a1*v*v*S/4+b2*S*S/2+b*S**3/3+q*S**4/4
    P=sp.Poly(P,w,z).as_expr()
    Z=(w,z); O={w:0,z:0}
    H=sp.hessian(P,Z)
    assert H.subs(O)==sp.diag(L*L,M*M)
    assert all(sp.diff(P,e).subs(O)==0 for e in Z)
    D=lambda inds:sp.simplify(sp.diff(P,*[Z[i] for i in inds]).subs(O))
    t3={a:D(a) for a in itertools.product(range(2),repeat=3)}
    t4={a:D(a) for a in itertools.product(range(2),repeat=4)}
    n=-s; cos2=c*c-n*n; sin2=2*c*n; cos4=1-2*sin2*sin2
    # Literal NST Eqs. (A45)-(A51), vertices c=-grad(V), NST angle -theta_RMT.
    source={
       "111":-3*L*L/v*c**3-3*v*a*c*n*n+2*b*n**3,
       "222":-3*M*M/v*n**3-3*v*a*n*c*c-2*b*c**3,
       "112":n/(2*v)*(v*v*a-2*L*L-M*M+(3*v*v*a-2*L*L-M*M)*cos2)-2*b*c*n*n,
       "122":c/(2*v)*(v*v*a-L*L-2*M*M-(3*v*v*a-L*L-2*M*M)*cos2)+2*b*n*c*c}
    s4={
       "1111":-3*c**4/(2*v*v)*(L*L+M*M+(L*L-M*M)*cos2)-6*(q*n*n+a*c*c)*n*n,
       "2222":-3*n**4/(2*v*v)*(L*L+M*M+(L*L-M*M)*cos2)-6*(q*c*c+a*n*n)*c*c,
       "1122":-3*sin2**2/(8*v*v)*(L*L+M*M+(L*L-M*M)*cos2)-(3*q+a+3*(a-q)*cos4)/4}
    for dic, tt in ((source,t3),(s4,t4)):
        assert all(sp.simplify(x+tt[tuple(int(a)-1 for a in k)])==0 for k,x in dic.items())
    m=asmp
    return {"mass":list(map(m,(L,M))),"s":m(s),"c":m(c),"v":m(v),"a2":m(a),"lam":m(lam),"a1":m(a1),"b1":m(-a1*v*v/4),
        "g3":{k:m(x) for k,x in t3.items()},"g4":{k:m(x) for k,x in t4.items()},
        "nst_signed_vertices":{k:m(y) for k,y in (source|s4).items()}}


def asmp(x): return mp.mpf(str(sp.N(x,mp.mp.dps+4)))


def a_ms(m2,mu):
    return m2*(1+mp.log(mu*mu/m2))/(16*mp.pi**2)


def roots(k,a,b):
    if k==0: return []
    w=a-b-k
    d=mp.sqrt(mp.mpc(w*w-4*k*b))
    return [(-w+d)/(2*k),(-w-d)/(2*k)]


def b_ms(k,a,b,mu,method):
    rr=roots(k,a,b)
    real=sorted(mp.re(t) for t in rr if mp.im(t)==0 and 0<mp.re(t)<1)
    cuts=[mp.mpf(0)]+real+[mp.mpf(1)]
    # timelike D-i0: + imaginary part of -log(D-i0).
    D=lambda u:(1-u)*b+u*a-k*u*(1-u)
    negwidth=sum(y-x for x,y in zip(cuts,cuts[1:]) if D((x+y)/2)<0)
    if method=="integral":
        # Factor at real roots to avoid cancellation near logarithmic endpoints.
        def fun(u):
            d=D(u) if len(real)!=2 else k*(u-real[0])*(u-real[1])
            return -mp.log(abs(d)/mu**2) if d else mp.mpf(0)
        re=mp.quad(fun,cuts)
    elif method=="roots":
        if k==0:
            if a==b: re=-mp.log(a/mu**2)
            else:
                F=lambda u:u*(mp.log(u/mu**2)-1)
                re=-(F(a)-F(b))/(a-b)
        else:
            F=lambda u:u*(mp.log(u)-1) if u else 0
            re=-mp.log(k/mu**2)-sum(mp.re(F(1-r)-F(-r)) for r in rr)
    else: raise ValueError(method)
    return (re+1j*mp.pi*negwidth)/(16*mp.pi**2)


def a28_explicit(ew,mu):
    GF,W,Z=(mp.mpf(ew[k]) for k in ("GF","MW","MZ"))
    g0sq=8*W*W*GF/mp.sqrt(2)
    bracket=4*mp.log(mu**2/W**2)+(mp.mpf(7)/2*Z**2/(Z**2-W**2)-2)*mp.log(W**2/Z**2)+6
    return {"g0sq":g0sq,"v_GF_tree_different_from_fixture_GeV":2*W/mp.sqrt(g0sq),
            "explicit_bracket":bracket,"signed_explicit_bracket_contribution":-g0sq*bracket/(16*mp.pi**2),
            "PiW_momentum_difference_uncomputed":True}


def calculate(t,mu):
    m2=[m*m for m in t["mass"]]; g=t["g3"]; q=t["g4"]
    A=[a_ms(z,mu) for z in m2]
    t_cw=[-sum(g[i,a,a]*A[a] for a in range(2))/2 for i in range(2)]
    vev_shift=[-t_cw[i]/m2[i] for i in range(2)]
    def CW(x,y):
        z=(x,y)
        H=mp.matrix([[m2[i]*(i==j)+sum(g[i,j,a]*z[a] for a in range(2))
              +sum(q[i,j,a,b]*z[a]*z[b] for a,b in itertools.product(range(2),repeat=2))/2
              for j in range(2)] for i in range(2)])
        w=mp.eigsy(H,eigvals_only=True)
        return sum(L**2*(mp.log(L/mu**2)-mp.mpf("1.5")) for L in w)/(64*mp.pi**2)

    rows=[]
    errors=[]
    for i in range(2):
        B={}
        for j,k in ((0,0),(0,1),(1,1)):
           z=b_ms(m2[i],m2[j],m2[k],mu,"integral")
           zc=b_ms(m2[i],m2[j],m2[k],mu,"roots")
           errors.append(abs(z-zc))
           B[(j,k)]=z
           B[(k,j)]=z
        bubble=sum(g[i,j,k]**2*B[j,k] for j,k in itertools.product(range(2),repeat=2))/2
        seagull=sum(q[i,i,j,j]*A[j] for j in range(2))/2
        reducible=-sum(g[i,i,j]*vev_shift[j] for j in range(2))
        zero=seagull+sum(g[i,j,k]**2*b_ms(mp.mpf(0),m2[j],m2[k],mu,"roots") for j,k in itertools.product(range(2),repeat=2))/2
        fn=(lambda x:CW(x,0)) if i==0 else (lambda x:CW(0,x))
        dx=mp.mpf("1e-13")
        cw1=(fn(-2*dx)-8*fn(-dx)+8*fn(dx)-fn(2*dx))/(12*dx)
        cw2=(-fn(-2*dx)+16*fn(-dx)-30*fn(0)+16*fn(dx)-fn(2*dx))/(12*dx*dx)
        e1=abs(cw1-t_cw[i])/max(1,abs(t_cw[i]))
        e2=abs(cw2+mp.re(zero))/max(1,abs(zero))
        rows.append({"state":i+1,"pole_p2":m2[i],"bubbles":bubble,"quartic_seagull":seagull,
                     "1PI":bubble+seagull,"1PR_tadpole":reducible,"Pi_partial":bubble+seagull+reducible,
                     "CW_gradient_scaled_error":e1,"CW_curvature_scaled_error":e2,
                     "B0_integral_by_index":{f"{j+1}{k+1}":B[j,k] for j,k in ((0,0),(0,1),(1,1))}})
    P,Q=(mp.re(x["Pi_partial"]) for x in rows)
    c,s,v=t["c"],t["s"],t["v"]
    X=c*c*P+s*s*Q;Y=s*s*P+c*c*Q;A1=2*c*s*(P-Q)/v
    return {"mu_bar":mu,"tadpole_cw_gradient_mass_eigenbasis":t_cw,"local_CW_minimum_shift_not_EW_v_matching":vev_shift,
            "rows":rows,"max_B0_integral_root_difference":max(errors),
            "diagonal_scalar_pole_block_response":{"dcH":-X/2,"db2":Y,"da1":A1,"db1":-v*v*A1/4,"dlambda":X/(2*v*v)}}


def pack(x):
    if isinstance(x,dict):return {str(k):pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [pack(v) for v in x]
    if isinstance(x,mp.mpc):return {"re":pack(x.real),"im":pack(x.imag)}
    if isinstance(x,mp.mpf):return mp.nstr(x,61)
    return x


def result():
    r=load();t=tree(r["benchmark"])
    scans=[calculate(t,mp.mpf(q)) for q in r["mu_bar_GeV"]]
    a=a28_explicit(r["nst_A20_reference"],mp.mpf(r["mu_bar_GeV"][0]))
    f=lambda k:hashlib.sha256((ROOT/k).read_bytes()).hexdigest()
    return {"scope":"massive physical scalar internal h1,h2 subset only; not full pole-to-MS parameters",
            "inputs":r,"phase_angle_NST":-t["s"],"nst_signed_vertices_by_index":t["nst_signed_vertices"],
            "partial_mass_diagrams":scans,"a28_explicit_reference_not_full_EW":a,
            "linear_response_EW_dv_over_v_coefficient":{"dcH":0,"db2":-t["a2"]*t["v"]**2,
                "da1":-t["a1"],"db1":t["b1"],"dlambda":-2*t["lam"]},
            "runtime":{"python":platform.python_version(),"mpmath":mp.__version__,"sympy":sp.__version__,"digits":mp.mp.dps},
            "sha256":{k:f(k) for k in ("inputs.json","pole.py","test_pole.py")}}


def main():
    p=argparse.ArgumentParser();p.add_argument("out",type=Path);a=p.parse_args()
    mp.mp.dps=70
    r=result();z=pack(r)
    a.out.write_text(json.dumps(z,indent=2,sort_keys=True,allow_nan=False)+"\n")
    for d in r["partial_mass_diagrams"]:
        print("mu",mp.nstr(d["mu_bar"],8),"max B0 error",mp.nstr(d["max_B0_integral_root_difference"],7))
        for o in d["rows"]: print(o["state"],*[mp.nstr(o[k],12) for k in ("1PI","1PR_tadpole","Pi_partial","CW_gradient_scaled_error","CW_curvature_scaled_error")])
    print("A28 explicit only",z["a28_explicit_reference_not_full_EW"]["signed_explicit_bracket_contribution"])


if __name__=="__main__":main()
