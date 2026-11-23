"""Version-pinned real xSM tree benchmark and scalar tensor input contract."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as S

ROOT=Path(__file__).resolve().parent
f=S.symbols("r1 r2 r3 r4 s",real=True)
r1,r2,r3,r4,s=f
origin=dict.fromkeys(f,0)
cH,lam,a1,a2,b1,b2,b3,b4=S.symbols("cH lam a1 a2 b1 b2 b3 b4",real=True)
pars=(cH,lam,a1,a2,b1,b2,b3,b4)
rho=sum(h*h for h in f[:4])/2
V=cH*rho+lam*rho**2+(a1*s+a2*s*s)*rho/2+b1*s+b2*s*s/2+b3*s**3/3+b4*s**4/4


def inputs():
    return json.loads((ROOT/"inputs.json").read_text())


def setup():
    x=inputs()
    m1,m2,st,v,aa2,bb3,bb4=[S.Rational(x[k]) for k in ("mh1_GeV","mh2_GeV","sin_theta","v_h_GeV_tree_reference","a2","b3_GeV","b4")]
    ct=S.sqrt(1-st*st)
    R=S.Matrix([[ct,st],[-st,ct]])
    D=S.diag(m1*m1,m2*m2)
    # Construct two literally printed sets independently from the Hessian.
    def params(h,hS,pp):
        return dict(zip(pars,(-h/2,h/(2*v*v),pp,aa2,-pp*v*v/4,hS-aa2*v*v/2,bb3,bb4)))
    Hh=m1*m1*ct*ct+m2*m2*st*st
    Hs=m2*m2*ct*ct+m1*m1*st*st
    p=2*st*ct*(m1*m1-m2*m2)/v
    return R,D,v,params(Hh,Hs,p),params(Hs,Hh,-p)


def tensor_rows(p):
    axes={
      "G1_s":(s,), "G2_HH":(r1,r1), "G2_SS":(s,s),
      "G3_HHS":(r1,r1,s),"G3_SSS":(s,s,s),
      "G4_HHHH":(r1,r1,r1,r1),"G4_HHKK":(r1,r1,r2,r2),
      "G4_HHSS":(r1,r1,s,s),"G4_SSSS":(s,s,s,s),
      "G3_HKS_zero":(r1,r2,s)}
    return {k:S.diff(V,*a).subs(origin).subs(p) for k,a in axes.items()}


def engineering3(p,T):
    return dict(zip(pars,[p[k]*w for k,w in zip(pars,(1,T,S.sqrt(T),T,1/S.sqrt(T),1,S.sqrt(T),T))]))


def compute():
    R,D,v,p,q=setup()
    vv={r1:0,r2:0,r3:0,r4:v,s:0}
    G=S.hessian(V,(r4,s)).subs(vv)
    H=G.subs(p)
    P=G.subs(q)
    assert (R*H*R.T-D).applyfunc(S.simplify).is_zero_matrix
    assert (R*P*R.T-S.diag(D[1,1],D[0,0])).applyfunc(S.simplify).is_zero_matrix
    assert (R*H*R.T-R*P*R.T).applyfunc(S.simplify).is_zero_matrix is False
    assert all(S.simplify(S.diff(V,x).subs(vv).subs(p))==0 for x in f)
    u,w=S.symbols("mass_state_1 mass_state_2")
    hfluct,sfluct=R.T*S.Matrix([u,w])
    mass_cubic=S.diff(V.subs(p).subs({r1:0,r2:0,r3:0,r4:v+hfluct,s:sfluct},simultaneous=True),u,u,w).subs({u:0,w:0}).simplify()
    st,ct=R[0,1],R[0,0]
    printed=2*p[b3]*st**2*ct+p[a1]*ct*(ct**2-2*st**2)/2+(2*ct**2-st**2)*st*v*p[a2]-6*p[lam]*st*ct**2*v
    assert S.simplify(mass_cubic-printed)==0
    bareH,bareS=(S.diff(V,k,k).subs(dict.fromkeys(f[:4],0)).subs(p) for k in (r1,s))
    roots=sorted(S.solve(bareH,s),key=float)
    s_stationary=S.solve(S.diff(bareS,s),s)[0]
    minS=bareS.subs(s,s_stationary)
    assert p[a2]>0 and p[b4]>0 and minS>0
    raw=inputs()
    m3=[]
    for T0 in raw["T_engineering_map_GeV"]:
        T=S.Rational(T0)
        p3=engineering3(p,T)
        scaled=V.subs(p).subs({x:S.sqrt(T)*x for x in f},simultaneous=True)/T
        assert S.expand(V.subs(p3)-scaled)==0
        m3.append(dict(T=T,scalar_coefficients={str(k):z for k,z in p3.items()},tensor_entries=tensor_rows(p3)))
    sha=lambda f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()
    return dict(inputs=raw,scope="paper tree input only and engineering 3D normalization, not physical-to-MSbar loop conversion, thermal DR or DRalgo execution",
      runtime={"python":platform.python_version(),"sympy":S.__version__},
      input_and_source_sha256={k:sha(k) for k in ("inputs.json","benchmark.py","test_benchmark.py")},
      journal_4D={str(k):z for k,z in p.items()},arxiv_v1_literal_4D={str(k):z for k,z in q.items()},
      journal_neutral_hessian=H,arxiv_v1_literal_neutral_hessian=P,
      journal_diag_with_eq_2_3=R*H*R.T,arxiv_v1_diag_with_eq_2_3=R*P*R.T,
      physically_light_state_doublet_purity_journal=ct*ct,physically_light_state_doublet_purity_arxiv_literal=st*st,
      lambda211_tree_from_potential_GeV=mass_cubic,tree_grad_tensors=tensor_rows(p),
      H0_bare_positive_domain={"doublet_curvature_at_s0":bareH.subs(s,0),"singlet_curvature_at_s0":bareS.subs(s,0),
                              "doublet_zeroes_s":roots,"singlet_min_curvature":minS,"singlet_min_at_s":s_stationary},
      engineering_3D_not_matched=m3)


def pack(x):
    if isinstance(x,S.MatrixBase):return [[pack(y) for y in row] for row in x.tolist()]
    if isinstance(x,dict):return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [pack(v) for v in x]
    if isinstance(x,S.Expr):return str(S.N(S.simplify(x),45))
    return x


def main():
    p=argparse.ArgumentParser();p.add_argument("out",type=Path);a=p.parse_args()
    r=pack(compute())
    a.out.write_text(json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+"\n")
    for k in ("journal_4D","arxiv_v1_literal_4D","journal_diag_with_eq_2_3","arxiv_v1_diag_with_eq_2_3","lambda211_tree_from_potential_GeV","H0_bare_positive_domain"):
        print(k,r[k])


if __name__=="__main__":main()
