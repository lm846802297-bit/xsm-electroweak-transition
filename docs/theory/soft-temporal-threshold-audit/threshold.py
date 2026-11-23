"""One-loop static temporal threshold; no full tadpole or transition."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parent
WEIGHTS = (3,1,8)


def points():
    return json.loads((ROOT/"fixtures.json").read_text())


def matching(row, printed=False):
    T,g,gp,gs,yt,a1,a2=(mp.mpf(row[k]) for k in ("T","g","gp","gs","yt","a1_4","a2_4"))
    mass2=(mp.mpf(11)/6*g*g*T*T,mp.mpf(11)/6*gp*gp*T*T,2*gs*gs*T*T)
    if row.get("equal_debye"): assert mass2[0]==mass2[1]
    den=(4*mp.pi)**2
    x=tuple(g2*mp.sqrt(T)*a1/((1 if printed else 2)*den) for g2 in (g*g,gp*gp,0))
    y=tuple(g2*T*a2/(2*den) for g2 in (g*g,gp*gp,0))
    h=(g*g*T/4,gp*gp*T/4,-2*T*gs*gs*yt*yt/den)
    return {"mass2":mass2,"x":x,"y":y,"h":h,"k":g*gp*T/2}


def m2at(c,s):
    return tuple(m+2*x*s+2*y*s*s for m,x,y in zip(c["mass2"],c["x"],c["y"]))


def shifted(c,d):
    return dict(c,mass2=m2at(c,d),x=tuple(x+2*y*d for x,y in zip(c["x"],c["y"])))


def determinant(c,s,rho):
    a,b,z=(m+2*rho*h for m,h in zip(m2at(c,s),c["h"]))
    center=(a+b)/2
    gap=mp.sqrt(((a-b)/2)**2+(c["k"]*rho)**2)
    eig=(a,a,center+gap,center-gap)+(z,)*8
    if min(eig)<=0: raise ValueError("nonpositive temporal squared mass")
    return -sum(m*mp.sqrt(m) for m in eig)/(12*mp.pi)


def analytical(c,s):
    m=tuple(mp.sqrt(q) for q in m2at(c,s))
    u=tuple(x+2*y*s for x,y in zip(c["x"],c["y"]))
    y,h=c["y"],c["h"]
    W=lambda f:sum(n*f(M,U,Y,H) for n,M,U,Y,H in zip(WEIGHTS,m,u,y,h))
    pi=mp.pi
    A={}
    A["b1"]=-W(lambda m,u,y,h:m*u)/(4*pi)
    A["b2_y"]=-W(lambda m,u,y,h:m*y)/(2*pi)
    A["b2_x2"]=-W(lambda m,u,y,h:u*u/m)/(4*pi)
    A["cH"]=-W(lambda m,u,y,h:m*h)/(4*pi)
    A["a1"]=-W(lambda m,u,y,h:h*u/m)/(2*pi)
    A["a2_y"]=-W(lambda m,u,y,h:h*y/m)/(2*pi)
    A["a2_x2"]=W(lambda m,u,y,h:h*u*u/m**3)/(4*pi)
    A["lambda_neutral_mix"]=-c["k"]**2/(8*pi*(m[0]+m[1]))
    A["lambda"]=-W(lambda m,u,y,h:h*h/m)/(8*pi)+A["lambda_neutral_mix"]
    return dict(A,b2=A["b2_y"]+A["b2_x2"],a2=A["a2_y"]+A["a2_x2"])


def integrals(c,s):
    # Differentiate 1/2 Tr log of the Gaussian inverse propagator before
    # integrating. At rho=0, the propagator is diagonal. Only rho^2
    # can have a two-off-diagonal insertion at the derivatives tested.
    ms=m2at(c,s)
    ds=tuple(2*x+4*y*s for x,y in zip(c["x"],c["y"]))
    dss=tuple(4*y for y in c["y"])
    dr=tuple(2*h for h in c["h"])
    pi=mp.pi
    def integrand(p,key):
        P=p*p
        G=tuple(1/(P+q) for q in ms)
        # P*(G-1/P), convergent, no subtraction of nearly equal floating numbers.
        P_Gsub=tuple(-q*g for q,g in zip(ms,G))
        W=lambda f:sum(n*f(g,t,a,b,c) for n,g,t,a,b,c in zip(WEIGHTS,G,P_Gsub,ds,dss,dr))
        return {
           "b1": lambda : W(lambda g,t,a,b,c:t*a)/2,
           "b2_y": lambda : W(lambda g,t,a,b,c:t*b)/2,
           "b2_x2": lambda : -P*W(lambda g,t,a,b,c:(g*a)**2)/2,
           "cH": lambda : W(lambda g,t,a,b,c:t*c)/2,
           "a1": lambda : -P*W(lambda g,t,a,b,c:g*g*c*a),
           "a2_y": lambda : -P*W(lambda g,t,a,b,c:g*g*c*b)/2,
           "a2_x2": lambda : P*W(lambda g,t,a,b,c:g**3*c*a*a),
           "lambda_neutral_mix": lambda : -P*G[0]*G[1]*c["k"]**2/2,
           "lambda": lambda : -P*(W(lambda g,t,a,b,c:(g*c)**2)+2*G[0]*G[1]*c["k"]**2)/4,
        }[key]()/(2*pi*pi)
    keys=("b1","b2_y","b2_x2","cH","a1","a2_y","a2_x2","lambda_neutral_mix","lambda")
    R={k:mp.quad(lambda p:integrand(p,k),[0,mp.sqrt(ms[0]),mp.inf]) for k in keys}
    return dict(R,b2=R["b2_y"]+R["b2_x2"],a2=R["a2_y"]+R["a2_x2"])


def rel(a,b): return max(abs(a[k]-b[k])/max(abs(a[k]),abs(b[k]),mp.mpf(1)) for k in a)


def pack(x):
    if isinstance(x,dict): return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [pack(v) for v in x]
    if isinstance(x,mp.mpf): return mp.nstr(x,62)
    return x


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument("out",type=Path)
    p=cli.parse_args()
    mp.mp.dps=65
    records=[]
    for row in points():
        inp=matching(row); printed=matching(row,printed=True)
        s,d,rho=(mp.mpf(row[k]) for k in ("s","shift","rho"))
        a=analytical(inp,s); q=integrals(inp,s); a0=analytical(inp,0); b0=analytical(printed,0)
        v=determinant(inp,s,rho)
        e1=rel(a,q)
        e2=rel(a,analytical(shifted(inp,d),s-d))
        e3=abs(v-determinant(shifted(inp,d),s-d,rho))/max(1,abs(v))
        assert max(e1,e2,e3)<mp.mpf("1e-59"),(e1,e2,e3)
        records.append(dict(fixture=row,matched_3d=inp,analytical_background=a,quadrature_background=q,
              background_quadrature_max_scaled_error=e1,shift_derivatives_max_scaled_error=e2,
              shift_determinant_scaled_error=e3,origin_static_candidate=a0,origin_2026_printed_x=b0,
              origin_candidate_minus_printed={k:a0[k]-b0[k] for k in a0},
              origin_determinant=determinant(inp,0,0),nonzero_rho_determinant=v))
        print(row["label"],"quadrature/shift errors",*[mp.nstr(x,6) for x in (e1,e2,e3)])
        print("origin candidate minus printed:",{k:mp.nstr(a0[k]-b0[k],11) for k in ("b1","b2","a1","a2")})
    hash_file=lambda f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()
    result=pack(dict(scope="only static one-loop temporal-gauge threshold, not full EFT",
      runtime={"python":platform.python_version(),"mpmath":mp.__version__,"decimal_precision":65},
      files_sha256={k:hash_file(k) for k in ("threshold.py","fixtures.json","test_threshold.py")},records=records))
    p.out.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":main()
