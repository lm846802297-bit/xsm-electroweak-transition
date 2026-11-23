"""NST common vacuum tadpole, diagram and independent effective-potential methods."""
import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp

HERE=Path(__file__).resolve().parent
KINDS=('top','vector','goldstone','physical')


def fixture(cfg,bench=None):
    p=cfg['prior']
    path=HERE/p['path']
    data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=p['sha256']:
        raise ValueError('prior input hash mismatch')
    old=json.loads(data)
    for k in ('GF_GeV_minus2','MW_GeV','MZ_GeV','Mt_GeV','Mh1_GeV','Mh2_GeV','sin_theta','mu_bar_GeV'):
        if old[k]!=cfg[k]:
            raise ValueError(f'prior pole input differs: {k}')
    if bench:
        b=bench.read_bytes()
        if hashlib.sha256(b).hexdigest()!=cfg['sources']['journal_benchmark']['audited_input_sha256']:
            raise ValueError('journal input hash mismatch')
        b=json.loads(b)
        for k in ('a2','b3_GeV','b4','sin_theta'):
            if b[k]!=cfg[k]:
                raise ValueError(f'journal benchmark differs: {k}')


def ew_model(cfg):
    gf,w,z,t,m1,m2,sj,a2,b3,b4=map(lambda k:mp.mpf(cfg[k]),(
        'GF_GeV_minus2','MW_GeV','MZ_GeV','Mt_GeV','Mh1_GeV','Mh2_GeV','sin_theta','a2','b3_GeV','b4'))
    if cfg['gauge_xi']!='1':
        raise ValueError('Only xi=1 is in scope')
    v=(mp.sqrt(2)*gf)**(-mp.mpf('.5'))
    g2=4*w*w/(v*v)
    gz2=4*z*z/(v*v)
    s=-sj
    c=mp.sqrt(1-s*s)
    m1,m2=m1*m1,m2*m2
    a1=2*c*s*(m2-m1)/v
    lam=(c*c*m1+s*s*m2)/(2*v*v)
    b2=s*s*m1+c*c*m2-a2*v*v/2
    return dict(v=v,g2=g2,gz2=gz2,mw=w*w,mz=z*z,mt=t*t,mh=(m1,m2),c=c,s=s,a1=a1,a2=a2,
                b1=-a1*v*v/4,b2=b2,b3=b3,b4=b4,lam=lam)


def nst_diagrams(m,mu):
    v,g2,gz2,w,z,t,c,s,a,b=(m[k] for k in ('v','g2','gz2','mw','mz','mt','c','s','a2','b3'))
    h1,h2=m['mh']
    co=c*c-s*s
    c111=-3*h1*c**3/v-3*v*a*c*s*s+2*b*s**3
    c222=-3*h2*s**3/v-3*v*a*s*c*c-2*b*c**3
    c112=s*(v*v*a-2*h1-h2+(3*v*v*a-2*h1-h2)*co)/(2*v)-2*b*c*s*s
    c122=c*(v*v*a-h1-2*h2-(3*v*v*a-h1-2*h2)*co)/(2*v)+2*b*s*c*c
    den=16*mp.pi**2
    A=lambda q:q*(1+2*mp.log(mu/mp.sqrt(q)))/den
    DA=lambda q:3*A(q)-2*q/den  # (D-1)A through MS finite part
    return {k:xs for k,xs in zip(KINDS,[
        [-12*r*t*A(t)/v for r in (c,s)],
        [v*r*(2*g2*DA(w)+gz2*DA(z))/4 for r in (c,s)],
        [h*r*(2*A(w)+A(z))/(2*v) for h,r in zip((h1,h2),(c,s))],
        [-(c111*A(h1)+c122*A(h2))/2,-(c112*A(h1)+c222*A(h2))/2]])}


def Hessian(m,h,S):
    v,a1,a2,b2,b3,b4,l=(m[k] for k in ('v','a1','a2','b2','b3','b4','lam'))
    p=v+h
    hS=p*(a1/2+a2*S)
    return mp.matrix([[-l*v*v+3*l*p*p+(a1*S+a2*S*S)/2,hS],
                      [hS,b2+2*b3*S+3*b4*S*S+a2*p*p/2]])


def CW(m,mu,h,S):
    v,a1,a2,l=(m[k] for k in ('v','a1','a2','lam'))
    p=v+h
    w,z,t=(m[k]*(p/v)**2 for k in ('mw','mz','mt'))
    G=-l*v*v+l*p*p+(a1*S+a2*S*S)/2
    vals=mp.eigsy(Hessian(m,h,S),eigvals_only=True)
    f=lambda q,d:q*q*(mp.log(q/(mu*mu))-d)/(64*mp.pi**2)
    sc=mp.mpf('1.5')
    vec=mp.mpf(5)/6
    return {'top':-12*f(t,sc),'vector':6*f(w,vec)+3*f(z,vec),
            'goldstone':2*(f(w+G,sc)-f(w,sc))+f(z+G,sc)-f(z,sc),
            'physical':sum(f(q,sc) for q in vals)}


def potential_gradient(m,mu,step):
    # Differentiate the whole field-dependent determinant in the interaction basis.
    data=[]
    for j in range(2):
        a=[CW(m,mu,d*step if j==0 else 0,d*step if j==1 else 0) for d in (-2,-1,1,2)]
        data.append({k:(a[0][k]-8*a[1][k]+8*a[2][k]-a[3][k])/(12*step) for k in KINDS})
    return data


def response(m,T):
    v=m['v']
    # Use separate tree mass gradients for W/Z/top.
    w=-(m['g2']*v/2)*T
    z=-(m['gz2']*v/2)*T
    top=-T/v
    dg=(w-w)/m['mw']
    dgp=dg+(z-m['mz']*w/m['mw'])/(m['mz']-m['mw'])
    return dict(T_h_GeV=T,PiW_tad_GeV2=w,PiZ_tad_GeV2=z,Sigma_top_tad=top,
                W_input_A30=-w/m['mw'],top_A30=2*top,A30_tad_sum=dg-w/m['mw']+2*top,
                A28_tad=dg,A29_tad=dgp,pole_matching_relative_dv2=w/m['mw']-dg,
                local_minimum_relative_dv=T/v)


def result(cfg):
    m=ew_model(cfg)
    c,s,h=m['c'],m['s'],m['mh']
    H=Hessian(m,0,0)
    rows=[]
    for scale in cfg['mu_bar_GeV']:
        mu=mp.mpf(scale)
        di=nst_diagrams(m,mu)
        trs=[]
        for f in map(mp.mpf,cfg['fd_steps_GeV']):
            hx,sx=potential_gradient(m,mu,f)
            pot={k:[-c*hx[k]+s*sx[k],-s*hx[k]-c*sx[k]] for k in KINDS}
            err=max(abs(di[k][i]-pot[k][i])/max(1,abs(di[k][i])) for k in KINDS for i in range(2))
            dh=mp.lu_solve(H,-mp.matrix([sum(hx.values()),sum(sx.values())]))
            trs.append(dict(step_GeV=f,tadpoles_from_potential=pot,maximum_scaled_error=err,
                            local_linear_response_doublet_GeV=dh[0],local_linear_response_singlet_GeV=dh[1]))
        ti=[sum(di[k][i] for k in KINDS) for i in range(2)]
        T=c*ti[0]/h[0]+s*ti[1]/h[1]
        ans=response(m,T)
        rows.append(dict(mu_bar_GeV=scale,diagram_tadpoles_GeV3=di,diagram_totals_GeV3=ti,
                         response=ans,independent_potential=trs))
    return dict(inputs=cfg,dps=mp.mp.dps,tree=m,rows=rows)


def serial(x):
    if isinstance(x,dict):
        return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):
        return [serial(y) for y in x]
    if isinstance(x,mp.mpf):
        return mp.nstr(x,50)
    return x


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('output',type=Path)
    p.add_argument('--inputs',type=Path,default=HERE/'inputs.json')
    p.add_argument('--check-benchmark',type=Path)
    o=p.parse_args()
    mp.mp.dps=65
    cfg=json.loads(o.inputs.read_bytes())
    fixture(cfg,o.check_benchmark)
    o.output.write_text(json.dumps(serial(result(cfg)),indent=2)+'\n')
