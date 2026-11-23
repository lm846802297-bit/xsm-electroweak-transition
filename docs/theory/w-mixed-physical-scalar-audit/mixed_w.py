"""A58 mixed physical Higgs/charged vector and Goldstone bubbles, xi=1."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import mpmath as mp

HERE=Path(__file__).resolve().parent


def read_exact(path,sha):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha:
        raise ValueError('Pinned dependency mismatch: '+path.name)
    return raw


def fixture(cfg):
    p=cfg['prior']
    old=json.loads(read_exact(HERE/p['input_path'],p['input_sha256']))
    for k in ('GF_GeV_minus2','MW_GeV','MZ_GeV','Mt_GeV','Mh1_GeV','Mh2_GeV','sin_theta','gauge_xi','mu_bar_GeV'):
        if old[k]!=cfg[k]:
            raise ValueError('Prior input mismatch: '+k)
    source=HERE/p['pv_path']
    read_exact(source,p['pv_sha256'])
    spec=importlib.util.spec_from_file_location('pinned_analytic_pv',source)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.pv_moments


def shifted(a,b,s,mu2):
    """Direct D-dimensional shift and IBP; stripped 1/(16 pi^2)."""
    c=b-a-s
    disc=c*c-4*s*a
    roots=[(-c-mp.sqrt(disc))/(2*s),(-c+mp.sqrt(disc))/(2*s)] if s and disc>=0 else []
    real=sorted({r for r in roots if 0<r<1})
    edges=[0,*real,1] if real else [0,mp.mpf('.5'),1]
    def delta(x):
        return s*(x-roots[0])*(x-roots[-1]) if roots else a*(1-x)+b*x-s*x*(1-x)
    def log(q):
        return mp.log(mu2/abs(q))+(mp.j*mp.pi if q<0 else 0) if q else mp.mpf(0)
    b0=mp.quad(lambda x:log(delta(x)),edges)
    b00=mp.quad(lambda x:delta(x)*(1+log(delta(x)))/2,edges)
    return b0,b00


def zero_moments(a,b,mu2):
    if a==b:
        L=mp.log(mu2/a)
        return L,-L/2
    w=b-a
    J0=lambda q:q*(mp.log(q)-1)
    J1=lambda q:q*q*(mp.log(q)/2-mp.mpf('.25'))
    K0=J0(b)-J0(a)
    K1=J1(b)-J1(a)
    return mp.log(mu2)-K0/w,-mp.log(mu2)/2+(K1-a*K0)/(w*w)


def reduced(a,b,s,mu2,pv):
    """A41/A42 tensor reduction and a genuinely separate A44 zero check."""
    if s:
        b0,b1=pv(mp.log(mu2/a),b/a,s/a)
    else:
        b0,b1=zero_moments(a,b,mu2)
    A_b=b*(1+mp.log(mu2/b))
    n=2*a*b0+A_b-(b-a-s)*b1
    pole=3*(a+b)/2-s/2
    bareD4=n/6
    valid=bareD4+pole/9
    A44=(A_b+a*b0)/4+(a+b)/8 if not s else None
    return {'B0':b0,'B1':b1,'B00':valid,'B00_early_D4':bareD4,
            'B00_zero_A44':A44,'B00_evanescent':pole/9}


def output(cfg,pv):
    if cfg['gauge_xi']!='1':
        raise ValueError('xi=1 only')
    G,W,m1,m2,st=map(lambda k:mp.mpf(cfg[k]),('GF_GeV_minus2','MW_GeV','Mh1_GeV','Mh2_GeV','sin_theta'))
    b=W*W
    g2=8*b*G/mp.sqrt(2)
    loop=16*mp.pi**2
    weights=(1-st*st,st*st)
    rows=[]
    for scale in cfg['mu_bar_GeV']:
        mu2=mp.mpf(scale)**2
        points=[]
        for s in (mp.mpf(0),b):
            states=[]
            for m,w in zip((m1,m2),weights):
                a=m*m
                c=g2*w/loop
                ib0,ib00=shifted(a,b,s,mu2)
                red=reduced(a,b,s,mu2,pv)
                p=c*(b*red['B0']-red['B00'])
                p0=c*(b*red['B0']-red['B00_early_D4'])
                d=c*(b*ib0-ib00)
                states.append(dict(mass_GeV=m,weight=w,direct_B0=ib0,direct_B00_GeV2=ib00,
                                   reduced_stripped=red,Pi_direct_GeV2=d,Pi_PV_GeV2=p,
                                   Pi_PV_early_D4_GeV2=p0,Pi_difference_GeV2=d-p,
                                   Pi_uv_pole_GeV2=c*(b-(a+b)/4+s/12)))
            points.append(dict(s_GeV2=s,states=states,Pi_GeV2=sum(v['Pi_PV_GeV2'] for v in states),
                               Pi_early_D4_GeV2=sum(v['Pi_PV_early_D4_GeV2'] for v in states)))
        z,p=points
        dg=mp.re(p['Pi_GeV2']-z['Pi_GeV2'])/b
        early=mp.re(p['Pi_early_D4_GeV2']-z['Pi_early_D4_GeV2'])/b
        rows.append(dict(mu_bar_GeV=scale,points=points,A28_W_bubble_relative=dg,
                         A28_wrong_early_D4=early,A30_W_input_relative_only=dg-mp.re(p['Pi_GeV2'])/b,
                         A28_log_scale_slope=g2/(12*loop)))
    return dict(inputs=cfg,dps=mp.mp.dps,g02=g2,rows=rows)


def serial(x):
    if isinstance(x,dict):
        return {k:serial(v) for k,v in x.items()}
    if isinstance(x,list):
        return [serial(v) for v in x]
    if isinstance(x,mp.mpc):
        return [serial(mp.re(x)),serial(mp.im(x))]
    if isinstance(x,mp.mpf):
        return mp.nstr(x,50)
    return x


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('output',type=Path)
    p.add_argument('--inputs',type=Path,default=HERE/'inputs.json')
    o=p.parse_args()
    mp.mp.dps=65
    c=json.loads(o.inputs.read_bytes())
    o.output.write_text(json.dumps(serial(output(c,fixture(c))),indent=2)+'\n')
