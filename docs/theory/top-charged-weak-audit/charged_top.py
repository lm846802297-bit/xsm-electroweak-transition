"""NST A61 on-shell charged W and Goldstone corrections with massless bottom."""
import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp

HERE=Path(__file__).resolve().parent


def moments(L,b):
    """Finite B0 and B1, without 1/(16 pi^2), from closed B0 and A42."""
    z=1-b
    edge=z*mp.log(abs(z)) if z else mp.mpf(0)
    b0=L+2-edge-b*mp.log(b)+mp.j*mp.pi*max(z,0)
    b1=(-z*b0-b*(1+L-mp.log(b)))/2
    return b0,b1


def direct_moment(L,b):
    """Integrate the shifted linear fermion numerator with the causal branch."""
    z=1-b
    def f(x):
        d=x*(x-z)
        if not d:
            return mp.mpc(0)
        return x*(L-mp.log(abs(d))+(mp.j*mp.pi if d<0 else 0))
    ends=[0,z,1] if 0<z<1 else [0,mp.mpf('.5'),1]
    return mp.quad(f,ends)


def project(L,b,g2,yt2,v2):
    b0,b1=moments(L,b)
    i=direct_moment(L,b)
    common=v2/(16*mp.pi**2)
    pvW=-common*g2/4*(2*b1+1)
    pvG=-common*yt2/2*b1
    dirW=common*g2/4*(2*i-1)
    dirG=common*yt2/2*i
    return {'B0_finite_without_loop_factor':b0,'B1_finite_without_loop_factor':b1,
            'B1_pole_without_loop_factor':-mp.mpf(1)/2,'B1_minus_shifted_projection':b1+i,
            'W_PV':pvW,'Goldstone_PV':pvG,'W_trace':dirW,'Goldstone_trace':dirG,
            'W_direct_minus_PV':dirW-pvW,'Goldstone_direct_minus_PV':dirG-pvG,
            'W_evanescent_finite':-common*g2/4,'sum':pvW+pvG,
            'A30_charged_relative':2*mp.re(pvW+pvG),
            'A30_charged_relative_early_D4':2*mp.re(pvW+pvG+common*g2/4)}


def result(cfg):
    W,Z,M,G,V=map(lambda k:mp.mpf(cfg[k]),('MW_GeV','MZ_GeV','Mt_GeV','GF_GeV_minus2','Vtb_abs'))
    if cfg['mb_GeV']!='0' or cfg['gauge_xi']!='1':
        raise ValueError('Only the massless-bottom xi=1 block is implemented.')
    g2=8*W**2*G/mp.sqrt(2)
    yt2=g2*M**2/(2*W**2)
    b=(W/M)**2
    w_width=V**2*g2*(M**2-W**2)*(2+M**2/W**2)/4 * (M**2-W**2)/(16*mp.pi*M**3)
    rows=[]
    for mu in cfg['mu_bar_GeV']:
        L=2*mp.log(mp.mpf(mu)/M)
        p=project(L,b,g2,yt2,V**2)
        rows.append({'mu_bar_GeV':mu,'L':L,**p,
                     'A30_charged_additive':yt2*p['A30_charged_relative'],
                     'bW_width_from_optical_GeV':2*M*mp.im(p['sum']),
                     'bW_width_from_spin_sum_GeV':w_width})
    return {'dps':mp.mp.dps,'inputs':cfg,'g02':g2,'yt02':yt2,'mw2_over_mt2':b,'rows':rows}


def serialized(o):
    if isinstance(o,dict):
        return {k:serialized(v) for k,v in o.items()}
    if isinstance(o,list):
        return [serialized(i) for i in o]
    if isinstance(o,mp.mpc):
        return [serialized(mp.re(o)),serialized(mp.im(o))]
    if isinstance(o,mp.mpf):
        return mp.nstr(o,50)
    return o


def verify_prior(cfg,path):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=cfg['comparison_fixture']['sha256']:
        raise ValueError('Pinned historical fixture SHA-256 mismatch.')
    data=json.loads(raw)
    for key in ('GF_GeV_minus2','MW_GeV','MZ_GeV','Mt_GeV','mu_bar_GeV'):
        if cfg[key]!=data[key]:
            raise ValueError('Changed benchmark input: '+key)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('output',type=Path)
    p.add_argument('--inputs',type=Path,default=HERE/'inputs.json')
    p.add_argument('--check-prior',type=Path,help='verify the immutable comparison input when available')
    opt=p.parse_args()
    mp.mp.dps=65
    cfg=json.loads(opt.inputs.read_bytes())
    if opt.check_prior:
        verify_prior(cfg,opt.check_prior)
    opt.output.write_text(json.dumps(serialized(result(cfg)),indent=2)+'\n')
