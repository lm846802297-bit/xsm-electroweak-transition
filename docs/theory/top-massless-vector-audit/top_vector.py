"""NST A61 on-shell normalization for the gluon and photon exchange diagrams."""
import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp

HERE=Path(__file__).resolve().parent


def pv(L):
    """Strip 1/(16 pi^2), then multiply D and D-2 before MS-bar subtraction."""
    b0=L+2
    b1=-L/2-mp.mpf(3)/2
    p0,p1=mp.mpf(1),-mp.mpf(1)/2
    early=4*b0+2*b1
    evanescent=-2*(p0+p1)
    return {'B0_finite':b0,'B1_finite':b1,'B0_pole':p0,'B1_pole':p1,
            'combination_pole':4*p0+2*p1,'combination_finite':early+evanescent,
            'combination_early_D4':early,'evanescent':evanescent}


def trace(L, below_ratio=None):
    """Integrate projected D-dimensional vector numerator, no B0/B1 identity.

    At the pole external spinors put /p=M. The optional below-threshold
    parameter varies only the denominator, for an endpoint limiting check.
    It is not an off-shell mass-shift calculation.
    """
    r=mp.mpf(1) if below_ratio is None else below_ratio
    if not 0 <= r <= 1:
        raise ValueError("Only the real approach from below the cut is used.")
    def integrand(x):
        # Exact on-shell delta/m^2=(1-x)^2. Quadrature has no x=1 node.
        u=1-x
        if u==0:
            return mp.mpf(0)
        delta=(u*u if r==1 else u*(1-r*x))
        return (4-2*x)*(L-mp.log(delta))-2*u
    return mp.quad(integrand,[0,mp.mpf('0.5'),1])


def output(cfg):
    mw,mz,mt,gf=map(lambda k:mp.mpf(cfg[k]),('MW_GeV','MZ_GeV','Mt_GeV','GF_GeV_minus2'))
    gs2,cf,qt=map(lambda k:mp.mpf(cfg[k]),('gs_squared_fixed_NST_A2','Cf','Qt'))
    g2=8*mw**2*gf/mp.sqrt(2)
    sw2=1-(mw/mz)**2
    e2=g2*sw2
    yt2=g2*mt**2/(2*mw**2)
    rows=[]
    for mu in cfg['mu_bar_GeV']:
        L=2*mp.log(mp.mpf(mu)/mt)
        a=pv(L)
        exact=4+3*L
        num=trace(L)
        qcd=-cf*gs2*a['combination_finite']/(16*mp.pi**2)
        qed=-qt**2*e2*a['combination_finite']/(16*mp.pi**2)
        coeff=-2*(cf*gs2+qt**2*e2)/(16*mp.pi**2)
        rows.append({'mu_bar_GeV':mu,'L':L,**a,'trace_finite':num,
                     'trace_minus_A61':num-a['combination_finite'],'analytic_finite':exact,
                     'Sigma_sum_gluon':qcd,'Sigma_sum_photon':qed,
                     'A30_yt2_relative_gluon_only':2*qcd,'A30_yt2_relative_photon_only':2*qed,
                     'A30_yt2_relative_gluon_plus_photon_only':2*(qcd+qed),
                     'A30_yt2_additive_gluon_plus_photon_only':yt2*2*(qcd+qed),
                     'A30_yt2_relative_early_D4':coeff*a['combination_early_D4'],
                     'A30_yt2_relative_UV_pole_coeff':3*coeff})
    return {'inputs':cfg,'digits':mp.mp.dps,'g02':g2,'sW02':sw2,'e02':e2,'yt02':yt2,
            'alpha_s_fixed_NST':gs2/(4*mp.pi),'one_loop_zero_GeV':mt*mp.exp(-mp.mpf(2)/3),
            'rows':rows}


def fixture_check(cfg, wpath):
    raw=wpath.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=cfg['src_WZ_inputs']['sha256']:
        raise ValueError('Exact archived W/Z inputs hash mismatch.')
    old=json.loads(raw)
    for k in ('GF_GeV_minus2','MW_GeV','MZ_GeV','Mt_GeV','mu_bar_GeV'):
        if old[k]!=cfg[k]:
            raise ValueError('Archived W/Z fixture mismatch: '+k)


def serialize(x):
    if isinstance(x,dict):
        return {k:serialize(v) for k,v in x.items()}
    if isinstance(x,list):
        return [serialize(v) for v in x]
    if isinstance(x,mp.mpf):
        return mp.nstr(x,50)
    return x


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('output',type=Path)
    p.add_argument('--inputs',type=Path,default=HERE/'inputs.json')
    p.add_argument('--wz-inputs',type=Path)
    args=p.parse_args()
    mp.mp.dps=65
    cfg=json.loads(args.inputs.read_bytes())
    fixture_check(cfg,args.wz_inputs or HERE/cfg['src_WZ_inputs']['path'])
    args.output.write_text(json.dumps(serialize(output(cfg)),indent=2)+'\n')
