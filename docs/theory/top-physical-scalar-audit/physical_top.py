"""NST A61 CP-even scalar pole terms, separate PV and shifted numerator."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import mpmath as mp

HERE=Path(__file__).resolve().parent


def checked(path,digest):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=digest:
        raise ValueError('Archived fixture hash mismatch: '+path.name)
    return raw


def verify_fixture(cfg,benchmark=None):
    prior=cfg['prior']
    old=json.loads(checked(HERE/prior['input_path'],prior['input_sha256']))
    for k in ('GF_GeV_minus2','MW_GeV','MZ_GeV','Mt_GeV','mu_bar_GeV'):
        if cfg[k]!=old[k]:
            raise ValueError('Changed common NST input: '+k)
    if benchmark:
        record=json.loads(checked(benchmark,cfg['journal_benchmark']['audited_input_sha256']))
        for a,b in [('Mh1_GeV','mh1_GeV'),('Mh2_GeV','mh2_GeV'),('sin_theta','sin_theta')]:
            if cfg[a]!=record[b]:
                raise ValueError('Changed journal benchmark input: '+a)
    src=HERE/prior['pv_source_path']
    checked(src,prior['pv_source_sha256'])
    spec=importlib.util.spec_from_file_location('xsm_pinned_neutral_pv',src)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.pv_moments


def direct(L,b,r=1):
    """Direct signed 1+x spinor projection. Physical scalar block uses r=1."""
    c=b-1-r
    roots=[]
    disc=c*c-4*r
    if disc>=0:
        mid=-c/(2*r)
        half=mp.sqrt(disc)/(2*r)
        roots=[mid-half,mid+half]
    interior=sorted({z for z in roots if 0<z<1})
    parts=[0,*interior,1] if interior else [0,mp.mpf('0.5'),1]
    def f(x):
        delta=r*(x-roots[0])*(x-roots[-1]) if roots else (1-x)*(1-r*x)+b*x
        if not delta:
            return mp.mpc(0)
        return (1+x)*(L-mp.log(abs(delta))+(mp.j*mp.pi if delta<0 else 0))
    return mp.quad(f,parts)


def mix(light,heavy,w):
    return (1-w)*light+w*heavy


def heavy_asym(L,b):
    return mp.mpf('1.5')*(L-mp.log(b))+mp.mpf('1.25')


def heavy_bound(b):
    return 2*((1+1/b)*mp.log1p(1/b)+mp.log(b)/b)


def result(cfg,pv):
    w,z,m,G,h1,h2,st=map(lambda k:mp.mpf(cfg[k]),
                         ('MW_GeV','MZ_GeV','Mt_GeV','GF_GeV_minus2','Mh1_GeV','Mh2_GeV','sin_theta'))
    s2=st*st
    y2=(8*w*w*G/mp.sqrt(2))*m*m/(2*w*w)
    norm=y2/(16*mp.pi**2)  # relative correction is twice the selected self-energy sum.
    rows=[]
    for mu in cfg['mu_bar_GeV']:
        L=2*mp.log(mp.mpf(mu)/m)
        slots=[]
        for mass in (h1,h2):
            b=(mass/m)**2
            b0,b1=pv(L,b,mp.mpf(1))
            reduced=b0-b1
            shifted=direct(L,b)
            slots.append(dict(mass_GeV=mass,B0_stripped=b0,B1_stripped=b1,PV=reduced,
                              direct=shifted,diff=shifted-reduced,unweighted_A30=norm*reduced))
        light,heavy=(s['unweighted_A30'] for s in slots)
        rel=mix(light,heavy,s2)
        rows.append(dict(mu_bar_GeV=mu,L=L,scalars=slots,
                         weighted_light_A30=(1-s2)*light,weighted_heavy_A30=s2*heavy,
                         A30_relative_1PI=rel,A30_additive_1PI=y2*rel,
                         A30_change_from_zero_mixing=s2*(heavy-light),A30_zero_mixing=light))
    L=2*mp.log(z/m)
    a0=mp.re(rows[0]['scalars'][0]['direct'])
    limits=[]
    for val in map(mp.mpf,cfg['heavy_mass_controls_GeV']):
        b=(val/m)**2
        d=mp.re(direct(L,b))
        a=heavy_asym(L,b)
        limits.append(dict(heavy_mass_GeV=val,direct_kernel=d,asymptotic_kernel=a,
                           asymptotic_minus_direct=a-d,bound=heavy_bound(b),
                           A30_difference_fixed_mixing=s2*norm*(d-a0),
                           A30_difference_scaling_mixing=s2*(h2/val)**2*norm*(d-a0)))
    return dict(inputs=cfg,dps=mp.mp.dps,yt02=y2,sin2_theta=s2,rows=rows,
                heavy_limits_mu_bar_MZ=limits)


def ser(x):
    if isinstance(x,dict):
        return {k:ser(v) for k,v in x.items()}
    if isinstance(x,list):
        return [ser(v) for v in x]
    if isinstance(x,mp.mpc):
        return [ser(mp.re(x)),ser(mp.im(x))]
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
    pv=verify_fixture(cfg,o.check_benchmark)
    o.output.write_text(json.dumps(ser(result(cfg,pv)),indent=2)+'\n')
