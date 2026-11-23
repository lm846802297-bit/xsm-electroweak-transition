"""Independent shifted numerator and PV routes for the neutral NST top block."""
import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp

HERE=Path(__file__).resolve().parent
HALF=lambda:mp.mpf('0.5')


def delta_data(b,r):
    """b=MZ^2/Mt^2; r=k^2/Mt^2. Finite kernels require b,r>0."""
    if b<=0 or r<=0:
        raise ValueError('Positive mass ratio and timelike invariant required')
    c=b-1-r
    disc=c*c-4*r
    bounds=[]
    if disc>=0:
        w=mp.sqrt(disc)
        bounds=[(-c-w)/(2*r),(-c+w)/(2*r)]
    return c,disc,bounds


def pv_moments(L,b,r):
    """Strip loop factor. Closed scalar logarithm and NST A42 vector reduction."""
    c,disc,roots=delta_data(b,r)
    if disc<0:
        h=mp.sqrt(-disc)/(2*r)
        def primitive(t):
            return t*mp.log(t*t+h*h)-2*t+2*h*mp.atan(t/h)
        v=c/(2*r)
        logint=mp.log(r)+primitive(v+1)-primitive(v)
        cut=mp.mpf(0)
    else:
        def xlog(t):
            return t*mp.log(abs(t)) if t else mp.mpf(0)
        lo,hi=roots
        logint=mp.log(r)+sum(xlog(1-a)+xlog(a)-1 for a in roots)
        cut=max(mp.mpf(0),min(1,hi)-max(0,lo))
    b0=L-logint+mp.j*mp.pi*cut
    b1=((b-1-r)*b0+(1+L)-b*(1+L-mp.log(b)))/(2*r)
    return b0,b1


def signed_trace_integrals(L,b,r,cL,cR,yt2):
    """Shift the spinor numerator directly; no scalar/vector PV moments."""
    c,disc,roots=delta_data(b,r)
    div=[mp.mpf(0)]+sorted({w for w in roots if 0<w<1})+[mp.mpf(1)]
    if len(div)==2:
        div.insert(1,HALF())
    def k(x):
        q=r*x*x+c*x+1 if disc<0 else r*(x-roots[0])*(x-roots[-1])
        if not q:  # endpoint of one subinterval
            return mp.mpc(0)
        return L-mp.log(abs(q))+(mp.j*mp.pi if q<0 else 0)
    def z(x):
        # D-vector numerator times the pole before subtraction.
        same=cL*cL+cR*cR
        flip=cL*cR
        return (-4*flip+same*x)*k(x)+2*flip-same*x
    def g(x):
        # Two pseudoscalar vertices give the projected numerator x-1.
        return yt2/2*(x-1)*k(x)
    return mp.quad(z,div),mp.quad(g,div)


def fixture(cfg,path):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=cfg['fixture']['sha256']:
        raise ValueError('Prior charged block input SHA mismatch')
    prior=json.loads(raw)
    for k in ('GF_GeV_minus2','MW_GeV','MZ_GeV','Mt_GeV','mu_bar_GeV','gauge_xi'):
        if cfg[k]!=prior[k]:
            raise ValueError('Changed prior fixture value: '+k)


def output(cfg):
    if cfg['gauge_xi']!='1':
        raise ValueError('Only xi=1 supported')
    mw,mz,mt,gf=map(lambda k:mp.mpf(cfg[k]),('MW_GeV','MZ_GeV','Mt_GeV','GF_GeV_minus2'))
    g2=8*mw**2*gf/mp.sqrt(2)
    sw2=1-(mw/mz)**2
    gp2=g2*sw2/(1-sw2)
    gz=mp.sqrt(g2+gp2)
    yt2=g2*mt**2/(2*mw**2)
    cL,cR=gz*(HALF()-2*sw2/3),gz*(-2*sw2/3)
    # Use printed A61 coefficients, independently of the chiral-vertex route.
    a=gp2*(3*g2-gp2)/(9*(g2+gp2))
    v=(9*g2*g2-6*g2*gp2+17*gp2*gp2)/(72*(g2+gp2))
    loop=16*mp.pi**2
    b=(mz/mt)**2
    eps=-2*a-v
    rows=[]
    for mu in cfg['mu_bar_GeV']:
        L=2*mp.log(mp.mpf(mu)/mt)
        b0,b1=pv_moments(L,b,mp.mpf(1))
        zz=(4*a*b0-2*v*b1+eps)/loop
        gg=-yt2*(b0+b1)/(2*loop)
        zd,gd=signed_trace_integrals(L,b,mp.mpf(1),cL,cR,yt2)
        total=zz+gg
        rows.append({'mu_bar_GeV':mu,'L':L,'B0_without_loop_factor':b0,'B1_without_loop_factor':b1,
                     'Z_PV':zz,'G0_PV':gg,'Z_direct':zd/loop,'G0_direct':gd/loop,
                     'Z_direct_minus_PV':zd/loop-zz,'G0_direct_minus_PV':gd/loop-gg,
                     'Z_evanescent':eps/loop,'sum':total,'A30_relative_neutral_only':2*mp.re(total),
                     'A30_additive_neutral_only':yt2*2*mp.re(total),
                     'A30_relative_early_D4':2*mp.re(total-eps/loop)})
    return dict(inputs=cfg,dps=mp.mp.dps,g02=g2,gprime02=gp2,yt02=yt2,cL=cL,cR=cR,
                A61_mass_coefficient=a,A61_vector_coefficient=v,MZ2_over_Mt2=b,rows=rows)


def serialize(x):
    if isinstance(x,dict):
        return {k:serialize(v) for k,v in x.items()}
    if isinstance(x,list):
        return [serialize(v) for v in x]
    if isinstance(x,mp.mpc):
        return [serialize(mp.re(x)),serialize(mp.im(x))]
    if isinstance(x,mp.mpf):
        return mp.nstr(x,50)
    return x


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('output',type=Path)
    p.add_argument('--inputs',type=Path,default=HERE/'inputs.json')
    p.add_argument('--prior-inputs',type=Path)
    o=p.parse_args()
    mp.mp.dps=65
    cfg=json.loads(o.inputs.read_bytes())
    fixture(cfg,o.prior_inputs or HERE/cfg['fixture']['path'])
    o.output.write_text(json.dumps(serialize(output(cfg)),indent=2)+'\n')
