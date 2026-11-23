"""NST massless W bubbles and composition of archived closed-fermion W/Z inputs."""
import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
WTB = HERE.parent/'ew-w-top-bottom-audit'/'native-results.json'
Z = HERE.parent/'ew-z-fermion-audit'/'native-results.json'


def signed(value):
    if isinstance(value,list):
        return mp.mpc(*value)
    return mp.mpf(value)


def read_pinned(path, sha):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha:
        raise ValueError(f'Dependency hash does not match: {path.name}')
    return json.loads(raw)


def direct_w(s, mu2, g2, multiplicity):
    if s == 0:
        return mp.mpf(0)
    def integrand(x):
        u=x*(1-x)
        if u==0:
            return mp.mpf(0)
        return 2*u*s*(mp.log(mu2/(u*s))+mp.j*mp.pi)
    return g2*multiplicity/(16*mp.pi**2)*mp.quad(integrand,[0,mp.mpf("0.5"),1])


def pv_w(s, mu2, g2, multiplicity, early_d4=False):
    if s == 0:
        return mp.mpf(0)
    b0=mp.log(mu2/s)+2+mp.j*mp.pi
    return g2*multiplicity/(16*mp.pi**2)*s*(b0/3-(0 if early_d4 else mp.mpf(1)/9))


def projections(mw2, mz2, w0, wp, zp):
    """NST A28/A29 self-energy terms only; arguments may be complex at poles."""
    sw, sz = mp.re(wp), mp.re(zp)
    dg=(sw-w0)/mw2
    dgp=dg+(sz-mz2*sw/mw2)/(mz2-mw2)
    return {'A28_dg2_relative':dg, 'A29_dgprime2_relative':dgp,
            'dv2_relative_from_W_input':sw/mw2-dg}


def run(cfg, wtb_path, z_path):
    deps=cfg['archive_dependencies']
    wtb=read_pinned(wtb_path,deps['wtb_native_sha256'])
    zz=read_pinned(z_path,deps['z_native_sha256'])
    mw2, mz2=(mp.mpf(cfg[k])**2 for k in ('MW_GeV','MZ_GeV'))
    g2=8*mw2*mp.mpf(cfg['GF_GeV_minus2'])/mp.sqrt(2)
    gp2=g2*(mz2-mw2)/mw2
    mult=cfg['W_massless_quark_families']*cfg['quark_colors']+cfg['W_massless_lepton_families']
    rows=[]
    for mu in cfg['mu_bar_GeV']:
        w=next(r for r in wtb['rows'] if r['mu_bar']==mu)
        z=next(r for r in zz['rows'] if r['scale_GeV']==mu and r['point']=='pole')
        if abs(signed(w['g0sq_GF_tree'])-g2)>mp.mpf('1e-51') or abs(signed(z['s_GeV2'])-mz2)>mp.mpf('1e-46'):
            raise ValueError('Mismatched archived EW inputs')
        mu2=mp.mpf(mu)**2
        w_light=direct_w(mw2,mu2,g2,mult)
        wp_tb=signed(w['PiW_tb_MSfinite_at_pole_GeV2'])
        w0=signed(w['PiW_tb_MSfinite_at_zero_GeV2'])
        zp=signed(z['all_GeV2'])
        mix=projections(mw2,mz2,w0,wp_tb+w_light,zp)
        p_light=g2*mult*mw2/(48*mp.pi**2)
        poles=projections(mw2,mz2,signed(w['PiW_tb_UV_pole_coeff_at_zero_GeV2']),
                          signed(w['PiW_tb_UV_pole_coeff_at_pole_GeV2'])+p_light,
                          sum(signed(i['pole_coeff_GeV2']) for i in z['groups'].values()))
        rows.append({'mu_bar_GeV':mu,'W_light_zero_GeV2':direct_w(mp.mpf(0),mu2,g2,mult),
                     'W_light_direct_pole_GeV2':w_light,'W_light_PV_pole_GeV2':pv_w(mw2,mu2,g2,mult),
                     'W_light_PV_early_D4_GeV2':pv_w(mw2,mu2,g2,mult,True),
                     'W_light_direct_minus_PV_GeV2':w_light-pv_w(mw2,mu2,g2,mult),
                     'W_light_width_from_Impi_GeV':mp.im(w_light)/mp.sqrt(mw2),
                     'W_light_tree_width_GeV':g2*mult*mp.sqrt(mw2)/(48*mp.pi),
                     'W_full_closed_fermion_zero_GeV2':w0,'W_full_closed_fermion_pole_GeV2':wp_tb+w_light,
                     'Z_full_closed_fermion_pole_GeV2':zp,**mix,
                     'A28_dg2_additive':g2*mix['A28_dg2_relative'],
                     'A29_dgprime2_additive':gp2*mix['A29_dgprime2_relative'],
                     'pole_coefficients_in_self_energy_projections':poles})
    return {'working_decimal_digits':mp.mp.dps,'inputs':cfg,'g02':g2,'gprime02':gp2,'W_massless_multiplicity':mult,'rows':rows}


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
    p.add_argument('--wtb',type=Path,default=WTB)
    p.add_argument('--z',type=Path,default=Z)
    args=p.parse_args()
    mp.mp.dps=65
    args.output.write_text(json.dumps(serialize(run(json.loads(args.inputs.read_bytes()),args.wtb,args.z)),indent=2)+'\n')
