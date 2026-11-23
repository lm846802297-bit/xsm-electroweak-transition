"""NST Z closed-fermion block; all functions expect m^2 and s=k^2 in GeV^2."""
import argparse
import json
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent


def u(x):
    return x * (1 - x)


def direct(m2, s, mu2, v, a):
    """Trace + shifted D-dimensional parameter reduction, no PV identities."""
    if m2 < 0 or s < 0 or (m2 > 0 and s >= 4*m2):
        raise ValueError("Use nonnegative s and stay below the massive cut.")
    if m2 == s == 0:
        return mp.mpf(0)  # scaleless in dimensional regularization
    def integrand(x):
        delta = m2 - u(x)*s
        if not delta:
            return mp.mpf(0)  # endpoint of massless weighted log
        log = mp.log(mu2 / abs(delta))
        if delta < 0:
            log += 1j*mp.pi  # log(mu^2/(Delta-i0))
        return 2 * (v*v*u(x)*s - a*a*delta)*log
    return mp.quad(integrand, [0, mp.mpf('0.5'), 1])


def b0_closed(m2, s, mu2):
    """Analytic finite B0, here defined at zero and below the massive cut."""
    if not m2:
        if not s:
            return mp.mpf(0)
        return mp.log(mu2/s) + 2 + 1j*mp.pi
    if s == 0:
        return mp.log(mu2/m2)
    if not 0 < s < 4*m2:
        raise ValueError('Massive comparison intentionally limited to 0<s<4m^2.')
    r = mp.sqrt(4*m2/s-1)
    return mp.log(mu2/m2) + 2 - 2*r*mp.atan(1/r)


def pv(m2, s, mu2, q, t3, sw2, early_d4=False):
    """A59/A60 with A41/A42/A43; retain the numerator pole of 2B00."""
    if m2 == s == 0:
        return mp.mpf(0)
    af = 4*q*sw2*(2*q*sw2-2*t3)
    A = m2*(1+mp.log(mu2/m2)) if m2 else mp.mpf(0)
    B = b0_closed(m2,s,mu2)
    n = (2*m2-s/2)*B + A
    n_pole = 3*m2-s/2
    b00_twice = n/3 + (0 if early_d4 else 2*n_pole/9)
    return ((1+af)*(b00_twice-A-(m2-s/2)*B) + af*m2*B)/2


def pole_reduced(m2, s, v, a):
    return (v*v+a*a)*s/3 - 2*a*a*m2


def group_values(cfg, s, mu2, mw2, mz2, g02):
    sw2 = 1-mw2/mz2
    gz2 = g02/(1-sw2)
    result = {}
    for f in cfg["fermions"]:
        key=f["group"]
        m2=mp.mpf(cfg[f["mass_key"]])**2 if "mass_key" in f else mp.mpf(0)
        q, t3 = [mp.mpf(f[k]) for k in ("charge","T3")]
        v,a = t3-2*q*sw2,t3
        scalar = f["Nc"]*f["copies"]*gz2/(16*mp.pi**2)
        traced = scalar*direct(m2,s,mu2,v,a)
        reduced = scalar*pv(m2,s,mu2,q,t3,sw2)
        early = scalar*pv(m2,s,mu2,q,t3,sw2,early_d4=True)
        treewidth = f["Nc"]*f["copies"]*gz2*(v*v+a*a)*mp.sqrt(mz2)/(48*mp.pi) if not m2 else mp.mpf(0)
        result[key] = {"direct_GeV2": traced, "pv_GeV2": reduced, "difference_GeV2": traced-reduced,
                       "early_d4_GeV2": early, "pole_coeff_GeV2": scalar*pole_reduced(m2,s,v,a),
                       "massless_tree_partial_width_GeV": treewidth}
    return result


def results(cfg):
    mw2,mz2 = [mp.mpf(cfg[k])**2 for k in ('MW_GeV','MZ_GeV')]
    g02 = 8*mw2*mp.mpf(cfg['GF_GeV_minus2'])/mp.sqrt(2)
    gp02 = g02*(mz2-mw2)/mw2
    rows=[]
    for scale in cfg['mu_bar_GeV']:
        mu2=mp.mpf(scale)**2
        for label,s in [('zero',mp.mpf(0)),('pole',mz2)]:
            groups=group_values(cfg,s,mu2,mw2,mz2,g02)
            top=groups['top']['direct_GeV2']
            light=sum(v['direct_GeV2'] for k,v in groups.items() if k!='top')
            width=sum(v['massless_tree_partial_width_GeV'] for v in groups.values())
            rows.append({
                'scale_GeV':scale,'point':label,'s_GeV2':s,'groups':groups,
                'top_GeV2':top,'massless_sum_GeV2':light,'all_GeV2':top+light,
                'A29_Z_relative_top_only':mp.re(top)/(mz2-mw2) if label=='pole' else None,
                'A29_Z_relative_massless_only':mp.re(light)/(mz2-mw2) if label=='pole' else None,
                'A29_Z_relative_all_fermions':mp.re(top+light)/(mz2-mw2) if label=='pole' else None,
                'A29_Z_delta_gprime2_top_only':gp02*mp.re(top)/(mz2-mw2) if label=='pole' else None,
                'A29_Z_delta_gprime2_all_fermions':gp02*mp.re(top+light)/(mz2-mw2) if label=='pole' else None,
                'massless_optical_residue_GeV':mp.im(light)/mp.sqrt(mz2)-width if label=='pole' else None,
                'massless_tree_total_width_GeV':width if label=='pole' else None})
    return {'dps':mp.mp.dps,'g02':g02,'gprime02':gp02,'sW02':1-mw2/mz2,'rows':rows}


def serial(x):
    if isinstance(x,dict):
        return {k:serial(v) for k,v in x.items()}
    if isinstance(x,list):
        return [serial(i) for i in x]
    if isinstance(x,(mp.mpf,mp.mpc)):
        return [mp.nstr(mp.re(x), 52),mp.nstr(mp.im(x), 52)] if isinstance(x,mp.mpc) else mp.nstr(x,52)
    return x


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    parser.add_argument('--inputs',type=Path,default=HERE/'inputs.json')
    opts=parser.parse_args()
    mp.mp.dps=65
    opts.output.write_text(json.dumps(serial(results(json.loads(opts.inputs.read_text()))),indent=2)+'\n')
