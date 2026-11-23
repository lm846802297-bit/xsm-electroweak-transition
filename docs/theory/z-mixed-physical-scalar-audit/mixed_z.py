"""NST A59 mixed h_i-Z/neutral-Goldstone self-energy loops, Feynman gauge."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

from mpmath import mp

HERE = Path(__file__).resolve().parent


def pinned_bytes(spec):
    p = HERE / spec["path"]
    data = p.read_bytes()
    if hashlib.sha256(data).hexdigest() != spec["sha256"]:
        raise ValueError("Pinned input mismatch: " + p.name)
    return p, data


def backend(spec):
    path, _ = pinned_bytes(spec)
    src = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(src)
    src.loader.exec_module(mod)
    return mod


def shifted(a, b, s, mu2):
    """Absolute finite loop via D-dimensional shift and transverse IBP, stripped."""
    mid = mp.mpf("0.5")
    sp = (mp.sqrt(a) + mp.sqrt(b))**2
    sm = (mp.sqrt(a) - mp.sqrt(b))**2
    if s > sp:
        term = mp.sqrt((s - sp)*(s - sm))
        l, r = (s+a-b-term)/(2*s), (s+a-b+term)/(2*s)
        edges = [0, l, r, 1]

        def delta(x):
            return s*(x-l)*(x-r)
    else:
        edges = [0, mid, 1]

        def delta(x):
            return (1-x)*a + x*b - s*x*(1-x)

    def integrate(x):
        q = delta(x)
        if not q:  # value at isolated endpoints does not affect the quadrature
            return mp.mpf(0)
        log = mp.log(mu2/abs(q)) + (mp.j*mp.pi if q < 0 else 0)
        return b*log - q*(1+log)/2
    return mp.quad(integrate, edges)


def unitary_imag(a, b, s, gz2w):
    """Imaginary part via physical h_i+Z phase space and polarization sum."""
    sp = (mp.sqrt(a) + mp.sqrt(b))**2
    sm = (mp.sqrt(a) - mp.sqrt(b))**2
    if s <= sp:
        return mp.mpf(0)
    beta = mp.sqrt((s-sp)*(s-sm))/s
    p2 = s*beta**2/4
    return gz2w*b*beta*(1+p2/(3*b))/(16*mp.pi)


def second_direct(a, b, s, c):
    def integrand(x):
        u = x*(1-x)
        q = (1-x)*a + x*b - s*u
        return u*u*(b/q**2+1/(2*q))
    return c*mp.quad(integrand, [0, mp.mpf("0.5"), 1])


def second_spectral(a, b, s, c):
    """Above-cut density with t=threshold/u; two subtractions, no UV fit."""
    sp = (mp.sqrt(a) + mp.sqrt(b))**2
    sm = (mp.sqrt(a) - mp.sqrt(b))**2

    def integrand(u):
        beta = mp.sqrt((1-u)*(1-u*sm/sp))
        return (b*u*beta+sp*beta**3/12)*sp/(sp-s*u)**3
    return 2*c*mp.quad(integrand, [0, mp.mpf("0.5"), 1])


def serial(x):
    if isinstance(x, dict):
        return {key: serial(value) for key, value in x.items()}
    if isinstance(x, list):
        return [serial(value) for value in x]
    if isinstance(x, mp.mpc):
        return [serial(x.real), serial(x.imag)]
    if isinstance(x, mp.mpf):
        return mp.nstr(x, 59)
    return x


def run():
    conf = json.loads((HERE / "inputs.json").read_text())
    inp = json.loads(pinned_bytes(conf["benchmark"])[1])
    if inp["gauge_xi"] != "1":
        raise ValueError("This projection requires Feynman gauge")
    reduce = backend(conf["a41_a42_pv_backend"]).reduced
    analytic_pv = backend(conf["a41_b0_b1_backend"]).pv_moments
    with mp.workdps(conf["working_dps"]):
        gf, mw, mz, m1, m2, st = [mp.mpf(inp[k]) for k in
            ("GF_GeV_minus2", "MW_GeV", "MZ_GeV", "Mh1_GeV", "Mh2_GeV", "sin_theta")]
        b, w2 = mz*mz, mw*mw
        g2, gz2 = [4*mp.sqrt(2)*gf*x for x in (w2,b)]
        gp2 = gz2-g2
        L = 16*mp.pi**2
        states = [(m1,1-st*st),(m2,st*st)]
        rows = []
        for scale in inp["mu_bar_GeV"]:
            mu2 = mp.mpf(scale)**2
            points = []
            for s in (mp.mpf(0),b):
                values = []
                for m, w in states:
                    a = m*m
                    c = gz2*w/L
                    pv = reduce(a,b,s,mu2,analytic_pv)
                    pi = c * (b*pv["B0"] - pv["B00"])
                    direct = c * shifted(a,b,s,mu2)
                    wrong = c * (b*pv["B0"] - pv["B00_early_D4"])
                    values.append({
                        "m_h_GeV": m, "weight": w, "direct_GeV2": [mp.re(direct),mp.im(direct)],
                        "pv_GeV2": [mp.re(pi),mp.im(pi)], "direct_minus_pv_abs_GeV2": abs(direct-pi),
                        "B0_stripped": pv["B0"], "B00_stripped_GeV2": pv["B00"],
                        "omitted_B00_evanescent_stripped_GeV2": pv["B00_evanescent"],
                        "wrong_early_D4_real_GeV2": mp.re(wrong),
                        "MS_log_slope_GeV2": c*(b-(a+b)/4+s/12)})
                points.append({"s_GeV2":s,"states":values,
                               "Pi_Z_real_GeV2":sum(v["pv_GeV2"][0] for v in values),
                               "Pi_Z_imag_GeV2":sum(v["pv_GeV2"][1] for v in values),
                               "wrong_early_D4_real_GeV2":sum(v["wrong_early_D4_real_GeV2"] for v in values),
                               "MS_log_slope_GeV2":sum(v["MS_log_slope_GeV2"] for v in values)})
            p = points[1]
            rel = p["Pi_Z_real_GeV2"]/(b-w2)
            rows.append({"mu_bar_GeV": scale, "points": points, "A29_Z_only_relative": rel,
                         "A29_Z_only_additive_gp2": g2/w2*p["Pi_Z_real_GeV2"],
                         "A29_Z_only_wrong_early_D4_relative":p["wrong_early_D4_real_GeV2"]/(b-w2)})
        spectra = []
        for m,w in states:
            a=m*m
            c=gz2*w/L
            sp=(m+mz)**2
            probe=mp.mpf(conf["above_threshold_s_ratio"])*sp
            sd=second_direct(a,b,b,c)
            ss=second_spectral(a,b,b,c)
            cd=mp.im(c*shifted(a,b,probe,b))
            cs=unitary_imag(a,b,probe,gz2*w)
            spectra.append({"m_h_GeV":m,"threshold_s_GeV2":sp,
                            "above_cut_probe_s_GeV2":probe,
                            "physical_pole_imag_GeV2":unitary_imag(a,b,b,gz2*w),
                            "sample_cut_shifted_imag_GeV2":cd,"sample_cut_unitary_imag_GeV2":cs,
                            "cut_difference_abs_GeV2":abs(cd-cs),
                            "Pi_second_derivative_shifted_GeV_minus2":sd,
                            "Pi_second_derivative_spectral_GeV_minus2":ss,
                            "Pi_second_derivative_difference_abs_GeV_minus2":abs(sd-ss)})
        return serial({"configuration":conf,"inherited_physics":inp,
                       "g02":g2,"gZ02":gz2,"gp02":gp2,"rows":rows,"spectral_checks":spectra})


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("native-results.json")
    n = run()
    path.write_text(json.dumps(n, indent=2)+"\n")
    for row in n["rows"]:
        print(row["mu_bar_GeV"], row["points"][0]["Pi_Z_real_GeV2"],
              row["points"][-1]["Pi_Z_real_GeV2"], row["A29_Z_only_relative"])
