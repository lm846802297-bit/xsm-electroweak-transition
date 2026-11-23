"""NST (A58) h_i-containing W seagulls and matched immutable bubble sum."""
import hashlib
import json
import sys
from pathlib import Path
from decimal import Decimal

from mpmath import mp

HERE = Path(__file__).resolve().parent


def locked_json(spec):
    data = (HERE / spec["path"]).read_bytes()
    if hashlib.sha256(data).hexdigest() != spec["sha256"]:
        raise ValueError(f"Prior input hash mismatch: {spec['path']}")
    return json.loads(data)


def a36(m, mu):
    return m**2 * (1 + 2 * mp.log(mu / m)) / (16 * mp.pi**2)


def proper_time(m, mu, ratio, guard):
    """Euclidean coincident propagator, then convert to the sign of A36."""
    with mp.extradps(guard):
        cutoff2 = (m * ratio)**2
        z = m**2 / cutoff2
        raw = cutoff2 * mp.exp(-z) - m**2 * mp.e1(z)
        msbar_ct = cutoff2 + m**2 * (mp.log(mu**2 / cutoff2) + mp.euler)
        converted = -(raw - msbar_ct) / (16 * mp.pi**2)
        without_gamma = converted - mp.euler * m**2 / (16 * mp.pi**2)
        bound = m**2 * z / (32 * mp.pi**2)
    return +converted, +without_gamma, +bound


def higgs_mass_hessian_half(g2, signed_doublet_axis):
    """The coincidence contraction uses only diagonal Hessian entries."""
    return [(g2 / 4) * item**2 for item in signed_doublet_axis]


def number(value):
    return mp.nstr(value, 59)


def run():
    conf = json.loads((HERE / "inputs.json").read_text())
    inp = locked_json(conf["prior_inputs"])
    bub = locked_json(conf["prior_bubbles"])
    if inp != bub["inputs"]:
        raise ValueError("Prior input and bubble fixture disagree")
    if inp["gauge_xi"] != "1":
        raise ValueError("This partial audit fixes Feynman gauge")
    with mp.workdps(conf["working_dps"]):
        mw = mp.mpf(inp["MW_GeV"])
        g2 = 4 * mp.sqrt(2) * mp.mpf(inp["GF_GeV_minus2"]) * mw**2
        s = mp.mpf(inp["sin_theta"])
        masses = [mp.mpf(inp[k]) for k in ("Mh1_GeV", "Mh2_GeV")]
        weights = [1 - s**2, s**2]
        slope = sum(g2 * w * m**2 / (64 * mp.pi**2) for w, m in zip(weights, masses))
        rows = []
        if len(inp["mu_bar_GeV"]) != len(bub["rows"]):
            raise ValueError("Scale count mismatch")
        for mu_text, prior in zip(inp["mu_bar_GeV"], bub["rows"]):
            if mp.mpf(prior["mu_bar_GeV"]) != mp.mpf(mu_text):
                raise ValueError("Scale mismatch")
            mu = mp.mpf(mu_text)
            subtotal = mp.mpf(0)
            states = []
            for idx, (w, m) in enumerate(zip(weights, masses)):
                eq = g2 * w * a36(m, mu) / 4
                subtotal += eq
                cutoffs = []
                for rv in conf["cutoff_to_mass_ratios"]:
                    # Reparse decimal fixtures with guards; the last bound gap
                    # is O((m/Lambda)^4), not resolvable in the 65-dps main sum.
                    with mp.extradps(conf["cutoff_guard_dps"]):
                        mh = mp.mpf(inp[("Mh1_GeV", "Mh2_GeV")[idx]])
                        sh = mp.mpf(inp["sin_theta"])
                        muh = mp.mpf(mu_text)
                        g2h = 4 * mp.sqrt(2) * mp.mpf(inp["GF_GeV_minus2"]) * mp.mpf(inp["MW_GeV"])**2
                        hc = higgs_mass_hessian_half(g2h, [mp.sqrt(1 - sh**2), -sh])[idx]
                        eqh = g2h * ((1 - sh**2) if idx == 0 else sh**2) * a36(mh, muh) / 4
                        converted, bad, bound = proper_time(mh, muh, mp.mpf(rv), conf["cutoff_guard_dps"])
                        cutoffs.append({
                            "Lambda_over_m": rv, "Pi_GeV2": number(hc * converted),
                            "signed_difference_GeV2": number(hc * converted - eqh),
                            "rigorous_error_bound_GeV2": number(hc * bound),
                            "Pi_wrong_no_gamma_conversion_GeV2": number(hc * bad),
                        })
                states.append({"mass_GeV": str(m), "w": str(w),
                               "Pi_A58_GeV2": number(eq), "cutoffs": cutoffs})
            zero, pole = prior["points"]
            if mp.mpf(zero["s_GeV2"]) != 0 or Decimal(pole["s_GeV2"]) != Decimal(inp["MW_GeV"])**2:
                raise ValueError("Prior momentum points do not match")
            if mp.mpf(pole["Pi_GeV2"][1]) != 0:
                raise ValueError("Prior physical pole unexpectedly has an absorptive part")
            b0, bp = mp.mpf(zero["Pi_GeV2"]), mp.mpf(pole["Pi_GeV2"][0])
            mix0 = g2 * a36(masses[0], mu) / 4
            rows.append({
                "mu_bar_GeV": mu_text, "states": states,
                "seagull": {"Pi_W_0_GeV2": number(subtotal), "Pi_W_pole_GeV2": number(subtotal),
                            "A28_relative": "0.0",
                            "A30_W_input_relative_only": number(-subtotal / mw**2),
                            "zero_mix_Pi_GeV2": number(mix0),
                            "mixing_increment_GeV2": number(subtotal - mix0)},
                "physical_h_internal_W_1PI_bubbles_plus_seagulls": {
                    "Pi_W_0_GeV2": number(b0 + subtotal), "Pi_W_pole_GeV2": number(bp + subtotal),
                    "A28_relative": number((bp - b0) / mw**2),
                    "A30_W_input_relative_only": number(-(b0 + subtotal) / mw**2),
                    "prior_rounding_floor_for_totals_GeV2": "1e-48",
                },
            })
        return {"configuration": conf, "inherited_physics": inp,
                "g02": number(g2), "dPi_seagull_dln_mu2_GeV2": number(slope), "rows": rows}


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("native-results.json")
    result = run()
    out.write_text(json.dumps(result, indent=2) + "\n")
    for row in result["rows"]:
        print(row["mu_bar_GeV"], row["seagull"]["Pi_W_0_GeV2"],
              row["seagull"]["A30_W_input_relative_only"])
