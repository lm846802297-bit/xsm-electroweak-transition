import unittest
from mpmath import mp
from seagull import a36, higgs_mass_hessian_half, proper_time, run


class Seagull(unittest.TestCase):
    def setUp(self):
        self.ctx = mp.workdps(65)
        self.ctx.__enter__()
        self.addCleanup(self.ctx.__exit__, None, None, None)

    def test_signed_alternate_regulator_and_cutoff_bound(self):
        for row in run()["rows"]:
            for state in row["states"]:
                for c in state["cutoffs"]:
                    diff = mp.mpf(c["signed_difference_GeV2"])
                    bound = mp.mpf(c["rigorous_error_bound_GeV2"])
                    self.assertTrue(0 < diff <= bound)
                self.assertLess(diff, mp.mpf("1e-48"))

    def test_conversion_negative_control(self):
        m, mu = mp.mpf(350), mp.mpf("91.1876")
        cv, bad, _ = proper_time(m, mu, mp.mpf("1e26"), 65)
        self.assertLess(abs(cv - a36(m, mu)), mp.mpf("1e-48"))
        self.assertLess(abs((cv - bad) - mp.euler * m**2 / (16 * mp.pi**2)), mp.mpf("1e-60"))
        self.assertGreater(abs(cv - bad), 400)

    def test_doublet_hessian_zero_mix_and_degeneracy(self):
        self.assertEqual(higgs_mass_hessian_half(4, [1, 0]), [1, 0])
        for s in map(mp.mpf, ("0", "0.1", "0.9", "1")):
            h = higgs_mass_hessian_half(mp.mpf(4), [mp.sqrt(1 - s*s), -s])
            self.assertLess(abs(sum(k * a36(mp.mpf(350), mp.mpf(91)) for k in h)
                                - a36(mp.mpf(350), mp.mpf(91))), mp.mpf("1e-60"))

    def test_analytic_scale_slope(self):
        data = run()
        r0, r1 = data["rows"][:2]
        diff = mp.mpf(r1["seagull"]["Pi_W_0_GeV2"]) - mp.mpf(r0["seagull"]["Pi_W_0_GeV2"])
        prediction = 2 * mp.log(mp.mpf(r1["mu_bar_GeV"]) / mp.mpf(r0["mu_bar_GeV"])) * mp.mpf(
            data["dPi_seagull_dln_mu2_GeV2"])
        self.assertLess(abs(prediction - diff), mp.mpf("1e-56"))

    def test_w_projectors_and_prior_sum(self):
        data = run()
        mw2 = mp.mpf(data["inherited_physics"]["MW_GeV"])**2
        for row in data["rows"]:
            sea = row["seagull"]
            comb = row["physical_h_internal_W_1PI_bubbles_plus_seagulls"]
            self.assertEqual(sea["Pi_W_0_GeV2"], sea["Pi_W_pole_GeV2"])
            self.assertEqual(mp.mpf(sea["A28_relative"]), 0)
            self.assertLess(abs(mp.mpf(sea["A30_W_input_relative_only"])
                                + mp.mpf(sea["Pi_W_0_GeV2"]) / mw2), mp.mpf("1e-58"))
            self.assertLess(abs(mp.mpf(comb["A30_W_input_relative_only"])
                                + mp.mpf(comb["Pi_W_0_GeV2"]) / mw2), mp.mpf("1e-58"))
        lead = data["rows"][0]["physical_h_internal_W_1PI_bubbles_plus_seagulls"]
        self.assertLess(abs(mp.mpf(lead["A28_relative"])
                            - mp.mpf("0.00023886236871387302539808445219619518053578296701534")),
                        mp.mpf("1e-51"))


if __name__ == "__main__":
    unittest.main()
