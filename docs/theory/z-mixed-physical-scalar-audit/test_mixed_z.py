import unittest

from mpmath import mp
from mixed_z import run


class MixedZ(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = run()

    def test_absolute_reduced_vs_shifted_and_real_pole(self):
        with mp.workdps(65):
            for r in self.data["rows"]:
                for p in r["points"]:
                    self.assertEqual(mp.mpf(p["Pi_Z_imag_GeV2"]), 0)
                    for v in p["states"]:
                        self.assertLess(mp.mpf(v["direct_minus_pv_abs_GeV2"]), mp.mpf("1e-48"))

    def test_independent_positive_unitarity_cut_and_dispersion(self):
        with mp.workdps(65):
            for sp in self.data["spectral_checks"]:
                self.assertEqual(mp.mpf(sp["physical_pole_imag_GeV2"]), 0)
                self.assertGreater(mp.mpf(sp["sample_cut_unitary_imag_GeV2"]), 0)
                self.assertLess(mp.mpf(sp["cut_difference_abs_GeV2"]), mp.mpf("1e-48"))
                self.assertGreater(mp.mpf(sp["Pi_second_derivative_spectral_GeV_minus2"]), 0)
                self.assertLess(mp.mpf(sp["Pi_second_derivative_difference_abs_GeV_minus2"]), mp.mpf("1e-48"))

    def test_a29_sign_and_no_w_input(self):
        with mp.workdps(65):
            d = self.data
            mz = mp.mpf(d["inherited_physics"]["MZ_GeV"])
            mw = mp.mpf(d["inherited_physics"]["MW_GeV"])
            gp2 = mp.mpf(d["gp02"])
            for r in d["rows"]:
                shot = mp.mpf(r["points"][1]["Pi_Z_real_GeV2"])
                rel = mp.mpf(r["A29_Z_only_relative"])
                add = mp.mpf(r["A29_Z_only_additive_gp2"])
                self.assertLess(abs(shot/(mz*mz-mw*mw)-rel), mp.mpf("1e-55"))
                self.assertLess(abs(add-rel*gp2), mp.mpf("1e-55"))
            self.assertLess(mp.mpf(d["rows"][0]["A29_Z_only_relative"]), 0)

    def test_finite_evanescent_term_has_positive_shift_when_omitted(self):
        with mp.workdps(65):
            d = self.data
            gz2 = mp.mpf(d["gZ02"])
            mz = mp.mpf(d["inherited_physics"]["MZ_GeV"])
            for r in d["rows"]:
                for p in r["points"]:
                    s = mp.mpf(p["s_GeV2"])
                    for v in p["states"]:
                        m = mp.mpf(v["m_h_GeV"])
                        w = mp.mpf(v["weight"])
                        e = (3*(m*m+mz*mz)-s)/18
                        self.assertGreater(e, 0)
                        self.assertLess(abs(mp.mpf(v["omitted_B00_evanescent_stripped_GeV2"])-e), mp.mpf("1e-52"))
                        err = mp.mpf(v["wrong_early_D4_real_GeV2"])-mp.mpf(v["pv_GeV2"][0])
                        self.assertLess(abs(err-gz2*w*e/(16*mp.pi**2)), mp.mpf("1e-54"))

    def test_ms_scale_slope(self):
        with mp.workdps(65):
            r0,r1,r2 = self.data["rows"]
            for ri in (r1,r2):
                dl = 2*mp.log(mp.mpf(ri["mu_bar_GeV"])/mp.mpf(r0["mu_bar_GeV"]))
                for pi,p0 in zip(ri["points"],r0["points"]):
                    change = mp.mpf(pi["Pi_Z_real_GeV2"])-mp.mpf(p0["Pi_Z_real_GeV2"])
                    slope = mp.mpf(p0["MS_log_slope_GeV2"])
                    self.assertLess(abs(change-slope*dl),mp.mpf("1e-54"))


if __name__ == "__main__":
    unittest.main()
