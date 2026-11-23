import unittest
import mpmath as mp
from pole import result, load, b_ms


class ScalarVacuumPole(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=70
        cls.r=result()

    def test_real_pole_and_heavy_threshold(self):
        for row in self.r["partial_mass_diagrams"]:
            l,h=row["rows"]
            self.assertEqual(l["Pi_partial"].imag,0)
            self.assertGreater(h["Pi_partial"].imag,0)
            self.assertEqual(h["B0_integral_by_index"]["12"].imag,0)
            self.assertGreater(h["B0_integral_by_index"]["11"].imag,0)

    def test_independent_bubble_and_static_CW(self):
        for r in self.r["partial_mass_diagrams"]:
            self.assertLess(r["max_B0_integral_root_difference"],mp.mpf("1e-57"))
            for x in r["rows"]:
                self.assertLess(x["CW_gradient_scaled_error"],mp.mpf("1e-40"))
                self.assertLess(x["CW_curvature_scaled_error"],mp.mpf("1e-35"))

    def test_tadpoles_are_numerically_nonzero(self):
        x=self.r["partial_mass_diagrams"][0]["rows"]
        self.assertTrue(all(abs(i["1PR_tadpole"])>1 for i in x))
        for i in x:
            self.assertEqual(i["Pi_partial"],i["1PI"]+i["1PR_tadpole"])

    def test_a28_is_sign_sentinel_not_full_gauge_match(self):
        x=self.r["a28_explicit_reference_not_full_EW"]
        self.assertLess(x["signed_explicit_bracket_contribution"],0)
        self.assertTrue(x["PiW_momentum_difference_uncomputed"])

    def test_nst_angle_convention(self):
        self.assertEqual(self.r["phase_angle_NST"],mp.mpf("-0.1"))


if __name__=="__main__":unittest.main()
