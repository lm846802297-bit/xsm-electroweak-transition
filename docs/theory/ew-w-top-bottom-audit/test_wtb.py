import unittest
import mpmath as mp
from wtb import output, load, params, loopfac


class WtbSubset(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.r=output()["rows"]

    def test_below_threshold_domain(self):
        with self.assertRaises(ValueError): params(1,1,1)
        with self.assertRaises(ValueError): params(0,0,1)
        self.assertGreater(load()["Nc"],0)

    def test_independent_reductions_and_small_momentum(self):
        for r in self.r:
            self.assertLess(r["max_two_methods_abs_error"],mp.mpf("1e-59"))
            self.assertLess(r["zero_analytic_check_error"],mp.mpf("1e-59"))
            self.assertLess(r["B0_closed_form_max_error"],mp.mpf("1e-59"))
            self.assertLess(r["zero_slope_analytic_check_error"],mp.mpf("1e-42"))

    def test_evanescent_times_pole(self):
        w=mp.mpf(load()["MW_GeV"])
        for r in self.r:
            expect=r["g0sq_GF_tree"]*load()["Nc"]*loopfac()/9
            found=r["A28_tb_momentum_difference_only"]-r["A28_tb_wrong_if_D4_before_UV_subtraction"]
            # At equal external W scales the zero and on-shell DR terms differ.
            self.assertLess(abs(found+expect),mp.mpf("1e-59"))
            u=r["PiW_tb_UV_pole_coeff_at_pole_GeV2"]-r["PiW_tb_UV_pole_coeff_at_zero_GeV2"]
            self.assertLess(abs(u/(w*w)-3*expect),mp.mpf("1e-59"))

    def test_input_response_cancellation(self):
        w=mp.mpf(load()["MW_GeV"])
        for r in self.r:
            a=r["PiW_tb_MSfinite_at_zero_GeV2"]/(w*w)
            self.assertLess(abs(2*r["dv_over_v_from_tb_W_input_only"]-a),mp.mpf("1e-59"))
            self.assertLess(abs(r["A30_relative_y_t_squared_from_tb_W_input_only_no_top_self_energy"]+a),mp.mpf("1e-59"))


if __name__=="__main__":unittest.main()
