import unittest
import mpmath as mp
from screened import load, result, evaluate, scaled


class StaticRingAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.r=result()

    def test_positive_domain(self):
        m=load()["mode_data"][0]
        with self.assertRaises(ValueError): evaluate(m,0,"unknown")
        with self.assertRaises(ValueError): evaluate(dict(m,P0="-3600"),0,"field")
        with self.assertRaises(ValueError): evaluate(dict(m,X0="0"),0,"field")

    def test_zero_insertion_and_double_count(self):
        zero,constant,variable=self.r["records"]
        self.assertEqual(zero["totals"]["ring_only"],0)
        for c in (zero,constant,variable):
            for v in c["modes"].values():
                self.assertLess(scaled(v["wrong_add_without_subtract"]-v["screened_once"],v["bare_static"]),mp.mpf("1e-60"))
        self.assertEqual(constant["totals"]["wrong_omit_Pi_prime"],constant["totals"]["screened_once"])

    def test_slope_term_matches_chain_rule(self):
        T=mp.mpf(load()["T"])
        for R in self.r["records"]:
            for d in R["modes"].values():
                n=d["multiplicity"]
                ref=-n*T*d["Mscreened"]*d["Pi_prime"]/(8*mp.pi)
                self.assertLess(scaled(d["screened_once"]-d["wrong_omit_Pi_prime"],ref),mp.mpf("1e-60"))

    def test_explicit_insertion_removal(self):
        for o in self.r["overlap_control"].values():
            self.assertLess(scaled(o["coefficient"],o["derivative_at_zero"]),mp.mpf("1e-60"))
            rows=o["rows"]
            for a,b in zip(rows,rows[1:]):
                f=a["subtract_insertion_remainder"]/b["subtract_insertion_remainder"]
                self.assertLess(abs(f-4),mp.mpf(".01"))

    def test_independent_routes(self):
        for c in self.r["records"]:
            self.assertLess(c["scaled_radial"],mp.mpf("1e-59"))
            self.assertLess(c["scaled_fivepoint"],mp.mpf("1e-43"))


if __name__=="__main__": unittest.main()
