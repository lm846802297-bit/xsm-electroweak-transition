import unittest
import mpmath as mp
from threshold import analytical, integrals, determinant, shifted, points, matching, rel


class TemporalThresholdTests(unittest.TestCase):
    def setUp(self): mp.mp.dps = 65

    def test_second_stage_matches_momentum_integrals(self):
        for row in points():
            c, s = matching(row), mp.mpf(row["s"])
            self.assertLess(rel(analytical(c,s),integrals(c,s)),mp.mpf("1e-59"))

    def test_field_shift_and_matrix_degeneracy(self):
        for row in points():
            c, s, d = matching(row), mp.mpf(row["s"]), mp.mpf(row["shift"])
            self.assertLess(rel(analytical(c,s),analytical(shifted(c,d),s-d)),mp.mpf("1e-59"))
            for rho in (mp.mpf(0),mp.mpf(row["rho"])):
                v = determinant(c,s,rho)
                w = determinant(shifted(c,d),s-d,rho)
                self.assertLess(abs(v-w)/max(1,abs(v)),mp.mpf("1e-59"))

    def test_determinant_derivatives_include_neutral_mixing(self):
        for row in points():
            c, s = matching(row), mp.mpf(row["s"])
            v=determinant
            a=analytical(c,s)
            direct={"b1":mp.diff(lambda z:v(c,z,0),s),
                    "b2":mp.diff(lambda z:v(c,z,0),s,2),
                    "cH":mp.diff(lambda rho:v(c,s,rho),0),
                    "a1":2*mp.diff(lambda z,rho:v(c,z,rho),(s,0),(1,1)),
                    "a2":mp.diff(lambda z,rho:v(c,z,rho),(s,0),(2,1)),
                    "lambda":mp.diff(lambda rho:v(c,s,rho),0,2)/2}
            self.assertLess(rel(direct,a),mp.mpf("1e-59"))

    def test_factor_two_is_not_full_curvature_factor(self):
        c = matching(points()[0]); p = matching(points()[0],printed=True)
        a,b=analytical(c,0),analytical(p,0)
        for k,q in {"b1":2,"b2_x2":4,"a1":2,"a2_x2":4,"b2_y":1,"a2_y":1,"lambda":1}.items():
            self.assertLess(abs(a[k]*q-b[k]),mp.mpf("1e-59"))
        self.assertNotEqual(a["b2"]*2,b["b2"])

if __name__=="__main__":
    unittest.main()
