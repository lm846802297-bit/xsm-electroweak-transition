import unittest
import mpmath as mp
from tadpole import radial, bessel, highT2, bound_R, split, example, load


class ScalarTadpole(unittest.TestCase):
    def test_massless_and_monotone(self):
        self.assertEqual(radial(0)[0],1/12)
        self.assertEqual(highT2(0),1/12)
        f=[radial(x)[0] for x in load()["x_scan"]]
        self.assertTrue(all(a>b>0 for a,b in zip(f,f[1:])))
        for x in (-0.2,float("nan"),float("inf")):
            with self.assertRaises(ValueError): radial(x)

    def test_three_routes_and_ms_scale_cancel(self):
        for x in load()["x_scan"]:
            f=radial(x)[0]; B,b,_=bessel(x)
            self.assertLessEqual(abs(B-f),b+1e-14)
            lo,hi=bound_R(x)
            self.assertLessEqual(lo-1e-14,f-highT2(x))
            self.assertLessEqual(f-highT2(x),hi+1e-14)
            for mu in (.8,1,4*mp.pi):
                low,high=split(x,mu)["thermal_reconstructed_interval"]
                self.assertLessEqual(low-1e-14,f)
                self.assertLessEqual(f,high+1e-14)

    def test_high_precision_integral(self):
        with mp.workdps(55):
            for x in map(mp.mpf,("0.05","0.5","3")):
                def f(q):
                    E=mp.sqrt(q*q+x*x)
                    return q*q/E/mp.expm1(E)/(2*mp.pi**2)
                independent=mp.quad(f,[0,1,3,10,mp.inf])
                self.assertLess(abs(independent-radial(x)[0]),mp.mpf("2e-15"))

    def test_exact_cancellation_requires_equal_masses(self):
        a,b,c=[example(x) for x in load()["examples"]]
        self.assertEqual(a["thermal_scalar_total"],0)
        self.assertEqual(b["massless_eq_2_8_local"],0)
        self.assertGreater(b["thermal_scalar_total"],0)
        self.assertLess(c["thermal_scalar_total"],c["massless_eq_2_8_local"])


if __name__=="__main__": unittest.main()
