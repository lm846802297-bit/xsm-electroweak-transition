import json
import unittest

import mpmath as mp

from mixed_w import HERE, fixture, output, reduced, shifted


class MixedWTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.cfg=json.loads((HERE/'inputs.json').read_bytes())
        cls.pv=staticmethod(fixture(cls.cfg))
        cls.a=output(cls.cfg,cls.pv)
        cls.b=mp.mpf(cls.cfg['MW_GeV'])**2
        cls.tol=mp.mpf('1e-50')

    def test_direct_and_PV_both_momenta_and_A44_zero(self):
        for row in self.a['rows']:
            for p in row['points']:
                for a in p['states']:
                    d=a['reduced_stripped']
                    self.assertLess(abs(a['Pi_difference_GeV2']),self.tol)
                    self.assertLess(abs(d['B0']-a['direct_B0']),self.tol)
                    self.assertLess(abs(d['B00']-a['direct_B00_GeV2']),self.tol)
                    self.assertEqual(mp.im(a['Pi_PV_GeV2']),0)
                    if not p['s_GeV2']:
                        self.assertLess(abs(d['B00_zero_A44']-d['B00']),self.tol)

    def test_log_scale_slope_and_W_only_identity(self):
        first=self.a['rows'][0]
        for row in self.a['rows']:
            log=2*mp.log(mp.mpf(row['mu_bar_GeV'])/mp.mpf(first['mu_bar_GeV']))
            r=row['A28_W_bubble_relative']
            self.assertLess(abs(r-first['A28_W_bubble_relative']-row['A28_log_scale_slope']*log),self.tol)
            self.assertLess(abs(row['A30_W_input_relative_only']+mp.re(row['points'][0]['Pi_GeV2'])/self.b),self.tol)
            for now,ref in zip(row['points'],first['points']):
                for a,b in zip(now['states'],ref['states']):
                    self.assertLess(abs(a['Pi_PV_GeV2']-b['Pi_PV_GeV2']-a['Pi_uv_pole_GeV2']*log),self.tol)

    def test_premature_dimension_negative_control(self):
        exact_delta=self.a['g02']/(18*16*mp.pi**2)
        for row in self.a['rows']:
            true=row['A28_W_bubble_relative']
            early=row['A28_wrong_early_D4']
            self.assertGreater(true,early)
            self.assertLess(abs(true-early-exact_delta),self.tol)

    def test_equal_mass_cut(self):
        a=mp.mpf(25)
        for r in map(mp.mpf,('3.999','4','4.001','10')):
            s=a*r
            b0,b00=shifted(a,a,s,a)
            red=reduced(a,a,s,a,self.pv)
            self.assertLess(abs(b0-red['B0']),self.tol)
            self.assertLess(abs(b00-red['B00']),self.tol)
            beta=mp.sqrt(1-4/r) if r>4 else 0
            self.assertLess(abs(mp.im(b0)-mp.pi*beta),self.tol)
            self.assertLess(abs(mp.im(b00)+mp.pi*s*beta**3/12),self.tol)

    def test_precision_and_source_benchmark(self):
        with mp.workdps(45):
            coarse=output(self.cfg,self.pv)
        for a,b in zip(self.a['rows'],coarse['rows']):
            self.assertLess(abs(a['A28_W_bubble_relative']-b['A28_W_bubble_relative']),mp.mpf('1e-40'))
        a=self.a['rows'][0]['A28_W_bubble_relative']
        self.assertLess(abs(a-mp.mpf('0.000238862368713873025398')),mp.mpf('1e-23'))


if __name__=='__main__':
    unittest.main()
