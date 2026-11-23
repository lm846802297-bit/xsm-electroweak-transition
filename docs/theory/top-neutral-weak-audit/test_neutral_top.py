import json
import unittest

import mpmath as mp

from neutral_top import HERE, fixture, output, pv_moments, signed_trace_integrals


class NeutralTopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.cfg=json.loads((HERE/'inputs.json').read_bytes())
        cls.a=output(cls.cfg)
        cls.tol=mp.mpf('1e-55')

    def test_full_signed_projection_both_routes_and_couplings(self):
        a=self.a
        cL,cR=a['cL'],a['cR']
        self.assertLess(abs(a['A61_mass_coefficient']+cL*cR),self.tol)
        self.assertLess(abs(a['A61_vector_coefficient']-(cL*cL+cR*cR)/2),self.tol)
        for row in a['rows']:
            self.assertLess(abs(row['Z_direct_minus_PV']),self.tol)
            self.assertLess(abs(row['G0_direct_minus_PV']),self.tol)
            self.assertEqual(mp.im(row['Z_PV']),0)
            self.assertEqual(mp.im(row['G0_PV']),0)

    def test_evanescent_finite_shift(self):
        cL,cR=self.a['cL'],self.a['cR']
        eps=(2*cL*cR-(cL*cL+cR*cR)/2)/(16*mp.pi**2)
        for row in self.a['rows']:
            self.assertLess(abs(row['Z_evanescent']-eps),self.tol)
            self.assertLess(abs(row['A30_relative_neutral_only']-row['A30_relative_early_D4']-2*eps),self.tol)

    def test_equal_mass_threshold_and_causal_cut(self):
        L=mp.mpf('-.5')
        for r in map(mp.mpf,('3.999','4','4.001','10')):
            b0,b1=pv_moments(L,mp.mpf(1),r)
            expect=mp.pi*mp.sqrt(1-4/r) if r>4 else 0
            self.assertLess(abs(mp.im(b0)-expect),self.tol)
            self.assertLess(abs(b1+b0/2),self.tol)
            if r==4:
                self.assertLess(abs(b0-(L+2)),self.tol)
            cL,cR,y=self.a['cL'],self.a['cR'],self.a['yt02']
            zd,gd=signed_trace_integrals(L,mp.mpf(1),r,cL,cR,y)
            self.assertLess(abs(zd+4*cL*cR*b0+(cL*cL+cR*cR)*b1-2*cL*cR+(cL*cL+cR*cR)/2),self.tol)
            self.assertLess(abs(gd+y*(b0+b1)/2),self.tol)
        # At the pseudothreshold the double zero is outside [0,1].
        self.assertEqual(mp.im(pv_moments(L,mp.mpf('.25'),mp.mpf('.25'))[0]),0)

    def test_precision_benchmark_and_prior_input(self):
        fixture(self.cfg,HERE/self.cfg['fixture']['path'])
        with mp.workdps(45):
            coarse=output(self.cfg)
        for a,b in zip(self.a['rows'],coarse['rows']):
            for key in ('Z_PV','G0_PV','Z_direct','G0_direct','A30_relative_neutral_only'):
                self.assertLess(abs(a[key]-b[key]),mp.mpf('1e-40'))
        self.assertLess(abs(self.a['rows'][0]['A30_relative_neutral_only']-mp.mpf('0.00016375310025433244827')),mp.mpf('1e-22'))


if __name__=='__main__':
    unittest.main()
