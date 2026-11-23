import json
import unittest

import mpmath as mp

from zferm import HERE, b0_closed, direct, pole_reduced, pv, results


class NeutralFermionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.cfg=json.loads((HERE/'inputs.json').read_text())
        cls.result=results(cls.cfg)
        cls.z=mp.mpf(cls.cfg['MZ_GeV'])
        cls.m2=mp.mpf(cls.cfg['Mt_GeV'])**2
        cls.s=cls.z**2
        cls.tol=mp.mpf('1e-60')

    def test_independent_pv_at_all_scales_and_momenta(self):
        for row in self.result['rows']:
            for group,item in row['groups'].items():
                with self.subTest(scale=row['scale_GeV'], point=row['point'], group=group):
                    self.assertLess(abs(item['difference_GeV2']),self.tol)

    def test_zero_momentum_vector_ward_and_axial_a44(self):
        a=mp.mpf('0.5')
        mu2=self.s
        L=mp.log(mu2/self.m2)
        # A44: 2 B00=(2/D)(A+m²B0). Its evanescent finite term is m²/2.
        b00_twice_a44=(self.m2*(1+L)+self.m2*L)/2+self.m2/2
        other=-self.m2*(1+L)-2*self.m2*L
        pure_axial_pv=a*a*(b00_twice_a44+other)
        self.assertEqual(direct(self.m2,mp.mpf(0),mu2,1,0),0)
        self.assertLess(abs(direct(self.m2,mp.mpf(0),mu2,0,a)-pure_axial_pv),self.tol)
        self.assertLess(abs(pure_axial_pv+2*a*a*self.m2*L),self.tol)

    def test_massless_closed_form_and_optical_sign(self):
        for row in self.result['rows'][1::2]:
            mu2=mp.mpf(row['scale_GeV'])**2
            for group, item in row['groups'].items():
                if group=='top':
                    self.assertEqual(mp.im(item['direct_GeV2']),0)
                    continue
                expected=item['pole_coeff_GeV2']*(mp.log(mu2/self.s)+mp.mpf(5)/3+1j*mp.pi)
                pi=item['direct_GeV2']
                self.assertLess(abs(pi-expected),self.tol)
                self.assertGreater(mp.im(pi),0)
                self.assertLess(abs(mp.im(pi)/self.z-item['massless_tree_partial_width_GeV']),self.tol)

    def test_scale_difference_is_its_uv_pole_coefficient(self):
        for i in (0,1):
            row1,row2=self.result['rows'][i],self.result['rows'][i+2]
            log=2*mp.log(mp.mpf(row2['scale_GeV'])/mp.mpf(row1['scale_GeV']))
            for key,a in row1['groups'].items():
                b=row2['groups'][key]
                self.assertLess(abs(b['direct_GeV2']-a['direct_GeV2']-log*a['pole_coeff_GeV2']),self.tol)

    def test_early_d4_loses_the_predicted_finite_pole_term(self):
        q,t3=mp.mpf(2)/3,mp.mpf('0.5')
        sw2=self.result['sW02']
        af=4*q*sw2*(2*q*sw2-2*t3)
        for s in (mp.mpf(0),self.s):
            discrepancy=pv(self.m2,s,self.s,q,t3,sw2)-pv(self.m2,s,self.s,q,t3,sw2,True)
            self.assertLess(abs(discrepancy-(1+af)*(3*self.m2-s/2)/9),self.tol)
        self.assertEqual(pole_reduced(self.m2,0,1,0),0)
        with self.assertRaises(ValueError):
            b0_closed(self.m2,4*self.m2,self.s)
        with self.assertRaises(ValueError):
            direct(self.m2,4*self.m2,self.s,q,t3)


if __name__ == '__main__':
    unittest.main()
