import json
import unittest

import mpmath as mp

from charged_top import HERE, direct_moment, moments, result


class ChargedTopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.inp=json.loads((HERE/'inputs.json').read_bytes())
        cls.a=result(cls.inp)
        cls.tol=mp.mpf('1e-55')

    def test_full_complex_projection_and_evanescent_term(self):
        for r in self.a['rows']:
            for key in ('W_direct_minus_PV','Goldstone_direct_minus_PV','B1_minus_shifted_projection'):
                self.assertLess(abs(r[key]),self.tol)
            delta=r['A30_charged_relative_early_D4']-r['A30_charged_relative']
            self.assertLess(abs(delta-self.a['g02']/(32*mp.pi**2)),self.tol)

    def test_equal_invariant_threshold_and_cut(self):
        L=mp.mpf('-0.5')
        self.assertEqual(mp.re(moments(L,mp.mpf(1))[1]),-L/2-mp.mpf('.5'))
        for b in map(mp.mpf,('1.01','1','0.99999','0.80')):
            pv=moments(L,b)[1]
            self.assertLess(abs(pv+direct_moment(L,b)),self.tol)
            self.assertLess(abs(mp.im(pv)+(mp.pi/2)*max(1-b,0)**2),self.tol)

    def test_optical_versus_physical_W_spin_sum(self):
        for r in self.a['rows']:
            self.assertGreater(mp.im(r['W_PV']),0)
            self.assertGreater(mp.im(r['Goldstone_PV']),0)
            self.assertLess(abs(r['bW_width_from_optical_GeV']-r['bW_width_from_spin_sum_GeV']),self.tol)

    def test_precision_and_reference_scale(self):
        with mp.workdps(45):
            coarse=result(self.inp)
        for r,a in zip(coarse['rows'],self.a['rows']):
            for key in ('W_PV','Goldstone_PV','A30_charged_relative','bW_width_from_optical_GeV'):
                self.assertLess(abs(r[key]-a[key]),mp.mpf('1e-40'))
        self.assertLess(abs(self.a['rows'][0]['A30_charged_relative']-mp.mpf('0.004217220461587993106')),mp.mpf('1e-22'))


if __name__=='__main__':
    unittest.main()
