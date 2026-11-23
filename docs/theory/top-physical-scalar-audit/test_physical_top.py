import json
import unittest

import mpmath as mp

from physical_top import HERE, direct, mix, result, verify_fixture


class PhysicalTopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.cfg=json.loads((HERE/'inputs.json').read_bytes())
        cls.pv=staticmethod(verify_fixture(cls.cfg))
        cls.a=result(cls.cfg,cls.pv)
        cls.tol=mp.mpf('1e-55')

    def test_signed_pv_independent_direct_and_real_physical_poles(self):
        for row in self.a['rows']:
            for h in row['scalars']:
                self.assertLess(abs(h['diff']),self.tol)
                self.assertEqual(mp.im(h['PV']),0)
        self.assertLess(abs(self.a['rows'][0]['A30_relative_1PI']-mp.mpf('-0.006569593417412465597034')),mp.mpf('1e-23'))

    def test_zero_mix_equal_mass_and_common_log(self):
        for row in self.a['rows']:
            a,b=[h['unweighted_A30'] for h in row['scalars']]
            self.assertEqual(mix(a,b,0),a)
            self.assertLess(abs(mix(a,a,self.a['sin2_theta'])-a),self.tol)
            self.assertLess(abs(row['A30_relative_1PI']-a-row['A30_change_from_zero_mixing']),self.tol)
            self.assertLess(abs(row['A30_change_from_zero_mixing']-self.a['rows'][0]['A30_change_from_zero_mixing']),self.tol)

    def test_threshold_with_positive_scalar_weight(self):
        L=mp.mpf('-.5')
        for r in map(mp.mpf,('3.999','4','4.001')):
            b0,b1=self.pv(L,mp.mpf(1),r)
            shifted=direct(L,mp.mpf(1),r)
            self.assertLess(abs(b0-b1-shifted),self.tol)
            cut=mp.mpf('1.5')*mp.pi*mp.sqrt(1-4/r) if r>4 else 0
            self.assertLess(abs(mp.im(shifted)-cut),self.tol)
            if r==4:
                self.assertLess(abs(shifted-mp.mpf('1.5')*(L+2)),self.tol)

    def test_asymptotic_bound_and_mixing_family(self):
        factors=self.a['heavy_limits_mu_bar_MZ']
        fixed=[]; suppressed=[]
        for r in factors:
            self.assertGreater(r['asymptotic_minus_direct'],0)
            self.assertLessEqual(r['asymptotic_minus_direct'],r['bound'])
            fixed.append(abs(r['A30_difference_fixed_mixing']))
            suppressed.append(abs(r['A30_difference_scaling_mixing']))
        self.assertEqual(fixed,sorted(fixed))
        self.assertEqual(suppressed,sorted(suppressed,reverse=True))
        self.assertLess(suppressed[-1],suppressed[0]/1000)

    def test_precision(self):
        with mp.workdps(45):
            coarse=result(self.cfg,self.pv)
        for a,b in zip(self.a['rows'],coarse['rows']):
            self.assertLess(abs(a['A30_relative_1PI']-b['A30_relative_1PI']),mp.mpf('1e-40'))


if __name__=='__main__':
    unittest.main()
