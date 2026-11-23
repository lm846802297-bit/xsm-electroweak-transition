import json
import unittest

import mpmath as mp

from shared_tadpole import HERE, Hessian, fixture, response, result


class SharedTadpoleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.cfg=json.loads((HERE/'inputs.json').read_bytes())
        fixture(cls.cfg)
        cls.a=result(cls.cfg)
        cls.tol=mp.mpf('1e-55')

    def test_tree_masses_and_signed_higgs_rotation(self):
        m=self.a['tree']
        R=mp.matrix([[m['c'],-m['s']],[m['s'],m['c']]])
        D=R*Hessian(m,0,0)*R.T
        expected=mp.diag(m['mh'])
        self.assertLess(max(abs(x-y) for x,y in zip(D,expected)),self.tol)
        self.assertLess(m['s'],0)
        self.assertGreater(m['b1'],0)

    def test_all_signed_diagrams_vs_CW_and_local_response(self):
        for r in self.a['rows']:
            a,b=r['independent_potential']
            self.assertLess(a['maximum_scaled_error'],mp.mpf('1e-36'))
            self.assertLess(b['maximum_scaled_error'],a['maximum_scaled_error'])
            self.assertLess(abs(b['local_linear_response_doublet_GeV']-r['response']['T_h_GeV']),mp.mpf('1e-36'))
        ref=self.a['rows'][0]['diagram_totals_GeV3']
        expected=('-389125.35836266758283616','82611.96159663132390740')
        for a,b in zip(ref,expected):
            self.assertLess(abs(a-mp.mpf(b)),mp.mpf('1e-16'))

    def test_common_projection_nulls_and_vev_signs(self):
        for r in self.a['rows']:
            p=r['response']
            for k in ('A28_tad','A29_tad','A30_tad_sum'):
                self.assertLess(abs(p[k]),self.tol)
            self.assertLess(abs(p['local_minimum_relative_dv']*2+p['pole_matching_relative_dv2']),self.tol)
        for t in (mp.mpf('-10.2'),mp.mpf('0'),mp.mpf('5.49')):
            p=response(self.a['tree'],t)
            self.assertLess(abs(p['W_input_A30']+p['top_A30']),self.tol)
        self.assertGreater(self.a['rows'][0]['response']['top_A30'],0)
        self.assertLess(self.a['rows'][0]['response']['W_input_A30'],0)

    def test_precision(self):
        with mp.workdps(50):
            coarse=result(self.cfg)
        for a,b in zip(self.a['rows'],coarse['rows']):
            self.assertLess(abs(a['response']['T_h_GeV']-b['response']['T_h_GeV']),mp.mpf('1e-43'))


if __name__=='__main__':
    unittest.main()
