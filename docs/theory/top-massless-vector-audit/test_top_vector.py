import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import mpmath as mp

from top_vector import HERE, fixture_check, output, pv, trace


class TopVectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.cfg=json.loads((HERE/'inputs.json').read_bytes())
        cls.a=output(cls.cfg)
        cls.eps=mp.mpf('1e-60')

    def test_pv_trace_analytic_and_early_D4(self):
        for r in self.a['rows']:
            self.assertLess(abs(r['trace_minus_A61']),self.eps)
            self.assertLess(abs(r['combination_finite']-r['analytic_finite']),self.eps)
            self.assertEqual(r['combination_pole'],3)
            self.assertEqual(r['evanescent'],-1)
            self.assertLess(abs(r['combination_early_D4']-r['combination_finite']-1),self.eps)

    def test_pole_scale_standard_QCD_and_QED(self):
        r=self.a['rows'][-1]
        q=mp.mpf(self.cfg['Qt'])
        exp_qcd=-mp.mpf(4)/3*self.a['alpha_s_fixed_NST']/mp.pi
        exp_qed=-q*q*self.a['e02']/(4*mp.pi**2)
        self.assertLess(abs(r['Sigma_sum_gluon']-exp_qcd),self.eps)
        self.assertLess(abs(r['Sigma_sum_photon']-exp_qed),self.eps)
        self.assertLess(abs(r['A30_yt2_relative_gluon_plus_photon_only']-2*(exp_qcd+exp_qed)),self.eps)

    def test_scale_shift_and_zero(self):
        first=self.a['rows'][0]
        log=2*mp.log(mp.mpf(self.a['one_loop_zero_GeV'])/mp.mpf(self.cfg['Mt_GeV']))
        self.assertLess(abs(pv(log)['combination_finite']),self.eps)
        for r in self.a['rows']:
            shift=r['A30_yt2_relative_gluon_plus_photon_only']-first['A30_yt2_relative_gluon_plus_photon_only']
            expect=(r['L']-first['L'])*r['A30_yt2_relative_UV_pole_coeff']
            self.assertLess(abs(shift-expect),self.eps)

    def test_threshold_limit_from_below_real_and_convergent(self):
        target=trace(mp.mpf(0))
        errors=[]
        for n in (4,8,12):
            value=trace(mp.mpf(0),1-mp.mpf(10)**(-n))
            self.assertEqual(mp.im(value),0)
            errors.append(abs(value-target))
        self.assertGreater(errors[0],errors[1])
        self.assertGreater(errors[1],errors[2])
        self.assertLess(errors[2],mp.mpf('1e-10'))
        with self.assertRaises(ValueError):
            trace(mp.mpf(0), mp.mpf('1.01'))

    def test_archived_input_hash_fail_closed(self):
        path=HERE/self.cfg['src_WZ_inputs']['path']
        fixture_check(self.cfg,path)
        with TemporaryDirectory() as tmp:
            bad=Path(tmp)/'copy.json'
            bad.write_bytes(path.read_bytes()+b' ')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):
                fixture_check(self.cfg,bad)


if __name__=='__main__':
    unittest.main()
