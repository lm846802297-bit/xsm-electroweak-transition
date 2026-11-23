from fractions import Fraction
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import mpmath as mp

from fermion_gauge import HERE, WTB, Z, projections, read_pinned, run


class FermionGaugeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=65
        cls.cfg=json.loads((HERE/'inputs.json').read_bytes())
        cls.out=run(cls.cfg,WTB,Z)
        cls.rows=cls.out['rows']
        cls.w2,cls.z2=(mp.mpf(cls.cfg[k])**2 for k in ('MW_GeV','MZ_GeV'))
        cls.eps=mp.mpf('1e-48')

    def test_trace_pv_scaleless_width_and_early_D4(self):
        g2=self.out['g02']
        n=self.out['W_massless_multiplicity']
        for r in self.rows:
            self.assertEqual(r['W_light_zero_GeV2'],0)
            self.assertLess(abs(r['W_light_direct_minus_PV_GeV2']),mp.mpf('1e-60'))
            self.assertGreater(mp.im(r['W_light_direct_pole_GeV2']),0)
            self.assertLess(abs(r['W_light_tree_width_GeV']-r['W_light_width_from_Impi_GeV']),self.eps)
            extra=r['W_light_PV_early_D4_GeV2']-r['W_light_PV_pole_GeV2']
            self.assertLess(abs(extra-g2*n*self.w2/(144*mp.pi**2)),self.eps)

    def test_coupling_projection_and_archived_benchmark(self):
        r=self.rows[0]
        direct=(mp.re(r['Z_full_closed_fermion_pole_GeV2'])-mp.re(r['W_full_closed_fermion_pole_GeV2']))/(self.z2-self.w2)-r['W_full_closed_fermion_zero_GeV2']/self.w2
        self.assertLess(abs(direct-r['A29_dgprime2_relative']),self.eps)
        self.assertLess(abs(r['A28_dg2_relative']-mp.mpf('0.01322393860552000651')),mp.mpf('1e-20'))
        self.assertLess(abs(r['A29_dgprime2_relative']-mp.mpf('0.05007569924190788834')),mp.mpf('1e-20'))

    def test_tadpole_null_directions_and_nonzero_vev_direction(self):
        a=projections(self.w2,self.z2,mp.mpf(12),mp.mpc(5,3),mp.mpc(8,2))
        for frac in [mp.mpf('0.004'),mp.mpf('-0.002')]:
            b=projections(self.w2,self.z2,12+frac*self.w2,mp.mpc(5,3)+frac*self.w2,mp.mpc(8,2)+frac*self.z2)
            self.assertLess(abs(b['A28_dg2_relative']-a['A28_dg2_relative']),self.eps)
            self.assertLess(abs(b['A29_dgprime2_relative']-a['A29_dgprime2_relative']),self.eps)
            self.assertLess(abs(b['dv2_relative_from_W_input']-a['dv2_relative_from_W_input']-frac),self.eps)

    def test_scale_log_and_fermion_rg_from_Weyl_charges(self):
        b2=Fraction(2,3)*3*(3+1)*Fraction(1,2)
        fields=[(6,Fraction(1,6)),(3,Fraction(2,3)),(3,Fraction(-1,3)),(2,Fraction(-1,2)),(1,Fraction(-1))]
        by=Fraction(2,3)*3*sum(n*y*y for n,y in fields)
        for key,beta,coupling in [('A28_dg2_relative',b2,'g02'),('A29_dgprime2_relative',by,'gprime02')]:
            expected=mp.mpf(beta.numerator)/beta.denominator*self.out[coupling]/(16*mp.pi**2)
            first=self.rows[0]
            for r in self.rows:
                self.assertLess(abs(r['pole_coefficients_in_self_energy_projections'][key]-expected),self.eps)
                log=2*mp.log(mp.mpf(r['mu_bar_GeV'])/mp.mpf(first['mu_bar_GeV']))
                self.assertLess(abs(r[key]-first[key]-expected*log),self.eps)

    def test_dependencies_fail_closed(self):
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/'native-results.json'
            p.write_bytes(WTB.read_bytes()+b' ')
            with self.assertRaisesRegex(ValueError,'Dependency hash'):
                read_pinned(p,self.cfg['archive_dependencies']['wtb_native_sha256'])


if __name__=='__main__':
    unittest.main()
