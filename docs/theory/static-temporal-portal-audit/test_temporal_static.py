import math
import unittest
import mpmath as mp
from temporal_static import symbolic, epsilon_limit, holonomy_quad, holonomy_bessel, leading, next_mass_slope

class StaticWardAudit(unittest.TestCase):
 def test_symbolic(self):
  s=symbolic()
  for gauge in ('SU2','U1'):
   self.assertEqual(s[gauge]['shift_candidate_difference'],'0')
   self.assertEqual(s[gauge]['RMT_x_over_derived_x'],'2')
   self.assertEqual(s[gauge]['RMT_y_over_derived_y'],'1')
   self.assertEqual(s[gauge]['Schicho_a1_squared_over_derived'],'1')
  self.assertEqual(s['SU2_U1_mixed_trace'],'0')
 def test_evanescent_loop(self):
  with mp.workdps(70):
   lim=epsilon_limit(mp.mpf('1e-25'))
   self.assertLess(abs(lim/(mp.mpf(1)/(8*mp.pi**2))-1),mp.mpf('1e-22'))
 def test_independent_massive_holonomy(self):
  for r in (0.02,0.1,0.5,1.,2.,4.):
   q=holonomy_quad(r);b=holonomy_bessel(r)
   self.assertAlmostEqual(q,b,places=12)
  h=.02
  self.assertLess(abs((holonomy_quad(h)-leading())/h**2 /next_mass_slope()-1),1e-4)
