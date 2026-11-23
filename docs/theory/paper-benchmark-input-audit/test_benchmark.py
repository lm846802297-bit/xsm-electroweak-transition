import unittest
import sympy as S
from benchmark import setup, compute, tensor_rows, engineering3, a1,a2,b1,b3,b4,cH,b2,lam


class PaperInputContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=compute()

    def test_fixed_rotation_and_mislabelling(self):
        r=self.r
        self.assertEqual(r["physically_light_state_doublet_purity_journal"],S.Rational(".99"))
        self.assertEqual(r["physically_light_state_doublet_purity_arxiv_literal"],S.Rational(".01"))
        self.assertEqual(S.simplify(r["journal_diag_with_eq_2_3"][0,0]),S.Rational("125.1")**2)
        self.assertEqual(S.simplify(r["arxiv_v1_diag_with_eq_2_3"][0,0]),350**2)

    def test_tensor_multiplicities(self):
        _,_,_,p,_=setup()
        expected=lambda d:{"G1_s":d[b1],"G2_HH":d[cH],"G2_SS":d[b2],
               "G3_HHS":d[a1]/2,"G3_SSS":2*d[b3],"G4_HHHH":6*d[lam],
               "G4_HHKK":2*d[lam],"G4_HHSS":d[a2],"G4_SSSS":6*d[b4],"G3_HKS_zero":0}
        for q in (p,engineering3(p,S.Rational("82.12"))):
            self.assertEqual(tensor_rows(q),expected(q))

    def test_bare_positivity_domain(self):
        z=self.r["H0_bare_positive_domain"]
        self.assertLess(z["doublet_zeroes_s"][0],0)
        self.assertGreater(z["doublet_zeroes_s"][1],0)
        self.assertLess(z["doublet_curvature_at_s0"],0)
        _,_,_,p,_=setup()
        self.assertEqual(z["singlet_min_curvature"],p[b2]-p[b3]**2/(3*p[b4]))
        self.assertGreater(z["singlet_min_curvature"],0)

    def test_portal_sign_and_branch(self):
        _,_,_,p,q=setup()
        self.assertLess(p[a1],0)
        self.assertEqual(S.simplify(q[a1]+p[a1]),0)
        self.assertEqual(p[b3],q[b3])  # not a full S-reflection of the tree action


if __name__=="__main__":unittest.main()
