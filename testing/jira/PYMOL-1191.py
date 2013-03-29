'''
PYMOL-1191
MOE 2012 surfaces
'''

from pymol import cmd, testing, stored

@testing.requires('incentive')
class TestPYMOL1191(testing.PyMOLTestCase):

    def testLoadMOE2012(self):
        cmd.load(self.datafile('benzene.moe.gz'))

        v = cmd.get_names()
        self.assertTrue('benzene.system' in v)
        self.assertTrue('benzene.Quick_Surface_Ligand_' in v)
