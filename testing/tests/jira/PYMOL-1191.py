'''
PYMOL-1191
MOE 2012 surfaces
'''

from pymol import cmd, testing, stored

@testing.requires('incentive')
@testing.requires_version('1.7.4')
class TestPYMOL1191(testing.PyMOLTestCase):

    def testLoadMOE2012(self):
        cmd.load(self.datafile('benzene.moe.gz'))

        v = cmd.get_names()
        self.assertTrue('benzene.system' not in v)
        self.assertTrue('benzene.chain_1' in v)
        self.assertTrue('benzene.Quick_Surface_Ligand_' in v)
