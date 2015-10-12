'''
Alignment fails with atoms w/o coordinates in a state
'''

import pymol
from pymol import cmd, testing

@testing.requires_version('1.7.7')
class Test2597(testing.PyMOLTestCase):

    def testAlignMissingCoords(self):
        filename = self.datafile('1t46-frag.pdb')

        cmd.load(filename, 'm1')
        cmd.remove('m1 & resi 600-605')
        cmd.load(filename, 'm1', state=2)

        cmd.load(filename, 'm2')
        cmd.remove('m2 & resi 620-625')
        cmd.load(filename, 'm2', state=2)

        cmd.align('m1', 'm2', object='aln', mobile_state=1, target_state=1, cycles=0)

        self.assertIn('aln', cmd.get_names())
